#!/usr/bin/env python3
"""Gera graficos e tabelas a partir de dados/medicoes_<maquina>.csv (e perfil_<maquina>.csv).
Uso: python3 scripts/plotar.py            (processa todas as maquinas encontradas em dados/)
Saidas em graficos/: tempo, speedup (com curva ideal e Lei de Amdahl), eficiencia,
escalabilidade fraca, escalonamento, reducao vs atomic vs critical, overhead de comunicacao MPI,
e dados/resumo_<maquina>.csv (tabela de speedup/eficiencia) + dados/perfil_resumo_<maquina>.md.
Estatistica: mediana das repeticoes (a execucao de aquecimento ja e descartada em medir.sh).
"""
import glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("graficos", exist_ok=True)

def amdahl(f, p):
    return 1.0 / ((1.0 - f) + f / p)

def fracao_paralela(maq):
    """f = t_laco / (t_init + t_laco + t_soma), do perfil do baseline instrumentado."""
    arq = f"dados/perfil_{maq}.csv"
    if not os.path.exists(arq):
        return None, None
    pf = pd.read_csv(arq)
    g = pf.groupby("N")[["t_init", "t_laco", "t_soma"]].median()
    g["total"] = g.sum(axis=1)
    g["f"] = g["t_laco"] / g["total"]
    with open(f"dados/perfil_resumo_{maq}.md", "w") as fh:
        fh.write("| N | t_init (s) | t_laco (s) | t_soma (s) | % laco (f) | speedup max (Amdahl, p->inf) |\n|---|---|---|---|---|---|\n")
        for N, r in g.iterrows():
            fh.write(f"| {N} | {r.t_init:.6f} | {r.t_laco:.6f} | {r.t_soma:.6f} | {100*r.f:.3f}% | {1/(1-r.f):.0f} |\n")
    return g["f"], g

for arq in sorted(glob.glob("dados/medicoes_*.csv")):
    maq = os.path.basename(arq)[len("medicoes_"):-len(".csv")]
    df = pd.read_csv(arq)
    df["variante"] = df["variante"].fillna("-")
    med = (df.groupby(["tipo", "versao", "variante", "N", "T", "workers"])
             [["tempo", "t_soma", "comm_max", "comm_media"]].median().reset_index())
    f_por_N, _ = fracao_paralela(maq)
    pref = f"graficos/{maq}"

    # ---- escalabilidade forte -------------------------------------------------
    forte = med[med.tipo == "forte"]
    tseq = forte[forte.versao == "seq"].set_index("N")["tempo"]
    Ns = sorted(tseq.index)
    linhas = []
    for v in ("omp", "mpi"):
        for N in Ns:
            s = forte[(forte.versao == v) & (forte.N == N)].sort_values("workers")
            for _, r in s.iterrows():
                sp = tseq[N] / r.tempo
                linhas.append(dict(versao=v, N=N, workers=int(r.workers), tempo=r.tempo,
                                   speedup=sp, eficiencia=sp / r.workers))
    res = pd.DataFrame(linhas)
    if not res.empty:
        res.to_csv(f"dados/resumo_{maq}.csv", index=False)

        # tempo
        fig, ax = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4), squeeze=False)
        for j, N in enumerate(Ns):
            a = ax[0, j]
            for v, mk in (("omp", "o-"), ("mpi", "s-")):
                s = res[(res.versao == v) & (res.N == N)]
                a.plot(s.workers, s.tempo, mk, label={"omp": "OpenMP (threads)", "mpi": "MPI (processos)"}[v])
            ws = sorted(res.workers.unique())
            a.plot(ws, [tseq[N] / w for w in ws], "k--", label="ideal (T_seq/p)")
            a.axhline(tseq[N], color="gray", ls=":", label="baseline sequencial")
            a.set_xscale("log", base=2); a.set_yscale("log")
            a.set_xlabel("threads / processos"); a.set_ylabel("tempo do laco (s)"); a.set_title(f"N={N}")
            a.grid(True, which="both", alpha=.3)
        ax[0, 0].legend(fontsize=8)
        fig.suptitle(f"Tempo — escalabilidade forte ({maq})"); fig.tight_layout()
        fig.savefig(f"{pref}_tempo.png", dpi=140); plt.close(fig)

        # speedup e eficiencia
        for nome, col, ylabel in (("speedup", "speedup", "speedup = T_seq / T(p)"),
                                  ("eficiencia", "eficiencia", "eficiencia = speedup / p")):
            fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
            for a, v, tit in ((ax[0], "omp", "OpenMP"), (ax[1], "mpi", "MPI")):
                for N in Ns:
                    s = res[(res.versao == v) & (res.N == N)]
                    a.plot(s.workers, s[col], "o-", label=f"N={N}")
                ws = sorted(res.workers.unique())
                if col == "speedup":
                    a.plot(ws, ws, "k--", label="ideal")
                    if f_por_N is not None:
                        Nref = max(Ns)
                        if Nref in f_por_N.index:
                            a.plot(ws, [amdahl(f_por_N[Nref], w) for w in ws], ":", color="red",
                                   label=f"Amdahl (f={f_por_N[Nref]:.4f}, N={Nref})")
                else:
                    a.axhline(1.0, color="k", ls="--", label="ideal")
                a.set_xscale("log", base=2); a.set_xlabel("threads / processos"); a.set_ylabel(ylabel)
                a.set_title(tit); a.grid(True, alpha=.3); a.legend(fontsize=8)
            fig.suptitle(f"{col.capitalize()} — escalabilidade forte ({maq})"); fig.tight_layout()
            fig.savefig(f"{pref}_{nome}.png", dpi=140); plt.close(fig)

    # ---- escalabilidade fraca ---------------------------------------------------
    fraca = med[med.tipo == "fraca"]
    if not fraca.empty:
        fig, ax = plt.subplots(figsize=(6, 4.2))
        for v, mk in (("omp", "o-"), ("mpi", "s-")):
            s = fraca[fraca.versao == v].sort_values("workers")
            if s.empty: continue
            t1 = s[s.workers == 1]["tempo"].iloc[0]
            ax.plot(s.workers, t1 / s.tempo, mk, label={"omp": "OpenMP", "mpi": "MPI"}[v])
        ax.axhline(1.0, color="k", ls="--", label="ideal")
        ax.set_xscale("log", base=2); ax.set_ylim(0, 1.2)
        ax.set_xlabel("threads / processos (N = N0*sqrt(p))"); ax.set_ylabel("eficiencia fraca = T(1)/T(p)")
        ax.set_title(f"Escalabilidade fraca ({maq})"); ax.grid(True, alpha=.3); ax.legend()
        fig.tight_layout(); fig.savefig(f"{pref}_fraca.png", dpi=140); plt.close(fig)

    # ---- escalonamento ----------------------------------------------------------
    sc = med[med.tipo == "sched"]
    if not sc.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(sc.variante, sc.tempo)
        ax.set_ylabel("tempo do laco (s)")
        ax.set_title(f"schedule — N={int(sc.N.iloc[0])}, {int(sc.workers.iloc[0])} threads ({maq})")
        ax.grid(True, axis="y", alpha=.3); fig.tight_layout()
        fig.savefig(f"{pref}_schedule.png", dpi=140); plt.close(fig)

    # ---- reduction x atomic x critical ------------------------------------------
    rd = med[med.tipo == "reducao"]
    if not rd.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        rot = rd.versao.map({"omp": "reduction", "omp_atomic": "atomic", "omp_critical": "critical"})
        ax.bar(rot, rd.t_soma)
        ax.set_yscale("log"); ax.set_ylabel("tempo da soma (s, escala log)")
        ax.set_title(f"Protecao da soma — N={int(rd.N.iloc[0])} ({maq})")
        ax.grid(True, axis="y", alpha=.3); fig.tight_layout()
        fig.savefig(f"{pref}_reducao.png", dpi=140); plt.close(fig)

    # ---- overhead de comunicacao MPI --------------------------------------------
    mp = forte[forte.versao == "mpi"]
    if not mp.empty:
        fig, ax = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4), squeeze=False)
        for j, N in enumerate(Ns):
            s = mp[mp.N == N].sort_values("workers")
            frac = 100 * s.comm_media / s.tempo
            ax[0, j].bar(s.workers.astype(str), frac)
            ax[0, j].set_title(f"N={N}"); ax[0, j].set_xlabel("processos")
            ax[0, j].set_ylabel("% do tempo em comunicacao (media)"); ax[0, j].grid(True, axis="y", alpha=.3)
        fig.suptitle(f"Overhead de comunicacao MPI (troca de bordas) ({maq})"); fig.tight_layout()
        fig.savefig(f"{pref}_comunicacao.png", dpi=140); plt.close(fig)

    print(f"[{maq}] graficos em graficos/{maq}_*.png ; tabelas em dados/resumo_{maq}.csv")
