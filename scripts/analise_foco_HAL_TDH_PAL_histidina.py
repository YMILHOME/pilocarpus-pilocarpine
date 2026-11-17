# Script drafted with AI assistance (Claude, Anthropic); reviewed and validated by the authors.

"""
Análise focada: HAL, TDH (TAL), PAL e via de biossíntese da histidina (HISN1–8)
================================================================================
Pilocarpus microphyllus – transcriptomas de raiz, caule, folha e ráquis.

Para cada enzima e tecido:
  * eggNOG-mapper  -> busca pelo número EC (coluna EC)
  * Mercator4      -> busca pelo nome do bin
  * origem         -> grupo taxonômico do eggNOG (max_annot_lvl): planta x microrganismo

Saídas (pasta resultados_foco/):
  tabela_foco.csv            contagens por enzima x tecido (eggNOG, Mercator, ambos, união, planta)
  candidatos_foco.csv        todos os transcritos candidatos, com origem e anotação (p/ filogenia)
  fig_foco_heatmap.(png|svg) heatmap: enzimas x tecidos, "eggNOG | Mercator" em cada célula
  fig_foco_origem.(png|svg)  origem taxonômica dos candidatos: raiz x tecidos aéreos

Uso no Colab:  !python analise_foco_HAL_TDH_PAL_histidina.py
"""

import io
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================== CONFIGURAÇÃO ==============================
# "aba" = aba do Excel com a anotação COMPLETA do eggNOG
AMOSTRAS = {
    "Root":   {"eggnog": "root_eggnog.xlsx",   "aba": "total_RnaSeq_anot", "mercator": "root_mercator.txt"},
    "Stem":   {"eggnog": "stem_eggnog.xlsx",   "aba": "Raw_Stem",          "mercator": "stem_mercator.txt"},
    "Leaf":   {"eggnog": "leaf_eggnog.xlsx",   "aba": "filtro1",           "mercator": "leaf_mercator.txt"},
    "Rachis": {"eggnog": "rachis_eggnog.xlsx", "aba": "Filtro1",           "mercator": "rachis_mercator.txt"},
}
TECIDO_MISTO = "Root"   # tecido que ganha a coluna extra "plant only"

# grupo, rótulo, ECs (eggNOG), KOs do KEGG (eggNOG), regex do bin (Mercator; None = sem bin)
# O eggNOG conta o transcrito se tiver o EC OU o KO (a IGP sintase, por ex., vem sem EC).
HB = r"histidine biosynthesis\."
ALVOS = [
    ("Histidine biosynthesis", "HISN1 · ATP-PRT",          ["2.4.2.17"],             ["K00765"],
     HB + r"ATP phosphoribosyl ?transferase"),
    ("Histidine biosynthesis", "HISN2 · PRA-PH/CH",        ["3.6.1.31", "3.5.4.19"], ["K11755", "K01523", "K01496"],
     HB + r"bifunctional phosphoribosyl-ATP"),
    ("Histidine biosynthesis", "HISN3 · ProFAR isomerase", ["5.3.1.16"],             ["K01814", "K24017"],
     HB + r"phosphoribosyl-formimino"),
    ("Histidine biosynthesis", "HISN4 · IGP synthase",     ["4.3.2.10"],             ["K01663", "K02501", "K02500"],
     HB + r"imidazoleglycerol-phosphate synthase"),
    ("Histidine biosynthesis", "HISN5 · IGP dehydratase",  ["4.2.1.19"],             ["K01693", "K01089"],
     HB + r"imidazoleglycerol-phosphate dehydratase"),
    ("Histidine biosynthesis", "HISN6 · HolP aminotransf.", ["2.6.1.9"],             ["K00817"],
     HB + r"histidinol-phosphate aminotransferase"),
    ("Histidine biosynthesis", "HISN7 · HolP phosphatase", ["3.1.3.15"],             ["K04486", "K05602", "K01089"],
     HB + r"histidinol-phosphate phosphatase"),
    ("Histidine biosynthesis", "HISN8 · histidinol DH",    ["1.1.1.23"],             ["K00013"],
     HB + r"histidinol dehydrogenase"),
    ("Ammonia-lyases",         "HAL",                      ["4.3.1.3"],              ["K01745"], None),
    ("Ammonia-lyases",         "PAL",                      ["4.3.1.24", "4.3.1.25"], ["K10775", "K13064"],
     r"phenylalanine ammonia.lyase \*?\(PAL\)"),
    ("Ammonia-lyases",         "TDH (TAL)",                ["4.3.1.19"],             ["K01754", "K17989"],
     r"threonine degradation\.threonine dehydratase"),
]
# Anotações do eggNOG sem EC que também contam (só HAL: 1 transcrito "Histidine ammonia-lyase")
TEXTO_EGGNOG = {"HAL": r"histidine ammonia.lyase"}

USAR_KO = True   # False = só EC

PASTA_SAIDA = Path("resultados_foco")
CORES = {"navy": "#13284B", "orange": "#F07D26", "gray": "#6B7A6E", "green": "#3A8D3F",
         "light": "#B9C2CC"}
# ==========================================================================


def normaliza_id(x):
    return str(x).strip().strip("'").lower()


def le_eggnog(caminho, aba=None):
    if str(caminho).lower().endswith((".xlsx", ".xls")):
        bruto = pd.read_excel(caminho, sheet_name=aba if aba else 0, header=None, dtype=str).fillna("")
        i_cab = next(i for i, r in bruto.iterrows() if any(str(v).strip().upper() == "EC" for v in r))
        df = bruto.iloc[i_cab + 1:].copy()
        df.columns = [str(c) for c in bruto.iloc[i_cab]]
        df = df[~df.iloc[:, 0].astype(str).str.startswith("##")]
        df = df[df.iloc[:, 0].astype(str).str.strip() != ""].replace("", "-")
    else:
        linhas = [l for l in open(caminho, encoding="utf-8-sig", errors="replace") if not l.startswith("##")]
        sep = max(["\t", ";", ","], key=linhas[0].count)
        df = pd.read_csv(io.StringIO("".join(linhas)), sep=sep, dtype=str).fillna("-")
    df.columns = [c.strip().lstrip("#").strip() for c in df.columns]
    col_id = next((c for c in df.columns if c.lower() in ("query", "query_name")), df.columns[0])
    df["id"] = df[col_id].map(normaliza_id)
    for c in ("EC", "KEGG_ko", "Description", "Preferred_name", "max_annot_lvl", "seed_ortholog"):
        df[c] = df[c].fillna("-").astype(str) if c in df.columns else "-"
    df["ec_set"] = df["EC"].map(lambda c: {e.strip() for e in c.split(",")})
    df["ko_set"] = df["KEGG_ko"].map(lambda c: {k.strip().replace("ko:", "") for k in c.split(",")})
    return df


def le_mercator(caminho):
    df = pd.read_csv(caminho, sep="\t", dtype=str, quotechar="'").fillna("")
    df.columns = [c.strip().strip("'").upper() for c in df.columns]
    for c in ("NAME", "IDENTIFIER"):
        df[c] = df[c].str.strip().str.strip("'")
    df = df[df["IDENTIFIER"] != ""].copy()
    df["id"] = df["IDENTIFIER"].map(normaliza_id)
    return df


def grupo_taxonomico(nivel):
    n = str(nivel).split("|")[-1].strip().lower()
    if n in ("streptophyta", "viridiplantae", "embryophyta", "magnoliopsida", "chlorophyta"):
        return "Plant"
    if n in ("", "-", "nan"):
        return "Not in eggNOG"
    if n == "eukaryota":
        return "Eukaryota (unresolved)"
    if n in ("metazoa", "fungi", "opisthokonta", "amoebozoa", "peronosporales", "pythiales",
             "oomycota", "stramenopiles", "alveolata", "apicomplexa", "nematoda", "arthropoda"):
        return "Other eukaryotes"
    return "Bacteria/Archaea"


def coleta():
    regs = []
    for tecido, arq in AMOSTRAS.items():
        print(f"Lendo {tecido}...")
        egg = le_eggnog(arq["eggnog"], arq.get("aba")).set_index("id", drop=False)
        merc = le_mercator(arq["mercator"])
        for grupo, rot, ecs, kos, rx in ALVOS:
            sel = egg["ec_set"].map(lambda s: any(e in s for e in ecs))
            if USAR_KO:
                sel |= egg["ko_set"].map(lambda s: any(k in s for k in kos))
            chave = rot.split(" ")[0]
            if chave in TEXTO_EGGNOG:
                sel |= egg["Description"].str.contains(TEXTO_EGGNOG[chave], case=False, regex=True)
            ids_e = set(egg.index[sel])
            if rx:
                m = merc[merc["NAME"].str.contains(rx, flags=re.I, regex=True)]
                ids_m = set(m["id"])
                bin_m = dict(zip(m["id"], m["NAME"]))
            else:
                ids_m, bin_m = set(), {}
            for i in ids_e | ids_m:
                info = egg.loc[i] if i in egg.index else None
                if isinstance(info, pd.DataFrame):
                    info = info.iloc[0]
                regs.append({
                    "grupo": grupo, "enzima": rot, "tecido": tecido, "transcrito": i,
                    "eggnog": i in ids_e, "mercator": i in ids_m,
                    "origem": grupo_taxonomico(info["max_annot_lvl"] if info is not None else ""),
                    "max_annot_lvl": info["max_annot_lvl"] if info is not None else "",
                    "EC": info["EC"] if info is not None else "",
                    "KEGG_ko": info["KEGG_ko"] if info is not None else "",
                    "seed_ortholog": info["seed_ortholog"] if info is not None else "",
                    "eggnog_description": info["Description"] if info is not None else "",
                    "mercator_bin": bin_m.get(i, ""),
                })
    return pd.DataFrame(regs)


def tabela(c):
    enz = [a[1] for a in ALVOS]
    linhas = []
    for (e, t), g in c.groupby(["enzima", "tecido"]):
        linhas.append({"enzima": e, "tecido": t,
                       "eggNOG": int(g.eggnog.sum()), "Mercator": int(g.mercator.sum()),
                       "ambos": int((g.eggnog & g.mercator).sum()), "uniao": len(g),
                       "planta_eggNOG": int((g.eggnog & (g.origem == "Plant")).sum()),
                       "planta_Mercator": int((g.mercator & (g.origem == "Plant")).sum()),
                       "planta_uniao": int((g.origem == "Plant").sum())})
    t = pd.DataFrame(linhas)
    idx = pd.MultiIndex.from_product([enz, list(AMOSTRAS)], names=["enzima", "tecido"])
    return t.set_index(["enzima", "tecido"]).reindex(idx, fill_value=0).reset_index()


def fig_heatmap(t):
    enz = [a[1] for a in ALVOS]
    cols = []
    for tec in AMOSTRAS:
        cols.append((tec, False))
        if tec == TECIDO_MISTO:
            cols.append((tec, True))
    rot_cols = [f"{tec}\n(plant only)" if p else tec for tec, p in cols]
    U = np.zeros((len(enz), len(cols)))
    txt = [["" for _ in cols] for _ in enz]
    for r, e in enumerate(enz):
        for c, (tec, p) in enumerate(cols):
            row = t[(t.enzima == e) & (t.tecido == tec)].iloc[0]
            eg, me = (row.planta_eggNOG, row.planta_Mercator) if p else (row.eggNOG, row.Mercator)
            U[r, c] = row.planta_uniao if p else row.uniao
            sem_bin = [a for a in ALVOS if a[1] == e][0][4] is None
            txt[r][c] = f"{eg} | –" if sem_bin else f"{eg} | {me}"
    vmax = max(1, U.max())
    fig, ax = plt.subplots(figsize=(1.35 * len(cols) + 3.6, 0.48 * len(enz) + 1.9))
    im = ax.imshow(np.log1p(U), cmap="Blues", vmin=0, vmax=np.log1p(vmax), aspect="auto")
    for r in range(len(enz)):
        for c in range(len(cols)):
            ax.text(c, r, txt[r][c], ha="center", va="center", fontsize=9, fontweight="bold",
                    color="white" if np.log1p(U[r, c]) > np.log1p(vmax) * 0.6 else CORES["navy"])
    ax.set_xticks(range(len(cols)), rot_cols, fontsize=9)
    ax.set_yticks(range(len(enz)), enz, fontsize=9)
    n_his = sum(1 for a in ALVOS if a[0] == "Histidine biosynthesis")
    ax.axhline(n_his - 0.5, color=CORES["orange"], lw=2.2)
    ax.tick_params(length=0)
    ax.text(-0.33, (n_his - 1) / 2, "Histidine\nbiosynthesis", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=9, color=CORES["orange"], fontweight="bold")
    ax.text(-0.33, n_his + (len(enz) - n_his - 1) / 2, "Ammonia-\nlyases",
            transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=9,
            color=CORES["orange"], fontweight="bold")
    ax.set_title("Candidate transcripts per tissue  (eggNOG-mapper | Mercator4)",
                 color=CORES["navy"], fontweight="bold", fontsize=11)
    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    ticks = [v for v in (0, 1, 2, 5, 10, 20, 50, 100) if v <= vmax]
    cb.set_ticks(np.log1p(ticks), labels=[str(v) for v in ticks])
    cb.set_label("Transcripts (either tool)")
    fig.text(0.01, 0.005, "– : Mercator4 has no HAL bin. 'Plant only' = transcripts assigned to "
             "Streptophyta by eggNOG-mapper.", fontsize=7.5, color=CORES["gray"])
    salva(fig, "fig_foco_heatmap")


def fig_origem(c):
    ordem = ["Plant", "Eukaryota (unresolved)", "Other eukaryotes", "Bacteria/Archaea", "Not in eggNOG"]
    cores = [CORES["green"], CORES["light"], CORES["gray"], CORES["orange"], "#E3E6EA"]
    c = c.copy()
    c["bloco"] = c["enzima"].where(c["grupo"] != "Histidine biosynthesis", "Histidine pathway (HISN1–8)")
    c["local"] = np.where(c["tecido"] == TECIDO_MISTO, TECIDO_MISTO, "Aerial (stem+leaf+rachis)")
    blocos = ["Histidine pathway (HISN1–8)"] + [a[1] for a in ALVOS if a[0] == "Ammonia-lyases"]
    tab = (c.groupby(["bloco", "local"])["origem"].value_counts().unstack(fill_value=0)
           .reindex(columns=ordem, fill_value=0))
    idx = [(b, l) for b in blocos for l in (TECIDO_MISTO, "Aerial (stem+leaf+rachis)")]
    tab = tab.reindex(idx, fill_value=0)
    tab.to_csv(PASTA_SAIDA / "origem_foco.csv")
    tab = tab[tab.sum(axis=1) > 0]   # tira combinações sem transcritos
    pct = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0) * 100
    fig, ax = plt.subplots(figsize=(7.8, 0.42 * len(tab) + 1.6))
    rot = [f"{b} · {l} (n={int(tab.loc[(b, l)].sum())})" for b, l in tab.index]
    esq = np.zeros(len(tab))
    for col, cor in zip(ordem, cores):
        v = pct[col].fillna(0).values
        ax.barh(rot, v, left=esq, color=cor, label=col, height=0.7)
        for i, x in enumerate(v):
            if x >= 8:
                ax.text(esq[i] + x / 2, i, f"{x:.0f}%", ha="center", va="center", fontsize=8,
                        color=CORES["navy"] if col in ("Eukaryota (unresolved)", "Not in eggNOG") else "white")
        esq += v
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of candidate transcripts (taxonomy from eggNOG-mapper)")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    salva(fig, "fig_foco_origem")


def salva(fig, nome):
    for ext in ("png", "svg"):
        fig.savefig(PASTA_SAIDA / f"{nome}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def main():
    PASTA_SAIDA.mkdir(exist_ok=True)
    c = coleta()
    t = tabela(c)
    c.sort_values(["grupo", "enzima", "tecido", "origem"]).to_csv(PASTA_SAIDA / "candidatos_foco.csv", index=False)
    t.to_csv(PASTA_SAIDA / "tabela_foco.csv", index=False)
    fig_heatmap(t)
    fig_origem(c)
    print(t.to_string(index=False))
    print(f"\nArquivos em: {PASTA_SAIDA.resolve()}")


if __name__ == "__main__":
    main()
