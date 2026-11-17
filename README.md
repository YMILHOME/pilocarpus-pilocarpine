# In search of candidate enzymes for pilocarpine biosynthesis in *Pilocarpus microphyllus*

Code and methods for the poster presented at the **Brazilian Symposium on Bioinformatics (BSB 2026)**, Ribeirão Preto, SP.

**Yasmim Aquino Milhome**, Bruno Costa Evangelista Souza, Alison Pelri Albuquerque Menezes, Acácio Raphael Aguiar Moraes, Fernando Araripe Gonçalves Torres3, Natalia Florencio Martins.
Universidade Federal do Ceará (UFC) · Embrapa Agroindústria Tropical · Universidade de Brasília (UnB)

> **Data availability.** *P. microphyllus* assemblies, predicted proteins, annotation tables and alignments are part of unpublished work and are **not** included. They are available from the authors upon reasonable request. Raw reads: NCBI BioProject **PRJNA1220646**.

---

## Overview

Pilocarpine is an imidazole alkaloid whose biosynthetic pathway in jaborandi (*Pilocarpus microphyllus*) is still unresolved. The working hypothesis is that the imidazole ring derives from **histidine** and the aliphatic moiety from **threonine**, with ammonia-lyases as candidate entry enzymes.

This repository contains the workflow used to:

1. identify candidate transcripts for **histidine biosynthesis (HISN1–8)**, **histidine ammonia-lyase (HAL)**, **phenylalanine ammonia-lyase (PAL)** and **threonine dehydratase (TDH)** in root, stem, leaf and rachis transcriptomes;
2. compare annotations from **eggNOG-mapper** and **Mercator4**;
3. assign the **taxonomic origin** (plant vs. non-plant) of each candidate;
4. place candidates in **maximum-likelihood phylogenies** (methods described below; trees not included).

## Pipeline

```
Raw reads (PRJNA1220646)
  └─ rnaSPAdes ............ de novo assembly (per tissue)
      └─ TransDecoder ..... ORF / protein prediction
          ├─ eggNOG-mapper  functional annotation (EC, KEGG KO, max_annot_lvl)
          └─ Mercator4 .... plant-specific functional bins
              └─ candidate screening (EC / KO / Mercator bin)  → scripts/
                  └─ MUSCLE → trimAl → IQ-TREE (BIC model, 1,000 UFBoot)
```

Most steps were run on the European Galaxy server, **[usegalaxy.eu](https://usegalaxy.eu)** (The Galaxy Community, 2026). Mercator4 was run on its web server, alignments were built with MUSCLE in MEGA 12, trimmed with trimAl and analysed with IQ-TREE on their web servers, and candidate screening and figures were done with the scripts in this repository.

| Step | Tool | Platform | Version | Main parameters |
|---|---|---|---|---|
| Assembly | rnaSPAdes | usegalaxy.eu | `TODO` | default |
| ORF prediction | TransDecoder | usegalaxy.eu | `TODO` | `TODO` |
| Annotation | eggNOG-mapper | usegalaxy.eu | `TODO` | `TODO` |
| Plant annotation | Mercator4 | Mercator4 web server | `TODO` | default |
| Candidate screening | Apps Script + Python | Google Sheets / Colab | — | see `scripts/` |
| Alignment | MUSCLE | MEGA 12 | MEGA 12 | default |
| Trimming | trimAl | trimAl web server (`TODO` URL) | 1.5.1 | `TODO` |
| Phylogeny | IQ-TREE | W-IQ-TREE web server | 2.4 | model selection by BIC, 1,000 UFBoot |

The Galaxy workflow (tools and parameters, without data) is available at: `TODO: link to the published Galaxy workflow`.

## Target enzymes

| Enzyme | EC | KEGG KO | Mercator4 bin (regex) |
|---|---|---|---|
| HISN1 · ATP phosphoribosyltransferase | 2.4.2.17 | K00765 | `histidine biosynthesis.ATP phosphoribosyltransferase` |
| HISN2 · PRA-PH / PRA-CH | 3.6.1.31, 3.5.4.19 | K11755, K01523, K01496 | `histidine biosynthesis.bifunctional phosphoribosyl-ATP…` |
| HISN3 · ProFAR isomerase | 5.3.1.16 | K01814, K24017 | `histidine biosynthesis.phosphoribosyl-formimino…` |
| HISN4 · IGP synthase | 4.3.2.10 | K01663, K02501, K02500 | `histidine biosynthesis.imidazoleglycerol-phosphate synthase` |
| HISN5 · IGP dehydratase | 4.2.1.19 | K01693, K01089 | `histidine biosynthesis.imidazoleglycerol-phosphate dehydratase` |
| HISN6 · histidinol-phosphate aminotransferase | 2.6.1.9 | K00817 | `histidine biosynthesis.histidinol-phosphate aminotransferase` |
| HISN7 · histidinol-phosphate phosphatase | 3.1.3.15 | K04486, K05602, K01089 | `histidine biosynthesis.histidinol-phosphate phosphatase` |
| HISN8 · histidinol dehydrogenase | 1.1.1.23 | K00013 | `histidine biosynthesis.histidinol dehydrogenase` |
| HAL · histidine ammonia-lyase | 4.3.1.3 | K01745 | — (no Mercator4 bin) |
| PAL · phenylalanine ammonia-lyase | 4.3.1.24, 4.3.1.25 | K10775, K13064 | `phenylalanine ammonia lyase (PAL)` |
| TDH · threonine dehydratase | 4.3.1.19 | K01754, K17989 | `threonine degradation.threonine dehydratase` |

A transcript is counted by eggNOG-mapper if it carries the **EC or the KO**; by Mercator4 if it falls in the bin. HAL also includes transcripts whose eggNOG description is "histidine ammonia-lyase" without an EC.

Taxonomic origin is taken from eggNOG-mapper `max_annot_lvl` and grouped as *Plant* (Streptophyta/Viridiplantae), *Bacteria/Archaea*, *Other eukaryotes* or *Eukaryota (unresolved)*. This is an annotation-based assignment, not a contamination test.


## Usage

Requirements: Python ≥ 3.9, `pandas`, `numpy`, `matplotlib`, `openpyxl` (all available on Google Colab).

1. Edit the `CONFIGURAÇÃO` block at the top of `analise_foco_HAL_TDH_PAL_histidina.py` with your file names (eggNOG-mapper output as `.xlsx`, `.csv` or `.tsv`; Mercator4 results as `.txt`) and, for Excel files, the sheet with the full annotation.
2. Run:
   ```bash
   python scripts/analise_foco_HAL_TDH_PAL_histidina.py
   python scripts/figuras_poster.py
   ```
3. Outputs are written to `resultados_foco/`: count tables, candidate list, taxonomic origin and figures (PNG, SVG, PDF).

## References

- Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics (Oxford, England), 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348
- Thomas KF Wong, Nhan Ly-Trong, Huaiyan Ren, Piyumal Demotte, Hector Banos, Andrew J Roger, Edward Susko, Chris Bielow, Nicola De Maio, Nick Goldman, Matthew W Hahn, Mario dos Reis, Le Sy Vinh, Gavin Huttley, Robert Lanfear, Bui Quang Minh (2026) IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models. Mol. Biol. Eviol.. https://doi.org/10.1093/molbev/msag117
- The Galaxy Community. Galaxy for accessible, reproducible, and collaborative data analyses: 2026 update. Nucleic Acids Research, 2026; gkag469, https://doi.org/10.1093/nar/gkag469

## Acknowledgments

We acknowledge the European Galaxy server (usegalaxy.eu), operated by the University of Freiburg as part of de.NBI / ELIXIR-DE, which was used for most of the analyses.

## Citation

Milhome, Y. A. *et al.* In search of candidate enzymes for pilocarpine biosynthesis in *Pilocarpus microphyllus* Stapf. Poster, Brazilian Symposium on Bioinformatics (BSB 2026), Ribeirão Preto, SP.

## Contact

Yasmim Aquino Milhome · ymilhome@gmail.com
