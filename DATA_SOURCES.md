# Data Sources

## Protein-DNA structural cognate sources

RCSB PDB entries are used as structural cognate protein-DNA cases. A PDB
structure is not treated as quantitative specificity ground truth in benchmark
v0.2.

| source_type | accession | url | retrieval_date | notes |
|---|---|---|---|---|
| PDB | 8TAC | https://www.rcsb.org/structure/8TAC | 2026-08-29 | Designed DNA binding protein |
| PDB | 1D66 | https://www.rcsb.org/structure/1D66 | 2026-08-29 | GAL4 protein-DNA complex |
| PDB | 1TUP | https://www.rcsb.org/structure/1TUP | 2026-08-29 | p53-DNA complex |
| PDB | 4KPY | https://www.rcsb.org/structure/4KPY | 2026-08-29 | TtAgo/DNA complex |
| PDB | 8XIH | https://www.rcsb.org/structure/8XIH | 2026-08-29 | protein-DNA complex |
| PDB | 1O3Q | https://www.rcsb.org/structure/1O3Q | 2026-08-29 | CAP-DNA complex |
| PDB | 1O3R | https://www.rcsb.org/structure/1O3R | 2026-08-29 | CAP-DNA complex |
| PDB | 1O3S | https://www.rcsb.org/structure/1O3S | 2026-08-29 | CAP-DNA complex |
| PDB | 1O3T | https://www.rcsb.org/structure/1O3T | 2026-08-29 | CAP-DNA complex |
| PDB | 4JBM | https://www.rcsb.org/structure/4JBM | 2026-08-29 | murine DNA binding protein complex |
| PDB | 1BNZ | https://www.rcsb.org/structure/1BNZ | 2026-08-29 | SSO7D protein/DNA complex |
| PDB | 9C0F | https://www.rcsb.org/structure/9C0F | 2026-08-29 | piggyBat transposase protein-DNA complex |
| PDB | 1WVL | https://www.rcsb.org/structure/1WVL | 2026-08-29 | Sac7d-GCN4 with DNA decamer |
| PDB | 4E54 | https://www.rcsb.org/structure/4E54 | 2026-08-29 | UV-DDB complex with DNA |
| PDB | 4E5Z | https://www.rcsb.org/structure/4E5Z | 2026-08-29 | UV-DDB complex with DNA |
| PDB | 1B3T | https://www.rcsb.org/structure/1B3T | 2026-08-29 | EBNA-1 protein/DNA complex |

## Genome source

| source_type | accession | url | retrieval_date | notes |
|---|---|---|---|---|
| Genome FASTA | GRCh38 chr22 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | 2026-08-29 | Genome candidate retrieval demo chromosome |

## Experimental specificity pilot sources

Layer C v0.2 used a small JASPAR CORE pilot. Scores in
`data/processed/experimental_specificity_small.csv` are PWM log2-odds values
derived from JASPAR position frequency matrices. They are not raw PBM
enrichment scores and are not cross-assay normalized.

| source_type | accession | protein | UniProt | URL | retrieval_date | notes |
|---|---|---|---|---|---|---|
| JASPAR CORE PFM | MA0139.1 | CTCF | P49711 | https://jaspar.elixir.no/matrix/MA0139.1/ | 2026-08-29 | PFM-derived PWM score pilot |
| JASPAR CORE PFM | MA0106.3 | TP53 | P04637 | https://jaspar.elixir.no/matrix/MA0106.3/ | 2026-08-29 | PFM-derived PWM score pilot |
| JASPAR CORE PFM | MA0493.1 | Klf1 | P46099 | https://jaspar.elixir.no/matrix/MA0493.1/ | 2026-08-29 | PFM-derived PWM score pilot |
| JASPAR CORE PFM | MA0035.4 | GATA1 | P15976 | https://jaspar.elixir.no/matrix/MA0035.4/ | 2026-08-29 | PFM-derived PWM score pilot |
| JASPAR CORE PFM | MA0079.5 | SP1 | P08047 | https://jaspar.elixir.no/matrix/MA0079.5/ | 2026-08-29 | PFM-derived PWM score pilot |
| UniProt FASTA | P49711 | CTCF | P49711 | https://rest.uniprot.org/uniprotkb/P49711.fasta | 2026-08-29 | Protein sequence provenance |
| UniProt FASTA | P04637 | TP53 | P04637 | https://rest.uniprot.org/uniprotkb/P04637.fasta | 2026-08-29 | Protein sequence provenance |
| UniProt FASTA | P46099 | Klf1 | P46099 | https://rest.uniprot.org/uniprotkb/P46099.fasta | 2026-08-29 | Protein sequence provenance |
| UniProt FASTA | P15976 | GATA1 | P15976 | https://rest.uniprot.org/uniprotkb/P15976.fasta | 2026-08-29 | Protein sequence provenance |
| UniProt FASTA | P08047 | SP1 | P08047 | https://rest.uniprot.org/uniprotkb/P08047.fasta | 2026-08-29 | Protein sequence provenance |

## Designed DBP uPBM experimental specificity sources

Layer C v0.3 uses GEO GSE237017 processed and raw uPBM supplementary files for designed DNA-binding proteins. These files provide experimental 7-mer PBM E-scores, median intensities, and z-scores for DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48.

| source_type | accession | URL | retrieval_date | local_manifest | notes |
|---|---|---|---|---|---|
| GEO Series | GSE237017 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237017 | 2026-09-01 | `metadata/v0_3/gse237017_series_metadata.json` | Series metadata for computational design of sequence-specific DNA-binding proteins |
| GEO Family SOFT | GSE237017_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237017/soft/GSE237017_family.soft.gz | 2026-09-01 | `metadata/v0_3/gse237017_samples.csv` | Programmatically parsed GSM sample metadata |
| GEO Supplementary Files | 12 processed 7-mer files and 12 raw spot-data files | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237017/suppl/ | 2026-09-01 | `metadata/v0_3/gse237017_file_manifest.csv` | Raw files preserved under `data/raw/gse237017/` with SHA256 and file size |
| Nature Supplementary Workbook | 41594_2025_1669_MOESM3_ESM.xlsx | https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM3_ESM.xlsx | 2026-09-01 | `metadata/v0_3/designed_dbp_design_source_manifest.csv` | Official supplementary tables used for designed DBP protein sequences and intended target DNA sequences |
| Nature Source Data Fig. 4 | 41594_2025_1669_MOESM12_ESM.xlsx | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM12_ESM.xlsx | 2026-09-01 | `metadata/v0_3_1/paper_source_data_manifest.csv` | Source data used to cross-check DBP48 sequence C context |
| Nature Source Data Extended Data Fig. 8 | 41594_2025_1669_MOESM20_ESM.xls | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM20_ESM.xls | 2026-09-01 | `metadata/v0_3_1/paper_source_data_manifest.csv` | Source data used to reproduce published uPBM motif percentile values |

The v0.3 benchmark output is `data/processed/v0_3/designed_dbp_upbm_v0_3.parquet`. The v0.3.1 corrected outputs are under `data/processed/v0_3_1/`, including oriented rows and RC-class units. PBM E-score is treated as a processed experimental 7-mer specificity/enrichment score for per-protein ranking, not as Kd or absolute cross-protein affinity.

## Literature anchor

The designed binder example 8TAC is documented on the RCSB page and linked to the 2025 Nature article on computational design of sequence-specific DNA-binding proteins.
