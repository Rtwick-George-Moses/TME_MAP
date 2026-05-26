# TCGA T Cell Exhaustion & TME Suppression Explorer

A Streamlit app that analyzes how the tumor microenvironment (TME) suppresses T cell immune responses across 15 TCGA cancer types. It measures which immune checkpoint pathways are being actively exploited by tumors, how strongly each is upregulated compared to healthy tissue, and which pathways the TME co-activates together — directly informing immunotherapy strategy.

---

## What it does

For a given cancer type, this tool answers three questions:

1. **Which inhibitory pathways is the TME engaging?** — By measuring the expression of 23 TME ligands that bind 16 T cell inhibitory receptors, compared against GTEx normal tissue baselines.

2. **How much is each pathway upregulated?** — All values are reported as log₂ fold-change over normal tissue (e.g., log₂FC = 3 means the ligand is 8× higher in tumor than in healthy tissue of the same organ).

3. **Which pathways does the TME co-activate together?** — A network graph shows Spearman correlations between pathway activation scores across patients. A thick blue edge between PD-1 and TIGIT means tumors that upregulate PD-L1/PD-L2 also upregulate CD155/CD112 — supporting combination anti-PD-1 + anti-TIGIT therapy.

---

## Receptors & Ligands Tracked

**16 inhibitory receptors** on T cells, across 6 families:

| Family | Receptors |
|--------|-----------|
| Classical Checkpoints | PD-1, CTLA-4, LAG-3, TIM-3, TIGIT |
| Co-inhibitory Ig Superfamily | BTLA, CD160, VISTA, CD200R |
| SLAM Family | 2B4 (CD244) |
| NK-related Receptors | KLRG1, NKG2A, CD96 |
| Metabolic Exhaustion Markers | CD39 |
| Additional Inhibitory Receptors | LAIR-1, ILT2 (LILRB1) |

**23 TME ligands** including PD-L1, PD-L2, B7-1/B7-2, FGL1, Galectin-3, Galectin-9, HMGB1, CD155, CD112, HVEM, HLA-E, HLA-G, E-cadherin, COL1A1, CD48, CD73, A2A receptor, VSIG3, CD200, IDO1, B7-H3, B7-H4.

---

## Data Sources

- **Tumor expression:** [TCGA](https://portal.gdc.cancer.gov/) via GDC API — STAR-Counts TPM from RNA-seq, 15 cancer types, up to 500 samples each.
- **Normal tissue baseline:** [GTEx](https://gtexportal.org/) V2 API — median TPM per gene in matched healthy tissue (e.g., GTEx Lung for TCGA-LUAD, GTEx Breast for TCGA-BRCA).
- **Demographics:** GDC cases endpoint — race, ethnicity, gender, AJCC pathologic stage.

### TCGA → GTEx Tissue Mapping

| TCGA Project | Cancer Type | GTEx Normal Tissue |
|---|---|---|
| TCGA-BRCA | Breast Invasive Carcinoma | Breast Mammary Tissue |
| TCGA-LUAD | Lung Adenocarcinoma | Lung |
| TCGA-LUSC | Lung Squamous Cell Carcinoma | Lung |
| TCGA-SKCM | Skin Cutaneous Melanoma | Skin Sun Exposed Lower leg |
| TCGA-COAD | Colon Adenocarcinoma | Colon Transverse |
| TCGA-BLCA | Bladder Urothelial Carcinoma | Bladder |
| TCGA-HNSC | Head & Neck Squamous Cell Carcinoma | Minor Salivary Gland |
| TCGA-KIRC | Kidney Renal Clear Cell Carcinoma | Kidney Cortex |
| TCGA-LIHC | Liver Hepatocellular Carcinoma | Liver |
| TCGA-PRAD | Prostate Adenocarcinoma | Prostate |
| TCGA-GBM | Glioblastoma Multiforme | Brain Cortex |
| TCGA-OV | Ovarian Serous Cystadenocarcinoma | Ovary |
| TCGA-UCEC | Uterine Corpus Endometrial Carcinoma | Uterus |
| TCGA-PAAD | Pancreatic Adenocarcinoma | Pancreas |
| TCGA-STAD | Stomach Adenocarcinoma | Stomach |

---

## Installation & Usage

### Requirements

- Python 3.9+
- ~500 MB disk space for the SQLite database (if using offline mode)

### Setup

```bash
git clone <repo-url>
cd tcga_app
pip install -r requirements.txt
```

### Option A: Offline mode (recommended)

Pre-download all TCGA data into a local SQLite database:

```bash
python download_to_db.py                          # all 15 projects, 300 samples each
python download_to_db.py --max-samples 500         # more samples
python download_to_db.py --projects TCGA-SKCM TCGA-BRCA  # specific projects
```

Then run the app:

```bash
streamlit run Explorer.py
```

The Explorer checks for `tcga_data.db` on startup. If found, data loads instantly. GTEx baselines are fetched on first load per project and cached permanently by Streamlit.

### Option B: Live mode

Run without downloading — the app fetches from GDC API on the fly:

```bash
streamlit run Explorer.py
```

This is slower (each project takes 30-60 seconds to download) but requires no setup beyond `pip install`.

---

## Project Structure

```
tcga_app/
├── Explorer.py              # Main Streamlit app (network graph, charts, filtering)
├── config.py                # Receptor/ligand definitions, Ensembl IDs, families, mappings
├── download_to_db.py        # Offline TCGA data downloader → SQLite
├── requirements.txt         # Python dependencies
├── README.md                # This file
└── pages/
    └── Literature_Overview.py  # Methodology docs, receptor profiles, references
```

---

## Methodology

### Normalization Pipeline

1. **TPM** (within-sample): corrects for gene length and sequencing depth
2. **Upper Quartile** (between-sample): scales each sample's 75th percentile to the cohort median, correcting for technical batch variation
3. **log₂(x + 1)** transform for correlation analysis

### Core Metric: log₂ Fold-Change

All charts report **log₂(tumor_TPM / GTEx_normal_median_TPM)** — the standard differential expression metric. Key values:

| log₂FC | Fold-change | Interpretation |
|--------|-------------|----------------|
| 0 | 1× | Same as normal tissue |
| 1 | 2× | Doubled |
| 3 | 8× | Strongly upregulated |
| 5 | 32× | Very strongly upregulated |
| 10 | 1024× | Massively upregulated |

### Ligand Activation Score

For each receptor, we compute a per-patient score:

1. Identify which TME ligands bind that receptor
2. Back-transform from log₂ to linear TPM
3. Sum all ligands in linear space (biologically correct — each molecule independently engages the receptor)
4. log₂ fold-change vs the corresponding GTEx normal tissue sum

### Co-activation Network

For each pair of receptors, Spearman rank correlation is computed between their ligand activation scores across patients:

- **Blue edges (+ρ):** Co-activation — TME upregulates both pathways together → supports combination therapy
- **Red edges (−ρ):** Inverse correlation — TME tends to use one OR the other → may define distinct patient subgroups
- **Dashed edges:** The pair also shares at least one ligand
- **Black dotted edges:** Identical ligand set (trivial ρ=1, e.g., BTLA & CD160 both bind only HVEM)
- **Dashed orange edges:** Shared ligand but no significant correlation at current threshold
- **Node size:** Total log₂FC of that receptor's ligands vs normal tissue

---

## App Features

### Explorer Page

- **Co-activation network graph** — interactive Plotly graph with hover details, edge types, family-colored nodes
- **📊 Pathway Activation** — ridgeline/violin density plot showing per-patient ligand activation distributions for each receptor
- **🎯 Ligand Breakdown** — stacked bar chart: log₂FC per ligand per receptor, showing which ligand dominates each pathway
- **🔴 TME Ligand Activity** — horizontal bar chart of all 23 ligands ranked by log₂FC, with Q1–Q3 error bars
- **🔥 Correlation Matrix** — heatmap of ligand-activation Spearman ρ between all receptor pairs
- **📋 Edge Table** — sortable dataframe of all network edges with ρ, p-values, ligand lists
- **👥 Demographics** — sample counts by race and AJCC stage
- **🧪 Receptor–Ligand Pairs** — per-pair analysis with receptor log₂FC, ligand log₂FC, suppressive score, and receptor-ligand correlation

### Sidebar Filters

- **TCGA Project** — select cancer type
- **Population (Race)** — filter by race/ethnicity
- **Receptor Families** — multiselect to include/exclude specific families
- **Stage** — segmented control for AJCC pathologic stage (I/II/III/IV)
- **Network parameters** — |ρ| threshold, p-value cutoff, max samples

### Literature Overview Page

- Full methodology documentation
- TCGA → GTEx tissue mapping table
- Receptor–ligand mapping table (all 16×23 pairs)
- Notes on non-classical interactions (HVEM bidirectional signaling, CD39/CD73/A2A adenosine pathway)
- T cell exhaustion background
- Detailed receptor profiles with mechanism descriptions, Ensembl links, and PubMed references

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| **log₂FC over raw TPM** | Compresses dynamic range — constitutive genes (collagen at 1000 TPM) and immune-specific genes (PD-L1 at 15 TPM) become visually comparable |
| **GTEx baselines over cohort median** | Measures tumor-specific upregulation vs healthy tissue, not just relative position within a tumor cohort |
| **Sum ligands in linear space** | Summing log values = multiplying, which is biologically wrong. Each ligand molecule independently engages the receptor |
| **Both positive and negative correlations** | Positive = co-activation (combination therapy). Negative = mutually exclusive strategies (patient subgroup stratification) |
| **Skip identical ligand pairs** | BTLA & CD160 both bind only HVEM → ρ=1 trivially, shown as dotted line instead |
| **Violin over bar charts** | Reveals bimodal distributions (cold vs hot tumors) that bars with error bars would hide |
| **Upper Quartile normalization** | TPM is within-sample only. UQ corrects between-sample technical variation, same approach as GDC's own FPKM-UQ |

---

## Limitations

- **Bulk RNA-seq only** — cannot distinguish which cell type expresses a gene. PD-L1 mRNA could come from tumor cells, macrophages, or DCs.
- **Transcript ≠ protein** — post-transcriptional regulation (ubiquitination, glycosylation, receptor internalization) can decouple mRNA from surface protein.
- **GTEx tissue matching is approximate** — e.g., TCGA-HNSC uses GTEx "Minor Salivary Gland" as the closest available match for head & neck tissue.
- **No spatial information** — co-expression across bulk samples doesn't prove co-localization within the tumor.
- **TCGA demographics skew white/American** — population-level findings should be validated in more diverse cohorts.

---

## References

1. Thorsson V et al. (2018). "The Immune Landscape of Cancer." *Immunity* 48(4):812-830. [PMID: 29628290](https://pubmed.ncbi.nlm.nih.gov/29628290)
2. Ayers M et al. (2017). "IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade." *J Clin Invest* 127(8):2930-2940. [PMID: 28650338](https://pubmed.ncbi.nlm.nih.gov/28650338)
3. GTEx Consortium (2020). "The GTEx Consortium atlas of genetic regulatory effects across human tissues." *Science* 369(6509):1318-1330. [PMID: 32913098](https://pubmed.ncbi.nlm.nih.gov/32913098)
4. Bullard JH et al. (2010). "Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments." *BMC Bioinformatics* 11:94. [PMID: 20167110](https://pubmed.ncbi.nlm.nih.gov/20167110)
5. Robinson MD & Oshlack A (2010). "A scaling normalization method for differential expression analysis of RNA-seq data." *Genome Biol* 11:R25. [PMID: 20196867](https://pubmed.ncbi.nlm.nih.gov/20196867)