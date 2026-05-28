"""
Literature Overview — T Cell Exhaustion Receptors
==================================================
Methodology, receptor-ligand mapping, and detailed receptor profiles.
"""

import streamlit as st
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from config import RECEPTORS, LIGANDS, FAMILY_COLORS, FAMILY_LIST, TCGA_TO_GTEX_TISSUE, TCGA_PROJECTS

st.set_page_config(
    page_title="Literature Overview — T Cell Exhaustion Receptors",
    page_icon="📚",
    layout="wide",
)

st.markdown("""
<style>
.main .block-container { padding-top: 2rem; max-width: 900px; }
.receptor-card {
    border-left: 4px solid; padding: 1rem 1.2rem; margin-bottom: 1.5rem;
    border-radius: 0 8px 8px 0;
    background: rgba(128,128,128,0.04);
}
.ref-item { font-size: 0.85rem; line-height: 1.6; margin-bottom: 0.3rem; }
.ref-item a { color: #667eea; text-decoration: none; font-weight: 500; }
.ref-item a:hover { text-decoration: underline; color: #5a6fd6; }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════
# HEADER
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("# 📚 Literature Overview")
st.markdown("### T Cell Exhaustion, TME Suppression & Methodology")

st.markdown(
    f"This page documents the methodology behind the Explorer and provides "
    f"detailed profiles for each of the **{len(RECEPTORS)} receptors** and "
    f"**{len(LIGANDS)} TME ligands** tracked across **{len(TCGA_PROJECTS)} TCGA cancer types**."
)

st.markdown("---")


# ═══════════════════════════════════════════════════════════════════════════
# METHODOLOGY (placed first, before receptor profiles)
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("## Methodology")

st.markdown("### What This Tool Measures")
st.markdown(
    "This tool answers the question: **for a given cancer type, which T cell "
    "inhibitory pathways is the tumor microenvironment (TME) most actively "
    "engaging, and which pathways does it co-activate together?**"
)
st.markdown(
    "It does this by measuring the expression of TME ligands (the molecules "
    "tumor/stromal/myeloid cells produce to shut down T cells) and comparing "
    "them to healthy tissue. Everything is expressed as **% above normal** — "
    "how much more of each suppressive ligand the tumor produces compared to "
    "the same organ in a healthy person."
)

st.markdown("### Data Sources")
st.markdown(
    f"**Tumor data:** The Cancer Genome Atlas (TCGA) via the NCI Genomic Data "
    f"Commons (GDC) API. **All {len(TCGA_PROJECTS)} TCGA cancer types** are supported, "
    f"including 10 rare cancers and 2 hematologic malignancies. STAR-Counts gene expression "
    f"quantification files provide Transcripts Per Million (TPM) values from the STAR aligner "
    f"against GRCh38/hg38 with GENCODE v36 annotations."
)
st.markdown(
    "**Normal tissue baseline:** The Genotype-Tissue Expression (GTEx) project "
    "provides RNA-seq expression from **54 non-diseased tissue sites** across "
    "~1000 healthy individuals. For each TCGA cancer type, we query the GTEx "
    "API for the median TPM of each gene in the corresponding normal tissue. "
    "The tissue mappings are:"
)

# Show TCGA → GTEx mapping table
mapping_rows = []
for proj, name in sorted(TCGA_TO_GTEX_TISSUE.items()):
    mapping_rows.append({
        "TCGA Project": proj,
        "Cancer Type": TCGA_PROJECTS.get(proj, ""),
        "GTEx Normal Tissue": name.replace("_", " "),
    })
mapping_df = pd.DataFrame(mapping_rows)
st.dataframe(mapping_df, use_container_width=True, hide_index=True, height=400)

st.markdown("#### GTEx Mapping Notes")
st.markdown(
    "Not all TCGA cancer types have an exact match in GTEx. Approximate mappings are used "
    "where necessary:"
)
st.markdown(
    "- **TCGA-READ** (Rectum) → Colon Transverse — rectum not available in GTEx; colon is anatomically closest.\n"
    "- **TCGA-MESO** (Mesothelioma) → Lung — pleural mesothelioma originates in the lung lining.\n"
    "- **TCGA-SARC** (Sarcoma) → Adipose Subcutaneous — soft tissue sarcomas arise from mesenchymal tissue; adipose is the closest available.\n"
    "- **TCGA-THYM** (Thymoma) → Blood Vessel — no thymus tissue in GTEx; vascular endothelium is a rough proxy.\n"
    "- **TCGA-UVM** (Uveal Melanoma) → Skin Sun Exposed — no eye tissue in GTEx; melanocyte lineage shared with skin.\n"
    "- **TCGA-CHOL** (Cholangiocarcinoma) → Liver — bile duct epithelium originates from hepatic tissue.\n"
    "- **TCGA-PCPG** (Pheochromocytoma) → Adrenal Gland — correct anatomical origin.\n"
    "- **TCGA-LAML** (Acute Myeloid Leukemia) → Whole Blood — hematologic cancer; myeloid origin.\n"
    "- **TCGA-DLBC** (Diffuse Large B-Cell Lymphoma) → EBV-transformed Lymphocytes — B-cell baseline "
    "matching the cell of origin. Shows what malignant transformation upregulates vs normal B cells."
)
st.markdown(
    "⚠ **Hematologic cancers** (TCGA-LAML, TCGA-DLBC) lack a traditional solid tumor microenvironment. "
    "The receptor-ligand co-activation framework was designed for solid tumors where TME cells "
    "express suppressive ligands to shut down infiltrating T cells. In blood cancers, the malignant "
    "cells themselves ARE immune cells — they constitutively express HLA-E, HVEM, CD48, B7-1/B7-2, "
    "PD-L1, and other molecules at high levels because that's normal immune cell biology, not active "
    "suppression."
)
st.markdown(
    "**No good GTEx baseline exists for these cancers.** GTEx does not include lymph node tissue. "
    "The closest options are all problematic:"
)
st.markdown(
    "- **Spleen** — mostly red pulp (blood filtration), not germinal centers where B cells live. "
    "Massively underestimates normal B cell gene expression.\n"
    "- **EBV-transformed lymphocytes** — B cell lines in culture, no stromal component. "
    "Structural genes (collagen, cadherins) show as massively upregulated because the baseline has zero stroma.\n"
    "- **Whole Blood** — contains myeloid precursors but not in the tissue context of bone marrow. "
    "Used for LAML as the least-bad option."
)
st.markdown(
    "Results for LAML and DLBC should be interpreted as **relative expression patterns within the "
    "cohort** rather than true upregulation above normal tissue. The co-activation network (which "
    "uses rank-based Spearman correlation) is still informative for identifying which pathways "
    "co-vary across patients, even if the absolute log₂FC values are inflated."
)

st.markdown("### Why Transcript Levels Are a Reasonable Approximation")
st.markdown(
    "This tool uses mRNA transcript abundance (TPM) as a proxy for protein-level "
    "activity. This is an approximation, not a direct measurement."
)
st.markdown(
    "**Transcription correlates with protein for most immune checkpoints.** "
    "Multiple studies have validated that for PD-L1, PD-L2, CTLA-4 ligands, "
    "galectin-9, CD155, and others, mRNA expression correlates with protein "
    "expression measured by IHC or flow cytometry at the bulk tissue level. "
    "The TCGA Pan-Cancer Immune Landscape study (Thorsson et al., 2018, "
    "*Immunity*) demonstrated that RNA-based immune signatures are prognostically "
    "informative and concordant with protein-level immune phenotyping."
)
st.markdown(
    "**TCGA provides the largest uniformly processed cohort available.** "
    "With 11,000+ tumors across 33 cancer types, all processed through the same "
    "STAR pipeline, TCGA offers unmatched statistical power. No proteomic dataset "
    "comes close for immune checkpoint molecules."
)
st.markdown(
    "**Immune checkpoint therapies were developed using transcriptomic biomarkers.** "
    "TMB, IFN-γ signatures, and T cell-inflamed gene expression profiles used in "
    "clinical trials (KEYNOTE, CheckMate) are all RNA-based. The FDA-approved "
    "companion diagnostics for pembrolizumab include gene expression scores."
)
st.markdown("**Known limitations:**")
st.markdown(
    "Post-transcriptional regulation can decouple mRNA from protein. PD-L1's protein "
    "half-life is regulated by ubiquitination (CMTM6/CMTM4), glycosylation affects "
    "antibody detection, and exosomal PD-L1 acts at a distance. CTLA-4 is constitutively "
    "internalized and recycled. Bulk RNA-seq cannot distinguish which cell type expresses "
    "a gene — PD-L1 mRNA could come from tumor cells, macrophages, or DCs."
)

st.markdown("### Normalization Pipeline")
st.markdown(
    "**Within-sample (TPM):** Corrects for gene length and sequencing depth within "
    "each sample. Makes genes comparable within a single patient."
)
st.markdown(
    "**Between-sample (Upper Quartile):** Corrects for technical variation between "
    "samples. Each sample is scaled so its 75th percentile expression equals the "
    "median 75th percentile across all samples. Same approach as GDC's FPKM-UQ."
)
st.markdown(
    "**Log₂(x + 1) transform:** Compresses dynamic range for correlation analysis."
)

st.markdown("### log₂ Fold-Change — The Core Metric")
st.markdown(
    "All charts and node sizes use **log₂ fold-change over normal tissue** "
    "(log₂FC) as the primary unit. This is the standard metric in differential "
    "expression analysis (DESeq2, limma, edgeR all report log₂FC). The formula is:"
)
st.markdown(
    "    *log₂FC = log₂( (mean_tumor_TPM + 0.1) / GTEx_normal_median_TPM )*"
)
st.markdown(
    "The +0.1 pseudocount prevents log(0) for unexpressed genes. Only positive "
    "values are shown (upregulation in tumor vs normal). Key interpretation:"
)
st.markdown(
    "| log₂FC | Fold-change | Meaning |\n"
    "|--------|-------------|----------|\n"
    "| 0 | 1× | Same as normal tissue |\n"
    "| 1 | 2× | Doubled in tumor |\n"
    "| 2 | 4× | Quadrupled |\n"
    "| 3 | 8× | 8-fold higher |\n"
    "| 5 | 32× | Strongly upregulated |\n"
    "| 10 | 1024× | Massively upregulated |"
)
st.markdown(
    "**Why log₂FC, not raw percentage?** Raw percentage above normal produces "
    "misleading visuals when comparing genes with very different baseline expression. "
    "A ligand going from 0.5 TPM (normal skin) to 500 TPM (melanoma) is a 100,000% "
    "increase that visually crushes everything else. On the log₂ scale, that's ~10 — "
    "a large but readable value alongside a ligand at log₂FC = 3 (8× increase). "
    "The log scale correctly compresses the dynamic range so that constitutive "
    "structural genes (collagen, HLA-E) and immune-specific genes (PD-L1, galectin-9) "
    "are visually comparable."
)
st.markdown(
    "**For ligands with near-zero normal expression** (GTEx baseline < 0.1 TPM), "
    "the baseline is floored at 0.1 to avoid division by zero. This means a tumor "
    "expressing 10 TPM of a normally-silent gene shows as log₂(10/0.1) ≈ 6.6."
)
st.markdown(
    "**GTEx baseline is required.** The tool will not run without GTEx data — there is no "
    "fallback to cohort median. This ensures all log₂FC values are measured against healthy "
    "tissue, not relative to the tumor cohort."
)

st.markdown("### Ligand Activation Score & Co-activation Network")
st.markdown(
    "For each inhibitory receptor, we compute a per-patient **Ligand Activation Score**:"
)
st.markdown(
    "**Step 1 — Identify ligands.** Each receptor has known TME binding partners "
    "(see Receptor–Ligand Mapping table below)."
)
st.markdown(
    "**Step 2 — Sum in linear space.** For each patient, back-transform ligand "
    "expression from log₂ to linear TPM, then sum all ligands for a receptor. "
    "This gives the total ligand burden available to engage that receptor. "
    "Summing must be in linear space (not log) because adding log values would "
    "multiply rather than add the underlying quantities."
)
st.markdown(
    "**Step 3 — Correlate across patients.** For each pair of receptors, Spearman "
    "rank correlation is computed between their ligand activation scores. Both "
    "positive and negative correlations are shown if they pass the |ρ| threshold."
)

st.markdown("### Network Edge Types")
st.markdown(
    "| Edge Style | Color | Meaning |\n"
    "|---|---|---|\n"
    "| **Solid** | Blue | Positive co-activation (+ρ): TME upregulates both pathways together |\n"
    "| **Solid** | Red | Inverse correlation (−ρ): when one pathway is active, the other tends to be low |\n"
    "| **Dashed** | Blue/Red | Same as above, but the pair also shares at least one ligand |\n"
    "| **Dotted** | Black | Identical ligand set (e.g., BTLA & CD160 both bind only HVEM) — ρ=1 trivially |\n"
    "| **Dashed** | Orange | Shared ligand but no significant |ρ| at current threshold |"
)
st.markdown(
    "Edge thickness is proportional to |ρ|. Node size reflects total log₂ fold-change "
    "of that receptor's ligands compared to normal tissue — larger nodes receive more "
    "suppressive pressure from the TME."
)
st.markdown(
    "**Positive edges (blue)** reveal co-suppressive strategies: the TME engages both "
    "pathways simultaneously, suggesting combination therapy may be needed. "
    "**Negative edges (red)** reveal mutually exclusive strategies: the TME tends to use "
    "one OR the other, suggesting either monotherapy might be sufficient for different "
    "patient subgroups."
)

st.markdown("### Ligand Activation Distribution (Ridgeline Chart)")
st.markdown(
    "The **Pathway Activation** tab shows a ridgeline density plot for each receptor. "
    "Each row contains per-ligand KDE density curves (colored) overlaid on a gray "
    "histogram of per-patient total activation. The x-axis is log₂ fold-change vs "
    "GTEx normal tissue."
)
st.markdown(
    "**Gray histogram bars** show the per-patient combined ligand activation — one value "
    "per patient computed as log₂((L1 + L2 + ...) / (L1_normal + L2_normal + ...)). "
    "Hover shows the exact patient count per bin."
)
st.markdown(
    "**Colored density curves** show each ligand's individual distribution. Where multiple "
    "curves overlap, those ligands co-vary across patients. Where one curve is shifted far "
    "right while another sits near zero, the pathway is driven by a single dominant ligand."
)
st.markdown(
    "**Why density instead of a bar chart?** A bar chart with error bars assumes the "
    "distribution is roughly symmetric around the mean. In reality, ligand activation "
    "is often bimodal: some tumors are 'cold' (low ligand expression, near normal) "
    "while others are 'hot' (high expression, strong suppression). The ridgeline "
    "reveals this structure."
)

st.markdown("### Ligand Breakdown (Compound Bar Chart)")
st.markdown(
    "The **Ligand Breakdown** tab shows a horizontal bar chart with compound y-axis labels "
    "(Receptor · Ligand). Each bar is one receptor-ligand pair, showing the mean log₂FC "
    "with Q1–Q3 error bars for patient-level variability."
)
st.markdown(
    "Bars are grouped by receptor (most activated at top) and sorted by log₂FC within each "
    "group. This format avoids the thin-bar problem of grouped bar charts where receptors "
    "with 1 ligand get the same width allocation as receptors with 3."
)

st.markdown("### Clinical Interpretation")
st.markdown(
    "A thick blue edge between PD-1 and TIGIT means the TME co-activates both pathways — "
    "supporting combination immunotherapy (anti-PD-1 + anti-TIGIT). A thick red edge "
    "between two receptors means they are inversely activated — the TME tends to use "
    "one or the other, which may define distinct patient subgroups. A large isolated "
    "node means that pathway is heavily activated but independent — "
    "monotherapy may suffice. Stage filtering reveals how the suppressive network "
    "evolves with disease progression."
)

st.markdown("### Sample Size Independence")
st.markdown(
    "TCGA cohorts vary significantly in size — TCGA-BRCA has ~1095 patients while "
    "TCGA-PAAD has ~178. A natural question is whether the metrics in this tool are "
    "biased by sample count. Here is why they are not:"
)
st.markdown(
    "**log₂ fold-change** is computed as log₂(mean_tumor_TPM / GTEx_normal_TPM). "
    "The mean is already a per-patient average — it produces one value regardless "
    "of whether 178 or 1095 patients contributed. Larger cohorts give a more precise "
    "estimate of the true mean, but the value itself does not inflate with n."
)
st.markdown(
    "**Spearman correlation (ρ)** is a rank-based statistic. The ρ value measures "
    "the strength of monotonic association between two variables, independent of "
    "sample size. What does depend on n is statistical significance: with 1095 patients, "
    "even a tiny ρ = 0.08 can achieve p < 0.05, while with 178 patients you need "
    "ρ ≈ 0.15. This is why we filter on **effect size (|ρ| threshold)** rather than "
    "relying solely on p-values. The default |ρ| ≥ 0.5 ensures that only biologically "
    "meaningful correlations appear as edges, regardless of cohort size."
)
st.markdown(
    "**Q1/Q3 error bars and violin shapes** are quantiles and kernel density estimates "
    "of the per-patient distribution. More patients produce smoother, more precise "
    "estimates of the true distribution shape, but the quartile values themselves "
    "do not scale with n."
)
st.markdown(
    "**Ligand activation scores** are computed per-patient (sum of ligand TPMs for "
    "each individual), then correlated across patients. The per-patient score is "
    "independent of how many other patients exist in the cohort."
)
st.markdown(
    "**Population comparison chart** is the one area where sample size matters most. "
    "If White patients have n=800 and Black patients have n=50, the White mean is much "
    "more precisely estimated. We show sample counts in the axis labels (e.g., "
    "'White (n=800)') so readers can assess the reliability of each group's estimate. "
    "These are exploratory observations — formal comparisons would require matched "
    "cohorts and statistical testing with appropriate correction for confounders."
)
st.markdown(
    "**The ridgeline histogram bars** show raw frequency counts within each receptor's "
    "row. Since every patient has every gene measured, all receptors within a project "
    "share the same n — the bar heights reflect the shape of the distribution, not "
    "differences in sample availability. No normalization is needed."
)

st.markdown("### Methodology References")
method_refs = [
    ("Thorsson V et al. (2018). \"The Immune Landscape of Cancer.\" <i>Immunity</i> 48(4):812-830.", "29628290"),
    ("Ayers M et al. (2017). \"IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade.\" <i>J Clin Invest</i> 127(8):2930-2940.", "28650338"),
    ("GTEx Consortium (2020). \"The GTEx Consortium atlas of genetic regulatory effects across human tissues.\" <i>Science</i> 369(6509):1318-1330.", "32913098"),
    ("Bullard JH et al. (2010). \"Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments.\" <i>BMC Bioinformatics</i> 11:94.", "20167110"),
    ("Robinson MD & Oshlack A (2010). \"A scaling normalization method for differential expression analysis of RNA-seq data.\" <i>Genome Biol</i> 11:R25.", "20196867"),
]
for i, (text, pmid) in enumerate(method_refs, 1):
    st.markdown(
        f'<div class="ref-item">{i}. {text} '
        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}" target="_blank">PubMed: {pmid}</a></div>',
        unsafe_allow_html=True,
    )

st.markdown("---")


# ═══════════════════════════════════════════════════════════════════════════
# RECEPTOR–LIGAND MAPPING TABLE
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("## Receptor–Ligand Mapping")
st.markdown(
    "Each T cell inhibitory receptor is engaged by one or more ligands expressed "
    "in the tumor microenvironment. This table shows every receptor–ligand pair "
    "tracked in the Explorer."
)

st.markdown("#### ⚠ Non-classical interactions")
st.markdown(
    "Most pairs in this tool follow the classical model: a receptor on the T cell "
    "binds a ligand secreted or surface-expressed by tumor/stromal/myeloid cells in "
    "the TME. However, a few entries are non-classical and should be interpreted "
    "with that context:"
)
st.markdown(
    "**HVEM (TNFRSF14) → BTLA / CD160:** HVEM is a bidirectional signaling molecule. "
    "It is expressed on T cells themselves as well as on tumor cells and APCs. When "
    "HVEM on a tumor cell engages BTLA or CD160 on a T cell, it delivers an inhibitory "
    "signal — this is the interaction we measure. However, HVEM on the T cell can also "
    "deliver co-stimulatory signals in the opposite direction (via LIGHT/LTα), so its "
    "role is context-dependent. Bulk RNA-seq cannot distinguish which cell type expresses "
    "HVEM in a given sample."
)
st.markdown(
    "**CEACAM1 → TIM-3:** CEACAM1 functions as a cis-ligand for TIM-3, meaning both "
    "molecules are co-expressed on the same T cell surface and interact in cis rather "
    "than in trans (across cells). This cis-interaction is required for TIM-3's inhibitory "
    "signaling. We track CEACAM1 as a TIM-3 'ligand' because its expression in the bulk "
    "tumor sample reflects the availability of CEACAM1 to stabilize TIM-3 signaling, "
    "but the interaction is mechanistically different from a TME-derived ligand like PD-L1."
)
st.markdown(
    "**A2A receptor (ADORA2A) / CD39 (ENTPD1) / CD73 (NT5E) — the adenosine pathway:** "
    "A2A is the actual inhibitory receptor on T cells — a GPCR that binds extracellular "
    "adenosine and suppresses T cell function via cAMP elevation. CD39 and CD73 are the "
    "enzymes that generate adenosine: CD39 (on exhausted T cells and Tregs) converts "
    "ATP → AMP, and CD73 (on tumor and stromal cells) converts AMP → adenosine. We model "
    "A2A as the receptor and CD39+CD73 as its 'ligands' because their expression indicates "
    "adenosine is being generated to engage A2A. This is a three-step metabolic cascade "
    "rather than a direct receptor-ligand binding event — CD39 and CD73 don't bind A2A "
    "directly, but their co-expression is the best RNA-seq proxy for adenosine pathway "
    "activation."
)
st.markdown(
    "**LSECtin (CLEC4G) → LAG-3:** LSECtin was originally characterized in liver sinusoidal "
    "endothelial cells but is also expressed on melanoma and other tumor cells. It binds "
    "LAG-3 and inhibits IFN-γ production from T cells. Its expression is tissue-dependent "
    "and may be very low outside liver-related cancers."
)

rl_rows = []
for receptor_gene, rinfo in RECEPTORS.items():
    # Find ligands for this receptor
    receptor_ligands = [
        (lg, linfo) for lg, linfo in LIGANDS.items()
        if receptor_gene in linfo.get("receptors", [])
    ]
    if receptor_ligands:
        for lg, linfo in receptor_ligands:
            rl_rows.append({
                "Receptor": rinfo["label"],
                "Receptor Gene": receptor_gene,
                "Receptor Ensembl": rinfo["ensembl"],
                "Ligand": linfo["label"],
                "Ligand Gene": lg,
                "Ligand Ensembl": linfo["ensembl"],
                "Family": rinfo["family"],
            })
    else:
        rl_rows.append({
            "Receptor": rinfo["label"],
            "Receptor Gene": receptor_gene,
            "Receptor Ensembl": rinfo["ensembl"],
            "Ligand": "—",
            "Ligand Gene": "—",
            "Ligand Ensembl": "—",
            "Family": rinfo["family"],
        })

# Also add general TME suppressive molecules not tied to a specific receptor
for lg, linfo in LIGANDS.items():
    if not linfo["receptors"]:
        rl_rows.append({
            "Receptor": "TME (general)",
            "Receptor Gene": "—",
            "Receptor Ensembl": "—",
            "Ligand": linfo["label"],
            "Ligand Gene": lg,
            "Ligand Ensembl": linfo["ensembl"],
            "Family": "General TME Suppression",
        })

rl_df = pd.DataFrame(rl_rows)
st.dataframe(rl_df, use_container_width=True, height=500, hide_index=True)

st.caption(
    f"{len(RECEPTORS)} receptors × {len(LIGANDS)} ligands tracked. "
    "Receptors are on the T cell surface. Ligands are expressed by "
    "tumor cells, stromal cells, or myeloid cells in the TME."
)

st.markdown("---")


# ═══════════════════════════════════════════════════════════════════════════
# BACKGROUND ON T CELL EXHAUSTION
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("## What is T Cell Exhaustion?")
st.markdown(
    "T cell exhaustion is a state of progressive dysfunction that develops when T cells "
    "are exposed to persistent antigen stimulation — as occurs in chronic viral infections "
    "and in the tumor microenvironment. Exhausted T cells lose their ability to kill target "
    "cells and produce effector cytokines (IFN-γ, TNF-α, IL-2) in a hierarchical manner: "
    "IL-2 production is lost first, followed by TNF-α, and finally IFN-γ and cytotoxic "
    "capacity."
)
st.markdown(
    "The hallmark of exhaustion is the sustained, high-level co-expression of multiple "
    "inhibitory receptors on the T cell surface. No single receptor defines exhaustion — "
    "rather, it's the combinatorial co-expression pattern (e.g., PD-1⁺ TIM-3⁺ LAG-3⁺ TIGIT⁺) "
    "that marks the dysfunctional state. This is why co-activation network analysis is "
    "informative: it reveals which suppressive pathways the TME engages simultaneously."
)
st.markdown(
    "Anti-PD-1 therapy works primarily by reinvigorating the \"progenitor exhausted\" (Tpex) "
    "subset rather than fully reversing terminal exhaustion. Understanding which additional "
    "pathways the TME co-activates alongside PD-1 is critical for designing effective "
    "combination immunotherapy."
)

st.markdown("---")


# ═══════════════════════════════════════════════════════════════════════════
# RECEPTOR OVERVIEW TABLE
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("## Receptor Reference")

families = {}
for gene, info in RECEPTORS.items():
    families.setdefault(info["family"], []).append((gene, info))

overview_rows = []
for gene, info in RECEPTORS.items():
    overview_rows.append({
        "Gene Symbol": gene,
        "Common Name": info["label"],
        "Ensembl ID": info["ensembl"],
        "Chromosome": info.get("chromosome", "—"),
        "Family": info["family"],
    })

overview_df = pd.DataFrame(overview_rows)
st.dataframe(overview_df, use_container_width=True, height=400, hide_index=True)

st.markdown("---")


# ═══════════════════════════════════════════════════════════════════════════
# DETAILED RECEPTOR PROFILES
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("## Detailed Receptor Profiles")

for fam_name in FAMILY_LIST:
    if fam_name not in families:
        continue

    color = FAMILY_COLORS.get(fam_name, "#888888")
    st.markdown(f"### {fam_name}")

    for gene, info in families[fam_name]:
        ensembl_url = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={info['ensembl']}"

        st.markdown(
            f'<div class="receptor-card" style="border-left-color: {color};">',
            unsafe_allow_html=True,
        )

        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown(f"#### {info['label']}  (`{gene}`)")
        with col2:
            st.markdown(f"[Ensembl ↗]({ensembl_url})")

        # Gene info
        st.markdown(
            f"**Ensembl:** `{info['ensembl']}` · "
            f"**Chromosome:** {info.get('chromosome', '—')} · "
            f"**Family:** {info['family']}"
        )

        # Ligands for this receptor
        receptor_ligs = [
            LIGANDS[lg]["label"] for lg in LIGANDS
            if gene in LIGANDS[lg].get("receptors", [])
        ]
        if receptor_ligs:
            st.markdown(f"**TME Ligands:** {', '.join(receptor_ligs)}")

        # Mechanism
        st.markdown("**How it shuts down the T cell:**")
        st.markdown(info.get("mechanism", "_No mechanism description available._"))

        # References
        refs = info.get("references", [])
        if refs:
            st.markdown("**Key references:**")
            for i, ref in enumerate(refs, 1):
                pmid = None
                if "PMID:" in ref:
                    pmid_part = ref.split("PMID:")[-1].strip()
                    pmid = "".join(c for c in pmid_part if c.isdigit())

                if pmid:
                    cite_text = ref.split("PMID:")[0].strip().rstrip(".")
                    st.markdown(
                        f'<div class="ref-item">{i}. {cite_text}. '
                        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}" target="_blank">'
                        f'PubMed: {pmid}</a></div>',
                        unsafe_allow_html=True,
                    )
                else:
                    st.markdown(f'<div class="ref-item">{i}. {ref}</div>', unsafe_allow_html=True)

        st.markdown("</div>", unsafe_allow_html=True)
        st.markdown("")

    st.markdown("---")

st.markdown(
    "**Population differences:** Variability across racial/ethnic groups may reflect "
    "differences in tumor immunogenicity, germline immune gene polymorphisms, or HLA "
    "diversity. These are exploratory observations from TCGA bulk RNA-seq and should be "
    "validated with single-cell data and controlled cohort studies."
)

st.markdown("---")
st.caption(
    f"Literature Overview · {len(RECEPTORS)} receptors · {len(LIGANDS)} ligands · "
    f"{len(TCGA_PROJECTS)} cancer types · References link to PubMed"
)