"""
Shared configuration, receptor definitions, and constants.
"""

import plotly.express as px

# ══════════════════════════════════════════════════════════════════════════════
# T CELL EXHAUSTION / INHIBITORY RECEPTORS
# Only receptors expressed ON the T cell surface (+ key exhaustion TFs)
# ══════════════════════════════════════════════════════════════════════════════

RECEPTORS = {
    # ── Classical checkpoint receptors ────────────────────────────────────
    "PDCD1": {
        "label": "PD-1",
        "family": "Classical Checkpoints",
        "ensembl": "ENSG00000188389",
        "chromosome": "2q37.3",
        "desc": "Primary exhaustion marker; ITIM/ITSM-bearing receptor",
        "mechanism": (
            "PD-1 (Programmed Cell Death Protein 1) is the canonical T cell exhaustion receptor. "
            "When engaged by its ligands PD-L1 or PD-L2 on tumor cells or APCs, PD-1's cytoplasmic "
            "ITIM and ITSM motifs recruit phosphatases SHP-1 and SHP-2, which dephosphorylate key "
            "signaling molecules downstream of the TCR and CD28. This dampens PI3K/AKT and "
            "RAS/MAPK/ERK signaling, reducing IL-2 production, proliferation, and cytotoxic "
            "effector function. Sustained PD-1 expression under chronic antigen stimulation "
            "drives the epigenetic exhaustion program via TOX."
        ),
        "references": [
            "Ishida Y et al. (1992) EMBO J 11(11):3887-95. PMID: 1396582",
            "Sharpe AH & Pauken KE (2018) Nat Rev Immunol 18:153-167. PMID: 28990585",
            "Patsoukis N et al. (2020) Sci Adv 6(38):eabd7929. PMID: 32948597",
        ],
    },
    "CTLA4": {
        "label": "CTLA-4",
        "family": "Classical Checkpoints",
        "ensembl": "ENSG00000163599",
        "chromosome": "2q33.2",
        "desc": "CD28 family; competes for B7 ligands",
        "mechanism": (
            "CTLA-4 (Cytotoxic T-Lymphocyte Associated Protein 4) outcompetes the co-stimulatory "
            "receptor CD28 for binding to B7-1 (CD80) and B7-2 (CD86) on antigen-presenting cells, "
            "due to its ~20-fold higher affinity. This removes essential co-stimulatory signals "
            "needed for full T cell activation. Additionally, CTLA-4 recruits phosphatase PP2A "
            "to its cytoplasmic tail, directly inhibiting AKT signaling. CTLA-4 also physically "
            "removes B7 ligands from APC surfaces via trans-endocytosis, reducing co-stimulation "
            "for other T cells in the microenvironment."
        ),
        "references": [
            "Brunet JF et al. (1987) Nature 328:267-270. PMID: 3496540",
            "Leach DR et al. (1996) Science 271:1734-1736. PMID: 8596936",
            "Qureshi OS et al. (2011) Science 332:600-603. PMID: 21474713",
        ],
    },
    "LAG3": {
        "label": "LAG-3",
        "family": "Classical Checkpoints",
        "ensembl": "ENSG00000089692",
        "chromosome": "12p13.31",
        "desc": "MHC-II binding; synergizes with PD-1",
        "mechanism": (
            "LAG-3 (Lymphocyte Activation Gene 3) is structurally homologous to CD4 and binds "
            "MHC class II with higher affinity than CD4, competitively blocking CD4-MHC-II "
            "interaction. Its intracellular KIEELE motif is required for inhibitory function. "
            "LAG-3 also binds FGL1 (fibrinogen-like protein 1) and galectin-3 in the tumor "
            "microenvironment. It synergizes strongly with PD-1: dual LAG-3/PD-1 blockade shows "
            "superior anti-tumor responses compared to single-agent therapy. LAG-3 is a hallmark "
            "of the terminally exhausted T cell state."
        ),
        "references": [
            "Triebel F et al. (1990) J Exp Med 171:1393-1405. PMID: 1692078",
            "Woo SR et al. (2012) Cancer Res 72:917-927. PMID: 22186141",
            "Wang J et al. (2019) Cell 176:334-349. PMID: 30580966",
        ],
    },
    "HAVCR2": {
        "label": "TIM-3",
        "family": "Classical Checkpoints",
        "ensembl": "ENSG00000135077",
        "chromosome": "5q33.3",
        "desc": "Galectin-9/CEACAM1 receptor; marks terminal exhaustion",
        "mechanism": (
            "TIM-3 (T cell Immunoglobulin and Mucin domain-containing protein 3, encoded by "
            "HAVCR2) binds multiple ligands: galectin-9, phosphatidylserine, CEACAM1, and HMGB1. "
            "Galectin-9 binding triggers calcium influx-dependent cell death in Th1 cells. "
            "CEACAM1 forms a heterodimer with TIM-3 in cis, which is required for TIM-3's "
            "inhibitory function — without CEACAM1, TIM-3 cannot signal. TIM-3 marks the most "
            "terminally exhausted CD8+ T cells in tumors and co-expression with PD-1 defines "
            "the dysfunctional population least responsive to anti-PD-1 therapy."
        ),
        "references": [
            "Monney L et al. (2002) Nature 415:536-541. PMID: 11823861",
            "Zhu C et al. (2005) Nat Immunol 6:1245-1252. PMID: 16286920",
            "Huang YH et al. (2015) Nature 517:386-390. PMID: 25363763",
        ],
    },
    "TIGIT": {
        "label": "TIGIT",
        "family": "Classical Checkpoints",
        "ensembl": "ENSG00000181847",
        "chromosome": "3q13.31",
        "desc": "ITIM-bearing; CD155/CD112 receptor",
        "mechanism": (
            "TIGIT (T cell Immunoreceptor with Ig and ITIM domains) competes with the "
            "co-stimulatory receptor CD226 (DNAM-1) for binding to CD155 (PVR) and CD112 "
            "(nectin-2) on tumor cells and APCs. TIGIT binds CD155 with much higher affinity "
            "than CD226, displacing the activating signal. Its cytoplasmic ITIM recruits SHIP-1 "
            "phosphatase, inhibiting PI3K and MAPK signaling. TIGIT also directly inhibits "
            "NK cell cytotoxicity and promotes the generation of tolerogenic dendritic cells. "
            "TIGIT+ PD-1+ CD8 T cells represent a severely exhausted population in tumors."
        ),
        "references": [
            "Yu X et al. (2009) Nat Immunol 10:48-57. PMID: 19011627",
            "Johnston RJ et al. (2014) Cancer Cell 26:923-937. PMID: 25465800",
            "Chauvin JM et al. (2015) J Clin Invest 125:2046-2058. PMID: 25866972",
        ],
    },

    # ── Co-inhibitory Ig superfamily ─────────────────────────────────────
    "BTLA": {
        "label": "BTLA",
        "family": "Co-inhibitory Ig Superfamily",
        "ensembl": "ENSG00000186265",
        "chromosome": "3q13.2",
        "desc": "HVEM receptor; ITIM/ITSM signaling",
        "mechanism": (
            "BTLA (B and T Lymphocyte Attenuator) binds HVEM (TNFRSF14), a TNF receptor "
            "superfamily member — an unusual interaction between Ig superfamily and TNFR "
            "superfamily members. BTLA's ITIM and ITSM motifs recruit SHP-1 and SHP-2, "
            "attenuating TCR signaling. Unlike most exhaustion receptors which are upregulated, "
            "BTLA is constitutively expressed on naive T cells and downregulated upon "
            "differentiation to effector cells, but re-expressed during exhaustion. BTLA "
            "particularly restrains CD8+ T cell persistence in the tumor microenvironment."
        ),
        "references": [
            "Watanabe N et al. (2003) Nat Immunol 4:670-679. PMID: 12796776",
            "Sedy JR et al. (2005) Nat Immunol 6:90-98. PMID: 15568026",
            "Derré L et al. (2010) J Clin Invest 120:157-167. PMID: 20038811",
        ],
    },
    "CD160": {
        "label": "CD160",
        "family": "Co-inhibitory Ig Superfamily",
        "ensembl": "ENSG00000117281",
        "chromosome": "1q21.1",
        "desc": "GPI-anchored; HVEM receptor; marks severe exhaustion",
        "mechanism": (
            "CD160 is a GPI-anchored glycoprotein that, like BTLA, binds HVEM. CD160 also "
            "interacts with classical and non-classical MHC-I molecules. On CD8+ T cells, "
            "CD160 engagement delivers an inhibitory signal that reduces IFN-γ and TNF-α "
            "production. CD160 co-expression with PD-1 and 2B4 defines the most severely "
            "exhausted T cell subset in chronic viral infections (HCV, HIV) and in the "
            "tumor microenvironment. Blocking CD160 partially restores T cell function."
        ),
        "references": [
            "Cai G et al. (2008) Nat Immunol 9:176-185. PMID: 18193050",
            "Bengsch B et al. (2010) PLoS Pathog 6:e1000947. PMID: 20548953",
            "Viganò S et al. (2014) Blood 124:2657-2665. PMID: 25212332",
        ],
    },
    "CD244": {
        "label": "2B4 (CD244)",
        "family": "SLAM Family",
        "ensembl": "ENSG00000122223",
        "chromosome": "1q23.3",
        "desc": "SLAM family; CD48 receptor; inhibitory without SAP",
        "mechanism": (
            "2B4 (CD244/SLAMF4) binds CD48, a ubiquitously expressed SLAM family member. "
            "In healthy T cells with functional SAP adaptor, 2B4 delivers co-stimulatory signals. "
            "However, in exhausted T cells where SAP is downregulated, 2B4 switches to an "
            "inhibitory receptor by recruiting EAT-2 and phosphatases instead of SAP-associated "
            "signaling molecules. This dual function makes 2B4 context-dependent: activating "
            "in acute immune responses, inhibitory in chronic exhaustion. Co-expression of 2B4 "
            "with PD-1, CD160, and KLRG1 marks severely exhausted virus-specific CD8+ T cells."
        ),
        "references": [
            "Brown MH et al. (1998) J Exp Med 188:2083-2090. PMID: 9841922",
            "Wherry EJ et al. (2007) Immunity 27:670-684. PMID: 17950003",
            "Bengsch B et al. (2010) PLoS Pathog 6:e1000947. PMID: 20548953",
        ],
    },
    "VSIR": {
        "label": "VISTA",
        "family": "Co-inhibitory Ig Superfamily",
        "ensembl": "ENSG00000107738",
        "chromosome": "10q22.1",
        "desc": "B7 family; PSGL-1 receptor; suppresses early T cell activation",
        "mechanism": (
            "VISTA (V-domain Ig Suppressor of T cell Activation, encoded by VSIR/C10orf54) "
            "functions both as a ligand on APCs and as a receptor on T cells. As a receptor, "
            "VISTA engagement by VSIG-3 at acidic pH (common in the tumor microenvironment) "
            "suppresses T cell activation at early stages. VISTA inhibits ERK and AKT "
            "phosphorylation and promotes the expression of Foxp3 (Treg differentiation). "
            "VISTA is unique in that its inhibitory function is pH-dependent, making it "
            "particularly active in the acidic tumor microenvironment."
        ),
        "references": [
            "Wang L et al. (2011) J Exp Med 208:577-592. PMID: 21383057",
            "Johnston RJ et al. (2019) Nature 574:565-570. PMID: 31645726",
            "Yuan L et al. (2021) Front Immunol 12:660714. PMID: 33995384",
        ],
    },
    "CD200R1": {
        "label": "CD200R",
        "family": "Co-inhibitory Ig Superfamily",
        "ensembl": "ENSG00000163606",
        "chromosome": "3q13.2",
        "desc": "CD200 receptor; delivers inhibitory signal via DOK2",
        "mechanism": (
            "CD200R1 binds CD200, which is widely overexpressed on tumor cells. Upon engagement, "
            "CD200R1's cytoplasmic NPxY motif recruits the adaptor protein DOK2, which then "
            "recruits RasGAP and SHIP, inhibiting RAS/MAPK and PI3K signaling respectively. "
            "This suppresses T cell proliferation and cytokine production. The CD200/CD200R "
            "axis is particularly important in glioblastoma, ovarian cancer, and CLL, where "
            "CD200 overexpression creates an immunosuppressive barrier."
        ),
        "references": [
            "Wright GJ et al. (2003) Immunity 18:391-402. PMID: 12648455",
            "Rygiel TP et al. (2012) Eur J Immunol 42:1174-1184. PMID: 22539293",
            "Coles SJ et al. (2012) Leukemia 26:2146-2148. PMID: 22425898",
        ],
    },

    # ── NK receptor family on T cells ────────────────────────────────────
    "KLRG1": {
        "label": "KLRG1",
        "family": "NK-related Receptors",
        "ensembl": "ENSG00000139187",
        "chromosome": "12p13.31",
        "desc": "E-cadherin receptor; marks senescent/exhausted T cells",
        "mechanism": (
            "KLRG1 (Killer cell Lectin-like Receptor G1) binds E-cadherin, N-cadherin, and "
            "R-cadherin on epithelial and tumor cells. Its cytoplasmic ITIM recruits SHIP-1 "
            "and SHP-2, inhibiting PI3K/AKT and TCR signaling. KLRG1 marks short-lived "
            "effector CD8 T cells and is associated with both terminal differentiation and "
            "exhaustion. KLRG1+ T cells have reduced proliferative capacity and shortened "
            "telomeres. Co-expression with PD-1 and 2B4 identifies an exhausted population "
            "distinct from the KLRG1- memory precursor pool."
        ),
        "references": [
            "Voehringer D et al. (2002) Eur J Immunol 32:3049-3056. PMID: 12385025",
            "Joshi NS et al. (2007) Immunity 27:281-295. PMID: 17723218",
            "Li L et al. (2018) J Clin Invest 128:1329-1340. PMID: 29324452",
        ],
    },
    "KLRC1": {
        "label": "NKG2A",
        "family": "NK-related Receptors",
        "ensembl": "ENSG00000134545",
        "chromosome": "12p13.2",
        "desc": "HLA-E receptor; tumor-specific exhaustion marker",
        "mechanism": (
            "NKG2A (encoded by KLRC1) forms a heterodimer with CD94 and binds HLA-E, a "
            "non-classical MHC-I molecule that presents leader peptides from other HLA molecules. "
            "NKG2A's cytoplasmic ITIM recruits SHP-1, directly inhibiting cytotoxic granule "
            "release and cytokine production. NKG2A is specifically upregulated on tumor-"
            "infiltrating CD8+ T cells but not on virus-specific exhausted T cells, making it "
            "a tumor-specific exhaustion marker. Monalizumab (anti-NKG2A) is in clinical trials "
            "for combination checkpoint blockade."
        ),
        "references": [
            "Braud VM et al. (1998) Nature 391:795-799. PMID: 9486650",
            "André P et al. (2018) Cell 175:1731-1743. PMID: 30503213",
            "Chen X et al. (2022) Mol Oncol 2(1):e111. PMID: 35036984",
        ],
    },
    "CD96": {
        "label": "CD96 (TACTILE)",
        "family": "NK-related Receptors",
        "ensembl": "ENSG00000153283",
        "chromosome": "3q13.13",
        "desc": "CD155 receptor; competes with TIGIT",
        "mechanism": (
            "CD96 (TACTILE) belongs to the same receptor family as TIGIT and CD226, all "
            "competing for CD155 (PVR) binding. CD96 has intermediate affinity for CD155 "
            "(between TIGIT and CD226) and delivers an inhibitory signal via its cytoplasmic "
            "ITIM-like motif. CD96 blockade enhances NK cell and T cell anti-tumor responses "
            "in preclinical models. In the exhaustion context, CD96 co-expression with TIGIT "
            "creates redundant inhibition through the CD155 axis, which may explain resistance "
            "to single-agent TIGIT blockade."
        ),
        "references": [
            "Fuchs A et al. (2004) J Immunol 172:3994-3998. PMID: 15034010",
            "Chan CJ et al. (2014) Nat Immunol 15:431-438. PMID: 24658051",
            "Blake SJ et al. (2016) Clin Cancer Res 22:5183-5190. PMID: 27267849",
        ],
    },

    # ── Adenosine pathway receptor on T cells ──────────────────────────────
    "ADORA2A": {
        "label": "A2A Receptor",
        "family": "Metabolic Exhaustion Markers",
        "ensembl": "ENSG00000128271",
        "chromosome": "22q11.23",
        "desc": "Adenosine receptor; GPCR; suppresses T cell function via cAMP",
        "mechanism": (
            "The A2A receptor (ADORA2A) is a G-protein coupled receptor on T cells that binds "
            "extracellular adenosine — the end product of the CD39/CD73 ectonucleotidase pathway. "
            "In the tumor microenvironment, CD39 (on exhausted T cells and Tregs) hydrolyzes ATP "
            "to AMP, and CD73 (on tumor and stromal cells) converts AMP to adenosine. When "
            "adenosine engages A2A, it activates adenylyl cyclase, elevating intracellular cAMP "
            "and suppressing TCR signaling, cytokine production (IFN-γ, TNF-α, IL-2), and "
            "cytotoxic function. A2A signaling also promotes Treg differentiation and inhibits "
            "NK cell cytotoxicity. A2A receptor blockade is under clinical investigation as "
            "an immunotherapy strategy (e.g., ciforadenant/CPI-444)."
        ),
        "references": [
            "Ohta A et al. (2006) Proc Natl Acad Sci USA 103:13132-13137. PMID: 16916931",
            "Vigano S et al. (2019) Front Immunol 10:925. PMID: 31114584",
            "Fong L et al. (2020) Cancer Discov 10:40-51. PMID: 31732494",
        ],
    },

    # ── Additional inhibitory receptors ──────────────────────────────────
    "LAIR1": {
        "label": "LAIR-1",
        "family": "Additional Inhibitory Receptors",
        "ensembl": "ENSG00000167613",
        "chromosome": "19q13.42",
        "desc": "Collagen receptor; ITIM-bearing; inhibits T cell activation",
        "mechanism": (
            "LAIR-1 (Leukocyte Associated Immunoglobulin-like Receptor 1) binds collagen and "
            "collagen-domain containing proteins (including C1q complement) in the extracellular "
            "matrix. The tumor stroma is collagen-rich, providing abundant LAIR-1 ligands. "
            "Upon collagen binding, LAIR-1's two cytoplasmic ITIMs recruit SHP-1 and SHP-2, "
            "inhibiting TCR signaling and cytokine production. LAIR-1 is upregulated on "
            "T cells in the tumor microenvironment and contributes to collagen-mediated "
            "immune suppression."
        ),
        "references": [
            "Meyaard L et al. (1997) Immunity 7:283-290. PMID: 9285412",
            "Lebbink RJ et al. (2006) J Exp Med 203:1419-1425. PMID: 16754720",
            "Peng DH et al. (2020) Nat Commun 11:4520. PMID: 32908152",
        ],
    },
    "LILRB1": {
        "label": "ILT2 (LILRB1)",
        "family": "Additional Inhibitory Receptors",
        "ensembl": "ENSG00000104972",
        "chromosome": "19q13.42",
        "desc": "MHC-I receptor; inhibits CD8 T cell cytotoxicity",
        "mechanism": (
            "LILRB1 (ILT2/CD85j/LIR-1) binds a broad range of classical and non-classical "
            "MHC-I molecules, including HLA-G which is overexpressed by many tumors. LILRB1's "
            "four cytoplasmic ITIMs recruit SHP-1, directly inhibiting TCR proximal signaling "
            "and cytotoxic granule polarization. On CD8+ T cells, LILRB1 inhibits target cell "
            "killing even when the TCR is engaged. HLA-G-expressing tumors exploit LILRB1 "
            "to evade CD8+ T cell-mediated killing. LILRB1 expression increases on TILs "
            "compared to peripheral T cells."
        ),
        "references": [
            "Colonna M et al. (1997) J Exp Med 186:1809-1818. PMID: 9382880",
            "Shiroishi M et al. (2003) Proc Natl Acad Sci 100:8856-8861. PMID: 12853576",
            "Barkal AA et al. (2018) Nat Immunol 19:76-84. PMID: 29180808",
        ],
    },
}

GENE_ENSEMBL = {gene: info["ensembl"] for gene, info in RECEPTORS.items()}

TCGA_PROJECTS = {
    "TCGA-BRCA": "Breast Invasive Carcinoma",
    "TCGA-LUAD": "Lung Adenocarcinoma",
    "TCGA-LUSC": "Lung Squamous Cell Carcinoma",
    "TCGA-SKCM": "Skin Cutaneous Melanoma",
    "TCGA-COAD": "Colon Adenocarcinoma",
    "TCGA-BLCA": "Bladder Urothelial Carcinoma",
    "TCGA-HNSC": "Head & Neck Squamous Cell Carcinoma",
    "TCGA-KIRC": "Kidney Renal Clear Cell Carcinoma",
    "TCGA-LIHC": "Liver Hepatocellular Carcinoma",
    "TCGA-PRAD": "Prostate Adenocarcinoma",
    "TCGA-GBM":  "Glioblastoma Multiforme",
    "TCGA-OV":   "Ovarian Serous Cystadenocarcinoma",
    "TCGA-UCEC": "Uterine Corpus Endometrial Carcinoma",
    "TCGA-PAAD": "Pancreatic Adenocarcinoma",
    "TCGA-STAD": "Stomach Adenocarcinoma",
}

# TCGA project → corresponding GTEx normal tissue (for baseline)
# GTEx tissue IDs from https://gtexportal.org/api/v2/dataset/tissueSiteDetail
TCGA_TO_GTEX_TISSUE = {
    "TCGA-BRCA": "Breast_Mammary_Tissue",
    "TCGA-LUAD": "Lung",
    "TCGA-LUSC": "Lung",
    "TCGA-SKCM": "Skin_Sun_Exposed_Lower_leg",
    "TCGA-COAD": "Colon_Transverse",
    "TCGA-BLCA": "Bladder",
    "TCGA-HNSC": "Minor_Salivary_Gland",       # closest available
    "TCGA-KIRC": "Kidney_Cortex",
    "TCGA-LIHC": "Liver",
    "TCGA-PRAD": "Prostate",
    "TCGA-GBM":  "Brain_Cortex",
    "TCGA-OV":   "Ovary",
    "TCGA-UCEC": "Uterus",
    "TCGA-PAAD": "Pancreas",
    "TCGA-STAD": "Stomach",
}

GTEX_API_BASE = "https://gtexportal.org/api/v2"

RACE_MAP = {
    "white": "White",
    "black or african american": "Black / African American",
    "asian": "Asian",
    "american indian or alaska native": "Native American / Alaska Native",
    "native hawaiian or other pacific islander": "Pacific Islander",
    "not reported": "Not Reported",
}
POPULATION_GROUPS = ["All"] + list(RACE_MAP.values())

GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_FILES_ENDPOINT = f"{GDC_API_BASE}/files"
GDC_CASES_ENDPOINT = f"{GDC_API_BASE}/cases"
GDC_DATA_ENDPOINT = f"{GDC_API_BASE}/data"

# Plotly colors for families
_PX_COLORS = px.colors.qualitative.D3
FAMILY_LIST = sorted(set(info["family"] for info in RECEPTORS.values()))
FAMILY_COLORS = {f: _PX_COLORS[i % len(_PX_COLORS)] for i, f in enumerate(FAMILY_LIST)}


# ══════════════════════════════════════════════════════════════════════════════
# RECEPTOR → LIGAND MAPPING (TME-expressed ligands/enzymes that engage each receptor)
# Gene symbol, Ensembl ID, common name, and which receptor(s) it activates
# ══════════════════════════════════════════════════════════════════════════════

LIGANDS = {
    # PD-1 ligands
    "CD274":    {"ensembl": "ENSG00000120217", "label": "PD-L1",        "receptors": ["PDCD1"]},
    "PDCD1LG2": {"ensembl": "ENSG00000197646", "label": "PD-L2",       "receptors": ["PDCD1"]},
    # CTLA-4 ligands
    "CD80":     {"ensembl": "ENSG00000121594", "label": "B7-1 (CD80)",  "receptors": ["CTLA4"]},
    "CD86":     {"ensembl": "ENSG00000114013", "label": "B7-2 (CD86)",  "receptors": ["CTLA4"]},
    # LAG-3 ligands
    "FGL1":     {"ensembl": "ENSG00000104760", "label": "FGL1",         "receptors": ["LAG3"]},
    "LGALS3":   {"ensembl": "ENSG00000131981", "label": "Galectin-3",   "receptors": ["LAG3"]},
    "CLEC4G":   {"ensembl": "ENSG00000182566", "label": "LSECtin (CLEC4G)", "receptors": ["LAG3"]},
    # TIM-3 ligands
    "LGALS9":   {"ensembl": "ENSG00000168961", "label": "Galectin-9",   "receptors": ["HAVCR2"]},
    "HMGB1":    {"ensembl": "ENSG00000189403", "label": "HMGB1",        "receptors": ["HAVCR2"]},
    "CEACAM1":  {"ensembl": "ENSG00000079385", "label": "CEACAM1",      "receptors": ["HAVCR2"]},
    # TIGIT / CD96 ligands
    "PVR":      {"ensembl": "ENSG00000073008", "label": "CD155 (PVR)",  "receptors": ["TIGIT", "CD96"]},
    "PVRL2":    {"ensembl": "ENSG00000130202", "label": "CD112 (Nectin-2)", "receptors": ["TIGIT"]},
    "NECTIN3":  {"ensembl": "ENSG00000078043", "label": "CD113 (Nectin-3)", "receptors": ["TIGIT"]},
    "NECTIN1":  {"ensembl": "ENSG00000110400", "label": "CD111 (Nectin-1)", "receptors": ["CD96"]},
    # BTLA / CD160 ligand
    "TNFRSF14": {"ensembl": "ENSG00000157873", "label": "HVEM",         "receptors": ["BTLA", "CD160"]},
    # 2B4 (CD244) ligand
    "CD48":     {"ensembl": "ENSG00000117091", "label": "CD48",         "receptors": ["CD244"]},
    # VISTA ligand
    "VSIG3":    {"ensembl": "ENSG00000155659", "label": "VSIG3 (IGSF11)", "receptors": ["VSIR"]},
    # CD200R ligand
    "CD200":    {"ensembl": "ENSG00000091972", "label": "CD200",        "receptors": ["CD200R1"]},
    # KLRG1 ligands (cadherins)
    "CDH1":     {"ensembl": "ENSG00000039068", "label": "E-cadherin",   "receptors": ["KLRG1"]},
    "CDH2":     {"ensembl": "ENSG00000170558", "label": "N-cadherin",   "receptors": ["KLRG1"]},
    "CDH4":     {"ensembl": "ENSG00000179242", "label": "R-cadherin",   "receptors": ["KLRG1"]},
    # NKG2A ligand
    "HLA_E":    {"ensembl": "ENSG00000204592", "label": "HLA-E",        "receptors": ["KLRC1"]},
    # LAIR-1 ligand (collagen — use COL1A1 as representative)
    "COL1A1":   {"ensembl": "ENSG00000108821", "label": "Collagen-I (COL1A1)", "receptors": ["LAIR1"]},
    # LILRB1 ligand
    "HLA_G":    {"ensembl": "ENSG00000204632", "label": "HLA-G",        "receptors": ["LILRB1"]},
    # Adenosine pathway: CD39 + CD73 generate adenosine → engages A2A receptor
    "ENTPD1":   {"ensembl": "ENSG00000138185", "label": "CD39 (ENTPD1)", "receptors": ["ADORA2A"]},
    "NT5E":     {"ensembl": "ENSG00000135318", "label": "CD73 (NT5E)",  "receptors": ["ADORA2A"]},
    # IDO1 — immunosuppressive enzyme in TME (tryptophan catabolism)
    "IDO1":     {"ensembl": "ENSG00000131203", "label": "IDO1",          "receptors": []},
    # CD276 (B7-H3) — broad TME suppressive molecule
    "CD276":    {"ensembl": "ENSG00000103855", "label": "B7-H3 (CD276)", "receptors": []},
    # VTCN1 (B7-H4) — another TME suppressive ligand
    "VTCN1":    {"ensembl": "ENSG00000134258", "label": "B7-H4",         "receptors": []},
}

LIGAND_ENSEMBL = {gene: info["ensembl"] for gene, info in LIGANDS.items()}

# Combined gene list for downloading both receptors + ligands
ALL_GENE_ENSEMBL = {**GENE_ENSEMBL, **LIGAND_ENSEMBL}