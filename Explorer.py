"""
TCGA T Cell Exhaustion Receptor Explorer — Main Page
=====================================================
Network explorer + TME suppressive ligand activity.
Reads from local SQLite database (tcga_data.db) if available,
otherwise falls back to live GDC API calls.
Run download_to_db.py first for offline use.
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from itertools import combinations
import requests
import sqlite3
import io
import time
import os

from config import (
    RECEPTORS, GENE_ENSEMBL, LIGANDS, LIGAND_ENSEMBL, ALL_GENE_ENSEMBL,
    TCGA_PROJECTS, RACE_MAP, POPULATION_GROUPS,
    FAMILY_COLORS, GDC_CASES_ENDPOINT, GDC_FILES_ENDPOINT, GDC_DATA_ENDPOINT,
    TCGA_TO_GTEX_TISSUE, GTEX_API_BASE,
)

st.set_page_config(page_title="TCGA T Cell Exhaustion Explorer", page_icon="🧬",
                   layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
.main .block-container { padding-top: 2rem; max-width: 1200px; }
.gdc-badge { display: inline-block; padding: 0.25rem 0.75rem; border-radius: 20px;
    font-size: 0.75rem; font-weight: 600; }
.gdc-live { background: #d4edda; color: #155724; border: 1px solid #c3e6cb; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# DATABASE PATH
# ══════════════════════════════════════════════════════════════════════════════

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tcga_data.db")
DB_AVAILABLE = os.path.exists(DB_PATH)


# ══════════════════════════════════════════════════════════════════════════════
# SQLITE OFFLINE LOADER
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_data(ttl=None, show_spinner=False)
def load_from_db(project_id):
    """
    Load demographics and expression from local SQLite.
    Returns (expression_df, demographics_df).
    """
    conn = sqlite3.connect(DB_PATH)

    check = pd.read_sql_query(
        "SELECT COUNT(*) as n FROM expression WHERE project_id = ?",
        conn, params=(project_id,)
    )
    if check["n"].iloc[0] == 0:
        conn.close()
        raise GDCError(
            f"Project {project_id} not found in local database ({DB_PATH}). "
            f"Re-run: python download_to_db.py --projects {project_id}"
        )

    demo_df = pd.read_sql_query(
        "SELECT * FROM demographics WHERE project_id = ?",
        conn, params=(project_id,)
    )
    if not demo_df.empty:
        demo_df = demo_df.set_index("case_id")
        demo_df["race_label"] = demo_df["race"].map(RACE_MAP).fillna("Other / Unknown")

        def simplify_stage(s):
            if not isinstance(s, str) or s.lower() in ("not reported", "not available", ""):
                return "Not Reported"
            s_upper = s.upper().strip()
            for stage_num, patterns in [
                ("Stage IV",  ["STAGE IV", "IV"]),
                ("Stage III", ["STAGE III", "III"]),
                ("Stage II",  ["STAGE II", "II"]),
                ("Stage I",   ["STAGE I", "STAGE 1", " I"]),
            ]:
                for pat in patterns:
                    if pat in s_upper:
                        if pat == " I" and any(s_upper.endswith(x) for x in ["II", "III", "IV"]):
                            continue
                        return stage_num
            if "STAGE 0" in s_upper or "STAGE X" in s_upper or "IS" in s_upper:
                return "Stage 0/IS"
            return "Not Reported"
        demo_df["stage"] = demo_df["stage_raw"].apply(simplify_stage)

    expr_long = pd.read_sql_query(
        "SELECT case_id, gene, log2_tpm FROM expression WHERE project_id = ?",
        conn, params=(project_id,)
    )
    conn.close()

    if expr_long.empty:
        raise GDCError(f"No expression data for {project_id} in database.")

    expr_df = expr_long.pivot(index="case_id", columns="gene", values="log2_tpm")
    expr_df.index.name = "case_id"
    return expr_df, demo_df


# ══════════════════════════════════════════════════════════════════════════════
# NORMALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def upper_quartile_normalize(df):
    """
    Upper-quartile (UQ) between-sample normalization.
    TPM is already within-sample normalized (gene length + library size).
    UQ normalization corrects for remaining cross-sample technical variation
    by dividing each sample's values by its 75th percentile, then scaling
    to a common reference (median of all 75th percentiles).

    Input: log2(TPM+1) matrix (samples × genes)
    Output: UQ-normalized log2 matrix
    """
    # Back-transform to linear TPM
    linear = np.power(2, df) - 1

    # Compute 75th percentile per sample (row), ignoring zeros
    uq_per_sample = linear.apply(
        lambda row: np.percentile(row[row > 0], 75) if (row > 0).sum() > 0 else 1.0,
        axis=1,
    )

    # Scale factor: sample UQ / median UQ
    median_uq = uq_per_sample.median()
    scale_factors = median_uq / uq_per_sample.replace(0, median_uq)

    # Apply scaling and re-log
    normalized = linear.multiply(scale_factors, axis=0)
    return np.log2(normalized + 1)


# ══════════════════════════════════════════════════════════════════════════════
# GDC API
# ══════════════════════════════════════════════════════════════════════════════

class GDCError(Exception):
    pass

def _gdc_post(endpoint, payload, retries=3, timeout=60):
    last_err = None
    for attempt in range(retries):
        try:
            resp = requests.post(endpoint, json=payload,
                                 headers={"Content-Type": "application/json"}, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.Timeout as e:
            last_err = e; time.sleep(2 ** attempt) if attempt < retries - 1 else None
        except requests.exceptions.ConnectionError as e:
            last_err = e; time.sleep(2 ** attempt) if attempt < retries - 1 else None
        except requests.exceptions.HTTPError as e:
            raise GDCError(f"GDC HTTP {e.response.status_code}: {e.response.text[:500]}")
        except Exception as e:
            raise GDCError(f"Unexpected GDC error: {e}")
    raise GDCError(f"GDC unreachable after {retries} tries ({endpoint}). Last: {last_err}")


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_case_demographics(project_id):
    all_cases, bs, off = [], 500, 0
    while True:
        data = _gdc_post(GDC_CASES_ENDPOINT, {
            "filters": {"op": "and", "content": [
                {"op": "=", "content": {"field": "project.project_id", "value": project_id}},
                {"op": "=", "content": {"field": "files.analysis.workflow_type", "value": "STAR - Counts"}},
                {"op": "=", "content": {"field": "files.data_type", "value": "Gene Expression Quantification"}},
            ]},
            "fields": ("case_id,submitter_id,demographic.race,demographic.ethnicity,"
                       "demographic.gender,diagnoses.age_at_diagnosis,"
                       "diagnoses.ajcc_pathologic_stage,diagnoses.tumor_stage"),
            "size": bs, "from": off,
        })
        if "data" not in data: raise GDCError(f"Bad /cases response")
        hits = data["data"]["hits"]
        if not hits: break
        for h in hits:
            demo = h.get("demographic", {}); diag = h.get("diagnoses", [{}])
            d0 = diag[0] if diag else {}
            # Get stage — try ajcc_pathologic_stage first, fall back to tumor_stage
            raw_stage = d0.get("ajcc_pathologic_stage") or d0.get("tumor_stage") or "Not Reported"
            all_cases.append({"case_id": h["case_id"], "submitter_id": h.get("submitter_id",""),
                "race": (demo.get("race") or "not reported").lower(),
                "ethnicity": (demo.get("ethnicity") or "not reported").lower(),
                "gender": (demo.get("gender") or "not reported").lower(),
                "age_at_diagnosis_days": d0.get("age_at_diagnosis"),
                "stage_raw": raw_stage})
        off += bs
        if off >= data["data"]["pagination"]["total"]: break
    if not all_cases: raise GDCError(f"No cases for {project_id}")
    df = pd.DataFrame(all_cases).set_index("case_id")
    df["race_label"] = df["race"].map(RACE_MAP).fillna("Other / Unknown")

    # Normalize stage to simplified groups (Stage I, II, III, IV, Not Reported)
    def simplify_stage(s):
        if not isinstance(s, str) or s.lower() in ("not reported", "not available", ""):
            return "Not Reported"
        s_upper = s.upper().strip()
        # Match Roman numerals — extract the first stage number
        for stage_num, patterns in [
            ("Stage IV",  ["STAGE IV", "IV"]),
            ("Stage III", ["STAGE III", "III"]),
            ("Stage II",  ["STAGE II", "II"]),
            ("Stage I",   ["STAGE I", "STAGE 1", " I"]),
        ]:
            for pat in patterns:
                if pat in s_upper:
                    # Avoid false match: "Stage III" matching "III" in "Stage IIIA" is fine
                    # But "Stage I" shouldn't match "Stage III" — check it's not followed by I or V
                    if pat == " I" and any(s_upper.endswith(x) for x in ["II", "III", "IV"]):
                        continue
                    return stage_num
        # Handle "Stage 0" / "Stage X" / other
        if "STAGE 0" in s_upper or "STAGE X" in s_upper or "IS" in s_upper:
            return "Stage 0/IS"
        return "Not Reported"

    df["stage"] = df["stage_raw"].apply(simplify_stage)
    return df


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_expression_file_ids(project_id):
    all_files, bs, off = [], 500, 0
    while True:
        data = _gdc_post(GDC_FILES_ENDPOINT, {
            "filters": {"op": "and", "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
                {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
                {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
                {"op": "=", "content": {"field": "data_format", "value": "TSV"}},
                {"op": "=", "content": {"field": "access", "value": "open"}},
            ]},
            "fields": "file_id,file_name,cases.case_id", "size": bs, "from": off,
        })
        if "data" not in data: raise GDCError("Bad /files response")
        hits = data["data"]["hits"]
        if not hits: break
        for h in hits:
            for c in h.get("cases", []):
                all_files.append({"file_id": h["file_id"], "case_id": c["case_id"]})
        off += bs
        if off >= data["data"]["pagination"]["total"]: break
    if not all_files: raise GDCError(f"No files for {project_id}")
    return pd.DataFrame(all_files)


@st.cache_data(ttl=3600, show_spinner=False)
def _download_batch(_batch_ids_tuple):
    """Download a single batch of expression files from GDC. Cached per batch."""
    import tarfile
    fids = list(_batch_ids_tuple)
    target_ids = {eid.split(".")[0]: gene for gene, eid in ALL_GENE_ENSEMBL.items()}
    try:
        resp = requests.post(GDC_DATA_ENDPOINT, json={"ids": fids},
                             headers={"Content-Type": "application/json"}, timeout=120)
        resp.raise_for_status()
    except requests.exceptions.Timeout:
        raise GDCError(f"GDC /data timed out for batch of {len(fids)}.")
    except requests.exceptions.ConnectionError as e:
        raise GDCError(f"Cannot connect to GDC: {e}")
    except requests.exceptions.HTTPError as e:
        raise GDCError(f"GDC HTTP {e.response.status_code}")

    records = {}
    try:
        with tarfile.open(fileobj=io.BytesIO(resp.content), mode="r:gz") as tar:
            for member in tar.getmembers():
                if not member.name.endswith(".tsv"): continue
                parts = member.name.split("/")
                fid = parts[0] if len(parts) > 1 else member.name
                f = tar.extractfile(member)
                if not f: continue
                try: tsv = pd.read_csv(f, sep="\t", comment="#")
                except: continue
                if "gene_id" not in tsv.columns: continue
                tpm_col = next((c for c in ["tpm_unstranded","fpkm_unstranded","unstranded"] if c in tsv.columns), None)
                if not tpm_col: continue
                tsv["gid"] = tsv["gene_id"].str.split(".").str[0]
                sub = tsv.loc[tsv["gid"].isin(target_ids), ["gid", tpm_col]].copy()
                sub["sym"] = sub["gid"].map(target_ids)
                if not sub.empty:
                    records[fid] = {row["sym"]: np.log2(float(row[tpm_col]) + 1) for _, row in sub.iterrows()}
    except Exception as e:
        raise GDCError(f"Failed to parse batch: {e}")
    return records


def download_expression_batched(file_ids, max_files=200, batch_size=25, progress_bar=None):
    """
    Download expression data in batches with progress updates.
    Each batch is independently cached via _download_batch.
    """
    fids = file_ids[:max_files]
    total = len(fids)

    # Split into batches
    batches = [fids[i:i+batch_size] for i in range(0, total, batch_size)]

    all_records = {}
    downloaded = 0

    for bi, batch in enumerate(batches):
        # Update progress: map batch progress to 25%–75% of overall bar
        if progress_bar is not None:
            pct = 25 + int(50 * downloaded / total)
            progress_bar.progress(pct, text=f"Downloading expression profiles: {downloaded}/{total} complete...")

        batch_records = _download_batch(tuple(batch))
        all_records.update(batch_records)
        downloaded += len(batch)

    if progress_bar is not None:
        progress_bar.progress(75, text=f"Downloaded {downloaded}/{total} profiles. Extracted {len(all_records)} samples.")

    if not all_records:
        raise GDCError(f"No expression data extracted from {total} files.")

    expr = pd.DataFrame.from_dict(all_records, orient="index")
    expr.index.name = "file_id"
    for g in ALL_GENE_ENSEMBL:
        if g not in expr.columns:
            expr[g] = np.nan
    return expr


def get_full_dataset(project_id, max_samples=200, progress_bar=None):
    """
    Load data for a project. Tries local SQLite DB first,
    falls back to live GDC API if DB unavailable.
    Returns (receptor_df, ligand_df, demo_df, source_label).
    """
    def _update(pct, text):
        if progress_bar is not None:
            progress_bar.progress(pct, text=text)

    # ── Try SQLite first ─────────────────────────────────────────────────
    if DB_AVAILABLE:
        try:
            _update(10, f"Loading {project_id} from local database...")
            expr, demo = load_from_db(project_id)
            _update(50, f"Loaded {len(expr)} samples. Filling missing values...")

            for c in expr.columns:
                expr[c] = expr[c].fillna(expr[c].median())
            _update(70, "Applying upper-quartile normalization...")

            expr = upper_quartile_normalize(expr)
            _update(90, "Splitting receptor and ligand matrices...")

            rcols = [c for c in expr.columns if c in RECEPTORS]
            lcols = [c for c in expr.columns if c in LIGANDS]

            if len(expr) < 10:
                raise GDCError(f"Only {len(expr)} samples in DB for {project_id}.")

            _update(100, "Done (from local database).")
            return expr[rcols], expr[lcols], demo, "LOCAL DB"

        except GDCError:
            _update(15, f"{project_id} not in local DB. Trying GDC API...")

    # ── Fall back to GDC API ─────────────────────────────────────────────
    _update(5, "Querying GDC for patient demographics & stage...")
    demo = fetch_case_demographics(project_id)
    _update(15, f"Found {len(demo)} cases. Locating expression files...")

    files = fetch_expression_file_ids(project_id).drop_duplicates(subset="case_id", keep="first")
    n_files = min(len(files), max_samples)
    _update(20, f"Found {len(files)} files. Preparing to download {n_files} profiles...")

    fids = files["file_id"].tolist()[:max_samples]
    expr = download_expression_batched(fids, max_files=max_samples, batch_size=25, progress_bar=progress_bar)
    _update(78, f"Extracted {len(expr)} samples. Mapping to cases...")

    f2c = dict(zip(files["file_id"], files["case_id"]))
    expr["case_id"] = expr.index.map(lambda x: f2c.get(x))
    expr = expr.dropna(subset=["case_id"]).set_index("case_id")
    expr = expr.dropna(how="all")
    _update(82, f"{len(expr)} valid samples. Filling missing values...")

    if len(expr) < 10:
        raise GDCError(f"Only {len(expr)} valid samples for {project_id}.")

    for c in expr.columns:
        expr[c] = expr[c].fillna(expr[c].median())
    _update(88, "Applying upper-quartile normalization...")

    expr = upper_quartile_normalize(expr)
    _update(95, "Splitting receptor and ligand matrices...")

    rcols = [c for c in expr.columns if c in RECEPTORS]
    lcols = [c for c in expr.columns if c in LIGANDS]
    _update(100, "Done (from GDC API).")

    return expr[rcols], expr[lcols], demo, "GDC API"


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# GTEx BASELINE (fetched on-demand, cached indefinitely)
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_data(ttl=None, show_spinner=False)
def fetch_gtex_baseline(project_id):
    """
    Fetch median gene expression from GTEx V2 API for the healthy tissue
    matching this TCGA cancer type. Cached indefinitely after first fetch.

    Returns (pd.Series of gene→median_tpm, gtex_tissue_name).
    Returns empty Series if GTEx is unreachable or no mapping exists.
    """
    gtex_tissue = TCGA_TO_GTEX_TISSUE.get(project_id)
    if not gtex_tissue:
        return pd.Series(dtype=float), ""

    baselines = {}
    for gene_symbol, ensembl_id in ALL_GENE_ENSEMBL.items():
        try:
            resp = requests.get(
                f"{GTEX_API_BASE}/expression/medianGeneExpression",
                params={
                    "datasetId": "gtex_v8",
                    "gencodeId": ensembl_id,
                    "tissueSiteDetailId": gtex_tissue,
                },
                timeout=15,
            )
            if resp.status_code == 200:
                data = resp.json()
                items = data.get("data", [])
                if items:
                    baselines[gene_symbol] = items[0].get("median", 0)
                else:
                    baselines[gene_symbol] = 0
            else:
                baselines[gene_symbol] = 0
        except Exception:
            baselines[gene_symbol] = 0
        time.sleep(0.05)  # rate limit politeness

    return pd.Series(baselines), gtex_tissue


def compute_tissue_baseline(ligand_df, gtex_baseline=None):
    """
    Get per-ligand tissue baseline in linear TPM.

    If GTEx baseline is available (from DB), use it — this is the median
    expression in healthy normal tissue (e.g., GTEx Lung for TCGA-LUAD).

    If not available, fall back to cohort median (less ideal but functional).

    Returns a Series: gene → baseline linear TPM.
    """
    baselines = {}
    for lg in ligand_df.columns:
        if lg in LIGANDS:
            if gtex_baseline is not None and not gtex_baseline.empty and lg in gtex_baseline.index:
                # GTEx provides median TPM directly (already linear)
                baselines[lg] = max(0, gtex_baseline[lg])
            else:
                # Fallback: cohort median
                linear = np.power(2, ligand_df[lg]) - 1
                baselines[lg] = linear.median()
    return pd.Series(baselines)

def compute_ligand_activation_scores(ligand_df, active_receptors=None):
    """
    For each receptor, compute per-patient Ligand Activation Score.
    Only includes receptors in active_receptors (or all if None).
    """
    receptors_to_use = active_receptors if active_receptors is not None else RECEPTORS
    scores = {}
    for receptor_gene, rinfo in receptors_to_use.items():
        lig_genes = [lg for lg, linfo in LIGANDS.items()
                     if receptor_gene in linfo.get("receptors", []) and lg in ligand_df.columns]
        if lig_genes:
            # Back-transform each ligand to linear TPM, sum, re-log
            linear_sum = (np.power(2, ligand_df[lig_genes]) - 1).sum(axis=1)
            scores[receptor_gene] = np.log2(linear_sum + 1)
    return pd.DataFrame(scores)


def compute_coexpression(ligand_df, thresh=0.3, p_thresh=0.05):
    """
    Compute co-activation of inhibitory pathways via LIGAND correlation.

    For each pair of receptors (A, B):
      - Ligand Activation Score_A = Σ expression of all A's ligands per patient
      - Ligand Activation Score_B = Σ expression of all B's ligands per patient
      - Edge weight = |Spearman(Score_A, Score_B)|

    Skips pairs where the ligand gene sets are identical (would give trivial ρ=1.0).
    """
    activation = compute_ligand_activation_scores(ligand_df)
    receptors_with_ligands = activation.columns.tolist()

    # Pre-compute which ligand genes map to each receptor (for labels + overlap check)
    receptor_ligand_genes = {}
    receptor_ligand_labels = {}
    for rg in receptors_with_ligands:
        genes = [lg for lg, linfo in LIGANDS.items()
                 if rg in linfo.get("receptors", []) and lg in ligand_df.columns]
        receptor_ligand_genes[rg] = set(genes)
        receptor_ligand_labels[rg] = ", ".join(
            LIGANDS[lg]["label"] for lg in genes
        ) or "none"

    edges = []
    for g1, g2 in combinations(receptors_with_ligands, 2):
        # Skip if ligand gene sets are identical (e.g., BTLA & CD160 both use only HVEM)
        if receptor_ligand_genes[g1] == receptor_ligand_genes[g2]:
            continue

        valid = activation[[g1, g2]].dropna()
        if len(valid) < 10:
            continue
        rho, pv = stats.spearmanr(valid[g1], valid[g2])
        if abs(rho) >= thresh and pv < p_thresh:
            edges.append({
                "source": g1, "target": g2,
                "weight": abs(rho), "rho": rho, "pval": pv,
                "co_prob": (abs(rho) + 1) / 2,
                "ligands_A": receptor_ligand_labels[g1],
                "ligands_B": receptor_ligand_labels[g2],
                "mean_score_A": activation[g1].mean(),
                "mean_score_B": activation[g2].mean(),
            })
    return pd.DataFrame(edges)


def find_shared_ligand_pairs(ligand_df):
    """
    Find receptor pairs that share one or more ligands.
    Returns list of dicts with shared ligand info.
    These get drawn as dashed edges (distinct from correlation edges).
    """
    # Build receptor → ligand gene set mapping
    receptor_ligs = {}
    for rg in RECEPTORS:
        genes = [lg for lg, linfo in LIGANDS.items()
                 if rg in linfo.get("receptors", []) and lg in ligand_df.columns]
        if genes:
            receptor_ligs[rg] = set(genes)

    shared_pairs = []
    for g1, g2 in combinations(receptor_ligs.keys(), 2):
        overlap = receptor_ligs[g1] & receptor_ligs[g2]
        if overlap:
            shared_labels = ", ".join(LIGANDS[lg]["label"] for lg in overlap)
            shared_pairs.append({
                "source": g1, "target": g2,
                "shared_ligands": shared_labels,
                "shared_genes": overlap,
                "is_identical": receptor_ligs[g1] == receptor_ligs[g2],
            })
    return shared_pairs


def build_graph(edge_df, shared_ligand_pairs=None, active_receptors=None):
    """Build network. Only includes receptors in active_receptors (or all if None)."""
    G = nx.Graph()
    receptors_to_use = active_receptors if active_receptors is not None else RECEPTORS
    for g, info in receptors_to_use.items():
        G.add_node(g, **info)

    # Track which pairs share ligands and whether it's a perfect overlap
    shared_map = {}  # (source, target) → {"shared_ligands": ..., "is_identical": bool}
    if shared_ligand_pairs:
        for sp in shared_ligand_pairs:
            key = tuple(sorted([sp["source"], sp["target"]]))
            shared_map[key] = sp

    # Add correlation edges, tagging them if they also share ligands
    for _, r in edge_df.iterrows():
        s, t = r["source"], r["target"]
        key = tuple(sorted([s, t]))
        sp = shared_map.pop(key, None)  # consume from shared_map if present

        if sp:
            # Has both correlation AND shared ligands
            edge_type = "both"
            shared_ligands = sp["shared_ligands"]
        else:
            edge_type = "correlation"
            shared_ligands = ""

        G.add_edge(s, t,
                    weight=r["weight"], rho=r["rho"], pval=r["pval"], co_prob=r["co_prob"],
                    ligands_A=r.get("ligands_A", ""), ligands_B=r.get("ligands_B", ""),
                    mean_score_A=r.get("mean_score_A", 0), mean_score_B=r.get("mean_score_B", 0),
                    edge_type=edge_type, shared_ligands=shared_ligands)

    # Remaining shared_map entries have no correlation edge — add as shared-only
    for key, sp in shared_map.items():
        s, t = sp["source"], sp["target"]
        # Determine style: identical overlap = dotted black, partial = dashed orange
        etype = "shared_identical" if sp["is_identical"] else "shared_partial"
        G.add_edge(s, t,
                    weight=0, rho=0, pval=1, co_prob=0.5,
                    ligands_A="", ligands_B="",
                    mean_score_A=0, mean_score_B=0,
                    edge_type=etype,
                    shared_ligands=sp["shared_ligands"])
    return G

def hierarchical_layout(G):
    """
    Radial tree layout grouped by family.
    Outer ring for members, with generous angular spacing to reduce overlap.
    """
    fams = {}
    for n, d in G.nodes(data=True):
        fams.setdefault(d.get("family", "Other"), []).append(n)

    pos = {}
    nf = len(fams)
    # Give each family an equal angular slice
    slice_angle = (2 * np.pi) / nf if nf > 0 else 2 * np.pi

    for fi, (fam, mems) in enumerate(sorted(fams.items())):
        center_angle = (slice_angle * fi) - np.pi / 2
        nm = len(mems)

        if nm == 1:
            pos[mems[0]] = (5.0 * np.cos(center_angle), 5.0 * np.sin(center_angle))
        else:
            # Spread members within the family's angular slice (use 70% of slice)
            member_spread = slice_angle * 0.7
            for mi, m in enumerate(sorted(mems)):
                offset = member_spread * (mi - (nm - 1) / 2) / max(nm - 1, 1)
                a = center_angle + offset
                # Alternate radius slightly to reduce label overlap
                r = 5.0 + (0.6 if mi % 2 == 1 else 0)
                pos[m] = (r * np.cos(a), r * np.sin(a))

    return pos


def compute_tme_suppression(receptor_df, ligand_df, gtex_baseline=None):
    """
    For each receptor, compute log₂FC for the receptor itself and each of its ligands.
    Also compute receptor-ligand correlation (Spearman).
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    rows = []
    for gene, rinfo in RECEPTORS.items():
        if gene not in receptor_df.columns:
            continue
        matching_ligands = []
        for lgene, linfo in LIGANDS.items():
            if gene in linfo.get("receptors", []) and lgene in ligand_df.columns:
                matching_ligands.append((lgene, linfo))

        r_linear = np.power(2, receptor_df[gene]) - 1
        r_bl = gtex_baseline[gene] if (gtex_baseline is not None and not gtex_baseline.empty and gene in gtex_baseline.index) else max(r_linear.median(), 0.1)
        r_bl = max(r_bl, 0.1)
        r_fc = max(0, np.log2((r_linear.mean() + 0.1) / r_bl))

        for lgene, linfo in matching_ligands:
            l_linear = np.power(2, ligand_df[lgene]) - 1
            l_bl = max(baseline.get(lgene, 0), 0.1)
            l_fc = max(0, np.log2((l_linear.mean() + 0.1) / l_bl))

            valid = pd.concat([receptor_df[gene], ligand_df[lgene]], axis=1).dropna()
            if len(valid) >= 10:
                rho, pval = stats.spearmanr(valid.iloc[:,0], valid.iloc[:,1])
            else:
                rho, pval = np.nan, np.nan

            rows.append({
                "Receptor": rinfo["label"],
                "Receptor Gene": gene,
                "Ligand": linfo["label"],
                "Ligand Gene": lgene,
                "Receptor log₂FC": r_fc,
                "Ligand log₂FC": l_fc,
                "Suppressive Score": r_fc + l_fc,
                "R-L Correlation (ρ)": rho,
                "R-L p-value": pval,
            })

    for lgene, linfo in LIGANDS.items():
        if not linfo["receptors"] and lgene in ligand_df.columns:
            l_linear = np.power(2, ligand_df[lgene]) - 1
            l_bl = max(baseline.get(lgene, 0), 0.1)
            l_fc = max(0, np.log2((l_linear.mean() + 0.1) / l_bl))
            rows.append({
                "Receptor": "TME (general)",
                "Receptor Gene": "—",
                "Ligand": linfo["label"],
                "Ligand Gene": lgene,
                "Receptor log₂FC": np.nan,
                "Ligand log₂FC": l_fc,
                "Suppressive Score": l_fc,
                "R-L Correlation (ρ)": np.nan,
                "R-L p-value": np.nan,
            })

    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def create_network(G, pos, title="", activation_scores=None, gtex_baseline=None, ligand_df=None):
    """
    Network graph. Node size = potential activation (total ligand log₂FC vs normal tissue).
    """
    fig = go.Figure()
    MIN_W, MAX_W = 1.0, 14.0
    aw = [d["weight"] for _,_,d in G.edges(data=True)]
    wn = min(aw) if aw else 0; wx = max(aw) if aw else 1; wr = wx-wn if wx>wn else 1

    # Compute log₂FC per receptor for node sizing
    _baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline) if ligand_df is not None else {}
    mean_log2fc = {}
    for rg in RECEPTORS:
        lgs = [lg for lg, li in LIGANDS.items()
               if rg in li.get("receptors", []) and ligand_df is not None and lg in ligand_df.columns]
        if lgs and ligand_df is not None:
            total_fc = 0
            for lg in lgs:
                lin = np.power(2, ligand_df[lg]) - 1
                bl = max(_baseline.get(lg, 0), 0.1)
                total_fc += max(0, np.log2((lin.mean() + 0.1) / bl))
            mean_log2fc[rg] = total_fc

    fv = list(mean_log2fc.values()) if mean_log2fc else [0]
    fmin = min(fv); fmax = max(fv); frng = fmax - fmin if fmax > fmin else 1
    MIN_NODE, MAX_NODE = 15, 55
    def fc_to_size(gene):
        v = mean_log2fc.get(gene, 0)
        return MIN_NODE + ((v - fmin) / frng) * (MAX_NODE - MIN_NODE)

    # Categorize edges
    all_edges = list(G.edges(data=True))
    corr_edges = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "correlation"]
    both_edges = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "both"]
    shared_identical = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "shared_identical"]
    shared_partial = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "shared_partial"]

    # Weight scaling uses only edges that have a correlation
    weighted_edges = corr_edges + both_edges
    aw = [d["weight"] for _,_,d in weighted_edges]
    wn = min(aw) if aw else 0; wx = max(aw) if aw else 1; wr = wx-wn if wx>wn else 1

    def _hover_corr(lu, lv, d):
        """Build hover text for a correlation edge."""
        h = (f"<b>{lu} ↔ {lv}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Co-activation ρ:</b> {d['rho']:.4f}<br>"
             f"<b>Weight:</b> {d['weight']:.4f}<br>"
             f"<b>p-value:</b> {d['pval']:.2e}<br>"
             f"<b>{lu} ligands:</b> {d.get('ligands_A','')}<br>"
             f"<b>{lu} mean ligand score:</b> {d.get('mean_score_A',0):.3f}<br>"
             f"<b>{lv} ligands:</b> {d.get('ligands_B','')}<br>"
             f"<b>{lv} mean ligand score:</b> {d.get('mean_score_B',0):.3f}")
        shared = d.get("shared_ligands", "")
        if shared:
            h += f"<br><b>⚠ Also shares ligand:</b> {shared}"
        return h

    def _add_edge_traces(fig, u, v, d, line_dict, hover_text):
        x0,y0=pos[u]; x1,y1=pos[v]
        fig.add_trace(go.Scatter(x=[x0,x1,None],y=[y0,y1,None],mode="lines",
            line=line_dict,hoverinfo="text",hovertext=[hover_text,hover_text,None],showlegend=False))
        fig.add_trace(go.Scatter(x=[(x0+x1)/2],y=[(y0+y1)/2],mode="markers",
            marker=dict(size=max(24, line_dict.get("width",3)*2.5),color="rgba(0,0,0,0)"),
            hoverinfo="text",hovertext=hover_text,showlegend=False))

    # ── 1. Shared-identical edges: black dotted (same ligand set, ρ=1 trivially) ──
    for u,v,d in shared_identical:
        lu=RECEPTORS[u]["label"]; lv=RECEPTORS[v]["label"]
        shared = d.get("shared_ligands","?")
        h = (f"<b>{lu} ↔ {lv}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Identical ligand set:</b> {shared}<br>"
             f"Both receptors bind exactly the same TME ligand(s).<br>"
             f"Co-activation is guaranteed (ρ = 1.0 by definition).")
        _add_edge_traces(fig, u, v, d,
            line_dict=dict(width=3, color="rgba(0,0,0,0.45)", dash="dot"), hover_text=h)

    # ── 2. Shared-partial edges (no correlation): dashed orange ───────────
    for u,v,d in shared_partial:
        lu=RECEPTORS[u]["label"]; lv=RECEPTORS[v]["label"]
        shared = d.get("shared_ligands","?")
        h = (f"<b>{lu} ↔ {lv}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Shared ligand:</b> {shared}<br>"
             f"These receptors share some (not all) TME ligands.<br>"
             f"No significant independent co-activation at current threshold.")
        _add_edge_traces(fig, u, v, d,
            line_dict=dict(width=3, color="rgba(255,165,0,0.5)", dash="dash"), hover_text=h)

    # ── 3. Correlation-only edges: blue if ρ>0, red if ρ<0, thickness ∝ |ρ| ──
    for u,v,d in corr_edges:
        lu=RECEPTORS[u]["label"]; lv=RECEPTORS[v]["label"]
        w=d["weight"]; rho=d["rho"]; t=(w-wn)/wr; lp=MIN_W+(t**1.5)*(MAX_W-MIN_W); op=0.35+0.55*t
        ec = f"rgba(99,110,250,{op:.2f})" if rho >= 0 else f"rgba(239,85,59,{op:.2f})"
        h = _hover_corr(lu, lv, d)
        _add_edge_traces(fig, u, v, d,
            line_dict=dict(width=lp, color=ec), hover_text=h)

    # ── 4. Both (correlation + shared ligands): dashed, blue/red ∝ ρ ──────
    for u,v,d in both_edges:
        lu=RECEPTORS[u]["label"]; lv=RECEPTORS[v]["label"]
        w=d["weight"]; rho=d["rho"]; t=(w-wn)/wr; lp=MIN_W+(t**1.5)*(MAX_W-MIN_W); op=0.35+0.55*t
        ec = f"rgba(99,110,250,{op:.2f})" if rho >= 0 else f"rgba(239,85,59,{op:.2f})"
        h = _hover_corr(lu, lv, d)
        _add_edge_traces(fig, u, v, d,
            line_dict=dict(width=lp, color=ec, dash="dash"), hover_text=h)

    # ── Legend entries for edge types ──────────────────────────────────────
    has_pos = any(d["rho"] >= 0 for _,_,d in corr_edges + both_edges) if (corr_edges or both_edges) else False
    has_neg = any(d["rho"] < 0 for _,_,d in corr_edges + both_edges) if (corr_edges or both_edges) else False
    if has_pos:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(99,110,250,0.7)"),
            name="Co-activation (+ρ)",showlegend=True))
    if has_neg:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(239,85,59,0.7)"),
            name="Inverse correlation (−ρ)",showlegend=True))
    if both_edges:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(99,110,250,0.7)",dash="dash"),
            name="+ shared ligand (dashed)",showlegend=True))
    if shared_identical:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(0,0,0,0.45)",dash="dot"),
            name="Identical ligand (ρ=1)",showlegend=True))
    if shared_partial:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(255,165,0,0.5)",dash="dash"),
            name="Shared ligand (no sig. ρ)",showlegend=True))
    bf={}
    for n,d in G.nodes(data=True):
        f=d.get("family","Other"); bf.setdefault(f,{"x":[],"y":[],"text":[],"hover":[],"sizes":[]})
        x,y=pos[n]; deg=G.degree(n); lb=d.get("label",n); ds=d.get("desc","")
        aw2=np.mean([d2["weight"] for _,_,d2 in G.edges(n,data=True)]) if deg>0 else 0
        fc = mean_log2fc.get(n, 0)
        node_sz = fc_to_size(n)

        lig_list = [LIGANDS[lg]["label"] for lg in LIGANDS
                    if n in LIGANDS[lg].get("receptors",[]) and ligand_df is not None and lg in ligand_df.columns]

        bf[f]["x"].append(x); bf[f]["y"].append(y); bf[f]["text"].append(lb)
        bf[f]["sizes"].append(node_sz)
        bf[f]["hover"].append(
            f"<b>{lb}</b> ({n})<br>{ds}<br>Family: {f}<br>"
            f"━━━━━━━━━━━━━━━<br>"
            f"<b>Potential activation (log₂FC):</b> {fc:.2f}<br>"
            f"<b>Ligands:</b> {', '.join(lig_list) if lig_list else 'none'}<br>"
            f"<b>Connections:</b> {deg}<br>"
            f"<b>Avg edge weight:</b> {aw2:.4f}"
        )
    for f,v in bf.items():
        fig.add_trace(go.Scatter(x=v["x"],y=v["y"],mode="markers+text",
            marker=dict(size=v["sizes"],color=FAMILY_COLORS.get(f),line=dict(width=2,color="white")),
            text=v["text"],textposition="top center",textfont=dict(size=11, family="Arial"),
            hoverinfo="text",hovertext=v["hover"],name=f,legendgroup=f))
    fig.update_layout(title=dict(text=title,x=0.5),
        xaxis=dict(showgrid=False,zeroline=False,showticklabels=False),
        yaxis=dict(showgrid=False,zeroline=False,showticklabels=False,scaleanchor="x",scaleratio=1),
        hovermode="closest",height=800,margin=dict(l=60,r=60,t=60,b=60),template="plotly",dragmode="pan")
    return fig

def create_barplot(df, project_id, gtex_baseline=None, ligand_df=None):
    """
    Per-receptor ligand activation distribution showing individual ligand contributions.
    Each receptor row has overlaid violins colored by ligand.
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline) if ligand_df is not None else {}

    # Build per-patient, per-ligand, per-receptor data
    rows = []
    receptor_medians = {}
    for gene, rinfo in RECEPTORS.items():
        if gene not in df.columns:
            continue
        lig_genes = [(lg, linfo) for lg, linfo in LIGANDS.items()
                     if gene in linfo.get("receptors", []) and ligand_df is not None and lg in ligand_df.columns]
        if not lig_genes or ligand_df is None:
            continue

        # Per-ligand log₂FC for each patient
        for lg, linfo in lig_genes:
            linear = np.power(2, ligand_df[lg]) - 1
            bl = max(baseline.get(lg, 0), 0.1)
            patient_fc = np.log2((linear + 0.1) / bl)
            patient_fc = patient_fc.loc[patient_fc.index.isin(df.index)]

            for val in patient_fc.values:
                rows.append({
                    "Receptor": rinfo["label"],
                    "Ligand": linfo["label"],
                    "Family": rinfo["family"],
                    "log₂FC": val,
                })

        # Also compute total for sorting
        linear_sum = pd.Series(0.0, index=ligand_df.index)
        baseline_sum = 0.0
        for lg, linfo in lig_genes:
            linear_sum += np.power(2, ligand_df[lg]) - 1
            baseline_sum += max(baseline.get(lg, 0), 0.1)
        total_fc = np.log2((linear_sum + 0.1) / max(baseline_sum, 0.1))
        total_fc = total_fc.loc[total_fc.index.isin(df.index)]
        receptor_medians[rinfo["label"]] = total_fc.median()

    if not rows:
        fig = go.Figure()
        fig.update_layout(title="No receptor data available")
        return fig

    rdf = pd.DataFrame(rows)

    # Sort receptors by median total activation descending
    receptor_order = sorted(receptor_medians.keys(), key=lambda r: receptor_medians[r], reverse=True)

    # Consistent ligand colors
    all_ligands = rdf["Ligand"].unique().tolist()
    colors = px.colors.qualitative.Plotly + px.colors.qualitative.D3
    ligand_colors = {lig: colors[i % len(colors)] for i, lig in enumerate(all_ligands)}

    fig = go.Figure()

    # Track which ligands we've added to legend
    legend_shown = set()

    for receptor_name in reversed(receptor_order):
        sub = rdf[rdf["Receptor"] == receptor_name]
        ligands_in_receptor = sub["Ligand"].unique()

        for lig_name in ligands_in_receptor:
            lig_sub = sub[sub["Ligand"] == lig_name]
            show_legend = lig_name not in legend_shown
            legend_shown.add(lig_name)

            fig.add_trace(go.Violin(
                x=lig_sub["log₂FC"],
                y=[receptor_name] * len(lig_sub),
                orientation="h",
                side="positive",
                width=1.5,
                line_color=ligand_colors[lig_name],
                fillcolor=ligand_colors[lig_name],
                opacity=0.5,
                meanline_visible=True,
                name=lig_name,
                legendgroup=lig_name,
                showlegend=show_legend,
                hovertemplate=(
                    f"<b>{receptor_name}</b> ← {lig_name}<br>"
                    "log₂FC: %{x:.2f}<extra></extra>"
                ),
                scalemode="width",
            ))

    fig.update_layout(
        title=f"Ligand Activation Distribution per Receptor — {TCGA_PROJECTS[project_id]}",
        xaxis_title="log₂(Ligand TPM / Normal)",
        yaxis=dict(categoryorder="array", categoryarray=receptor_order),
        height=max(600, len(receptor_order) * 55),
        margin=dict(l=140, r=20, t=50, b=50),
        template="plotly",
        violingap=0.05,
        violinmode="overlay",
        legend_title="TME Ligand",
    )

    # Zero reference line
    fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5,
                  annotation_text="Normal", annotation_position="top")

    return fig

def create_corrmatrix(ligand_df):
    """Spearman correlation of ligand activation scores (not receptor expression)."""
    activation = compute_ligand_activation_scores(ligand_df)
    if activation.empty:
        return go.Figure()
    lb=[RECEPTORS[g]["label"] for g in activation.columns]; corr=activation.corr(method="spearman")
    fig=go.Figure(go.Heatmap(z=corr.values,x=lb,y=lb,colorscale="RdBu_r",zmid=0,zmin=-1,zmax=1,
        hovertemplate="<b>%{x} × %{y}</b><br>Ligand-activation ρ=%{z:.3f}<extra></extra>",colorbar=dict(title="ρ")))
    fig.update_layout(title="Ligand-Activation Correlation Matrix (TME suppressive pressure)",
        xaxis=dict(tickangle=45,tickfont=dict(size=9)),yaxis=dict(autorange="reversed",tickfont=dict(size=9)),
        height=700,margin=dict(l=140,r=20,t=50,b=140),template="plotly")
    return fig


def create_tme_heatmap(tme_df):
    """Heatmap: rows=receptors, cols=ligands, values=suppressive score."""
    pivot = tme_df.pivot_table(index="Receptor", columns="Ligand", values="Suppressive Score", aggfunc="first")
    pivot = pivot.dropna(how="all").dropna(axis=1, how="all")
    if pivot.empty:
        return None
    fig = go.Figure(go.Heatmap(
        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale="YlOrRd", hovertemplate="<b>%{y} ← %{x}</b><br>Score: %{z:.3f}<extra></extra>",
        colorbar=dict(title="Score"),
    ))
    fig.update_layout(title="TME Suppressive Score (log₂FC Receptor + Ligand)",
        xaxis=dict(tickangle=45, tickfont=dict(size=10), title="Ligand (TME)"),
        yaxis=dict(tickfont=dict(size=10), title="Receptor (T cell)", autorange="reversed"),
        height=600, margin=dict(l=140, r=20, t=50, b=120), template="plotly")
    return fig


def create_ligand_barplot(ligand_df, gtex_baseline=None):
    """Bar plot of ligand log₂ fold-change over normal tissue with Q1–Q3 variability."""
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    gc = [c for c in ligand_df.columns if c in LIGANDS]

    rows = []
    for g in gc:
        linear = np.power(2, ligand_df[g]) - 1
        bl = max(baseline.get(g, 0), 0.1)
        mean_tpm = linear.mean()
        log2fc_mean = max(0, np.log2((mean_tpm + 0.1) / bl))

        # Per-patient log₂FC for variability
        patient_fc = np.log2((linear + 0.1) / bl).clip(lower=0)
        q1 = patient_fc.quantile(0.25)
        q3 = patient_fc.quantile(0.75)

        rows.append({
            "gene": g, "label": LIGANDS[g]["label"],
            "log2fc": log2fc_mean, "q1": q1, "q3": q3,
            "iqr_lo": max(0, log2fc_mean - q1),
            "iqr_hi": max(0, q3 - log2fc_mean),
            "fold": mean_tpm / bl, "mean_tpm": mean_tpm, "baseline": bl,
        })

    df = pd.DataFrame(rows).sort_values("log2fc", ascending=True)

    fig = go.Figure(go.Bar(
        x=df["log2fc"].values, y=df["label"].values, orientation="h",
        marker=dict(color="indianred"),
        error_x=dict(
            type="data", symmetric=False,
            array=df["iqr_hi"].values, arrayminus=df["iqr_lo"].values,
            thickness=1.5, width=3,
        ),
        customdata=np.column_stack([
            df["fold"].values, df["mean_tpm"].values, df["baseline"].values,
            df["q1"].values, df["q3"].values,
        ]),
        hovertemplate=(
            "<b>%{y}</b><br>"
            "log₂FC: %{x:.2f} (= %{customdata[0]:.1f}× normal)<br>"
            "Q1: %{customdata[3]:.2f} | Q3: %{customdata[4]:.2f}<br>"
            "Tumor: %{customdata[1]:.1f} TPM | Normal: %{customdata[2]:.1f} TPM"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        title="TME Ligand Upregulation vs Normal Tissue",
        xaxis_title="log₂(Tumor / Normal)",
        height=600, margin=dict(l=180, r=20, t=50, b=40), template="plotly",
    )
    return fig


def create_receptor_activation_chart(ligand_df, project_id, gtex_baseline=None, active_receptors=None):
    """
    Two charts using log₂ fold-change over normal tissue.
    Only includes receptors in active_receptors.
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    baseline_source = "GTEx normal tissue" if (gtex_baseline is not None and not gtex_baseline.empty) else "cohort median"
    receptors_to_use = active_receptors if active_receptors is not None else RECEPTORS

    rows = []
    for receptor_gene, rinfo in receptors_to_use.items():
        lig_genes = [
            (lg, linfo) for lg, linfo in LIGANDS.items()
            if receptor_gene in linfo.get("receptors", []) and lg in ligand_df.columns
        ]
        if not lig_genes:
            continue
        for lg, linfo in lig_genes:
            linear_vals = np.power(2, ligand_df[lg]) - 1
            bl = max(baseline.get(lg, 0), 0.1)

            # log₂FC per patient
            log2fc_vals = np.log2((linear_vals + 0.1) / bl).clip(lower=0)

            mean_tpm = linear_vals.mean()
            mean_fc = log2fc_vals.mean()
            q1 = log2fc_vals.quantile(0.25)
            median_fc = log2fc_vals.median()
            q3 = log2fc_vals.quantile(0.75)

            rows.append({
                "Receptor": rinfo["label"],
                "Ligand": linfo["label"],
                "Mean log₂FC": mean_fc,
                "Median log₂FC": median_fc,
                "Q1": q1,
                "Q3": q3,
                "IQR_lower": max(0, mean_fc - q1),
                "IQR_upper": max(0, q3 - mean_fc),
                "Mean TPM": mean_tpm,
                "Normal TPM": bl,
                "Fold Change": mean_tpm / bl,
            })

    if not rows:
        return None, None

    df = pd.DataFrame(rows)

    # Sort receptors by total log₂FC
    receptor_totals = df.groupby("Receptor")["Mean log₂FC"].sum().sort_values(ascending=False)
    receptor_order = receptor_totals.index.tolist()
    df["_sort"] = df["Receptor"].map({r: i for i, r in enumerate(receptor_order)})
    df = df.sort_values("_sort")

    # ── Chart 1: Stacked bar (log₂FC) ────────────────────────────
    stacked_fig = px.bar(
        df, x="Receptor", y="Mean log₂FC", color="Ligand",
        barmode="stack",
        hover_data={
            "Ligand": True,
            "Mean log₂FC": ":.2f",
            "Fold Change": ":.1f",
            "Mean TPM": ":.1f",
            "Normal TPM": ":.1f",
            "_sort": False,
        },
        title=f"Ligand Upregulation vs {baseline_source} — {TCGA_PROJECTS[project_id]}",
        template="plotly",
        category_orders={"Receptor": receptor_order},
    )
    stacked_fig.update_layout(
        xaxis_title="T Cell Receptor",
        yaxis_title="log₂(Tumor / Normal)",
        xaxis_tickangle=35,
        height=450,
        margin=dict(l=60, r=20, t=60, b=100),
        legend_title="TME Ligand",
    )

    # ── Chart 2: Grouped bar with IQR (log₂FC variability) ───────
    detail_fig = go.Figure()

    ligand_names = df["Ligand"].unique().tolist()
    colors = px.colors.qualitative.Plotly
    ligand_color_map = {lig: colors[i % len(colors)] for i, lig in enumerate(ligand_names)}

    for ligand_name in ligand_names:
        sub = df[df["Ligand"] == ligand_name].copy()
        sub = sub.set_index("Receptor").reindex(receptor_order).dropna(subset=["Mean log₂FC"]).reset_index()
        if sub.empty:
            continue

        detail_fig.add_trace(go.Bar(
            name=ligand_name,
            x=sub["Receptor"],
            y=sub["Mean log₂FC"],
            marker_color=ligand_color_map[ligand_name],
            error_y=dict(
                type="data",
                symmetric=False,
                array=sub["IQR_upper"].values,
                arrayminus=sub["IQR_lower"].values,
                thickness=1.5,
                width=4,
            ),
            hovertemplate=(
                "<b>%{x}</b> ← %{fullData.name}<br>"
                "Mean log₂FC: %{y:.2f} (= %{customdata[0]:.1f}× normal)<br>"
                "Q1: %{customdata[1]:.2f} | Median: %{customdata[2]:.2f} | Q3: %{customdata[3]:.2f}<br>"
                "Tumor: %{customdata[4]:.1f} TPM | Normal: %{customdata[5]:.1f} TPM"
                "<extra></extra>"
            ),
            customdata=sub[["Fold Change", "Q1", "Median log₂FC", "Q3", "Mean TPM", "Normal TPM"]].values,
        ))

    detail_fig.update_layout(
        barmode="group",
        title=f"Ligand Upregulation Variability (Q1–Q3) — {TCGA_PROJECTS[project_id]}",
        xaxis_title="T Cell Receptor",
        yaxis_title="log₂(Tumor / Normal)",
        xaxis_tickangle=35,
        height=500,
        margin=dict(l=60, r=20, t=60, b=100),
        legend_title="TME Ligand",
        template="plotly",
    )

    return stacked_fig, detail_fig


# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("#### Filters")
    population = st.selectbox("Population (Race)", POPULATION_GROUPS, index=0)

    st.markdown("---")
    st.markdown("#### Network Parameters")
    corr_threshold = st.slider("Min |ρ| threshold", 0.1, 0.8, 0.50, 0.05,
        help="Minimum absolute Spearman ρ between ligand activation scores to draw an edge.")
    p_threshold = st.select_slider("p-value", [0.001, 0.005, 0.01, 0.05], value=0.05)
    st.markdown("---")
    db_status = f"**Local DB found** (`tcga_data.db`)" if DB_AVAILABLE else "No local DB — using GDC API"
    st.caption(
        f"{db_status}. "
        f"{len(RECEPTORS)} receptors · {len(LIGANDS)} ligands. "
        "TPM from STAR-Counts, UQ-normalized, log₂(x+1)."
    )


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

st.subheader("T Cell Exhaustion & TME Suppression Explorer", divider="gray")
st.caption(
    "Identifying which immune checkpoint receptors are co-activated by the tumor microenvironment. "
    "If the TME consistently upregulates the same groups of inhibitory ligands together, "
    "those receptor pathways fire as a unit — and blocking just one may not be enough. "
    "This tool finds those groupings across TCGA cancer types so we can design combination "
    "immunotherapies that target the right clusters together."
)

# Project + family selectors side by side
col_proj, col_fam = st.columns([1, 3])
with col_proj:
    project_id = st.selectbox(
        "Cancer Type",
        list(TCGA_PROJECTS.keys()),
        format_func=lambda x: f"{x} — {TCGA_PROJECTS[x]}",
        index=list(TCGA_PROJECTS.keys()).index("TCGA-SKCM"),
        help="TCGA cancer cohort to analyze. Each project contains RNA-seq from hundreds of tumors.",
    )
    st.caption(f"{TCGA_PROJECTS[project_id]}")

with col_fam:
    all_families = sorted(set(info["family"] for info in RECEPTORS.values()))
    selected_families = st.multiselect(
        "Receptor Families",
        options=all_families,
        default=all_families,
        help="Which receptor families to include in the network and charts. Deselect to focus on specific pathway groups.",
    )
    if not selected_families:
        active_families = set(all_families)
    else:
        active_families = set(selected_families)

    n_active_r = sum(1 for info in RECEPTORS.values() if info["family"] in active_families)
    n_active_l = sum(1 for linfo in LIGANDS.values()
                     if not linfo["receptors"] or any(
                         RECEPTORS[r]["family"] in active_families
                         for r in linfo["receptors"] if r in RECEPTORS))
    st.caption(f"{n_active_r} receptors · {n_active_l} ligands selected")

# Build filtered receptor/ligand sets
ACTIVE_RECEPTORS = {g: info for g, info in RECEPTORS.items() if info["family"] in active_families}
ACTIVE_LIGANDS = {}
for lg, linfo in LIGANDS.items():
    if not linfo["receptors"]:
        ACTIVE_LIGANDS[lg] = linfo
    elif any(r in ACTIVE_RECEPTORS for r in linfo["receptors"]):
        ACTIVE_LIGANDS[lg] = linfo

n_active_r = len(ACTIVE_RECEPTORS)
n_active_l = len(ACTIVE_LIGANDS)

# Load data — get all samples (no max cap)
progress_bar = st.progress(0, text="Loading data...")
try:
    receptor_df, ligand_df, demo_df, data_source = get_full_dataset(
        project_id, max_samples=99999, progress_bar=progress_bar
    )
    progress_bar.empty()
except GDCError as e:
    progress_bar.empty()
    st.error(
        f"**Data loading failed**\n\n{e}\n\n"
        f"**To use offline mode**, run the downloader first:\n"
        f"```\npython download_to_db.py\n```\n"
        f"This saves all TCGA projects to `tcga_data.db` for instant offline access."
    )
    st.stop()
except Exception as e:
    progress_bar.empty()
    st.error(f"**Error:** `{type(e).__name__}: {e}`"); st.stop()

# Fetch GTEx baseline (cached after first call, independent of TCGA data source)
with st.spinner("Fetching GTEx normal tissue baselines (first time only)..."):
    gtex_baseline, gtex_tissue = fetch_gtex_baseline(project_id)

# ── Stage segmented control (placed prominently at top of page) ──────────────
# Build stage options from the data
STAGE_ORDER = ["All Stages", "Stage I", "Stage II", "Stage III", "Stage IV", "Not Reported"]
if "stage" in demo_df.columns:
    available_stages = demo_df["stage"].value_counts()
    stage_options = ["All Stages"]
    for s in ["Stage I", "Stage II", "Stage III", "Stage IV", "Stage 0/IS", "Not Reported"]:
        if s in available_stages.index and available_stages[s] >= 1:
            n = available_stages[s]
            stage_options.append(s)
else:
    stage_options = ["All Stages"]

st.markdown(f"##### {TCGA_PROJECTS[project_id]} · {population}")

# Segmented control for stage
selected_stage = st.segmented_control(
    "Cancer Stage (AJCC Pathologic)",
    options=stage_options,
    default="All Stages",
    help="Filter patients by AJCC pathologic stage at diagnosis. Stage data from GDC clinical records.",
)
if selected_stage is None:
    selected_stage = "All Stages"

bl_status = f"GTEx baseline ({gtex_tissue}) ✓" if not gtex_baseline.empty else "cohort baseline"
if data_source == "LOCAL DB":
    db_size = os.path.getsize(DB_PATH) / (1024*1024)
    st.markdown(f'<span class="gdc-badge gdc-live">● OFFLINE — DB ({db_size:.0f} MB) · {bl_status} · UQ-norm</span>', unsafe_allow_html=True)
else:
    st.markdown(f'<span class="gdc-badge gdc-live">● LIVE — GDC API · {bl_status} · UQ-norm</span>', unsafe_allow_html=True)

# ── Apply filters: population + stage ────────────────────────────────────────
filter_idx = demo_df.index  # start with all

# Population filter
if population != "All" and "race_label" in demo_df.columns:
    filter_idx = filter_idx.intersection(demo_df[demo_df["race_label"] == population].index)

# Stage filter
if selected_stage != "All Stages" and "stage" in demo_df.columns:
    filter_idx = filter_idx.intersection(demo_df[demo_df["stage"] == selected_stage].index)

# Apply to expression data
r_filt = receptor_df.loc[receptor_df.index.intersection(filter_idx)]
l_filt = ligand_df.loc[ligand_df.index.intersection(filter_idx)]

if len(r_filt) < 10:
    n_found = len(r_filt)
    filter_desc = []
    if population != "All": filter_desc.append(f"race='{population}'")
    if selected_stage != "All Stages": filter_desc.append(f"stage='{selected_stage}'")
    st.warning(f"Only {n_found} samples for {' + '.join(filter_desc)} — showing all samples.")
    r_filt, l_filt = receptor_df, ligand_df

# Show active filter summary
filter_parts = []
if population != "All": filter_parts.append(f"**{population}**")
if selected_stage != "All Stages": filter_parts.append(f"**{selected_stage}**")
if n_active_r < len(RECEPTORS): filter_parts.append(f"**{n_active_r} receptors**")
if filter_parts:
    st.caption(f"Filtered to: {' · '.join(filter_parts)} — {len(r_filt)} samples")

# ── Apply family filter to columns ───────────────────────────────────────
active_rcols = [c for c in r_filt.columns if c in ACTIVE_RECEPTORS]
active_lcols = [c for c in l_filt.columns if c in ACTIVE_LIGANDS]
r_filt = r_filt[active_rcols]
l_filt = l_filt[active_lcols]

edge_df = compute_coexpression(l_filt, corr_threshold, p_threshold)
shared_pairs = find_shared_ligand_pairs(l_filt)
G = build_graph(edge_df, shared_ligand_pairs=shared_pairs, active_receptors=ACTIVE_RECEPTORS)
pos = hierarchical_layout(G)

# Compute ligand activation scores for metrics
activation_scores = compute_ligand_activation_scores(l_filt, active_receptors=ACTIVE_RECEPTORS)

c1,c2,c3,c4 = st.columns(4)
c1.metric("Samples", len(r_filt))
c2.metric("Receptors", n_active_r)
c3.metric("Ligands", n_active_l)
avg_rho = f"{edge_df['rho'].abs().mean():.3f}" if len(edge_df) > 0 else "—"
c4.metric("Mean |ρ| (ligand)", avg_rho)

stage_label = selected_stage if selected_stage != "All Stages" else "All Stages"
net_fig = create_network(G, pos,
    title=f"TME Co-Activation — {project_id} · {population} · {stage_label}",
    activation_scores=activation_scores,
    gtex_baseline=gtex_baseline,
    ligand_df=l_filt)
st.plotly_chart(net_fig, use_container_width=True, config={"scrollZoom": True})
st.caption(
    "**Node size = potential activation (ligand log₂FC vs normal).** "
    "**Blue** = co-activation (+ρ). **Red** = inverse correlation (−ρ). "
    "**Dashed** = also shares a ligand. **Black dotted** = identical ligand set. "
    "**Dashed orange** = shared ligand, no significant ρ."
)

# ── Tabs ─────────────────────────────────────────────────────────────────────
tab0, tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📊 Pathway Activation", "🎯 Ligand Breakdown", "🔴 TME Ligand Activity",
    "🔥 Correlation Matrix", "📋 Edge Table", "👥 Demographics", "🧪 Receptor–Ligand Pairs"
])

with tab0:
    st.markdown("### Per-Receptor Ligand Activation Distribution")
    st.markdown(
        "Each row is a receptor. Overlaid violins show the **per-patient distribution of each ligand's "
        "log₂FC** vs GTEx normal tissue. Different colors = different ligands. "
        "**Further right = more upregulated in tumor → stronger activation of that receptor pathway.** "
        "Near zero = similar to healthy tissue. Left of the dashed line = below normal. "
        "Where violins overlap, those ligands co-vary across patients. "
        "Wider = more patients at that level. Hover for ligand identity."
    )
    st.plotly_chart(create_barplot(r_filt, project_id, gtex_baseline=gtex_baseline, ligand_df=l_filt), use_container_width=True)

with tab1:
    st.markdown("### Ligand Breakdown by Receptor")
    baseline_desc = "**GTEx healthy tissue**" if (not gtex_baseline.empty) else "**cohort median**"
    st.markdown(
        f"Stacked chart shows **log₂ fold-change** of each ligand over "
        f"{baseline_desc} baseline. A value of 1 = 2× normal, 3 = 8× normal. "
        "Receptors sorted by potential activation."
    )
    stacked_chart, detail_chart = create_receptor_activation_chart(l_filt, project_id, gtex_baseline=gtex_baseline, active_receptors=ACTIVE_RECEPTORS)
    if stacked_chart:
        st.plotly_chart(stacked_chart, use_container_width=True)
    else:
        st.info("No receptor-ligand activation data available.")

with tab2:
    st.markdown("### TME Suppressive Ligand Landscape")
    st.markdown(
        "These are the ligands and enzymes expressed **by tumor/stromal/myeloid cells** in the TME "
        "that engage the T cell inhibitory receptors. Higher expression = more suppressive microenvironment."
    )
    st.plotly_chart(create_ligand_barplot(l_filt, gtex_baseline=gtex_baseline), use_container_width=True)

    # Heatmap
    tme_df = compute_tme_suppression(r_filt, l_filt, gtex_baseline=gtex_baseline)
    if not tme_df.empty:
        hm = create_tme_heatmap(tme_df)
        if hm: st.plotly_chart(hm, use_container_width=True)

with tab3:
    st.plotly_chart(create_corrmatrix(l_filt), use_container_width=True)

with tab4:
    if len(edge_df) > 0:
        d = edge_df.copy()
        d["Receptor A"] = d["source"].map(lambda g: RECEPTORS[g]["label"])
        d["Receptor B"] = d["target"].map(lambda g: RECEPTORS[g]["label"])
        d = d[["Receptor A","ligands_A","mean_score_A","Receptor B","ligands_B","mean_score_B","weight","rho","pval"]]
        d.columns = ["Receptor A","A's Ligands","A Mean Lig. Score",
                      "Receptor B","B's Ligands","B Mean Lig. Score",
                      "Weight (ρ)","Co-activation ρ","p-value"]
        st.dataframe(d.sort_values("Weight (ρ)", ascending=False), use_container_width=True, height=500)
    else:
        st.info("No edges at current thresholds.")

with tab5:
    if not demo_df.empty:
        col_a, col_b = st.columns(2)
        with col_a:
            if "race_label" in demo_df.columns:
                rc = demo_df["race_label"].value_counts().reset_index(); rc.columns=["Race","Count"]
                st.plotly_chart(px.bar(rc,x="Race",y="Count",title=f"Samples by Race — {project_id}",
                    template="plotly").update_layout(height=400,xaxis_tickangle=30), use_container_width=True)
        with col_b:
            if "stage" in demo_df.columns:
                sc = demo_df["stage"].value_counts().reset_index(); sc.columns=["Stage","Count"]
                # Sort by stage order
                stage_order = ["Stage 0/IS","Stage I","Stage II","Stage III","Stage IV","Not Reported"]
                sc["_sort"] = sc["Stage"].apply(lambda x: stage_order.index(x) if x in stage_order else 99)
                sc = sc.sort_values("_sort").drop(columns="_sort")
                st.plotly_chart(px.bar(sc,x="Stage",y="Count",title=f"Samples by Stage — {project_id}",
                    template="plotly", color="Stage").update_layout(height=400,showlegend=False),
                    use_container_width=True)
    else:
        st.info("No demographics.")

with tab6:
    st.markdown("### Receptor–Ligand Pair Analysis")
    st.markdown(
        "For each inhibitory receptor on T cells, we identify its known TME ligand(s) and "
        "compute: mean expression of both, a combined suppressive score, and the Spearman "
        "correlation between receptor and ligand expression across patient samples."
    )
    tme_df2 = compute_tme_suppression(r_filt, l_filt, gtex_baseline=gtex_baseline)
    if not tme_df2.empty:
        display_cols = ["Receptor","Ligand","Receptor log₂FC","Ligand log₂FC",
                        "Suppressive Score","R-L Correlation (ρ)","R-L p-value"]
        st.dataframe(
            tme_df2[display_cols].sort_values("Suppressive Score", ascending=False),
            use_container_width=True, height=500,
        )
    else:
        st.info("No receptor-ligand pairs computed.")

# ── Population comparison ────────────────────────────────────────────────────
st.markdown("---")
st.markdown("## Population-Level Variability")
st.caption(
    "Compares mean TME ligand upregulation (log₂FC vs GTEx normal) across racial/ethnic groups. "
    "Differences may reflect variation in tumor immunogenicity, stromal composition, "
    "or cohort demographics — these are exploratory observations, not controlled comparisons."
)
if not demo_df.empty and "race_label" in demo_df.columns:
    gcols = [c for c in l_filt.columns if c in LIGANDS]
    pm_rows = []
    for rl in demo_df["race_label"].unique():
        rc2 = demo_df[demo_df["race_label"]==rl].index
        le = ligand_df.loc[ligand_df.index.intersection(rc2), gcols]
        if len(le) < 5:
            continue
        # For each receptor, sum its ligands' log₂FC
        for rg, rinfo in ACTIVE_RECEPTORS.items():
            lig_genes = [lg for lg, linfo in LIGANDS.items()
                         if rg in linfo.get("receptors", []) and lg in le.columns]
            if not lig_genes:
                continue
            total_fc = 0
            for lg in lig_genes:
                linear = np.power(2, le[lg]) - 1
                bl = gtex_baseline[lg] if (gtex_baseline is not None and not gtex_baseline.empty and lg in gtex_baseline.index) else max(linear.median(), 0.1)
                bl = max(bl, 0.1)
                total_fc += max(0, np.log2((linear.mean() + 0.1) / bl))
            pm_rows.append({
                "Population": f"{rl} (n={len(le)})",
                "Receptor": rinfo["label"],
                "log₂(Tumor / Normal)": total_fc,
            })
    if pm_rows:
        pm_df = pd.DataFrame(pm_rows)
        pf = px.bar(pm_df, x="Receptor", y="log₂(Tumor / Normal)", color="Population",
                    barmode="group",
                    title=f"Ligand Activation per Receptor by Race — {project_id}",
                    template="plotly")
        pf.update_layout(height=550, xaxis_tickangle=45)
        st.plotly_chart(pf, use_container_width=True)

st.markdown("---")
st.caption(f"TCGA Explorer · {len(RECEPTORS)} receptors + {len(LIGANDS)} ligands · UQ-normalized · GDC API")