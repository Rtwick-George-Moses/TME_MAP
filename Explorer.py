"""
TCGA T Cell Exhaustion & TME Suppression Explorer
===================================================
Receptor-normalized analysis using LRAR (Ligand log₂FC − Receptor log₂FC).
This cancels pipeline batch effects between TCGA (STAR) and GTEx (RSEM).

Reads from local SQLite database (tcga_data.db) if available,
otherwise falls back to live GDC API calls.
Run download_to_db.py first for offline use.
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
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
.main .block-container { padding-top: 2rem; max-width: 1400px; }
</style>
""", unsafe_allow_html=True)

BG = dict(plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)")

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


def download_expression_batched(file_ids, max_files=0, batch_size=25, progress_bar=None):
    """
    Download expression data in batches with progress updates.
    Each batch is independently cached via _download_batch.
    """
    fids = file_ids if max_files <= 0 else file_ids[:max_files]
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


def get_full_dataset(project_id, max_samples=0, progress_bar=None):
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
            _update(80, "Splitting receptor and ligand matrices...")

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
    n_files = len(files) if max_samples <= 0 else min(len(files), max_samples)
    _update(20, f"Found {len(files)} files. Preparing to download {n_files} profiles...")

    fids = files["file_id"].tolist() if max_samples <= 0 else files["file_id"].tolist()[:max_samples]
    expr = download_expression_batched(fids, max_files=0, batch_size=25, progress_bar=progress_bar)
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
    _update(92, "Splitting receptor and ligand matrices...")

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
    Get per-ligand tissue baseline in linear TPM from GTEx.
    Requires GTEx baseline — no fallback to cohort median.
    """
    baselines = {}
    for lg in ligand_df.columns:
        if lg in LIGANDS:
            if gtex_baseline is not None and not gtex_baseline.empty and lg in gtex_baseline.index:
                baselines[lg] = max(0, gtex_baseline[lg])
            else:
                baselines[lg] = 0.1  # floor for genes missing from GTEx
    return pd.Series(baselines)

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

# ══════════════════════════════════════════════════════════════════════════════
# CORE: Receptor-Normalized Activation
# ══════════════════════════════════════════════════════════════════════════════

def compute_receptor_log2fc(receptor_df, gtex_baseline):
    """Compute log₂FC for each receptor gene vs GTEx."""
    fc = {}
    for gene in receptor_df.columns:
        if gene not in RECEPTORS:
            continue
        linear = np.power(2, receptor_df[gene]) - 1
        bl = gtex_baseline[gene] if (gtex_baseline is not None and not gtex_baseline.empty and gene in gtex_baseline.index) else 0.1
        bl = max(bl, 0.1)
        fc[gene] = np.log2((linear.mean() + 0.1) / bl)
    return fc


def compute_normalized_scores(receptor_df, ligand_df, gtex_baseline):
    """
    For each receptor-ligand pair, compute:
      - Ligand log₂FC vs GTEx
      - Receptor log₂FC vs GTEx
      - LRAR = Ligand log₂FC − Receptor log₂FC

    Each ligand independently compared against the receptor's full upregulation.
    Positive LRAR = ligand outpaces receptor (TME actively suppressing)
    Zero = balanced
    Negative = receptor outpaces ligand (pathway is ligand-limited)
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    r_fc = compute_receptor_log2fc(receptor_df, gtex_baseline)

    rows = []
    for gene, rinfo in RECEPTORS.items():
        if gene not in receptor_df.columns or gene not in r_fc:
            continue
        receptor_fc = r_fc[gene]

        for lgene, linfo in LIGANDS.items():
            if gene not in linfo.get("receptors", []) or lgene not in ligand_df.columns:
                continue

            l_linear = np.power(2, ligand_df[lgene]) - 1
            l_bl = max(baseline.get(lgene, 0), 0.1)
            ligand_fc = np.log2((l_linear.mean() + 0.1) / l_bl)

            lrar = ligand_fc - receptor_fc

            rows.append({
                "Receptor": rinfo["label"],
                "Receptor Gene": gene,
                "Ligand": linfo["label"],
                "Ligand Gene": lgene,
                "Receptor log₂FC": receptor_fc,
                "Ligand log₂FC": ligand_fc,
                "LRAR": lrar,
                "Family": rinfo["family"],
            })

    return pd.DataFrame(rows)


def compute_per_patient_lrar(receptor_df, ligand_df, gtex_baseline):
    """
    Per-patient LRAR for each receptor-ligand pair.
    LRAR = Ligand_log₂FC − Receptor_log₂FC (per patient).
    Each ligand is independently compared against the receptor's full upregulation.
    Returns dict: (receptor_gene, ligand_gene) → Series of per-patient LRAR values.
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    results = {}

    for gene, rinfo in RECEPTORS.items():
        if gene not in receptor_df.columns:
            continue
        r_linear = np.power(2, receptor_df[gene]) - 1
        r_bl = gtex_baseline[gene] if (gtex_baseline is not None and not gtex_baseline.empty and gene in gtex_baseline.index) else 0.1
        r_bl = max(r_bl, 0.1)
        r_fc_per_patient = np.log2((r_linear + 0.1) / r_bl)

        for lgene, linfo in LIGANDS.items():
            if gene not in linfo.get("receptors", []) or lgene not in ligand_df.columns:
                continue
            l_linear = np.power(2, ligand_df[lgene]) - 1
            l_bl = max(baseline.get(lgene, 0), 0.1)
            l_fc_per_patient = np.log2((l_linear + 0.1) / l_bl)

            lrar = l_fc_per_patient - r_fc_per_patient
            results[(gene, lgene)] = lrar

    return results


def compute_normalized_activation_scores(receptor_df, ligand_df, gtex_baseline, active_receptors=None):
    """
    Per-patient normalized activation score per receptor.
    Sum of (ligand_FC − receptor_FC) across all ligands for that receptor.
    """
    baseline = compute_tissue_baseline(ligand_df, gtex_baseline=gtex_baseline)
    receptors_to_use = active_receptors if active_receptors is not None else RECEPTORS
    scores = {}

    for gene, rinfo in receptors_to_use.items():
        if gene not in receptor_df.columns:
            continue
        r_linear = np.power(2, receptor_df[gene]) - 1
        r_bl = gtex_baseline[gene] if (gtex_baseline is not None and not gtex_baseline.empty and gene in gtex_baseline.index) else 0.1
        r_bl = max(r_bl, 0.1)
        r_fc = np.log2((r_linear + 0.1) / r_bl)

        lig_genes = [lg for lg, li in LIGANDS.items()
                     if gene in li.get("receptors", []) and lg in ligand_df.columns]
        if not lig_genes:
            continue

        total_lrar = pd.Series(0.0, index=ligand_df.index)
        for lg in lig_genes:
            l_linear = np.power(2, ligand_df[lg]) - 1
            l_bl = max(baseline.get(lg, 0), 0.1)
            l_fc = np.log2((l_linear + 0.1) / l_bl)
            total_lrar += l_fc - r_fc

        scores[gene] = total_lrar

    return pd.DataFrame(scores)


def compute_normalized_coexpression(receptor_df, ligand_df, gtex_baseline, thresh=0.3, p_thresh=0.05):
    """Spearman correlation of normalized activation scores."""
    activation = compute_normalized_activation_scores(receptor_df, ligand_df, gtex_baseline)
    receptors_with_scores = activation.columns.tolist()

    receptor_ligand_genes = {}
    receptor_ligand_labels = {}
    for rg in receptors_with_scores:
        genes = [lg for lg, li in LIGANDS.items()
                 if rg in li.get("receptors", []) and lg in ligand_df.columns]
        receptor_ligand_genes[rg] = set(genes)
        receptor_ligand_labels[rg] = ", ".join(LIGANDS[lg]["label"] for lg in genes) or "none"

    edges = []
    for g1, g2 in combinations(receptors_with_scores, 2):
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
                "ligands_A": receptor_ligand_labels[g1],
                "ligands_B": receptor_ligand_labels[g2],
                "mean_score_A": activation[g1].mean(),
                "mean_score_B": activation[g2].mean(),
            })
    return pd.DataFrame(edges)


# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION (all charts use LRAR instead of raw log₂FC)
# ══════════════════════════════════════════════════════════════════════════════

def create_normalized_ridgeline(receptor_df, ligand_df, gtex_baseline, project_id, df_filt):
    """Ridgeline density of per-patient LRAR per receptor per ligand."""
    patient_lrar = compute_per_patient_lrar(receptor_df, ligand_df, gtex_baseline)

    rows = []
    receptor_medians = {}
    for (rg, lg), vals in patient_lrar.items():
        rinfo = RECEPTORS[rg]
        linfo = LIGANDS[lg]
        vals_filt = vals.loc[vals.index.isin(df_filt.index)].dropna()
        for v in vals_filt.values:
            rows.append({"Receptor": rinfo["label"], "Ligand": linfo["label"], "LRAR": v})

        # Track median for sorting
        receptor_medians.setdefault(rinfo["label"], []).append(vals_filt.median())

    if not rows:
        fig = go.Figure()
        fig.update_layout(title="No data available")
        return fig

    # Average median across ligands for sorting
    for k in receptor_medians:
        receptor_medians[k] = np.mean(receptor_medians[k])

    rdf = pd.DataFrame(rows)
    receptor_order = sorted(receptor_medians.keys(), key=lambda r: receptor_medians[r], reverse=True)

    all_ligands = rdf["Ligand"].unique().tolist()
    colors_list = px.colors.qualitative.Plotly + px.colors.qualitative.D3
    ligand_colors = {lig: colors_list[i % len(colors_list)] for i, lig in enumerate(all_ligands)}

    n_receptors = len(receptor_order)
    fig = make_subplots(rows=n_receptors, cols=1, shared_xaxes=True,
                        vertical_spacing=0.0, subplot_titles=None)

    legend_shown = set()

    for ri, receptor_name in enumerate(receptor_order):
        row_idx = ri + 1
        sub = rdf[rdf["Receptor"] == receptor_name]

        # Histogram — per-patient combined LRAR: log₂(ΣL_tumor / ΣL_normal) − log₂(R_tumor / R_normal)
        receptor_gene = next((g for g, info in RECEPTORS.items() if info["label"] == receptor_name), None)
        if receptor_gene and receptor_gene in receptor_df.columns:
            lig_genes_for_hist = [lg for lg, li in LIGANDS.items()
                                  if receptor_gene in li.get("receptors", []) and lg in ligand_df.columns]
            if lig_genes_for_hist:
                # Sum ligands in linear TPM per patient
                tumor_sum = pd.Series(0.0, index=ligand_df.index)
                normal_sum = 0.0
                for lg in lig_genes_for_hist:
                    tumor_sum += np.power(2, ligand_df[lg]) - 1
                    bl = gtex_baseline[lg] if (gtex_baseline is not None and not gtex_baseline.empty and lg in gtex_baseline.index) else 0.1
                    normal_sum += max(bl, 0.1)

                ligand_combined_fc = np.log2((tumor_sum + 0.1) / max(normal_sum, 0.1))

                # Receptor FC per patient
                r_linear = np.power(2, receptor_df[receptor_gene]) - 1
                r_bl = gtex_baseline[receptor_gene] if (gtex_baseline is not None and not gtex_baseline.empty and receptor_gene in gtex_baseline.index) else 0.1
                r_bl = max(r_bl, 0.1)
                receptor_fc = np.log2((r_linear + 0.1) / r_bl)

                # Combined LRAR per patient
                combined = (ligand_combined_fc - receptor_fc).loc[ligand_combined_fc.index.isin(df_filt.index)].dropna().values

                if len(combined) > 0:
                    counts, bin_edges = np.histogram(combined, bins=25)
                    counts_norm = counts / counts.max() * 0.6 if counts.max() > 0 else counts
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    bin_width = bin_edges[1] - bin_edges[0]
                    fig.add_trace(go.Bar(
                        x=bin_centers, y=counts_norm, width=bin_width * 0.9,
                        marker=dict(color="rgba(180,180,180,0.35)", line=dict(width=0.5, color="rgba(150,150,150,0.3)")),
                        showlegend=(ri == 0), name="Patient distribution", legendgroup="hist",
                        customdata=counts,
                        hovertemplate=f"<b>{receptor_name}</b> (combined)<br>LRAR: %{{x:.2f}}<br>Patients: %{{customdata}}<extra></extra>",
                    ), row=row_idx, col=1)

        # KDE per ligand
        for lig_name in sub["Ligand"].unique():
            lig_sub = sub[sub["Ligand"] == lig_name]
            show_legend = lig_name not in legend_shown
            legend_shown.add(lig_name)

            vals = lig_sub["LRAR"].dropna().values
            if len(vals) < 3:
                continue
            try:
                kde = gaussian_kde(vals, bw_method=0.3)
            except Exception:
                continue

            x_grid = np.linspace(vals.min() - 0.5, vals.max() + 0.5, 200)
            y_kde = kde(x_grid)
            y_kde = y_kde / y_kde.max() if y_kde.max() > 0 else y_kde

            color = ligand_colors[lig_name]
            import plotly.colors as pc
            rgb = pc.hex_to_rgb(color) if color.startswith("#") else (99, 110, 250)

            fig.add_trace(go.Scatter(
                x=x_grid, y=y_kde, mode="lines", fill="tozeroy",
                fillcolor=f"rgba({rgb[0]},{rgb[1]},{rgb[2]},0.35)",
                line=dict(color=f"rgba({rgb[0]},{rgb[1]},{rgb[2]},0.85)", width=1.5),
                name=lig_name, legendgroup=lig_name, showlegend=show_legend,
                hovertemplate=f"<b>{receptor_name}</b> ← {lig_name}<br>LRAR: %{{x:.2f}}<extra></extra>",
            ), row=row_idx, col=1)

        fig.update_yaxes(title_text=receptor_name, title_font=dict(size=10),
                         title_standoff=5, showticklabels=False, showgrid=False,
                         zeroline=False, row=row_idx, col=1)
        fig.add_vline(x=0, line_dash="dot", line_color="gray", opacity=0.5, row=row_idx, col=1)

    fig.update_xaxes(title_text="LRAR (Ligand log₂FC − Receptor log₂FC)", row=n_receptors, col=1)

    fig.update_layout(
        title=f"Normalized Ligand Activation Ridgeline — {TCGA_PROJECTS[project_id]}",
        height=max(600, n_receptors * 70), margin=dict(l=120, r=20, t=50, b=50),
        template="plotly", legend_title="TME Ligand", showlegend=True, barmode="overlay", **BG,
    )
    return fig


def create_normalized_breakdown(norm_df, project_id):
    """Horizontal bar chart of LRAR per receptor×ligand."""
    all_ligands = norm_df["Ligand"].unique().tolist()
    colors_list = px.colors.qualitative.Plotly + px.colors.qualitative.D3
    lig_color_map = {lig: colors_list[i % len(colors_list)] for i, lig in enumerate(all_ligands)}

    # Sort receptors by mean LRAR
    receptor_mean = norm_df.groupby("Receptor")["LRAR"].mean().sort_values(ascending=False)
    receptor_order = receptor_mean.index.tolist()

    y_labels, x_vals, colors, customdata = [], [], [], []
    for rname in receptor_order:
        rsub = norm_df[norm_df["Receptor"] == rname].sort_values("LRAR", ascending=False)
        for _, row in rsub.iterrows():
            y_labels.append(f"{rname}  ·  {row['Ligand']}")
            x_vals.append(row["LRAR"])
            colors.append(lig_color_map.get(row["Ligand"], "#888"))
            customdata.append([row["Ligand log₂FC"], row["Receptor log₂FC"]])

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=y_labels, x=x_vals, orientation="h", marker=dict(color=colors),
        customdata=customdata, showlegend=False,
        hovertemplate=(
            "<b>%{y}</b><br>"
            "LRAR: %{x:.2f}<br>"
            "Ligand log₂FC: %{customdata[0]:.2f} | Receptor log₂FC: %{customdata[1]:.2f}"
            "<extra></extra>"
        ),
    ))

    # Legend
    for lig in all_ligands:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=lig_color_map[lig]), name=lig, showlegend=True))

    fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5,
                  annotation_text="Balanced (0)", annotation_position="top")

    fig.update_layout(
        title=f"LRAR per Receptor–Ligand Pair — {TCGA_PROJECTS[project_id]}",
        xaxis_title="LRAR (Ligand log₂FC − Receptor log₂FC)",
        yaxis=dict(categoryorder="array", categoryarray=list(reversed(y_labels)), tickfont=dict(size=9)),
        height=max(500, len(y_labels) * 24), margin=dict(l=200, r=20, t=60, b=50),
        legend_title="TME Ligand", template="plotly", **BG,
    )
    return fig


def create_normalized_corrmatrix(receptor_df, ligand_df, gtex_baseline):
    """Correlation matrix of normalized activation scores."""
    activation = compute_normalized_activation_scores(receptor_df, ligand_df, gtex_baseline)
    if activation.empty:
        return go.Figure()
    lb = [RECEPTORS[g]["label"] for g in activation.columns]
    corr = activation.corr(method="spearman")
    fig = go.Figure(go.Heatmap(
        z=corr.values, x=lb, y=lb, colorscale="RdBu_r", zmid=0, zmin=-1, zmax=1,
        hovertemplate="<b>%{x} × %{y}</b><br>ρ=%{z:.3f}<extra></extra>",
        colorbar=dict(title="ρ"),
    ))
    fig.update_layout(
        title="Normalized Co-Activation Correlation Matrix",
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(autorange="reversed", tickfont=dict(size=9)),
        height=700, margin=dict(l=140, r=20, t=50, b=140), template="plotly", **BG,
    )
    return fig


def create_normalized_network(edge_df, norm_df, ligand_df, active_receptors=None):
    """Network graph with node size = mean LRAR, edges = normalized correlation."""

    receptors_to_use = active_receptors if active_receptors is not None else RECEPTORS
    G = nx.Graph()
    for g, info in receptors_to_use.items():
        G.add_node(g, **info)

    # Node sizing from mean LRAR
    receptor_lrar = {}
    for rg, rinfo in receptors_to_use.items():
        rsub = norm_df[norm_df["Receptor Gene"] == rg]
        if not rsub.empty:
            receptor_lrar[rg] = rsub["LRAR"].mean()

    lv = list(receptor_lrar.values()) if receptor_lrar else [0]
    lmin, lmax = min(lv), max(lv)
    lrng = lmax - lmin if lmax > lmin else 1
    MIN_NODE, MAX_NODE = 15, 55

    def lrar_to_size(gene):
        v = receptor_lrar.get(gene, 0)
        return MIN_NODE + ((v - lmin) / lrng) * (MAX_NODE - MIN_NODE)

    # Build shared ligand map
    shared_pairs = find_shared_ligand_pairs(ligand_df)
    shared_map = {}
    if shared_pairs:
        for sp in shared_pairs:
            key = tuple(sorted([sp["source"], sp["target"]]))
            shared_map[key] = sp

    # Add correlation edges, tagging shared ligands
    for _, r in edge_df.iterrows():
        s, t = r["source"], r["target"]
        if s not in receptors_to_use or t not in receptors_to_use:
            continue
        key = tuple(sorted([s, t]))
        sp = shared_map.pop(key, None)

        if sp:
            edge_type = "both"
            shared_ligands = sp["shared_ligands"]
        else:
            edge_type = "correlation"
            shared_ligands = ""

        G.add_edge(s, t, weight=r["weight"], rho=r["rho"], pval=r["pval"],
                   ligands_A=r.get("ligands_A", ""), ligands_B=r.get("ligands_B", ""),
                   mean_score_A=r.get("mean_score_A", 0), mean_score_B=r.get("mean_score_B", 0),
                   edge_type=edge_type, shared_ligands=shared_ligands)

    # Remaining shared_map entries: no correlation edge
    for key, sp in shared_map.items():
        s, t = sp["source"], sp["target"]
        if s not in receptors_to_use or t not in receptors_to_use:
            continue
        etype = "shared_identical" if sp["is_identical"] else "shared_partial"
        G.add_edge(s, t, weight=0, rho=0, pval=1,
                   ligands_A="", ligands_B="",
                   mean_score_A=0, mean_score_B=0,
                   edge_type=etype, shared_ligands=sp["shared_ligands"])

    pos = hierarchical_layout(G)
    fig = go.Figure()

    # Categorize edges
    all_edges = list(G.edges(data=True))
    corr_edges = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "correlation"]
    both_edges = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "both"]
    shared_identical = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "shared_identical"]
    shared_partial = [(u,v,d) for u,v,d in all_edges if d.get("edge_type") == "shared_partial"]

    weighted = corr_edges + both_edges
    MIN_W, MAX_W = 1.0, 14.0
    aw = [d["weight"] for _,_,d in weighted]
    wn = min(aw) if aw else 0; wx = max(aw) if aw else 1; wr = wx - wn if wx > wn else 1

    def _hover_corr(lu, lv_l, d, u_gene, v_gene):
        # Use the same mean LRAR as the nodes
        u_lrar = receptor_lrar.get(u_gene, 0)
        v_lrar = receptor_lrar.get(v_gene, 0)
        h = (f"<b>{lu} ↔ {lv_l}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Co-activation ρ:</b> {d['rho']:.4f}<br>"
             f"<b>Weight:</b> {d['weight']:.4f}<br>"
             f"<b>p-value:</b> {d['pval']:.2e}<br>"
             f"<b>{lu} ligands:</b> {d.get('ligands_A','')}<br>"
             f"<b>{lu} mean LRAR:</b> {u_lrar:.2f}<br>"
             f"<b>{lv_l} ligands:</b> {d.get('ligands_B','')}<br>"
             f"<b>{lv_l} mean LRAR:</b> {v_lrar:.2f}")
        shared = d.get("shared_ligands", "")
        if shared:
            h += f"<br><b>⚠ Also shares ligand:</b> {shared}"
        return h

    def _add_edge(fig, u, v, d, line_dict, hover_text):
        x0,y0=pos[u]; x1,y1=pos[v]
        fig.add_trace(go.Scatter(x=[x0,x1,None],y=[y0,y1,None],mode="lines",
            line=line_dict,hoverinfo="text",hovertext=[hover_text,hover_text,None],showlegend=False))
        fig.add_trace(go.Scatter(x=[(x0+x1)/2],y=[(y0+y1)/2],mode="markers",
            marker=dict(size=max(24, line_dict.get("width",3)*2.5),color="rgba(0,0,0,0)"),
            hoverinfo="text",hovertext=hover_text,showlegend=False))

    # 1. Shared-identical: black dotted
    for u,v,d in shared_identical:
        lu=RECEPTORS[u]["label"]; lv_l=RECEPTORS[v]["label"]
        shared = d.get("shared_ligands","?")
        h = (f"<b>{lu} ↔ {lv_l}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Identical ligand set:</b> {shared}<br>"
             f"Both receptors bind exactly the same TME ligand(s).<br>"
             f"Co-activation is guaranteed (ρ = 1.0 by definition).")
        _add_edge(fig, u, v, d, line_dict=dict(width=3, color="rgba(0,0,0,0.45)", dash="dot"), hover_text=h)

    # 2. Shared-partial: dashed orange
    for u,v,d in shared_partial:
        lu=RECEPTORS[u]["label"]; lv_l=RECEPTORS[v]["label"]
        shared = d.get("shared_ligands","?")
        h = (f"<b>{lu} ↔ {lv_l}</b><br>━━━━━━━━━━━━━━━<br>"
             f"<b>Shared ligand:</b> {shared}<br>"
             f"These receptors share some (not all) TME ligands.<br>"
             f"No significant independent co-activation at current threshold.")
        _add_edge(fig, u, v, d, line_dict=dict(width=3, color="rgba(255,165,0,0.5)", dash="dash"), hover_text=h)

    # 3. Correlation-only: solid blue/red
    for u,v,d in corr_edges:
        lu=RECEPTORS[u]["label"]; lv_l=RECEPTORS[v]["label"]
        w=d["weight"]; rho=d["rho"]; t=(w-wn)/wr; lp=MIN_W+(t**1.5)*(MAX_W-MIN_W); op=0.35+0.55*t
        ec = f"rgba(99,110,250,{op:.2f})" if rho >= 0 else f"rgba(239,85,59,{op:.2f})"
        _add_edge(fig, u, v, d, line_dict=dict(width=lp, color=ec), hover_text=_hover_corr(lu, lv_l, d, u, v))

    # 4. Both (correlation + shared): dashed blue/red
    for u,v,d in both_edges:
        lu=RECEPTORS[u]["label"]; lv_l=RECEPTORS[v]["label"]
        w=d["weight"]; rho=d["rho"]; t=(w-wn)/wr; lp=MIN_W+(t**1.5)*(MAX_W-MIN_W); op=0.35+0.55*t
        ec = f"rgba(99,110,250,{op:.2f})" if rho >= 0 else f"rgba(239,85,59,{op:.2f})"
        _add_edge(fig, u, v, d, line_dict=dict(width=lp, color=ec, dash="dash"), hover_text=_hover_corr(lu, lv_l, d, u, v))

    # Legend entries for edge types
    has_pos = any(d["rho"] >= 0 for _,_,d in corr_edges + both_edges) if weighted else False
    has_neg = any(d["rho"] < 0 for _,_,d in corr_edges + both_edges) if weighted else False
    if has_pos:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(99,110,250,0.7)"),name="Co-activation (+ρ)",showlegend=True))
    if has_neg:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(239,85,59,0.7)"),name="Inverse correlation (−ρ)",showlegend=True))
    if both_edges:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(99,110,250,0.7)",dash="dash"),name="+ shared ligand (dashed)",showlegend=True))
    if shared_identical:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(0,0,0,0.45)",dash="dot"),name="Identical ligand (ρ=1)",showlegend=True))
    if shared_partial:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode="lines",
            line=dict(width=3,color="rgba(255,165,0,0.5)",dash="dash"),name="Shared ligand (no sig. ρ)",showlegend=True))

    # Nodes
    bf = {}
    for n, d in G.nodes(data=True):
        f = d.get("family", "Other")
        bf.setdefault(f, {"x": [], "y": [], "text": [], "hover": [], "sizes": []})
        x, y = pos[n]; lb = d.get("label", n); ds = d.get("desc", "")
        lrar_val = receptor_lrar.get(n, 0)
        sz = lrar_to_size(n)
        deg = G.degree(n)
        avg_w = np.mean([d2["weight"] for _, _, d2 in G.edges(n, data=True)]) if deg > 0 else 0

        # List ligands and their LRAR values
        lig_lines = []
        for lg, linfo in LIGANDS.items():
            if n in linfo.get("receptors", []):
                pair_row = norm_df[(norm_df["Receptor Gene"] == n) & (norm_df["Ligand Gene"] == lg)]
                if not pair_row.empty:
                    pair_lrar = pair_row["LRAR"].iloc[0]
                    lig_lines.append(f"  {linfo['label']}: LRAR {pair_lrar:+.2f}")
                else:
                    lig_lines.append(f"  {linfo['label']}")
        lig_text = "<br>".join(lig_lines) if lig_lines else "none"

        bf[f]["x"].append(x); bf[f]["y"].append(y); bf[f]["text"].append(lb)
        bf[f]["sizes"].append(sz)
        bf[f]["hover"].append(
            f"<b>{lb}</b> ({n})<br>{ds}<br>Family: {f}<br>"
            f"━━━━━━━━━━━━━━━<br>"
            f"<b>Mean LRAR:</b> {lrar_val:.2f}<br>"
            f"<b>Ligands:</b><br>{lig_text}<br>"
            f"━━━━━━━━━━━━━━━<br>"
            f"<b>Connections:</b> {deg}<br>"
            f"<b>Avg edge weight:</b> {avg_w:.4f}"
        )

    for f, v in bf.items():
        fig.add_trace(go.Scatter(x=v["x"], y=v["y"], mode="markers+text",
            marker=dict(size=v["sizes"], color=FAMILY_COLORS.get(f), line=dict(width=2, color="white")),
            text=v["text"], textposition="top center", textfont=dict(size=11),
            hoverinfo="text", hovertext=v["hover"], name=f, legendgroup=f))

    fig.update_layout(
        title="Receptor-Normalized Co-Activation Network",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleanchor="x", scaleratio=1),
        hovermode="closest", height=800, margin=dict(l=60, r=60, t=60, b=60),
        template="plotly", dragmode="pan", **BG,
    )
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("#### Filters")
    population = st.selectbox("Population (Race)", POPULATION_GROUPS, index=0)
    st.markdown("---")
    db_status = f"**Local DB** (`tcga_data.db`)" if DB_AVAILABLE else "No local DB — GDC API"
    st.caption(f"{db_status}. {len(RECEPTORS)} receptors · {len(LIGANDS)} ligands.")

p_threshold = 0.05

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

st.subheader("T Cell Exhaustion & TME Suppression Explorer", divider="gray")
st.caption(
    "All ligand values have **their receptor's own upregulation subtracted** to produce "
    "the Ligand-to-Receptor Activation Ratio (LRAR = Ligand log₂FC − Receptor log₂FC). "
    "This normalizes out pipeline batch effects between TCGA and GTEx, because the "
    "systematic bias cancels in the subtraction. "
    "LRAR > 0 means ligand outpaces receptor → TME is actively suppressing. "
    "LRAR < 0 means receptor outpaces ligand → pathway is ligand-limited."
)

st.caption(
    "Identifying which immune checkpoint receptors are co-activated by the tumor microenvironment. "
    "If the TME consistently upregulates the same groups of inhibitory ligands together, "
    "those receptor pathways fire as a unit — and blocking just one may not be enough."
)

col_proj, col_fam = st.columns([1, 3])
with col_proj:
    project_id = st.selectbox("Cancer Type", list(TCGA_PROJECTS.keys()),
        format_func=lambda x: f"{x} — {TCGA_PROJECTS[x]}",
        index=list(TCGA_PROJECTS.keys()).index("TCGA-SKCM"),
        help="TCGA cancer cohort to analyze.")
    st.caption(f"{TCGA_PROJECTS[project_id]}")

with col_fam:
    all_families = sorted(set(info["family"] for info in RECEPTORS.values()))
    selected_families = st.multiselect("Receptor Families", options=all_families,
        default=all_families, help="Filter receptor families.")
    active_families = set(selected_families) if selected_families else set(all_families)
    n_active_r = sum(1 for info in RECEPTORS.values() if info["family"] in active_families)
    st.caption(f"{n_active_r} receptors selected")

ACTIVE_RECEPTORS = {g: info for g, info in RECEPTORS.items() if info["family"] in active_families}

# Load data
progress_bar = st.progress(0, text="Loading data...")
try:
    receptor_df, ligand_df, demo_df, data_source = get_full_dataset(project_id, max_samples=0, progress_bar=progress_bar)
    progress_bar.empty()
except (GDCError, Exception) as e:
    progress_bar.empty()
    st.error(f"**Data loading failed:** {e}")
    st.stop()

with st.spinner("Fetching GTEx baselines..."):
    gtex_baseline, gtex_tissue = fetch_gtex_baseline(project_id)

if gtex_baseline.empty:
    st.error("GTEx baseline unavailable. Check internet connection.")
    st.stop()

# ── Hematologic / poor-baseline warnings ─────────────────────────────────────
HEMATOLOGIC_PROJECTS = {"TCGA-LAML", "TCGA-DLBC"}
APPROXIMATE_BASELINE_PROJECTS = {
    "TCGA-THYM": "No thymus tissue in GTEx. Using Blood Vessel as a rough proxy — results are approximate.",
    "TCGA-UVM": "No eye/uveal tissue in GTEx. Using Skin (melanocyte lineage) — results reflect lineage not organ.",
    "TCGA-MESO": "No pleural tissue in GTEx. Using Lung — mesothelioma originates in the lung lining.",
    "TCGA-SARC": "No soft tissue in GTEx. Using Adipose Subcutaneous — sarcomas are mesenchymal, not adipose.",
}

if project_id in HEMATOLOGIC_PROJECTS:
    if project_id == "TCGA-DLBC":
        st.warning(
            f"⚠ **Diffuse Large B-Cell Lymphoma** presents as solid tumor masses in lymph nodes, "
            f"but the malignant cells are B cells (hematologic origin). TCGA collected tissue biopsies "
            f"from lymph node masses — so the data includes stroma, vasculature, and infiltrating T cells, "
            f"meaning the TME framework partially applies."
        )
    else:  # LAML
        st.warning(
            f"⚠ **Acute Myeloid Leukemia** is a blood cancer in the bone marrow with no solid tumor "
            f"microenvironment. The GTEx baseline (Whole Blood) lacks the marrow niche where AML resides, "
            f"and many checkpoint molecules are constitutively expressed by normal myeloid cells — so "
            f"log₂FC values will appear inflated. However, some overexpressed signals are clinically "
            f"relevant: Zhao et al. (2025, *Signal Transduction and Targeted Therapy*) showed that AML "
            f"cells drive NK cell exhaustion through overactivation of the NKG2A/HLA-E axis, suppressing "
            f"the PI3K-AKT pathway — and that blocking NKG2A reversed this exhaustion both in vitro and "
            f"in vivo. [DOI: 10.1038/s41392-025-02228-5](https://www.nature.com/articles/s41392-025-02228-5)"
        )
elif project_id in APPROXIMATE_BASELINE_PROJECTS:
    st.info(
        f"ℹ **Approximate GTEx baseline:** {APPROXIMATE_BASELINE_PROJECTS[project_id]}"
    )

# Stage filter
if "stage" in demo_df.columns:
    available_stages = demo_df["stage"].value_counts()
    stage_options = ["All Stages"] + [s for s in ["Stage I", "Stage II", "Stage III", "Stage IV", "Stage 0/IS", "Not Reported"]
                                       if s in available_stages.index and available_stages[s] >= 1]
else:
    stage_options = ["All Stages"]

st.markdown(f"##### {TCGA_PROJECTS[project_id]} · {population}")

col_stage, col_thresh = st.columns([1, 1], gap="small")
with col_stage:
    selected_stage = st.segmented_control("Cancer Stage (AJCC Pathologic)", options=stage_options, default="All Stages",
        help="Filter patients by AJCC pathologic stage at diagnosis.")
    if selected_stage is None:
        selected_stage = "All Stages"
with col_thresh:
    corr_threshold = st.slider("Min |ρ| correlation threshold", 0.1, 0.8, 0.60, 0.05,
        help="Minimum absolute Spearman ρ between LRAR-normalized activation scores to draw an edge. Edges also require p < 0.05.")

st.caption(f"● {data_source} · GTEx baseline ({gtex_tissue}) · LRAR-normalized")

# Apply filters
filter_idx = demo_df.index
if population != "All" and "race_label" in demo_df.columns:
    filter_idx = filter_idx.intersection(demo_df[demo_df["race_label"] == population].index)
if selected_stage != "All Stages" and "stage" in demo_df.columns:
    filter_idx = filter_idx.intersection(demo_df[demo_df["stage"] == selected_stage].index)

r_filt = receptor_df.loc[receptor_df.index.intersection(filter_idx)]
l_filt = ligand_df.loc[ligand_df.index.intersection(filter_idx)]

if len(r_filt) < 10:
    st.warning(f"Only {len(r_filt)} samples — showing all.")
    r_filt, l_filt = receptor_df, ligand_df

active_rcols = [c for c in r_filt.columns if c in ACTIVE_RECEPTORS]
active_lcols = [c for c in l_filt.columns if c in {lg for lg, li in LIGANDS.items()
                if not li["receptors"] or any(r in ACTIVE_RECEPTORS for r in li["receptors"])}]
r_filt = r_filt[active_rcols]
l_filt = l_filt[active_lcols]

# Compute normalized scores
norm_df = compute_normalized_scores(r_filt, l_filt, gtex_baseline)
edge_df = compute_normalized_coexpression(r_filt, l_filt, gtex_baseline, corr_threshold, p_threshold)

# Metrics
c1, c2, c3, c4 = st.columns(4)
c1.metric("Samples", len(r_filt))
c2.metric("Receptors", len(ACTIVE_RECEPTORS))
c3.metric("Mean LRAR", f"{norm_df['LRAR'].mean():.2f}" if not norm_df.empty else "—")
avg_rho = f"{edge_df['rho'].abs().mean():.3f}" if len(edge_df) > 0 else "—"
c4.metric("Mean |ρ|", avg_rho)

# Network
net_fig = create_normalized_network(edge_df, norm_df, l_filt, active_receptors=ACTIVE_RECEPTORS)
st.plotly_chart(net_fig, use_container_width=True, config={"scrollZoom": True})
st.caption(
    "**Node size = mean LRAR** (higher = TME ligand supply outpaces receptor). "
    "**Blue** = co-activation (+ρ). **Red** = inverse (−ρ). "
    "**Dashed gray line at LRAR=0** = balanced."
)

# Tabs
tab0, tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Normalized Ridgeline", "🎯 LRAR Breakdown", "🔥 Correlation Matrix",
    "📋 Edge Table", "👥 Demographics", "🧪 Receptor–Ligand Pairs"
])

with tab0:
    st.markdown("### Per-Receptor LRAR Distribution")
    st.markdown(
        f"**n = {len(r_filt)} patients.** Each row is a receptor. Values show "
        "**Ligand log₂FC − Receptor log₂FC** per patient. "
        "**LRAR > 0** = ligand outpaces receptor (TME actively suppressing). "
        "**LRAR < 0** = receptor outpaces ligand (pathway is ligand-limited)."
    )
    st.plotly_chart(create_normalized_ridgeline(r_filt, l_filt, gtex_baseline, project_id, r_filt),
                    use_container_width=True)
    with st.expander("How to read this chart"):
        st.markdown(
            "Same ridgeline format as the main Explorer, but values are LRAR not raw log₂FC.\n\n"
            "- **Dashed line at 0** = balanced (ligand and receptor equally upregulated).\n"
            "- **Curves shifted right of 0** = TME is producing more ligand than the receptor can absorb.\n"
            "- **Curves shifted left of 0** = receptor is upregulated but ligand supply is limited.\n"
            "- Pipeline batch effects are normalized out because both numerator and denominator "
            "are affected equally by the TCGA vs GTEx quantification difference."
        )

with tab1:
    st.markdown("### LRAR Breakdown by Receptor")
    st.markdown("Each bar = Ligand log₂FC − Receptor log₂FC. Dashed line at 0 = balanced.")
    if not norm_df.empty:
        st.plotly_chart(create_normalized_breakdown(norm_df, project_id), use_container_width=True)
    with st.expander("How to read this chart"):
        st.markdown(
            "- **Bar > 0** = ligand is more upregulated than its receptor → TME is actively engaging this pathway.\n"
            "- **Bar < 0** = receptor is more upregulated than ligand → pathway may be receptor-limited.\n"
            "- **Bar ≈ 0** = proportional upregulation.\n"
            "- Hover shows both the raw ligand and receptor log₂FC values."
        )

with tab2:
    st.plotly_chart(create_normalized_corrmatrix(r_filt, l_filt, gtex_baseline), use_container_width=True)
    with st.expander("How to read this chart"):
        st.markdown(
            "Same as main Explorer correlation matrix, but using LRAR-normalized activation scores. "
            "Blue = co-activation, Red = inverse. "
            "Pipeline batch effects are factored out."
        )

with tab3:
    if len(edge_df) > 0:
        d = edge_df.copy()
        d["Receptor A"] = d["source"].map(lambda g: RECEPTORS[g]["label"])
        d["Receptor B"] = d["target"].map(lambda g: RECEPTORS[g]["label"])
        d = d[["Receptor A", "ligands_A", "mean_score_A", "Receptor B", "ligands_B", "mean_score_B", "weight", "rho", "pval"]]
        d.columns = ["Receptor A", "A's Ligands", "A Mean LRAR",
                      "Receptor B", "B's Ligands", "B Mean LRAR",
                      "Weight (|ρ|)", "ρ", "p-value"]
        st.dataframe(d.sort_values("Weight (|ρ|)", ascending=False), use_container_width=True, height=500)
    else:
        st.info("No edges at current thresholds.")

with tab4:
    if not demo_df.empty:
        col_a, col_b = st.columns(2)
        with col_a:
            if "race_label" in demo_df.columns:
                rc = demo_df["race_label"].value_counts().reset_index(); rc.columns = ["Race", "Count"]
                st.plotly_chart(px.bar(rc, x="Race", y="Count", title=f"Samples by Race — {project_id}",
                    template="plotly").update_layout(height=400, xaxis_tickangle=30, **BG), use_container_width=True)
        with col_b:
            if "stage" in demo_df.columns:
                sc = demo_df["stage"].value_counts().reset_index(); sc.columns = ["Stage", "Count"]
                stage_order = ["Stage 0/IS", "Stage I", "Stage II", "Stage III", "Stage IV", "Not Reported"]
                sc["_sort"] = sc["Stage"].apply(lambda x: stage_order.index(x) if x in stage_order else 99)
                sc = sc.sort_values("_sort").drop(columns="_sort")
                st.plotly_chart(px.bar(sc, x="Stage", y="Count", title=f"Samples by Stage — {project_id}",
                    template="plotly", color="Stage").update_layout(height=400, showlegend=False, **BG),
                    use_container_width=True)

with tab5:
    st.markdown("### Receptor–Ligand LRAR Table")
    if not norm_df.empty:
        display_cols = ["Receptor", "Ligand", "Receptor log₂FC", "Ligand log₂FC", "LRAR"]
        st.dataframe(norm_df[display_cols].sort_values("LRAR", ascending=False),
                     use_container_width=True, height=500)
        with st.expander("How to read this table"):
            st.markdown(
                "- **Receptor log₂FC** = receptor upregulation in tumor vs GTEx normal.\n"
                "- **Ligand log₂FC** = ligand upregulation.\n"
                "- **LRAR** = Ligand − Receptor. Positive = ligand outpaces receptor; negative = receptor outpaces ligand.\n"
                "- Both log₂FC values may be inflated due to pipeline differences — "
                "the LRAR normalizes this out."
            )

# Population comparison
st.markdown("---")
st.markdown("## Population-Level LRAR")
st.caption("Mean LRAR per receptor across racial/ethnic groups.")
if not demo_df.empty and "race_label" in demo_df.columns:
    pm_rows = []
    for rl in demo_df["race_label"].unique():
        rc2 = demo_df[demo_df["race_label"] == rl].index
        r_pop = receptor_df.loc[receptor_df.index.intersection(rc2)]
        l_pop = ligand_df.loc[ligand_df.index.intersection(rc2)]
        if len(r_pop) < 5:
            continue
        pop_norm = compute_normalized_scores(r_pop, l_pop, gtex_baseline)
        for rg, rinfo in ACTIVE_RECEPTORS.items():
            rsub = pop_norm[pop_norm["Receptor Gene"] == rg]
            if not rsub.empty:
                pm_rows.append({
                    "Population": f"{rl} (n={len(r_pop)})",
                    "Receptor": rinfo["label"],
                    "LRAR": rsub["LRAR"].mean(),
                })
    if pm_rows:
        pm_df = pd.DataFrame(pm_rows)
        pf = px.bar(pm_df, x="Receptor", y="LRAR", color="Population", barmode="group",
                    title=f"LRAR per Receptor by Race — {project_id}", template="plotly")
        pf.update_layout(height=550, xaxis_tickangle=45, **BG)
        pf.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5)
        st.plotly_chart(pf, use_container_width=True)

st.markdown("---")
st.caption(f"TCGA Explorer · LRAR = Ligand log₂FC − Receptor log₂FC · {len(RECEPTORS)} receptors · {len(LIGANDS)} ligands")