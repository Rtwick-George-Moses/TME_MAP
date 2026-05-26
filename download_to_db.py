#!/usr/bin/env python3
"""
TCGA Data Downloader → SQLite
===============================
Downloads expression data (receptors + ligands) and demographics
for all TCGA projects from the GDC API and stores them in a local
SQLite database for offline use.

Usage:
    python download_to_db.py                   # all projects, default 300 samples each
    python download_to_db.py --max-samples 500 # more samples per project
    python download_to_db.py --projects TCGA-SKCM TCGA-BRCA  # specific projects only
    python download_to_db.py --db my_data.db   # custom output path

The resulting .db file is used by the Explorer app when present.
"""

import argparse
import sqlite3
import sys
import os
import time
import io
import tarfile
import json

import requests
import pandas as pd
import numpy as np

# Add parent dir so we can import config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RECEPTORS, LIGANDS, ALL_GENE_ENSEMBL, TCGA_PROJECTS,
    RACE_MAP, GDC_CASES_ENDPOINT, GDC_FILES_ENDPOINT, GDC_DATA_ENDPOINT,
)

DEFAULT_DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tcga_data.db")
BATCH_SIZE = 25


def gdc_post(endpoint, payload, retries=3, timeout=120):
    for attempt in range(retries):
        try:
            resp = requests.post(endpoint, json=payload,
                                 headers={"Content-Type": "application/json"}, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            if attempt < retries - 1:
                wait = 2 ** (attempt + 1)
                print(f"    Retry {attempt+1}/{retries} in {wait}s: {e}")
                time.sleep(wait)
            else:
                raise RuntimeError(f"GDC API failed after {retries} attempts: {e}")


def fetch_demographics(project_id):
    print(f"  Fetching demographics for {project_id}...")
    all_cases, bs, off = [], 500, 0
    while True:
        data = gdc_post(GDC_CASES_ENDPOINT, {
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
        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        for h in hits:
            demo = h.get("demographic", {})
            diag = h.get("diagnoses", [{}])
            d0 = diag[0] if diag else {}
            raw_stage = d0.get("ajcc_pathologic_stage") or d0.get("tumor_stage") or "Not Reported"
            all_cases.append({
                "case_id": h["case_id"],
                "submitter_id": h.get("submitter_id", ""),
                "race": (demo.get("race") or "not reported").lower(),
                "ethnicity": (demo.get("ethnicity") or "not reported").lower(),
                "gender": (demo.get("gender") or "not reported").lower(),
                "age_at_diagnosis_days": d0.get("age_at_diagnosis"),
                "stage_raw": raw_stage,
            })
        off += bs
        total = data["data"]["pagination"]["total"]
        print(f"    {min(off, total)}/{total} cases...")
        if off >= total:
            break
    print(f"    → {len(all_cases)} cases")
    return pd.DataFrame(all_cases) if all_cases else pd.DataFrame()


def fetch_file_ids(project_id):
    print(f"  Fetching file IDs for {project_id}...")
    all_files, bs, off = [], 500, 0
    while True:
        data = gdc_post(GDC_FILES_ENDPOINT, {
            "filters": {"op": "and", "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
                {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
                {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
                {"op": "=", "content": {"field": "data_format", "value": "TSV"}},
                {"op": "=", "content": {"field": "access", "value": "open"}},
            ]},
            "fields": "file_id,file_name,cases.case_id",
            "size": bs, "from": off,
        })
        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        for h in hits:
            for c in h.get("cases", []):
                all_files.append({"file_id": h["file_id"], "case_id": c["case_id"]})
        off += bs
        total = data["data"]["pagination"]["total"]
        if off >= total:
            break
    print(f"    → {len(all_files)} files")
    return pd.DataFrame(all_files) if all_files else pd.DataFrame()


def download_expression_batch(file_ids):
    """Download a batch and return dict of {file_id: {gene: log2_tpm}}."""
    target_ids = {eid.split(".")[0]: gene for gene, eid in ALL_GENE_ENSEMBL.items()}
    resp = requests.post(GDC_DATA_ENDPOINT, json={"ids": file_ids},
                         headers={"Content-Type": "application/json"}, timeout=300)
    resp.raise_for_status()

    records = {}
    with tarfile.open(fileobj=io.BytesIO(resp.content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".tsv"):
                continue
            parts = member.name.split("/")
            fid = parts[0] if len(parts) > 1 else member.name
            f = tar.extractfile(member)
            if not f:
                continue
            try:
                tsv = pd.read_csv(f, sep="\t", comment="#")
            except Exception:
                continue
            if "gene_id" not in tsv.columns:
                continue
            tpm_col = next((c for c in ["tpm_unstranded", "fpkm_unstranded", "unstranded"]
                            if c in tsv.columns), None)
            if not tpm_col:
                continue
            tsv["gid"] = tsv["gene_id"].str.split(".").str[0]
            sub = tsv.loc[tsv["gid"].isin(target_ids), ["gid", tpm_col]].copy()
            sub["sym"] = sub["gid"].map(target_ids)
            if not sub.empty:
                records[fid] = {row["sym"]: np.log2(float(row[tpm_col]) + 1)
                                for _, row in sub.iterrows()}
    return records


def download_all_expression(file_ids, max_samples):
    """Download expression in batches with progress."""
    fids = file_ids[:max_samples]
    total = len(fids)
    batches = [fids[i:i+BATCH_SIZE] for i in range(0, total, BATCH_SIZE)]

    all_records = {}
    for bi, batch in enumerate(batches):
        done = bi * BATCH_SIZE
        print(f"    Downloading batch {bi+1}/{len(batches)} ({done}/{total} complete)...")
        try:
            batch_records = download_expression_batch(batch)
            all_records.update(batch_records)
        except Exception as e:
            print(f"    ⚠ Batch {bi+1} failed: {e}. Skipping.")
            continue

    print(f"    → Extracted {len(all_records)} samples from {total} files")
    return all_records


def init_db(db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute("""
        CREATE TABLE IF NOT EXISTS demographics (
            project_id TEXT NOT NULL,
            case_id TEXT NOT NULL,
            submitter_id TEXT,
            race TEXT,
            ethnicity TEXT,
            gender TEXT,
            age_at_diagnosis_days INTEGER,
            stage_raw TEXT,
            PRIMARY KEY (project_id, case_id)
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS expression (
            project_id TEXT NOT NULL,
            case_id TEXT NOT NULL,
            gene TEXT NOT NULL,
            log2_tpm REAL,
            PRIMARY KEY (project_id, case_id, gene)
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)

    conn.commit()
    return conn


def save_project(conn, project_id, demo_df, expression_records, file_to_case):
    c = conn.cursor()

    # Clear existing data for this project
    c.execute("DELETE FROM demographics WHERE project_id = ?", (project_id,))
    c.execute("DELETE FROM expression WHERE project_id = ?", (project_id,))
    c.execute("DELETE FROM gtex_baseline WHERE project_id = ?", (project_id,))

    # Insert demographics
    if not demo_df.empty:
        for _, row in demo_df.iterrows():
            c.execute("""
                INSERT OR REPLACE INTO demographics
                (project_id, case_id, submitter_id, race, ethnicity, gender, age_at_diagnosis_days, stage_raw)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (project_id, row["case_id"], row.get("submitter_id", ""),
                  row.get("race", ""), row.get("ethnicity", ""),
                  row.get("gender", ""), row.get("age_at_diagnosis_days"),
                  row.get("stage_raw", "Not Reported")))

    # Insert expression (long format: one row per case × gene)
    expr_rows = []
    for file_id, gene_vals in expression_records.items():
        case_id = file_to_case.get(file_id)
        if not case_id:
            continue
        for gene, val in gene_vals.items():
            expr_rows.append((project_id, case_id, gene, val))

    c.executemany("""
        INSERT OR REPLACE INTO expression (project_id, case_id, gene, log2_tpm)
        VALUES (?, ?, ?, ?)
    """, expr_rows)

    conn.commit()
    print(f"    → Saved {len(demo_df)} demographics + {len(expr_rows)} expression rows")


def main():
    parser = argparse.ArgumentParser(description="Download TCGA data to SQLite for offline use")
    parser.add_argument("--db", default=DEFAULT_DB, help=f"Output database path (default: {DEFAULT_DB})")
    parser.add_argument("--max-samples", type=int, default=300, help="Max samples per project (default: 300)")
    parser.add_argument("--projects", nargs="*", default=None,
                        help="Specific projects to download (default: all)")
    args = parser.parse_args()

    projects = args.projects or list(TCGA_PROJECTS.keys())
    # Validate project names
    for p in projects:
        if p not in TCGA_PROJECTS:
            print(f"ERROR: Unknown project '{p}'. Valid: {', '.join(TCGA_PROJECTS.keys())}")
            sys.exit(1)

    print(f"═══════════════════════════════════════════════════════")
    print(f"TCGA Data Downloader → SQLite")
    print(f"Database: {args.db}")
    print(f"Projects: {len(projects)}")
    print(f"Max samples per project: {args.max_samples}")
    print(f"Genes tracked: {len(ALL_GENE_ENSEMBL)} ({len(RECEPTORS)} receptors + {len(LIGANDS)} ligands)")
    print(f"═══════════════════════════════════════════════════════\n")

    conn = init_db(args.db)

    success = 0
    failed = []

    for pi, project_id in enumerate(projects):
        print(f"\n[{pi+1}/{len(projects)}] {project_id} — {TCGA_PROJECTS[project_id]}")
        print(f"{'─' * 50}")

        try:
            # Demographics
            demo_df = fetch_demographics(project_id)
            if demo_df.empty:
                print(f"  ⚠ No cases found. Skipping.")
                failed.append((project_id, "No cases"))
                continue

            # File IDs
            files_df = fetch_file_ids(project_id)
            if files_df.empty:
                print(f"  ⚠ No expression files found. Skipping.")
                failed.append((project_id, "No files"))
                continue

            files_df = files_df.drop_duplicates(subset="case_id", keep="first")
            file_to_case = dict(zip(files_df["file_id"], files_df["case_id"]))

            # Expression
            fids = files_df["file_id"].tolist()
            expression_records = download_all_expression(fids, args.max_samples)

            if not expression_records:
                print(f"  ⚠ No expression data extracted. Skipping.")
                failed.append((project_id, "No expression"))
                continue

            # Save to DB
            save_project(conn, project_id, demo_df, expression_records, file_to_case)
            success += 1

        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            failed.append((project_id, str(e)))
            continue

    # Save metadata
    c = conn.cursor()
    c.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
              ("download_time", time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())))
    c.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
              ("max_samples", str(args.max_samples)))
    c.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
              ("gene_count", str(len(ALL_GENE_ENSEMBL))))
    c.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
              ("projects_downloaded", json.dumps([p for p in projects if p not in [f[0] for f in failed]])))
    conn.commit()
    conn.close()

    # Summary
    print(f"\n{'═' * 50}")
    print(f"DONE: {success}/{len(projects)} projects downloaded")
    if failed:
        print(f"FAILED ({len(failed)}):")
        for p, reason in failed:
            print(f"  {p}: {reason}")
    db_size = os.path.getsize(args.db) / (1024 * 1024)
    print(f"Database: {args.db} ({db_size:.1f} MB)")
    print(f"{'═' * 50}")


if __name__ == "__main__":
    main()