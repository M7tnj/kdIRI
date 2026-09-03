from __future__ import annotations

import io
import os
import re
import sys
import json
import time
import gzip
import shutil
import platform
import warnings
import urllib.request
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Dict, Any

import numpy as np
import pandas as pd
import requests
import networkx as nx
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, to_tree
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib_venn import venn2, venn3


try:
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/NotoSansSC-Regular.ttf')
except Exception:
    pass

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Noto Sans SC']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="networkx")


# Configuration

RESULTS_DIR = Path("/workspaces/kdIRI/results/")
FIGURES_DIR = Path("/workspaces/kdIRI/figures/")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

STRING_API   = "https://string-db.org/api"
STRING_SPECIES = 10116              
STRING_SCORE   = 400               
STRING_CHUNK   = 250                
HUB6_VALIDATED = ["BRCA1", "RAD51", "MCM7", "EXO1", "RFC3", "TP53"]

SEED = 20240901
np.random.seed(SEED)


# Utility helpers

def save_fig(fig: plt.Figure, name: str, *, width_in: float = 6,
             height_in: float = 5) -> Path:
    path = FIGURES_DIR / name
    fig.set_size_inches(width_in, height_in)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")
    return path


def safe_save_df(df: pd.DataFrame, name: str) -> Path:
    path = RESULTS_DIR / name
    df.to_csv(path, index=False)
    print(f"  Saved: {path}")
    return path


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")


# STEP 12 STRING PPI network


STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download"
STRING_VERSION = "12.0"
STRING_CACHE_DIR = RESULTS_DIR / "STRING_cache"
STRING_CACHE_DIR.mkdir(parents=True, exist_ok=True)


def _download_with_retry(url: str, dest: Path,
                          retries: int = 3, timeout: int = 300) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    for attempt in range(1, retries + 1):
        try:
            log(f"  Downloading {url} (attempt {attempt}/{retries})")
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        fh.write(chunk)
                tmp.rename(dest)
            return True
        except Exception as e:
            log(f"  Download error: {e}")
            time.sleep(5 * attempt)
    return False


def load_string_info() -> pd.DataFrame:
    url = f"{STRING_DOWNLOAD_BASE}/protein.info.v{STRING_VERSION}/10116.protein.info.v{STRING_VERSION}.txt.gz"
    dest = STRING_CACHE_DIR / f"10116.protein.info.v{STRING_VERSION}.txt.gz"
    if not _download_with_retry(url, dest):
        raise RuntimeError(f"Failed to download STRING protein.info from {url}")
    df = pd.read_csv(dest, sep="\t", dtype=str)

    sid_col = "#string_protein_id" if "#string_protein_id" in df.columns else df.columns[0]
    df = df.rename(columns={sid_col: "STRING_id",
                             "preferred_name": "GeneSymbol"})
    return df[["STRING_id", "GeneSymbol"]].dropna().drop_duplicates()


def load_string_links() -> pd.DataFrame:
    url = f"{STRING_DOWNLOAD_BASE}/protein.links.v{STRING_VERSION}/10116.protein.links.v{STRING_VERSION}.txt.gz"
    dest = STRING_CACHE_DIR / f"10116.protein.links.v{STRING_VERSION}.txt.gz"
    if not _download_with_retry(url, dest):
        raise RuntimeError(f"Failed to download STRING protein.links from {url}")
    df = pd.read_csv(dest, sep=r"\s+", dtype=str)
    df = df.rename(columns={"protein1": "fromId",
                             "protein2": "toId",
                             "combined_score": "combined_score"})
    df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")
    return df[["fromId", "toId", "combined_score"]].dropna().drop_duplicates()


def map_symbols_to_string_ids(symbols: List[str],
                              info_df: pd.DataFrame) -> pd.DataFrame:
    info_lower = info_df.copy()
    info_lower["sym_lower"] = info_lower["GeneSymbol"].str.lower()
    sym_to_sid = dict(zip(info_lower["sym_lower"], info_lower["STRING_id"]))
    sym_to_str = dict(zip(info_lower["sym_lower"], info_lower["GeneSymbol"]))

    rows = []
    for s in symbols:
        key = s.strip().lower()
        if key in sym_to_sid:
            rows.append({
                "GeneSymbol_input": s,
                "STRING_id":        sym_to_sid[key],
                "GeneSymbol_STRING": sym_to_str[key],
            })
    return pd.DataFrame(rows)


def filter_string_edges(string_ids: List[str],
                         links_df: pd.DataFrame,
                         min_score: int = STRING_SCORE) -> pd.DataFrame:

    id_set = set(string_ids)
    edges = links_df[
        links_df["fromId"].isin(id_set)
        & links_df["toId"].isin(id_set)
        & (links_df["combined_score"] >= min_score)
    ].copy()
    if edges.empty:
        return edges

    edges["__key"] = edges.apply(
        lambda r: "|".join(sorted([str(r["fromId"]), str(r["toId"])])),
        axis=1
    )
    edges = edges.drop_duplicates(subset="__key").drop(columns="__key")
    return edges.reset_index(drop=True)


def step12_string_ppi() -> Tuple[pd.DataFrame, pd.DataFrame, nx.Graph,
                                  List[str], List[str]]:
   
    log("STEP 12 — STRING PPI network (bulk-file method)")

    up_csv = RESULTS_DIR / "Upregulated_DEGs_for_STRING.csv"
    if not up_csv.exists():
        raise FileNotFoundError(
            f"Required R output not found: {up_csv}\n"
        )

    up_df = pd.read_csv(up_csv)

    sym_col = next((c for c in ["GeneSymbol", "SYMBOL", "Gene"]
                    if c in up_df.columns), None)
    if sym_col is None:
        raise ValueError(
            f"No gene-symbol column found in {up_csv} "
            f"(columns: {list(up_df.columns)})"
        )

    symbols = up_df[sym_col].dropna().astype(str).str.strip()
    symbols = symbols[symbols != ""].unique().tolist()
    log(f"  Up-regulated gene symbols submitted to STRING: {len(symbols)}")

    info_df  = load_string_info()
    links_df = load_string_links()
    log(f"  STRING info:  {len(info_df):,} proteins")
    log(f"  STRING links: {len(links_df):,} interactions")

    mapped = map_symbols_to_string_ids(symbols, info_df)
    if mapped.empty:
        raise RuntimeError("No up-regulated gene symbols mapped to STRING ids.")
    mapped.to_csv(RESULTS_DIR / "STRING_mapped_upregulated.csv", index=False)
    log(f"  STRING mapped proteins: {len(mapped)} / {len(symbols)}")

    string_ids = mapped["STRING_id"].astype(str).unique().tolist()
    edges_df = filter_string_edges(string_ids, links_df, STRING_SCORE)
    if edges_df.empty:
        raise RuntimeError(
            "STRING returned no edges among the input proteins. "
        )
    log(f"  STRING edges among input proteins (score >= {STRING_SCORE}): "
        f"{len(edges_df)}")

    id2sym = dict(zip(mapped["STRING_id"], mapped["GeneSymbol_STRING"]))
    edges_df = edges_df.copy()
    edges_df["fromSymbol"] = edges_df["fromId"].map(id2sym)
    edges_df["toSymbol"]   = edges_df["toId"].map(id2sym)
    safe_save_df(edges_df, "STRING_edge_list_upregulated.csv")

    G = nx.Graph()
    for _, r in edges_df.iterrows():
        a = r["fromSymbol"] if pd.notna(r["fromSymbol"]) else r["fromId"]
        b = r["toSymbol"]   if pd.notna(r["toSymbol"])   else r["toId"]
        G.add_edge(a, b, weight=float(r["combined_score"]))
    log(f"  NetworkX graph: {G.number_of_nodes()} vertices, "
        f"{G.number_of_edges()} edges")

    log("  Computing centralities …")
    degree_cent = dict(G.degree())
    btw_cent    = nx.betweenness_centrality(G, normalized=True)
    clo_cent    = nx.closeness_centrality(G)
    try:
        eig_cent = nx.eigenvector_centrality(G, max_iter=1000, tol=1e-6)
    except Exception:
        eig_cent = {n: np.nan for n in G.nodes}

    centralities = pd.DataFrame({
        "Gene":        list(G.nodes),
        "Degree":      [degree_cent[n] for n in G.nodes],
        "Betweenness": [btw_cent.get(n, 0.0)    for n in G.nodes],
        "Closeness":   [clo_cent.get(n, 0.0)    for n in G.nodes],
        "Eigenvector": [eig_cent.get(n, np.nan) for n in G.nodes],
    })
    centralities = centralities.sort_values("Degree", ascending=False) \
                               .reset_index(drop=True)
    safe_save_df(centralities, "hub_gene_centrality_metrics.csv")

    hub90 = centralities.head(90)["Gene"].tolist()
    hub10 = centralities.head(10)["Gene"].tolist()
    (RESULTS_DIR / "hub_genes_90.txt").write_text("\n".join(hub90))
    (RESULTS_DIR / "hub_genes_10.txt").write_text("\n".join(hub10))

    (RESULTS_DIR / "hub_genes_6_validated.txt").write_text("\n".join(HUB6_VALIDATED))

    log(f"  hub90 saved ({len(hub90)} genes); hub10 saved ({len(hub10)} genes)")
    return centralities, edges_df, G, hub10, hub90



# STEP 13 hubgene supplementary bar plots

def step13_hub_bar_plots(centralities: pd.DataFrame) -> None:
    log("STEP 13 — Hub-gene bar plots (4 metrics x top-20/30)")

    metric_cols = {
        "Degree":      "#E41A1C",
        "Betweenness": "#377EB8",
        "Closeness":   "#4DAF4A",
        "Eigenvector": "#984EA3",
    }

    for metric, color in metric_cols.items():
        for n in (20, 30):
            df = centralities.dropna(subset=[metric]) \
                             .nlargest(n, metric) \
                             .copy()
            df = df.iloc[::-1]
            df["Gene"] = pd.Categorical(df["Gene"], categories=df["Gene"], ordered=True)

            fig, ax = plt.subplots()
            ax.barh(df["Gene"], df[metric], color=color,
                    edgecolor="black", linewidth=0.5)
            ax.set_xlabel(metric)
            ax.set_ylabel(None)
            ax.set_title(f"Top-{n} hub genes by {metric}")
            ax.tick_params(axis='y', labelsize=7)
            save_fig(fig, f"HubFig_top{n}_{metric}.png",
                     width_in=6, height_in=max(4, n * 0.25))



# 14| HUB OF HUB refinement

def step14_hub10_subnetwork(hub10: List[str]) -> None:
    log("STEP 14 — Hub-10 sub-network refinement")
    if not hub10:
        log("  hub10 is empty — skipping.")
        return
    
    info_df  = load_string_info()
    links_df = load_string_links()
    mapped   = map_symbols_to_string_ids(hub10, info_df)
    if mapped.empty:
        log("  No hub10 genes mapped to STRING ids — skipping edge query.")
        empty = pd.DataFrame(columns=["fromId", "toId", "combined_score",
                                       "fromSymbol", "toSymbol"])
        safe_save_df(empty, "STRING_edge_list_hub10.csv")
        return
    string_ids = mapped["STRING_id"].astype(str).unique().tolist()
    edges = filter_string_edges(string_ids, links_df, STRING_SCORE)
    id2sym = dict(zip(mapped["STRING_id"], mapped["GeneSymbol_STRING"]))
    if not edges.empty:
        edges = edges.copy()
        edges["fromSymbol"] = edges["fromId"].map(id2sym)
        edges["toSymbol"]   = edges["toId"].map(id2sym)
    safe_save_df(edges, "STRING_edge_list_hub10.csv")
    log(f"  Hub-10 sub-network edges: {len(edges)}")



# 15| miRNA target prediction Venn

def predict_venn(mirdb_file: Optional[str] = None,
                 mirwalk_file: Optional[str] = None,
                 targetscan_file: Optional[str] = None,
                 out_csv: str = "miRNA_intersection.csv",
                 out_png: str = "Fig4D_miRNA_Venn.png") -> List[str]:
    log("STEP 15 — miRNA target-prediction Venn (template)")

    def _load(p: Optional[str]) -> List[str]:
        if not p or not Path(p).exists():
            return []
        df = pd.read_csv(p)
        if "miRNA" not in df.columns:
            raise ValueError(f"{p} must contain a 'miRNA' column")
        return df["miRNA"].dropna().astype(str).unique().tolist()

    sets = {
        "miRDB":      _load(mirdb_file),
        "miRWalk":    _load(mirwalk_file),
        "TargetScan": _load(targetscan_file),
    }
    non_empty = {k: v for k, v in sets.items() if len(v) > 0}

    fig, ax = plt.subplots()
    if len(non_empty) == 0:
        ax.text(0.5, 0.5, "No miRNA database CSVs provided\n(template only)",
                ha="center", va="center", fontsize=11)
        ax.axis("off")
    elif len(non_empty) == 1:
        ax.text(0.5, 0.5,
                f"{list(non_empty)[0]}: {len(next(iter(non_empty.values())))} miRNAs",
                ha="center", va="center")
        ax.axis("off")
    elif len(non_empty) == 2:
        venn2(list(non_empty.values()), set_labels=list(non_empty.keys()), ax=ax)
    else:
        venn3(list(non_empty.values()), set_labels=list(non_empty.keys()), ax=ax)
    ax.set_title("miRNA target prediction — database overlap")
    save_fig(fig, out_png)

    common = list(reduce_intersect(list(non_empty.values())))
    safe_save_df(pd.DataFrame({"miRNA": common}), out_csv)
    log(f"  Common miRNAs across provided DBs: {len(common)}")
    return common


def reduce_intersect(lists: List[List[str]]) -> List[str]:
    if not lists:
        return []
    out = set(lists[0])
    for nxt in lists[1:]:
        out &= set(nxt)
    return sorted(out)


def step15_mirna_table() -> None:
    table = pd.DataFrame({
        "Gene": ["BRCA1", "BRCA1", "BRCA1",
                 "RAD51", "RAD51",
                 "EXO1",
                 "MCM7",
                 "RFC3",
                 "TP53"],
        "miRNA": ["rno-miR-103-3p", "rno-miR-107-3p", "rno-miR-101b-3p",
                  "rno-miR-92a-2-5p", "rno-miR-935",
                  "rno-miR-103-3p",
                  "rno-miR-107-3p",
                  "rno-miR-935",
                  "rno-miR-103-3p"],
        "Predicted_by": ["miRDB, miRWalk, TargetScan"] * 9,
        "Evidence_type": ["predicted (3-database overlap)"] * 9,
    })
    safe_save_df(table, "gene_miRNA_prediction_table.csv")


# STEP 16| Spearman miRNA-mRNA correlation template
def run_mirna_mrna_spearman(ct_df: pd.DataFrame,
                            gene_cols: List[str],
                            mirna_cols: List[str],
                            out_csv: str = "miRNA_mRNA_spearman.csv",
                            out_png: str = "Fig_miRNA_mRNA_corr_heatmap.png") -> pd.DataFrame:
    log("STEP 16 — Spearman miRNA-mRNA correlation (template)")

    missing = [c for c in gene_cols + mirna_cols if c not in ct_df.columns]
    if missing:
        raise ValueError(f"ct_df is missing columns: {missing}")

    rows = []
    for g in gene_cols:
        for m in mirna_cols:
            rho, p = spearmanr(ct_df[g], ct_df[m], nan_policy="omit")
            rows.append({"Gene": g, "miRNA": m,
                         "rho": float(rho), "p.value": float(p)})
    res = pd.DataFrame(rows)
    res["adj.p.value"] = multipletests(res["p.value"].fillna(1.0),
                                       method="fdr_bh")[1]
    safe_save_df(res, out_csv)

    wide = res.pivot(index="Gene", columns="miRNA", values="rho")
    fig, ax = plt.subplots()
    im = ax.imshow(wide.values, aspect="auto",
                   cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(wide.columns)))
    ax.set_xticklabels(wide.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(wide.index)))
    ax.set_yticklabels(wide.index, fontsize=8)
    ax.set_title("Spearman correlation: mRNA vs miRNA")
    fig.colorbar(im, ax=ax, label="Spearman ρ")
    save_fig(fig, out_png, width_in=max(5, len(wide.columns) * 0.6),
             height_in=max(4, len(wide.index) * 0.5))
    return res


# STEP 17: Independent validation GSE9943 (Reviewer 1 Q1)
def download_geo_series_matrix(gse: str, dest: Path) -> Optional[Path]:
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/matrix/"
    url = base + f"{gse}_series_matrix.txt.gz"
    log(f"  Downloading {url}")
    try:
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 14):
                    fh.write(chunk)
        return dest
    except Exception as e:
        log(f"  Download failed: {e}")
        return None


def parse_geo_series_matrix(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, str]:

    opener = gzip.open if path.suffix == ".gz" else open
    sample_cols, expr_rows, platform = [], [], ""
    series_matrix_start = False
    pheno_lines: List[str] = []

    with opener(path, "rt", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Series_platform_id"):
                platform = line.split("\t", 1)[-1].strip('"').split("-")[0]
                platform = line.split("\t", 1)[-1].strip('"')
            elif line.startswith("!Sample_"):
                pheno_lines.append(line)
            elif line.startswith("!series_matrix_table_begin"):
                series_matrix_start = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                series_matrix_start = False
                break
            elif series_matrix_start:
                if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                    sample_cols = [c.strip('"') for c in line.split("\t")[1:]]
                elif line and not line.startswith("!"):
                    parts = line.split("\t")
                    expr_rows.append(parts)

    if not expr_rows:
        raise ValueError(f"No expression matrix found in {path}")

    ids = [r[0].strip('"') for r in expr_rows]
    values = np.array(
        [[float(v) if v not in ("", "NA", "null") else np.nan
          for v in r[1:]] for r in expr_rows],
        dtype=float
    )
    expr_df = pd.DataFrame(values, index=ids, columns=sample_cols)

    pheno_records: Dict[str, Dict[str, str]] = {}
    for line in pheno_lines:
        key, *vals = line.split("\t")
        key = key.lstrip("!").strip()
        vals = [v.strip('"') for v in vals]
        for i, col in enumerate(sample_cols):
            pheno_records.setdefault(col, {})[key] = vals[i] if i < len(vals) else ""
    pheno_df = pd.DataFrame.from_dict(pheno_records, orient="index")
    pheno_df.index.name = "GSM"

    return expr_df, pheno_df, platform


def simple_limma(expr: pd.DataFrame, groups: List[str]) -> pd.DataFrame:
    groups = pd.Series(groups, index=expr.columns)
    levels = list(dict.fromkeys(groups))
    if len(levels) < 2:
        raise ValueError("Need at least 2 group levels for DE analysis")

    ctrl_samples = [c for c, g in groups.items() if g == levels[0]]
    case_samples = [c for c, g in groups.items() if g == levels[1]]
    log(f"  DE: {levels[0]}={len(ctrl_samples)} vs {levels[1]}={len(case_samples)}")

    e_ctrl = expr[ctrl_samples].values.astype(float)
    e_case = expr[case_samples].values.astype(float)
    n_ctrl, n_case = e_ctrl.shape[1], e_case.shape[1]

    with np.errstate(invalid="ignore"):
        import warnings as _w
        with _w.catch_warnings():
            _w.filterwarnings("ignore", message="Mean of empty slice")
            _w.filterwarnings("ignore", message="Degrees of freedom")
            mean_ctrl = np.nanmean(e_ctrl, axis=1)
            mean_case = np.nanmean(e_case, axis=1)
            log_fc = mean_case - mean_ctrl
            var_ctrl = np.nanvar(e_ctrl, axis=1, ddof=1)
            var_case = np.nanvar(e_case, axis=1, ddof=1)
    se = np.sqrt(var_ctrl / n_ctrl + var_case / n_case)
    se = np.where(se == 0, np.nan, se)
    t_stat = log_fc / se

    num = (var_ctrl / n_ctrl + var_case / n_case) ** 2
    den = ((var_ctrl / n_ctrl) ** 2 / max(n_ctrl - 1, 1)
           + (var_case / n_case) ** 2 / max(n_case - 1, 1))
    df = num / np.where(den == 0, np.nan, den)

    with np.errstate(invalid="ignore"):
        p_vals = 2 * stats.t.sf(np.abs(t_stat), df=df)

    p_vals = np.where(np.isnan(p_vals), 1.0, p_vals)
    adj = multipletests(p_vals, method="fdr_bh")[1]

    return pd.DataFrame({
        "ProbeID":   expr.index,
        "AveExpr":   (mean_ctrl + mean_case) / 2,
        "logFC":     log_fc,
        "t":         t_stat,
        "P.Value":   p_vals,
        "adj.P.Val": adj,
    }).set_index("ProbeID")


def auto_extract_groups_gse9943(pheno: pd.DataFrame) -> List[str]:
    cand_cols = [c for c in pheno.columns
                 if c in ("Sample_title", "Sample_source_name_ch1",
                          "Sample_characteristics_ch1",
                          "Sample_characteristics_ch1.1",
                          "Sample_characteristics_ch1.2")]
    if not cand_cols:
        cand_cols = [c for c in pheno.columns
                     if c.startswith("Sample_characteristics")
                     or c == "Sample_title"]

    for col in cand_cols:
        vals = pheno[col].astype(str).str.lower()
        is_ctrl = vals.str.contains(r"control|sham", regex=True, na=False)
        is_case = vals.str.contains(r"i/r|iri|ischemia|aki|reperfusion",
                                     regex=True, na=False)
        if is_ctrl.any() and is_case.any() and not (is_ctrl & is_case).any():
            return ["Control" if c else "IRI"
                    for c in is_ctrl.tolist()]

    n = len(pheno)
    log(f"  Could not auto-detect groups; splitting {n} samples in half "
        f"(verify manually).")
    return ["Control"] * (n // 2) + ["IRI"] * (n - n // 2)


def step17_gse9943_validation() -> None:
    log("STEP 17 — GSE9943 independent validation")
    disc_csv = RESULTS_DIR / "DEG_full_all_genes.csv"
    if not disc_csv.exists():
        log(f"  {disc_csv} not found — skipping validation.")
        return
    disc = pd.read_csv(disc_csv)
    disc_syms = disc.get("GeneSymbol")
    if disc_syms is None:
        log("  DEG_full_all_genes.csv has no 'GeneSymbol' column — skipping.")
        return
    disc_syms = disc_syms.dropna().astype(str)
    disc_syms = disc_syms[disc_syms != ""].unique().tolist()
    log(f"  Discovery gene symbols: {len(disc_syms)}")

    gse_id = "GSE9943"
    dest = RESULTS_DIR / f"{gse_id}_series_matrix.txt.gz"
    if not dest.exists():
        ok = download_geo_series_matrix(gse_id, dest)
        if ok is None:
            log("  GSE9943 download failed — skipping validation.")
            return

    expr_val, pheno_val, platform = parse_geo_series_matrix(dest)
    log(f"  GSE9943: {expr_val.shape[0]} probes × {expr_val.shape[1]} samples "
        f"(platform {platform})")

    groups_val = auto_extract_groups_gse9943(pheno_val)
    log(f"  GSE9943 groups: {dict(pd.Series(groups_val).value_counts())}")

    deg_val = simple_limma(expr_val, groups_val)
    deg_val_filt = deg_val[(deg_val["logFC"].abs() > 2)
                           & (deg_val["P.Value"] < 0.05)].copy()
    safe_save_df(deg_val.reset_index(),       "GSE9943_DEG_full.csv")
    safe_save_df(deg_val_filt.reset_index(), "GSE9943_DEG_filtered.csv")

    sym_col = None
    for c in pheno_val.columns:
        if c.lower().startswith("sample_gene_platform") or "symbol" in c.lower():
            sym_col = c
            break

    val_syms: List[str] = []
    if "GeneSymbol" in deg_val_filt.columns:
        val_syms = deg_val_filt["GeneSymbol"].dropna().astype(str).unique().tolist()
    if not val_syms:
        log("  No probe→symbol mapping in series matrix — "
            "attempting GPL annotation download.")
        val_syms = _try_map_via_gpl(deg_val_filt.index.tolist(), platform)

    if not val_syms:
        log("  Could not map GSE9943 probes to symbols — "
            "overlap will be computed on probe IDs (less informative).")
        val_syms = deg_val_filt.index.astype(str).tolist()

    common = sorted(set(disc_syms) & set(val_syms))
    safe_save_df(pd.DataFrame({"Validated_gene": common}),
                 "GSE58438_vs_GSE9943_overlap_genes.csv")

    fig, ax = plt.subplots()
    venn2(subsets=(len(set(disc_syms) - set(val_syms)),
                   len(set(val_syms) - set(disc_syms)),
                   len(common)),
          set_labels=("GSE58438 (discovery)", "GSE9943 (validation)"),
          ax=ax)
    ax.set_title("GSE58438 vs GSE9943 DEG overlap")
    save_fig(fig, "Fig_GSE9943_validation_venn.png",
             width_in=6, height_in=5)


def _try_map_via_gpl(probe_ids: List[str], platform: str) -> List[str]:
    if not platform or not platform.startswith("GPL"):
        return []
    gpl_num = platform[3:]
    base = (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
            f"GPL{int(gpl_num)//1000}nnn/{platform}/annot/")
    urls = [base + f"{platform}.annot.gz",
            base + f"{platform}_annot.txt.gz"]
    dest = RESULTS_DIR / f"{platform}.annot.gz"
    if not dest.exists():
        ok = False
        for url in urls:
            log(f"  Downloading {url}")
            try:
                with requests.get(url, stream=True, timeout=180) as r:
                    if r.status_code != 200:
                        log(f"  HTTP {r.status_code}")
                        continue
                    with open(dest, "wb") as fh:
                        for chunk in r.iter_content(1 << 14):
                            fh.write(chunk)
                ok = True
                break
            except Exception as e:
                log(f"  Download error: {e}")
                continue
        if not ok:
            return []

    try:
        with gzip.open(dest, "rt", errors="ignore") as fh:
            lines = fh.readlines()
        header_idx = None
        for i, ln in enumerate(lines):
            if ln.startswith("!platform_table_begin"):
                header_idx = i + 1
                break
        if header_idx is None:
            for i, ln in enumerate(lines):
                if ln.startswith("ID\t"):
                    header_idx = i
                    break
        if header_idx is None:
            return []
        annot = pd.read_csv(io.StringIO("".join(lines[header_idx:])),
                            sep="\t", dtype=str, low_memory=False,
                            on_bad_lines="skip")
        if "ID" not in annot.columns:
            return []
        sym_col = next((c for c in ("Gene symbol", "Gene Symbol", "Gene_symbol",
                                     "Symbol", "SYMBOL", "GeneSymbol")
                        if c in annot.columns), None)
        if sym_col is None:
            sym_col = next((c for c in annot.columns if "symbol" in c.lower()),
                           None)
        if sym_col is None:
            return []
        annot[sym_col] = annot[sym_col].replace({"---": np.nan, "": np.nan})
        mapping = dict(zip(annot["ID"].astype(str),
                           annot[sym_col].astype(str)))
        syms = [mapping.get(p, np.nan) for p in probe_ids]
        syms = [s for s in syms if isinstance(s, str) and s.strip()
                and s.strip().lower() != "nan"]
        return list(dict.fromkeys(syms))
    except Exception as e:
        log(f"  GPL annotation parse failed: {e}")
        return []

# STEP 18 WGCNA
def load_expression_matrix() -> Tuple[pd.DataFrame, pd.Series]:

    saved = RESULTS_DIR / "GSE58438_expr_matrix.csv"
    if saved.exists():
        log(f"  Reading saved expression matrix from {saved}")
        expr = pd.read_csv(saved, index_col=0)
        sample_marker = "GSM1411057"
        if sample_marker in expr.columns:
            expr = expr.T
        elif sample_marker not in expr.index:
            log("  Warning: cannot find expected sample IDs in either axis.")
    else:
        log("  GSE58438_expr_matrix.csv not found in results/ — "
            "downloading from GEO.")
        dest = RESULTS_DIR / "GSE58438_series_matrix.txt.gz"
        if not dest.exists():
            ok = download_geo_series_matrix("GSE58438", dest)
            if ok is None:
                raise RuntimeError("GSE58438 download failed; cannot run WGCNA.")
        expr, pheno, _ = parse_geo_series_matrix(dest)
        expr = expr.T 
        expr.to_csv(saved)

    samples = ["GSM1411057", "GSM1411058", "GSM1411059", "GSM1411060",
               "GSM1411061", "GSM1411067", "GSM1411068", "GSM1411069",
               "GSM1411070"]
    available = [s for s in samples if s in expr.index]
    if len(available) != len(samples):
        log(f"  Warning: only {len(available)}/9 paper samples found in matrix.")
    expr = expr.loc[available]
    trait = pd.Series([0] * 5 + [1] * 4, index=available, name="AKI")
    log(f"  Expression matrix: {expr.shape[0]} samples × {expr.shape[1]} genes")
    return expr, trait


def good_samples_genes(datExpr: pd.DataFrame, frac_thresh: float = 0.95
                       ) -> pd.DataFrame:
    gfrac = datExpr.notna().mean(axis=0)
    good_genes = gfrac[gfrac >= frac_thresh].index
    datExpr = datExpr[good_genes].copy()

    means = datExpr.mean(axis=0, skipna=True)
    datExpr = datExpr.fillna(means)
    return datExpr


def pick_soft_threshold(datExpr: pd.DataFrame,
                        powers: Iterable[int] = range(1, 21)
                        ) -> pd.DataFrame:
    A = datExpr.corr().abs().values
    np.fill_diagonal(A, 0)
    records = []
    for p in powers:
        adj = A ** p
        k = adj.sum(axis=1) 
        hist, edges = np.histogram(k, bins=30)
        centers = (edges[:-1] + edges[1:]) / 2
        nonzero = hist > 0
        if nonzero.sum() < 5:
            records.append({"power": p, "SFT.R.sq": np.nan,
                            "slope": np.nan, "mean_connectivity": k.mean()})
            continue
        log_k  = np.log10(centers[nonzero])
        log_pk = np.log10(hist[nonzero] / hist.sum())
        slope, intercept, r, *_ = stats.linregress(log_k, log_pk)
        records.append({
            "power": p,
            "SFT.R.sq": r ** 2,
            "slope":   slope,
            "mean_connectivity": float(k.mean()),
        })
    return pd.DataFrame(records)


def build_tom(adj: np.ndarray) -> np.ndarray:
    n = adj.shape[0]
    k = adj.sum(axis=1) 
    L = adj @ adj  
    np.fill_diagonal(L, 0)
    min_k = np.minimum.outer(k, k)
    denom = min_k + 1 - adj
    tom = (L + adj) / denom
    np.fill_diagonal(tom, 1.0)
    return tom


def blockwise_modules(datExpr: pd.DataFrame, power: int = 6,
                      min_module_size: int = 30,
                      merge_cut_height: float = 0.25) -> Tuple[np.ndarray, pd.DataFrame]:
    log(f"  Building adjacency with power={power} …")
    A = datExpr.corr().abs().values ** power
    np.fill_diagonal(A, 0)

    log("  Computing TOM …")
    tom = build_tom(A)
    diss_tom = 1.0 - tom
    np.fill_diagonal(diss_tom, 0.0)

    diss_tom = (diss_tom + diss_tom.T) / 2
    diss_tom = np.clip(diss_tom, 0.0, 1.0)

    log("  Hierarchical clustering (average linkage) …")
    condensed = squareform(diss_tom, checks=False)
    Z = linkage(condensed, method="average")

    def _filter_small(labels_arr: np.ndarray) -> np.ndarray:
        out = labels_arr.copy()
        for lab in set(out):
            mask = out == lab
            if mask.sum() < min_module_size:
                out[mask] = 0
        return out

    best_cut: Optional[Tuple[float, np.ndarray]] = None
    initial_labels = None
    for t in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
        cand = fcluster(Z, t=t, criterion="distance")
        filtered = _filter_small(cand)
        n_modules = len(set(filtered) - {0})
        log(f"  t={t}: {n_modules} non-grey modules "
            f"(largest raw: {pd.Series(cand).value_counts().max()})")
        if 3 <= n_modules <= 15:
            best_cut = (t, cand)
            break
    if best_cut is not None:
        initial_labels = best_cut[1]
        log(f"  Using cut t={best_cut[0]}")
    else:
        initial_labels = fcluster(Z, t=8, criterion="maxclust")
        log("  No cut height gave 3-15 modules — forcing maxclust k=8")

    unique_labels = sorted(set(initial_labels))
    relabel = {lab: i + 1 for i, lab in enumerate(unique_labels)}
    labels = np.array([relabel[lab] for lab in initial_labels])

    labels = _filter_small(labels)
    def _module_eigengene(gene_indices: np.ndarray) -> np.ndarray:
        sub = datExpr.iloc[:, gene_indices].values
        sub = StandardScaler().fit_transform(sub)
        pca = PCA(n_components=1, random_state=SEED)
        me = pca.fit_transform(sub).ravel()

        if np.corrcoef(me, sub.mean(axis=1))[0, 1] < 0:
            me = -me
        return me

    module_ids = sorted(set(labels) - {0})
    if not module_ids:
        return labels, pd.DataFrame()

    MEs = {}
    for m in module_ids:
        idx = np.where(labels == m)[0]
        MEs[f"ME{m}"] = _module_eigengene(idx)
    MEs_df = pd.DataFrame(MEs, index=datExpr.index)

    if MEs_df.shape[1] > 1:
        me_corr = MEs_df.corr().abs().values
        np.fill_diagonal(me_corr, 0)
        merge = (1 - me_corr) < merge_cut_height
        parent = list(range(len(MEs_df.columns)))

        def _find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def _union(a: int, b: int) -> None:
            ra, rb = _find(a), _find(b)
            if ra != rb:
                parent[rb] = ra

        for i in range(len(MEs_df.columns)):
            for j in range(i + 1, len(MEs_df.columns)):
                if merge[i, j]:
                    _union(i, j)

        col_to_root = {i: _find(i) for i in range(len(MEs_df.columns))}
        root_me_cols: Dict[int, List[int]] = {}
        for col_idx, root in col_to_root.items():
            root_me_cols.setdefault(root, []).append(col_idx)

        new_labels = labels.copy()
        merged_mes: Dict[str, np.ndarray] = {}
        for root, members in root_me_cols.items():
            m_root_id = int(MEs_df.columns[root].replace("ME", ""))

            for j in members:
                m_j = int(MEs_df.columns[j].replace("ME", ""))
                new_labels[labels == m_j] = m_root_id

            merged_mes[f"ME{m_root_id}"] = MEs_df.iloc[:, members].mean(axis=1).values
        labels = new_labels
        MEs_df = pd.DataFrame(merged_mes, index=datExpr.index)

    return labels, MEs_df


def wgcna_module_colors(labels: np.ndarray) -> List[str]:
    import matplotlib.colors as mcolors
    standard = ["grey", "turquoise", "blue", "brown", "yellow", "green",
                "red", "black", "pink", "magenta", "purple", "greenyellow",
                "tan", "salmon", "cyan", "midnightblue", "lightcyan",
                "grey60", "lightgreen", "lightyellow", "royalblue",
                "darkred", "darkgreen", "darkturquoise", "darkgrey",
                "orange", "darkorange", "white", "skyblue", "saddlebrown",
                "steelblue", "paleturquoise", "violet", "darkolivegreen",
                "darkmagenta", "sienna3", "yellowgreen", "skyblue3",
                "plum1", "orangered4", "mediumpurple3", "ivory"]
    color_map = {0: "grey"}
    color_map.update({i + 1: standard[(i + 1) % len(standard)]
                      for i in range(200)})
    return [color_map.get(int(l), "white") for l in labels]


def step18_wgcna() -> None:
    log("STEP 18 — WGCNA module discovery")

    datExpr_full, trait = load_expression_matrix()
    datExpr = good_samples_genes(datExpr_full)
    log(f"  After goodSamplesGenes filter: {datExpr.shape}")

    if datExpr.shape[1] > 5000:
        variances = datExpr.var(axis=0, skipna=True)
        top_vars = variances.nlargest(5000).index
        datExpr = datExpr[top_vars]
        log(f"  Restricted to top 5,000 variable genes for performance.")

    variances = datExpr.var(axis=0, skipna=True)
    datExpr = datExpr.loc[:, variances > 1e-8]

    log("  pickSoftThreshold …")
    sft = pick_soft_threshold(datExpr, powers=range(1, 21))
    safe_save_df(sft, "WGCNA_soft_threshold_table.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.scatter(sft["power"], sft["SFT.R.sq"], color="red")
    for _, r in sft.iterrows():
        ax1.annotate(int(r["power"]), (r["power"], r["SFT.R.sq"]))
    ax1.set_xlabel("Soft Threshold (power)")
    ax1.set_ylabel("Scale Free Topology Model Fit (signed R²)")
    ax1.set_title("Soft-threshold selection")
    ax1.axhline(0.85, color="grey", linestyle="--", linewidth=0.5)
    ax2.plot(sft["power"], sft["mean_connectivity"], marker="o")
    ax2.set_xlabel("Soft Threshold (power)")
    ax2.set_ylabel("Mean Connectivity")
    ax2.set_title("Mean connectivity decay")
    save_fig(fig, "WGCNA_soft_threshold.png", width_in=10, height_in=4)

    candidates = sft[sft["SFT.R.sq"] > 0.85]
    chosen_power = int(candidates["power"].iloc[0]) if not candidates.empty else 6
    log(f"  Chosen power = {chosen_power}")

    labels, MEs = blockwise_modules(datExpr, power=chosen_power,
                                     min_module_size=30,
                                     merge_cut_height=0.25)
    colors = wgcna_module_colors(labels)
    safe_save_df(pd.DataFrame({"ProbeID": datExpr.columns,
                                "Module": labels,
                                "ModuleColor": colors}),
                 "WGCNA_module_assignments.csv")

    if not datExpr.empty:
        A = datExpr.corr().abs().values ** chosen_power
        np.fill_diagonal(A, 0)
        tom = build_tom(A)
        diss = 1.0 - tom
        np.fill_diagonal(diss, 0.0)
        diss = (diss + diss.T) / 2
        diss = np.clip(diss, 0, 1)
        Z = linkage(squareform(diss, checks=False), method="average")

        fig, ax = plt.subplots(figsize=(10, 6))
        dendrogram(Z, no_labels=True, color_threshold=0.99, ax=ax,
                   truncate_mode="lastp", p=50, show_leaf_counts=True)
        ax.set_title("WGCNA dendrogram (last 50 merges) with module colors")
        save_fig(fig, "WGCNA_dendro_colors.png", width_in=10, height_in=6)

    if not MEs.empty:
        n = len(datExpr)
        rows = []
        for me in MEs.columns:
            r, p = pearsonr(MEs[me], trait.loc[MEs.index])
            rows.append({"Module": me, "cor": r, "p.value": p})
        mtr = pd.DataFrame(rows)
        mtr["adj.p.value"] = multipletests(mtr["p.value"].fillna(1.0),
                                           method="fdr_bh")[1]
        safe_save_df(mtr, "WGCNA_module_trait_correlations.csv")

        fig, ax = plt.subplots(figsize=(3, max(4, 0.4 * len(MEs.columns))))
        cmap = plt.cm.RdBu_r
        im = ax.imshow(mtr[["cor"]].values.reshape(-1, 1),
                       cmap=cmap, vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks([0]); ax.set_xticklabels(["AKI vs Control"], rotation=45,
                                                ha="right")
        ax.set_yticks(range(len(mtr))); ax.set_yticklabels(mtr["Module"])
        for i, (_, r) in enumerate(mtr.iterrows()):
            txt = f"{r['cor']:.2f}\n({r['p.value']:.1g})"
            ax.text(0, i, txt, ha="center", va="center", fontsize=7)
        ax.set_title("Module-trait relationships")
        fig.colorbar(im, ax=ax, label="Pearson r")
        save_fig(fig, "WGCNA_module_trait_heatmap.png", width_in=3.5,
                 height_in=max(4, 0.5 * len(MEs.columns) + 1))

        sizes = pd.Series(labels).value_counts()
        sizes = sizes[sizes.index != 0]
        if not sizes.empty:
            module_of_interest = int(sizes.idxmax())
            mod_color = wgcna_module_colors([module_of_interest])[0]
            log(f"  Largest module: {module_of_interest} ({mod_color}), "
                f"{sizes.max()} genes")
            mod_genes_mask = labels == module_of_interest
            mod_data = datExpr.loc[:, mod_genes_mask]
            me_col = f"ME{module_of_interest}"
            if me_col in MEs.columns:
                gmm = mod_data.apply(lambda g: pearsonr(g, MEs[me_col])[0])
                hub10_wgcna = gmm.nlargest(10).index.tolist()
                (RESULTS_DIR / "WGCNA_hub10.txt").write_text(
                    "\n".join(hub10_wgcna))
                log(f"  Top-10 hub genes (module {mod_color}): "
                    f"{', '.join(hub10_wgcna[:5])} …")


# 19 Exploratory ROC helper
def exploratory_roc(df: pd.DataFrame, marker: str,
                    group_col: str = "group",
                    group_levels: Tuple[str, str] = ("Sham", "IRI")
                    ) -> pd.DataFrame:
    df = df.copy()
    df["group_bin"] = df[group_col].map({group_levels[0]: 0,
                                          group_levels[1]: 1})
    sub = df[[marker, "group_bin"]].dropna()
    if sub.empty or sub["group_bin"].nunique() < 2:
        return pd.DataFrame()

    fpr, tpr, thr = roc_curve(sub["group_bin"], sub[marker])
    auc_val = auc(fpr, tpr)

    j = tpr - fpr
    best_idx = int(np.argmax(j))
    best_thr = thr[best_idx]
    best_sens = tpr[best_idx]
    best_spec = 1 - fpr[best_idx]
    lr_pos = best_sens / max(1 - best_spec, 1e-6)

    rng = np.random.default_rng(SEED)
    boot_aucs = []
    for _ in range(1000):
        idx = rng.integers(0, len(sub), len(sub))
        if sub["group_bin"].iloc[idx].nunique() < 2:
            continue
        try:
            f_, t_, _ = roc_curve(sub["group_bin"].iloc[idx],
                                  sub[marker].iloc[idx])
            boot_aucs.append(auc(f_, t_))
        except Exception:
            continue
    if boot_aucs:
        ci_low, ci_high = np.percentile(boot_aucs, [2.5, 97.5])
    else:
        ci_low = ci_high = np.nan

    cases = sub.loc[sub["group_bin"] == 1, marker]
    ctrls = sub.loc[sub["group_bin"] == 0, marker]
    try:
        _, p_val = stats.mannwhitneyu(cases, ctrls, alternative="two-sided")
    except Exception:
        p_val = np.nan

    return pd.DataFrame([{
        "marker":            marker,
        "n":                 int(len(sub)),
        "AUC":               float(auc_val),
        "AUC_lower95":       float(ci_low),
        "AUC_upper95":       float(ci_high),
        "p_value":           float(p_val),
        "cutoff":            float(best_thr),
        "sensitivity":       float(best_sens),
        "specificity":       float(best_spec),
        "likelihood_ratio":  float(lr_pos),
    }])


def step19_roc_template() -> None:
    log("STEP 19 — Exploratory ROC helper (template; no data)")
    src = (
        "def exploratory_roc(df, marker, group_col='group',\n"
        "                    group_levels=('Sham','IRI')):\n"
        "    # See proj_steps12_20.py for full implementation.\n"
        "    # Returns: marker, n, AUC, AUC_lower95, AUC_upper95, p_value,\n"
        "    #         cutoff, sensitivity, specificity, likelihood_ratio.\n"
        "    # NOTE: exploratory only — no independent validation cohort.\n"
        "    pass\n"
    )
    (RESULTS_DIR / "ROC_exploratory_helper.py").write_text(src)
    log("  Template saved to results/ROC_exploratory_helper.py")


# 20 Session info
def step20_session_info() -> None:
    log("STEP 20 — Session info + manifest")
    info = {
        "python_version":  sys.version,
        "platform":        platform.platform(),
        "numpy":           np.__version__,
        "pandas":          pd.__version__,
        "networkx":        nx.__version__,
        "matplotlib":      matplotlib.__version__,
        "scipy":           __import__("scipy").__version__,
        "sklearn":         __import__("sklearn").__version__,
        "requests":        requests.__version__,
        "statsmodels":     __import__("statsmodels").__version__,
    }
    with open(RESULTS_DIR / "python_sessionInfo.txt", "w") as fh:
        for k, v in info.items():
            fh.write(f"{k}: {v}\n")

    pngs = sorted(p.name for p in FIGURES_DIR.glob("*.png"))
    csvs = sorted(p.name for p in RESULTS_DIR.glob("*.csv"))
    txts = sorted(p.name for p in RESULTS_DIR.glob("*.txt"))
    others = sorted(p.name for p in RESULTS_DIR.iterdir()
                    if p.is_file() and p.suffix not in (".csv", ".txt"))
    with open(RESULTS_DIR / "MANIFEST.txt", "w") as fh:
        fh.write("# Figures\n")
        fh.writelines(f"- {p}\n" for p in pngs)
        fh.write("\n# Result tables (CSV)\n")
        fh.writelines(f"- {p}\n" for p in csvs)
        fh.write("\n# Text artifacts\n")
        fh.writelines(f"- {p}\n" for p in txts)
        fh.write("\n# Other\n")
        fh.writelines(f"- {p}\n" for p in others)
    log("  Pipeline finished.")


# Main
def main() -> None:
    print("=" * 78)
    print(" GSE58438 Renal IRI pipeline — Python continuation (steps 12–20)")
    print("=" * 78)

    centralities, edges_df, G, hub10, hub90 = step12_string_ppi()
    step13_hub_bar_plots(centralities)

    step14_hub10_subnetwork(hub10)
    step15_mirna_table()
    predict_venn()
    log("STEP 16 — Spearman template (no qPCR data; function defined)")
    rng = np.random.default_rng(SEED)
    demo_ct = pd.DataFrame({
        "BRCA1": rng.normal(0, 1, 12),
        "RAD51": rng.normal(0, 1, 12),
        "rno-miR-103-3p": -rng.normal(0, 1, 12),  # anti-correlated
        "rno-miR-107-3p": -rng.normal(0, 1, 12),
    })
    run_mirna_mrna_spearman(demo_ct,
                            gene_cols=["BRCA1", "RAD51"],
                            mirna_cols=["rno-miR-103-3p", "rno-miR-107-3p"])
    step17_gse9943_validation()
    step18_wgcna()
    step19_roc_template()
    step20_session_info()

    print("\n=== Python pipeline (steps 12-20) finished ===")
    print(f"PNG figures -> {FIGURES_DIR.resolve()}")
    print(f"CSV tables  -> {RESULTS_DIR.resolve()}")


if __name__ == "__main__":
    main()
