# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, norm
from pathlib import Path
from typing import Tuple, Dict, List
import math

def load_type_data(path: str | Path) -> pd.DataFrame:
    """Load altExonUsage_devel_type.gz and preprocess fields"""

    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Filter column with "Neuron"
    neurons = [col for col in df.columns if "Neuron" in col]
    type_neurons = df[["Exon", "Gene"] + neurons].copy()

    # Combine Exon + Gene into GE (Exon::Gene)
    type_neurons["GE"] = type_neurons["Exon"] + "::" + type_neurons["Gene"]
    type_neurons = type_neurons.set_index("GE")
    type_neurons = type_neurons.drop(columns=["Exon", "Gene"])

    type_neurons = type_neurons.dropna(axis=0, how='any')

    return type_neurons

def compute_spearman_corr(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Spearman correlation matrix & Clean up row and column names"""
    
    corr, p = spearmanr(df, axis=0, nan_policy='omit')
    df = pd.DataFrame(corr, index=df.columns, columns=df.columns)

    return df

def preprocess_type(df: pd.DataFrame) -> pd.DataFrame:
    """Rename "Hippocampus" -> "HIPP" & "VisCortex" -> "VIS" """
    rename = (df.index
        .str.replace("Neuron", "", regex=False)
        .str.replace("Hippocampus", "HIPP", regex=False)
        .str.replace("VisCortex", "VIS", regex=False)
    )

    df.index = rename
    df.columns = rename

    return df

def spearman_corr_labels(CN: pd.DataFrame) -> pd.DataFrame:
    """Extract metadata to create labels for Fig.4a"""
    labels = pd.DataFrame({"CT": CN.index})
    labels[["Age", "Region", "Type"]] = labels["CT"].str.split("_", expand=True)

    return labels

def load_subtype_data(path: str | Path) -> pd.DataFrame:
    """Load altExonUsage_devel_subtype.gz and preprocess fields"""

    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Filter column with "Inh","Excite","Granule","NIPC"
    patterns = ["Inh", "Excite", "Granule", "NIPC"]
    df_cols = [c for c in df.columns if any(p in c for p in patterns)]
    subtype = df[["Exon", "Gene"] + df_cols].copy()

    # Combine Exon + Gene into GE (Exon::Gene)
    subtype["GE"] = subtype["Exon"] + "::" + subtype["Gene"]
    subtype = subtype.set_index("GE")
    subtype = subtype.drop(columns=["Exon", "Gene"])

    threshold = 0.95
    subtype = subtype.loc[:, subtype.isna().mean() < threshold]  # filter columns
    subtype = subtype.loc[subtype.isna().mean(axis=1) <= threshold, :]  # filter rows

    return subtype

def preprocess_subtype(df: pd.DataFrame) -> pd.DataFrame:
    """Extract upper triangule matrix and assign subtype labels"""
    
    def get_upper_tri(patterns: List[str]) -> np.array:
        """Get upper triangle matrix for selected patterns"""
        cols = [c for c in df.columns if any(p in c for p in patterns)]
        rows = [r for r in df.index if any(p in r for p in patterns)]
        ut_mat = df.loc[rows, cols]

        return ut_mat.values[np.triu_indices_from(ut_mat, k =1)]
    
    # Separate out Excite, Inhib, and Inh
    CE = get_upper_tri(["Excite", "Granule"])
    CR = get_upper_tri(["Inh"])
    CI = get_upper_tri(["Inhib"])

    # Combine into one dataframe
    df = pd.DataFrame({
        "CorValue": np.concatenate([CE, CR, CI]),
        "Type": ["Excite"]*len(CE) + ["Inh_wCR"]*len(CR) + ["Inhib"]*len(CI)
    })

    # Set order
    df["Type"] = pd.Categorical(df["Type"], categories=["Excite","Inh_wCR","Inhib"], ordered=True)

    return df

def compute_spearman_ci(x, y):
    """Compute Spearman correlation and confidence interval"""
    alpha = 0.05
    r, p = spearmanr(x, y, nan_policy="omit")
    n = np.sum(~(np.isnan(x) | np.isnan(y)))

    if n <= 3:
        return r, np.nan, np.nan
    
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf(1 - alpha/2)
    low = np.tanh(z - z_crit * se)
    high = ow = np.tanh(z + z_crit * se)

    return r, low, high



