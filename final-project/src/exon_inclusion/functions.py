# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path
from typing import Tuple, Dict, List

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
    CN = pd.DataFrame(corr, index=df.columns, columns=df.columns)

    return CN

def preprocess_type(CN: pd.DataFrame) -> pd.DataFrame:
    """Rename "Hippocampus" -> "HIPP" & "VisCortex" -> "VIS" """
    rename = (CN.index
        .str.replace("Neuron", "", regex=False)
        .str.replace("Hippocampus", "HIPP", regex=False)
        .str.replace("VisCortex", "VIS", regex=False)
    )

    CN.index = rename
    CN.columns = rename

    return CN

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

    return subtype

def preprocess_subtype(df: pd.DataFrame) -> pd.DataFrame:
    """ """
    threshold = 0.95

    # Filter missing data
    df = df.loc[:, df.isna().mean() <= threshold]
    df = df.loc[df.isna().mean(axis=1) <= threshold, :]

    return df