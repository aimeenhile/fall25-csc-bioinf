# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path
from typing import Tuple, Dict

REGIONS = ["VIS", "HIPP"]

def load_type_neurons(path: str | Path) -> pd.DataFrame:
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

    return type_neurons

def compute_spearman_corr(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Spearman correlation matrix"""

    corr, _ = spearmanr(df, axis=0, nan_policy='omit')
    CN = pd.DataFrame(corr, index=df.columns, columns=df.columns)

    return CN

def clean_names(CN: pd.DataFrame) -> pd.DataFrame:
    """Clean up row and column names"""

    # Rename "Hippocampus" -> "HIPP" & "VisCortex" -> "VIS"
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


