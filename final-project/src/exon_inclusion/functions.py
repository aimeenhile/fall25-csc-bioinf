# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
import scipy as sp
from pathlib import Path
from typing import Tuple, Dict

REGIONS = ["VIS", "HIPP"]

def load_type_neurons(path: str | Path):
    """Load altExonUsage_devel_type.gz and preprocess fields"""

    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Filter column with "Neuron"
    neurons = [col for col in df.columns if "Neuron" in col]
    type_neurons = df[["Exon", "Gene"] + neurons].copy()

    # Combine Exon + Gene into GE
    type_neurons["GE"] = type_neurons["Exon"] + "::" + type_neurons["Gene"]
    type_neurons = type_neurons.set_index("GE")
    type_neurons = type_neurons.drop(columns=["Exon", "Gene"])

    return type_neurons

def get_correlation():
    """ """