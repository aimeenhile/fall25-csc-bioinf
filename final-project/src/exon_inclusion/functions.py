# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
import scipy as sp
from pathlib import Path
from typing import Tuple, Dict

REGIONS = ["VIS", "HIPP"]

def load_exon_usage(path: str | Path):
    """Load altExonUsage_devel_type.gz and preprocess fields"""

    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Filter only neuron columns
    neuron_cols = [c for c in df.columns if 'Neuron' in c]
    df_neurons = df[['Exon', 'Gene'] + neuron_cols].copy()

    # Combine Exon and Gene into a single row identifier
    df_neurons['Exon_Gene'] = df_neurons['Exon'] + '::' + df_neurons['Gene']
    df_neurons = df_neurons.drop(columns=['Exon', 'Gene'])

    # Set combined identifier as index
    df_neurons.set_index('Exon_Gene', inplace=True)

    # Clean column names
    new_cols = df_neurons.columns.str.replace('Neuron', '', regex=False)
    new_cols = new_cols.str.replace('Hippocampus', 'HIPP', regex=False)
    new_cols = new_cols.str.replace('VisCortex', 'VIS', regex=False)
    df_neurons.columns = new_cols

    return df_neurons