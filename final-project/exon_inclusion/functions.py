# exon_inclusion/functions.py

from __future__ import annotations
import pandas as pd
import numpy as np
import scipy as sp
from pathlib import Path
from typing import Tuple, Dict

def load_exon_usage(path: str | Path) -> pd.DataFrame:
    """Load altExonUsage_devel_type.gz and preprocess fields"""

    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Split Exon into chr + start + end + strand 
