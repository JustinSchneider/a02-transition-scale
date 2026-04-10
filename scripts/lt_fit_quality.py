#!/usr/bin/env python
"""Compute fit-quality summary for resolved LITTLE THINGS galaxies.

Reads the NB12 fit results CSV and prints median reduced chi-squared
for the resolved subsample (excluding DDO 50, which pegs at Rt lower bound).
"""

import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path(__file__).resolve().parent.parent / "results"
df = pd.read_csv(results_dir / "NB12_lt_rt_fits.csv")

# Resolved galaxies, excluding DDO 50 (Rt pegs at lower bound)
resolved = df[(df["resolved"] == True) & (~df["galaxy_id"].str.contains("DDO_50"))]

chi2r = resolved["chi2_r"]
print(f"Resolved LITTLE THINGS (excl. DDO 50): N = {len(resolved)}")
print(f"  Median chi2_r = {np.median(chi2r):.2f}")
print(f"  Range: {chi2r.min():.2f} to {chi2r.max():.2f}")
print(f"  Values: {sorted(chi2r.round(2).tolist())}")
