"""Generate Fig. 8: DDO 154 three-pipeline RT fit comparison.

Single panel showing observed rotation curves (data points with error bars)
and RT model fits from SPARC, THINGS, and LITTLE THINGS, with vertical
lines marking each pipeline's fitted taper radius Rt.

Output: manuscript/figures/fig08_ddo154_three_pipeline.{pdf,png}
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Project imports
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.database import get_engine, get_session, query_profiles_as_dataframe
from src.physics import rt_model_velocity

# ---------------------------------------------------------------------------
# Publication style (matches NB18)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# Colorblind-safe palette (from NB18)
C_SPARC = "#999999"
C_THINGS = "#0072B2"
C_LT = "#D55E00"

# ---------------------------------------------------------------------------
# Load DDO 154 radial profiles from all three pipelines
# ---------------------------------------------------------------------------
engine = get_engine()
session = get_session(engine)

sparc = query_profiles_as_dataframe(session, "DDO154")
things = query_profiles_as_dataframe(session, "DDO154_THINGS")
lt = query_profiles_as_dataframe(session, "DDO_154_LITTLE_THINGS")

session.close()

# ---------------------------------------------------------------------------
# Best-fit RT parameters (from results CSVs, verified in plan)
# ---------------------------------------------------------------------------
pipelines = [
    {
        "name": "SPARC",
        "data": sparc,
        "omega": 15.524704,
        "Rt": 3.250104,
        "color": C_SPARC,
        "marker": "o",
        "zorder_data": 1,
        "zorder_fit": 2,
    },
    {
        "name": "THINGS",
        "data": things,
        "omega": 12.094323,
        "Rt": 4.821186,
        "color": C_THINGS,
        "marker": "s",
        "zorder_data": 3,
        "zorder_fit": 4,
    },
    {
        "name": "LITTLE THINGS",
        "data": lt,
        "omega": 34.163320,
        "Rt": 1.104753,
        "color": C_LT,
        "marker": "^",
        "zorder_data": 5,
        "zorder_fit": 6,
    },
]

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))

for p in pipelines:
    df = p["data"]
    r = df["radius_kpc"].values
    v_obs = df["v_obs"].values
    v_err = df["v_err"].values
    v_bary = df["v_baryon_total"].values

    N = len(r)

    # Observed data (error bars)
    ax.errorbar(
        r, v_obs, yerr=v_err,
        fmt=p["marker"], ms=3, color=p["color"],
        alpha=0.45, elinewidth=0.6, capsize=0,
        label=f'{p["name"]} ($N={N}$)',
        zorder=p["zorder_data"],
    )

    # RT model fit (smooth curve evaluated on dense grid)
    r_fine = np.linspace(r.min(), r.max(), 300)
    v_bary_fine = np.interp(r_fine, r, v_bary)
    v_model_fine = rt_model_velocity(r_fine, v_bary_fine, p["omega"], p["Rt"])
    ax.plot(
        r_fine, v_model_fine, "-",
        color=p["color"], linewidth=1.8,
        zorder=p["zorder_fit"],
    )

    # Taper radius Rt (vertical line)
    ax.axvline(
        p["Rt"], color=p["color"], linestyle=":", linewidth=1.0,
        alpha=0.7, zorder=p["zorder_fit"] - 1,
    )

ax.set_xlabel(r"$R$ [kpc]")
ax.set_ylabel(r"$V$ [km s$^{-1}$]")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend(loc="lower right")

plt.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig_dir = ROOT / "manuscript" / "figures"
fig_dir.mkdir(parents=True, exist_ok=True)

name = "fig08_ddo154_three_pipeline"
fig.savefig(fig_dir / f"{name}.pdf")
fig.savefig(fig_dir / f"{name}.png")
print(f"Saved: {name}.pdf, {name}.png")
