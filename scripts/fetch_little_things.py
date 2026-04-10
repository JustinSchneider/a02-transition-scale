"""
Paper 3 Block 4: Fetch LITTLE THINGS rotation curves and galaxy metadata.

Downloads the Oh et al. (2015) catalog from VizieR (J/AJ/149/180) and saves
three tables as CSV files in data/raw/LITTLE_THINGS/:

  - galaxies.csv:  Galaxy properties and mass modeling results (~26 rows)
  - rotdmbar.csv:  Total rotation curves incl. baryonic contribution (~1700 rows)
  - rotdm.csv:     DM-only rotation curves, baryonic subtracted (~1700 rows)

Requires: astroquery, pandas

Usage:
    python scripts/fetch_little_things.py
"""

import os
import sys

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "raw", "LITTLE_THINGS")
CATALOG = "J/AJ/149/180"
TABLE_GALAXIES = f"{CATALOG}/galaxies"
TABLE_ROTCURVES = f"{CATALOG}/rotdmbar"
TABLE_ROTDM = f"{CATALOG}/rotdm"

# DDO 154 appears in SPARC, THINGS, and LITTLE THINGS
# VizieR uses underscores in names (e.g., "DDO_154")
DDO154_PATTERN = "DDO154"


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------
def fetch_vizier_tables():
    """Fetch galaxies, rotdmbar, and rotdm tables from VizieR."""
    try:
        from astroquery.vizier import Vizier
    except ImportError:
        print("ERROR: astroquery is not installed. Run: pip install astroquery")
        sys.exit(1)

    print(f"Fetching catalog {CATALOG} from VizieR...")
    print("(this may take a moment)\n")

    try:
        v = Vizier(row_limit=-1)
        catalogs = v.get_catalogs(CATALOG)
    except Exception as e:
        print(f"ERROR: Failed to fetch from VizieR: {e}")
        sys.exit(1)

    # Verify expected tables exist
    available_keys = list(catalogs.keys())
    print(f"Available tables: {available_keys}\n")

    for tbl_key in [TABLE_GALAXIES, TABLE_ROTCURVES, TABLE_ROTDM]:
        if tbl_key not in available_keys:
            print(f"ERROR: Expected table '{tbl_key}' not found.")
            print(f"Available: {available_keys}")
            sys.exit(1)

    # Convert astropy tables to pandas, handling masked columns
    df_meta = catalogs[TABLE_GALAXIES].filled().to_pandas()
    df_curves = catalogs[TABLE_ROTCURVES].filled().to_pandas()
    df_dm = catalogs[TABLE_ROTDM].filled().to_pandas()

    return df_meta, df_curves, df_dm


def save_tables(df_meta, df_curves, df_dm):
    """Write DataFrames to CSV in the output directory."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    meta_path = os.path.join(OUTPUT_DIR, "galaxies.csv")
    curves_path = os.path.join(OUTPUT_DIR, "rotdmbar.csv")
    dm_path = os.path.join(OUTPUT_DIR, "rotdm.csv")

    df_meta.to_csv(meta_path, index=False)
    df_curves.to_csv(curves_path, index=False)
    df_dm.to_csv(dm_path, index=False)

    print(f"Saved: {meta_path}")
    print(f"Saved: {curves_path}")
    print(f"Saved: {dm_path}\n")


def print_summary(df_meta, df_curves, df_dm):
    """Print summary statistics and flag cross-sample galaxies."""
    print("=" * 65)
    print("  LITTLE THINGS Data Summary (Oh et al. 2015)")
    print("=" * 65)

    # Galaxy metadata
    print(f"\n--- galaxies table ({len(df_meta)} rows) ---\n")
    print(f"  Columns: {list(df_meta.columns)}\n")

    # Rotation curves
    print(f"--- rotdmbar table ({len(df_curves)} rows) ---\n")
    print(f"  Columns: {list(df_curves.columns)}")

    # DM-only rotation curves
    print(f"\n--- rotdm table ({len(df_dm)} rows) ---\n")
    print(f"  Columns: {list(df_dm.columns)}")

    # Identify the galaxy name column
    name_col = None
    for candidate in ["Name", "Galaxy", "name", "galaxy"]:
        if candidate in df_meta.columns:
            name_col = candidate
            break

    if name_col is None:
        print("\n  WARNING: Could not identify galaxy name column in metadata.")
        print(f"  Available columns: {list(df_meta.columns)}")
        return

    galaxies = df_meta[name_col].astype(str).tolist()
    print(f"\n  Galaxy names ({name_col} column):")
    for g in galaxies:
        print(f"    {g}")

    # Count unique galaxies in rotation curve table
    curve_name_col = None
    for candidate in ["Name", "Galaxy", "name", "galaxy"]:
        if candidate in df_curves.columns:
            curve_name_col = candidate
            break

    if curve_name_col:
        n_curve_galaxies = df_curves[curve_name_col].nunique()
        print(f"\n  Unique galaxies in rotation curves: {n_curve_galaxies}")

    # Flag DDO 154
    print("\n" + "-" * 65)
    normalized = [g.strip().replace(" ", "").replace("_", "").upper() for g in galaxies]
    ddo154_found = any(n == "DDO154" for n in normalized)

    if ddo154_found:
        print("  *** DDO 154 FOUND — cross-sample galaxy (SPARC + THINGS + LITTLE THINGS)")
        if curve_name_col:
            ddo154_rows = df_curves[
                df_curves[curve_name_col].str.strip().str.replace(" ", "")
                .str.replace("_", "").str.upper() == "DDO154"
            ]
            print(f"  *** DDO 154 rotation curve points: {len(ddo154_rows)}")
    else:
        print("  DDO 154 not found in galaxy list (unexpected — check name formatting)")

    print("-" * 65)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 65)
    print("  Paper 3 Block 4: LITTLE THINGS Data Acquisition")
    print("  Source: Oh et al. (2015), AJ, 149, 180")
    print(f"  VizieR catalog: {CATALOG}")
    print("=" * 65)
    print()

    df_meta, df_curves, df_dm = fetch_vizier_tables()
    save_tables(df_meta, df_curves, df_dm)
    print_summary(df_meta, df_curves, df_dm)

    print("\nDone. Data saved to data/raw/LITTLE_THINGS/")


if __name__ == "__main__":
    main()
