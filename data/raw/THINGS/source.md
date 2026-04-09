# THINGS Data Source

## Rotation Curves (Curves/)

**Primary source:**
de Blok, W. J. G. et al. (2008). "High-Resolution Rotation Curves and Galaxy
Mass Models from THINGS." AJ, 136, 2648.
https://doi.org/10.1088/0004-6256/136/6/2648

19 galaxies with tilted-ring rotation curves derived from THINGS HI velocity
fields. Files are named `{GALAXY}.curve.02`.

**Column format** (from README.TXT):

| Col | Content                              | Units  | Use   |
|-----|--------------------------------------|--------|-------|
| 1   | Radius                               | arcsec | Yes   |
| 2   | Rotation velocity                    | km/s   | Yes   |
| 3   | Scatter in V_rot along ring          | km/s   | No    |
| 4   | Formal error in tilted-ring fit      | km/s   | No    |
| 5   | V_rot from approaching side          | km/s   | No    |
| 6   | V_rot from receding side             | km/s   | No    |
| 7   | Combined error (preferred)           | km/s   | Yes   |
| 8   | Position angle                       | deg    | Note  |
| 9   | Inclination                          | deg    | Note  |
| 10  | Systemic velocity                    | km/s   | Note  |

For RT fitting, use columns 1, 2, 7. Radius must be converted from arcsec to
kpc using the galaxy distance.

## Mass Models (MassModels/)

**Primary source:** Same as above (de Blok et al. 2008, Tables 4-7).

ROTMAS output files containing baryonic velocity decomposition. Files are named
`{GALAXY}.{HALO}.{ML}.REV.dat` where:
- HALO = ISO (pseudo-isothermal) or NFW
- ML = fix (fixed M/L) or free (fitted M/L)
- `.Kr.dat` suffix = Kroupa IMF; no suffix = diet-Salpeter IMF

**Column format** (from file headers):

| Col | Content  | Units | Use  |
|-----|----------|-------|------|
| 1   | Radius   | kpc   | Yes  |
| 2   | V_gas    | km/s  | Yes  |
| 3   | V_disk   | km/s  | Yes  |
| 4   | V_bulge  | km/s  | Yes  |
| 5   | V_obs    | km/s  | Yes  |
| 6   | err V_obs| km/s  | Yes  |
| 7   | V_halo   | km/s  | No   |
| 8   | V_total  | km/s  | No   |

The `.fix.REV.dat` files provide V_disk and V_bulge at M/L = 1.0 (same
convention as SPARC). Our analysis applies Y_disk = 0.5 and Y_bulge = 0.7
identically to the SPARC pipeline.

**Coverage note:** NGC3627 and NGC4826 have rotation curves but no mass model
files. These two galaxies are THINGS-only (not in SPARC), so baryonic
decomposition must be derived independently for them in Block 2.

## Parent Survey

Walter, F. et al. (2008). "THINGS: The HI Nearby Galaxy Survey." AJ, 136, 2563.
https://doi.org/10.1088/0004-6256/136/6/2563

34 galaxies observed with the VLA in HI 21-cm. de Blok et al. (2008) derived
rotation curves for 19 of these.

## Data Obtained Via

Garcia-Arroyo, G. et al. (2026). "Data-Driven Modeling of Rotation Curves with
Artificial Neural Networks." Physics of the Dark Universe, 52, 102240.
https://doi.org/10.1016/j.dark.2026.102240 (arXiv: 2404.05833)

The rotation curve and mass model files were originally provided by
Prof. W. J. G. de Blok and redistributed in the supplementary data of
Garcia-Arroyo et al. (2026).

## Galaxy Inventory (19 galaxies)

| Galaxy   | RC | MassModel | In SPARC | Notes                          |
|----------|----|-----------|----------|--------------------------------|
| DDO154   | Y  | Y         | Yes      |                                |
| IC2574   | Y  | Y         | Yes      |                                |
| NGC925   | Y  | Y         | No       | THINGS-only                    |
| NGC2366  | Y  | Y         | Yes      |                                |
| NGC2403  | Y  | Y         | Yes      | Has bulge; also has 1disk variant |
| NGC2841  | Y  | Y         | Yes      | Has bulge                      |
| NGC2903  | Y  | Y         | Yes      | Has bulge; also has allcurve variant |
| NGC2976  | Y  | Y         | Yes      |                                |
| NGC3031  | Y  | Y         | No       | THINGS-only; has bulge         |
| NGC3198  | Y  | Y         | Yes      | Has bulge; also has 1disk variant |
| NGC3521  | Y  | Y         | Yes      |                                |
| NGC3621  | Y  | Y         | No       | THINGS-only                    |
| NGC3627  | Y  | **N**     | No       | THINGS-only; needs baryonic decomp |
| NGC4736  | Y  | Y         | No       | THINGS-only; has bulge         |
| NGC4826  | Y  | **N**     | No       | THINGS-only; needs baryonic decomp |
| NGC5055  | Y  | Y         | Yes      | Has bulge                      |
| NGC6946  | Y  | Y         | Yes      | Has bulge                      |
| NGC7331  | Y  | Y         | Yes      | Has bulge                      |
| NGC7793  | Y  | Y         | Yes      | Also has short variant         |
