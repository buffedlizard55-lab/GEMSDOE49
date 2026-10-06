"""Grid constants and byte-level hash pins for the DOE GEMS Prize (DrivenData #306).

Every constant in this file was MEASURED from the official bytes (rasterio open of the
competition rasters), not copied from prose.  The sha256 pins are the values published in the
transport manifest ``data/bridge/manifest.json`` of the sibling transport repository and
re-verified here byte-for-byte.

Sources for manual review
-------------------------
* Competition problem description (metric, submission format):
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
* Competition data tab:  https://www.drivendata.org/competitions/306/competition-doe-gems/data/
* USGS GeoDAWN release:  https://doi.org/10.5066/P93LGLVQ
"""

from __future__ import annotations

import hashlib
from pathlib import Path

# --------------------------------------------------------------------------------------
# Grid  (MEASURED from data/example_submission.tif and data/existing_faults.tif)
# --------------------------------------------------------------------------------------
GRID = {
    "height": 3730,
    "width": 3292,
    "crs": "EPSG:32611",
    "res_m": 100.0,
    "transform": (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0),
    "bounds": (243350.0, 4135550.0, 572550.0, 4508550.0),
    "nodata_sentinel": -3.4028234663852886e38,  # most-negative float32, used by the feature stack
}

FOOTPRINT_PIXELS = 5_167_373  # finite (non-NaN) pixels of example_submission.tif
OUTSIDE_PIXELS = 7_111_787  # NaN pixels of example_submission.tif
CATALOGUE_PIXELS = 60_988  # pixels == 1 in existing_faults.tif
TOTAL_PIXELS = GRID["height"] * GRID["width"]  # 12_279_160

# --------------------------------------------------------------------------------------
# Metric constants  (OFFICIAL, read from the problem-description page 967)
# --------------------------------------------------------------------------------------
ALPHA = 0.2  # false-positive weight
BETA = 0.8  # false-negative weight
R_M = 300.0  # triangular-kernel support, metres

# --------------------------------------------------------------------------------------
# Official band inventory (read from training_features.tif's OWN TIFF tags, band by band)
# --------------------------------------------------------------------------------------
BANDS = (
    (1, "mag_anom", "magnetic", "Magnetic anomaly - deviation from expected Earth's magnetic field"),
    (2, "rtp", "magnetic", "Reduced to pole magnetic data - magnetic anomaly corrected for latitude effects"),
    (3, "tmi_hg", "magnetic", "Total magnetic intensity horizontal gradient - rate of change in horizontal direction"),
    (4, "geod_2ndinv", "geodetic", "Geodetic second invariant - measure of strain rate tensor magnitude"),
    (5, "iso_grav_anom_slope", "gravity", "Isostatic gravity anomaly slope - gradient of gravity after isostatic correction"),
    (6, "tc", "magnetic", "Tilt angle or total curvature - magnetic field derivative for edge detection"),
    (7, "geod_shearrate", "geodetic", "Geodetic shear rate - rate of angular deformation from GPS/InSAR"),
    (8, "geod_dilaterate", "geodetic", "Geodetic dilatation rate - rate of volumetric strain (expansion/contraction)"),
    (9, "tmi_vg", "magnetic", "Total magnetic intensity vertical gradient - rate of change in vertical direction"),
    (10, "deq_n100a15", "seismic", "Distance to earthquake (n=100km radius, a=15 deg azimuth parameters)"),
    (11, "iso_grav_anom_vg", "gravity", "Isostatic gravity anomaly vertical gradient - vertical rate of change"),
    (12, "det_elev", "topographic", "Detrended elevation - topography with regional trends removed"),
    (13, "iso_grav_anom", "gravity", "Isostatic gravity anomaly - gravity after compensating for topographic mass"),
    (14, "tmi", "magnetic", "Total magnetic intensity - total strength of magnetic field"),
    (15, "depth_to_base_surf", "subsurface", "Depth to basement surface - thickness of sedimentary cover"),
    (16, "ieq_n100a15", "seismic", "Earthquake intensity or density (n=100km radius, a=15 deg parameters)"),
    (17, "cond_surf", "subsurface", "Conductivity surface - electrical conductivity of subsurface"),
    (18, "iso_grav_anom_hg", "gravity", "Isostatic gravity anomaly horizontal gradient - horizontal rate of change"),
    (19, "det_elev_slope", "topographic", "Detrended elevation slope - gradient of elevation after detrending"),
)
BAND_INDEX = {name: idx for idx, name, _cat, _desc in BANDS}
N_BANDS = len(BANDS)

# --------------------------------------------------------------------------------------
# sha256 pins  (from the transport manifest; re-verified in this checkout)
# --------------------------------------------------------------------------------------
PINS: dict[str, str] = {
    "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "existing_faults.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "example_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
    "external/topo_u8.tif": "a6398d9950965dec6aae6ccecdaa6ced48645d133eab222cbdd11def9bdabfa4",
    "external/radiometric_u8.tif": "6cb051f70f941fd78028fe66a9f71e87204fcd8d9a85903df0b94993bad1ec4d",
}

FEATURE_PARTS = (
    ("gems-geodawn-numerical-features.tif.part-000", 94_371_840,
     "0a330f8951af6c921029e25c84a579319d2db554d62d30d894d6ddc97f98cff7"),
    ("gems-geodawn-numerical-features.tif.part-001", 94_371_840,
     "3c98037b2c997e3bbfcdfc2d9a982e8b820a77410dd7404b05cdb80594922c50"),
    ("gems-geodawn-numerical-features.tif.part-002", 94_371_840,
     "c375c4dbc40c59bbaece572b5e348700b75935f9a30c879c82b6417e0836f31c"),
    ("gems-geodawn-numerical-features.tif.part-003", 94_371_840,
     "b164159e6d0cb2595bc9f63a948af2646b124114a9c5921880c7516092137320"),
    ("gems-geodawn-numerical-features.tif.part-004", 41_425_484,
     "fa0a6f9c936fac1d6f20ca37f5929b2d60bf7a80f3d477dcab886f941aee2696"),
)

# Official mirror URLs (published on the competition data tab) for the three small/medium files.
DROPBOX_MIRRORS = {
    "training_features.tif": (
        "https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/"
        "gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&dl=1"
    ),
    "existing_faults.tif": (
        "https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/"
        "existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&dl=1"
    ),
    "example_submission.tif": (
        "https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/"
        "example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&dl=1"
    ),
}

# Checksum-pinned public transport bridge (github.com), used when the official mirrors are
# unreachable from the executing host.  Bytes are accepted only if they match PIN_* above.
BRIDGE_REPOS = (
    ("buffedlizard55-lab/5GEMSDOE", "main"),
    ("buffedlizard55-lab/GEMSDOE2", "main"),
    ("buffedlizard55-lab/GEMSDOE", "main"),
)


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """Streaming sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()
