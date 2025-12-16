"""Convert PyPSA-Eur hydro inflow series to country-level CSV.

Snakemake runs this lightweight helper to transform the hydro reservoir
inflow profile from PyPSA-Eur into a table indexed by time, with columns
named after country codes. The resulting CSV is saved in the intermediate
folder.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import sys

import pandas as pd
import xarray as xr

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent

if TYPE_CHECKING:
    from snakemake.script import Snakemake


def main(snakemake: "Snakemake") -> None:
    inflow = xr.open_dataarray(snakemake.input.profile)
    df = inflow.to_pandas()
    df.index.name = "time"
    df.columns = pd.Index(df.columns).astype(str).str.upper()

    out_path = Path(snakemake.output.csv)
    ensure_parent(out_path)
    df.to_csv(out_path)


if __name__ == "__main__":
    main(snakemake)
