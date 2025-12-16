"""Aggregate PyPSA-Eur biomass potentials to country-level CSV for EnergyScope.

The script sums per-region biomass potential columns, converts them to GWh, and
stores the results keyed by ISO country code. Snakemake uses the output when
building EnergyScope resource assumptions.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import sys

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Aggregate biomass potentials by country
def aggregate_biomass_potentials(biomass_file: Path) -> pd.DataFrame:
    frame = pd.read_csv(biomass_file)
    id_column = "name" if "name" in frame.columns else frame.columns[0]
    frame = frame.rename(columns={id_column: "region"})
    if "geometry" in frame.columns:
        frame = frame.drop(columns="geometry")

    frame["region"] = frame["region"].astype(str).str.strip()
    frame["country"] = frame["region"].str.replace(" ", "", regex=False).str[:2]

    numeric_cols = [c for c in frame.columns if c not in {"region", "country"}]
    frame[numeric_cols] = frame[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    grouped = frame.groupby("country")[numeric_cols].sum()
    grouped.index.name = "country"
    return grouped * 1000.0 # potentials in GWh


def main(snakemake: "Snakemake") -> None:
    biomass_path = Path(snakemake.input.biomass)
    output_path = Path(snakemake.output.potentials)
    potentials = aggregate_biomass_potentials(biomass_path)
    ensure_parent(output_path)
    potentials.to_csv(output_path)


if __name__ == "__main__":
    main(snakemake)  
