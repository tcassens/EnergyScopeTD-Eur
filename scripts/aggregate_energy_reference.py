"""Aggregate per-cluster energy reference tables to country totals.

PyPSA-Eur exports per-cluster energy references, this script sums them across
clusters using ISO mappings so EnergyScope can split demand by sector. Both
general energy totals and industrial-only totals are joined into one CSV.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import sys

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent, is_valid_country, country_key

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Aggregate energy reference tables by country
def aggregate_by_country(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        return pd.DataFrame()
    df.index = df.index.astype(str).str.strip()
    numeric = df.apply(pd.to_numeric, errors="coerce")
    codes = numeric.index.map(country_key)
    mask = codes.map(is_valid_country)
    aggregated = numeric.loc[mask].groupby(codes[mask]).sum()
    aggregated.index.name = "country"
    return aggregated


def main(snakemake: "Snakemake") -> None:
    energy_totals = Path(snakemake.input.energy_totals)
    industrial_totals = Path(snakemake.input.industrial_totals)
    target = Path(snakemake.output.csv)

    energy = aggregate_by_country(energy_totals)
    industrial = aggregate_by_country(industrial_totals)

    combined = energy.join(industrial, how="outer", lsuffix="_energy", rsuffix="_industrial").fillna(0.0)
    ensure_parent(target)
    combined.to_csv(target)


if __name__ == "__main__":
    main(snakemake)
