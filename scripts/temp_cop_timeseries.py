"""Build country-aggregated temperature and COP time series.

Snakemake calls this script to produce population-weighted soil/air temperature
profiles and coefficient-of-performance (COP) profiles for each heat system /
heat source combination. The input data are outputs from PyPSA-Eur,
combined with population layouts to aggregate clusters into country-level time
series. The resulting CSVs are saved in the intermediate folder.
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


def aggregate_population_weighted(
    dataarray: xr.DataArray, pop_layout: pd.DataFrame
) -> pd.DataFrame:
    node_population = pop_layout["total"].astype(float)
    node_countries = pop_layout["ct"].astype(str)
    df = dataarray.to_pandas().reindex(columns=node_population.index)
    weighted = df.mul(node_population, axis=1)
    country_weighted = weighted.T.groupby(node_countries).sum().T
    country_population = node_population.groupby(node_countries).sum()
    aggregated = country_weighted.divide(country_population, axis=1)
    aggregated.index.name = "time"
    return aggregated


def slugify(value: str) -> str:
    return (
        str(value)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
    )


def export_cop_profiles(
    cop_profiles: xr.DataArray,
    pop_layout: pd.DataFrame,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for heat_system in cop_profiles.coords["heat_system"].values:
        for heat_source in cop_profiles.coords["heat_source"].values:
            subset = cop_profiles.sel(
                heat_system=heat_system, heat_source=heat_source
            ).transpose("time", "name")
            aggregated = aggregate_population_weighted(subset, pop_layout)
            filename = (
                f"{slugify(heat_system)}_{slugify(heat_source)}_cop_time_series.csv"
            )
            aggregated.to_csv(out_dir / filename)


def main(snakemake: "Snakemake") -> None:
    temp_soil = xr.open_dataarray(snakemake.input.temp_soil_nc)
    temp_air = xr.open_dataarray(snakemake.input.temp_air_nc)
    cop_profiles = xr.open_dataarray(snakemake.input.cop_profiles)
    pop_layout = pd.read_csv(snakemake.input.pop_csv, index_col=0)
    out_path_soil = Path(snakemake.output.soil_temp_time_series)
    out_path_air = Path(snakemake.output.air_temp_time_series)
    cop_dir = Path(snakemake.output.cop_time_series_dir)

    soil_temp_weighted = aggregate_population_weighted(temp_soil, pop_layout)
    air_temp_weighted = aggregate_population_weighted(temp_air, pop_layout)
    ensure_parent(out_path_soil)
    ensure_parent(out_path_air)
    soil_temp_weighted.to_csv(out_path_soil)
    air_temp_weighted.to_csv(out_path_air)

    export_cop_profiles(cop_profiles, pop_layout, cop_dir)


if __name__ == "__main__":
    main(snakemake)
