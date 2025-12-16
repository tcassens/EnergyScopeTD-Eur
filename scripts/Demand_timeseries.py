"""Generate normalized demand time series per country and carrier for EnergyScope.

This script reads load profiles from a PyPSA-Eur network, aggregates them to
country level, normalizes each carrier profile, and pairs the result with
transport demand intensities derived from the mobility dataset. One CSV per
country plus a totals file are written for the EnergyScope conversion workflow.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Dict, Tuple

import sys

import numpy as np
import pandas as pd
import xarray as xr

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent, is_valid_country, country_key

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Normalize a pandas Series to sum to 1, handling edge cases
def normalized(series: pd.Series) -> pd.Series:
    cleaned = series.astype(float).clip(lower=0.0)
    total = cleaned.sum()
    if total <= 0 or not np.isfinite(total):
        return pd.Series(0.0, index=cleaned.index)
    return cleaned / total

# Aggregate transport data by country
def aggregate_transport(transport_df: pd.DataFrame) -> pd.DataFrame:
    renamed = transport_df.copy()
    renamed.columns = renamed.columns.map(country_key)
    renamed = renamed.loc[:, renamed.columns.map(is_valid_country)]
    aggregated = renamed.T.groupby(level=0).sum().T
    aggregated.index = pd.RangeIndex(len(aggregated))
    return aggregated

# Load country-specific demand profiles from PyPSA-Eur dataset
def load_country_profiles(
    ds: xr.Dataset,
) -> Tuple[pd.Index, Dict[str, pd.DataFrame], Dict[str, pd.Series]]:
    loads_ts = ds["loads_t_p_set"].to_pandas()
    loads_ts.columns = loads_ts.columns.astype(str).str.strip()

    all_load_ids = pd.Index(map(str, ds.indexes["loads_i"])).str.strip()

    if "loads_p_set" in ds.data_vars:
        static_vals = ds["loads_p_set"].to_pandas().astype(float)
        static_vals.index = static_vals.index.astype(str).str.strip()
    else:
        static_vals = pd.Series(0.0, index=all_load_ids)

    missing = all_load_ids.difference(loads_ts.columns)
    if not missing.empty:
        static_panel = pd.DataFrame(
            {name: static_vals[name] for name in missing},
            index=loads_ts.index,
        )
        loads_panel = pd.concat([loads_ts, static_panel], axis=1)
    else:
        loads_panel = loads_ts

    load_to_carrier = ds["loads_carrier"].to_pandas().astype(str)
    load_to_carrier.index = load_to_carrier.index.astype(str).str.strip()

    load_to_bus = ds["loads_bus"].to_pandas().astype(str)
    load_to_bus.index = load_to_bus.index.astype(str).str.strip()

    bus_to_country = ds["buses_country"].to_pandas().astype(str).str.strip()
    bus_to_country.index = bus_to_country.index.astype(str).str.strip()

    load_to_country = load_to_bus.map(bus_to_country)

    k_country = load_to_country.reindex(loads_panel.columns)
    k_carrier = load_to_carrier.reindex(loads_panel.columns)
    valid = k_country.notna() & k_carrier.notna()

    grouped = (
        loads_panel.loc[:, valid]
        .T.groupby([k_country[valid], k_carrier[valid]])
        .sum()
        .T
    )

    snapshots = pd.Index(ds.indexes["snapshots"])
    country_profiles: Dict[str, pd.DataFrame] = {}
    country_totals: Dict[str, pd.Series] = {}

    for country in grouped.columns.get_level_values(0).unique():
        if not is_valid_country(country):
            continue
        panel = grouped.xs(country, axis=1, level=0).fillna(0.0)
        panel.index = snapshots
        normalized_panel = panel.apply(normalized, axis=0)
        country_profiles[country] = normalized_panel
        country_totals[country] = panel.sum(axis=0)

    return snapshots, country_profiles, country_totals

# Write country-specific CSV files including transport demand
def write_country_csvs(
    snapshots: pd.Index,
    country_profiles: Dict[str, pd.DataFrame],
    country_totals: Dict[str, pd.Series],
    transport_csv: Path,
    output_dir: Path,
) -> None:
    ensure_parent(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    transport = aggregate_transport(pd.read_csv(transport_csv, index_col=0).astype(float))

    totals_frame = pd.DataFrame.from_dict(country_totals, orient="index")
    totals_frame.index.name = "country"

    all_countries = sorted(set(country_profiles.keys()) | set(transport.columns))
    totals_frame = (
        pd.DataFrame(index=all_countries)
        .join(totals_frame, how="left")
        .fillna(0.0)
    )


    for country, profile in sorted(country_profiles.items()):
        zero = pd.Series(0.0, index=snapshots)
        transport_ts = transport[country] if country in transport.columns else zero
        transport_profile = normalized(transport_ts)

        result = profile.copy()
        result["passenger_mobility"] = transport_profile
        result["freight_mobility"] = transport_profile

        out_file = output_dir / f"Demand_timeseries_{country}.csv"
        ensure_parent(out_file)
        result.to_csv(out_file)

    totals_path = output_dir / "Demand_timeseries_totals.csv"
    ensure_parent(totals_path)
    totals_frame.to_csv(totals_path)


def main(snakemake: "Snakemake") -> None:
    network_path = Path(snakemake.input.network)
    transport_csv = Path(snakemake.input.transport_csv)
    output_dir = Path(snakemake.output.Demand_timeseries_dir)

    with xr.open_dataset(network_path) as ds:
        snapshots, country_profiles, country_totals = load_country_profiles(ds)

    write_country_csvs(
        snapshots,
        country_profiles,
        country_totals,
        transport_csv,
        output_dir,
    )


if __name__ == "__main__":
    main(snakemake)
