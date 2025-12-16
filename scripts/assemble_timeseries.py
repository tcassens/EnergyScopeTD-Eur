"""Combine demand, mobility and RET profiles into EnergyScope-ready time series.

For each country available in the demand and RET directories, the script
computes weighted electricity/heat profiles, attaches passenger/freight demand
shares, aggregates renewable capacity factors into the five technology groups,
and writes the consolidated `time_series_<country>.csv` file expected by
EnergyScope typical days selction.
"""


from collections.abc import Mapping as MappingABC
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, Mapping

import sys

import numpy as np
import pandas as pd
import xarray as xr

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import (
    ensure_parent,
    is_valid_country,
    load_indexed_csv,
)

if TYPE_CHECKING:
    from snakemake.script import Snakemake


TECH_GROUPS: Mapping[str, Iterable[str]] = {
    "PV": ["solar", "solar rooftop"],
    "Wind_onshore": ["onwind"],
    "Wind_offshore": ["offwind-ac", "offwind-dc", "offwind-float"],
    "Hydro_river": ["ror"],
    "Solar": [
        "rural solar thermal",
        "urban central solar thermal",
        "urban decentral solar thermal",
    ],
}

ELECTRICITY_COMPONENTS = [
    "electricity",
    "industry electricity",
    "agriculture electricity",
]

HEAT_COMPONENTS = [
    "rural heat",
    "urban central heat",
    "urban decentral heat",
]

TIMESERIES_COLUMNS = [
    "{PERIODS}",
    "Electricity (%_elec)",
    "Space Heating (%_sh)",
    "Passanger mobility (%_pass)",
    "Freight mobility (%_freight)",
    "PV",
    "Wind_onshore",
    "Wind_offshore",
    "Hydro_river",
    "Solar",
]


def weighted_mean_profile(
    profile: pd.DataFrame,
    totals: pd.Series,
    columns: Iterable[str],
) -> pd.Series:
    index = profile.index
    available = [col for col in columns if col in profile.columns]
    if not available:
        return pd.Series(0.0, index=index)

    weighted_sum = pd.Series(0.0, index=index)
    total_weight = 0.0
    for column in available:
        weight = float(totals.get(column, 0.0))
        if weight <= 0 or not np.isfinite(weight):
            continue
        weighted_sum = weighted_sum.add(profile[column] * weight, fill_value=0.0)
        total_weight += weight

    if total_weight <= 0:
        return profile[available].mean(axis=1).fillna(0.0)

    combined = weighted_sum / total_weight
    return combined.fillna(0.0)


def load_demand_profile(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "snapshots" in df.columns:
        df = df.set_index("snapshots")
    df.index = pd.Index(df.index).astype(int)
    return df


def load_ret_profile(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "snapshot" in df.columns:
        df = df.set_index("snapshot")
    df.index = pd.Index(df.index).astype(int)
    return df


def country_from_path(path: Path) -> str:
    stem = path.stem
    parts = stem.split("_")
    if not parts:
        return ""
    return parts[-1].upper()


def build_country_map(paths) -> Dict[str, Path]:
    if isinstance(paths, MappingABC):
        iterator = paths.items()
    else:
        iterator = ((None, entry) for entry in paths)
    mapping: Dict[str, Path] = {}
    for key, raw in iterator:
        path = Path(raw)
        if key is None or str(key).strip() == "":
            country = country_from_path(path)
        else:
            country = str(key).strip().upper()
        if not is_valid_country(country):
            continue
        mapping[country] = path
    return mapping


def country_carrier_capacities(ds: xr.Dataset) -> Dict[str, pd.Series]:
    carrier = ds["generators_carrier"].to_pandas().astype(str).str.strip()
    carrier.index = carrier.index.astype(str).str.strip()

    p_nom = ds["generators_p_nom"].to_pandas().astype(float)
    p_nom.index = p_nom.index.astype(str).str.strip()

    if "generators_p_nom_max" in ds:
        p_nom_max = ds["generators_p_nom_max"].to_pandas().astype(float)
        p_nom_max.index = p_nom_max.index.astype(str).str.strip()
    else:
        p_nom_max = pd.Series(np.nan, index=p_nom.index)

    extendable = ds["generators_p_nom_extendable"].to_pandas().astype(bool)
    extendable.index = extendable.index.astype(str).str.strip()
    capacities = p_nom_max.where(extendable, p_nom).fillna(p_nom)

    gen_to_bus = ds["generators_bus"].to_pandas().astype(str).str.strip()
    gen_to_bus.index = gen_to_bus.index.astype(str).str.strip()

    bus_to_country = ds["buses_country"].to_pandas().astype(str).str.strip()
    bus_to_country.index = bus_to_country.index.astype(str).str.strip()

    gen_to_country = gen_to_bus.map(bus_to_country)
    df = pd.DataFrame(
        {
            "country": gen_to_country,
            "carrier": carrier,
            "capacity": capacities,
        }
    ).dropna()

    grouped = df.groupby(["country", "carrier"])["capacity"].sum()
    results: Dict[str, pd.Series] = {}
    for country in grouped.index.get_level_values(0).unique():
        if not is_valid_country(country):
            continue
        series = grouped.xs(country)
        series.name = "capacity"
        results[country] = series
    return results


def aggregate_capacity_factors(
    ret_profiles: pd.DataFrame,
    capacities: pd.Series,
) -> pd.DataFrame:
    if capacities is None:
        capacities = pd.Series(dtype=float)

    result: Dict[str, pd.Series] = {}
    for group, carriers in TECH_GROUPS.items():
        available = [c for c in carriers if c in ret_profiles.columns]
        if not available:
            result[group] = pd.Series(0.0, index=ret_profiles.index)
            continue
        weights = capacities.reindex(available).fillna(0.0)
        total_weight = weights.sum()
        data = ret_profiles.loc[:, available]
        if total_weight <= 0 or not np.isfinite(total_weight):
            result[group] = data.mean(axis=1).fillna(0.0)
            continue
        weighted = data.mul(weights, axis=1).sum(axis=1) / total_weight
        result[group] = weighted.fillna(0.0)

    aggregated = pd.DataFrame(result).reindex(ret_profiles.index)
    aggregated.index.name = "snapshot"
    return aggregated


def assemble_country_timeseries(
    country: str,
    demand_df: pd.DataFrame,
    totals: pd.Series,
    cf_df: pd.DataFrame,
) -> pd.DataFrame:
    electricity_profile = weighted_mean_profile(
        demand_df,
        totals,
        ELECTRICITY_COMPONENTS,
    )
    heat_profile = weighted_mean_profile(
        demand_df,
        totals,
        HEAT_COMPONENTS,
    )

    passenger_profile = demand_df.get("passenger_mobility")
    if passenger_profile is None:
        passenger_profile = pd.Series(0.0, index=demand_df.index)
    freight_profile = demand_df.get("freight_mobility")
    if freight_profile is None:
        freight_profile = pd.Series(0.0, index=demand_df.index)

    cf_grouped = cf_df.reindex(demand_df.index, fill_value=0.0)
    cf_grouped = cf_grouped.reindex(columns=list(TECH_GROUPS.keys()), fill_value=0.0)

    timeseries = pd.DataFrame(
        {
            "{PERIODS}": np.arange(1, len(demand_df.index) + 1, dtype=int),
            "Electricity (%_elec)": electricity_profile.values,
            "Space Heating (%_sh)": heat_profile.values,
            "Passanger mobility (%_pass)": passenger_profile.values,
            "Freight mobility (%_freight)": freight_profile.values,
            "PV": cf_grouped["PV"].values,
            "Wind_onshore": cf_grouped["Wind_onshore"].values,
            "Wind_offshore": cf_grouped["Wind_offshore"].values,
            "Hydro_river": cf_grouped["Hydro_river"].values,
            "Solar": cf_grouped["Solar"].values,
        }
    )
    timeseries.columns = TIMESERIES_COLUMNS
    return timeseries


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def main(snakemake: "Snakemake") -> None:
    network_path = Path(snakemake.input.network)
    ret_dir = Path(snakemake.input.RET_timeseries_dir)
    demand_dir = Path(snakemake.input.Demand_timeseries_dir)

    ret_paths = build_country_map(sorted(ret_dir.glob("RET_timeseries_*.csv")))
    demand_paths = build_country_map(sorted(demand_dir.glob("Demand_timeseries_*.csv")))

    output_dir = Path(snakemake.output.assembled_dir)
    timeseries_dir = ensure_output_dir(output_dir)

    if not demand_paths:
        raise ValueError("No demand time series provided.")
    totals_path = demand_dir / "Demand_timeseries_totals.csv"
    totals_frame = load_indexed_csv(totals_path)

    with xr.open_dataset(network_path) as ds:
        capacities = country_carrier_capacities(ds)

    for country in sorted(demand_paths.keys()):
        if country not in ret_paths:
            continue
        demand_df = load_demand_profile(demand_paths[country])
        totals = (
            totals_frame.loc[country].fillna(0.0)
            if country in totals_frame.index
            else pd.Series(dtype=float)
        )

        ret_df = load_ret_profile(ret_paths[country])
        cf_df = aggregate_capacity_factors(ret_df, capacities.get(country))
        cf_df = cf_df.reindex(demand_df.index, fill_value=0.0)

        assembled = assemble_country_timeseries(country, demand_df, totals, cf_df)
        out_file = timeseries_dir / f"time_series_{country}.csv"
        ensure_parent(out_file)
        assembled.to_csv(out_file, index=False)


if __name__ == "__main__":
    main(snakemake)
