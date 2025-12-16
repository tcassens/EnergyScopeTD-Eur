"""Produce RET (Renewable Energy Technology) time series per country.

Given a PyPSA-Eur network file, the script computes technology-specific capacity
factor profiles aggregated by country, weighting each generator by its installed
or extendable capacity. Snakemake uses the resulting CSVs as RET inputs for the
EnergyScope conversion workflow.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Dict

import sys

import numpy as np
import pandas as pd
import xarray as xr

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent, is_valid_country

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Compute a capacity-factor profile weighted by installed/potential capacity
def weighted_capacity_factor(
    cf_panel: pd.DataFrame,
    capacities: pd.Series,
    mask: pd.Series,
) -> pd.Series:
    ids = mask[mask].index.intersection(cf_panel.columns)
    if ids.empty:
        return pd.Series(0.0, index=cf_panel.index)
    data = cf_panel.loc[:, ids]
    weights = capacities.loc[ids]
    if np.isinf(weights).any(): # Fallback to simple average if any capacity is infinite
        return data.mean(axis=1)
    weights = weights.fillna(0.0)
    total_capacity = weights.sum()
    if total_capacity <= 0 or not np.isfinite(total_capacity): # Fallback to simple average if no valid capacity
        return data.mean(axis=1)
    return data.multiply(weights, axis=1).sum(axis=1) / total_capacity # Weighted average

# Compute capacity factors per country without aggregating carriers
def country_capacity_factors(ds: xr.Dataset) -> Dict[str, pd.DataFrame]:
    cf_panel = ds["generators_t_p_max_pu"].to_pandas()
    cf_panel.columns = cf_panel.columns.astype(str).str.strip()

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

    carriers_available = carrier.dropna().unique()

    results: Dict[str, pd.DataFrame] = {}
    for country in gen_to_country.dropna().unique():
        if not is_valid_country(country):
            continue
        mask_base = gen_to_country.eq(country)
        factors: Dict[str, pd.Series] = {}
        for tech in sorted(carriers_available):
            tech_mask = mask_base & carrier.eq(tech)
            if not tech_mask.any():
                continue
            factors[tech] = weighted_capacity_factor(
                cf_panel,
                capacities,
                tech_mask,
            )
        if not factors:
            continue
        country_profile = pd.DataFrame(factors)
        country_profile.index.name = "snapshot"
        results[country] = country_profile
    return results

# Write country-specific CSV files
def write_country_csvs(country_profiles: Dict[str, pd.DataFrame], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for country, profile in sorted(country_profiles.items()):
        out_file = output_dir / f"RET_timeseries_{country}.csv"
        ensure_parent(out_file)
        profile.to_csv(out_file)


def main(snakemake: "Snakemake") -> None:
    output_dir = Path(snakemake.output.RET_timeseries_dir)
    with xr.open_dataset(snakemake.input.network) as ds:
        cf_per_country_dict = country_capacity_factors(ds)
    write_country_csvs(cf_per_country_dict, output_dir)

if __name__ == "__main__":
    main(snakemake) 
