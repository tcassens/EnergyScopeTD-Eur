"""Aggregate PyPSA capacities by country and carrier.

This Snakemake script reads a PyPSA network dataset and derives three tables:
installable generator capacities, currently installed generator capacities, and
installed storage capacities. Each table is aggregated per country and carrier,
expressed in GW, and written to CSV files whose locations are provided through
the Snakemake `input`/`output` mappings.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import sys

import numpy as np
import pandas as pd
import xarray as xr

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent, clean_strings

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Prepare mappings from dataset
def prepare_mappings(ds: xr.Dataset) -> dict[str, pd.Series]:
    bus_country = clean_strings(ds["buses_country"].to_pandas())
    gen_bus = clean_strings(ds["generators_bus"].to_pandas())
    gen_country = gen_bus.map(bus_country)
    gen_carrier = clean_strings(ds["generators_carrier"].to_pandas())
    storage_bus = clean_strings(ds["storage_units_bus"].to_pandas())
    storage_country = storage_bus.map(bus_country)
    storage_carrier = clean_strings(ds["storage_units_carrier"].to_pandas())
    countries = pd.Index(sorted(bus_country.dropna().unique()))
    return {
        "countries": countries,
        "bus_country": bus_country,
        "gen_country": gen_country,
        "gen_carrier": gen_carrier,
        "storage_country": storage_country,
        "storage_carrier": storage_carrier,
    }

# Aggregate installable capacities
def aggregate_installable_capacities(ds: xr.Dataset, maps: dict[str, pd.Series]) -> pd.DataFrame:
    generators = ds["generators_p_nom_max"].to_pandas().astype(float)
    gen_country = maps["gen_country"]
    gen_carrier = maps["gen_carrier"]
    valid = gen_country.notna() & gen_carrier.notna()
    grouped = (
        generators[valid]
        .groupby([gen_country[valid], gen_carrier[valid]])
        .sum()
        .unstack(fill_value=0.0)
        .sort_index()
    )
    return grouped / 1000.0

# Aggregate storage capacities
def aggregate_storage_capacities(ds: xr.Dataset, maps: dict[str, pd.Series]) -> pd.DataFrame:
    storage_nom = ds["storage_units_p_nom"].to_pandas().astype(float)
    storage_country = maps["storage_country"]
    storage_carrier = maps["storage_carrier"]
    valid = storage_country.notna() & storage_carrier.notna()
    if not valid.any():
        return pd.DataFrame(index=maps["countries"], dtype=float)
    grouped = (
        storage_nom[valid]
        .groupby([storage_country[valid], storage_carrier[valid]])
        .sum()
        .unstack(fill_value=0.0)
        .sort_index()
    )
    return grouped / 1000.0

# Aggregate installed capacities
def aggregate_installed_capacities(ds: xr.Dataset, maps: dict[str, pd.Series]) -> pd.DataFrame:
    generators = ds["generators_p_nom"].to_pandas().astype(float)
    gen_country = maps["gen_country"]
    gen_carrier = maps["gen_carrier"]
    valid = gen_country.notna() & gen_carrier.notna()
    grouped = (
        generators[valid]
        .groupby([gen_country[valid], gen_carrier[valid]])
        .sum()
        .unstack(fill_value=0.0)
        .sort_index()
    )
    return grouped / 1000.0


def main(snakemake: "Snakemake") -> None:
    network_path = Path(snakemake.input.network)
    outputs = {
        "installable": Path(snakemake.output.installable),
        "installed": Path(snakemake.output.installed),
        "storage_installed": Path(snakemake.output.storage_installed),
    }

    with xr.open_dataset(network_path) as ds:
        maps = prepare_mappings(ds)
        countries = maps["countries"]
        installable = aggregate_installable_capacities(ds, maps).reindex(countries, fill_value=0.0)
        installed = aggregate_installed_capacities(ds, maps).reindex(countries, fill_value=0.0)
        storage = aggregate_storage_capacities(ds, maps).reindex(countries, fill_value=0.0)
        dfs = {"installable": installable, "installed": installed, "storage_installed": storage}

    for kind, frame in dfs.items():
        target = outputs[kind]
        ensure_parent(target)
        frame.to_csv(target)

if __name__ == "__main__":
    main(snakemake) 
