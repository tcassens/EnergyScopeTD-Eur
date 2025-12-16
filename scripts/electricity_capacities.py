"""Calculate inter-country electricity exchange capacities from a PyPSA network.

For each country represented in the network, the script sums AC line capacities
and DC link capacities to produce total import/export potentials in MW. The
resulting table is exported via Snakemake for the EnergyScope conversion step.
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

# Load and compute inter-country electricity exchange capacities
def load_electricity_capacities(network_file: Path) -> pd.DataFrame:
    with xr.open_dataset(network_file) as ds:
        bus_country = clean_strings(ds["buses_country"].to_pandas()).replace({"": np.nan})
        countries = pd.Index(sorted(bus_country.dropna().unique()))
        totals = pd.DataFrame(0.0, index=countries, columns=["import_mw", "export_mw"])

        line_from = clean_strings(ds["lines_bus0"].to_pandas())
        line_to = clean_strings(ds["lines_bus1"].to_pandas())
        line_cap = ds["lines_s_nom"].to_pandas().astype(float).fillna(0.0)
        for line_id, capacity in line_cap.items():
            if capacity <= 0:
                continue
            origin_bus = line_from.get(line_id)
            target_bus = line_to.get(line_id)
            origin_country = bus_country.get(origin_bus)
            target_country = bus_country.get(target_bus)
            if pd.isna(origin_country) or pd.isna(target_country) or origin_country == target_country:
                continue
            totals.loc[origin_country, ["import_mw", "export_mw"]] += capacity
            totals.loc[target_country, ["import_mw", "export_mw"]] += capacity

        link_from = clean_strings(ds["links_bus0"].to_pandas())
        link_to = clean_strings(ds["links_bus1"].to_pandas())
        link_car = clean_strings(ds["links_carrier"].to_pandas())
        link_cap = ds["links_p_nom"].to_pandas().astype(float).fillna(0.0)
        dc_mask = link_car.str.lower().fillna("").eq("dc")
        for link_id in dc_mask[dc_mask].index:
            capacity = float(link_cap.get(link_id, 0.0))
            if capacity <= 0:
                continue
            exporter = bus_country.get(link_from.get(link_id))
            importer = bus_country.get(link_to.get(link_id))
            if pd.isna(exporter) or pd.isna(importer) or exporter == importer:
                continue
            totals.loc[exporter, "export_mw"] += capacity
            totals.loc[importer, "import_mw"] += capacity

    totals.columns = ["import_MW", "export_MW"]
    totals.index.name = "country"
    return totals


def main(snakemake: "Snakemake") -> None:
    network_path = Path(snakemake.input.network)
    output_path = Path(snakemake.output.caps)
    capacities = load_electricity_capacities(network_path)
    ensure_parent(output_path)
    capacities.to_csv(output_path)


if __name__ == "__main__":
    main(snakemake) 
