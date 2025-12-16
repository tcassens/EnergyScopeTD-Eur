"""Map biomass potential columns to EnergyScope resource names.

This script takes the country-level biomass potentials CSV and a
mapping of source columns to EnergyScope resource names (from config.yaml)
and produces the aggregated resource CSV used downstream in the conversion.
"""


from pathlib import Path
from typing import TYPE_CHECKING, Mapping

import sys

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import ensure_parent, load_indexed_csv

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Clean up mapping to ensure all keys and values are strings
def _normalise_mapping(mapping: Mapping[str, str] | None) -> dict[str, str] | None:
    if not mapping:
        return None
    return {str(source): str(resource) for source, resource in mapping.items()}

# Aggregate biomass potentials according to the provided mapping
def assign_biomass(
    potentials: pd.DataFrame,
    mapping: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    mapping_dict = dict(mapping or {})
    data: dict[str, pd.Series] = {}
    for source_col, resource in mapping_dict.items():
        if source_col not in potentials.columns:
            continue
        series = potentials[source_col]
        if resource in data:
            data[resource] = data[resource].add(series, fill_value=0.0)
        else:
            data[resource] = series
    if not data:
        return pd.DataFrame(index=potentials.index)
    df = pd.DataFrame(data).reindex(potentials.index).fillna(0.0)
    df.index.name = potentials.index.name
    return df


def main(snakemake: "Snakemake") -> None:
    potentials_path = Path(snakemake.input.potentials)
    resources_path = Path(snakemake.output.resources)
    mapping_param = getattr(snakemake.params, "mapping", None)

    potentials = load_indexed_csv(potentials_path)
    mapping = _normalise_mapping(mapping_param)
    resources = assign_biomass(potentials, mapping=mapping)

    ensure_parent(resources_path)
    resources.to_csv(resources_path)


if __name__ == "__main__":
    main(snakemake) 
