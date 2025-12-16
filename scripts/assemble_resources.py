"""Generate EnergyScope resource CSVs by merging biomass data and grid capacities.

The script blends per-country biomass potentials with calculated electricity
import/export limits (converted from MW to GWh) and fills an EnergyScope
resources template. Each country receives its own `resources_<country>.csv`
ready for the EnergyScope conversion workflow.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import sys

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import load_indexed_csv, load_template

if TYPE_CHECKING:
    from snakemake.script import Snakemake

MW_TO_GWH = 8760/1000

# Build country-specific resource tables
def build_country_resource_tables(
    template: pd.DataFrame,
    biomass_resources: pd.DataFrame,
    electricity_caps: pd.DataFrame,
    electricity_resource: str = "ELECTRICITY",
    export_resource: str = "ELEC_EXPORT",
) -> dict[str, pd.DataFrame]:
    resources_by_country: dict[str, pd.DataFrame] = {}
    biomass_index = pd.Index(biomass_resources.index) if not biomass_resources.empty else pd.Index([])
    electricity_index = pd.Index(electricity_caps.index) if not electricity_caps.empty else pd.Index([])
    countries = biomass_index.union(electricity_index)

    for country in countries:
        country_df = template.copy(deep=True)
        if not biomass_resources.empty and country in biomass_resources.index:
            for resource in biomass_resources.columns:
                value = float(biomass_resources.loc[country, resource])
                mask = country_df["resource"].eq(resource)
                if mask.any():
                    country_df.loc[mask, "avail"] = value
        import_gwh = export_gwh = 0.0
        if not electricity_caps.empty and country in electricity_caps.index:
            import_gwh = float(electricity_caps.loc[country, "import_MW"]) * MW_TO_GWH
            export_gwh = float(electricity_caps.loc[country, "export_MW"]) * MW_TO_GWH
        mask_imp = country_df["resource"].eq(electricity_resource)
        if mask_imp.any():
            country_df.loc[mask_imp, "avail"] = import_gwh
        mask_exp = country_df["resource"].eq(export_resource)
        if mask_exp.any():
            country_df.loc[mask_exp, "avail"] = export_gwh
        resources_by_country[country] = country_df.reset_index(drop=True)
    return resources_by_country

# Format resource table to EnergyScope expectations
def energyscope_resource_format(resources: pd.DataFrame) -> pd.DataFrame:
    formatted = resources.copy(deep=True)
    rename_map = {"resource": "parameter name"}
    if "category" in formatted.columns:
        rename_map["category"] = "Category"
    if "subcategory" in formatted.columns:
        rename_map["subcategory"] = "Subcategory"
    formatted = formatted.rename(columns=rename_map)

    column_order = [
        "Category",
        "Subcategory",
        "parameter name",
        "avail",
        "gwp_op",
        "c_op",
    ]
    existing = [col for col in column_order if col in formatted.columns]
    remaining = [col for col in formatted.columns if col not in existing]
    if existing:
        formatted = formatted.loc[:, existing + remaining]
    return formatted


def main(snakemake: "Snakemake") -> None:
    biomass_resources = load_indexed_csv(
        Path(snakemake.input.biomass_resources),
        non_numeric_columns=None,
        coerce_numeric=True,
        fillna_value=0.0,
    )
    electricity_caps = load_indexed_csv(
        Path(snakemake.input.electricity_caps),
        non_numeric_columns=None,
        coerce_numeric=True,
        fillna_value=0.0,
    )

    template_config = getattr(snakemake.params, "template", None)
    template = load_template(template_config, required_column= "resource")

    out_dir = Path(snakemake.output.resources_dir)
    electricity_resource = getattr(snakemake.params, "electricity_resource", "ELECTRICITY")
    export_resource = getattr(snakemake.params, "export_resource", "ELEC_EXPORT")
    out_dir.mkdir(parents=True, exist_ok=True)

    resources_by_country = build_country_resource_tables(
        template=template,
        biomass_resources=biomass_resources,
        electricity_caps=electricity_caps,
        electricity_resource=electricity_resource,
        export_resource=export_resource,
    )

    for country, table in resources_by_country.items():
        formatted = energyscope_resource_format(table)
        formatted.to_csv(out_dir / f"resources_{country}.csv", index=False)

if __name__ == "__main__":
    main(snakemake)
