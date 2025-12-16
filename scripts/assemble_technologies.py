"""Build country-specific EnergyScope technology tables from PyPSA capacities.

The script reads installable potentials, installed capacities, and storage
limits, injects those values into an EnergyScope technology template (filling
`f_max` per available technology) and exports formatted `technologies_<country>.csv`
files that the conversion workflow consumes.
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

from conversion_energyscope.scripts._helpers import clean_strings, ensure_parent, load_template, load_indexed_csv

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Maps from PyPSA technology columns to technology names
WIND_COLUMN_MAP = {
    "onwind": "wind_onshore_GW",
    "offwind-ac": "wind_offshore_ac_GW",
    "offwind-dc": "wind_offshore_dc_GW",
    "offwind-float": "wind_offshore_float_GW",
}

# Maps from PyPSA technology columns to technology names
PV_COLUMN_MAP = {
    "solar": "pv_ground_GW",
    "solar rooftop": "pv_rooftop_GW",
    "solar-hsat": "pv_hsat_GW",
}

# Build country-specific technology tables
def build_country_technology_tables(
        template: pd.DataFrame,
        installable: pd.DataFrame,
        installed: pd.DataFrame,
        storage: pd.DataFrame,
) -> dict[str, pd.DataFrame]:

    countries = installable.index
    country_tables: dict[str, pd.DataFrame] = {}

    wind_offshore_installable = [c for c in WIND_COLUMN_MAP if c != "onwind" and c in installable.columns]
    pv_installable_cols = [c for c in PV_COLUMN_MAP if c in installable.columns]

    # Helper to get value or zero
    def value_from(frame: pd.DataFrame, country: str, column: str) -> float:
        if column not in frame.columns:
            return 0.0
        value = frame.at[country, column]
        return float(value) if pd.notna(value) else 0.0

    for country in countries:
        country_df = template.copy(deep=True)

        mask_w_on = country_df["technology"] == "WIND_ONSHORE"
        if mask_w_on.any():
            country_df.loc[mask_w_on, "f_max"] = round(value_from(installable, country, "onwind"), 3)

        mask_w_off = country_df["technology"] == "WIND_OFFSHORE"
        if mask_w_off.any():
            installable_sum = float(installable.loc[country, wind_offshore_installable].sum()) if wind_offshore_installable else 0.0
            country_df.loc[mask_w_off, "f_max"] = round(installable_sum, 3)

        mask_pv = country_df["technology"] == "PV"
        if mask_pv.any():
            installable_sum = float(installable.loc[country, pv_installable_cols].sum()) if pv_installable_cols else 0.0
            country_df.loc[mask_pv, "f_max"] = round(installable_sum, 3)

        mask_phs = country_df["technology"] == "PHS"
        if mask_phs.any():
            phs_capacity = value_from(storage, country, "PHS")
            country_df.loc[mask_phs, "f_max"] = round(phs_capacity, 3)

        mask_river = country_df["technology"] == "HYDRO_RIVER"
        if mask_river.any():
            river_capacity = value_from(installed, country, "ror")
            country_df.loc[mask_river, "f_max"] = round(river_capacity, 3)

        mask_nuc = country_df["technology"] == "NUCLEAR"
        if mask_nuc.any():
            nuclear_capacity = value_from(installed, country, "nuclear")
            country_df.loc[mask_nuc, "f_max"] = round(nuclear_capacity, 3)

        country_tables[country] = country_df.reset_index(drop=True)

    return country_tables

# Format technology table to EnergyScope expectations
def energyscope_technology_format(technologies: pd.DataFrame) -> pd.DataFrame:
    formatted = technologies.copy(deep=True)
    rename_map = {
        "category": "Category",
        "subcategory": "Subcategory",
        "tech_name": "Technologies name",
        "technology": "Technologies param",
        "comment": "Comment",
    }
    formatted = formatted.rename(columns=rename_map)

    for col in ["Category", "Subcategory", "Technologies name", "Comment"]:
        if col not in formatted.columns:
            formatted[col] = ""

    column_order = [
        "Category",
        "Subcategory",
        "Technologies name",
        "Technologies param",
        "c_inv",
        "c_maint",
        "gwp_constr",
        "lifetime",
        "c_p",
        "fmin_perc",
        "fmax_perc",
        "f_min",
        "f_max",
        "Comment",
    ]
    existing = [c for c in column_order if c in formatted.columns]
    remaining = [c for c in formatted.columns if c not in existing]
    return formatted.loc[:, existing + remaining]


def main(snakemake: "Snakemake") -> None:
    installable = load_indexed_csv(
        Path(snakemake.input.installable),
        non_numeric_columns=None,
        coerce_numeric=True,
        fillna_value=0.0,
    )
    installed = load_indexed_csv(
        Path(snakemake.input.installed),
        non_numeric_columns=None,
        coerce_numeric=True,
        fillna_value=0.0,
    )
    storage = load_indexed_csv(
        Path(snakemake.input.storage_installed),
        non_numeric_columns=None,
        coerce_numeric=True,
        fillna_value=0.0,
    )
    
    template_config = getattr(snakemake.params, "template", None)
    template = load_template(template_config, required_column="technology")

    country_tables = build_country_technology_tables(
        template=template,
        installable=installable,
        installed=installed,
        storage=storage,
    )

    out_dir = Path(snakemake.output.technologies_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for country, table in country_tables.items():
        formatted = energyscope_technology_format(table)
        formatted.to_csv(out_dir / f"technologies_{country}.csv", index=False)

if __name__ == "__main__":
    main(snakemake) 
