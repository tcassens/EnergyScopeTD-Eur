"""Translate PyPSA demand totals into EnergyScope demand sector CSVs.

Demand totals, mobility tables and energy-reference splits are combined to
derive sector-level demands for electricity, high/low-temperature heat, and
mobility. The values overwrite an EnergyScope template and yield
`demands_<country>.csv` files for the conversion workflow.
"""


from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, Optional

import sys

import numpy as np
import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from conversion_energyscope.scripts._helpers import (
    ensure_parent,
    is_valid_country,
    load_indexed_csv,
    country_key,
    load_template,
)

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Mapping of demand components
HEAT_COMPONENTS = ["rural heat", "urban central heat", "urban decentral heat"]
HEAT_HIGH_CARRIERS = ["gas for industry", "solid biomass for industry", "H2 for industry"]
SECTOR_COLUMNS = ["HOUSEHOLDS", "SERVICES", "INDUSTRY", "TRANSPORTATION"]

ELECTRICITY_REFERENCE = ["electricity residential", "electricity services"]
SPACE_REFERENCE = ["total residential space", "total services space"]
WATER_REFERENCE = ["total residential water", "total services water"]

# Distribute total according to reference values
def distribute_with_reference(total: float, references: Iterable[float]) -> np.ndarray:
    ref_array = np.array(list(references), dtype=float)
    n = len(ref_array)
    if n == 0:
        return np.zeros(0, dtype=float)

    finite_mask = np.isfinite(ref_array) & (ref_array >= 0)
    if finite_mask.sum() == 0:
        return np.full(n, total / n if n else 0.0, dtype=float)

    result = np.zeros(n, dtype=float)
    if finite_mask.sum() == 1:
        idx = np.argmax(finite_mask)
        result[idx] = total
        return result

    ref_filtered = ref_array[finite_mask]
    ref_sum = ref_filtered.sum()
    if ref_sum <= 0:
        return np.full(n, total / n if n else 0.0, dtype=float)
    result[finite_mask] = total * ref_filtered / ref_sum
    return result

# Load mobility totals from CSV
def load_mobility_totals(path: Path) -> Dict[str, Dict[str, float]]:
    df = pd.read_csv(path)
    if df.empty:
        return {}
    df["country"] = df["country"].astype(str).str.upper()
    df = df.set_index("country")
    return df.to_dict(orient="index")

# Get reference values from energy reference table
def get_reference_values(
    energy_reference: Optional[pd.Series],
) -> Dict[str, np.ndarray]:
    if energy_reference is None:
        return {
            "electricity": np.array([np.nan, np.nan], dtype=float),
            "space": np.array([np.nan, np.nan], dtype=float),
            "water": np.array([np.nan, np.nan], dtype=float),
        }

    electricity_ref = np.array(
        [float(energy_reference.get(key, np.nan)) for key in ELECTRICITY_REFERENCE],
        dtype=float,
    )
    space_ref = np.array(
        [float(energy_reference.get(key, np.nan)) for key in SPACE_REFERENCE],
        dtype=float,
    )
    water_ref = np.array(
        [float(energy_reference.get(key, np.nan)) for key in WATER_REFERENCE],
        dtype=float,
    )
    return {
        "electricity": electricity_ref,
        "space": space_ref,
        "water": water_ref,
    }

# Compute sector totals from overall demand totals
def compute_sector_totals(
    totals: pd.Series,
    energy_reference: Optional[pd.Series],
    mobility_info: Optional[Dict[str, float]],
) -> Dict[str, Dict[str, float]]:
    building_electricity = float(totals.get("electricity", 0.0))
    industry_electricity = float(totals.get("industry electricity", 0.0))
    heat_total = float(totals.reindex(HEAT_COMPONENTS).fillna(0.0).sum())
    low_temp_industry = float(totals.get("low-temperature heat for industry", 0.0))
    high_temp_heat = float(totals.reindex(HEAT_HIGH_CARRIERS).fillna(0.0).sum())

    references = get_reference_values(energy_reference)
    space_totals = distribute_with_reference(
        heat_total,
        [np.nansum(references["space"]), np.nansum(references["water"])],
    )
    space_heat_mwh = float(space_totals[0]) if space_totals.size >= 1 else 0.0
    hot_water_mwh = float(space_totals[1]) if space_totals.size >= 2 else 0.0

    elec_split = distribute_with_reference(
        building_electricity,
        references["electricity"],
    )
    elec_res_mwh = float(elec_split[0]) if elec_split.size >= 1 else building_electricity / 2.0
    elec_serv_mwh = float(elec_split[1]) if elec_split.size >= 2 else building_electricity - elec_res_mwh

    space_sector_split = distribute_with_reference(space_heat_mwh, references["space"])
    space_heat_households = float(space_sector_split[0]) if space_sector_split.size >= 1 else space_heat_mwh / 2.0
    space_heat_services = float(space_sector_split[1]) if space_sector_split.size >= 2 else space_heat_mwh - space_heat_households
    space_heat_industry = low_temp_industry / 2.0

    hot_water_sector_split = distribute_with_reference(hot_water_mwh, references["water"])
    hot_water_households = float(hot_water_sector_split[0]) if hot_water_sector_split.size >= 1 else hot_water_mwh / 2.0
    hot_water_services = float(hot_water_sector_split[1]) if hot_water_sector_split.size >= 2 else hot_water_mwh - hot_water_households
    hot_water_industry = low_temp_industry / 2.0

    passenger_total_km = 0.0
    freight_total_km = 0.0
    if mobility_info is not None:
        passenger_total_km = float(mobility_info.get("mio_pkm_wo_aviation", 0.0))
        freight_total_km = float(mobility_info.get("mio_tkm_wo_aviation", 0.0)) 

    totals_dict: Dict[str, Dict[str, float]] = {
        "ELECTRICITY": {
            "HOUSEHOLDS": elec_res_mwh / 1000.0,
            "SERVICES": elec_serv_mwh / 1000.0,
            "INDUSTRY": industry_electricity / 1000.0,
            "TRANSPORTATION": 0.0,
        },
        "HEAT_HIGH_T": {
            "HOUSEHOLDS": 0.0,
            "SERVICES": 0.0,
            "INDUSTRY": high_temp_heat / 1000.0,
            "TRANSPORTATION": 0.0,
        },
        "HEAT_LOW_T_SH": {
            "HOUSEHOLDS": space_heat_households / 1000.0,
            "SERVICES": space_heat_services / 1000.0,
            "INDUSTRY": space_heat_industry / 1000.0,
            "TRANSPORTATION": 0.0,
        },
        "HEAT_LOW_T_HW": {
            "HOUSEHOLDS": hot_water_households / 1000.0,
            "SERVICES": hot_water_services / 1000.0,
            "INDUSTRY": hot_water_industry / 1000.0,
            "TRANSPORTATION": 0.0,
        },
        "MOBILITY_PASSENGER": {
            "HOUSEHOLDS": 0.0,
            "SERVICES": 0.0,
            "INDUSTRY": 0.0,
            "TRANSPORTATION": passenger_total_km,
        },
        "MOBILITY_FREIGHT": {
            "HOUSEHOLDS": 0.0,
            "SERVICES": 0.0,
            "INDUSTRY": 0.0,
            "TRANSPORTATION": freight_total_km,
        },
    }

    return totals_dict

# Generate overrides for demand template
def template_overrides_from_totals(totals: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    overrides: Dict[str, Dict[str, float]] = {}

    def get_row(name: str) -> Dict[str, float]:
        row = totals.get(name, {})
        return {sector: float(row.get(sector, 0.0)) for sector in SECTOR_COLUMNS}

    overrides["ELECTRICITY"] = get_row("ELECTRICITY")
    overrides["HEAT_HIGH_T"] = get_row("HEAT_HIGH_T")
    overrides["HEAT_LOW_T_SH"] = get_row("HEAT_LOW_T_SH")
    overrides["HEAT_LOW_T_HW"] = get_row("HEAT_LOW_T_HW")
    overrides["MOBILITY_PASSENGER"] = get_row("MOBILITY_PASSENGER")
    overrides["MOBILITY_FREIGHT"] = get_row("MOBILITY_FREIGHT")
    return overrides

# Apply overrides to demand template
def apply_demand_template(
    template: pd.DataFrame, overrides: Dict[str, Dict[str, float]]
) -> pd.DataFrame:
    if "parameter name" not in template.columns:
        raise ValueError("Demand template must contain a 'parameter name' column")

    table = template.copy(deep=True)
    table = table.set_index("parameter name", drop=False)

    for parameter, sector_values in overrides.items():
        if parameter not in table.index:
            new_row = {column: 0.0 if column in SECTOR_COLUMNS else "" for column in table.columns}
            new_row["parameter name"] = parameter
            table.loc[parameter] = new_row
        for sector in SECTOR_COLUMNS:
            if sector not in table.columns:
                continue
            table.at[parameter, sector] = float(sector_values.get(sector, 0.0))

    return table.reset_index(drop=True)


def main(snakemake: "Snakemake") -> None:
    demand_timeseries_dir = Path(snakemake.input.Demand_timeseries_dir)
    demand_totals_path = demand_timeseries_dir / "Demand_timeseries_totals.csv"
    mobility_path = Path(snakemake.input.mobility_totals)
    energy_reference_path = Path(snakemake.input.energy_reference)
    output_dir = Path(snakemake.output.demands_dir)
    template_config = getattr(snakemake.params, "template", None)
    if template_config is None:
        raise ValueError("No demand template configured")

    demand_template = load_template(template_config, required_column="parameter name")

    demand_totals = load_indexed_csv(demand_totals_path)
    mobility_totals = load_mobility_totals(mobility_path)
    energy_reference = load_indexed_csv(energy_reference_path)

    output_dir.mkdir(parents=True, exist_ok=True)

    for country in sorted(demand_totals.index):
        if not is_valid_country(country):
            continue
        totals = demand_totals.loc[country]
        ref = energy_reference.loc[country] if country in energy_reference.index else None
        mobility_info = mobility_totals.get(country)
        sector_totals = compute_sector_totals(totals, ref, mobility_info)
        overrides = template_overrides_from_totals(sector_totals)
        table = apply_demand_template(demand_template, overrides)

        out_file = output_dir / f"demands_{country}.csv"
        ensure_parent(out_file)
        table.to_csv(out_file, index=False)


if __name__ == "__main__":
    main(snakemake)
