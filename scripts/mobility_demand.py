"""Assemble mobility demand (pkm/tkm/aviation energy) for all countries.

The script mixes direct readings from national JRC-IDEES files (EU-27) with
derived estimates for non-EU countries based on Eurostat energy balance
statistics. Results are aligned into the columns expected by EnergyScope
(passenger-km, tonne-km) and additional information for aviation-only 
components and ktoe. Output is written as a single CSV that Snakemake passes 
downstream in the conversion workflow.
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

from conversion_energyscope.scripts._helpers import ensure_parent, clean_strings

if TYPE_CHECKING:
    from snakemake.script import Snakemake


# EU-27 codes used in the transport demand dataset
EU27_COUNTRIES = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR", "HR",
    "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO", "SE", "SI", "SK",
]

# Extract transport data for EU-27 countries directly from JRC-IDEES dataset
def get_transport_EU_27(year: int, base_eu: str | Path) -> Dict[str, pd.DataFrame]:
    path = Path(base_eu)

    if not path.exists():
        raise FileNotFoundError(f"Transport path does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"Expected transport directory, got file: {path}")

    data: dict[str, dict[str, float]] = {}

    for country in EU27_COUNTRIES:
        country_file = path / country / f"JRC-IDEES-2021_Transport_{country}.xlsx"
        if not country_file.exists():
            raise FileNotFoundError(f"Missing transport file for {country}: {country_file}")

        raw = pd.read_excel(country_file, sheet_name="Transport", header=None)

        try:
            year_col = raw.columns[(raw.iloc[0] == year)][0]
        except IndexError as exc:
            raise KeyError(f"{year=} not found in {country_file}") from exc

        try:
            code_col = raw.columns[(raw.iloc[0] == "Code")][0]
        except IndexError as exc:
            raise KeyError(f"'Code' column not found in {country_file}") from exc

        codes = (
            raw[raw[code_col].notna()]
            .set_index(code_col)[year_col]
            .astype(float)
        )

        try:
            total_pkm = codes[f"Activity.Mpkm.{country}.Tr"]
            avia_pkm = codes[f"Activity.Mpkm.{country}.Tr.Avia.Passenger"]
            total_tkm = codes[f"Activity.Mtkm.{country}.Tr"]
            avia_tkm = codes[f"Activity.Mtkm.{country}.Tr.Avia.Freight"]
            ktoe_aviation = codes[f"FEC.ktoe.{country}.Tr.Avia.Freight"] + codes[f"FEC.ktoe.{country}.Tr.Avia.Passenger"]

        except KeyError as exc:
            raise KeyError(f"Transport data missing for country '{country}' in {country_file}") from exc

        data[country] = {
            "mio_pkm": total_pkm,
            "mio_pkm_aviation": avia_pkm,
            "mio_pkm_wo_aviation": total_pkm - avia_pkm,
            "mio_tkm": total_tkm,
            "mio_tkm_aviation": avia_tkm,
            "mio_tkm_wo_aviation": total_tkm - avia_tkm,
            "ktoe_aviation": ktoe_aviation,
        }
    return data

# Extract transport data for non-EU countries using EU averages and national energy balances (Eurostat, Swiss Federal Office of Energy)
def get_transport_non_EU(year: int, base_eu27: str | Path, base_non_eu: str | Path) -> Dict[str, pd.DataFrame]:
    country_transport: Dict[str, pd.DataFrame] = {}
    year_label = str(year)

    eu_df = pd.read_excel(base_eu27, sheet_name="Transport", header=None)

    year_col = eu_df.columns[(eu_df.iloc[0] == year)].tolist()[0]
    code_col = eu_df.columns[(eu_df.iloc[0] == "Code")].tolist()[0]
    eu_codes = (
        eu_df[eu_df[code_col].notna()]
        .set_index(code_col)[year_col]
        .astype(float)
    )

    road_pass_energy = eu_codes["FEC.ktoe.EU27.Tr.Road.Passenger"]
    road_pass_activity = eu_codes["Activity.Mpkm.EU27.Tr.Road.Passenger"]
    road_freight_energy = eu_codes["FEC.ktoe.EU27.Tr.Road.Freight"]
    road_freight_activity = eu_codes["Activity.Mtkm.EU27.Tr.Road.Freight"]

    rail_pass_energy = eu_codes["FEC.ktoe.EU27.Tr.Rail.Passenger"]
    rail_pass_activity = eu_codes["Activity.Mpkm.EU27.Tr.Rail.Passenger"]
    rail_freight_energy = eu_codes["FEC.ktoe.EU27.Tr.Rail.Freight"]
    rail_freight_activity = eu_codes["Activity.Mtkm.EU27.Tr.Rail.Freight"]

    avia_pass_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Passenger"]
    avia_pass_activity = eu_codes["Activity.Mpkm.EU27.Tr.Avia.Passenger"]
    avia_pass_dom_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Passenger.Domestic"]
    avia_pass_int_energy = (
        eu_codes["FEC.ktoe.EU27.Tr.Avia.Passenger.IntraEEAwUK"]
        + eu_codes["FEC.ktoe.EU27.Tr.Avia.Passenger.ExtraEEAwUK"]
    )

    avia_freight_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight"]
    avia_freight_activity = eu_codes["Activity.Mtkm.EU27.Tr.Avia.Freight"]
    avia_freight_dom_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight.Domestic"]
    avia_freight_int_energy = (
        eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight.IntraEEAwUK"]
        + eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight.ExtraEEAwUK"]
    )

    passenger_nonavi_intensity = (
        road_pass_energy + rail_pass_energy
    ) / (road_pass_activity + rail_pass_activity)
    freight_nonavi_intensity = (
        road_freight_energy + rail_freight_energy
    ) / (road_freight_activity + rail_freight_activity)
    passenger_avi_intensity = avia_pass_energy / avia_pass_activity
    freight_avi_intensity = avia_freight_energy / avia_freight_activity

    road_passenger_share = road_pass_energy / (road_pass_energy + road_freight_energy)
    rail_passenger_share = rail_pass_energy / (rail_pass_energy + rail_freight_energy)
    dom_avi_passenger_share = avia_pass_dom_energy / (
        avia_pass_dom_energy + avia_freight_dom_energy
    )
    int_avi_passenger_share = avia_pass_int_energy / (
        avia_pass_int_energy + avia_freight_int_energy
    )

    def build_transport_dataframe(
        road_total_ktoe: float,
        rail_total_ktoe: float,
        dom_avi_ktoe: float,
        int_avi_ktoe: float,
        navigation_ktoe: float,
    ) -> pd.DataFrame:
        road_passenger_energy = road_total_ktoe * road_passenger_share
        road_freight_energy = road_total_ktoe - road_passenger_energy

        rail_passenger_energy = rail_total_ktoe * rail_passenger_share
        rail_freight_energy = rail_total_ktoe - rail_passenger_energy

        dom_avi_passenger_energy = dom_avi_ktoe * dom_avi_passenger_share
        dom_avi_freight_energy = dom_avi_ktoe - dom_avi_passenger_energy

        int_avi_passenger_energy = int_avi_ktoe * int_avi_passenger_share
        int_avi_freight_energy = int_avi_ktoe - int_avi_passenger_energy

        passenger_nonavi_energy = road_passenger_energy + rail_passenger_energy
        freight_nonavi_energy = (
            road_freight_energy + rail_freight_energy + navigation_ktoe
        )

        passenger_nonavi_mpkm = (
            passenger_nonavi_energy / passenger_nonavi_intensity
            if passenger_nonavi_intensity > 0
            else 0.0
        )
        passenger_avi_mpkm = (
            (dom_avi_passenger_energy + int_avi_passenger_energy) / passenger_avi_intensity
            if passenger_avi_intensity > 0
            else 0.0
        )

        freight_nonavi_mtkm = (
            freight_nonavi_energy / freight_nonavi_intensity
            if freight_nonavi_intensity > 0
            else 0.0
        )
        freight_avi_mtkm = (
            (dom_avi_freight_energy + int_avi_freight_energy) / freight_avi_intensity
            if freight_avi_intensity > 0
            else 0.0
        )

        total_mio_pkm = passenger_nonavi_mpkm + passenger_avi_mpkm
        total_mio_tkm = freight_nonavi_mtkm + freight_avi_mtkm
        avia_energy_ktoe = dom_avi_ktoe + int_avi_ktoe

        return pd.DataFrame(
            {
                "mio_pkm": [total_mio_pkm],
                "mio_pkm_aviation": [passenger_avi_mpkm],
                "mio_pkm_wo_aviation": [total_mio_pkm - passenger_avi_mpkm],
                "mio_tkm": [total_mio_tkm],
                "mio_tkm_aviation": [freight_avi_mtkm],
                "mio_tkm_wo_aviation": [total_mio_tkm - freight_avi_mtkm],
                "ktoe_aviation": [avia_energy_ktoe],
            }
        )

    # Helper function to read transport flows from national energy balance sheets
    def read_transport_flows(country_code: str, base_non_eu: str | Path) -> Dict[str, float]:
        file_country = "UK" if country_code == "GB" else country_code
        data_file = base_non_eu / (
            f"{file_country}-Energy-balance-sheets-April-2023-edition.xlsb"
        )
        raw = pd.read_excel(data_file, sheet_name=year_label, engine="pyxlsb")

        flow_col = raw.columns[1]
        sub_flow_col = raw.columns[2]
        code_col = raw.columns[7]
        value_col = raw.columns[8]

        table = (
            raw[[flow_col, sub_flow_col, code_col, value_col]]
            .rename(
                columns={
                    flow_col: "flow",
                    sub_flow_col: "sub_flow",
                    code_col: "code",
                    value_col: "value",
                }
            )
        )

        series = pd.to_numeric(
            table[table["code"].notna()].set_index("code")["value"],
            errors="coerce",
        ).fillna(0.0)

        codes = {
            "road_total_ktoe": "FC_TRA_ROAD_E",
            "rail_total_ktoe": "FC_TRA_RAIL_E",
            "dom_avi_ktoe": "FC_TRA_DAVI_E",
            "int_avi_ktoe": "INTAVI",
            "dom_nav_ktoe": "FC_TRA_DNAVI_E",
            "int_nav_ktoe": "FC_TRA_INAVI_E",
        }

        return {
            name: float(series.get(code, 0.0))
            for name, code in codes.items()
        }

    for country in ["AL", "BA", "CH", "GB", "ME", "MK", "NO", "RS", "XK"]:
        if country == "CH":
            PJ_TO_KTOE = 23.88458966275
            swiss_file = Path("data/switzerland-new_format-all_years.csv") 
            swiss_df = pd.read_csv(swiss_file, index_col="item")

            swiss_row = swiss_df[year_label]
            road_total_ktoe = swiss_row["total road"] * PJ_TO_KTOE
            rail_total_ktoe = swiss_row["total rail"] * PJ_TO_KTOE
            dom_avi_ktoe = swiss_row.get("total domestic aviation", 0.0) * PJ_TO_KTOE
            int_avi_ktoe = swiss_row.get("total international aviation", 0.0) * PJ_TO_KTOE
            navigation_ktoe = (
                swiss_row.get("total domestic navigation", 0.0)
                + swiss_row.get("total international navigation", 0.0)
            ) * PJ_TO_KTOE

            country_transport[country] = build_transport_dataframe(
                road_total_ktoe,
                rail_total_ktoe,
                dom_avi_ktoe,
                int_avi_ktoe,
                navigation_ktoe,
            )
        else:
            flows = read_transport_flows(country, base_non_eu=base_non_eu)
            navigation_ktoe = flows["dom_nav_ktoe"] + flows["int_nav_ktoe"]

            country_transport[country] = build_transport_dataframe(
                flows["road_total_ktoe"],
                flows["rail_total_ktoe"],
                flows["dom_avi_ktoe"],
                flows["int_avi_ktoe"],
                navigation_ktoe,
            )

    return country_transport

def get_transport_all_countries(year: int, base_eu: str | Path, base_non_eu: str | Path, base_eu27: str | Path) -> Dict[str, pd.DataFrame]:
    eu_data = get_transport_EU_27(year, base_eu)
    non_eu_data = get_transport_non_EU(year, base_eu27, base_non_eu)

    desired_columns = [
        "mio_pkm",
        "mio_pkm_aviation",
        "mio_pkm_wo_aviation",
        "mio_tkm",
        "mio_tkm_aviation",
        "mio_tkm_wo_aviation",
        "ktoe_aviation",
    ]

    # Helper to align dataframes to desired columns
    def ensure_columns(dataset: Dict[str, pd.DataFrame | Dict[str, float]]) -> Dict[str, pd.Series]:
        aligned: Dict[str, pd.Series] = {}
        for country, df in dataset.items():
            if isinstance(df, pd.DataFrame):
                row = df.iloc[0]
            elif isinstance(df, dict):
                row = pd.Series(df)
            else:
                raise TypeError(f"Unexpected transport data type for '{country}': {type(df)}")
            aligned[country] = row.reindex(desired_columns)
        return aligned

    merged = {**ensure_columns(eu_data), **ensure_columns(non_eu_data)}
    df_mobility_per_country = pd.DataFrame.from_dict(merged, orient="index")
    df_mobility_per_country.index.name = "country"
    return df_mobility_per_country

def main(snakemake: "Snakemake") -> None:
    input_path_eu_countries = Path(snakemake.input.mobility_eu_countries)
    input_path_eu27 = Path(snakemake.input.mobility_eu27)
    input_path_non_eu = Path(snakemake.input.mobility_non_eu_energy_demands)
    output_path = Path(snakemake.output.mobility_demands)

    mobility = get_transport_all_countries(
        year=2019,
        base_eu=input_path_eu_countries,
        base_non_eu=input_path_non_eu,
        base_eu27=input_path_eu27,
    )
    ensure_parent(output_path)
    mobility.to_csv(output_path)

if __name__ == "__main__":
    main(snakemake) 
