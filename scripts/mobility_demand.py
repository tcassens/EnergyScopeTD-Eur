"""Assemble mobility demand (pkm/tkm/aviation energy) for all countries.

The script reads EU-27 transport activity and energy from JRC-IDEES and
computes non-EU transport activity from aggregated PyPSA-Eur `energy_totals.csv`
using EU-27 energy-per-activity intensities.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Dict

import sys

import pandas as pd

WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "EnergyScopeTD-Eur"
if str(WORKFLOW_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_ROOT))

from scripts._helpers import ensure_parent

if TYPE_CHECKING:
    from snakemake.script import Snakemake


# EU-27 codes used in the transport demand dataset
EU27_COUNTRIES = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR", "HR",
    "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO", "SE", "SI", "SK",
]
NON_EU_COUNTRIES = ["AL", "BA", "CH", "GB", "ME", "MK", "NO", "RS", "XK"]
TWH_TO_KTOE = 1000.0 / 11.63

# Extract transport data for EU-27 countries directly from JRC-IDEES dataset
def get_transport_EU_27(year: int, base_eu: str | Path) -> Dict[str, pd.DataFrame]:
    path = Path(base_eu)

    if not path.exists():
        raise FileNotFoundError(f"Transport path does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"Expected transport directory, got file: {path}")

    data: dict[str, dict[str, float]] = {}

    for country in EU27_COUNTRIES:
        country_dir = path / country
        candidates = sorted(country_dir.glob(f"JRC-IDEES-*_Transport_{country}.xlsx"))
        if not candidates:
            raise FileNotFoundError(
                f"Missing transport file for {country} under {country_dir}"
            )
        country_file = candidates[-1]

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

# Extract transport data for non-EU countries from aggregated PyPSA-Eur totals
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

    def code_value(*candidates: str) -> float:
        for key in candidates:
            if key in eu_codes.index:
                return float(eu_codes[key])
        raise KeyError(f"None of the code variants found: {candidates}")

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
        code_value(
            "FEC.ktoe.EU27.Tr.Avia.Passenger.IntraEEAwCHUK",
            "FEC.ktoe.EU27.Tr.Avia.Passenger.IntraEEAwUK",
            "FEC.ktoe.EU27.Tr.Avia.Passenger.IntraEU27",
        )
        + code_value(
            "FEC.ktoe.EU27.Tr.Avia.Passenger.ExtraEEAwCHUK",
            "FEC.ktoe.EU27.Tr.Avia.Passenger.ExtraEEAwUK",
            "FEC.ktoe.EU27.Tr.Avia.Passenger.ExtraEU27",
        )
    )

    avia_freight_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight"]
    avia_freight_activity = eu_codes["Activity.Mtkm.EU27.Tr.Avia.Freight"]
    avia_freight_dom_energy = eu_codes["FEC.ktoe.EU27.Tr.Avia.Freight.Domestic"]
    avia_freight_int_energy = (
        code_value(
            "FEC.ktoe.EU27.Tr.Avia.Freight.IntraEEAwCHUK",
            "FEC.ktoe.EU27.Tr.Avia.Freight.IntraEEAwUK",
            "FEC.ktoe.EU27.Tr.Avia.Freight.IntraEU27",
        )
        + code_value(
            "FEC.ktoe.EU27.Tr.Avia.Freight.ExtraEEAwCHUK",
            "FEC.ktoe.EU27.Tr.Avia.Freight.ExtraEEAwUK",
            "FEC.ktoe.EU27.Tr.Avia.Freight.ExtraEU27",
        )
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

    # Consume aggregated transport-relevant energy totals directly.
    # Required columns: country, year, road/rail/domestic+international aviation/navigation.
    base_non_eu_path = Path(base_non_eu)
    if not (base_non_eu_path.is_file() and base_non_eu_path.suffix == ".csv"):
        raise ValueError(
            f"Expected non-EU transport input as CSV file, got: {base_non_eu_path}"
        )

    totals = pd.read_csv(base_non_eu_path)
    required_cols = {
        "country",
        "year",
        "total road",
        "total rail",
        "total domestic aviation",
        "total international aviation",
        "total domestic navigation",
        "total international navigation",
    }
    missing_cols = sorted(required_cols.difference(totals.columns))
    if missing_cols:
        raise ValueError(
            f"Missing required columns in {base_non_eu_path}: {missing_cols}"
        )

    year_totals = totals[totals["year"] == year].set_index("country")
    for country in NON_EU_COUNTRIES:
        if country not in year_totals.index:
            raise KeyError(
                f"Country '{country}' missing in {base_non_eu_path} for year {year}"
            )
        row = year_totals.loc[country]
        road_total_ktoe = float(row["total road"]) * TWH_TO_KTOE
        rail_total_ktoe = float(row["total rail"]) * TWH_TO_KTOE
        dom_avi_ktoe = float(row["total domestic aviation"]) * TWH_TO_KTOE
        int_avi_ktoe = float(row["total international aviation"]) * TWH_TO_KTOE
        navigation_ktoe = (
            float(row["total domestic navigation"])
            + float(row["total international navigation"])
        ) * TWH_TO_KTOE

        country_transport[country] = build_transport_dataframe(
            road_total_ktoe,
            rail_total_ktoe,
            dom_avi_ktoe,
            int_avi_ktoe,
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
