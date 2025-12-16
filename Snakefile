from pathlib import Path

from snakemake.io import directory

# Load configuration
configfile: "EnergyScopeTD-Eur/config.yaml"

# Parameters
WEATHER_YEARS = config.get("weather_years", [2013])
CLUSTERS = config.get("clusters", [39])
HORIZONS = config.get("horizons", [2050])
COUNTRIES = config.get("countries", [])
TD_CONFIG = config.get("typical_days", {})

# Templates
RESOURCE_TEMPLATE = config.get("resource_template")
BIOMASS_MAPPING = config.get("biomass_mapping")
TECH_TEMPLATE = config.get("technology_template")
DEMAND_TEMPLATE = config.get("demand_template")

# PyPSA outputs
PYPSA_NET_PATH = config["paths"]["pypsa_network"]
BIOMASS_CSV = config["paths"]["biomass_csv"]
MOBILITY_EU27_CSV = config["paths"]["eu27_mobility_csv"]
MOBILITY_EU_COUNTRIES_CSV = config["paths"]["europe_mobility_csv"]
MOBILITY_NON_EU_ENERGY_DEMANDS_CSV = config["paths"]["europe_energy_damands_csv"]
TRANSPORT_DEMAND_CSV = config["paths"]["transport_demand_csv"]
ENERGY_TOTALS_CSV = config["paths"]["energy_totals_csv"]
TEMP_SOIL_NC = config["paths"]["temp_soil_nc"]
TEMP_AIR_NC = config["paths"]["temp_air_nc"]
POP_CSV = config["paths"]["pop_csv"]
COP_NC = config["paths"]["cop_profiles_nc"]
HYDRO_PROFILE_NC = config["paths"]["profile_hydro_nc"]
INDUSTRIAL_PRODUCTION_CSV = config["paths"]["industrial_production_csv"]

# Output directories
BASE_OUT_DIR = "EnergyScopeTD-Eur/output/Europe_{year}_{clusters}_{horizon}"
INTERMEDIATE_DIR = BASE_OUT_DIR + "/intermediate"
RESOURCES_DIR = BASE_OUT_DIR + "/resources"
TECH_DIR = BASE_OUT_DIR + "/technologies"
TIME_DIR = BASE_OUT_DIR + "/timeseries"
DEMAND_DIR = BASE_OUT_DIR + "/demands"
DAT_DIR = BASE_OUT_DIR + "/dat_files"

# Rule to build all outputs
rule all:
    input:
        expand(
            DAT_DIR + "/ESTD_{country}_data_core.dat",
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
            country=COUNTRIES,
        ),
        expand(
            DAT_DIR + "/ESTD_{country}_12TD.dat",
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
            country=COUNTRIES,
        ),
        expand(
            INTERMEDIATE_DIR + "/hydro_inflow_timeseries.csv",
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
        ),
        expand(
            INTERMEDIATE_DIR + "/soil_temp_time_series.csv",    
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
        ),
        expand(
            INTERMEDIATE_DIR + "/air_temp_time_series.csv",    
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
        ),
        expand(
            INTERMEDIATE_DIR + "/cop_time_series",  
            year=WEATHER_YEARS,
            clusters=CLUSTERS,
            horizon=HORIZONS,
        ),

rule biomass_potentials:
    input:
        biomass=lambda wc: BIOMASS_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        potentials=INTERMEDIATE_DIR + "/biomass_potentials_per_country.csv"
    conda:
        "envs/resources.yaml"
    script:
        "scripts/biomass_potentials.py"


rule biomass_mapping:
    input:
        potentials=INTERMEDIATE_DIR + "/biomass_potentials_per_country.csv"
    output:
        resources=INTERMEDIATE_DIR + "/biomass_potentials_per_country_energyscope.csv"
    params:
        mapping=BIOMASS_MAPPING
    conda:
        "envs/resources.yaml"
    script:
        "scripts/biomass_mapping.py"


rule electricity_capacities:
    input:
        network=lambda wc: PYPSA_NET_PATH.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        caps=INTERMEDIATE_DIR + "/electricity_capacities_per_country.csv"
    conda:
        "envs/resources.yaml"
    script:
        "scripts/electricity_capacities.py"


rule assemble_resources:
    input:
        biomass_resources=INTERMEDIATE_DIR + "/biomass_potentials_per_country_energyscope.csv",
        electricity_caps=INTERMEDIATE_DIR + "/electricity_capacities_per_country.csv"
    output:
        resources_dir=directory(RESOURCES_DIR),
        resources=expand(RESOURCES_DIR + "/resources_{country}.csv", country=COUNTRIES, allow_missing=True)
    params:
        template=RESOURCE_TEMPLATE
    conda:
        "envs/resources.yaml"
    script:
        "scripts/assemble_resources.py"


rule aggregate_capacities:
    input:
        network=lambda wc: PYPSA_NET_PATH.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        installable=INTERMEDIATE_DIR + "/capacities_per_country_installable.csv",
        installed=INTERMEDIATE_DIR + "/capacities_per_country_installed.csv",
        storage_installed=INTERMEDIATE_DIR + "/capacities_per_country_storage_installed.csv"
    conda:
        "envs/technologies.yaml"
    script:
        "scripts/aggregate_capacities.py"


rule assemble_technologies:
    input:
        installable=INTERMEDIATE_DIR + "/capacities_per_country_installable.csv",
        installed=INTERMEDIATE_DIR + "/capacities_per_country_installed.csv",
        storage_installed=INTERMEDIATE_DIR + "/capacities_per_country_storage_installed.csv"
    output:
        technologies_dir=directory(TECH_DIR),
        technologies=expand(TECH_DIR + "/technologies_{country}.csv", country=COUNTRIES, allow_missing=True)
    params:
        template=TECH_TEMPLATE
    conda:
        "envs/technologies.yaml"
    script:
        "scripts/assemble_technologies.py"


rule mobility_demands:
    input:
        mobility_eu_countries=lambda wc: MOBILITY_EU_COUNTRIES_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        mobility_eu27=lambda wc: MOBILITY_EU27_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        mobility_non_eu_energy_demands=lambda wc: MOBILITY_NON_EU_ENERGY_DEMANDS_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        mobility_demands=INTERMEDIATE_DIR + "/mobilty_demand_per_country.csv"
    conda:
        "envs/demands.yaml"
    script:
        "scripts/mobility_demand.py"

rule RET_timeseries:
    input:
        network=lambda wc: PYPSA_NET_PATH.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        RET_timeseries_dir=directory(INTERMEDIATE_DIR + "/RET_timeseries_per_country")
    conda:
        "envs/demands.yaml"
    script:
        "scripts/RET_timeseries.py"

rule Demand_timeseries:
    input:
        network=lambda wc: PYPSA_NET_PATH.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        transport_csv=lambda wc: TRANSPORT_DEMAND_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        Demand_timeseries_dir=directory(INTERMEDIATE_DIR + "/Demand_timeseries_per_country")
    conda:
        "envs/demands.yaml"
    script:
        "scripts/Demand_timeseries.py"

rule assemble_timeseries:
    input:
        network=lambda wc: PYPSA_NET_PATH.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        RET_timeseries_dir=INTERMEDIATE_DIR + "/RET_timeseries_per_country",
        Demand_timeseries_dir=INTERMEDIATE_DIR + "/Demand_timeseries_per_country",
    output:
        assembled_dir=directory(TIME_DIR),
        technologies=expand(TIME_DIR + "/time_series_{country}.csv", country=COUNTRIES, allow_missing=True)
    conda:
        "envs/demands.yaml"
    script:
        "scripts/assemble_timeseries.py"

rule energy_and_non_energy_totals:
    input:
        energy_totals=lambda wc: ENERGY_TOTALS_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        industrial_totals=lambda wc: INDUSTRIAL_PRODUCTION_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
    output:
        csv=INTERMEDIATE_DIR + "/energy_and_non_energy_totals_per_country.csv"
    conda:
        "envs/demands.yaml"
    script:
        "scripts/aggregate_energy_reference.py"

rule assemble_demands:
    input:
        Demand_timeseries_dir=INTERMEDIATE_DIR + "/Demand_timeseries_per_country",
        mobility_totals=INTERMEDIATE_DIR + "/mobilty_demand_per_country.csv",
        energy_reference=INTERMEDIATE_DIR + "/energy_and_non_energy_totals_per_country.csv"
    output:
        demands_dir=directory(DEMAND_DIR),
        demands=expand(DEMAND_DIR + "/demands_{country}.csv", country=COUNTRIES, allow_missing=True)
    params:
        template=DEMAND_TEMPLATE
    conda:
        "envs/demands.yaml"
    script:
        "scripts/assemble_demands.py"

rule energyscope_dat:
    input:
        demand=DEMAND_DIR + "/demands_{country}.csv",
        resources=RESOURCES_DIR + "/resources_{country}.csv",
        technologies=TECH_DIR + "/technologies_{country}.csv",
    output:
        datfile=DAT_DIR + "/ESTD_{country}_data_core.dat",
    conda:
        "envs/demands.yaml"
    script:
        "scripts/csv_to_dat.py"


rule typical_days:
    input:
        timeseries=TIME_DIR + "/time_series_{country}.csv",
        mod="conversion_energyscope/core_files/TD_main.mod"
    output:
        td_dat=DAT_DIR + "/ESTD_{country}_12TD.dat"
    params:
        n_TDs=TD_CONFIG.get("n_TDs", 12),
        n_ts=TD_CONFIG.get("n_ts", 9),
        weights=TD_CONFIG.get("weights"),
        print_graphs=TD_CONFIG.get("print_graphs", False),
    conda:
        "envs/demands.yaml"
    script:
        "scripts/TDgenerator.py"

rule hydro_inflow_timeseries:
    input:
        profile=lambda wc: HYDRO_PROFILE_NC.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        )
    output:
        csv=INTERMEDIATE_DIR + "/hydro_inflow_timeseries.csv"
    conda:
        "envs/resources.yaml"
    script:
        "scripts/hydro_inflow_timeseries.py"


rule temperature_cop_timeseries:
    input:
        temp_soil_nc=lambda wc: TEMP_SOIL_NC.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        temp_air_nc=lambda wc: TEMP_AIR_NC.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        pop_csv=lambda wc: POP_CSV.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
        cop_profiles=lambda wc: COP_NC.format(
            year=wc.year,
            clusters=wc.clusters,
            horizon=wc.horizon,
        ),
    output:
        soil_temp_time_series=INTERMEDIATE_DIR + "/soil_temp_time_series.csv",
        air_temp_time_series=INTERMEDIATE_DIR + "/air_temp_time_series.csv",
        cop_time_series_dir=directory(INTERMEDIATE_DIR + "/cop_time_series"),
    conda:
        "envs/resources.yaml"
    script:
        "scripts/temp_cop_timeseries.py"
