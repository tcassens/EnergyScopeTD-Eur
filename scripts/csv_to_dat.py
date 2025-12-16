"""Convert EnergyScope CSV outputs plus additional data into AMPL .dat files.

Snakemake calls this script once per country to gather the core EnergyScope
tables, per-country demand/resource/technology CSVs, and shared JSON settings.
It then writes the `.dat` file required by the EnergyScope model by
writing sets and parameters with the expected layout.
"""


##############################################################################################################
# IMPORT PACKAGES
##############################################################################################################

from pathlib import Path
from typing import TYPE_CHECKING

import json
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from snakemake.script import Snakemake


csv_delimiter = ','
dat_space = '	'
CORE_DATA_DIR = Path(__file__).resolve().parents[1] / "core_files"
COUNTRY_FILE_PATTERNS = {
    'demand': ("demands", "demands_{country}.csv"),
    'resources': ("resources", "resources_{country}.csv"),
    'technologies': ("technologies", "technologies_{country}.csv"),
}
CORE_FILE_MAP = {
    'end_use_categories': "END_USES_CATEGORIES.csv",
    'layers_in_out': "Layers_in_out.csv",
    'storage_eff_in': "Storage_eff_in.csv",
    'storage_eff_out': "Storage_eff_out.csv",
    'storage_characteristics': "Storage_characteristics.csv",
}
JSON_FILENAME = "misc.json"
BASE_OUTPUT_TEMPLATE = "conversion_energyscope/output/Europe_{year}_{clusters}_{horizon}"


def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required CSV not found: {path}")
    return pd.read_csv(path, delimiter=csv_delimiter)


def _collect_input_tables(run_dir: Path, country: str, core_dir: Path) -> dict:
    data: dict[str, object] = {}
    for key, filename in CORE_FILE_MAP.items():
        data[key] = _load_csv(core_dir / filename)

    country_code = country.strip().upper()
    for key, (subdir, template) in COUNTRY_FILE_PATTERNS.items():
        file_path = run_dir / subdir / template.format(country=country_code)
        data[key] = _load_csv(file_path)

    json_path = core_dir / JSON_FILENAME
    if not json_path.exists():
        raise FileNotFoundError(f"Required JSON not found: {json_path}")
    with json_path.open() as json_fp:
        data['json'] = json.load(json_fp)

    return data



##############################################################################################################
# FUNCTIONS
##############################################################################################################

def convert_to_dat(run_dir: Path, country: str, output_path: Path, *, core_dir: Path = CORE_DATA_DIR):
    '''
    Convert CSV files + json file to dat file

    Parameters
    ----------
    run_dir : Path
        base directory containing per-country CSV folders (demands, resources, technologies)
    country : str
        two-letter country code identifying which per-country files to read
    output_path : Path
        path to the dat file that will be written
    core_dir : Path, optional
        directory containing shared EnergyScope core CSVs

    Returns
    -------
    None
    '''
    scenario_dir = Path(run_dir)
    core_dir = Path(core_dir)
    country_code = country.strip().upper()
    data = _collect_input_tables(scenario_dir, country_code, core_dir)
    
    # Write dat file
    datfile_path = Path(output_path)
    datfile_path.parent.mkdir(parents=True, exist_ok=True)
    with datfile_path.open('w') as datfile:
        datfile = write_SETS(datfile,data)
        datfile = write_PARAMS(datfile,data)

    return None

##################################################

def write_SET(datfile,setname,set_content):
    '''
    Write a SET in the dat file
    
    Parameters
    ----------
    datfile : file
        dat file to be written
    setname : str
        name of the SET
    set_content : list
        list of elements of the SET
    
    Returns
    -------
    datfile : file
        dat file with SET written
    '''

    datfile.write('set ' + setname + ' :=')
    for i in range(len(set_content)):
        datfile.write(dat_space+str(set_content[i]))
    datfile.write(';\n')

    return datfile

##################################################

def write_SET_entry(datfile,setname,entry,content):
    '''
    Write an entry of a SET in the dat file
    
    Parameters
    ----------
    datfile : file
        dat file to be written
    setname : str
        name of the SET
    entry : str
        name of the entry
    content : list
        list of elements of the SET entry
    
    Returns
    -------
    datfile : file
        dat file with SET entry written
    '''

    datfile.write('set ' + setname +'["'+entry+'"]'+ ' :=')
    for i in range(len(content)):
        datfile.write(dat_space+str(content[i]))
    datfile.write(';\n')

    return datfile

##################################################

def write_SETS(datfile, data):
    ''' 
    Write SETS section of dat file

    Parameters
    ----------
    datfile : file
        dat file to be written
    data : dict
        dictionary containing all dataframes and json file

    Returns
    -------
    datfile : file
        dat file with SETS section written
    '''

    # Write SETS section
    datfile.write('# SETS\n\n')

    # SECTORS
    sectors = list(data['demand'].columns[3:-1])
    datfile = write_SET(datfile,'SECTORS',sectors)

    # END_USES_INPUT
    eui = list(data['demand']['parameter name'])
    datfile = write_SET(datfile,'END_USES_INPUT',eui)

    # END_USES_CATEGORIES
    categories = data['end_use_categories']['END_USES_CATEGORIES'].unique().tolist()
    datfile = write_SET(datfile,'END_USES_CATEGORIES',categories)

    # RESOURCES
    resources = data['resources']['parameter name'].tolist()
    datfile = write_SET(datfile,'RESOURCES',resources)

    # RES_IMPORT_CONSTANT
    ric = data['json']['RES_IMPORT_CONSTANT']
    datfile = write_SET(datfile,'RES_IMPORT_CONSTANT',ric)
    
    # BIOFUELS
    biofuels = data['resources']['parameter name'][data['resources']['Subcategory']=='Biofuel'].tolist()
    datfile = write_SET(datfile,'BIOFUELS',biofuels)

    # RE_RESOURCES
    re_resources = data['resources']['parameter name'][data['resources']['Category']=='Renewable'].tolist()
    datfile = write_SET(datfile,'RE_RESOURCES',re_resources)

    # EXPORT
    export = data['resources']['parameter name'][data['resources']['Subcategory']=='Export'].tolist()
    datfile = write_SET(datfile,'EXPORT',export)

    ###################
    datfile.write('\n')
    ###################

    # END_USE_TYPES_OF_CATEGORIES
    for category in categories:
        entry = data['end_use_categories']['END_USES_TYPES_OF_CATEGORY'][data['end_use_categories']['END_USES_CATEGORIES']==category].tolist()
        datfile = write_SET_entry(datfile,'END_USES_TYPES_OF_CATEGORY',category, entry)

    ###################
    datfile.write('\n')
    ###################

    # TECHNOLOGIES_OF_END_USES_TYPE
    types = data['end_use_categories']['END_USES_TYPES_OF_CATEGORY'].to_list()
    for type in types:
        # Filter: only include conversion technologies, not storage technologies
        # Storage technologies have category: 'Electricity storage', 'Thermal storage', or 'Other storage'
        storage_categories = ['Electricity storage', 'Thermal storage', 'Other storage']
        mask = (data['technologies']['Subcategory']==type) & (~data['technologies']['Category'].isin(storage_categories))
        entry = data['technologies']['Technologies param'][mask].tolist()
        # Only write if entry is not empty
        if entry:
            datfile = write_SET_entry(datfile,'TECHNOLOGIES_OF_END_USES_TYPE',type, entry)
    
    ###################
    datfile.write('\n')
    ###################

    # STORAGE_TECH
    storage_tech = data['storage_characteristics']['param:'].tolist()
    datfile = write_SET(datfile,'STORAGE_TECH',storage_tech)

    # INFRASTRUCTURE
    infrastructure = data['technologies']['Technologies param'][data['technologies']['Subcategory']=='INFRASTRUCTURE'].tolist() 
    datfile = write_SET(datfile,'INFRASTRUCTURE',infrastructure)

    ##################
    datfile.write('\n')
    ###################

    ## EVs
    evs = data['json']['evs']
    # EVS_BATT
    evs_batt = evs['EVs_BATT']
    datfile = write_SET(datfile,'EVs_BATT',evs_batt)
    # V2G
    cars = evs['CAR']
    datfile = write_SET(datfile,'V2G',cars)

     # STORAGE_DAILY
    storage_daily = data['json']['STORAGE_DAILY']
    datfile = write_SET(datfile,'STORAGE_DAILY',storage_daily)

    ##################
    datfile.write('\n')
    ###################

    # STORAGE_OF_END_USES_TYPES
    # Exclude AMMONIA and METHANOL from storage sets
    excluded_types = ['AMMONIA', 'METHANOL']
    for type in types:
        if type not in excluded_types:
            mask=data['storage_eff_in'][type]!=0
            content = data['storage_eff_in']['param:'][mask].tolist()
            if content != []:
                datfile = write_SET_entry(datfile,'STORAGE_OF_END_USES_TYPES',type, content)

    ##################
    datfile.write('\n')
    ###################
    # TS_OF_DEC_TECH
    # Map each decentralized heating tech to its corresponding storage tech using naming convention
    low_t_dec_storages = data['storage_eff_in']['param:'][data['storage_eff_in']['HEAT_LOW_T_DECEN']>0].tolist()
    low_t_dec_techs = data['technologies']['Technologies param'][data['technologies']['Subcategory'] == 'HEAT_LOW_T_DECEN'].tolist()

    for dec_tech in low_t_dec_techs:
        if (dec_tech != 'DEC_SOLAR'): # 'DEC_SOLAR' is not storable
            # Find matching storage tech by extracting the suffix and matching it
            # e.g., DEC_HP_ELEC -> TS_DEC_HP_ELEC
            suffix = dec_tech.replace('DEC_', '')
            matching_storage = f'TS_DEC_{suffix}'
            
            # Only add if the matching storage tech exists in the available storages
            if matching_storage in low_t_dec_storages:
                datfile = write_SET_entry(datfile,'TS_OF_DEC_TECH', dec_tech, [matching_storage])
    
    # EVs_BATT_OF_V2G
    for i,car in enumerate(cars):
        datfile = write_SET_entry(datfile,'EVs_BATT_OF_V2G',car,[evs['EVs_BATT'][i]])
    
    ###################
    datfile.write('\n')
    ###################
    cogen = data['json']['COGEN']
    datfile = write_SET(datfile, 'COGEN', cogen)

    boilers = data['json']['BOILERS']
    datfile = write_SET(datfile, 'BOILERS', boilers)

    ###################
    datfile.write('\n')
    ###################
    return datfile

##################################################

def write_PARAMS(datfile, data):
    '''
    Write PARAMS section of dat file
    
    Parameters
    ----------
    datfile : file
        dat file to be written
    data : dict
        dictionary containing all dataframes and json file
    
    Returns
    -------
    datfile : file
        dat file with PARAMS section written
    '''

    # Write PARAMETERS section
    datfile.write('\n\n# PARAMETERS\n')

    # Write scalar parameters
    for entry in data['json']:
        if type(data['json'][entry])==float or type(data['json'][entry])==int or type(data['json'][entry])==str:
            datfile.write('param ' + entry + ' := ' + str(data['json'][entry]) + ';\n')
        
    datfile.write('\n')

    for entry in data['json']:
        if type(data['json'][entry])==dict:
            keys = list(data['json'][entry].keys())
            if entry == 'evs':
                for key in keys[1:]:
                    if key == 'EVs_BATT':
                        continue
                    
                    datfile.write('param ' + key + ' := \n')
                    for i in range(len(data['json'][entry][key])):
                        datfile.write(data['json'][entry][keys[0]][i] + ' ' + str(data['json'][entry][key][i]) + '\n')
                    datfile.write(';\n\n')
            elif entry == 'loss_network':
                datfile.write('param ' + entry + ' := \n')
                for key in keys:
                    datfile.write(key + ' ' + str(data['json'][entry][key]) + '\n')
                datfile.write(';\n\n')

            else:
                continue

    # Write matrices
    matrices = ['demand','layers_in_out', 'technologies', 'resources', 'storage_eff_in', 'storage_eff_out', 'storage_characteristics']
    drop_rows = [[],[], [],        [], [], [], []]
    drop_cols = [[-1,1,0],[], [-1,2,1,0], [1,0], [], [], []]
    names    = ['end_uses_demand_year','layers_in_out', '', '', 'storage_eff_in', 'storage_eff_out', '']

    for i in range(len(matrices)):
        
        data_matrix = data[matrices[i]]

        for drop_row in drop_rows[i]:
            data_matrix = data_matrix.drop(drop_row)

        for drop_col in drop_cols[i]:
            data_matrix = data_matrix.drop(data_matrix.columns[drop_col], axis=1)

        datfile.write(df_to_param_mat(data_matrix, names[i]))


    return datfile

##################################################

def df_to_scalar_param(datfile,df, name):
    '''
    Convert dataframe to scalar parameter

    Parameters
    ----------
    datfile : file
        dat file to be written
    df : dataframe
        dataframe to be converted
    name : str
        name of scalar parameter
    
    Returns
    -------
    dat : str
        string containing scalar parameter
    '''

    dat = 'param ' + name + ' :=\n'
    for i in df.index:
        dat += str(i) + dat_space + str(df.loc[i, 0]) + '\n'
    dat += ';\n'
    return dat

##################################################

def df_to_param_mat(df, name):
    '''
    Convert dataframe to parameter matrix

    Parameters
    ----------
    datfile : file
        dat file to be written
    df : dataframe
        dataframe to be converted
    name : str
        name of parameter matrix
    
    Returns
    -------
    dat : str
        string containing parameter matrix
    '''

    dat = 'param ' + name + ' : '
    for i in df.columns[1:]:
        dat += str(i) + dat_space
    dat += ':=\n'

    for i in df.index:
        row = df.loc[i].tolist()
        for j in range(len(row)):
            if row[j] == float('inf'):
                row[j] = "Infinity"
            dat += str(row[j])  + dat_space
        dat += '\n'
    dat += ';\n\n'
    return dat

##################################################





##############################################################################################################
# MAIN
##############################################################################################################

def main(snakemake):
    run_dir = Path(BASE_OUTPUT_TEMPLATE.format(**snakemake.wildcards))
    country = snakemake.wildcards.country.upper()
    out_path = Path(snakemake.output.datfile)
    convert_to_dat(run_dir, country, out_path, core_dir=CORE_DATA_DIR)

if __name__ == "__main__":
    main(snakemake)
