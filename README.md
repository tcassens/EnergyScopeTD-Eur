# EnergyScopeTD-Eur

Workflow to generate **EnergyScopeTD** input data (`.dat`) for **European countries**, using **PyPSA-Eur** preprocessing outputs as the main data backbone.

This repo is built as a **Snakemake workflow** that:
1) runs/uses a PyPSA-Eur preprocessing run, then  
2) converts the relevant outputs into **per-country EnergyScopeTD-ready** input files, including a **Typical Days** file.

> **Designed usage:** clone this repo into a PyPSA-Eur checkout and run the workflow from the **PyPSA-Eur root**.

---

## What it produces

For each country in `config.yaml -> countries`, and for each `{year, clusters, horizon}` combination, the workflow writes:

- `EnergyScopeTD-Eur/output/Europe_{year}_{clusters}_{horizon}/dat_files/ESTD_{country}_data_core.dat`
- `EnergyScopeTD-Eur/output/Europe_{year}_{clusters}_{horizon}/dat_files/ESTD_{country}_12TD.dat`

and several intermediate CSVs in:

- `.../intermediate/`
- `.../resources/`
- `.../technologies/`
- `.../demands/`
- `.../timeseries/`

---

## Repository layout

- `Snakefile` – main conversion workflow
- `config.yaml` – conversion settings + input paths + country list
- `config_pypsa_override.yaml` – convenience override config for running PyPSA-Eur preprocessing
- `scripts/` – conversion scripts
- `core_files/` – EnergyScopeTD core tables + `misc.json` required by `.dat` writer
- `envs/` – conda env definitions for Snakemake rules
- `solve_per_country.ipynb` – example notebook to run/inspect country data

---

## Requirements

- A working **PyPSA-Eur** setup (including its data and Snakemake dependencies)
- **Snakemake** (recommended via PyPSA-Eur environment)
- **AMPL + license** for Typical Day generation via `amplpy`  

---
## Quickstart

### 1. Clone the repository

Clone this repository into your **PyPSA-Eur root directory**:

```bash
git clone https://github.com/tcassens/EnergyScopeTD-Eur.git
```

All paths in the workflow assume it is referenced as `EnergyScopeTD-Eur/...`.

---

### 2. Configure the workflow

Edit the configuration files:

```text
EnergyScopeTD-Eur/config.yaml
```
```text
EnergyScopeTD-Eur/config_pypsa_override.yaml
```

**Key settings:**

- `weather_years`
- `clusters`
- `horizon`
- `countries` (ISO country codes)
- `typical_days` settings (number of TDs, weights)

---

### 3. Run PyPSA-Eur preprocessing (if needed)

For a faster and more reliable PyPSA-Eur workflow, you can optionally download the
preprocessed weather data in advance by running:

```bash
python EnergyScopeTD-Eur/scripts/retrieve_weather_data.py <year>

If the required PyPSA-Eur outputs do not yet exist, run:

```bash
snakemake   --configfile EnergyScopeTD-Eur/config_pypsa_override.yaml   --until prepare_sector_networks   --cores <insert_number_of_cores>
```

Make sure the selected year, clusters, and horizon match `config.yaml`.

---

### 4. Run the EnergyScopeTD conversion workflow

From the PyPSA-Eur root directory:

```bash
snakemake   --snakefile EnergyScopeTD-Eur/Snakefile   --cores <insert_number_of_cores>
```

This generates both core `.dat` files and Typical Day files.

---

## Citation and context

This workflow builds on and reuses data from:

- **PyPSA-Eur [1]**
- **EnergyScopeTD [2]**

If you use this workflow in academic work, please cite the corresponding upstream projects.

---

## References

[1] J. Hörsch, F. Hofmann, D. Schlachtberger, and T. Brown.  
*PyPSA-Eur: An Open Optimisation Model of the European Transmission System*.  
arXiv:1806.01613, 2018.  
URL: http://arxiv.org/abs/1806.01613  
GitHub: https://github.com/PyPSA/pypsa-eur

[2] G. Limpens, S. Moret, H. Jeanmart, and F. Maréchal.  
*EnergyScope TD: A novel open-source model for regional energy systems*.  
Applied Energy, 255:113729, 2019.  
DOI: https://doi.org/10.1016/j.apenergy.2019.113729  
URL: https://linkinghub.elsevier.com/retrieve/pii/S0306261919314163

---

## Author

**Till Cassens**  
ETH Zürich

---

## License

MIT License (see `LICENSE`)
