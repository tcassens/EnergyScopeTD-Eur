# Typical Day (TD) generator based on yearly data series
# Need to additionally install openpyxl

import sys
import os
from pathlib import Path
from typing import TYPE_CHECKING
from amplpy import AMPL, DataFrame
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Still ToDo:
# - add clustering errors checks as in Gauthier's paper
# - fix error in mod file (square root in formula)
# - add output plots


def run_td_generator(
    fn_in: Path,
    fn_out: Path,
    mod_file: Path,
    *,
    n_TDs: int = 30,
    n_ts: int = 15,
    weights=None,
    printGraphs: bool = True,
):
    """
    Execute the original TD logic with configurable inputs/outputs.
    All arguments mirror the legacy top-level variables.
    """
    fn_in = Path(fn_in).resolve()
    fn_out = Path(fn_out).resolve()
    mod_file = Path(mod_file).resolve()

    base_path = Path(os.path.dirname(os.path.abspath(__file__)))
    weights = weights or [1, 1, 1, 1, 1, 1, 1, 1, 1]

    n_days = 365
    n_hours = 24

    # Read time series (CSV input)
    df_ts = pd.read_csv(fn_in)
    for col in ['t', '{PERIODS}', 'PERIODS']:
        if col in df_ts.columns:
            df_ts = df_ts.drop(columns=[col])
    print(df_ts)

    ts_names_all = df_ts.columns.to_numpy()
    n_ts_all = ts_names_all.shape[0]
    ts_names = ts_names_all[0:n_ts]
    n_dimensions = n_ts * n_hours

    df_ts_norm = df_ts.iloc[:, 0:n_ts]
    ts_sums = df_ts_norm.sum(axis=0)
    ts_sums_safe = ts_sums.where(ts_sums != 0, 1)
    df_ts_norm = df_ts_norm / ts_sums_safe
    df_ts_norm = df_ts_norm * weights

    print('-----------------------')
    print('There are ' + str(n_ts_all) + ' time series in total. ' + str(n_ts) + ' of them are used for the clustering')
    print('This is their sum before normalization and weighing:')
    print(df_ts.sum(axis=0))
    print('This is their sum after normalization and weighing:')
    print(df_ts_norm.sum(axis=0))
    print('-----------------------')

    ts_dict = {}
    index = 0
    for ts_name in ts_names_all:
        ts_dict[ts_name] = {}
        values_temp = np.zeros((n_days, n_hours))
        for d in range(n_days):
            for h in range(n_hours):
                values_temp[d][h] = df_ts[ts_name].iloc[n_hours * d + h]

        ts_dict[ts_name]['Values'] = values_temp
        ts_dict[ts_name]['Values_sum'] = np.sum(ts_dict[ts_name]['Values'])
        if index < n_ts:
            ts_dict[ts_name]['Weights'] = weights[index]
        index += 1
        ts_dict[ts_name]['Duration_curve'] = np.sort(df_ts[ts_name])[::-1]

    os.chdir(os.path.dirname(__file__) or os.curdir)

    try:
        ampl = AMPL()
        ampl.read(str(mod_file))

        Nbr_TD = ampl.getParameter('Nbr_TD')
        Nbr_TD.setValues([n_TDs])
        print('The number of TDs is ' + str(Nbr_TD.value()))
        print('-----------------------')
        print('Writing AMPL problem and solving it. Hold on...')

        df = DataFrame('DAYS')
        days = [str(x) for x in range(1, n_days + 1)]
        df._setColumn('DAYS', days)
        ampl.setData(df, 'DAYS')

        df = DataFrame('DIMENSIONS')
        dimensions = [str(x) for x in range(1, n_dimensions + 1)]
        df._setColumn('DIMENSIONS', dimensions)
        ampl.setData(df, 'DIMENSIONS')

        df = DataFrame(('DAYS', 'DIMENSIONS'), 'Ndata')

        Ndata = np.zeros((n_days, n_dimensions))
        for d in range(n_days):
            for ts in range(n_ts):
                for h in range(n_hours):
                    r = d
                    c = n_hours * ts + h
                    Ndata[r][c] = df_ts_norm.iloc[n_hours * d + h, ts]

        df._setValues({
            (day, dimension): Ndata[i][j]
            for i, day in enumerate(days)
            for j, dimension in enumerate(dimensions)
        })
        ampl.setData(df)

        ampl.setOption('solver', 'cplex')
        ampl.setOption('cplex_options', 'outlev=1 timelimit=3600')

        ampl.solve()

        print('Objective: {}'.format(ampl.getObjective('Euclidean_distance').value()))

        Selected_TD_df = ampl.getVariable('Selected_TD').getValues().toPandas()
        Selected_TD_df = Selected_TD_df[(Selected_TD_df['Selected_TD.val']) > 0]
        tds_selected = []
        for i, row in Selected_TD_df.iterrows():
            tds_selected.append(int(i))

        tds_selected.sort()

        print('-----------------------')
        print('The selected TDs are: ')
        print(tds_selected)
        print('-----------------------')

        Cluster_matrix_df = ampl.getVariable('Cluster_matrix').getValues().toPandas()
        Cluster_matrix_df = Cluster_matrix_df[(Cluster_matrix_df['Cluster_matrix.val']) > 0]

        d_td_dict = {}
        for i, row in Cluster_matrix_df.iterrows():
            d_temp = int(i[1])
            td_temp = int(i[0])
            d_td_dict[d_temp] = td_temp

        print(d_td_dict)

        TD_dict = {}
        index = 1
        for td in tds_selected:
            TD_dict[td] = {}
            TD_dict[td]['Name'] = index
            index += 1
            for ts_name in ts_names_all:
                TD_dict[td][ts_name] = ts_dict[ts_name]['Values'][td - 1, :]
            TD_dict[td]['Days'] = []
            for d in days:
                if d_td_dict[int(d)] == td:
                    TD_dict[td]['Days'].append(int(d))
            TD_dict[td]['Number of Days'] = len(TD_dict[td]['Days'])

        counter = 0
        for td in tds_selected:
            counter += TD_dict[td]['Number of Days']
        assert (counter == 365)

        scaling_factors = {}
        for ts_name in ts_names_all:
            ts_sum_orig = ts_dict[ts_name]['Values_sum']
            ts_sum_new = 0
            for td in tds_selected:
                ts_sum_new += np.sum(TD_dict[td][ts_name]) * TD_dict[td]['Number of Days']
            scaling_factors[ts_name] = ts_sum_new / ts_sum_orig

        for ts_name in ts_names_all:
            print(ts_name + ' time series scaling factor: ' + str(scaling_factors[ts_name]))
        print('-----------------------')

        for td in tds_selected:
            for ts_name in ts_names_all:
                scaled = TD_dict[td][ts_name] / scaling_factors[ts_name]
                TD_dict[td][ts_name + '_scaled'] = np.nan_to_num(scaled, nan=0.0)

        for ts_name in ts_names_all:
            ts_sum_new = 0
            for td in tds_selected:
                ts_sum_new += np.sum(TD_dict[td][ts_name + '_scaled']) * TD_dict[td]['Number of Days']
            if ts_sum_new - ts_dict[ts_name]['Values_sum'] < 1e-6:
                print(ts_name + ' time series scaled correctly! Sum after scaling: ' + str(round(ts_sum_new, 2)))

        print('-----------------------')

        for ts_name in ts_names_all:
            values_temp = np.zeros((n_days * n_hours))
            i = 0
            for d in days:
                td = d_td_dict[int(d)]
                for h in range(n_hours):
                    values_temp[i] = TD_dict[td][ts_name + '_scaled'][h]
                    i += 1
            ts_dict[ts_name]['Values_TDs'] = values_temp

        fn_out.parent.mkdir(parents=True, exist_ok=True)
        f = open(fn_out, 'w')

        f.write('set T_H_TD :=\n')
        hours = [str(x) for x in range(1, n_hours + 1)]

        for d in days:
            for h in hours:
                t = n_hours * (int(d) - 1) + int(h)
                td = d_td_dict[int(d)]
                f.write('(\t' + str(t) + '\t,\t' + h + '\t,\t' + str(TD_dict[td]['Name']) + '\t)\n')

        f.write(';\n\n')

        ts_lookup = {name.lower(): name for name in ts_names_all}

        def require_series(label: str) -> str:
            ts_name = ts_lookup.get(label.lower())
            if ts_name is None:
                raise ValueError(f"Required time series '{label}' is missing from {ts_names_all}.")
            return ts_name

        def write_hour_td_param(param_name: str, ts_label: str) -> None:
            ts_name = require_series(ts_label)
            f.write(f'param {param_name} :\t')
            for td in range(1, n_TDs + 1):
                f.write(str(td) + '\t')
            f.write(':=\n')
            for h in hours:
                f.write(h + '\t')
                for td in tds_selected:
                    temp = round(TD_dict[td][ts_name + '_scaled'][int(h) - 1], 10)
                    f.write(str(temp) + '\t')
                f.write('\n')
            f.write(';\n\n')

        TIMESERIES_MAPPINGS = [
            ('electricity_time_series', 'Electricity (%_elec)'),
            ('heating_time_series', 'Space Heating (%_sh)'),
            ('mob_pass_time_series', 'Passanger mobility (%_pass)'),
            ('mob_freight_time_series', 'Freight mobility (%_freight)'),
        ]

        for param_name, series_label in TIMESERIES_MAPPINGS:
            write_hour_td_param(param_name, series_label)

        f.write('param c_p_t :=\n')

        def write_cpt_block(tech_label: str, ts_label: str) -> None:
            ts_name = require_series(ts_label)
            f.write(f'["{tech_label}",*,*]:\t')
            for td in range(1, n_TDs + 1):
                f.write(str(td) + '\t')
            f.write(':=\n')
            for h in hours:
                f.write(h + '\t')
                for td in tds_selected:
                    temp = round(TD_dict[td][ts_name + '_scaled'][int(h) - 1], 10)
                    f.write(str(temp) + '\t')
                f.write('\n')
            f.write('\n')

        CPT_MAPPINGS = [
            ('PV', 'PV'),
            ('WIND_ONSHORE', 'Wind_onshore'),
            ('WIND_OFFSHORE', 'Wind_offshore'),
            ('HYDRO_RIVER', 'Hydro_river'),
            ('DHN_SOLAR', 'Solar'),
            ('DEC_SOLAR', 'Solar'),
        ]

        for tech_label, series_label in CPT_MAPPINGS:
            write_cpt_block(tech_label, series_label)

        f.write(';\n\n')
        f.close()

        if printGraphs:
            if not os.path.exists(base_path / '../../Figures/TDs'):
                os.makedirs(base_path / '../../Figures/TDs')
            for ts_name in ts_names_all:
                x_data = [h for h in range(1, n_hours * n_days + 1)]
                y_data = ts_dict[ts_name]['Duration_curve']
                plt.plot(x_data, y_data, label='Reference')
                y2_data = np.sort(ts_dict[ts_name]['Values_TDs'])[::-1]
                plt.plot(x_data, y2_data, label='TDs')
                plt.title(ts_name + ' Duration Curve')
                plt.xlabel('Hours')
                plt.legend()
                plt.savefig(base_path / f'../../Figures/TDs/{ts_name}.pdf')
                plt.close()

    except Exception as e:
        print(e)
        raise


def main(snakemake: "Snakemake") -> None:
    """Snakemake entry point mirroring other workflow scripts."""

    run_td_generator(
        Path(snakemake.input.timeseries),
        Path(snakemake.output.td_dat),
        Path(snakemake.input.mod),
        n_TDs=int(snakemake.params.n_TDs),
        n_ts=int(snakemake.params.n_ts),
        weights=snakemake.params.weights,
        printGraphs=snakemake.params.print_graphs,
    )


if __name__ == "__main__":
    main(snakemake)
