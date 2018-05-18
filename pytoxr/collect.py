#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License.
"""
import glob
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys
import numpy as np
from scipy.stats import ttest_ind

def collect_data_in_folder(target_dir):
    """ Collect all ToxR data in a directory, and save in a subdirectory "collected".

    Parameters
    ----------
    target_dir : str
        Directory with .pda, .txt and .xls input files
        Also should contain subdirectories with parsed data from SoftMax Pro to excel/csv tables.

    """
    sys.stdout.write("Starting collect_data_in_folder.\n")
    sys.stdout.flush()

    collect_out_dir = os.path.join(target_dir, "collected")
    if not os.path.isdir(collect_out_dir):
        os.makedirs(collect_out_dir)
    collected_excel = os.path.join(collect_out_dir, "collected.xls")

    # target_dir = r"B:\xbackups\MarkTeese\Praktikum\2017 MEM & MEM PROTEINS\2017_data\old_toxr_results_sep_folders\group3"
    txt_file_list = glob.glob(os.path.join(target_dir, "*.txt"))


    n = 0
    analysed_txt_file_list = []
    exp_name_list = []
    df = pd.DataFrame()
    dfnu_list = []

    for txt_file in txt_file_list:
        print(txt_file)
        filename = os.path.basename(txt_file)
        if "excluded" in filename:
            # skip this file, do not include with other datapoints
            print("{} skipped (has 'excluded' in the filename)".format(filename))
            continue
        pda_file_path = txt_file[:-4] + ".pda"
        samples_file_path = txt_file[:-4] + ".xls"
        exp_name = os.path.basename(txt_file)[:-4][0:60]
        out_dir = os.path.join(os.path.dirname(txt_file), exp_name)
        # for compatibility purposes, currently only excel 2003 is recommended
        excel_format = ".xls"
        out_parsed_excel = os.path.join(out_dir, "{}_parsed{}".format(exp_name, excel_format))
        # check if pda and samples files exist
        if os.path.isfile(pda_file_path) and os.path.isfile(samples_file_path):
            pass
        else:
            # skip this text file
            continue

        # check if parsed excel file exists
        if not os.path.isfile(out_parsed_excel):
            print("{} skipped, no analysed data found".format(exp_name))
            continue

        # open dataframe normalised unique
        dfnu = pd.read_excel(out_parsed_excel, sheetname="norm_unique")
        # rename columns to have a prefix that shows the file name
        newcols = ["{}||{}".format(n, c) for c in dfnu.columns]
        dfnu.columns = newcols
        # add dataframe, filename and experiment names to lists
        dfnu_list.append(dfnu)
        analysed_txt_file_list.append(txt_file)
        exp_name_list.append(exp_name)

        n += 1

    # join all dfnu dataframes
    df = pd.concat(dfnu_list, axis=1)

    """
    Dataframe df should now look like this.

    Input files are labelled 0 to 8, for example. Columns are renamed to input_file_number||col_heading

                     0||order  0||error    0||OD1    0||OD2    0||OD3  ...    \
    P02724-0_GpA_wt    1        NaN        0.893932  0.856831  1.030202  ...
    AZ2                2        NaN        0.747473  0.753499  0.569780  ...
    P02724-0_GpA_G83A  3        NaN        0.878750  0.798697  1.018703  ...
    dTM                4        NaN        0.772738  1.039832  1.467613  ...
    Q6ZRP7-2_QSOX2_wt  5        NaN        1.000000  1.000000  1.000000  ...

                        8||MU2    8||MU3  8||MU_mean  8||MU_std  8||Vm
    P02724-0_GpA_wt    0.718398  0.820212  0.743628    0.417787   Vm01
    AZ2                1.070686 -1.716484  0.440827    16.463324  Vm17
    P02724-0_GpA_G83A  0.375756  0.337701  0.359586    0.538124   Vm25
    dTM               -0.029046  0.027206 -0.017234    0.347562   Vm02
    Q6ZRP7-2_QSOX2_wt  1.000000  1.000000  1.000000    1.000000   NaN
    """

    # copy across the sample order from the original samples file. since some norm-unique data is missing samples,
    # the orig sample order is most likely the max of the available sample orders
    order_cols = ["{}||order".format(n) for n in range(n)]
    df_order = df[order_cols].copy()
    order_ser = df_order.max(axis=1)
    df["order"] = order_ser
    # sort by the sample order, for the following mean dataframes and figures
    df.sort_values("order", inplace=True)

    ttest_standard_cols = ["{}||ttest_standard".format(n) for n in range(n)]
    #ttest_standard_columns = set(ttest_standard_cols).intersection(set(df.index))
    ttest_standard_columns = set(ttest_standard_cols).intersection(set(df.columns))

    # WEIRD NON-REPEATABLE DF.T EFFECT TROUBLESHOOTING
    # print("ttest_standard_cols", ttest_standard_cols)
    # print("ttest_standard_columns", ttest_standard_columns)
    # print("df.index", df.index)
    # print("df.columns", df.columns)

    # check that there are actually some ttest standards labelled
    if ttest_standard_columns != set():
        conduct_ttest = True
        df_ttest_standard = df[ttest_standard_cols].copy()
        df_ttest_standard = df_ttest_standard.dropna(how="all")
        if df_ttest_standard.shape[0] == 1:
            ttest_standard = df_ttest_standard.index[0]
        else:
            raise ValueError("seems to be more than one ttest_standard selected in the original sample excel files."
                             "df_ttest_standard.index = {}".format(df_ttest_standard.index))
    else:
        conduct_ttest = False

    data_type_list = ["{}||OD_mean", "{}||Vi_mean", "{}||MU_mean"]
    writer = pd.ExcelWriter(collected_excel)

    data_is_normalised = True
    dfs = pd.read_excel(samples_file_path)
    standard_subset = dfs.loc[dfs.standard == True]
    standard_sample_name = standard_subset["sample"].unique()[0]
    if data_is_normalised:
        norm_string = ", normalised to {}".format(standard_sample_name)
    else:
        norm_string = ""
    titles_dict = {"OD" : "OD600, all experiments", "Vi" : "Initial Velocity (Vi), all experiments", "MU" : "Miller Units, all experiments"}
    y_axis_label_dict = {"OD": "OD600{}".format(norm_string), "Vi": "Initial Velocity (Vi, U/min){}".format(norm_string), "MU": "Miller Units (U/min/OD600){}".format(norm_string)}


    for data_type in data_type_list:

        excel_tab_name = data_type.split("||")[1].split("_")[0]
        title = titles_dict[excel_tab_name]
        y_axis_label = y_axis_label_dict[excel_tab_name]

        scatter_path = os.path.join(collect_out_dir, "{}_scatter_plot.png".format(excel_tab_name))
        box_path = os.path.join(collect_out_dir, "{}_box_plot.png".format(excel_tab_name))

        selected_cols = [data_type.format(n) for n in range(n)]
        df_sel = df[selected_cols].copy()
        df_sel.columns = exp_name_list
        df_sel.to_excel(writer, sheet_name=excel_tab_name)

        """
        df_sel looks like this, and only shows the normalised miller units

        2010915_prak_samples  group1_day1  group1_day2  \
        P02724-0_GpA_wt    0.814460              0.951131     1.001340
        AZ2                2.377500              2.795199     5.627452
        P02724-0_GpA_G83A  0.233870              0.262983     0.389499
        dTM                0.473288              0.048931     0.027387
        Q6ZRP7-2_QSOX2_wt  1.000000              1.000000     1.000000
        """

        plt.close("all")
        ax = df_sel.plot(marker="o", linestyle="", alpha=0.5, figsize=(8, 12))
        ax.set_xticks(range(df_sel.shape[0]))
        ax.set_xticklabels(df_sel.index, rotation=90)
        #ax.set_ylim(0, 6)
        ax.set_xlim(-1, df_sel.shape[0] + 1)
        ax.set_title(title)
        ax.set_ylabel(y_axis_label)
        plt.legend()
        plt.tight_layout()
        plt.savefig(scatter_path, dpi=240)
        plt.savefig(scatter_path[:-4] + ".pdf")
        plt.close("all")

        # pandas describe function doesn't seem to have an axis, so this clumsy method is used to flip columns and index twice
        df_sel_describe = df_sel.T.describe().T
        #df_sel_describe["std"] = df_sel_describe["std"].fillna(0)
        df_sel_describe["sem"] = df_sel_describe["std"] / np.sqrt(df_sel_describe["count"])
        df_sel_describe["order"] = order_ser

        print("conduct_ttest = {}".format(conduct_ttest))

        sample_list = df_sel.index.tolist()

        for sample_name in sample_list:
            data = df_sel.loc[sample_name, :].dropna()
            # get the median. Note that the median of [1,2,3,40] is 2.5.
            df_sel_describe.loc[sample_name, "median"] = data.median()

            if conduct_ttest and sample_name != ttest_standard:
                # get the standard "wildtype" data for the TTEST
                ttest_standard_data = df_sel.loc[ttest_standard, :]
                #sample_list.remove(ttest_standard)

                if len(data) > 1 and len(ttest_standard_data) > 1:
                    p_value = ttest_ind(ttest_standard_data, data)[1]
                else:
                    p_value = np.nan
                df_sel_describe.loc[sample_name, "ttest_p_value"] = p_value

            df_sel_describe["sign_stars"] = df_sel_describe["ttest_p_value"].apply(return_p_value_stars)
            df_sel_describe["ttest_significant"] = df_sel_describe["ttest_p_value"] < 0.05
        df_sel_describe.to_excel(writer, sheet_name=excel_tab_name + "_describe")

        max_value = df_sel_describe["max"].max()
        min_value = df_sel_describe["min"].max()

        """
        df_describe now looks like this

                           count      mean       std       min       25%       50%       75%       max       sem  order    median  ttest_p_value sign_stars  ttest_significant
        P02724-0_GpA_wt     14.0  1.000000  0.000000  1.000000  1.000000  1.000000  1.000000  1.000000  0.000000    1.0  1.000000       0.000073        ***               True
        AZ2                 13.0  1.163821  0.314529  0.964444  1.037625  1.065475  1.173914  2.167819  0.087235    2.0  1.065475       0.918356        NaN              False
        P02724-0_GpA_G83A   12.0  1.246863  0.466768  0.840600  1.052239  1.132847  1.236926  2.672585  0.134744    3.0  1.132847       0.578728        NaN              False
        dTM                 13.0  1.372324  0.462053  1.015609  1.067898  1.100609  1.624869  2.466696  0.128151    4.0  1.100609       0.135864        NaN              False
        Q6ZRP7-2_QSOX2_wt   14.0  1.173372  0.137769  0.926538  1.134087  1.155864  1.232545  1.403153  0.036820    5.0  1.155864            NaN        NaN              False

        """
        # plot barchart with two averaging methods, mean and median
        for stat_method in ["mean", "median"]:

            bar_path = os.path.join(collect_out_dir, "{}_{}_bar_plot.png".format(excel_tab_name, stat_method))
            plt.close("all")
            # fig = plt.Figure(figsize=(5,10))
            ax = df_sel_describe[stat_method].plot(kind="bar", figsize=(8, 12), yerr=df_sel_describe['sem'])
            #
            # sign_ = df_sel_describe["sign_stars"].copy()
            # sign_ser.index = range(sign_ser.shape[0])
            # sign_ser.dropna(inplace=True)


            if conduct_ttest:
                df_sel_describe["#"] = range(1, df_sel_describe.shape[0] + 1)
                sign_df = df_sel_describe.dropna(subset=["sign_stars"])

                for sample in sign_df.index:
                    x = sign_df.loc[sample, "#"] - 1
                    y = sign_df.loc[sample, stat_method]
                    sign_star = sign_df.loc[sample, "sign_stars"]
                    ax.annotate(sign_star, [x - 0.18, y + 0.005], horizontalalignment='center')
            ax.set_title(title)
            ax.set_ylabel(y_axis_label)
            #ax.errorbar(range(df_sel_describe.shape[0]), df_sel_describe['mean'], yerr=df_sel_describe['sem'])

            plt.tight_layout()
            plt.savefig(bar_path, dpi=240)
            plt.savefig(bar_path[:-4] + ".pdf")
            plt.close("all")

            meanpointprops = dict(marker='o', markerfacecolor='black', markersize=2)  # markeredgecolor='0.75',
            fig, ax = plt.subplots(figsize=(7, 10))
            df_sel.T.plot(kind="box", whis=1500, showmeans=True, meanprops=meanpointprops, showfliers=True, ax=ax)
            ax.set_xticklabels(df_sel.index, rotation=90)
            ax.set_title(title)
            ax.set_ylabel(y_axis_label)

            if conduct_ttest:
                for sample in sign_df.index:
                    x = sign_df.loc[sample, "#"]
                    y = sign_df.loc[sample, "75%"]
                    sign_star = sign_df.loc[sample, "sign_stars"]
                    ax.annotate(sign_star, [x - 0.18, y + 0.005], horizontalalignment='center')
            if data_is_normalised:
                if max_value > 6:
                    ylim_min = ax.get_ylim()[0]
                    ax.set_ylim(ylim_min,6)
            plt.tight_layout()
            plt.savefig(box_path, dpi=240)
            plt.savefig(box_path[:-4] + ".pdf")
            plt.close("all")

    sys.stdout.write("\n\ncollect_data_in_folder is finished.\Data from {} files collected in total.\n-----------------------------------------------------------------\n".format(n+1))
    sys.stdout.flush()

    writer.close()

def return_p_value_stars(p):
    if float(p) < 0.001:
        return "***"
    elif float(p) < 0.01:
        return "**"
    elif float(p) < 0.05:
        return "*"
    else:
        return np.nan