#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License...
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_heatmap_from_softmax_endpoint(datafile, skiprows, skipfooter, startcol, endcol, vmin):
    """

    Parameters
    ----------
    datafile
    skiprows
    skipfooter
    startcol
    endcol

    Returns
    -------

    Usage
    -----
    datafile = r"D:\data\toxr_pipette_test_Brigitte_Thermomax.txt"
    create_heatmap_from_softmax_endpoint(datafile, skiprows=6, skipfooter=103, startcol=2, endcol=14, vmin=0.075)

    """
    plt.style.use('seaborn-whitegrid')
    heatmap_path = datafile[:-4] + "_heatmap.png"
    df = pd.read_csv(datafile, skiprows=skiprows, skipfooter=skipfooter, sep="\t", engine="python")
    df = df.iloc[:, startcol:endcol]
    plt.close("all")
    fig, ax = plt.subplots()
    ax = sns.heatmap(df, ax=ax, vmin=vmin)
    ax.get_figure()
    # rotate the y-ticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    # set the x-axis labels on the top
    ax.xaxis.tick_top()
    fig.savefig(heatmap_path, dpi=240)

def aaa(df_or_series):
    """ Function for use in debugging.
    Saves pandas Series or Dataframes to a user-defined csv file.
    """
     # convert any series to dataframe
    if isinstance(df_or_series, pd.Series):
        df_or_series = df_or_series.to_frame()
    csv_out = r"D:\data\000_aaa_temp_df_out.csv"
    df_or_series.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)