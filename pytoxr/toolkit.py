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
import csv

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


def extract_aa_pos_from_sample_names(df_with_sn_as_index, start_aa, end_aa, newcol_name):
    """Extract the amino acid position from sample names in the index of a dataframe.

    Note this cannot be used if the data contains double-mutants. Consider using a dictionary approach instead.

    Parameters
    ----------
    df_with_sn_as_index : dataframe with the sample name as the index
        Dataframe to be altered, will be returned with the new column
    start_aa : int, starting amino acid in range of amino acids mutated
        UniProt style, NOT python indexing style
    end_aa : int, end amino acid in range of amino acids mutated
        UniProt style where the number is the last aa to be included, NOT python indexing style
    newcol_name : name of the new column containing the index

    Returns
    -------
    df_with_sn_as_index : original dataframe, with the new column
    """

    for i in range(start_aa,end_aa+1):
        # for each sample name in the dataframe with unique samples
        for sn in df_with_sn_as_index.index:
            # if the amino acid number (233, 234, etc) is in the sample name
            if str(i) in sn:
                # add the amino acid number to a new column
                df_with_sn_as_index.loc[sn, newcol_name] = i



def calc_av_over_3aa_window(df, data_col, new_col, n_aa = 4):
    ''' Function to average the properties of an amino acid (conservation, disruption, etc) across a window, consisting
    of three amino acids on the same side of an alpha helix (i, i-4, i+4)
    INPUTS:
    df: dataframe
    data_col: string showing the column name, which contains the data to be averaged
    new_col: string for the new column name, containing the average
    n_aa: number of amino acids separating the aa. Default = 4 (e.g. i, i-4, i+4)
    OUTPUTS:
    The original dataframe, with an extra column.
    '''
    import numpy as np
    # iterate over each row of the dataframe
    for aa_pos in df.index:
        # define the central datapoint
        disrupt_centre = df.loc[aa_pos,data_col]
        # define the "left" N-term datapoint (e.g. i-4)
        if aa_pos-n_aa in df.index:
            disrupt_left = df.loc[aa_pos-n_aa,data_col]
        else:
            disrupt_left = np.nan
        # define the "right" C-term datapoint (e.g. i+4)
        if aa_pos+n_aa in df.index:
            disrupt_right = df.loc[aa_pos+n_aa,data_col]
        else:
            disrupt_right = np.nan
        # calc average
        av_disrupt = np.nanmean([disrupt_centre, disrupt_left, disrupt_right])
        #print(aa_pos, disrupt_centre, disrupt_left, disrupt_right)
        # add average to dataframe
        df.loc[aa_pos,new_col] = av_disrupt
    return df