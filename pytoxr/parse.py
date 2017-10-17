#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License.
"""
import glob
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import sys

def parse_softmax(txt_path, samples_path):
    """Parse the text file output from Sofmax Pro, using the standard A600 endpoint & A405 kinetic template.
    This method uses the initial velocities calculated by Softmax Pro.

    Parameters
    ----------
    txt_path : str
        Path to txt data file.
    samples_path : str
        Path to an excel file with all of the sample names.
        An example of this text file should be in the "examples" subfolder.

    Dataframes
    ----------
    df = raw softmax pro output
    dfd = dataframe for data.
        index = Vm#, columns = Vi1      Vi2     Vi3       MU1       MU2, etc
    dfnu = dataframe for normalised, unique data
        index = sample_name, columns = as above for dfd

    Returns
    -------
    dfd : pd.DataFrame
        Dataframe for data.
    """
    exp_name = os.path.basename(txt_path)[:-4][0:60]
    out_dir = os.path.join(os.path.dirname(txt_path), exp_name)
    # path to output heatmap of raw OD600 values
    out_OD600_heatmap = os.path.join(out_dir, "{}_OD600_heatmap.png".format(exp_name))
    # path to excel output file
    # for compatibility purposes, currently only excel 2003 is recommended
    excel_format = ".xls"
    out_parsed_excel = os.path.join(out_dir, "{}_parsed{}".format(exp_name, excel_format))
    # only the non-normalised data is saved as csv
    out_parsed_csv = os.path.join(out_dir, "{}_parsed.csv".format(exp_name))

    # create output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


    """
    FORMAT OF THE OD600 section of the SOFTMAX PRO FILE

~End
Plate:	OD600	1,3	PlateFormat	Endpoint	Absorbance	Raw	TRUE	1						1	600	1	12	96	1	8	None
	Temperature(¡C)	1	2	3	4	5	6	7	8	9	10	11	12
	24,10	0,1094	0,1322	0,121	0,117	0,1214	0,1239	0,1128	0,1219	0,1191	0,1152	0,1172	0,1164
		0,153	0,1564	0,1582	0,1518	0,1636	0,1528	0,1448	0,1651	0,1371	0,1484	0,1491	0,1509
		0,1194	0,1218	0,1266	0,12	0,1171	0,1252	0,1155	0,1227	0,123	0,1204	0,1221	0,1159
		0,1217	0,1237	0,1239	0,119	0,1217	0,1245	0,1168	0,1241	0,1207	0,1168	0,1203	0,1203
		0,1152	0,119	0,1402	0,1184	0,1443	0,1219	0,1193	0,1254	0,1206	0,1173	0,1167	0,1165
		0,1253	0,1313	0,1435	0,1232	0,1261	0,1298	0,1239	0,1315	0,1133	0,1193	0,1157	0,1178
		0,1136	0,1143	0,1359	0,1172	0,1373	0,1275	0,1159	0,1281	0,1224	0,1195	0,1168	0,1143
		0,1078	0,1206	0,1243	0,1139	0,1199	0,1229	0,1172	0,121	0,1206	0,0379	0,0382	0,0407
    """

    # go through each line, searching for the various components
    regex_search_string = "Plate:\sOD600"
    # regex_search_string = "Group:\s+Results Kinetic\s+1"
    OD600_plate_is_found = False
    OD600_plate_line = 0
    first_line_with_OD600_data_int = None
    first_line_with_OD600_data = None

    with open(txt_path, "r") as f:
        for n, line in enumerate(f):
            match = re.match(regex_search_string, line)
            if match:
                OD600_plate_is_found = True
                OD600_plate_line = n
                first_line_with_OD600_data_int = OD600_plate_line + 2
            if n == first_line_with_OD600_data_int:
                first_line_with_OD600_data = line
    total_n_lines = n

    table_start = first_line_with_OD600_data_int
    table_end = first_line_with_OD600_data_int + 8
    skipfooter = total_n_lines - table_end + 1

    if "," in first_line_with_OD600_data:
        decimal = ","
    elif "." in first_line_with_OD600_data:
        decimal = "."
    df = pd.read_table(txt_path, skiprows=table_start, sep='\s+', decimal=decimal, header=None, skipfooter=skipfooter, engine="python")

    """
    OD600 data now looks something like this.
    The temperature in row 0 puts everything out of order.

            0       1       2       3       4       5       6       7       8        9       10      11      12
    0  24.1000  0.1094  0.1322  0.1210  0.1170  0.1214  0.1239  0.1128  0.1219   0.1191  0.1152  0.1172  0.116
    1   0.1530  0.1564  0.1582  0.1518  0.1636  0.1528  0.1448  0.1651  0.1371   0.1484  0.1491  0.1509     Na
    2   0.1194  0.1218  0.1266  0.1200  0.1171  0.1252  0.1155  0.1227  0.1230   0.1204  0.1221  0.1159     Na
    3   0.1217  0.1237  0.1239  0.1190  0.1217  0.1245  0.1168  0.1241  0.1207   0.1168  0.1203  0.1203     Na
    4   0.1152  0.1190  0.1402  0.1184  0.1443  0.1219  0.1193  0.1254  0.1206   0.1173  0.1167  0.1165     Na
    5   0.1253  0.1313  0.1435  0.1232  0.1261  0.1298  0.1239  0.1315  0.1133   0.1193  0.1157  0.1178     Na
    6   0.1136  0.1143  0.1359  0.1172  0.1373  0.1275  0.1159  0.1281  0.1224
    7   0.1078  0.1206  0.1243  0.1139  0.1199  0.1229  0.1172  0.1210  0.1206   0.1195  0.1168  0.1143     Na
    """
    # fix up alignment in dataframe, rename columns and index
    df.loc[0, :] = list(df.loc[0, 1:]) + [0]
    df = df.loc[:, 0:11]
    df.index = list("ABCDEFGH")
    df.columns = range(1, 13)

    """
    OD600 data is ready for heatmap

               1       2       3       4       5       6       7       8       9   \
    A  0.1094  0.1322  0.1210  0.1170  0.1214  0.1239  0.1128  0.1219  0.1191
    B  0.1530  0.1564  0.1582  0.1518  0.1636  0.1528  0.1448  0.1651  0.1371
    C  0.1194  0.1218  0.1266  0.1200  0.1171  0.1252  0.1155  0.1227  0.1230
    D  0.1217  0.1237  0.1239  0.1190  0.1217  0.1245  0.1168  0.1241  0.1207
    E  0.1152  0.1190  0.1402  0.1184  0.1443  0.1219  0.1193  0.1254  0.1206
    F  0.1253  0.1313  0.1435  0.1232  0.1261  0.1298  0.1239  0.1315  0.1133
    G  0.1136  0.1143  0.1359  0.1172  0.1373  0.1275  0.1159  0.1281  0.1224
    H  0.1078  0.1206  0.1243  0.1139  0.1199  0.1229  0.1172  0.1210  0.1206

    """

    ax = sns.heatmap(df)
    fig = ax.get_figure()
    # rotate the y-ticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    # set the x-axis labels on the top
    ax.xaxis.tick_top()
    fig.savefig(out_OD600_heatmap, dpi=240 )

    """
    OD600 heatmap is finished.
    Now grab the kinetic data, which is presented in a table at the bottom of the text file.

    This is the FORMAT OF THE INPUT SOFTMAX PRO FILE, at region of processed kinetic data

    ..... <----------- raw kinetic data is above. we currently use the sofmax pro calculated initial velocity
    ~End
    Group:	Results Kinetic	1                                <----------- group_results_line. Use regex to find the line number associated with this string.
    Sample	Wells	Sample#	V/max	Vmax/OD600	Mean	SD
    Vm01	A1	1	43,605	1231,785	1379,232	132,621  <----------- first_line_with_data. Use this to determine if . or , is the decimal.
        A2		59,253	1488,774
        A3		50,733	1417,135
    Vm02	B1	2	26,652	870,965	801,480	75,053
        B2		39,200	811,590
        B3		36,816	721,885

    ...

    ~End The last line with ~End shows the location of the end of the file. This is variable. Best to use dropna() to drop empty rows at the end of the dataframe.
    """

    # go through each line, searching for the various components
    regex_search_string = "Group:\s+Results Kinetic\s+1"
    group_results_are_found = False
    group_results_line = 0
    first_line_with_data_int = None
    first_line_with_data = None
    start_of_table = None

    with open(txt_path, "r") as f:
        for n, line in enumerate(f):
            match = re.match(regex_search_string, line)
            if match:
                group_results_are_found = True
                group_results_line = n
                first_line_with_data_int = group_results_line + 2
                start_of_table = group_results_line+1
            if n == start_of_table:
                start_of_table_data = line
            if n == first_line_with_data_int:
                first_line_with_data = line

    if not group_results_are_found:
        raise ValueError("Softmax template is not recognised. Check input file names, and that 'Group: Results Kinetic' is located in the text file.")

    # check if the data is saved using a German or English decimal point
    if "," in first_line_with_data:
        dec = ","
    elif "." in first_line_with_data:
        dec = "."
    else:
        raise TypeError("VersaMax exported text file seems to have an unknown decimal format. Try re-exporting data.")

    """
    pandas versions seem to use a different skiprows
    
    need to check that the skiprows is correctly aligned
    
    Group:	Results Kinetic	1
    Sample	Wells	Sample#	V/max	Vmax/OD600	Mean	SD
    Vm01	A1	1	55,264	790,996	776,806	48,225
    """

    if "Sample" in start_of_table_data:
        skiprows = group_results_line
    elif "Group" in start_of_table_data:
        skiprows = group_results_line + 1
    else:
        raise ValueError("skiprows does not work. check txt file format.")
    # open csv as a pandas dataframe, starting from "Sample	Wells	Sample#	V/max	Vmax/OD600	Mean" ...etc
    df = pd.read_csv(txt_path, sep='\t', skiprows=skiprows, decimal=dec)
    # drop the last footer rows, which may vary in length
    df.dropna(subset=["V/max"], inplace=True)
    # fill the names downwards, so that they can be used for creating pivot tables later
    df.fillna(method="pad", inplace=True)

    # search for datapoints that could not be fitted
    if "NoFit" in df["V/max"].tolist():
        df_no_fit = df.loc[df["V/max"] == "NoFit"].copy()
        df_no_fit["Vm_Well"] = df_no_fit["Sample"] + "_well_" + df_no_fit["Wells"]
        sample_names_with_no_fit = df_no_fit["Sample"].unique()
        samples_with_no_fit = df_no_fit["Vm_Well"].dropna().unique()
        sys.stdout.write("'NoFit' found in data, indicating measurement was impossible, probably due to very high enzyme concentrations.\n"
                         "Samples affected: {}".format(samples_with_no_fit))

        df.replace("NoFit", np.nan, inplace=True)
        df["V/max"] = df["V/max"].astype(float)
        df["Vmax/OD600"] = df["Vmax/OD600"].astype(float)

        #df["V/max"] = df["V/max"].replace("NoFit", np.nan).astype(float)
        #df["Vmax/OD600"] = df["Vmax/OD600"].replace("NoFit", np.nan).astype(float)
        # drow any rows with no fit
        #df = df.loc[df["V/max"].notnull()]
    else:
        sample_names_with_no_fit = []

    # calculate OD from V and V/OD (not currently supplied by the template)
    df["OD"] = df["V/max"] / df["Vmax/OD600"]


    """
    df should look like this. Basically, it looks like the SoftMax pro output using Jan Kirrbach's template.

      Sample Wells  Sample#   V/max  Vmax/OD600      Mean       SD        OD
    0   Vm01    A1      1.0  43.605    1231.785  1379.232  132.621  0.035400
    1    NaN    A2      NaN  59.253    1488.774       NaN      NaN  0.039800
    2    NaN    A3      NaN  50.733    1417.135       NaN      NaN  0.035800
    3   Vm02    B1      2.0  26.652     870.965   801.480   75.053  0.030601
    4    NaN    B2      NaN  39.200     811.590       NaN      NaN  0.048300
    """

    Vi_replicates = ["Vi1","Vi2","Vi3"] * int(df.shape[0]/3)
    MU_replicates = ["MU1","MU2","MU3"] * int(df.shape[0]/3)
    OD_replicates = ["OD1","OD2","OD3"] * int(df.shape[0]/3)
    df["Vi_rep"] = Vi_replicates
    df["MU_rep"] = MU_replicates
    df["OD_rep"] = OD_replicates

    dfv = df.pivot_table(values="V/max", index = "Sample", columns="Vi_rep")
    dfd = df.pivot_table(values="Vmax/OD600", index = "Sample", columns="MU_rep")
    dfo = df.pivot_table(values="OD", index = "Sample", columns="OD_rep")
    dfd = pd.concat([dfv, dfd, dfo], axis=1)

    """
    dfd now contains the values for each replicate, one after another

               Vi1      Vi2     Vi3       MU1       MU2       MU3       OD1     OD2     OD3
    Sample
    Vm01    43.605   59.253  50.733  1231.785  1488.774  1417.135  0.035400  0.0398  0.0358
    Vm02    26.652   39.200  36.816   870.965   811.590   721.885  0.030601  0.0483  0.0510
    Vm03    36.813   66.878  56.552  1191.375  1706.078  1516.140  0.030900  0.0392  0.0373
    Vm04    42.418   46.726  45.454  1258.694  1115.172  1014.602  0.033700  0.0419  0.0448
    Vm05    70.692  101.895  88.205  2380.192  3078.387  2183.284  0.029700  0.0331  0.0404

    Vi is the initial velocity.
    MU is the miller units
    OD is the OD600 (cell density)

    now the mean and standard deviation is simply recalculated for the Vi, MU and OD replicates
    """

    dfd["Vi_mean"] = dfd.loc[:, "Vi1":"Vi3"].mean(axis=1)
    dfd["Vi_std"] = dfd.loc[:, "Vi1":"Vi3"].std(axis=1)
    dfd["MU_mean"] = dfd.loc[:, "MU1":"MU3"].mean(axis=1)
    dfd["MU_std"] = dfd.loc[:, "MU1":"MU3"].std(axis=1)
    dfd["OD_mean"] = dfd.loc[:, "OD1":"OD3"].mean(axis=1)
    dfd["OD_std"] = dfd.loc[:, "OD1":"OD3"].std(axis=1)


    # open excel file with sample names
    dfs = pd.read_excel(samples_path, index_col=0)

    dfs["frame"] = dfs.frame.dropna().astype(int).astype(str)
    dfs.fillna("", inplace=True)
    for i in dfs.index:
        dfs.loc[i, "sample_python"] = dfs.loc[i, "uniprot":"mut"].dropna().astype(str).str.cat()
        # try:
        #     dfs.loc[i, "sample_python"] = dfs.loc[i, "uniprot":"mut"].dropna().astype(str).str.cat()
        # except:
        #     print(i)
        #     print(dfs.loc[i, "uniprot":"mut"])
        #     print(dfs.loc[i, "uniprot":"mut"].dropna())
        #     dfs.loc[i, "sample_python"] = dfs.loc[i, "uniprot":"mut"].dropna().astype(str).str.cat()
    if dfs["sample"].tolist() != dfs.sample_python.tolist():
        raise ValueError("sample column in excel does not exactly match data in preceding columns. Try repeating fill-down.")
    sample_name_dict = dfs.set_index("Vm#")["sample"].to_dict()
    order_dict = dfs.set_index("Vm#")["order"].to_dict()
    include_data_dict = dfs.set_index("Vm#")["include_data"].to_dict()

    """
    dfs should now look like this
        Vm# uniprot  - frame  _ shortname __  mut              sample      notes    include_data       sample_python
    #                                                                             #
    1  Vm01  P02724  -     0  _       GpA  _   wt     P02724-0_GpA_wt        GpA  1         True     P02724-0_GpA_wt
    2  Vm02     ΔTM                                               ΔTM        ΔTM  2         True                 ΔTM
    3  Vm03  Q6ZRP7  -     2  _     QSOX2  _  L4A  Q6ZRP7-2_QSOX2_L4A  QSOX2_L4A  3         True  Q6ZRP7-2_QSOX2_L4A
    4  Vm04  Q6ZRP7  -     2  _     QSOX2  _  S8A  Q6ZRP7-2_QSOX2_S8A  QSOX2_S8A  4         True  Q6ZRP7-2_QSOX2_S8A
    5  Vm05  Q6ZRP7  -     2  _     QSOX2  _  S8Q  Q6ZRP7-2_QSOX2_S8Q  QSOX2_S8Q  5         True  Q6ZRP7-2_QSOX2_S8Q
    """

    dfd["Vm"] = dfd.index
    dfd["sample"] = dfd["Vm"].replace(sample_name_dict)
    dfd["order"] = dfd["Vm"].replace(order_dict)
    dfd["include_data"] = dfd["Vm"].replace(include_data_dict)

    # note that there has been an error in some samples
    dfd.loc[sample_names_with_no_fit, "error"] = "NoFit"
    dfd["error"] = dfd["error"].fillna("none")

    # change the column order
    col_order = ['sample', 'order','include_data', 'error', 'OD1', 'OD2', 'OD3','OD_mean', 'OD_std','Vi1', 'Vi2', 'Vi3', 'Vi_mean', 'Vi_std', 'MU1', 'MU2', 'MU3', 'MU_mean', 'MU_std', 'Vm']
    dfd = dfd.reindex(columns=col_order)
    dfd.sort_values("order", inplace=True)

    dfn = pd.DataFrame()

    if True in dfs.standard.tolist():
        standard_subset = dfs.loc[dfs.standard == True]
        standard_sample_name_list = standard_subset.sample_python.unique()
        if len(standard_sample_name_list) != 1:
            raise ValueError("Multiple standards are selected. Double-check excel file with sample names.")
        standard_sample_name = standard_sample_name_list[0]
        standard_Vm_numbers = standard_subset.index.tolist()

        standard_values_ser_or_df = dfd.loc[dfd["sample"] == standard_sample_name]
        if isinstance(standard_values_ser_or_df, pd.DataFrame):
            standard_mean_values_ser = standard_values_ser_or_df.mean(axis=0)
        elif isinstance(standard_values_ser_or_df, pd.Series):
            standard_mean_values_ser = standard_values_ser_or_df
        else:
            raise TypeError("Oops! 'standard_values_ser_or_df' is neither a dataframe nor a series")

        for i in range(dfd.shape[0]):
            row_series = dfd.iloc[i, :]
            Vm = dfd.index[i]
            dfn[Vm] = row_series / standard_mean_values_ser

        dfn = dfn.T
        for transferred_col in ["sample", "order", "Vm", "include_data"]:
            dfn[transferred_col] = dfd[transferred_col]
        dfn = dfn.reindex(columns=col_order)
        dfn.sort_values("order", inplace=True)

        """
        dfn still has the Vm numbers as the index, which allows duplicates (e.g. two wells of QSOX2_wt in this experiment)

                          sample  order       OD1       OD2       OD3   OD_mean
        Vm01     P02724-0_GpA_wt      1  0.801497  0.853056  0.767382    0.8074
        Vm17                 AZ2      2  0.880439  0.935936  0.531282  0.779956
        Vm25   P02724-0_GpA_G83A      3  0.185782   0.16348   1.00166  0.457173
        Vm02                 dTM      4   1.20262   1.11343   1.10826   1.13996
        Vm09   Q6ZRP7-2_QSOX2_wt      5   1.05084  0.940833  0.935106  0.973782
        Vm24   Q6ZRP7-2_QSOX2_wt      6  0.949156   1.05917   1.06489   1.02622
        Vm10  Q6ZRP7-2_QSOX2_C1S      7   1.42985   1.56866     1.648   1.55175
        Vm18  Q6ZRP7-2_QSOX2_V2A      8   1.40206   1.50719   1.54405   1.48644
        Vm26  Q6ZRP7-2_QSOX2_V3A      9   1.11713  0.116117   1.39081  0.869624
        Vm03  Q6ZRP7-2_QSOX2_L4A     10  0.931424  0.765291  0.812374  0.834106
        """

        unique_samples = dfn["sample"].unique()
        dfn.set_index("sample", inplace=True)

        # exclude datapoints according to manual annotation in the "include_data" column
        dfn_include_data = dfn.loc[dfn.include_data]

        if dfn_include_data.empty:
            # all data has been labelled for exclusion here
            sys.stdout.write()
            sys.stdout.write("\n'{}' skipped, all data labeled as FALSE for 'include_data' in samples excel file.".format(exp_name))
            # skip this protein
            return dfd

        dfnu = pd.DataFrame()

        for s in unique_samples:
            if s not in dfn_include_data.index:
                # skip this unique sample, as it was manually excluded via the "include_data" column
                continue
            row_df_or_ser = dfn_include_data.loc[s, :]
            if isinstance(row_df_or_ser, pd.Series):
                dfnu.loc[:, s] = row_df_or_ser
            elif isinstance(row_df_or_ser, pd.DataFrame):
                # take the order of the first duplicate row as the correct one
                order = row_df_or_ser.iloc[0, :]["order"]
                # drop the text rows, so the dtype is changed to float
                # otherwise the mean of two values for the standard  (which should give 1.0) returns np.nan
                # reason for strange behaviour in pandas is unknown
                row_df_or_ser = row_df_or_ser.drop(["Vm", "order"], axis=1)
                mean_of_duplicate_samples = row_df_or_ser.mean(axis=0)
                # add the single row to the unique dataframe
                dfnu.loc[:, s] = mean_of_duplicate_samples
                dfnu.loc["order", s] = order

        dfnu = dfnu.T

        dropped_cols = ["include_data", "error"]
        dfnu.drop(dropped_cols, axis=1, inplace=True)

        ttest_standard_array = dfs.loc[dfs.ttest_standard == True].sample_python.unique()
        if len(ttest_standard_array) == 1:
            ttest_standard = ttest_standard_array[0]
            if ttest_standard in dfnu.index:
                dfnu.loc[ttest_standard, "ttest_standard"] = True
        elif len(ttest_standard_array) > 1:
            raise ValueError("seems to be two different ttest standards selected")
        elif len(ttest_standard_array) == 0:
            print("no ttest_standard selected")

    """
    dfnu is now completely normalised to your standard (e.g. QSOX2).
    This includes all OD measurements, initial velocity, and Miller Units.
    The index is now unique. The values of duplicates show the mean for that experiment.
    The values of duplicate standards (e.g. Q6ZRP7-2_QSOX2_wt) are now 1.0.

                          order       OD1       OD2       OD3   OD_mean     OD_      Vi1        Vi2        Vi3    Vi_mean      Vi_std
    P02724-0_GpA_wt       1  0.928313  0.898545  0.556443     0.748  0.0632076  0.594462   0.683314   0.632995   0.637554    0.669004
    AZ2                   2  0.709755  0.673359  0.401063  0.556855          0  0.563122   0.530232   0.466218   0.516645           0
    P02724-0_GpA_G83A     3   5.40776   5.13045   3.05578   4.24278          0  0.778933   0.733438   0.644892   0.714644           0
    dTM                   4   3.29763   2.30533  0.788788   1.86296    2.09499 0.0715175  0.0668064  0.0597716  0.0656486  0.00467747
    Q6ZRP7-2_QSOX2_wt     5         1         1         1         1          1         1          1          1          1           1
    """

    dfd.index.name = "Vm"
    dfd.to_csv(out_parsed_csv)

    with pd.ExcelWriter(out_parsed_excel) as writer:
        dfd.to_excel(writer, sheet_name="original")
        dfn.to_excel(writer, sheet_name="normalised")
        dfnu.to_excel(writer, sheet_name="norm_unique")
        writer.close()

    sys.stdout.write("\n'{}' parsed successfully".format(exp_name))
    sys.stdout.flush()

    return dfd

def parse_all_data_files_in_folder(target_dir, reparse_existing=False):
    """Parse all ToxR SoftMax Pro txt files in a given directory.

    Parameters
    ----------
    target_dir : str
        Directory with .pda, .txt and .xls input files

    """
    txt_file_list = glob.glob(os.path.join(target_dir, "*.txt"))

    #print(txt_file_list)
    counter = 0

    for txt_file in txt_file_list:
        #print(txt_file)
        pda_file_path = txt_file[:-4] + ".pda"
        samples_file_path = txt_file[:-4] + ".xls"
        exp_name = os.path.basename(txt_file)[:-4][0:60]
        out_dir = os.path.join(os.path.dirname(txt_file), exp_name)
        excel_format = ".xls"
        out_parsed_excel = os.path.join(out_dir, "{}_parsed{}".format(exp_name, excel_format))


        if os.path.isfile(pda_file_path) and os.path.isfile(samples_file_path):
            if not reparse_existing:
                if os.path.isfile(out_parsed_excel):
                    # skip analysis of this file
                    sys.stdout.write("\n{} skipped, parsed output file already exists.".format(exp_name))
                    sys.stdout.flush()
                    continue

            dfd = parse_softmax(txt_file, samples_file_path)
            counter += 1

    sys.stdout.write("\n\nparse_all_data_files_in_folder is finished.\n{} file(s) parsed in total.\n-----------------------------------------------------------------\n".format(counter))
    sys.stdout.flush()
