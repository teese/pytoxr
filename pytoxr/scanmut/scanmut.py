import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import ast
import eccpy.tools as eccpytools
import eccpy as eccpy
from pytoxr.toolkit import extract_aa_pos_from_sample_names
from pytoxr.fit_to_sine.fit_to_sine import fit_sin_to_data

def run_scanmut(settings_file, show_barchart=False, save_pdf=True):
    """ Runs the scanmut algorithm on the selected input files in the excel settings file.

    Scanmut takes typical scanning mutagenesis data, creates a barchart showing the values for each amino acid,
    and if desired, fits the data to a sine curve to judge if the data follows alpha-helical periodicity.

    The correct excel settings file template is necessary. Typically, the template is saved in the folder containing the
    scripts in the scanmut puthon module.

    Parameters
    ----------
    settings_file : excel file
        Path to excel file containing the scanmut settings.
            "settings" tab:
                List of datafiles to be processed, individual parameters associated with each file
                Each datafile should at LEAST contain a sample name, measured value, and optionally a y-error value.
                Datafiles can also contain columns with the original amino acid, position, and final amino acid. If
                this info is not available directly, you can get it from the dictionary on the settings excel file, or
                by automatic recognition of the mutation, based on the sample name.
            "dictionary" tab:
                 Columns forming a dictionary that links long sample names, e.g. insulin_precursor_Y108C, to separate
                 columns containing the relevant information ("original amino acid" = Y, "amino acid position" = 108,
                "final amino acid" = c.
            "explanations" tab:
                Explanation of elements in the excel file. Not used in the python script.

    show_barchart : boolean
        If True, the output figure will be shown as a popup window, or in IPython/Jupyter.
        For Ipython/Jupyter it is recommended to precede scanmut with the magic command %matplotlib inline.

    Saved Files and Figures
    -------
    All output files are saved in a subfolder based on the location of the input xlsx file:
    D:\Path\To\Your\Input\SettingsFile\scanmut_output\

    barchart : png, pdf
        Shows individual scanning mutagenesis data, clustered for each amino acid position, with bars coloured by
        the final amino acid.
    fit_sin_to_data_output : png, pdf, csv
        see tlabassays.sin_disrupt_fit.fit_sin_to_data for details regarding the fit-to-sine output files

    Usage
    -------
    import tlabtools
    settings_file = r"D:\Path\To\Your\Input\SettingsFile\scanmut_output\scanmut_settings.xlsx"
    tlabtools.run_scanmut(settings_file)
    """
    # open the settings as a pandas dataframe
    df_settings = pd.read_excel(settings_file, sheet_name = "settings")
    # fill any empty values (np.nan) with empty strings
    df_settings.fillna("", inplace=True)
    # open the dictionary showing the samplename and associated aa position as a pandas dataframe, dfd
    dfd = pd.read_excel(settings_file, sheetname = "dictionary", index_col="sample name")

    # normalise infile path to suit the operating system
    df_settings.loc[:,"infile"] = df_settings.loc[:,"infile"].apply(lambda x: os.path.normpath(x))

    for row in df_settings[df_settings.run_scanmut == True].index:
        sys.stdout.write("\n{}\n".format(df_settings.loc[row,"short experiment title"]))

        ###############################################################################################
        #                                                                                             #
        #         Read input file with scanmut data, define columns with data, etc                    #
        #                                                                                             #
        ###############################################################################################

        # define the path to the input file
        input_datafile = df_settings.loc[row,"infile"]
        # set up the output folder as the folder with the settings file, and a chosen subfolder
        outpath, outfile = eccpy.settings.setup_output_folder(settings_file, "scanmut_output")
        # obtain experiment title
        exp_title = df_settings.loc[row, "short experiment title"]
        outfile_base = outfile + "_" + exp_title
        # open the input datafile containing the scanning mutagenesis data
        if input_datafile[-4:] == ".csv":
            df = pd.read_csv(input_datafile, index_col = 0)
        elif input_datafile[-4:] == "xlsx":
            df = pd.read_excel(input_datafile, index_col=0)
        elif input_datafile[-4:] == ".xls":
            df = pd.read_excel(input_datafile, index_col=0)
        # raise an error if the dataframe is empty, due to a faulty excel file
        if df.empty:
            raise ValueError("The input file does not seem to contain data. Please check the filename, sheet_name,"
                             "and data contents.\nFile affected:\n{}".format(input_datafile))
        # set the index according to user input
        AAA_array = np.array(df_settings)

        if df_settings.loc[row, "index column"] != "":
            if df_settings.loc[row, "index column"] in df.columns:
                df.set_index(df_settings.loc[row, "index column"], drop=False, inplace=True)
        # define column with the y-axis data
        y_col = df_settings.loc[row,"y-data column"]
        # define column with the error bar data
        yerr_col = df_settings.loc[row, "y-error column"]
        # define final amino acid types (e.g. P in F255P, A256P) to be ignored
        list_ignored_final_aa = df_settings.loc[row,"ignored final amino acids"]
        # define the name of the wildtype
        wildtype_sample_name = df_settings.loc[row, "wildtype sample name"]

        ###############################################################################################
        #                                                                                             #
        #      Extract the amino acid positions, either from sample names or other columns            #
        #                                                                                             #
        ###############################################################################################

        # extract the amino acid position from the sample names
        if df_settings.loc[row,"amino-acid position identification method"] == "automatic":
            start_aa_pos = df_settings.loc[row,"start position"]
            end_aa_pos = df_settings.loc[row,"end position"]
            extract_aa_pos_from_sample_names(df, start_aa=start_aa_pos, end_aa=end_aa_pos, newcol_name="aa_pos")
        elif df_settings.loc[row,"amino-acid position identification method"] == "dictionary":
            # use the dictionary to extract the relevant amino acid positions from the sample names
            df = pd.concat([df, dfd.loc[df.index, :]], axis=1)
        elif "column_dictionary" in df_settings.loc[row,"amino-acid position identification method"]:
            column_dictionary = ast.literal_eval(df_settings.loc[row, "amino-acid position identification method"].split("=")[1])
            # if the dictionary contains the final columns as keys, reverse keys and values
            if "original amino acid" in column_dictionary.keys():
                column_dictionary = {v: k for k, v in column_dictionary.items()}
            # if all the columns in the dictionary are found in the dataframe, use dictionary to rename the columns
            if set(column_dictionary.keys()).issubset(set(df.columns)):
                df.rename(columns=column_dictionary, inplace=True)

        ###############################################################################################
        #                                                                                             #
        #      Setup wildtype, drop empty rows, create list of data and gap positions                 #
        #                                                                                             #
        ###############################################################################################

        # set the wildtype data as the first column to be shown
        if wildtype_sample_name != "":
            # set the amino acid position of the wildtype as 0, so that it sorts correctly
            df.loc[wildtype_sample_name, "amino acid position"] = df["amino acid position"].min() - 1
            df.loc[wildtype_sample_name, "original amino acid"] = "wt"
            df.loc[wildtype_sample_name, "final amino acid"] = "wt"

        # drop any of the columns that do not have a particular amino acid position (for example, positive control samples)
        df.dropna(subset=["amino acid position"], inplace=True)
        # drop any columns with mutations to amino acids that are excluded
        for final_aa in list_ignored_final_aa:
            df = df.loc[df['final amino acid'] != final_aa]

        # confirm that all of the amino acid positions are integers
        df["amino acid position"] = df["amino acid position"].astype(int)
        # sort the samples by the amino acid position in the dataframe
        df.sort_values(by=["amino acid position", "final amino acid"], inplace=True)
        # create a list of unique amino acid positions (for which there could be multiple mutations)
        list_unique_pos = df["amino acid position"].dropna().unique()
        # find the min and max of the amino acid positions
        min_aa_pos = list_unique_pos.min()
        max_aa_pos = list_unique_pos.max()
        # use min and max to create a contiguous list of positions, including any gaps
        list_contiguous_aa_pos = list(range(min_aa_pos,max_aa_pos+1))
        # create a list of gaps
        list_aa_pos_gaps = list(set(list_contiguous_aa_pos) - set(list_unique_pos))
        list_aa_pos_gaps.sort()
        if list_aa_pos_gaps != []:
            sys.stdout.write((", gaps (positions without data): " + str(list_aa_pos_gaps)))

        # set the index as the amino acid position
        df.set_index("amino acid position", inplace=True, drop=False)

        # double-check for errors, where the orig and final amino acid are the same
        errorseries = df["original amino acid"] == df["final amino acid"]
        if True in errorseries.value_counts():
            # collect the rows with the error
            error_rows = df.loc[errorseries].loc[:, ["original amino acid","final amino acid"]]
            # check if there is a wildtype row
            wtrow = error_rows[error_rows['original amino acid'] == "wt"].index.tolist()
            # drop the wildtype row from the error list
            error_rows.drop(wtrow, inplace=True)
            if not error_rows.empty:
                print("Warning, possible error found. Mutant listed as having wildtype sequence.\n""File:\n{}\n\nSample:"
                      "\n{}\n\n".format(input_datafile, error_rows))

        ###############################################################################################
        #                                                                                             #
        #  Move data to dfc_mult, which contains average values for all mutations, for fit to sine    #
        #                                                                                             #
        ###############################################################################################

        # create a new dataframe to store the average values for all mutations
        dfc_mult = pd.DataFrame()
        # create an empty list to hold the positions with no data (no mutants)
        for pos in list_contiguous_aa_pos:
            if pos in df.index:
                # select the row (or rows) with mutations at that position
                selected_data_at_aa_pos = df.loc[pos,:]
                # if selected_data_at_aa_pos is a series, there is only one mutation at that position.
                if isinstance(selected_data_at_aa_pos,pd.Series):
                    dfc_mult.loc[pos,"mean_mult"] = selected_data_at_aa_pos[y_col]
                    dfc_mult.loc[pos,"std_mult"] = 0
                # if selected_data_at_aa_pos is a dataframe, there are multiple mutations at that position. Take the mean.
                elif isinstance(selected_data_at_aa_pos,pd.DataFrame):
                    dfc_mult.loc[pos,"mean_mult"] = selected_data_at_aa_pos[y_col].mean()
                    dfc_mult.loc[pos,"std_mult"] = selected_data_at_aa_pos[y_col].std()
            else:
                # the position is not in the index, a gap is counted, the values will automatically be np.nan
                dfc_mult.loc[pos,"mean_mult"] = np.nan
                dfc_mult.loc[pos,"std_mult"] = 0

        if df_settings.loc[row, "run_fit_to_sin"] == True:
            #####################################################################################################
            #                                                                                                   #
            #                                        RUN FIT TO SINE                                            #
            #                                                                                                   #
            #####################################################################################################
            # define the x-axis for the sine calculation as the contiguous index integer
            x_sin = dfc_mult.index
            # the y-axis is the mean of all mutations at that position (note, NOT the mean of all replicates)
            y_sin = eccpytools.normalise_0_1(dfc_mult["mean_mult"])[0]
            title = "Scan. Mut. {}".format(exp_title[:25])
            output_tuple = fit_sin_to_data(x_sin, y_sin, title, outpath, outfile_base + "_sin")

            title_shuf = "Scan. Mut. {}, {}".format(exp_title[:25], "randomised")
            y_sin_shuffled = y_sin.copy()
            np.random.shuffle(y_sin_shuffled)
            output_tuple = fit_sin_to_data(x_sin, y_sin_shuffled, title_shuf, outpath, outfile_base  + "_sin" + "_randomised")

            print("Run sindisrupt is finished, files saved to {}".format(outpath))

        #####################################################################################################
        #                                                                                                   #
        #               CREATE BAR CHART WITH MUTATIONS TO MULTIPLE AMINO ACIDS                             #
        #                                                                                                   #
        #####################################################################################################

        # create a list of the final amino acids (ie., orig aa mutated to final aa)
        list_final_aa = df["final amino acid"].unique().tolist()
        # create empty list to hold the original amino acid
        list_orig_aa = []
        # create empty list to hold orig pos
        # list_orig_pos = []
        fontsize = 6
        plt.rcParams['font.size'] = fontsize
        # define width of bars in barchart, which will depend on number of bars per aa position
        max_n_aa_per_position = df["amino acid position"].value_counts().max()
        bar_width = 0.8 / max_n_aa_per_position
        # create figure and axes objects with matplotlib
        plt.close("all")
        fig_bar, ax_bar = plt.subplots()
        # create a list of colours, starting with black if the wildtype is included
        t20 = eccpy.tools.setup_t20_colour_list()
        # find the mutation number for that position, to place the bars next to each other in the barchart
        print(len(df.index.unique()))
        print(df.shape[0])
        for aa_pos in list_contiguous_aa_pos:
            if aa_pos in list_aa_pos_gaps:
                # NOTE some problem exists here with non-unique positions
                if len(df.index.unique()) == df.shape[0]:
                    df.loc[aa_pos,"mut_number_at_that_pos"] = 0
                    df.loc[aa_pos,"mut_offset_from_centre"] = 0
                    list_orig_aa.append("ND")
            # else:
            #     # add to a list of the actual positions with data
            #     list_orig_pos.append(aa_pos)
            else:
                # select data at that amino acid position, resulting in a dataframe or series (if only 1 mutant)
                df_series_aa_pos = df.loc[aa_pos,:]
                # if it's a series, there is only one mutant
                if isinstance(df_series_aa_pos, pd.Series):
                    df.loc[aa_pos,"mut_number_at_that_pos"] = 0
                    df.loc[aa_pos,"mut_offset_from_centre"] = 0
                    # select orig aa, add to list
                    list_orig_aa.append(df.loc[aa_pos,"original amino acid"])
                # if it's a dataframe, there are multiple mutants at that position
                elif isinstance(df_series_aa_pos, pd.DataFrame):
                    # count the mutants
                    n_mut = df_series_aa_pos.shape[0]
                    # insert the mutant number as a range of the number of mutants
                    df.loc[aa_pos,"mut_number_at_that_pos"] = range(n_mut)
                    offset = df.loc[aa_pos,"mut_number_at_that_pos"].mean()
                    df.loc[aa_pos,"mut_offset_from_centre"] = df.loc[aa_pos,"mut_number_at_that_pos"] - offset
                    list_orig_aa.append(df.loc[aa_pos,"original amino acid"].tolist()[0])

        # change the wildtype colour to grey, if it is in the dataset
        if list_orig_aa[0] == "wt":
            t20 = ["0.5"] + t20

        # define the general position of the bars on the x-axis, starting with 0
        df["barchart_index"] = df["amino acid position"] - min_aa_pos
        # define the exact position of each bar, based on the mutant number and calculated offset
        df["barchart_index_with_offset"] = df["barchart_index"] + df["mut_offset_from_centre"]*bar_width

        three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
        'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
        'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
        'G':'GLY', 'P':'PRO', 'C':'CYS'}

        # iterate through each final aa
        for n, aa in enumerate(list_final_aa):
            # set colour for those bars
            colour = t20[n]
            # select the data with that final aa
            df_final_aa = df.loc[df["final amino acid"] == aa]
            # obtain the position of the bars on the x-axis
            x_index = df_final_aa["barchart_index_with_offset"].tolist()
            # obtain the y-values and y-err
            y_values = df_final_aa[y_col].tolist()
            if yerr_col != "":
                y_err = df_final_aa[yerr_col].tolist()
            else:
                y_err = 0
            # define the legend-label
            if aa == "wt":
                legend_label = None
            else:
                legend_label = three_letter[aa]
            # plot the data
            bar_object = ax_bar.bar(x_index, y_values, color=colour, width=bar_width,
                                    yerr=y_err, error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1),
                                   label=legend_label, align="center")

        # obtain the number of datapoints on the x-axis
        n_dp = len(list_contiguous_aa_pos)
        # set the xticks
        xticks = np.array(range(n_dp))
        ax_bar.set_xticks(xticks)
        # create a series of the amino acid positions for the x-axis labels. Add the wt, if it is available.
        if list_orig_aa[0] == "wt":
            index = ["wt"] + list_contiguous_aa_pos[1:]
        else:
            index = list_contiguous_aa_pos
        # create a new dataframe to hold the xticklabels. use the list_contiguous_aa_pos as the index
        df_xticklabel = pd.DataFrame({"aa_pos": index,"list_orig_aa":list_orig_aa}, index=index).astype("U11")
        # for n, pos in enumerate(list_unique_pos):
        #     if pos in df_xticklabel.index:
        #         df_xticklabel.loc[pos, "aa"] = list_orig_aa[n]
        # df_xticklabel.fillna("", inplace=True)
        # df_xticklabel.index = df_xticklabel.index.fillna("")
        # df_xticklabel["xlabels"] = df_xticklabel.index  # + df_xticklabel["aa"]
        df_xticklabel["xlabels"] = df_xticklabel["aa_pos"].astype("U11") + "\n" + df_xticklabel["list_orig_aa"]
        ax_bar.set_xticklabels(df_xticklabel["xlabels"])
        # set the x-axis limits to be comfortable
        ax_bar.set_xlim(-bar_width*4, n_dp + bar_width*1)
        ax_bar.set_xlim(-1, n_dp + bar_width * 1)
        # set the y-axis limits as 5% above and below the max and min values, including error bars
        yerr_max = df[yerr_col].max() if yerr_col != "" else 0
        # extract the y-axis minimum or max from the settings file. Otherwise, pad ~10% around the min and max datapoints.
        if df_settings.loc[row,"ymax"] != "":
            ymax = df_settings.loc[row,"ymax"]
        else:
            ymax = (df[y_col].max() + yerr_max) * 1.1
        if df_settings.loc[row,"ymin"] != "":
            ymin = df_settings.loc[row,"ymin"]
        else:
            ymin = (df[y_col].min() - yerr_max) * 0.9
        ax_bar.set_ylim(ymin, ymax)
        # set the legend across the top (number of columns = number of final amino acids).
        ax_bar.legend(borderaxespad=2, ncol=len(list_final_aa))
        # set the y-axis title
        ax_bar.set_ylabel(df_settings.loc[row,"y-axis label"])
        # set title
        ax_bar.set_title(exp_title)
        fig_basename = os.path.join(outpath, outfile)
        fig_bar.savefig(fig_basename + "_" + exp_title + "_scanmut_barchart.png", format='png', dpi=150)
        if save_pdf == True:
            fig_bar.savefig(fig_basename + "_" + exp_title + "_scanmut_barchart.pdf", format='pdf')
        if show_barchart == True:
            plt.show()
        if df_settings.loc[row, "run_fit_to_sin"] != True:
            sys.stdout.write(", finished")

    print("Scanning mutagenesis analysis is finished.\nFiles saved to {}".format(outpath))
