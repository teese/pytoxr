#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         Mark Teese
Created:        06.11.2015
Dependencies:   Python 3.x
                Numpy
                SciPy
                Pandas
Purpose:        Fit a sine wave to experimental data
Credits:        All sections by Mark Teese.
"""
from __future__ import unicode_literals
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import os
import pandas as pd
from korbinian.utils import create_colour_lists
import pytoxr.mathfunctions as mathfunctions

def fit_sin_to_data(x, y, title, output_folder, output_base_filename, fit_unbiased=True, fit_alphahelix=True,x_grid=True,savefig=True, savecsv=True, saveexcel=True):
    """ Fits a sinusoidal curve to data.

    1) fits orig data to unbiased sine curve
    2) fits orig data to sine curve biased for alpha helical periodicity of 3.6

    The output files are decided by the user, and can include a png figure, excel and/or csv files.

    Parameters
    ----------
    x : arraylike (e.g. numpy 1D array, pandas Series)
        Data of the x axis (e.g. aa position)
    y : arraylike (e.g. numpy 1D array, pandas Series)
        Array-like data of the y axis (e.g. disruption)
    title : string
        Title of the output figure
    output_folder : string
        Output path. E.g. r"D:\Data\fit_to_sine_output\"
    output_base_filename : string
        Start of filename. E.g. "20160607_experiment1"
    fit_unbiased : bool
        If True, the sinusoidal curve will be fit to the data with an unbiased periodicity.
    fit_alphahelix : bool
        If True, the sinusoidal curve will befit to the data with alpha-helical periodicity of 3.6 (constant b = 1.7)
    savefig : bool
        If True, output figures will be saved
    savecsv : bool
        If True, output csv files will be saved
    saveexcel : bool
        If True, output excel files will be saved

    Saved Files and Figures
    -----------------------
    basename.xlsx : excel file with the following tabs
        sine_constants : sine constants from fitted curves
        yvalues_fitted : x and y datapoints used to draw the curves
    basename_sine_constants.csv : csv
        sine constants from fitted curves
    basename_yvalues_fitted.csv : csv
        x and y datapoints used to draw the curves
    basename.png : line graph of original input data,
        fitted unbiased sine curve, and fitted sine curve with fixed periodicity

    Usage Example
    -------------
    import tlabtools as tools
    import pandas as pd
    '''
    LaPointe 2013 ToxR Ala scanning data of FtsB that was fitted to a sine curve in publication
    '''
    lit_data_xls = r"D:\PyCharmProjects\tlabtools\tlabtools\fit_to_sine\examples\P0A6S5.xlsx"
    output_folder = r"D:\DATABASES_spitfire\20160607_fit_to_sine_output"
    output_base_filename = "20160607_P0A6S5"
    # read data file
    dflit = pd.read_excel(lit_data_xls, sheetname="Sheet1")
    # define x
    ftsb_x = dflit.aa_position
    # define y (normalised between 0 and 1)
    ftsb_y = tools.normalise_0_1(dflit["Disruption(1-mutation)"])[0]
    # set title
    P0A6S5_title = "P0A6S5, FtsB (Senes Lab)"
    # run fit sin to data
    tools.fit_sin_to_data(ftsb_x,ftsb_y,P0A6S5_title,output_folder,output_base_filename,savefig=True)
    # REPEAT WITH RANDOMIZED/SHUFFLED DATA
    ftsb_y_shuffled = ftsb_y.copy()
    np.random.shuffle(ftsb_y_shuffled)
    output_base_filename_random = "20160607_P0A6S5_RANDOM"
    P0A6S5_title_random = "P0A6S5, FtsB (Senes Lab), RANDOM"
    # run fit sin to data
    tools.fit_sin_to_data(ftsb_x,ftsb_y_shuffled,P0A6S5_title_random,output_folder,output_base_filename_random,savefig=True)
    """
    if fit_unbiased == False and fit_alphahelix == False:
        raise ValueError("Either fit_unbiased or fit_alphahelix must be True")

    #global orig_data_df, df_sine_constants
    output_folder_pdf = os.path.join(output_folder, "pdf")
    output_basename = os.path.join(output_folder, output_base_filename)
    output_basename_pdf = os.path.join(output_folder_pdf, output_base_filename)

    # create output folders, if necessary
    list_paths = [output_folder, output_folder_pdf]
    for path in list_paths:
        if not os.path.exists(path):
            os.makedirs(path)

    # create a number of x data points for smooth sine curves
    x_smooth = np.linspace(x.min(), x.max(), 500)
    # create dataframe to save the sine constants
    df_sine_constants = pd.DataFrame(index = ["a","b","c","d","periodicity"])
    # create dataframe to save the fitted curves
    df_fitted_curve = pd.DataFrame(index=x_smooth)

    if fit_unbiased:
        #####################################################################################################
        #                                                                                                   #
        #              CALCULATE SINE FIT TO FREE HELIX (ALTHOUGH GUESS HAS PERIODICITY OF 3.6)             #
        #                                                                                                   #
        #####################################################################################################

        # guess the constants in the initial curve
        sine_constants_guess = [1.0,1.7,0.2,0]
        # fit a sine curve to the data using the leastsq method
        sine_constants, cov, infodict, mesg, ier = scipy.optimize.leastsq(mathfunctions.residuals,
                                                                          sine_constants_guess,
                                                                          args=(mathfunctions.sine,x,y),
                                                                          full_output=1)

        # print the periodicity of the curve that best fits the data (hopefully ~3.6, corresponding to an alpha-helix!)
        periodicity_fitted_curve = 2 * np.pi / sine_constants[1]
        print("periodicity of fitted curve = %0.2f" % periodicity_fitted_curve)
        # add the sine constants and periodicity to the output dataframe
        df_sine_constants["sine_constants_unbiased"] = sine_constants.tolist() + [periodicity_fitted_curve]
        # plot the fitted curve to the data
        yvalues_fitted = mathfunctions.sine(sine_constants, x_smooth)
        df_fitted_curve['yvalues_fitted_unbiased'] = yvalues_fitted

    if fit_alphahelix:
        #####################################################################################################
        #                                                                                                   #
        #                CALCULATE SINE FIT TO PERFECT ALPHA HELIX (FIXED PERIODICITY OF 3.6)               #
        #                                                                                                   #
        #####################################################################################################
        #RECALCULATE SINE CONSTANTS ASSUMING A PERFECT HELIX (a=0.2, b=1.7, periodicity = 3.6)
        sine_constants_guess_perfhelix = [1.0,0.5]
        # fit a sine curve to the data using the leastsq method
        sine_constants_A, cov_A, infodict_A, mesg_A, ier_A = scipy.optimize.leastsq(mathfunctions.residuals,
                                                                             sine_constants_guess_perfhelix,
                                                                             args=(mathfunctions.sine_perfect_helix,x,y),
                                                                             full_output=1)
        # add the sine constants to output dataframe
        df_sine_constants.loc["c":"d","sine_constants_alpha_helix"] = sine_constants_A
        # plot the fitted curve to the data
        yvalues_fitted_perfhelix = mathfunctions.sine_perfect_helix(sine_constants_A, x_smooth)
        df_fitted_curve['yvalues_fitted_fixed_periodicity'] = yvalues_fitted_perfhelix

    if x_grid:
        raise ValueError("ax not created?")
        #ax.grid(axis='x', which='major', alpha=0.5)

    if savefig:
        #####################################################################################################
        #                                                                                                   #
        #                                           SAVE OUTPUT FIGURES                                     #
        #                                                                                                   #
        #####################################################################################################
        # set some figure settings
        plt.rcParams["savefig.dpi"] = 120
        plt.rc('font', family='Arial')
        fontsize = 12
        # create custom list of colours for figures
        colour_lists = create_colour_lists()
        # create new fig
        fig, ax = plt.subplots()
        # plot the raw disruption data against the amino acid number
        ax.plot(x, y, color = colour_lists['TUM_colours']['TUMBlue'], label = "original data")
        if fit_unbiased:
            ax.plot(x_smooth, yvalues_fitted, linestyle = "--", color = colour_lists['TUM_colours']['TUM4'], label = "fit to sine (unbiased)")
            ax.annotate(s='periodicity = %0.1f' % periodicity_fitted_curve,
                        xy = (0.04,0.9),color = colour_lists['TUM_colours']['TUM4'],
                        fontsize=fontsize, xytext=None, xycoords='axes fraction',
                        alpha=0.75)

        if fit_alphahelix:
            ax.plot(x_smooth, yvalues_fitted_perfhelix, linestyle = "--", color = colour_lists['TUM_accents']['green'], label = "fit to sine (fixed periodicity)")
            ax.annotate(s='periodicity = 3.6',
                        xy = (0.04,0.8),color = colour_lists['TUM_accents']['green'],
                        fontsize=fontsize, xytext=None, xycoords='axes fraction',
                        alpha=0.75)

        raise ValueError("code needs to be rearranged if you want this to work. orig_aa not found.")
        set_ticks_and_labels_and_save(fig, ax, x, y, orig_aa, output_basename, output_basename_pdf, fontsize, title)

    if saveexcel:
        # save data to new excel file
        writer = pd.ExcelWriter(output_basename + ".xlsx")
        df_sine_constants.to_excel(writer, sheet_name = "sine_constants")
        df_fitted_curve.to_excel(writer, sheet_name='yvalues_fitted')
        writer.save()
        writer.close()
    if savecsv:
        df_sine_constants.to_csv(output_basename + "_sine_constants.csv")
        df_fitted_curve.to_csv(output_basename + "_yvalues_fitted.csv")

    # return sine_constants, periodicity_fitted_curve

def set_ticks_and_labels_and_save(fig, ax, x, y, orig_aa, output_basename, output_basename_pdf, fontsize, title):
    """ Set the ticks and other parameters in the output figure, and save as png and pdf

    Parameters
    ----------
    fig : figure object
    ax : ax (plot) object
    x : array-like data of the x axis (e.g. aa position)
    y : array-like data of the y axis (e.g. disruption)
    output_basename : base filepath of output figures
    output_basename_pdf : base filepath of output pdf figures
    fontsize : figure fontsize
    title : figure title

    Saved Files and Figures
    -------
    output_basename + ".png" : output figure, png
    output_basename_pdf + ".pdf" : output figure, pdf
    """

    """ Adjusts some of the figure parameters, and saves as a png and pdf using the basepath indicated
    """
    ax.set_xticks(x)
    ax.set_xticklabels(x, fontsize  = fontsize, rotation = 90)
    #use original aa as x axis label
    ax.set_xticklabels(orig_aa, fontsize  = fontsize, rotation = 0)
    ax.set_xlabel('orignal aa', fontsize=fontsize)
    ax.set_ylabel('average value', fontsize=fontsize)
    y_dist = max(y)-min(y)
    ylim_min = min(y) - y_dist*0.2 if min(y) - y_dist*0.2 > 0 else 0
    ylim_max = max(y) + y_dist*0.2
    #ax.set_ylim(min(y),max(y) + 0.4)
    ax.set_ylim(ylim_min, ylim_max + 0.4)
    ax.set_title(title, fontsize=fontsize)
    # show the legend(as previously labelled)
    ax.legend(fontsize=fontsize)
    #autoscale_x_axis
    plt.autoscale(enable=True, axis='x', tight=True)
    # save figure
    fig.savefig(output_basename + ".png")
    fig.savefig(output_basename_pdf + ".pdf", format="pdf")
