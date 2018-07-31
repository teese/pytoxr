#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License.
"""
import matplotlib.colors as colors
import numpy as np
# figures plotted for each experiment can be put here.

def create_colour_lists():
    '''
    Converts several lists of rgb colours to the python format (normalized to between 0 and 1)
    Returns a dictionary that contains dictionaries of palettes with named colours (eg. TUM blues)
    and also lists of unnamed colours (e.g. tableau20)
    (copied from korbinian package of Mark Teese (MIT licence))
    '''
    output_dict = {}

    matplotlib_150 = list(colors.cnames.values())
    output_dict['matplotlib_150'] = matplotlib_150

    #define colour dictionaries. TUM colours are based on the style guide.
    colour_dicts = {
                    'TUM_colours' : {
                                    'TUMBlue':(34,99,169),
                                    'TUM1':(100,160,200),
                                    'TUM2':(1,51,89),
                                    'TUM3':(42,110,177),
                                    'TUM4':(153,198,231),
                                    'TUM5':(0,82,147)
                                    },
                    'TUM_oranges': {
                        'TUM0': (202, 101, 10),
                        'TUM1': (213, 148, 96),
                        'TUM2': (102, 49, 5),
                        'TUM3': (220, 108, 11),
                        'TUM4': (247, 194, 148),
                        'TUM5': (160, 78, 8)
                    },
                    'TUM_accents' : {
                                    'green':(162,183,0),
                                    'orange':(227,114,34),
                                    'ivory':(218,215,203),
                                    }
                    }

    #convert the nested dicts to python 0 to 1 format
    for c_dict in colour_dicts:
        for c in colour_dicts[c_dict]:
            #define r, g, b as ints
            r, g, b = colour_dicts[c_dict][c]
            #normalise r, g, b and add to dict
            colour_dicts[c_dict][c] = (r / 255., g / 255., b / 255.)
        #add normalised colours to output dictionary
        output_dict[c_dict] = colour_dicts[c_dict]

    #define colour lists
    colour_lists = {
                    'tableau20' : [
                                 (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
                                    ],
                    'tableau20blind' : [
                                         (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
                                         (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
                                         (255, 188, 121), (207, 207, 207)
                                          ]
                    }

    # custom colour lists
    output_dict['HTML_list01'] = ['#808080', '#D59460', '#005293', '#A1B11A', '#9ECEEC', '#0076B8', '#454545', "#7b3294", "#c2a5cf", "#008837", "#a6dba0"]
    output_dict['BGO_rgb'] = [(0,93,151), (83,167,214),(158,206,236),(94,152,49),(105,190,158),(166,213,196),(238,119,4),(214,168,82)] + colour_lists["tableau20"]
    output_dict['BGO_arr'] = [np.array(x)/255 for x in output_dict['BGO_rgb']]
    output_dict['BGO_dark_arr'] = [darken_or_lighten(c, -0.33) for c in output_dict['BGO_arr']]

    #normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list])/255.
        colour_array_tup = tuple(map(tuple,colour_array))
        colour_lists[rgb_list] = colour_array_tup
        #add normalised colours to output dictionary
        output_dict[rgb_list] = colour_lists[rgb_list]
    #create a mixed blue/grey colour list, with greys in decreasing darkness
    TUM_colours_list_with_greys = []
    grey = 0.7
    for c in colour_dicts['TUM_colours'].values():
        TUM_colours_list_with_greys.append('%0.2f' % grey)
        TUM_colours_list_with_greys.append(c)
        grey -= 0.1
    output_dict['TUM_colours_list_with_greys'] = TUM_colours_list_with_greys

    return output_dict

def darken_or_lighten(col_array, fraction):
    col_array = col_array + fraction
    col_array[col_array > 1] = 1
    col_array[col_array < 0] = 0
    return col_array