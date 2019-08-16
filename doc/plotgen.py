# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from os.path import join, isdir, isfile, dirname, realpath
from os import mkdir
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random
from sphinx.util.logging import getLogger
from sphinx.util import status_iterator
import biotite as biotite
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta
import biotite.sequence.graphics as graphics
import gecos
import gecos.cli as gecli


DOC_PATH = dirname(realpath(__file__))
PB_EXAMPLE_FILE_NAME \
    = join(dirname(realpath(__file__)), "pb_alignment.fasta")


plot_generators = {}

def empty_func(*args, **kwargs):
    pass 

def generate(plot_directory_path):
    # Overwrite 'plt.show()' to prevent matplotlib blocking the execution
    show_func = plt.show()
    plt.show = empty_func
    # Create plot directory if not already existing
    if not isdir(plot_directory_path):
        mkdir(plot_directory_path)
    # Log progress in terminal
    logger = getLogger('sphinx-gallery')
    logger.info("generating plots...", color="white")
    iterator = status_iterator(
        plot_generators.items(),
        "generating plots...",
        length=len(plot_generators),
        stringify_func=lambda val : val[0]
    )
    # Run '@plot_generator' functions
    for plot_name, function in iterator:
        plot_path_name = join(plot_directory_path, plot_name) + ".png"
        # Generate plot only when not already existing
        if not isfile(plot_path_name):
            figure = function()
            figure.savefig(plot_path_name)
    # Replace 'plt.show()' with original function
    plt.show = show_func

def plot_generator(function):
    function_name = function.__name__
    if not function_name.startswith("plot_"):
        raise ValueError("Plotter functions must start with 'plot_'")
    plot_name = function_name[5:]
    plot_generators[plot_name] = function
    return function


@plot_generator
def plot_main_example_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--matrix", "BLOSUM62",
        "--lmin", "60",
        "--lmax", "75",
        "-f", scheme_file
    ])
    return show_alignment(scheme_file)

@plot_generator
def plot_example_space():
    figure = plt.figure(figsize=(8.0, 4.0))
    ax = figure.add_subplot(121)
    show_space(ax, 40)
    ax = figure.add_subplot(122)
    show_space(ax, 70)
    figure.tight_layout()
    return figure

@plot_generator
def plot_no_constraints_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=["--seed", "0", "-f", scheme_file])
    return show_alignment(scheme_file)

@plot_generator
def plot_no_green_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "1",
        "--amin", "0",
        "--lmin", "50",
        "--lmax", "80",
        "-f", scheme_file
    ])
    return show_alignment(scheme_file)

@plot_generator
def plot_high_saturation_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "1",
        "--smin", "30",
        "--lmin", "55",
        "--lmax", "75",
        "-f", scheme_file
    ])
    return show_alignment(scheme_file)

@plot_generator
def plot_constrained_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "-c", "A", "70", "0", "0",
        "-c", "W", "70", "-10", "-45",
        "--lmin", "60",
        "--lmax", "75",
        "-f", scheme_file
    ])
    return show_alignment(scheme_file)

@plot_generator
def plot_high_contrast_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--contrast", "2000",
        "--lmin", "60",
        "--lmax", "75",
        "-f", scheme_file
    ])
    return show_alignment(scheme_file)

@plot_generator
def plot_show_space():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--show-space",
        "--dry-run",
        "--smin", "30",
        "--lmin", "60",
        "--lmax", "70",
        "-f", scheme_file
    ])
    return plt.gcf()

@plot_generator
def plot_show_scheme():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--show-scheme",
        "--smin", "30",
        "--lmin", "60",
        "--lmax", "70",
        "-f", scheme_file
    ])
    return plt.gcf()

@plot_generator
def plot_show_example():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--show-example",
        "--smin", "30",
        "--lmin", "60",
        "--lmax", "70",
        "-f", scheme_file
    ])
    return plt.gcf()

@plot_generator
def plot_show_score():
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--seed", "0",
        "--show-score",
        "--smin", "30",
        "--lmin", "60",
        "--lmax", "70",
        "-f", scheme_file
    ])
    return plt.gcf()

###############################################################################
# The protein block substitution matrix, the protein block reference angles
# and the RSMDA calculation were taken from the PBxplore software,
# governed by the following license:
# 
# 
# The MIT License (MIT)
# 
# Copyright (c) 2013 Poulain, A. G. de Brevern
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

@plot_generator
def plot_pb_scheme_alignment():
    scheme_file = biotite.temp_file("json")
    mat_file = biotite.temp_file("mat")
    with open(mat_file, "w") as file:
        # PB substitution matrix, adapted from PBxplore
        file.write(
            """
                a     b     c     d     e     f     g     h     i     j     k     l     m     n     o     p
            a  516   -59   113  -105  -411  -177   -27  -361    47  -103  -644  -259  -599  -372  -124   -83
            b  -59   541  -146  -210  -155  -310   -97    90   182  -128   -30    29  -745  -242  -165    22
            c  113  -146   360   -14  -333  -240    49  -438  -269  -282  -688  -682  -608  -455  -147     6
            d -105  -210   -14   221     5  -131  -349  -278  -253  -173  -585  -670 -1573 -1048  -691  -497
            e -411  -155  -333     5   520   185   186   138  -378   -70  -112  -514 -1136  -469  -617  -632
            f -177  -310  -240  -131   185   459   -99   -45  -445    83  -214   -88  -547  -629  -406  -552
            g  -27   -97    49  -349   186   -99   665   -99   -89  -118  -409  -138  -124   172   128   254
            h -361    90  -438  -278   138   -45   -99   632  -205   316   192  -108  -712  -359    95  -399
            i   47   182  -269  -253  -378  -445   -89  -205   696   186     8    15  -709  -269  -169   226
            j -103  -128  -282  -173   -70    83  -118   316   186   768   196     5  -398  -340  -117  -104
            k -644   -30  -688  -585  -112  -214  -409   192     8   196   568   -65  -270  -231  -471  -382
            l -259    29  -682  -670  -514   -88  -138  -108    15     5   -65   533  -131     8   -11  -316
            m -599  -745  -608 -1573 -1136  -547  -124  -712  -709  -398  -270  -131   241    -4  -190  -155
            n -372  -242  -455 -1048  -469  -629   172  -359  -269  -340  -231     8    -4   703    88   146
            o -124  -165  -147  -691  -617  -406   128    95  -169  -117  -471   -11  -190    88   716    58
            p  -83    22     6  -497  -632  -552   254  -399   226  -104  -382  -316  -155   146    58   609
            """
        )
    gecli.main(args=[
        "--seed", "0",
        "--alphabet", "abcdefghijklmnop",
        "--matrix", mat_file,
        "--contrast", "300",
        "--lmin", "65",
        "--lmax", "70",
        "-f", scheme_file
    ])

    colors = graphics.load_color_scheme(scheme_file)["colors"]
    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.gca()

    pb_alphabet = seq.LetterAlphabet("abcdefghijklmnop")
    fasta_file = fasta.FastaFile()
    fasta_file.read(PB_EXAMPLE_FILE_NAME)
    seq_strings = list(fasta_file.values())
    sequences = [seq.GeneralSequence(pb_alphabet, seq_str.replace("-",""))
                 for seq_str in seq_strings]
    trace = align.Alignment.trace_from_strings(seq_strings)
    alignment = align.Alignment(sequences, trace, score=None)
    
    graphics.plot_alignment_type_based(
        ax, alignment, symbols_per_line=60, spacing=2, color_scheme=colors
    )

    fig.tight_layout()
    return fig


def show_alignment(scheme_file):
    colors = graphics.load_color_scheme(scheme_file)["colors"]
    fig = plt.figure(figsize=(8.0, 2.5))
    ax = fig.gca()
    gecli.show_example(ax, colors)
    fig.tight_layout()
    return fig

def show_space(ax, lightness):
    space = gecos.ColorSpace()
    gecli.show_space(ax, space, lightness)