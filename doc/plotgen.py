# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from os.path import join, isdir, isfile, dirname, realpath
from os import mkdir
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random
import biotite as biotite
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.graphics as graphics
import gecos
import gecos.cli as gecli


DOC_PATH = dirname(realpath(__file__))


plot_generators = {}

def generate(plot_directory_path):
    if not isdir(plot_directory_path):
        mkdir(plot_directory_path)
    for plot_name, function in plot_generators.items():
        plot_path_name = join(plot_directory_path, plot_name) + ".png"
        # Generate plot only when not already existing
        if not isfile(plot_path_name):
            figure = function()
            plt.savefig(plot_path_name)

def plot_generator(function):
    function_name = function.__name__
    if not function_name.startswith("plot_"):
        raise ValueError("Plotter functions must start with 'plot_'")
    plot_name = function_name[5:]
    plot_generators[plot_name] = function
    return function



@plot_generator
def plot_main_example_alignment():
    random.seed(4)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--matrix", "BLOSUM62",
        "--lmin", "60",
        "--lmax", "75",
        "-s", scheme_file
    ])
    show_alignment(scheme_file)

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
    random.seed(0)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=["-s", scheme_file])
    show_alignment(scheme_file)

@plot_generator
def plot_no_green_scheme_alignment():
    random.seed(0)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--amin", "0",
        "--lmin", "50",
        "--lmax", "80",
        "-s", scheme_file
    ])
    show_alignment(scheme_file)

@plot_generator
def plot_high_saturation_scheme_alignment():
    random.seed(1)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--smin", "30",
        "--lmin", "55",
        "--lmax", "75",
        "-s", scheme_file
    ])
    show_alignment(scheme_file)

@plot_generator
def plot_constrained_scheme_alignment():
    random.seed(0)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "-c", "A", "70", "0", "0",
        "-c", "W", "70", "-20", "-40",
        "--lmin", "60",
        "--lmax", "75",
        "-s", scheme_file
    ])
    show_alignment(scheme_file)

@plot_generator
def plot_high_contrast_scheme_alignment():
    random.seed(0)
    scheme_file = biotite.temp_file("json")
    gecli.main(args=[
        "--contrast", "100",
        "--lmin", "60",
        "--lmax", "75",
        "-s", scheme_file
    ])
    show_alignment(scheme_file)



def show_alignment(scheme_file):
    colors = graphics.load_color_scheme(scheme_file)["colors"]
    fasta_file = fasta.FastaFile()
    fasta_file.read(join(DOC_PATH, "example_alignment.fasta"))
    alignment = fasta.get_alignment(fasta_file)
    alignment = alignment[:60]

    fig = plt.figure(figsize=(8.0, 2.5))
    ax = fig.add_subplot(111)
    graphics.plot_alignment_type_based(
        ax, alignment, symbols_per_line=len(alignment), color_scheme=colors
    )
    fig.tight_layout()
    return fig

def show_space(ax, lightness):
    space = gecos.ColorSpace()
    # For performance reasons, filter correct lightness
    # before converting to RGB 
    l = space.lab[..., 0]
    space.remove(l != lightness)
    # Remove first dimension (length = 1)
    rgb_space = space.get_rgb_space()[lightness]
    rgb_space[np.isnan(rgb_space)] = 0.7
    ax.imshow(np.transpose(rgb_space, axes=(1,0,2)), origin="lower",
                extent=(-128, 127,-128, 127), aspect="equal")
    ax.set_xlabel("a*")
    ax.set_ylabel("b*")
    ax.set_title(f"L* = {lightness}")