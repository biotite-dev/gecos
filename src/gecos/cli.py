import os
import argparse
import copy
import sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import biotite.sequence as seq
import biotite.sequence.align as align
from .space import ColorSpace
from .optimizer import ColorOptimizer
from .colors import convert_lab_to_rgb
from .file import write_color_scheme


def handle_error(func):
    def wrapped(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except InputError as e:
            print(e, file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print("An unexpected error occured:\n", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)
    return wrapped


class InputError(Exception):
    pass


@handle_error
def main(args=None):
    parser = argparse.ArgumentParser(
        description="This program automatically generates a color scheme for "
                    "sequence alignments. "
                    "The visual differences of the colors correspond to the "
                    "distances in a given substitution matrix."
    )
    
    space_group  = parser.add_argument_group(
        title="Color space",
        description=None
    )
    matrix_group = parser.add_argument_group(
        title="Substitution matrix",
        description=None
    )
    opt_group = parser.add_argument_group(
        title="Visual distance optimization",
        description=None
    )
    output_group = parser.add_argument_group(
        title="Output files",
        description=None
    )
    vis_group    = parser.add_argument_group(
        title="Visualization",
        description=None
    )

    space_group.add_argument(
        "--lightness", "-l", type=int, required=True,
        help="The lightness (brightness) of the color space (1 - 99). "
             "This argument must be provided."
    )
    space_group.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Show only the customized color space and terminate the program."
    )
    space_group.add_argument("--smin", type=int,
                             help="All colors in the space must "
                                  "have the specified saturation at minimum "
                                  "(a^2 + b^2 >= smin^2)."
    )
    space_group.add_argument("--smax", type=int,
                             help="All colors in the space must "
                                  "have the specified saturation at maximum "
                                  "(a^2 + b^2 <= smax^2)."
    )

    matrix_group.add_argument(
        "--alphabet", "-a",
        help="A custom alphabet to generate the scheme for (e.g. 'ACGT'). "
             "By default an alphabet containing the 20 amino acids is used. "
    )
    matrix_group.add_argument(
        "--matrix", "-m", default="BLOSUM62",
        help="The substitution matrix to calculate the pairwise symbol "
             "distances from. "
             "Can be an NCBI substitution matrix name (e.g. 'BLOSUM62') or "
             "alternatively a substitution matrix file in NCBI format "
             "(e.g. 'blosum62.mat'). "
             "Default: 'BLOSUM62'"
    )
    
    output_group.add_argument(
        "--scheme-file", "-s",
        type=argparse.FileType(mode="w"), default=sys.stdout,
        help="Write the generated color scheme into the specified file. "
             "The scheme is saved as Biotite-compatible JSON file. "
             "By default the scheme is output to STDOUT."
    )
    output_group.add_argument(
        "--name", default="scheme",
        help="The name of the color scheme that is used in the JSON file".
    )
    output_group.add_argument(
        "--pot-file", "-p", type=argparse.FileType(mode="w"),
        help="Write the potentials during the color scheme optimization "
             "into the specified file. "
             "Each line corresponds to one optimization step."
    )

    vis_group.add_argument(
        "--show-space", action="store_true",
        help="Show the generated color space."
    )
    vis_group.add_argument(
        "--show-scheme", action="store_true",
        help="Show the distribution of alphabet symbol in the color space."
    )
    vis_group.add_argument(
        "--show-example", action="store_true",
        help="Show an example multiple sequence alignment "
             "with newly generated color space."
    )
    vis_group.add_argument(
        "--show-pot", action="store_true",
        help="Show a plot of the potential during the optimization process."
    )

    args = parser.parse_args(args=args)


    alphabet = parse_alphabet(args.alphabet)
    matrix = parse_matrix(args.matrix)
    
    space = create_space(args.lightness)
    adjust_saturation(space, args.smin, args.smax)

    optimizer = ColorOptimizer(matrix, space)
    temps      = [100, 80, 60, 40, 20, 10, 8,   6,   4,   2,   1  ]
    step_sizes = [10,  8,  6,  4,  2,  1,  0.8, 0.6, 0.4, 0.2, 0.1]
    for temp, step_size in zip(temps, step_sizes): 
        with ProcessPoolExecutor() as executor:
            futures = [
                executor.submit(
                    optimize, copy.deepcopy(optimizer), 1000, temp, step_size
                ) for _ in range(10)
            ]
            splitted_optimizers = [future.result() for future in futures]
        pot = [opt.get_result().potential for opt in splitted_optimizers]
        optimizer = splitted_optimizers[np.argmin(pot)]
    result = optimizer.get_result()

    write_scheme(args.scheme_file, result, args.name)
    if args.pot_file:
        write_potential(args.pot_file, result)

    if args.show_space:
        show_space(space)
    if args.show_scheme:
        show_scheme(space, result)
    if args.show_example:
        show_example(result)
    if args.show_pot:
        show_potential(result)


def parse_alphabet(alphabet_str):
    if alphabet_str is None:
        return seq.LetterAlphabet(
            seq.ProteinSequence.alphabet.get_symbols()[:20]
        )
    else:
        if " " in alphabet_str:
            raise InputError("Alphabet may not contain whitespaces")
        try:
            return seq.LetterAlphabet(input)
        except Exception:
            raise InputError("Invalid alphabet")

def parse_matrix(matrix_str):
    if os.path.isfile(matrix_str):
        with open(matrix_str) as f:
            matrix_dict = align.SubstitutionMatrix.dict_from_str(f.read())
            return align.SubstitutionMatrix(alphabet, alphabet, matrix_dict)
    else:
        # String is a NCBI matrix name
        upper_matrix_str = matrix_str.upper()
        if upper_matrix_str not in align.SubstitutionMatrix.list_db():
            raise InputError(
                f"'{matrix_str}' is neither a file "
                f"nor a valid NCBI substitution matrix"
            )
        return align.SubstitutionMatrix(alphabet, alphabet, upper_matrix_str)


def create_space(lightness):
    if lightness < 1 or lightness > 99:
        raise InputError("Lightness value must be between 1 and 99") 
    return ColorSpace(lightness)

def adjust_saturation(space, smin, smax):
    if smin is not None:
        space.remove((a**2 + b**2 < smin**2))
    if smax is not None:
        space.remove((a**2 + b**2 > smax**2))


def write_scheme(file, result, name):
    write_color_scheme(file, result, name)

def write_potential(file, result):
    np.savetxt(file, result.potentials)


def show_space(space):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    rgb_matrix = space.get_rgb_matrix()
    rgb_matrix[np.isnan(rgb_matrix)] = 0.7
    ax.imshow(np.transpose(rgb_matrix, axes=(1,0,2)), origin="lower",
                extent=(-128, 127,-128, 127), aspect="equal")
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    plt.show(block=False)

def show_scheme(space, result):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.matshow(space.space.T, extent=(-128, 127,-128, 127),
               origin="lower", cmap=ListedColormap([(0.7,0.7,0.7), (1,1,1)]))
    for symbol, pos, color in zip(result.alphabet, result.coord, result.rgb_colors):
        ax.text(pos[1], pos[2], symbol, color=color,
                ha="center", va="center", size=14, weight="heavy")
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    plt.show(block=False)

def show_example(result):
    pass

def show_potential(result):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.plot(result.potentials)
    ax.set_xlabel("Step")
    ax.set_ylabel("Potential")
    plt.show(block=False)


def optimize(optimizer, n_steps, temp, step_size):
    optimizer.optimize(n_steps, temp, step_size)
    return optimizer