from multiprocessing import Pool
from os.path import join, dirname, realpath, isfile
import itertools
import argparse
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta
import biotite.sequence.graphics as graphics
from .space import ColorSpace
from .optimizer import ColorOptimizer, DefaultScoreFunction
from .file import write_color_scheme


FIGURE_WIDTH = 8.0
EXAMPLE_FILE_NAME \
    = join(dirname(realpath(__file__)), "example_alignment.fasta")


def handle_error(func):
    def wrapped(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except InputError as e:
            print(e, file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print("An unexpected error occured:\n", file=sys.stderr)
            raise
    return wrapped


class InputError(Exception):
    pass


def _optimize(seed, alphabet, score_func, space, constraints, 
              nsteps, beta_start, beta_end, stepsize_start, stepsize_end):
    """
    Worker function used for parallel execution of simulated annealing.
    """
    np.random.seed(seed)
    optimizer = ColorOptimizer(
        alphabet, score_func, space, constraints
    )
    optimizer.optimize(
        nsteps, beta_start, beta_end, stepsize_start, stepsize_end
    )
    return optimizer.get_result()

@handle_error
def main(args=None, result_container=None, show_plots=True):
    parser = argparse.ArgumentParser(
        description="This program automatically generates a color scheme for "
                    "sequence alignments. "
                    "The perceptual color differences correspond to the "
                    "distances in a given substitution matrix."
                    "\n"
                    "The algorithm tries to find an optimal color scheme by "
                    "means of color differences by performing a "
                    "Metropolis-Monte-Carlo optimization in color space.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    space_group  = parser.add_argument_group(
        title="Color space",
        description=None
    )
    matrix_group = parser.add_argument_group(
        title="Substitution matrix and alphabet",
        description=None
    )
    opt_group = parser.add_argument_group(
        title="Optimization",
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
        "--dry-run", "-n", action="store_true",
        help="Show only the customized color space and terminate the program."
    )
    space_group.add_argument(
        "--smin", type=int,
        help="All colors in the space must "
            "have the specified saturation at minimum "
            "(a^2 + b^2 >= smin^2).",
        metavar="S"
    )
    space_group.add_argument(
        "--smax", type=int,
        help="All colors in the space must "
             "have the specified saturation at maximum "
             "(a^2 + b^2 <= smax^2).",
        metavar="S"
    )
    space_group.add_argument(
        "--lmin", type=int,
        help="All colors in the space must "
             "have the specified lightness at minimum.",
        metavar="L"
    )
    space_group.add_argument(
        "--lmax", type=int,
        help="All colors in the space must "
             "have the specified lightness at maximum.",
        metavar="L"
    )
    space_group.add_argument(
        "--amin", type=int,
        help="All colors in the space must "
             "have the specified 'a*' value at minimum.",
        metavar="A"
    )
    space_group.add_argument(
        "--amax", type=int,
        help="All colors in the space must "
             "have the specified 'a*' value at maximum.",
        metavar="A"
    )
    space_group.add_argument(
        "--bmin", type=int,
        help="All colors in the space must "
             "have the specified 'b*' value at minimum.",
        metavar="B"
    )
    space_group.add_argument(
        "--bmax", type=int,
        help="All colors in the space must "
            "have the specified 'b*' value at maximum.",
        metavar="B"
    )

    matrix_group.add_argument(
        "--alphabet", "-a",
        help="A custom alphabet to generate the scheme for (e.g. 'ACGT'). "
             "By default an alphabet containing the 20 amino acids is used. ",
        metavar="SYMBOLS"
    )
    matrix_group.add_argument(
        "--matrix", "-m", default="BLOSUM62",
        help="The substitution matrix to calculate the pairwise symbol "
             "distances from. "
             "Can be an NCBI substitution matrix name (e.g. 'BLOSUM62') or "
             "alternatively a substitution matrix file in NCBI format "
             "(e.g. 'blosum62.mat').",
        metavar="MATRIX"
    )

    opt_group.add_argument(
        "--contrast", default=700, type=int,
        help="The contrast factor controls how strongly the symbols are "
             "pushed to the edges of the color space. "
             "At the minimum value '0' contrast is not rewarded.",
        metavar="FACTOR"
    )
    opt_group.add_argument(
        "--constraint", "-c", nargs=4, action="append",
        help="Constrain a symbol to a fixed position in color space. "
             "This argument takes four values: a symbol, "
             "and the color channels 'l*', 'a*' and 'b*'"
             "(e.g. 'A 65 10 15')."
             "Can be repeated to add constraints for multiple symbols.",
        metavar=("SYMBOL", "L", "A", "B")
    )
    opt_group.add_argument(
        "--nsteps", default=20000, type=int,
        help="The optimization process performs simulated annealing in order "
             "to find the optimal conformation of the color scheme. "
             "This parameter sets the total amount of optimization steps. "
             "With a higher number of steps the quality of the optimization "
             "increases at the cost of a longer runtime.",
        metavar="NUMBER"
    )
    opt_group.add_argument(
        "--delta", default="CIEDE2000",
        choices=["CIE76", "CIEDE94", "CIEDE2000"],
        help="The formula to use for the calculation of perceptual "
             "differences.",
        metavar="FORMULA"
    )
    opt_group.add_argument(
        "--seed", type=float,
        help="Start seed used for seeding the parallel runs. "
             "By default the seed is chosen randomly.",
        metavar="NUMBER"
    )
    opt_group.add_argument(
        "--nruns", default=16, type=int,
        help="Number of optimizations to run. "
             "From these runs, the color scheme with the best score is "
             "selected.",
        metavar="NUMBER"
    )
    opt_group.add_argument(
        "--nthreads", type=int,
        help="Number of optimization runs to run in parallel. "
             "By default this is equal to the number of CPUs.",
        metavar="NUMBER"
    )
    opt_group.add_argument(
        "--beta-start", default=1, type=float,
        help="Inverse temperature at start of simulated annealing.",
        metavar="FLOAT",
    )
    opt_group.add_argument(
        "--beta-end", default=500, type=float,
        help="Inverse temperature at end of simulated annealing.",
        metavar="FLOAT",
    )   
    opt_group.add_argument(
        "--stepsize-start", default=10, type=float,
        help="Step size temperature at start of simulated annealing.",
        metavar="FLOAT",
    )
    opt_group.add_argument(
        "--stepsize-end", default=0.2, type=float,
        help="Step size temperature at end of simulated annealing.",
        metavar="FLOAT",
    )     
    
    output_group.add_argument(
        "--scheme-file", "-f",
        type=argparse.FileType(mode="w"), default=sys.stdout,
        help="Write the generated color scheme into the specified file. "
             "The scheme is saved as Biotite-compatible JSON file. "
             "By default the scheme is output to STDOUT.",
        metavar="FILE"
    )
    output_group.add_argument(
        "--name", default="scheme",
        help="The name of the color scheme that is used in the JSON file.",
        metavar="NAME"
    )
    output_group.add_argument(
        "--score-file", type=argparse.FileType(mode="w"),
        help="Write the scores during the color scheme optimization "
             "into the specified file. "
             "Each line corresponds to one optimization step.",
        metavar="FILE"
    )

    vis_group.add_argument(
        "--lightness", type=int,
        help="Set the lightness level for color space visualization. "
             "The plots generated via the '--show-space' and the "
             "'--show-scheme' options show the color space excerpt at the "
             "given lightness level. "
             "By default the mean of the minimum and maximum lightness of the "
             "space is chosen.",
        metavar="LIGHTNESS"
    )
    vis_group.add_argument(
        "--show-space", action="store_true",
        help="Show the generated color space. "
             "The lightness of the '--lightness' option is used here. "
             "The gray area cannot be displayed in RGB space at the given "
             "lightness."
    )
    vis_group.add_argument(
        "--show-scheme", action="store_true",
        help="Show the distribution of alphabet symbols in the color space."
             "The lightness of the '--lightness' option is used here. "
             "The gray area cannot be displayed in RGB space at the given "
             "lightness."
    )
    vis_group.add_argument(
        "--show-example", action="store_true",
        help="Show an example multiple sequence alignment "
             "with newly generated color space. "
             "Cannot be used in combination with a custom '--alphabet' value."
    )
    vis_group.add_argument(
        "--show-score", action="store_true",
        help="Show a plot of the score during the optimization process."
    )

    args = parser.parse_args(args=args)


    alphabet = parse_alphabet(args.alphabet)
    matrix = parse_matrix(args.matrix, alphabet)
    
    # This lightness is only used fotr visualization purposes
    if args.lightness is not None:
        lightness = args.lightness
    elif args.lmin is not None and args.lmax is not None:
        lightness = (args.lmin + args.lmax) // 2
    else:
        lightness = 50

    space = ColorSpace()
    adjust_saturation(space, args.smin, args.smax)
    adjust_l(space, args.lmin, args.lmax)
    adjust_a(space, args.amin, args.amax)
    adjust_b(space, args.bmin, args.bmax)
    if args.dry_run:
        figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
        ax = figure.gca()
        show_space(ax, space, lightness)
        figure.tight_layout()
        if show_plots:
            plt.show()
        return

    constraints = np.full((len(alphabet), 3), np.nan)
    if args.constraint is not None:
        for symbol, l, a, b in args.constraint:
            constraints[alphabet.encode(symbol)] = (l,a,b)
    
    score_func = DefaultScoreFunction(matrix, args.contrast, args.delta)   
    
    # Simulated annealing     
    if args.seed is not None:
        np.random.seed(int(args.seed))
    # Different random seed for each run
    seeds = np.random.randint(0, 1000000, size=args.nruns)

    with Pool(args.nthreads) as p:
        results = p.starmap(_optimize, zip(
            seeds,
            itertools.repeat(matrix.get_alphabet1()),
            itertools.repeat(score_func),
            itertools.repeat(space),
            itertools.repeat(constraints),
            itertools.repeat(args.nsteps), 
            itertools.repeat(args.beta_start),
            itertools.repeat(args.beta_end),
            itertools.repeat(args.stepsize_start),
            itertools.repeat(args.stepsize_end),
        ))
    best_result = min(results, key=lambda x: x.score)

    scores = np.array([result.scores for result in results])
    scores_mean = np.mean(scores, axis=0)
    scores_std = np.std(scores, axis=0)
    scores_min = np.min(scores, axis=0)
    scores_max = np.max(scores, axis=0)                
    score_results = np.array(
        [scores_mean, scores_std, scores_min, scores_max]
    )
    
    if args.score_file:
        write_score(
            args.score_file, score_results,
            header="avg(score) std(score) min(score) max(score)"
        )
    write_scheme(args.scheme_file, best_result, args.name)


    if args.show_space:
        figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
        ax = figure.gca()
        show_space(ax, space, lightness)
        figure.tight_layout()
    if args.show_scheme:
        figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
        ax = figure.gca()
        show_scheme(ax, space, best_result, lightness)
        figure.tight_layout()
    if args.show_example:
        # Check whether a custom non-amino-acid alphabet is used
        if args.alphabet is not None:
            raise InputError(
                "The example alignment can only be shown "
                "for the amino acid alphabet"
            )
        figure = plt.figure(figsize=(FIGURE_WIDTH, 2.5))
        ax = figure.gca()
        show_example(ax, best_result.rgb_colors)
        figure.tight_layout()
    if args.show_score:
        figure = plt.figure(figsize=(FIGURE_WIDTH, 6.0))
        ax = figure.gca()
        show_score(ax, best_result.scores)
        figure.tight_layout()
    if show_plots:
        plt.show()

    # In case someone wants to use the CLI results in a Python script
    # In any other case: Ignore the following two lines
    if result_container is not None:
        result_container.append(best_result)


def parse_alphabet(alphabet_str):
    if alphabet_str is None:
        return seq.LetterAlphabet(
            seq.ProteinSequence.alphabet.get_symbols()[:20]
        )
    else:
        if " " in alphabet_str:
            raise InputError("Alphabet may not contain whitespaces")
        try:
            return seq.LetterAlphabet(alphabet_str)
        except Exception:
            raise InputError("Invalid alphabet")


def parse_matrix(matrix_str, alphabet):
    if isfile(matrix_str):
        with open(matrix_str) as f:
            matrix_dict = align.SubstitutionMatrix.dict_from_str(f.read())
            return align.SubstitutionMatrix(alphabet, alphabet, matrix_dict)
    else:
        # String is a NCBI matrix name
        # For user convenience there is no case sensitivity
        # -> Find fitting matrix
        matrix_list = align.SubstitutionMatrix.list_db()
        upper_matrix_str = matrix_str.upper()
        upper_matrix_list = [
            m.upper() for m in align.SubstitutionMatrix.list_db()
        ]
        try:
            matrix_str = matrix_list[upper_matrix_list.index(upper_matrix_str)]
        except:
            raise InputError(
                f"'{matrix_str}' is neither a file "
                f"nor a valid NCBI substitution matrix"
            )
        return align.SubstitutionMatrix(alphabet, alphabet, matrix_str)


def adjust_saturation(space, smin, smax):
    lab = space.lab
    a = lab[..., 1]
    b = lab[... ,2]
    if smin is not None:
        space.remove(a**2 + b**2 < smin**2)
    if smax is not None:
        space.remove(a**2 + b**2 > smax**2)

def adjust_l(space, lmin, lmax):
    lab = space.lab
    l = lab[..., 0]
    if lmin is not None:
        space.remove(l < lmin)
    if lmax is not None:
        space.remove(l > lmax)

def adjust_a(space, amin, amax):
    lab = space.lab
    a = lab[..., 1]
    if amin is not None:
        space.remove(a < amin)
    if amax is not None:
        space.remove(a > amax)

def adjust_b(space, bmin, bmax):
    lab = space.lab
    b = lab[..., 2]
    if bmin is not None:
        space.remove(b < bmin)
    if bmax is not None:
        space.remove(b > bmax)


def write_scheme(file, result, name):
    write_color_scheme(file, result, name)

def write_score(file, results, header=""):    
    results = np.transpose(results)
    np.savetxt(file, results, fmt="%.2f", header=header)


def show_space(ax, space, lightness):
    # For performance reasons, filter correct lightness
    # before converting to RGB 
    l = space.lab[..., 0]
    space = copy.deepcopy(space)
    space.remove(l != lightness)
    # Remove first dimension (length = 1)
    rgb_space = space.get_rgb_space()[lightness]
    rgb_space[np.isnan(rgb_space)] = 0.7
    ax.imshow(np.transpose(rgb_space, axes=(1,0,2)), origin="lower",
                extent=(-128, 127,-128, 127), aspect="equal")
    ax.set_xlabel("a*")
    ax.set_ylabel("b*")
    ax.set_title(f"L* = {lightness}")

def show_scheme(ax, space, result, lightness):
    s = space.space
    ax.matshow(space.space[lightness].T, extent=(-128, 127,-128, 127),
               origin="lower", cmap=ListedColormap([(0.7,0.7,0.7), (1,1,1)]))
    for symbol, pos, color in zip(
        result.alphabet, result.lab_colors, result.rgb_colors
    ):
        ax.text(pos[1], pos[2], symbol, color=color,
                ha="center", va="center", size=14, weight="heavy")
    ax.xaxis.set_ticks_position("bottom")
    ax.set_xlabel("a*")
    ax.set_ylabel("b*")
    ax.set_title(f"L* = {lightness}")

def show_example(ax, colors):
    fasta_file = fasta.FastaFile.read(EXAMPLE_FILE_NAME)
    alignment = fasta.get_alignment(fasta_file)
    alignment = alignment[:60]

    graphics.plot_alignment_type_based(
        ax, alignment, spacing=2.0, symbols_per_line=len(alignment),
        color_scheme=colors
    )

def show_score(ax, scores):
    ax.plot(scores)
    ax.set_xlabel("Step")
    ax.set_ylabel("Score")