import os
import copy
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


alphabet = None
matrix = None
space = None
optimizer = None
result = None
name = None

def main():
    try:
        dialog()
    except Exception:
        print()
        print(
            "Oops, something unexpected happened :( . "
            "Here is the error traceback:"
        )
        print()
        raise

def dialog():
    global alphabet
    global matrix
    global space
    global optimizer
    global result
    global name
    
    print(
        "Hello, this is Gecos, your friendly color scheme generator "
        "for sequence alignments. "
        "I will generate a color scheme for you, so that the visual "
        "differences of the colors correspond to the distances in a "
        "given substitution matrix. "
    )
    print()
    
    print(
        "Do you would like to create a color scheme for protein alignments?"
    )
    alphabet = process_input(choose_alphabet, {"y":"yes", "n":"no"})
    if alphabet is None:
        print(
            "Then please provide a custom alphabet to generate the scheme for "
            "(e.g. 'ACGT'):"
        )
        alphabet = process_input(parse_alphabet)
    
    print(
        "And now I need either the name of an NCBI substitution matrix "
        "(e.g. 'BLOSUM62') or a custom file name (e.g. 'blosum62.mat'):"
    )
    matrix = process_input(select_matrix)
    
    accepted_lightness = False
    accepted_space = False
    while not accepted_space:
        if not accepted_lightness:
            print(
                "Which lightness should the colors in the scheme have "
                "(1 - 99)?"
            )
            space = process_input(create_space)
            accepted_lightness = True
        print(
            "This is the color space I generated. "
            "Do you like to change something?"
        )
        show_space(space)
        mode = process_input(options={
            "1":"accept", "2":"rebuild", "3":"limit saturation"
        })
        if mode == "accept":
            accepted_space = True
        elif mode == "rebuild":
            accepted_lightness = False
        elif mode == "limit saturation":
            print(
                "Please give the minimum and maximum saturation "
                "(distance from center) of the color space (e.g. '25 - 50'):"
            )
            process_input(adjust_saturation)

    print(
        "Now I will arrange the alphabet symbols for you. "
        "This may take a while."
    )
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
    print(
        f"This is the scheme I found. "
        f"It has a potential of {result.potential:.2f}."
    )
    show_results(space, result)
    show_potential(result)
    
    print("How should the scheme be named (e.g. 'awesome_scheme')?")
    name = process_input()
    print(
        "Where do you like to save the color scheme "
        "(e.g. 'scheme.json')?"
    )
    process_input(save_scheme)
    print(
        "I saved your color scheme. Good bye!"
    )


def process_input(parse_function=None, options=None):
    while True:
        if options is not None:
            options_str = "["
            for option, value in options.items():
                options_str += option + ": " + value + ", "
            # Remove terminal comma
            options_str = options_str[:-2]
            options_str += "]"
            print(options_str)
        
        print("> ", end="")
        inp = input().strip().expandtabs(4)
        if len(inp) == 0:
            continue

        if options is not None:
            try:
                inp = options[inp.lower()]
            except KeyError:
                print("Please choose a proper option:")
                continue

        if parse_function is not None:
            try:
                return parse_function(inp)
            except InputError as e:
                print(str(e) + ".")
                print("Try again:")
        else:
            return inp


def choose_alphabet(input):
    if input == "yes":
        return seq.LetterAlphabet(
            seq.ProteinSequence.alphabet.get_symbols()[:20]
        )
    else:
        return None

def parse_alphabet(input):
    if " " in input:
        raise InputError("Alphabet may not contain whitespaces")
    try:
        return seq.LetterAlphabet(input)
    except Exception:
        raise InputError("Invalid alphabet")

def select_matrix(input):
    if os.path.isfile(input):
        with open(input) as f:
            matrix_dict = align.SubstitutionMatrix.dict_from_str(f.read())
            return align.SubstitutionMatrix(alphabet, alphabet, matrix_dict)
    else:
        # Input is a NCBI matrix name
        input = input.upper()
        if input not in align.SubstitutionMatrix.list_db():
            raise InputError(
                f"'{input}' is neither a file "
                f"nor a valid NCBI substitution matrix"
            )
        return align.SubstitutionMatrix(alphabet, alphabet, input)

def create_space(input):
    try:
        lightness = int(input)
    except ValueError:
        raise InputError("Value is not an integer")
    if lightness < 1 or lightness > 99:
        raise InputError("Value must be between 1 and 99") 
    return ColorSpace(lightness)

def adjust_saturation(input):
    global space
    try:
        min, max = input.split("-")
        min, max = int(min), int(max)
    except ValueError:
        raise InputError("The range is invalid")
    a = space.lab[:,:,1]
    b = space.lab[:,:,2]
    space.remove((a**2 + b**2 < min**2))
    space.remove((a**2 + b**2 > max**2))

def save_scheme(input):
    write_color_scheme(input, result, name)


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

def show_results(space, result):
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


class InputError(Exception):
    pass