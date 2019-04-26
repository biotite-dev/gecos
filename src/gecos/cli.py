import os
import copy
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
    print(alphabet)
    
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
        figure = plt.figure()
        ax = figure.add_subplot(111)
        rgb_matrix = space.get_rgb_matrix()
        rgb_matrix[np.isnan(rgb_matrix)] = 0.7
        ax.imshow(np.transpose(rgb_matrix, axes=(1,0,2)), origin="lower",
                  extent=(-128, 127,-128, 127), aspect="equal")
        ax.set_xlabel("a")
        ax.set_ylabel("b")
        plt.show(block=False)

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
    except Exception as e:
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


class InputError(Exception):
    pass