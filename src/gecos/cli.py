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
    print(
        "Hello, this is Gecos, your friendly color scheme generator "
        "for sequence alignments."
    )
    print(
        "Do you would like to create a color scheme for protein alignments?"
    )
    alphabet = parse_input(choose_alphabet, {"y":"yes", "n":"no"})
    if alphabet is None:
        print(
            "Then please provide an alphabet to generate the scheme for "
            "(e.g. 'ACGT'):"
        )
        alphabet = parse_input(parse_alphabet)
    print(alphabet)
    print(
        "And now I need the name of a substitution matrix file "
        "in BLAST format (e.g. 'blosum62.mat'):"
    )
    parse_input()
    print(
        "Which lightness should the color scheme have (1 - 99)"
    )
    parse_input()


def parse_input(parse_function=None, options=None):
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



class InputError(Exception):
    pass