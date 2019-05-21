# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["write_color_scheme"]

import json
from matplotlib.colors import to_hex


def write_color_scheme(file, result, name=""):
    """
    Write a color scheme in the *Biotite* compatible JSON format into
    a file.
    
    Parameters
    ----------
    file : file-like object
        The file to write the scheme into.
    result : Optimizer.Result
        The result from the optimization.
        Contains the color scheme.
    name : str, optional
        Name of the scheme. Will be written to the value of the
        ``"name"`` key.
    """
    scheme = {}
    scheme["name"] = name
    symbols = result.alphabet.get_symbols()
    scheme["alphabet"] = symbols
    scheme["colors"] = {s : to_hex(c) for s, c
                        in zip(symbols, result.rgb_colors)}
    json.dump(scheme, file, indent=4)