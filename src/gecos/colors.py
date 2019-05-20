# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["rgb_to_lab", "lab_to_rgb"]

import warnings
import numpy as np
from skimage.color import rgb2lab, lab2rgb


def rgb_to_lab(color):
    return _convert_color(color, rgb2lab)

def lab_to_rgb(color):
    return _convert_color(color, lab2rgb)

def _convert_color(color, convert_func):
    if not isinstance(color, np.ndarray):
        color = np.array(color, dtype=float)
    else:
        color = color.astype(float, copy=False)
    original_shape = color.shape
    color = color.reshape([1,-1,3])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        conv_color = convert_func(color)
    conv_color[
        :,
        (conv_color[0] == 0).any(axis=-1) | (conv_color[0] == 1).any(axis=-1),
        :
    ] = np.nan
    conv_color = conv_color.reshape(original_shape)
    return conv_color