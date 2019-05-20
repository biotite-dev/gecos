# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["rgb_to_lab", "lab_to_rgb"]

import warnings
import numpy as np
from skimage.color import rgb2lab, lab2rgb


def rgb_to_lab(rgb):
    """
    Convert one or multiple *RGB* color(s) to *Lab* color(s).
    
    Parameters
    ----------
    rgb : array-like, shape=(..., 3)
        The color to be converted.
        The only restriction on the array shape is that the length of
        the last dimension must be 3.
        *RGB* values are in range (0,1).
    
    Returns
    -------
    lab : ndarray
        The converted color.
    """
    return _convert_color(rgb, rgb2lab)

def lab_to_rgb(lab):
    """
    Convert one or multiple *Lab* color(s) to *RGB* color(s).
    
    Parameters
    ----------
    lab : array-like, shape=(..., 3)
        The color to be converted.
        The only restriction on the array shape is that the length of
        the last dimension must be 3.
    
    Returns
    -------
    rgb : ndarray
        The converted color. *RGB* values are in range (0,1).
    """
    return _convert_color(lab, lab2rgb)

def _convert_color(color, convert_func):
    if not isinstance(color, np.ndarray):
        color = np.array(color, dtype=float)
    else:
        color = color.astype(float, copy=False)
    original_shape = color.shape
    # Temporarily reshape into a scikit-image-accepted shape
    color = color.reshape([1,-1,3])
    # Convert color(s) using scikit-image
    # Ignore warnings about colors being not convertible
    # - they will be set later to NaN anyway
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        conv_color = convert_func(color)
    # Color components that cannot be properly converted are clamped
    # to (0/1) -> Use this to detect invalid values
    conv_color[
        :,
        (conv_color[0] == 0).any(axis=-1) | (conv_color[0] == 1).any(axis=-1),
        :
    ] = np.nan
    conv_color = conv_color.reshape(original_shape)
    return conv_color