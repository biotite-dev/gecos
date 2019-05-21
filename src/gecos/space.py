# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["ColorSpace"]

from os.path import join, dirname, realpath
import itertools
import numpy as np
from .colors import lab_to_rgb


SPACE_FILE_NAME = join(dirname(realpath(__file__)), "space.npy")


class ColorSpace():
    """
    Create a color space, that spans the complete *RGB* gamut of the *Lab*
    color space.
    The *Lab* components are discretized into integer values.

    A color space describes the boundaries of the color scheme to be generated.
    It uses the *Lab* color space.
    In addition to the inherent limit to colors, that can also be displayed in
    *RGB* space, further parts of the space can be removed by calling
    :func:`remove()`.

    Parameters
    ----------
    file_name : str, optional
        The path of a custom file to load the precalculated *RGB* convertible
        space from.

    Attributes
    ----------
    shape : tuple of int
        The shape of the space, i.e. the amount of *Lab* values in each
        dimension.
    lab : ndarray, shape=(100, 256, 256, 3), dtype=int
        The complete discretized *Lab* color space in the *ab* range of
        ``-128`` to ``127``.
    space : shape=(100, 256, 256), dtype=bool
        The allowed part of the `lab` attribute, i.e. the part that is
        convertible into *RGB* and was not manually removed.
    """

    def __init__(self, file_name=None):
        if file_name is None:
            file_name = SPACE_FILE_NAME
        
        l = np.arange(100)
        a = b = np.arange(-128,128)
        self._lab = np.zeros((100,256,256,3), dtype=int)
        self._lab[:,:,:,0] = l[:, np.newaxis, np.newaxis]
        self._lab[:,:,:,1] = a[np.newaxis, :, np.newaxis]
        self._lab[:,:,:,2] = b[np.newaxis, np.newaxis, :]

        with open(file_name, "rb") as file:
            self._space = np.load(file)

    def remove(self, mask):
        """
        Remove a portion of the color space.
        
        Parameters
        ----------
        space : ndarray, shape=(100, 256, 256), dtype=bool
            The space is removed where this mask is true.
        
        Examples
        --------
        Remove space below a defined lightness:

        >>> L_MIN = 50
        >>> space = ColorSpace()
        >>> lab = space.lab
        >>> l = lab[..., 0]
        >>> space.remove(l < L_MIN)
        """
        self._space &= ~mask
    
    def get_rgb_space(self):
        """
        Convert the *Lab* colors of the space in *RGB* colors.
        
        Returns
        -------
        rgb : ndarray, shape=(100,256,256,3), dtype=float
            The *RGB* colors.
            Colors that cannot be displayed in *RGB* or were manually
            removed from the space are ``NaN``.
        """
        rgb = lab_to_rgb(self._lab)
        rgb[~self._space] = np.nan
        return rgb

    @property
    def space(self):
        return self._space.copy()

    @property
    def shape(self):
        return (100, 256, 256)
    
    @property
    def lab(self):
        return self._lab.copy()

    @staticmethod
    def _generate(file_name=None):
        """
        Precalculate which part of the *Lab* space can be displayed in
        *RGB* and save the result as boolean mask into a
        *NumPy* ``.npy`` file.
        """
        lab = np.zeros((100, 256, 256, 3), dtype=int)
        lab[:,:,:,0] = np.arange(100      )[:, np.newaxis, np.newaxis]
        lab[:,:,:,1] = np.arange(-128, 128)[np.newaxis, :, np.newaxis]
        lab[:,:,:,2] = np.arange(-128, 128)[np.newaxis, np.newaxis, :]
        
        rgb = lab_to_rgb(lab)
        space = np.ones((100, 256, 256), dtype=bool)
        space[np.isnan(rgb).any(axis=-1)] = False
        
        if file_name is None:
            file_name = SPACE_FILE_NAME
        with open(file_name, "wb") as file:
            np.save(file, space)