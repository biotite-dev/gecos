# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import itertools
import pytest
import numpy as np
import gecos


np.random.seed(0)
@pytest.mark.parametrize(
    "rgb_color, ndim",
    itertools.product(np.random.rand(100, 3), [1,2,3,4])
)
def test_color_conversion(rgb_color, ndim):
    if ndim != 1:
        new_shape = tuple([1] * (ndim-1) + [3])
        rgb_color = rgb_color.reshape(new_shape)
    conv_rgb_color = gecos.lab_to_rgb(gecos.rgb_to_lab(rgb_color))
    assert np.allclose(conv_rgb_color, rgb_color)