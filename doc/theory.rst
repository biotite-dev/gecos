.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Background
==========

Color spaces
------------

The perceptual difference of two colors is approximately
the euclidean distance.

Potential function
------------------

The potential function :math:`V_T` is compound of two terms:
a sum of harmonic potentials between each pair of symbols :math:`V_H`
and a linear *contrast potential* :math:`V_C`:

.. math:: V_T = V_H + V_C

For both terms it is required that the input substitution matrix :math:`S` is
converted into a distance matrix :math:`D`:

.. math:: D_{ij} = S_{ii} - S_{ij}

For any substitution matrix :math:`S_{ii}`, should be the maximum value in the
row/column :math:`i`,
otherwise a symbol would be more similar to another symbol than to itself

On the other side the matrix :math:`C` denotes the pairwise perceptual
differences of the *L\*a\*b\** color for each pair of symbols.
The perceptual difference is approximated as the euclidean distance:

.. math:: C_{ij} = \sqrt{(L^*_i - L^*_j)^2 + (a^*_i - a^*_j)^2 + (b^*_i - b^*_j)^2}

While :math:`D` is constant, :math:`C` changes with every step of the
optimization process.

In order to translate the *L\*a\*b\** color distances into values in the
distance matrix :math:`D`, a scale factor :math:`f_s` is introduced.
:math:`f_s` is the proportion of the average distance in :math:`D` to the
average difference in :math:`C`:

.. math:: f_s = \frac{\left< D \right>}{\left< C \right>}

As :math:`C` is variable, :math:`f_s` also dynamically changes.

Harmonic potentials
^^^^^^^^^^^^^^^^^^^

A harmonic potential applied all pairs of symbols with the equilibrium
position being the respective value from the distance matrix:

.. math:: V_H = \sum_{ij} \left( f_s C_{ij} - D_{ij} \right)^2

This adjusts the relative distance between all symbols.

Contrast potential
^^^^^^^^^^^^^^^^^^

The *contrast potential* rewards symbol conformations with a high contrast,
i.e. a high average perceptual difference between the symbols.
Like the harmonic potentials adjusts the relative distances between the
symbols, the *contrast potential* tries to maximize the absolute distances.
The *contrast factor* :math:`f_c` is used for weighting of this term:

.. math:: V_C = -f_c \left< C \right>

This term drives the symbols to the edges of the color
space, thereby increasing the contrast.

Optimzation
-----------