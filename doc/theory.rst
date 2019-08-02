.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Background
==========

Color space
-----------

*Gecos* utilizes *CIE L\*a\*b\**, a color space that is designed to be
perceptually linear:
Changes of *L\*a\*b\** component values scale approximately linearly with the
visually perceived change of the corresponding color.
The *L\*a\*b\** color space contains the following components:

   - **L\*** - The lightness of the color. ``0`` is completely black and
     ``100`` is completely white.
   - **a\*** - The green-red component. Green is in the negative direction,
     red is in the positive direction.
   - **b\*** - The blue-yellow component. Blue is in the negative direction,
     yellow is in the positive direction.

The values for *a\** and *b\** are theoretically not limited in either
direction, but only a subspace (*gamut*) can be displayed on devices and can
be converted into *RGB*.

Due to the perceptual linearity, the perceptual difference of two *L\*a\*b\**
colors is approximately the euclidean distance of the *L\*a\*b\** components,
according to the *CIE76* formula.
Newer standards apply several corrections to the euclidean distance to deal
with perceptual non-uniformities.
However, *Gecos* still uses the simple euclidean distance as it can be
calculated very fast, which is crucial in the optimization process.

*Gecos* uses the package *scikit-image* for all kinds of color space
conversions.

.. _score_function: 

Score function
--------------

The score function :math:`S_T` is a compound of two terms:
a sum of harmonic potentials between each pair of symbols :math:`S_H`
and a linear *contrast score* :math:`S_C`:

.. math:: S_T = S_H + S_C

Note that a low score is desirable in this context.

.. note::
   
   In the following, :math:`\left< X \right>` denotes the arithmetic
   mean of all matrix elements
   (:math:`\left< X \right> = \frac{1}{n} \sum_{ij} X_{ij}`).
   For a triangular matrix only the non-zero diagonals are taken into account.

Construction of distance matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For both score terms it is required that the input substitution matrix
:math:`M` is converted into a triangular distance matrix :math:`D'`:

.. math:: D'_{ij} = \left( (S_{ii} - S_{ij}) + (S_{jj} - S_{ji}) / 2 \right)

For any substitution matrix :math:`M`, :math:`M_{ii}` should be the maximum
value in the row/column :math:`i`,
otherwise a symbol would be more similar to another symbol than to itself.
Thus, :math:`D'_{ii} = 0`.

:math:`D'` is scaled, so that the average distance is :math:`1`.

.. math:: D = \frac {D'} {\left< D' \right>} 

Perceptual difference matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On the other side the triangular matrix :math:`C` denotes the pairwise
perceptual differences of the *L\*a\*b\** color for each pair of symbols.
The perceptual difference is approximated as the euclidean distance:

.. math:: C_{ij} = \sqrt{(L^*_i - L^*_j)^2 + (a^*_i - a^*_j)^2 + (b^*_i - b^*_j)^2}

While :math:`D` is constant, :math:`C` changes with every step of the
optimization process.

In order to translate the *L\*a\*b\** color distances into values in the
distance matrix :math:`D`, a scale factor :math:`f_s` is introduced.
:math:`f_s` is the proportion of the average distance in :math:`D` to the
average difference in :math:`C`:

.. math:: f_s
   = \frac{\left< D \right>}{\left< C \right>}
   = \frac{ \frac{1}{n} \sum_{ij} D } { \frac{1}{n} \sum_{ij} C }
   = \frac{ \sum_{ij} D_{ij} } { \sum_{ij} C_{ij} }

As :math:`C` is variable, :math:`f_s` also dynamically changes.

Harmonic potentials
^^^^^^^^^^^^^^^^^^^

A harmonic potential applied all pairs of symbols with the equilibrium
position being the respective value from the distance matrix:

.. math:: S_H = \sum_{ij} \left( f_s C_{ij} - D_{ij} \right)^2

This adjusts the relative distance between all symbols.

Contrast score
^^^^^^^^^^^^^^

The *contrast score* rewards symbol conformations with a high contrast,
i.e. a high average perceptual difference between the symbols.
Like the harmonic potentials adjusts the relative distances between the
symbols, the *contrast score* tries to maximize the absolute distances.
A reciprocal function based on the average color differences is used here.
The *contrast factor* :math:`f_c` is a user-supplied parameter for weighting
this term:

.. math:: S_C = \frac{f_c}{\left< C \right>} 

Without this term, there would be no score difference between conformations,
that use a small or a large part of the color space, as long as the relative
distances are equal.
This term drives the symbols to the edges of the color
space, thereby increasing the contrast.

Optimization
------------

The purpose of *Gecos* is to find colors that match distances derived from a
substitution matrix as well as possible.
This means, that the software tries to optimize the *L\*a\*b\** values for all
symbols, so that the score function described above is minimized.
The *L\*a\*b\** values can be described as vector :math:`\vec{x}` with
:math:`n_s \times 3` dimensions, where :math:`n_s` is the amount of symbols
in the alphabet (e.g. 20 for amino acids). 

The optimization is performed via Metropolis-Monte-Carlo:
Starting from a random initial conformation :math:`\vec{x}_0` with a
score of :math:`S_0 = S_T(\vec{x}_0)`, the following
steps are performed:

   1) Perform random modifications on :math:`\vec{x}_n`:
      
      :math:`\vec{x}_{n+1} = f_M(\vec{x}_n)`

      :math:`f_M` is a function that adds a random value within a user-defined
      radius to :math:`\vec{x}`.
   
   2) Calculate the score of the new conformation:
      
      :math:`S_{n+1} = S_T(\vec{x}_{n+1})`
   
   3) Decide, whether to accept the new conformation based on the difference
      to the score of the conformation prior to modification:

      :math:`\Delta S = S_{n+1} - S_{n}`

      If :math:`\Delta S \leq 0`, then accept the new conformation.
      
      If :math:`\Delta S > 0`, then accept the new conformation with a
      probability of :math:`p = e^{ \frac{\Delta S}{T} }` where :math:`T`
      is the user-supplied temperature parameter.
      In case the new conformation is not accepted, the new conformation
      is replaced with the conformation prior to modification:

      :math:`\vec{x}_{n+1} = \vec{x}_n`

These steps are repeated until an acceptable score has been reached.

The command line interface uses a special variant, where the temperature is
stepwise decreased (simulated annealing).