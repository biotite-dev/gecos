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
   Hence, the amount of possible :math:`ij` pairings is equal to
   :math:`n = \frac{n_s (n_s - 1)} {2}`. :math:`n_s` is the number of symbols
   in the alphabet.


Construction of distance matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For both score terms it is required that the input substitution matrix
:math:`M` is converted into a triangular distance matrix :math:`D'`:

.. math:: D'_{ij} = \left( (S_{ii} - S_{ij}) + (S_{jj} - S_{ji}) \right) / 2 

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
   = \frac{ 1 } { \frac{1}{n} \sum_{ij} C }
   = \frac{ n } { \sum_{ij} C_{ij} }

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
substitution matrix as good as possible.
This means, that the software tries to optimize the *L\*a\*b\** values for all
symbols, so that the score function described above is minimized.
The *L\*a\*b\** values can be described as vector :math:`\vec{x}` with
:math:`n_s \times 3` dimensions, where :math:`n_s` is the amount of symbols
in the alphabet (e.g. 20 for amino acids). 

The optimization is performed via *simulated annealing* (SA).
The SA algorithm is basically an improved Monte-Carlo
optimization algorithm, which means sampling the solution space of the 
optimization problem with a given temperature :math:`T` or
rather an inverse temperature :math:`\beta` where :math:`\beta \sim 1/T`.

The improvement of SA over a simple 
Monte-Carlo optimization algorithm is to perform the optimization with an 
initially high temperature, or low inverse temperature accordingly, which 
is continuously cooled down over the course of the algorithms runtime.
The idea here comes from the physical process of annealing of, e.g., 
steel where you can make the observation that a slowly 
cooled steel has superior material characteristics.

The cooling down is steered by an annealing schedule which in our case is 
the exponential schedule, so we have

.. math:: \beta(t) = \beta_0 \cdot \exp \left( \tau \cdot t \right).
     
Furthermore, as SA is usually employed for combinatorial 
optimization problems, so problems defined on discrete space, we also use 
an exponential schedule for the step size 
    
.. math:: \delta(n) = \delta_0 \cdot \exp \left( \gamma \cdot t \right).
    
The step size is used for perturbing the current solution in each step of the
SA algorithm to find a new candidate solution.
So the idea for using the schedule here is to start with relatively large 
step size :math:`\delta_{start}` and to chose the rate  according to an 
target step size :math:`\delta_{end}`.
An according rate is easily derived  by claiming
:math:`\delta(N_{max})=\delta_{end}` which leads to

.. math:: \gamma = \frac{1}{N_{max}}\log \left( \frac{\delta_{end}}{\delta_{start}} \right).
 

Monte-Carlo algorithm
^^^^^^^^^^^^^^^^^^^^^

Starting from a random initial conformation :math:`\vec{x}_0` with a
score of :math:`S_0 = S_T(\vec{x}_0)`, the following
steps are performed:

   1) Perform random modifications on :math:`\vec{x}_n`:
      
      :math:`\vec{x}_{n+1} = \vec{x}_n + \Delta(\vec{x}_n)`

      where :math:`\Delta(\vec{x}_n)` is a random perturbation calculated using
      the step size :math:`\delta(n)`. 
  
   2) Calculate the score of the new conformation:
      
      :math:`S_{n+1} = S_T(\vec{x}_{n+1})`
                
   
   3) Decide, whether to accept the new conformation based on the difference
      to the score of the conformation prior to modification:

      :math:`\Delta S = S_{n+1} - S_{n}`

      If :math:`\Delta S \leq 0`, then accept the new conformation.
      
      If :math:`\Delta S > 0`, then accept the new conformation with a
      probability of 
      :math:`p = exp \left( \beta(n) \cdot \Delta S \right)` where :math:`\beta(n)`
      is the inverse temperature according to the exponential annealing schedule.
      
      
      The initial inverse temperature :math:`\beta_0` as 
      well as the rate :math:`\tau`, specifying how
      fast the inverse temperature increases, can be user specified.
      In case the new conformation is not accepted, the new conformation
      is replaced with the conformation prior to modification:

      :math:`\vec{x}_{n+1} = \vec{x}_n`

These steps are repeated until an stop criterion is met, which is just a fixed
number of iterations in this case.


