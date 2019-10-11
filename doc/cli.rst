.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Command line interface
======================

The command line interface (CLI) ist the common way to use *Gecos*.
It offers a range of options to customize the way color schemes are
generated.


Color space
-----------

The color space, in which the color scheme is optimized can be limited.
The range of *L\*a\*b\** values, as well as the minimum and maximum saturation
can be set.
The saturation :math:`s` is the euclidean distance from gray at the same
lightness:

.. math:: s = \sqrt{{a^*}^2 + {b^*}^2}

The following table shows which CLI option to use in order to set the lower or
upper limit of a color space value.

+-------+-------------+-------------+
|       | Lower limit | Upper limit |
+=======+=============+=============+
| *L\** | ``--lmin``  | ``--lmax``  |
+-------+-------------+-------------+
| *a\** | ``--amin``  | ``--amax``  |
+-------+-------------+-------------+
| *b\** | ``--bmin``  | ``--bmax``  |
+-------+-------------+-------------+
| *s*   | ``--smin``  | ``--smax``  |
+-------+-------------+-------------+

If you would like to display only the customized color space and skip the color
scheme generation use the ``--dry-run``/``-n`` option.


Substitution matrix and alphabet
--------------------------------

An alphabet is defined as the set of symbols a sequence may contain.
For example the unambiguous nucleobase alphabet contains the four symbols
``A``, ``C``, ``G`` and ``T``.

By default *Gecos* creates a color scheme for the amino acid alphabet based on
the *BLOSUM62* substitution matrix.

To change the substitution matrix, either a valid NCBI substitution matrix name
(e.g. *PAM250*) or a path to a text file containing a substitution matrix in
NCBI format can be given to the ``--matrix``/``-m`` option.

To create a color scheme for another alphabet
(e.g. nucleobases or a structural alphabet), the ``--alphabet``/``-a`` option
needs to be set. The option takes a string, where each character is treated as
individual symbol (e.g. ``"ACGT"`` for nucleobases).
Therefore, multi-character symbols are not allowed.
For multi-character symbols refer to the :doc:`Python API </api>` instead.


Optimization
------------

Based on the color space, the substitution matrix and the alphabet,
the optimizer tries to find a color scheme that optimally matches the matrix.
In order to increase the quality of the scheme tha amount of optimization steps
(``--nsteps``) or the number of runs (``--nruns``) can be increased.
However, increasing these values also extends the runtime of the optimization.
Note that ``--nruns`` can take advantage of multiple cores.

The simulated annealing can be adjusted even more fine grained by setting
the initial reverse temperature (``--beta``) and the rate of its exponential
growth (``--rate``). The step size decreases in the course of the simulated
annealing also in an exponential manner, which can be parameterized via
``--step-size-start`` and ``--step-size-end``.
The seed for the random number generator used by the algorithm is set with
the ``--seed`` option.
However, these parameters address the more advanced users.

A color can be fixed for a certain symbol by giving a
(symbol, *L\**, *a\**, *b\**) tuple to the ``--constraint``/``-c`` option
(e.g. ``A 50 0 0``).
The option can be repeated to fix the color for multiple symbols.
In addition to the relative color difference the optimizer also tries to find
a scheme with a high contrast. The weighting of the importance of the contrast
is adjusted with ``--contrast``.

The formula used for calculation of perceptual color differences
(:math:`\Delta E`) can be set via the ``--delta`` option.
While the default formula ``'CIEDE2000'`` yields the most accurate results,
``'CIE76'`` features the fastest calculation.


Output files
------------

If no output files are specified, *Gecos* simply outputs the generated color
scheme in JSON format to *STDOUT*.
Alternatively, the scheme can be directly saved to a file specified via
``--scheme-file``/``-f``.
The name of the scheme inside the file is set with ``--name``.
Optionally, the score values throughout the simulation can be saved via
``--score-file`` for further analysis.

Visualization
-------------

*Gecos* provides some visualization tools for direct evaluation of a color
space or a generated color scheme:

+---------------------+----------------------------------------------------------------------------------+
| ``--show-space``    | Shows the color space, including the user-supplied limitations.                  |
+---------------------+----------------------------------------------------------------------------------+
| ``--show-scheme``   | Shows conformation of symbols in the selected color space.                       |
+---------------------+----------------------------------------------------------------------------------+
| ``--show-score``    | Shows the score values throughout the simulation.                                |
+---------------------+----------------------------------------------------------------------------------+
| ``--show-example``  | Shows an example multiple protein sequence alignment using the new color scheme. |
|                     | Cannot be combined with a custom alphabet.                                       |
+---------------------+----------------------------------------------------------------------------------+

``--show-space`` and ``--show-scheme`` show a 2D *a\*b\** cross section of the
3D color space.
The *L\** value of this section is set with ``--lightness``/``-l``.