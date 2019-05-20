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
|       | Lower limit | Upper Limit |
+=======+=============+=============+
| *L\** | ``--lmin``  | ``--lmax``  |
+-------+-------------+-------------+
| *a\** | ``--amin``  | ``--amax``  |
+-------+-------------+-------------+
| *b\** | ``--bmin``  | ``--bmax``  |
+-------+-------------+-------------+
| *s*   | ``--smin``  | ``--smax``  |
+-------+-------------+-------------+

If you would like to display only the customized color and skip the color
scheme generation use the ``--dry-run``/``-n`` option.


Substitution matrix and alphabet
--------------------------------

An alphabet is defined as the set of symbols a sequence may contain.
For the unambiguous nucleobase alphabet contains the four symbols
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
(``--nsteps``) or the number of parallel optimization starts (``--nparallel``)
can be increased.
However, increasing these values also extends the runtime of the optimization.
Note that ``--nparallel`` can take advantage of multiple cores.

A color can be fixed for a certain symbol by giving a
(symbol, *L\**, *a\**, *b\**) tuple to the ``--constraint``/``-c`` option
(e.g. ``A 50 0 0``).
The option can be repeated to fix the color for multiple symbols.



Output files
------------

Visualization
-------------