.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Background
==========

Color spaces
------------

The perceptual difference of two colors is approximately
the euclidean distance

Scoring function
----------------

At first the selected substitution matrix is converted into a distance
matrix.
Then a harmonic potential applied all pairs of symbols with the  being the
respective value from the distance matrix
An additional *contrast potential* drives the symbol to the edges of the color
space, thereby increasing the contrast.

Optimzation
-----------