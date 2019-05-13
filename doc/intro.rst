.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Getting started
===============

*Gecos* generates a color scheme by performing a Metropolis-Monte-Carlo
optimization in color space.
In short it means that the algorithm tries to assign colors to the symbols
(e.g. amino acids), whose pairwise perceptual differences is proportional to
the respective distances calculated from a substitution matrix.

There are dozens of different color spaces with *RGB* probably being the most
common one.
Despite its popularity, the *RGB* color space does not do well when it comes
to perceptual linearity:
Changing an *RGB* color value by a particular amount does not result in a
visual difference of the same amount. 
Due to this issue *Gecos* uses the *CIE L\*a\*b\** color space instead, that
behaves perceptually approximately linear.
The color space consists of three components:

   - **L\*** - The lightness of the color. ``0`` is completely black and
     ``100`` is completely white.
   - **a\*** - The green-red component. Green is in the negative direction,
     red is in the positive direction.
   - **b\*** - The blue-yellow component. Blue is in the negative direction,
     yellow is in the positive direction.

While values for *a\** and *b\** are not limited in either direction,
only a small space is displayable and hence can be converted into *RGB* colors.
Consequently the optimization process is also restricted to the displayable
subspace.
The following plots show the displayable *a\*b\** space at two different
*L\** levels.
The gray area are *L\*a\*b\** values that cannot be converted into *RGB* space.

.. image:: /plots/example_space.png


Installation
------------

In order to use *Gecos* you need to have Python (at least 3.6) installed.
Furthermore, the following Python packages are required:

   - **biotite**
   - **numpy**
   - **matplotlib**
   - **colormath**

If these prerequisites are met, *Gecos* is simply installed via

.. code-block:: console

   $ pip install gecos


Usage
-----