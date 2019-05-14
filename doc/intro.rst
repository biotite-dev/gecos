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
The gray area consists of *L\*a\*b\** values that cannot be converted into
*RGB* space.

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


Basic usage
-----------

The most simple invocation is simply

.. code-block:: console
   
   $ gecos

By default *Gecos* uses the *BLOSUM62* matrix to generate a color scheme, which
is printed to console.
Alternatively, You can save the color scheme into a file via the ``-s`` option.
The color scheme is printed in a *Biotite* compatible JSON format and will
look something like this:

.. code-block:: json
   
    {
        "name": "scheme",
        "alphabet": ["A","C","D","E","F","G","H","I","K","L"
                     "M","N","P","Q","R","S","T","V","W","Y"],
        "colors": {
            "A": "#7c7b8b",
            "C": "#17ebd9",
            "D": "#740365",
            "E": "#992651",
            "F": "#f3df8c",
            "G": "#140a1a",
            "H": "#b41308",
            "I": "#e8eafe",
            "K": "#fe83aa",
            "L": "#f0eee6",
            "M": "#fcdbce",
            "N": "#d0388b",
            "P": "#ba82fd",
            "Q": "#873429",
            "R": "#fe7878",
            "S": "#744759",
            "T": "#4c5e53",
            "V": "#afcbe0",
            "W": "#d5e70b",
            "Y": "#aa7e00"
        }
    }

The value of ``"name"`` is obviously the given name of the color scheme.
It can be adjusted with the ``--name`` option.
The ``"alphabet"`` maps to a list of symbols comprised by the alphabet the
scheme is intended for.
The most important field is ``"colors"``:
It maps to a dictionary, where a color is assigned to each symbol of the
alphabet.
Even though *Gecos* assigns a color to all symbols in ``"alphabet"``,
the format allows that ``"colors"`` assigns colors only to a subset of the
symbols in alphabet.

.. note::
   
   Although the format is compliant with the *Biotite* color scheme format,
   the *Biotite* amino acid alphabet contains additional symbols for the
   ambiguous amino acids and the stop codon.
   Hence incorporating a *Gecos* color scheme into *Biotite* requires that
   that the symbols ``"B"``, ``"Z"``, ``"X"`` and ``"*"`` are appended at the
   end of the ``"alphabet"`` value.
   Editing ``"colors"`` is not necessary.

As the color space was not restricted in any way, the generated color scheme
contains the whole lightness range -  from pitch-black to pure white.
Alignments visualized with this color scheme look accordingly:

.. image:: /plots/no_constraints_scheme_alignment.png

Although this scheme has a high contrast and the color differences are well
aligned with the substitution matrix, such a wide lightness range is seldom
intended.
To constrain the lightness range, you can give *Gecos* a minimum and a
maximum lightness level:

.. code-block:: console
   
   $ gecos --lmin 60 --lmax 75 -s a_color_scheme.json

.. image:: /plots/main_example_alignment.png

However, the minimum and the maximum lightness should not be too close, lest
the contrast will be quite low.

Color constraints
-----------------

The *a\** and *b\** components can be restrained in the same way, to create
a color scheme that is shifted into a certain hue.
This can, for example, be used to create a color scheme for red-green deficient
people.
For this purpose the green region will be removed, i.e. *a\** starts at
``0``.
In order to compensate for the lost contrast, the lightness range is increased:

.. code-block:: console

   $ gecos --amin 0 --lmin 50 --lmax 80 -s no_green_scheme.json

.. image:: /plots/no_green_scheme_alignment.png

Likewise the saturation range can be set.
The saturation is the euclidean distance of the *a\*b\** components to
gray (``0``, ``0``):

.. code-block:: console

   $ gecos --smin 30 --lmin 55 --lmax 75 -s saturated_scheme.json

.. image:: /plots/high_saturation_scheme_alignment.png

Last but not least, you can constrain a symbol to a specfic *L\*a\*b\** color
via the ``--constraint`` or ``-c`` option.
The optimization will not change the color of constrained symbols
In the following example, we want *alanine* to be gray and *tryptophane* to be
blue, both with a lightness of ``70``:

.. code-block:: console

   $ gecos -c A 70 0 0 -c W 70 -20 -40 --lmin 60 --lmax 75 -s constrained_scheme.json

.. image:: /plots/constrained_scheme_alignment.png

Adjusting the contrast
----------------------

*Gecos'* optimization process contains an additional potential that penalizes
low contrast color conformations, i.e. average low distances between the
symbols.
This behavior can be customized by setting the ``--contrast`` option.
When the value is ``0``, low contrast schemes are not penalized.
The higher the value, the more the symbols are driven to the edges of the
color space.
However, the magnitude of the value required to obtain a certain effect is
strongly dependent on the size of the alphabet and the substitution matrix.
Thus, a bit of experimentation is necessary to find an optimal value for this
option.
The following example creates a high contrast color scheme:

.. code-block:: console
   
   $ gecos --contrast 100 --lmin 60 --lmax 75 -s high_contrast_scheme.json

.. image:: /plots/high_contrast_scheme_alignment.png

.. note::
   
   Use the ``--contrast`` parameter with caution.
   Increasing the contrast parameter also means, that the substitution matrix
   based distance is weighted less strongly.
   Consequently, although a high contrast color scheme may look appealing,
   it also may not represent the similarity of symbols very well.

Color space and scheme preview
------------------------------

You do not need to create an alignment yourself in order to evaluate a newly
created color scheme.
*Gecos* provides some visualization capabilities by itself, so you can
directly discard a color scheme you do not like.

At first, you can output your selected color space with the ``--show-space``
option.
The additional ``--dry-run`` terminates the program after the color space
has been displayed:

.. code-block:: console
   
   $ gecos --show-space --dry-run --smin 30 --lmin 60 --lmax 70

.. image:: /plots/show_space.png

The plot is a 2D projection of the color space at a fixed lightness.
The lightness value in the plot is the average of the ``--lmin`` and the
``--lmax`` value.
The displayed lightness value can be customized with the ``--lightness``
option.
The *hole* in the center of the plot is causes by the saturation constraint.

The ``--show-scheme`` option shows the symbol conformation in color space
after the optimization.
Again the plot is a 2D projection at a fixed lightness.
The white area shows the allowed color space at the given lightness:

.. code-block:: console
   
   $ gecos --show-scheme --smin 30 --lmin 60 --lmax 70

.. image:: /plots/show_scheme.png

Some symbols might seem to be outside of the allowed space, but remember
that the white area is only the allowed space at the given lightness.

Finally, the ``--show-example`` options shows an example multiple protein
sequence alignment with the color scheme.

.. code-block:: console
   
   $ gecos --show-example --smin 30 --lmin 60 --lmax 70

.. image:: /plots/show_example.png

Custom matrices and alphabets
-----------------------------

While the default substitution matrix *Gecos* uses is *BLOSUM62* you can also
use a custom substitution matrix.
Either a valid NCBI substitution matrix name (e.g. `PAM250`) or a custom matrix
file in NCBI format can be supplied to the ``--matrix`` option.
Likewise, it is possible to generate a color scheme for a different alphabet
than the default amino acid alphabet, by setting the ``--alphabet`` option.

In order to demonstrate this the following example will generate a color scheme
for the *protein blocks* (PB) alphabet
(`de Brevern et al., 2000 <https://doi.org/10.1002/1097-0134(20001115)41:3\<271::AID-PROT10\>3.0.CO;2-Z>`_ ).
*protein blocks* consists of 16 symbols, where each one represents another
protein backbone conformation.
In a nutshell PBs can be used to encode a molecular 3D structure into a
sequence.