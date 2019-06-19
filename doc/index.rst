.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

Gecos - Generated Color Schemes
===============================

Multiple sequence alignments are often visualized by coloring the symbols
according to some kind of properties.
For example a color scheme for amino acids could use one color for
hydrophobic residues, another color for positively charged residues
and so forth.
Usually, such color schemes are created manually by experienced people
who have knowledge about the characteristics of the e.g. amino acids,
so they can assign equal or similar colors to amino acids that share
similar properties.

The *Gecos* software follows a different approach:
Instead of looking at specific, sometimes subjective properties,
it uses another source for estimating the similarity of symbols:
the substitution matrix itself.
Similar colors are assigned to high scoring pairs of symbols, low
scoring pairs get distant colors - in a completely automatic manner.
As a result the distance of two symbols in the substitution matrix corresponds
to the perceptual differences in the color scheme.

How about an example?
The following command line invocation creates a light color scheme.
An example alignment using the newly generated color scheme is displayed below.

.. code-block:: console
   
   $ gecos --matrix BLOSUM62 --lmin 60 --lmax 75 -s awesome_colors.json

.. image:: /plots/main_example_alignment.png


.. toctree::
   :maxdepth: 2
   :hidden:
   
   intro
   cli
   api
   theory
   logo

