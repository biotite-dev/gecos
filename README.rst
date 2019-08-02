Gecos - Generated Color Schemes for sequence alignment visualizations
=====================================================================

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

Installation
============

In order to use *Gecos* you need to have Python (at least 3.6) installed.
Furthermore, the following Python packages are required:

   - **biotite**
   - **numpy**
   - **matplotlib**
   - **scikit-image**

If these prerequisites are met, *Gecos* is simply installed via

.. code-block:: console

   $ pip install gecos