Gecos - Generated Color Schemes for sequence alignment visualizations
=====================================================================

Multiple sequence alignments are often visualized by coloring the symbols
according to some kind of properties.
For example a color scheme for amino acids could use one color for
hydrophobic residues, another color for positively charged residues
and so forth.
Usually, such color schemes are manually created by experienced people
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
   
   $ gecos --matrix BLOSUM62 --lmin 60 --lmax 75 -f awesome_colors.json

.. image:: https://raw.githubusercontent.com/biotite-dev/gecos/master/doc/static/assets/figures/main_example_alignment.png

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

Alternatively, *Gecos* can be installed via *Conda*:

.. code-block:: console

   $ conda install -c conda-forge gecos

Citation
--------

If you use *Gecos* in a scientific publication, please cite:

|

   P. Kunzmann, B. E. Mayer, K. Hamacher,
   "Substitution matrix based color schemes for sequence alignment visualization,"
   BMC Bioinformatics, vol. 21, pp. 209, 2020.
   doi: `10.1186/s12859-020-3526-6 <https://doi.org/10.1186/s12859-020-3526-6>`_