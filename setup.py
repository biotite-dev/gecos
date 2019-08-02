# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
import shlex
from src.gecos import __version__

long_description = """
Generated color schemes for sequence alignment visualizations
=============================================================

Multiple sequence alignments are often visualized by coloring the
symbols according to some kind of properties.
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
As a result the distance of two symbols in the substitution matrix
corresponds to the perceptual differences in the color scheme.
"""


class PyTestCommand(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to pytest")]

    def initialize_options(self):
        super().initialize_options()
        self.pytest_args = ''

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)


setup(
    name="gecos",
    version = __version__,
    description = \
        "Generated color schemes for sequence alignment visualizations",
    long_description = long_description,
    author = "Patrick Kunzmann",
    author_email = "padix.key@gmail.com",
    url = "",
    license = "BSD 3-Clause",
    classifiers = (
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ),
    
    zip_safe = False,
    packages = find_packages("src"),
    package_dir = {"" : "src"},
    entry_points = {
        "console_scripts": [
            "gecos = gecos.cli:main"
        ]
    },

    package_data = {
        "gecos" : ["space.npy", "example_alignment.fasta"],
    },
    
    install_requires = ["biotite",
                        "numpy",
                        "scikit-image"],
    python_requires = ">=3.6",
    
    cmdclass = {"test": PyTestCommand},
    tests_require = ["pytest"],
    
    command_options = {
        'build_sphinx':
            {"source_dir" : ("setup.py", "./doc"),
             "build_dir"  : ("setup.py", "./doc/_build"),
             "release"    : ("setup.py", __version__)}
    }
)
