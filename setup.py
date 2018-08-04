# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
import shlex
import glob
from os.path import join, abspath, dirname, normpath
import fnmatch
import os
from src.gecos import __version__

long_description = """
Generated color schemes.
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
    description = "Generated color schemes",
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
    
    install_requires = ["biotite",
                        "numpy",
                        "maplotlib"],
    python_requires = ">=3.5",
    
    cmdclass = {"test": PyTestCommand},
    tests_require = ["pytest"],
    
    command_options = {
        'build_sphinx':
            {"source_dir" : ("setup.py", "./doc"),
             "build_dir"  : ("setup.py", "./doc/_build"),
             "release"    : ("setup.py", __version__)}
    }
)
