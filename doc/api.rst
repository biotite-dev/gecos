.. This source code is part of the Gecos package and is distributed
   under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
   information.

.. currentmodule:: gecos

Python API
==========

Although the CLI can cope with most demands of *Gecos* users, there are some
aspects, where the CLI is not flexible enough.
This includes for example custom optimization processes or fine-tuned color
spaces.
This page explains the API of *Gecos'* underlying Python code. 


Colors
------

*Gecos* uses *scikit-image* to convert between *L\*\a\*b\** and *RGB* space.
However, the *scikit-image* function requires a specific shape for the
array of colors.
The following functions wrap the *scikit-image* functionality, generalizing
it for any shape.

.. autofunction:: rgb_to_lab

.. autofunction:: lab_to_rgb

In order to create a color scheme, the boundaries of the *L\*\a\*b\**
color space must be known.
For this purpose a :class:`ColorSpace` object is used.

.. autoclass:: ColorSpace

   .. automethod:: ColorSpace.remove

   .. automethod:: ColorSpace.get_rgb_space


Optimization
------------

A score function defines the objective of the color scheme generation.
It assigns a score to a given color conformation.
The color conformation is the entirety of the color component values for each
symbol, also called coordinates.
A more favorable color conformation gets a **lower** score than less favorable
one.

.. autoclass:: ScoreFunction

   .. automethod:: ScoreFunction.__call__()

The decision, whether a color conformation is favorable, is subject to the
actual implementation of the :class:`ScoreFunction`.
The mathematical details for the default score function are covered
:ref:`here <score_function>`.

.. autoclass:: DefaultScoreFunction

The :class:`ColorOptimizer` performs the actual optimization of the color
conformation.
This means, that it tries to adjust the coordinates, so that the value of the
score function gets minimal.

.. autoclass:: ColorOptimizer

   .. automethod:: ColorOptimizer.optimize()

   .. automethod:: ColorOptimizer.set_coordinates()

   .. automethod:: ColorOptimizer.get_result()

After optimizing the color conformation using
:func:`ColorOptimizer.optimize()`, the result of the optimization
is obtained via :func:`ColorOptimizer.get_result()`.

.. autoclass:: gecos::ColorOptimizer.Result


Color scheme output
-------------------

The color scheme from the :class:`ColorOptimizer.Result` can be written into
the *Biotite* compatible JSON format via :func:`write_color_scheme`.

.. autofunction:: write_color_scheme