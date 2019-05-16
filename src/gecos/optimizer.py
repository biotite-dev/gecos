# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["ColorOptimizer", "PotentialFunction", "DefaultPotentialFunction"]

from collections import namedtuple
import abc
import copy
import numpy as np
import numpy.random as random
import biotite.sequence as seq
import biotite.sequence.align as align
from .colors import lab_to_rgb


MIN_L = 0
MAX_L = 99
MIN_AB = -128
MAX_AB = 127


class ColorOptimizer(object):

    class Result(namedtuple("Result", ["alphabet", "trajectory", "potentials"])):

        @property
        def coord(self):
            return copy.deepcopy(self.trajectory[-1])
        
        @property
        def potential(self):
            return self.potentials[-1]

        @property
        def lab_colors(self):
            return copy.deepcopy(self.trajectory[-1])
        
        @property
        def rgb_colors(self):
            return lab_to_rgb(self.lab_colors.astype(int))
    
    def __init__(self, alphabet, potential_function, space, constraints=None):
        self._alphabet = alphabet
        self._n_symbols = len(alphabet)
        self._pot_func = potential_function
        self._space = space.space.copy()
        self._coord = None
        self._trajectory = []
        self._potentials = []

        if constraints is None:
            self._constraints = np.full((self._n_symbols, 3), np.nan)
        else:
            for constraint in constraints:
                if not np.isnan(constraint).any() and \
                   not self._is_allowed(constraint):
                    raise ValueError(
                        f"Constraint {constraint} is outside the allowed space"
                    )
            self._constraints = constraints.copy()

        ### Set initial conformation ###
        # Every symbol has the 'l', 'a' and 'b' coordinates
        # The coordinates are initially filled with values
        # that are guaranteed to be invalid (l cannot be -1)
        start_coord = np.full((self._n_symbols, 3), -1, dtype=float)
        # Chose start position from allowed positions at random
        for i in range(start_coord.shape[0]):
            while not self._is_allowed(start_coord[i]):
                drawn_coord = random.rand(3)
                drawn_coord[..., 0]  *= (MAX_L -MIN_L ) + MIN_L
                drawn_coord[..., 1:] *= (MAX_AB-MIN_AB) + MIN_AB
                start_coord[i] = drawn_coord
        self._apply_constraints(start_coord)
        self._set_coordinates(start_coord)

    def set_coordinates(self, coord):
        if coord.shape != (self._n_symbols, 3):
            raise ValueError(
                f"Given shape is {coord.shape}, "
                f"but expected shape is {(len(self._alphabet), 3)}"
            )
        for c in coord:
            if not self._is_allowed(c):
                raise ValueError(
                    f"Coordinates {c} are outside the allowed space"
                )
        coord = coord.copy()
        self._apply_constraints(coord)
        self._set_coordinates(coord)
    
    def _set_coordinates(self, coord, potential=None):
        self._coord = coord
        self._trajectory.append(coord)
        if potential is None:
            potential = self._pot_func(coord)
        self._potentials.append(potential)
    
    def optimize(self, n_steps, temp, step_size):
        for i in range(n_steps):
            pot = self._potentials[-1]
            new_coord = self._move(self._coord, step_size)
            new_pot = self._pot_func(new_coord)
            if new_pot < pot:
                self._set_coordinates(new_coord, new_pot)
            else:
                p = np.exp(-(new_pot-pot) / temp)
                if p > random.rand():
                    self._set_coordinates(new_coord, new_pot)
                else:
                    self._set_coordinates(self._coord, new_pot)

    def get_result(self):
        trajectory = np.array(self._trajectory)
        return ColorOptimizer.Result(
            alphabet = self._alphabet,
            trajectory = trajectory,
            potentials = np.array(self._potentials)
        )
    
    def _is_allowed(self, coord):
        if coord[0] < MIN_L  or coord[0] > MAX_L  or \
           coord[1] < MIN_AB or coord[1] > MAX_AB or \
           coord[2] < MIN_AB or coord[2] > MAX_AB:
                return False
        # Add sign to ensure the corresponding integer value
        # has an absolute value at least as high as the floating value
        # This ensures that no unallowed values
        # are classified as allowed
        return self._space[
            int(coord[0]) - MIN_L,
            int(coord[1]) - MIN_AB,
            int(coord[2]) - MIN_AB,
        ]
    
    def _move(self, coord, step):
        new_coord = coord + (random.rand(*coord.shape)-0.5) * 2 * step
        self._apply_constraints(new_coord)
        # Resample coordinates for alphabet symbols
        # when outside of the allowed area
        for i in range(new_coord.shape[0]):
            while not self._is_allowed(new_coord[i]):
                new_coord[i] = coord[i] + (random.rand(3)-0.5) * 2 * step
        return new_coord
    
    def _apply_constraints(self, coord):
        mask = ~(np.isnan(self._constraints).any(axis=-1))
        coord[mask] = self._constraints[mask]


class PotentialFunction(metaclass=abc.ABCMeta):

    def __init__(self, n_symbols):
        self._n_symbols = n_symbols

    @abc.abstractmethod
    def __call__(self, coord):
        if len(coord) != self._n_symbols:
            raise ValueError(
                f"Expected {self._n_symbols} coordinates, but got {len(coord)}"
            )


class DefaultPotentialFunction(PotentialFunction):

    def __init__(self, matrix, contrast=100):
        if not matrix.is_symmetric():
            raise ValueError("Substitution matrix must be symmetric")
        super().__init__(len(matrix.get_alphabet1()))
        self._matrix = self._calculate_distance_matrix(matrix)
        self._matrix_sum = np.sum(self._matrix)
        # Scale contrast factor internally
        # so the user does not need to type hight numbers
        self._contrast = contrast * 1000
    
    def __call__(self, coord):
        super().__call__(coord)
        dist = np.sqrt(
            np.sum(
                (coord[:, np.newaxis, :] - coord[np.newaxis, :, :])**2, axis=-1
            )
        )
        dist = np.tril(dist)
        dist_sum = np.sum(dist)
        # This factor translates visual distances
        # into substitution matrix distances
        scale_factor = self._matrix_sum / dist_sum
        # Harmonic potentials between each pair of symbols
        harmonic_pot = np.sum((dist*scale_factor - self._matrix)**2)
        # Contrast term: Favours conformations
        # with large absolute color differences
        # 'where=dist' includes all non-zeroes
        contrast_pot = self._contrast / dist_sum
        return harmonic_pot + contrast_pot
        
    @staticmethod
    def _calculate_distance_matrix(similarity_matrix):
        scores = similarity_matrix.score_matrix()
        diff_to_max = np.diag(scores) - scores
        distances = np.tril((diff_to_max + diff_to_max.T) / 2)
        # Side length of the triangular matrix
        n = len(scores) - 1
        norm_factor = n/2 * (n+1)
        # Scale, so that average distance is 1
        distances = distances / (np.sum(distances) / norm_factor)
        return distances
