from collections import namedtuple
import copy
import numpy as np
import numpy.random as random
import biotite.sequence as seq
import biotite.sequence.align as align
from .colors import convert_lab_to_rgb



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
            return convert_lab_to_rgb(self.lab_colors)
    
    def __init__(self, matrix, space, constraints=None, contrast=10):
        self._lightness = space.lightness
        self._space = space.space.copy()
        self._coord = None
        self._trajectory = []
        self._potentials = []
        
        if not matrix.is_symmetric():
            raise ValueError("Substitution matrix must be symmetric")
        self._alphabet = matrix.get_alphabet1()
        distance_matrix = self.calculate_distace_matrix(matrix)

        if constraints is None:
            self._constraints = np.full((len(self._alphabet), 2), np.nan)
        else:
            for constraint in constraints:
                if not np.isnan(constraint).any() and \
                   not self._is_allowed(constraint):
                    raise ValueError(
                        f"Constraint {constraint} is outside the allowed space"
                    )
            self._constraints = constraints.copy()
        
        ### Potential parameters ###
        # Under optimal conditions the distances in the simulation
        # should be equal to distances calculated from the
        # substitution matrix
        self._dist_opt = distance_matrix
        # The average optimal distance is used to relate the visual
        # distances to the distances in the substitution matrix
        self._mean_dist_opt = np.mean(self._dist_opt)
        # The contrast factor is used to force a distribution into the
        # edges of the color space
        self._contrast = contrast

        ### Set initial conformation ###
        MIN_AB = -128
        MAX_AB = 127
        # Every symbol has the 'a' and 'b' coordinates
        # The coordinates are initially filled with values
        # that are guaranteed to be invalid
        start_coord = np.full((len(self._alphabet), 2), MIN_AB-1, dtype=float)
        # Chose start position from allowed positions at random
        for i in range(start_coord.shape[0]):
            while not self._is_allowed(start_coord[i]):
                start_coord[i] = random.rand(2) * (MAX_AB-MIN_AB) + MIN_AB
        self._apply_constraints(start_coord)
        self._set_coordinates(start_coord)

    def set_coordinates(self, coord):
        if coord.shape != (len(self._alphabet), 2):
            raise ValueError(
                f"Given shape is {coord.shape}, "
                f"but expected shape is {(len(self._alphabet), 2)}"
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
            potential = self._potential_function(coord)
        self._potentials.append(potential)
    
    def set_seed(self, seed):
        random.seed(seed=seed)
    
    def optimize(self, n_steps, temp, step_size):
        for i in range(n_steps):
            pot = self._potentials[-1]
            new_coord = self._move(self._coord, step_size)
            new_pot = self._potential_function(new_coord)
            if new_pot < pot:
                self._set_coordinates(new_coord, new_pot)
            else:
                p = np.exp(-(new_pot-pot) / temp)
                if p > random.rand():
                    self._set_coordinates(new_coord, new_pot)
                else:
                    self._set_coordinates(self._coord, new_pot)

    def get_result(self):
        trajectory = np.zeros(
            (len(self._trajectory), len(self._alphabet), 3)
        )
        trajectory[:, :, 0 ] = self._lightness
        trajectory[:, :, 1:] = np.array(self._trajectory)
        return ColorOptimizer.Result(
            alphabet = self._alphabet,
            trajectory = trajectory,
            potentials = np.array(self._potentials)
        )
    
    def _is_allowed(self, coord):
        # Add sign to ensure the corresponding integer value
        # has an absolute value at least as high as the floating value
        # This ensures that no unallowed values
        # are classified as allowed
        return self._space[
            int(coord[0] + np.sign(coord[0])) + 128,
            int(coord[1] + np.sign(coord[1])) + 128
        ]
    
    def _move(self, coord, step):
        new_coord = coord + (random.rand(*coord.shape)-0.5) * 2 * step
        self._apply_constraints(new_coord)
        # Resample coordinates for alphabet symbols
        # when outside of the allowed area
        for i in range(new_coord.shape[0]):
            while not self._is_allowed(new_coord[i]):
                new_coord[i] = coord[i] + (random.rand(2)-0.5) * 2 * step
        return new_coord
    
    def _apply_constraints(self, coord):
        mask = ~(np.isnan(self._constraints).any(axis=-1))
        coord[mask] = self._constraints[mask]

    def _potential_function(self, coord):
        vis_dist = np.sqrt(
            np.sum(
                (coord[:, np.newaxis, :] - coord[np.newaxis, :, :])**2, axis=-1
            )
        )
        mean_vis_dist = np.mean(vis_dist)
        # This factor translates visual distances
        # into substitution matrix distances
        scale_factor = self._mean_dist_opt / mean_vis_dist
        # Harmonic potentials between each pair of symbols
        pot = (vis_dist*scale_factor - self._dist_opt)**2
        # Contrast term: Favours conformations
        # that take a large area of the color space
        pot += self._contrast * scale_factor
        return(np.sum(pot))
    
    @staticmethod
    def calculate_distace_matrix(similarity_matrix):
        scores = similarity_matrix.score_matrix()
        distances = np.max(scores, axis=0) - scores
        return distances