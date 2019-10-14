# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann, Benjamin Mayer"
__all__ = ["ColorOptimizer", "ScoreFunction", "DefaultScoreFunction"]

from collections import namedtuple
import abc
import copy
import skimage
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
    """
    Create an optimizer that tries to find an optimal color conformation
    within a given color space based on a score function.

    The optimizer tries to minimize the return value of the score
    function by adjusting the *Lab* values (coordinates) for each
    symbol in a given alphabet.

    The optimizer uses the random number generator from *NumPy*.
    Therefore, call :func:`numpy.random.seed()` to set the seed for the
    optimizer

    Parameters
    ----------
    alphabet : biotite.sequence.Alphabet
        The alphabet to calculate the color conformation for.
    score_function : ScoreFunction or callable
        The score function which should be minimized.
        When calling the object, its only parameter must be an array of
        coordinates with shape *(n, 3)*, where *n* is the length of the
        alphabet.
        Its return value must be a single float - the score.
    space : ColorSpace
        The color space that defines the allowed space for the
        coordinates.
    constraints : ndarray, shape=(n,3), dtype=float, optional
        An array whose non-NaN values are interpreted as constraints.
        Constrained values will be fixed during the optimization.
    """

    class Result(namedtuple("Result", ["alphabet", "trajectory", "scores"])):
        """
        The result of an optimization.
        Contains the final color scheme information as well as the
        course of the coordinates and the score during the optimization.

        Parameters
        ----------
        alphabet : biotite.sequence.Alphabet
            The alphabet the optimizer used.
        trajectory : ndarray, shape=(m,n,3), dtype=float
            The course of the coordinates during the simulation.
        scores : ndarray, shape=(m,), dtype=float
            The course of the score during the simulation.

        Attributes
        ----------
        alphabet : biotite.sequence.Alphabet
            The alphabet the optimizer used.
        trajectory : ndarray, shape=(m,n,3), dtype=float
            The course of coordinates during the simulation.
        lab_colors : ndarray, shape=(n,3), dtype=float
            The final *Lab* color conformation, i.e. the last element of
            `trajectory`.
        rgb_colors : ndarray, shape=(n,3), dtype=float
            The final color conformation converted into *RGB* colors.
        scores : ndarray, shape=(m,), dtype=float
            The course of the score during the simulation.
        score : float
            The final score, i.e. the last element of `scores`.
        
        """
        
        @property
        def score(self):
            return self.scores[-1]

        @property
        def lab_colors(self):
            return copy.deepcopy(self.trajectory[-1])
        
        @property
        def rgb_colors(self):
            return lab_to_rgb(self.lab_colors.astype(int))
    
    def __init__(self, alphabet, score_function, space, constraints=None):
        self._alphabet = alphabet
        self._n_symbols = len(alphabet)
        self._score_func = score_function
        self._space = space.space.copy()
        self._coord = None
        self._trajectory = []
        self._scores = []

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
        """
        Set the the coordinates of the current color conformation.
        Potential color constraints are applied on these.
        This coordinate changes will be tracked in the trajectory.
        
        Parameters
        ----------
        coord : ndarray, shape=(n,3), dtype=float
            The new coordinates.
        """
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
        
    def _set_coordinates(self, coord, score=None):
        self._coord = coord
        self._trajectory.append(coord)
        if score is None:
            score = self._score_func(coord)
        self._scores.append(score)
    
    def optimize(self, n_steps,
                 beta_start, rate_beta, stepsize_start, stepsize_end):
        r"""
        Perform a Simulated Annealing optimization on the current
        coordinate to minimize the score returned by the score function.
        
        This is basically a Monte-Carlo optimization where the
        temperature is varied according to a so called annealing
        schedule over the course of the optimization. 
        The algorithm is a heuristic thats motivated by the physical
        process of annealing.
        If we, e.g., cool steel than a slow cooling can yield a superior
        quality, whereas for a fast cooling the steel can become
        brittle.
        The same happens here within the search space for the given
        minimization task.                              
        
        Parameters
        ----------
        n_steps : int
            The number of Simulated-Annealing steps.
        beta_start : float
            The inverse start temperature, where the start temperature
            would be :math:`T_{start} = 1/(k_b \cdot \beta_{start})` with
            :math:`k_b` being the boltzmann constant.            
        rate: float
            The rate controlls how fast the inverse temperature is
            increased within the annealing schedule.
            Here the exponential schedule is chosen so we have
            :math:`\beta (t) = \beta_0 \cdot \exp(rate \cdot t)`.
        stepsize_start : float
            The radius in which the coordinates are randomly altered at
            the beginning of the simulated anneling algorithm.
            Like the inverse temperature the step size follows an
            exponential schedule, enabling the algorithm
            to do large perturbartions at the beginning of the algorithm
            run and increasingly smaller ones afterwards.
        stepsize_end : float
            The radius in which the coordinates are randomly altered at
            the end of the simulated annealing algorithm run.          
        """

        rate_stepsize = None
        beta = lambda i: beta_start*np.exp(rate_beta*i)

        #  Choose rate so that stepsize_end reached after n_steps
        #  derived from step_size(N_steps) = steps_end
        if stepsize_start == stepsize_end:
            rate_stepsize = 0
        else:
            if stepsize_end is None:
                rate_stepsize = -1
            else:            
                rate_stepsize = np.log(stepsize_end / stepsize_start) / n_steps
        step_size = lambda i: stepsize_start * np.exp(rate_stepsize * i)

        for i in range(n_steps):
        
            score = self._scores[-1]
            new_coord = self._move(self._coord, step_size(i))
            new_score = self._score_func(new_coord)
            
            if new_score < score:
                self._set_coordinates(new_coord, new_score)
                
            else:
                p_accept = np.exp( -beta(i) * (new_score-score))
                p = random.rand()
                
                if p <= p_accept:
                    self._set_coordinates(new_coord, new_score)
                else:
                    self._set_coordinates(self._coord, new_score)
                    
    def get_result(self):
        """
        Get the result of the optimization.

        Returns
        -------
        result : ColorOptimizer.Result
            The result.
        """
        trajectory = np.array(self._trajectory)
        self._scores = np.array(self._scores)
        return ColorOptimizer.Result(
            alphabet = self._alphabet,
            trajectory = trajectory,
            scores = self._scores
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


class ScoreFunction(metaclass=abc.ABCMeta):
    """
    Abstract base class for a score function.
    A score function calculates a score from a color conformation
    (coordinates).

    The score is calculated by calling the object with the coordinates
    as single argument.
    Hence, classes inheriting from this base class mut override the
    :func:`__call__()` method.

    Parameters
    ----------
    n_symbols : int
        The amount of symbols in the system.
        Equivalent to the length of the alphabet the color scheme is
        generated for.
        This value is used to check the shape of the coordinates when
        calling the score function.
    """

    def __init__(self, n_symbols):
        self._n_symbols = n_symbols

    @abc.abstractmethod
    def __call__(self, coord):
        """
        Calculate the score for the given coordinates.

        Parameters
        ----------
        coord : ndarray, shape=(n,3), dtype=float
            The coordinates.
        
        Returns
        -------
        score : float
            The score assigned to `coord`.
        """
        if len(coord) != self._n_symbols:
            raise ValueError(
                f"Expected {self._n_symbols} coordinates, but got {len(coord)}"
            )


class DefaultScoreFunction(ScoreFunction):
    """
    Create an instance of the default score function *Gecos* uses.

    The score function contains two terms:
    A sum of harmonic potentials between each pair of symbols, based on
    a substitution matrix, and *contrast score* that favors schemes with
    a high contrast.

    Parameters
    ----------
    matrix : biotite.sequence.align.SubstitutionMatrix
        A distance matrix is calculated from this score matrix.
        The equilibrium positions scale linearly with the values in the
        distance matrix.
    contrast : int, optional
        A weight for the *contrast score*.
    distance_formula : {'CIE76', 'CIEDE94', 'CIEDE2000'}, optional
        The formula to use for calculation of perceptual color
        difference.
        While ``'CIEDE2000'`` is the most accurate formula for the
        perceptual difference, ``'CIE76'`` features the fastest
        calculation.
    """

    def __init__(self, matrix, contrast=700, distance_formula="CIEDE2000"):
        if not matrix.is_symmetric():
            raise ValueError("Substitution matrix must be symmetric")
        super().__init__(len(matrix.get_alphabet1()))
        self._matrix = self._calculate_distance_matrix(matrix)
        self._n = DefaultScoreFunction._n_pairs(len(matrix.score_matrix()))
        self._contrast = contrast
        if distance_formula not in ["CIE76", "CIEDE94", "CIEDE2000"]:
            raise ValueError(
                f"Unknown color distance measure f'{distance_measure}'"
            )
        self._distance_formula = distance_formula
    
    def __call__(self, coord):
        super().__call__(coord)
        dist = DefaultScoreFunction._calculate_distance(
            coord, self._distance_formula
        )
        dist_sum = np.sum(dist)
        # This factor translates visual distances
        # into substitution matrix distances
        scale_factor = self._n / dist_sum
        # Harmonic potentials between each pair of symbols
        harmonic_score = np.sum((dist*scale_factor - self._matrix)**2)
        # Contrast term: Favours conformations
        # with large absolute color differences
        mean_dist = dist_sum / DefaultScoreFunction._n_pairs(len(dist))
        contrast_score = self._contrast / mean_dist
        return harmonic_score + contrast_score
    
    @staticmethod
    def _calculate_distance(coord, distance_formula):
        ind1, ind2 = np.tril_indices(len(coord), k=-1)
        flat_coord1 = coord[ind1]
        flat_coord2 = coord[ind2]
        if distance_formula == "CIEDE76":
            flat_dist = skimage.color.deltaE_ciede94(
                flat_coord1, flat_coord2
            )
        elif distance_formula == "CIEDE94":
            flat_dist = skimage.color.deltaE_cie76(
                flat_coord1, flat_coord2
            )
        else: #"CIEDE2000"
            flat_dist = skimage.color.deltaE_ciede2000(
                flat_coord1, flat_coord2
            )
        dist = np.zeros((len(coord),)*2)
        dist[ind1, ind2] = flat_dist
        return dist

    @staticmethod
    def _calculate_distance_matrix(similarity_matrix):
        scores = similarity_matrix.score_matrix()
        diff_to_max = np.diag(scores) - scores
        distances = np.tril((diff_to_max + diff_to_max.T) / 2)
        # Scale, so that average distance is 1
        n = DefaultScoreFunction._n_pairs(len(scores))
        distances /= (np.sum(distances) / n)
        return distances
    
    @staticmethod
    def _n_pairs(n_symbols):
        """
        Calculate the number of values in the lower triangle,
        excluding the main diagonal, of a
        matrix with a shape *(n_symbols, n_symbols)*.
        """
        return (n_symbols - 1) / 2 * n_symbols