# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import pytest
import numpy as np
from biotite.sequence import ProteinSequence, Alphabet
from biotite.sequence.align import SubstitutionMatrix
import gecos


MIN_COORD = np.array([ 0, -128, -128], dtype=float)
MAX_COORD = np.array([99,  127,  127], dtype=float)


def test_scale_factor():
    """
    The score for a random conformation of two symbols should be equal
    for all conformations at a contrast factor of 0, due to the scale
    factor.
    Due to scaling, the distance should always be the
    equilibrium distance of the harmonic potential -> score = 0.
    """
    N_RANDOM = 1000
    np.random.seed(0)

    # Create alphabet with two symbols
    # and a identity substution matrix for it
    alph = Alphabet(["A", "B"])
    matrix = SubstitutionMatrix(alph, alph, np.identity(len(alph)))
    # Important: contrast factor must be 0
    score_func = gecos.DefaultScoreFunction(matrix, contrast=0)

    space = gecos.ColorSpace()
    scores = []
    for _ in range(N_RANDOM):
        # Get score for random conformation
        random_coord = _draw_random(len(alph), space)
        scores.append(score_func(random_coord))
    
    assert scores == pytest.approx([0] * len(scores))


def test_optimization_score():
    """
    Assert that *n* randomly drawn conformations score worse than an
    optimization with *m* steps.
    """
    N_RANDOM = 10000
    N_STEPS = 1000
    np.random.seed(0)
    
    matrix = SubstitutionMatrix.std_protein_matrix()
    alphabet = ProteinSequence.alphabet
    space = gecos.ColorSpace()
    score_func = gecos.DefaultScoreFunction(matrix)
    
    optimizer = gecos.ColorOptimizer(alphabet, score_func, space)
    optimizer.optimize(N_STEPS)
    optimized_score = optimizer.get_result().score

    scores = []
    for _ in range(N_RANDOM):
        # Get score for random conformation
        random_coord = _draw_random(len(alphabet), space)
        scores.append(score_func(random_coord))
    
    assert (np.array(scores) > optimized_score).all()


def test_optimized_distances():
    """
    Assert that three symbols have an almost optimal conformation after
    the optimization.
    Three symbols can always be arranged optimally, as long as no
    distance is larger than the combined other two distances.
    """
    # Create alphabet with three symbols
    # and a toy substitution matrix for it
    alph = Alphabet(["A", "B", "C"])
    score_matrix = [
        [10,  5,  3],
        [ 5, 10,  2],
        [ 3,  2, 10]
    ]

    matrix = SubstitutionMatrix(alph, alph, np.array(score_matrix))
    # Contrast factor is 0 to optimize only for pairwise distances
    score_func = gecos.DefaultScoreFunction(
        matrix, contrast=0, distance_formula="CIE76"
    )
    ideal_distances = score_func._ideal_dist
    a_to_b_ref, a_to_c_ref, b_to_c_ref = ideal_distances

    np.random.seed(0)
    space = gecos.ColorSpace()
    optimizer = gecos.ColorOptimizer(alph, score_func, space)
    optimizer.optimize()
    result = optimizer.get_result()
    optimized_coord = result.lab_colors
    optimized_score = result.score
    start_score = result.scores[0]

    # Expect 0, since all three symbols can be arranged optimally,
    # i.e. all distances can be the equilibrium distance
    assert optimized_score < 0.01 * start_score
    
    # In an optimal conformation the normed distance
    # between each pair of symbols is equal to the respective distance
    # in the distance matrix
    a_to_b_test = _distance(optimized_coord[0], optimized_coord[1])
    a_to_c_test = _distance(optimized_coord[0], optimized_coord[2])
    b_to_c_test = _distance(optimized_coord[1], optimized_coord[2])
    mean_dist = np.mean([a_to_b_test, a_to_c_test, b_to_c_test])
    a_to_b_test /= mean_dist
    a_to_c_test /= mean_dist
    b_to_c_test /= mean_dist
    assert a_to_b_test == pytest.approx(a_to_b_ref, rel=0.1)
    assert a_to_c_test == pytest.approx(a_to_c_ref, rel=0.1)
    assert b_to_c_test == pytest.approx(b_to_c_ref, rel=0.1)


@pytest.mark.parametrize("constraint_seed", np.arange(10))
def test_constraints(constraint_seed):
    N_STEPS = 100
    alph = Alphabet([s for s in ProteinSequence.alphabet.get_symbols()[:20]])
    matrix = SubstitutionMatrix(alph, alph, "BLOSUM62")

    constraint_mask = np.random.choice([False, True], len(alph))
    # Arbitrarily set constraint colors to gray (50, 0, 0)
    constraint_pos = np.repeat(
        np.array([[50, 0, 0]], dtype=float),
        len(alph),
        axis=0
    )
    constraint_pos[~constraint_mask] = np.nan

    np.random.seed(constraint_seed)
    space = gecos.ColorSpace()
    score_func = gecos.DefaultScoreFunction(matrix)
    optimizer = gecos.ColorOptimizer(alph, score_func, space, constraint_pos)
    optimizer.optimize(N_STEPS)
    result = optimizer.get_result()
    opt_coord = result.lab_colors

    # No constrained colors should have changed,
    # but all non-constrained colors during the optimization
    assert opt_coord[constraint_mask].tolist() \
        == constraint_pos[constraint_mask].tolist()
    assert opt_coord[~constraint_mask].tolist() \
        != constraint_pos[~constraint_mask].tolist()


def _draw_random(n_symbols, space):
    space = space.space
    random_coord = np.zeros((n_symbols, 3))
    # Resample coordinates that are not in valid space until they
    # are in valid space
    resample_mask = np.ones(n_symbols, dtype=bool)
    while resample_mask.any():
        random_coord[resample_mask] = (
            np.random.rand(np.count_nonzero(resample_mask), 3)
            * (MAX_COORD - MIN_COORD) + MIN_COORD
        )
        resample_mask[_is_allowed(random_coord, space)] = False
    return random_coord


def _is_allowed(coord, space):
    mask = ((coord >= MIN_COORD) & (coord <= MAX_COORD)).all(axis=-1)
    # Only check values that are within valid index range
    ind = np.floor(coord[mask] - MIN_COORD).astype(int)
    mask[mask.copy()] = space[ind[..., 0], ind[..., 1], ind[..., 2]]
    return mask


def _distance(coord_a, coord_b):
    return np.sqrt(np.sum((coord_a - coord_b)**2))