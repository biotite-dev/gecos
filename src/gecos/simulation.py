import numpy as np
import numpy.random as random
import biotite.sequence as seq
import biotite.sequence.align as align


class Result():
    pass


def generate_color_scheme(matrix, space, constraints=None, n_steps=60000,
                          t_initial=10, t_final=0.01, ext_factor=1, seed=None):
    
    def is_allowed(coord):
        nonlocal space
        return space[int(coord[0])+128, int(coord[1])+128]
    
    def _move(coord):
        nonlocal space
        factor = 1.0
        new_coord = coord + (random.rand(*coord.shape)-0.5) * factor
        # Resample coordinates for alphabet symbols
        # when outside of the allowed area
        for i in range(new_coord.shape[0]):
            while not is_allowed(new_coord[i]):
                new_coord[i] = coord[i] + (random.rand(2)-0.5) * factor
        return new_coord


    def _potential_function(coord):
        nonlocal dist_opt
        nonlocal ext_factor
        dist = np.sqrt(
            np.sum(
                (coord[:, np.newaxis, :] - coord[np.newaxis, :, :])**2, axis=-1
            )
        )
        pot = (dist-dist_opt*ext_factor)**2
        return(np.sum(pot))
    

    result = Result()
    
    if seed is None:
        seed = random.randint(np.iinfo(np.int32).max, dtype=np.int32)
    random.seed(seed=seed)
    result.seed = seed

    space = space.space

    if matrix.get_alphabet1() != matrix.get_alphabet2():
        raise ValueError("The substiution matrix has unequal alphabets")
    alphabet = matrix.get_alphabet1()
    scores = matrix.score_matrix()
    dist_opt = np.max(scores, axis=0) - scores
    #mean_dist_opt = np.mean(dist_opt)

    temps = np.logspace(np.log10(t_initial), np.log10(t_final), n_steps)
    potentials = np.zeros(n_steps)
    coord = np.zeros((len(alphabet), 2), dtype=float)
    trajectory = np.zeros((n_steps, len(alphabet), 2), dtype=int)

    pot = _potential_function(coord)
    
    for i in range(n_steps):
        temp = temps[i]
        new_coord = _move(coord)
        new_pot = _potential_function(new_coord)
        if new_pot < pot:
            coord = new_coord
            pot = new_pot
        else:
            p = np.exp(-(new_pot-pot)/temp)
            if p > random.rand():
                coord = new_coord
                pot = new_pot
        trajectory[i] = coord.astype(int)
        potentials[i] = pot
    
    result.trajectory = trajectory
    result.potentials = potentials
    return result