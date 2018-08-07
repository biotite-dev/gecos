import numpy as np
import numpy.random as random
import biotite.sequence as seq
import biotite.sequence.align as align
from .colors import convert_lab_to_rgb


class Result():
    
    @property
    def lightness(self):
        return self._lightness
    
    @lightness.setter
    def lightness(self, val):
        self._lightness = val

    @property
    def alphabet(self):
        return self._alphabet
    
    @alphabet.setter
    def alphabet(self, val):
        self._alphabet= val

    @property
    def seed(self):
        return self._seed
    
    @seed.setter
    def seed(self, val):
        self._seed = val
    
    @property
    def potentials(self):
        return self._potentials
    
    @potentials.setter
    def potentials(self, val):
        self._potentials = val
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, val):
        self._positions = val
    
    @property
    def final_potential(self):
        return self._potentials[-1]
    
    @property
    def final_position(self):
        return self._positions[-1]
    
    @property
    def rgb_colors(self):
        ab = self.final_position
        lab = np.stack(
            (np.full(len(ab), self.lightness), ab[:,0], ab[:,1]), axis=-1
        )
        return convert_lab_to_rgb(lab)



def generate_color_scheme(matrix, space, constraints=None, n_steps=100000,
                          temp=(100, 0.1), step_size=(10,0.1),
                          ext_factor=1, seed=None):
    
    def is_allowed(coord):
        nonlocal space
        return space[int(coord[0])+128, int(coord[1])+128]
    
    def _move(coord, step):
        nonlocal space
        new_coord = coord + (random.rand(*coord.shape)-0.5) * 2 * step
        # Resample coordinates for alphabet symbols
        # when outside of the allowed area
        for i in range(new_coord.shape[0]):
            while not is_allowed(new_coord[i]):
                new_coord[i] = coord[i] + (random.rand(2)-0.5) * step
        return new_coord


    def _potential_function(coord):
        nonlocal dist_opt
        nonlocal mean_dist_opt
        nonlocal ext_factor
        vis_dist = np.sqrt(
            np.sum(
                (coord[:, np.newaxis, :] - coord[np.newaxis, :, :])**2, axis=-1
            )
        )
        mean_vis_dist = np.mean(vis_dist)
        scale_factor = mean_dist_opt / mean_vis_dist
        # Harmonic potential terms
        pot = (vis_dist*scale_factor - dist_opt)**2
        # extension term
        pot += ext_factor * scale_factor
        return(np.sum(pot))
    

    result = Result()

    result.lightness = space.l
    
    if seed is None:
        seed = random.randint(np.iinfo(np.int32).max, dtype=np.int32)
    random.seed(seed=seed)
    result.seed = seed

    a = space.a
    b = space.b
    space = space.space

    if matrix.get_alphabet1() != matrix.get_alphabet2():
        raise ValueError("The substiution matrix has unequal alphabets")
    alphabet = matrix.get_alphabet1()
    result.alphabet = alphabet
    scores = matrix.score_matrix()
    dist_opt = np.max(scores, axis=0) - scores
    mean_dist_opt = np.mean(dist_opt)

    temps = np.logspace(np.log10(temp[0]), np.log10(temp[1]), n_steps)
    steps = np.logspace(np.log10(
        step_size[0]), np.log10(step_size[1]), n_steps
    )

    potentials = np.zeros(n_steps)

    coord = np.zeros((len(alphabet), 2), dtype=float)
    # Find a suitable start position
    indices = np.nonzero(space)
    start_coord = (a[indices[0][0]], b[indices[1][0]])
    coord[:,:] = start_coord
    # Initial move to avoid error in _potential_function():
    # Without move, mean_vis_dist = 0 -> scale_factor = infinite
    coord = _move(coord, 1)
    trajectory = np.zeros((n_steps, len(alphabet), 2), dtype=int)

    pot = _potential_function(coord)
    
    for i in range(n_steps):
        temp = temps[i]
        new_coord = _move(coord, steps[i])
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
    
    result.positions = trajectory
    result.potentials = potentials
    return result