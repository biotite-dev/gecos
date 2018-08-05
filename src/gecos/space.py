import numpy as np
from .colors import convert_lab_to_rgb


class ColorSpace():

    def __init__(self, lightness):
        if lightness < 0 or lightness > 100:
            raise ValueError(f"Lightness value of {lightness} is invalid")
        self._lightness =lightness
        self._ab_values = np.arange(-128,128)
        self._space = np.ones((256,256), dtype=bool)
        self._lab = np.zeros((256,256,3))
        self._lab[:,:,0] = self._lightness
        self._lab[:,:,1] = np.repeat(self._ab_values[:, np.newaxis], 256, axis=-1)
        self._lab[:,:,2] = np.repeat(self._ab_values[np.newaxis, :], 256, axis= 0)
        rgb = convert_lab_to_rgb(self._lab)
        self._space[np.isnan(rgb).any(axis=-1)] = False

    def remove(self, mask):
        self._space &= mask
    
    def get_rgb_matrix(self):
        rgb = np.full((256,256,3), np.nan)
        for i in range(self._lab.shape[0]):
            for j in range(self._lab.shape[1]):
                if self._space[i,j]:
                    rgb[i,j] = convert_lab_to_rgb(self._lab[i,j])
        return rgb

    @property
    def shape(self):
        return (256, 256)
    
    @property
    def l(self):
        return self._lightness
    
    @property
    def a(self):
        return self._ab_values.copy()
    
    @property
    def b(self):
        return self._ab_values.copy()


        