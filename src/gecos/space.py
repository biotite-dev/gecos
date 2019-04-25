from os.path import join, dirname, realpath
import numpy as np
from .colors import convert_lab_to_rgb


SPACE_FILE_NAME = join(dirname(realpath(__file__)), "space.npy")


class ColorSpace():

    def __init__(self, lightness, file_name=None):
        if lightness < 0 or lightness > 100:
            raise ValueError(f"Lightness value of {lightness} is invalid")
        if file_name is None:
            file_name = SPACE_FILE_NAME
        
        self._lightness = lightness
        self._ab_values = np.arange(-128,128)
        self._lab = np.zeros((256,256,3))
        self._lab[:,:,0] = self._lightness
        self._lab[:,:,1] = np.repeat(self._ab_values[:, np.newaxis], 256, axis=-1)
        self._lab[:,:,2] = np.repeat(self._ab_values[np.newaxis, :], 256, axis= 0)

        with open(file_name, "rb") as file:
            self._space = np.load(file)[lightness]

    def remove(self, mask):
        self._space &= ~mask
    
    def get_rgb_matrix(self):
        rgb = np.full((256,256,3), np.nan)
        for i in range(self._lab.shape[0]):
            for j in range(self._lab.shape[1]):
                if self._space[i,j]:
                    rgb[i,j] = convert_lab_to_rgb(self._lab[i,j])
        return rgb

    @property
    def space(self):
        return self._space.copy()

    @property
    def shape(self):
        return (256, 256)
    
    @property
    def lab(self):
        return self._lab.copy()
    
    @property
    def lightness(self):
        return self._lightness

    @staticmethod
    def _generate(file_name=None):
        lab = np.zeros((100, 256, 256, 3), dtype=int)
        lab[:,:,:,0] = np.arange(100      )[:, np.newaxis, np.newaxis]
        lab[:,:,:,1] = np.arange(-128, 128)[np.newaxis, :, np.newaxis]
        lab[:,:,:,2] = np.arange(-128, 128)[np.newaxis, np.newaxis, :]
        
        rgb = convert_lab_to_rgb(lab)
        space = np.ones((100, 256, 256), dtype=bool)
        space[np.isnan(rgb).any(axis=-1)] = False
        
        if file_name is None:
            file_name = SPACE_FILE_NAME
        with open(file_name, "wb") as file:
            np.save(file, space)