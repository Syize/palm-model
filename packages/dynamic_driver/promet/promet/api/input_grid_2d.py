import numpy


class InputGrid2D:

    def __init__(
            self,
            x_mesh_flat: numpy.array,
            y_mesh_flat: numpy.array,
    ):
        assert x_mesh_flat.shape == y_mesh_flat.shape
        self.x_mesh_flat = x_mesh_flat
        self.y_mesh_flat = y_mesh_flat

    def __repr__(self):
        return self.__class__.__name__ + f'(size={self.x_mesh_flat.shape[0]})'

    def __str__(self):
        return self.__repr__()
