import numpy


class InputGrid3D:

    def __init__(
            self,
            x_mesh_flat: numpy.array,
            y_mesh_flat: numpy.array,
            z_mesh_flat: numpy.array,
    ):
        assert x_mesh_flat.shape == y_mesh_flat.shape == z_mesh_flat.shape
        self.x_mesh_flat = x_mesh_flat
        self.x_min = x_mesh_flat.min()
        self.x_max = x_mesh_flat.max()
        self.x_shift = self.x_min
        self.x_scale = 1.0 / ( self.x_max - self.x_min)
        self.y_mesh_flat = y_mesh_flat
        self.y_min = y_mesh_flat.min()
        self.y_max = y_mesh_flat.max()
        self.y_shift = self.y_min
        self.y_scale = 1.0 / ( self.y_max - self.y_min)
        self.z_mesh_flat = z_mesh_flat
        self.z_min = z_mesh_flat.min()
        self.z_max = z_mesh_flat.max()
        self.z_shift = self.z_min
        self.z_scale = 1.0 / ( self.z_max - self.z_min)

    def __add__(self, other):
        return InputGrid3D(
            x_mesh_flat=numpy.concatenate((self.x_mesh_flat, other.x_mesh_flat), axis=0),
            y_mesh_flat=numpy.concatenate((self.y_mesh_flat, other.y_mesh_flat), axis=0),
            z_mesh_flat=numpy.concatenate((self.z_mesh_flat, other.z_mesh_flat), axis=0),
        )

    def __repr__(self):
        return self.__class__.__name__ + f'(size={self.x_mesh_flat.shape[0]})'

    def __str__(self):
        return self.__repr__()
