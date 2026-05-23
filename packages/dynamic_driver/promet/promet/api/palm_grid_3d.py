import numpy


class PALMGrid3D:

    def __init__(
            self,
            x_palm: numpy.ndarray,
            y_palm: numpy.ndarray,
            z_palm: numpy.ndarray,
            origin_x: float = 0.0,
            origin_y: float = 0.0,
            origin_z: float = 0.0,
    ):
        self.x = x_palm + origin_x
        self.y = y_palm + origin_y
        self.z = z_palm + origin_z
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.origin_z = origin_z

        self.x_mesh, self.y_mesh, self.z_mesh = numpy.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.x_mesh_flat = self.x_mesh.flatten()
        self.y_mesh_flat = self.y_mesh.flatten()
        self.z_mesh_flat = self.z_mesh.flatten()

    def __repr__(self):
        return '{}({}, {})'.format(
            self.__class__.__name__,
            f'origin=({self.origin_x}, {self.origin_y}, {self.origin_z})',
            f'size=({self.x.shape[0]}, {self.y.shape[0]}, {self.z.shape[0]})',
        )

    def __str__(self):
        return self.__repr__()
