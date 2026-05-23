import numpy


class PALMGrid2D:

    def __init__(
            self,
            x_palm: numpy.ndarray,
            y_palm: numpy.ndarray,
            origin_x: float = 0.0,
            origin_y: float = 0.0,
    ):
        self.x = x_palm + origin_x
        self.y = y_palm + origin_y
        self.origin_x = origin_x
        self.origin_y = origin_y

        self.x_mesh, self.y_mesh = numpy.meshgrid(self.x, self.y, indexing='ij')
        self.x_mesh_flat = self.x_mesh.flatten()
        self.y_mesh_flat = self.y_mesh.flatten()
