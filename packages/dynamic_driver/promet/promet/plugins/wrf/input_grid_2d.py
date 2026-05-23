import numpy

from promet.api.input_grid_2d import InputGrid2D


class WRFInputGrid2D(InputGrid2D):

    @classmethod
    def centered(cls, x2d, y2d):
        x_flat = x2d.flatten(order='F')
        y_flat = y2d.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat,
            y_mesh_flat=y_flat,
        )

    @classmethod
    def staggered_x(cls, x2d, y2d, dx):
        # create flattened coordinates on u grid
        x2d_u = numpy.zeros(
            shape=(x2d.shape[0] + 1, x2d.shape[1]),
            dtype=x2d.dtype,
        )
        x2d_u[:x2d.shape[0], :] = x2d - dx * 0.5
        x2d_u[-1, :] = x2d[-1, :] + dx * 0.5
        x_flat_u = x2d_u.flatten(order='F')

        y2d_u = numpy.zeros(
            shape=(y2d.shape[0] + 1, y2d.shape[1]),
            dtype=y2d.dtype,
        )
        y2d_u[:y2d.shape[0], :] = y2d
        y2d_u[-1, :] = y2d[-1, :]
        y_flat_u = y2d_u.flatten(order='F')

        return cls(
            x_mesh_flat=x_flat_u,
            y_mesh_flat=y_flat_u,
        )

    @classmethod
    def staggered_y(cls, x2d, y2d, dy):
        # create flattened coordinates on v grid
        x2d_v = numpy.zeros(
            shape=(x2d.shape[0], x2d.shape[1] + 1),
            dtype=x2d.dtype,
        )
        x2d_v[:, :x2d.shape[1]] = x2d
        x2d_v[:, -1] = x2d[:, -1]
        x_flat_v = x2d_v.flatten(order='F')

        y2d_v = numpy.zeros(
            shape=(y2d.shape[0], y2d.shape[1] + 1),
            dtype=y2d.dtype,
        )
        y2d_v[:, :y2d.shape[1]] = y2d - dy * 0.5
        y2d_v[:, -1] = y2d[:, -1] + dy * 0.5
        y_flat_v = y2d_v.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat_v,
            y_mesh_flat=y_flat_v,
        )
