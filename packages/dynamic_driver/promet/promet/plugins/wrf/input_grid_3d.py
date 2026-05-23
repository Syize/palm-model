import numpy

from promet.api.input_grid_3d import InputGrid3D


class WRFInputGrid3D(InputGrid3D):

    @classmethod
    def centered(cls, x2d, y2d, z3d, terrain=None):
        if terrain is None:
            z3d_local = z3d
        else:
            z3d_local = numpy.concatenate((terrain[None, ...], z3d), axis=0)
        x2d_expanded = numpy.expand_dims(x2d, axis=0)
        x2d_tileed = numpy.tile(x2d_expanded, (z3d_local.shape[0], 1, 1))
        x_flat = x2d_tileed.flatten(order='F')
        y2d_expanded = numpy.expand_dims(y2d, axis=0)
        y2d_tileed = numpy.tile(y2d_expanded, (z3d_local.shape[0], 1, 1))
        y_flat = y2d_tileed.flatten(order='F')
        z_flat = z3d_local.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat,
            y_mesh_flat=y_flat,
            z_mesh_flat=z_flat,
        )

    @classmethod
    def centered_2d(cls, x2d, y2d, z2d):
        x_flat = x2d.flatten(order='F')
        y_flat = y2d.flatten(order='F')
        z_flat = z2d.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat,
            y_mesh_flat=y_flat,
            z_mesh_flat=z_flat,
        )

    @classmethod
    def staggered_x(cls, x2d, y2d, z3d, dx, terrain=None):
        if terrain is None:
            z3d_local = z3d
        else:
            z3d_local = numpy.concatenate((terrain[None, ...], z3d), axis=0)
        # create flattened coordinates on u grid
        ext_z = numpy.full(shape=(z3d_local.shape[0],), fill_value=0)
        x2d_u = numpy.zeros(
            shape=(x2d.shape[0], x2d.shape[1] + 1),
            dtype=x2d.dtype,
        )
        x2d_u[:, :x2d.shape[1]] = x2d - dx * 0.5
        x2d_u[:, -1] = x2d[:, -1] + dx * 0.5
        x_flat_u = (x2d_u[..., None] + ext_z).flatten(order='C')

        y2d_u = numpy.zeros(
            shape=(y2d.shape[0], y2d.shape[1] + 1),
            dtype=y2d.dtype,
        )
        y2d_u[:, :x2d.shape[1]] = y2d
        y2d_u[:, -1] = y2d[:, -1]
        y_flat_u = (y2d_u[..., None] + ext_z).flatten(order='C')

        z3d_u = numpy.zeros(
            shape=(z3d_local.shape[0], z3d_local.shape[1], z3d_local.shape[2] + 1),
            dtype=z3d_local.dtype,
        )
        z3d_u[:, :, :z3d_local.shape[2]] = z3d_local
        z3d_u[:, :, -1] = z3d_local[:, :, -1]
        z_flat_u = z3d_u.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat_u,
            y_mesh_flat=y_flat_u,
            z_mesh_flat=z_flat_u,
        )

    @classmethod
    def staggered_y(cls, x2d, y2d, z3d, dy, terrain=None):
        if terrain is None:
            z3d_local = z3d
        else:
            z3d_local = numpy.concatenate((terrain[None, ...], z3d), axis=0)
        # create flattened coordinates on v grid
        ext_z = numpy.full(shape=(z3d_local.shape[0],), fill_value=0)
        x2d_v = numpy.zeros(
            shape=(x2d.shape[0] + 1, x2d.shape[1]),
            dtype=x2d.dtype,
        )
        x2d_v[:x2d.shape[0], :] = x2d
        x2d_v[-1, :] = x2d[-1, :]
        x_flat_v = (x2d_v[..., None] + ext_z).flatten(order='C')

        y2d_v = numpy.zeros(
            shape=(y2d.shape[0] + 1, y2d.shape[1]),
            dtype=y2d.dtype,
        )
        y2d_v[:x2d.shape[0], :] = y2d - dy * 0.5
        y2d_v[-1, :] = y2d[-1, :] + dy * 0.5
        y_flat_v = (y2d_v[..., None] + ext_z).flatten(order='C')

        z3d_v = numpy.zeros(
            shape=(z3d_local.shape[0], z3d_local.shape[1] + 1, z3d_local.shape[2]),
            dtype=z3d_local.dtype,
        )
        z3d_v[:, :z3d_local.shape[1], :] = z3d_local
        z3d_v[:, -1, :] = z3d_local[:, -1, :]
        z_flat_v = z3d_v.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat_v,
            y_mesh_flat=y_flat_v,
            z_mesh_flat=z_flat_v,
        )

    @classmethod
    def staggered_z(cls, x2d, y2d, zw3d):
        # create flattened coordinates on w grid
        ext_zw = numpy.full(shape=(zw3d.shape[0],), fill_value=0)
        x_flat_w = (x2d[..., None] + ext_zw).flatten(order='C')
        y_flat_w = (y2d[..., None] + ext_zw).flatten(order='C')
        z_flat_w = zw3d.flatten(order='F')
        return cls(
            x_mesh_flat=x_flat_w,
            y_mesh_flat=y_flat_w,
            z_mesh_flat=z_flat_w,
        )


    @classmethod
    def centered_zsoil(cls, x2d, y2d, zsoil1d):
        ext_z = numpy.full(shape=(zsoil1d.shape[0],), fill_value=0)
        x_flat = (x2d[..., None] + ext_z).flatten(order='C')
        y_flat = (y2d[..., None] + ext_z).flatten(order='C')
        ext_h = numpy.full(shape=(x2d.shape[1], x2d.shape[0]), fill_value=0)
        zsoil_flat = (zsoil1d[..., None, None] + ext_h).flatten(order='F')
        return cls(
            x_mesh_flat=x_flat,
            y_mesh_flat=y_flat,
            z_mesh_flat=zsoil_flat,
        )