import numpy

from promet.api.input_grid_3d import InputGrid3D


class ICONInputGrid3D(InputGrid3D):

    @classmethod
    def centered(cls, x_cells, y_cells, z_cells):
        return cls(
            x_mesh_flat=x_cells,
            y_mesh_flat=y_cells,
            z_mesh_flat=z_cells,
        )

    @classmethod
    def centered_zsoil(cls, x_cells, y_cells, zsoil1d):
        ext_z = numpy.full(shape=(zsoil1d.shape[0],), fill_value=0)
        x_flat = (x_cells[..., None] + ext_z).flatten(order='C')
        y_flat = (y_cells[..., None] + ext_z).flatten(order='C')
        ext_h = numpy.full(shape=(x_cells.shape[0],), fill_value=0)
        zsoil_flat = (zsoil1d[..., None] + ext_h).flatten(order='F')
        return cls(
            x_mesh_flat=x_flat,
            y_mesh_flat=y_flat,
            z_mesh_flat=zsoil_flat,
        )
