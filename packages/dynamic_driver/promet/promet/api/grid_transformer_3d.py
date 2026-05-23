import numpy
from scipy.interpolate import griddata as gd

from promet.api.palm_grid_3d import PALMGrid3D
from promet.api.input_grid_3d import InputGrid3D


class GridTransformer3D:

    @staticmethod
    def transform(
            input_grid: InputGrid3D,
            input_array: numpy.array,
            output_grid: PALMGrid3D,
            normalize: bool = True,
    ) -> numpy.array:
        if normalize:
            out_array = gd(
                points=(
                    (input_grid.x_mesh_flat - input_grid.x_shift) * input_grid.x_scale,
                    (input_grid.y_mesh_flat - input_grid.y_shift) * input_grid.y_scale,
                    (input_grid.z_mesh_flat - input_grid.z_shift) * input_grid.z_scale,
                ),
                values=input_array,
                xi=(
                    (output_grid.x_mesh - input_grid.x_shift) * input_grid.x_scale,
                    (output_grid.y_mesh - input_grid.y_shift) * input_grid.y_scale,
                    (output_grid.z_mesh - input_grid.z_shift) * input_grid.z_scale,
                ),
                method='linear',
                fill_value=-9999.0,
            )
        else:
            out_array = gd(
                points=(input_grid.x_mesh_flat, input_grid.y_mesh_flat, input_grid.z_mesh_flat),
                values=input_array,
                xi=(output_grid.x_mesh, output_grid.y_mesh, output_grid.z_mesh),
                method='linear',
                fill_value=-9999.0,
            )
        arr_reversed = numpy.transpose(out_array, (2, 1, 0))
        return arr_reversed
