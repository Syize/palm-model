import numpy
from scipy.interpolate import griddata as gd

from promet.api.palm_grid_2d import PALMGrid2D
from promet.api.input_grid_2d import InputGrid2D


class GridTransformer2D:

    @staticmethod
    def transform(
            input_grid: InputGrid2D,
            input_array: numpy.array,
            output_grid: PALMGrid2D,
            method: str = 'linear',
    ) -> numpy.array:
        out_array = gd(
            points=(input_grid.x_mesh_flat, input_grid.y_mesh_flat),
            values=input_array,
            xi=(output_grid.x_mesh, output_grid.y_mesh),
            method=method,
            fill_value=-9999.0,
        )
        arr_reversed = numpy.transpose(out_array, (1, 0))
        return arr_reversed
