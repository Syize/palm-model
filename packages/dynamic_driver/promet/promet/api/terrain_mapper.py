import numpy
from copy import deepcopy
from pprint import pprint

from promet.api.grid_transformer_2d import GridTransformer2D
from promet.api.input_grid_2d import InputGrid2D
from promet.api.palm_grid_3d import PALMGrid3D
from promet.api.palm_grid_2d import PALMGrid2D


class TerrainMapper:

    @staticmethod
    def apply(
            input_grid: PALMGrid3D,
            input_array: numpy.array,
            input_terrain: numpy.array,
            output_terrain: numpy.array,
            mode: str,
            top_adjustment_index: int = -1,
    ):
        """
        - work on each lateral boundary individually, but pre-compute z shifting factors
        - only work with data up to z of max adjustment height.
        - get mesoskale terrain heigt for each n, (use original interpolator)
        - get palm terrain height for each n (read static driver)
        - compute tereign shift for all n simultaneously.
        """
        if mode in ['2d_left', '2d_right']:
            palm_grid_2d = PALMGrid2D(
                x_palm=input_grid.y - input_grid.origin_y,
                y_palm=input_grid.z - input_grid.origin_z,
                origin_x=input_grid.origin_y,
                origin_y=input_grid.origin_z,
            )
            sliced_array = input_array[:, :, 0]
        elif mode in ['2d_south', '2d_north']:
            palm_grid_2d = PALMGrid2D(
                x_palm=input_grid.x - input_grid.origin_x,
                y_palm=input_grid.z - input_grid.origin_z,
                origin_x=input_grid.origin_x,
                origin_y=input_grid.origin_z,
            )
            sliced_array = input_array[:, 0, :]
        else:
            raise NotImplementedError()

        z_mesh, x_mesh = numpy.meshgrid(palm_grid_2d.y, palm_grid_2d.x, indexing='ij')
        height_at_input_terrain = numpy.tile(input_terrain, (z_mesh.shape[0], 1))
        indices_below_input_terrain_height = numpy.argmax(z_mesh >= height_at_input_terrain, axis=0) - 1

        mask = numpy.arange(sliced_array.shape[0]) > indices_below_input_terrain_height[:, numpy.newaxis]
        mask = numpy.transpose(mask, (1, 0))
        fill_value = -9999.0
        sliced_array = numpy.where(mask, sliced_array, fill_value)

        if mode in ['2d_left', '2d_right']:
            input_array[:, :, 0] = sliced_array
        elif mode in ['2d_south', '2d_north']:
            input_array[:, 0, :] = sliced_array
        else:
            raise NotImplementedError()
        height_at_output_terrain = numpy.tile(output_terrain, (z_mesh.shape[0], 1))

        indices_below_output_terrain_height = numpy.argmax(z_mesh >= height_at_output_terrain, axis=0) - 1
        required_vertical_index_shift = indices_below_output_terrain_height - indices_below_input_terrain_height

        required_vertical_index_shift_2d = numpy.tile(required_vertical_index_shift, (z_mesh.shape[0], 1))

        z_index_1d = numpy.arange(z_mesh.shape[0])
        z_index_2d = numpy.tile(z_index_1d, (z_mesh.shape[1], 1))
        z_index_2d = numpy.transpose(z_index_2d, (1, 0))
        height_at_output_terrain = numpy.tile(indices_below_input_terrain_height, (z_mesh.shape[0], 1))

        if top_adjustment_index < 0:
            top_adjustment_index = z_mesh.shape[0] + top_adjustment_index
        factor_array = 1 - ((z_index_2d - height_at_output_terrain) / (top_adjustment_index - height_at_output_terrain))
        factor_array = numpy.where(factor_array > 1, 1, factor_array)
        factor_array = numpy.where(factor_array < 0, 0, factor_array)
        new_z_mesh = z_mesh * (1-factor_array) + (z_mesh + required_vertical_index_shift_2d*20) * factor_array
        new_z_mesh = numpy.transpose(new_z_mesh, (1, 0))


        palm_grid_2d_shifted = deepcopy(palm_grid_2d)
        palm_grid_2d_shifted.y_mesh_flat = new_z_mesh.flatten()

        sliced_array_flat = sliced_array.flatten(order='F')
        filter = numpy.where(sliced_array_flat != fill_value)

        input_grid_shifted = InputGrid2D(
            x_mesh_flat=palm_grid_2d_shifted.x_mesh_flat[filter],
            y_mesh_flat=palm_grid_2d_shifted.y_mesh_flat[filter],
        )
        output_array = GridTransformer2D.transform(
            input_grid=input_grid_shifted,
            input_array=sliced_array_flat[filter],
            output_grid=palm_grid_2d,
            method='nearest',
        )

        mask = numpy.arange(output_array.shape[0]) > indices_below_output_terrain_height[:, numpy.newaxis]
        mask = numpy.transpose(mask, (1, 0))
        fill_value = -9999.0
        output_array = numpy.where(mask, output_array, fill_value)

        if mode in ['2d_left', '2d_right']:
            input_array[:, :, 0] = output_array
        elif mode in ['2d_south', '2d_north']:
            input_array[:, 0, :] = output_array
        else:
            raise NotImplementedError()

        return input_array
