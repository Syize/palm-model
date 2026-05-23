import numpy

from promet.api.grid_transformer_2d import GridTransformer2D
from promet.api.input_grid_2d import InputGrid2D
from promet.api.logging import PrometException
from promet.api.logging import print_and_log_section
from promet.api.logging import print_and_log_step
from promet.api.namelist import Namelist
from promet.api.static_driver import StaticDriver
from promet.api.palm_grid_2d import PALMGrid2D
from promet.api.palm_grid_3d import PALMGrid3D


class PALMSetup:

    required_global_attributes = [
        'origin_lat',
        'origin_lon',
        'origin_time',
        'origin_x',
        'origin_y',
        'origin_z',
        'rotation_angle',
    ]

    required_dimensions = [
        'x',
        'xu',
        'y',
        'yv',
        'z',
        'zw',
        'zsoil',
    ]

    grid_modes_for_lateral_boundaries = ['2d_left', '2d_right', '2d_south', '2d_north']

    def __init__(
            self,
            namelist: Namelist,
            static_driver: StaticDriver,
    ):
        self.namelist = namelist
        self.static_driver = static_driver
        self._global_attributes = dict()
        self._dimensions = dict()
        self._palm_terrain_grids = dict()
        self._palm_terrain_arrays = dict()
        self._input_terrain_arrays_on_palm_grid = dict()

    @property
    def global_attributes(self):
        return self._global_attributes

    @property
    def dimensions(self):
        return self._dimensions

    @property
    def palm_terrain_grids(self):
        return self._palm_terrain_grids

    @property
    def palm_terrain_arrays(self):
        return self._palm_terrain_arrays

    @property
    def input_terrain_arrays_on_palm_grid(self):
        return self._input_terrain_arrays_on_palm_grid

    def load_global_attributes(
            self,
            verbose: bool = False,
    ) -> bool:
        print_and_log_section('Loading global attributes from PALM setup files...')
        for attribute in self.required_global_attributes:
            try:
                self._global_attributes[attribute] = self.static_driver.get_global_attribute(attribute)
                source = 'static driver'
            except AttributeError as e:
                self._global_attributes[attribute] = self.namelist[attribute]
                source = 'namelist'
            if verbose:
                print_and_log_step(
                    f'Found global attribute "{attribute}" = "{self._global_attributes[attribute]}" in {source}'
                )
        return False

    def get_dimension(
            self,
            dimension: str,
            verbose: bool = False,
    ):
        source = None
        if dimension == 'z':
            nz = self.namelist.parameters['nz']
            try:
                result_static = self.static_driver.get_dimension(dimension)
                if result_static['size'] == nz:
                    result = result_static
                    source = 'static_driver'
                else:
                    raise KeyError()
            except KeyError as e:
                result = self.namelist.get_dimension(dimension)
                source = 'namelist'
        else:
            result = self.static_driver.get_dimension(dimension)
            source = 'static_driver'
        if verbose:
            if source == 'namelist':
                print_and_log_step(
                    f'Created dimension "{dimension}" with size={result["size"]} from namelist info'
                )
            elif source == 'static_driver':
                print_and_log_step(
                    f'Found dimension "{dimension}" with size={result["size"]} in static driver'
                )
            else:
                raise PrometException(
                    f'Missing dimension "{dimension}"',
                    id='f1Ofeg',
                )
        return result

    def load_dimensions(
            self,
            zsoil: list,
            verbose: bool = False,
    ) -> bool:
        print_and_log_section('Loading dimensions from PALM setup files...')
        missing_dimensions = []
        for dimension in self.required_dimensions:
            try:
                self._dimensions[dimension] = self.get_dimension(dimension, verbose)
            except KeyError as e:
                missing_dimensions.append(dimension)
        print_and_log_section('Computing missing PALM dimensions...')
        dimension_affiliation = dict(
            xu='x',
            yv='y',
            zw='z',
        )
        for dimension in missing_dimensions:
            if dimension in dimension_affiliation:
                affiliated_dimension = dimension_affiliation[dimension]
                if affiliated_dimension in self._dimensions:
                    self._dimensions[dimension] = dict(
                        size=self._dimensions[affiliated_dimension]["size"] - 1,
                        long_name=f'{affiliated_dimension} coordinate of cell faces',
                        units='m',
                        axis=affiliated_dimension.upper(),
                        values=self._dimensions[affiliated_dimension]["values"][1:] -
                               0.5 * numpy.diff(self._dimensions[affiliated_dimension]["values"]),
                    )
                if dimension in self._dimensions:
                    assert self._dimensions[dimension]["size"] == len(self._dimensions[dimension]["values"])
                    if verbose:
                        print_and_log_step(
                            f'Computed dimension "{dimension}" with size={self._dimensions[dimension]["size"]} from dimension "{affiliated_dimension}"'
                        )
                else:
                    if verbose:
                        print_and_log_step(
                            f'Still missing dimension {dimension}'
                        )
            elif dimension == 'zsoil':
                self._dimensions[dimension] = dict(
                    size=len(zsoil),
                    long_name='depth below land surface',
                    units='m',
                    axis='Z',
                    values=numpy.array(zsoil),
                )
                print_and_log_step(
                    f'Loaded dimension "{dimension}" with size={self._dimensions[dimension]["size"]} from config'
                )
            else:
                raise PrometException(
                    f'Missing dimension "{dimension}"',
                    id='B7HMhA',
                )
        return False

    def get_grid_2d(
            self,
            x_dim: str,
            y_dim: str,
            mode: str = '2d',
    ) -> PALMGrid2D:
        assert x_dim in ['x', 'xu']
        assert y_dim in ['y', 'yv']
        assert mode in self.grid_modes_for_lateral_boundaries + ['2d_top', '2d']
        x_palm = self.dimensions[x_dim]["values"]
        y_palm = self.dimensions[y_dim]["values"]
        if mode == '2d_left':
            x_palm = self.dimensions[x_dim]["values"][:1]
        if mode == '2d_right':
            x_palm = self.dimensions[x_dim]["values"][-1:]
        if mode == '2d_south':
            y_palm = self.dimensions[y_dim]["values"][:1]
        if mode == '2d_north':
            y_palm = self.dimensions[y_dim]["values"][-1:]

        return PALMGrid2D(
            x_palm=x_palm,
            y_palm=y_palm,
            origin_x=self.global_attributes['origin_x'],
            origin_y=self.global_attributes['origin_y'],
        )

    def get_grid(
            self,
            x_dim: str,
            y_dim: str,
            z_dim: str,
            mode: str = '3d',
    ) -> PALMGrid3D:
        assert x_dim in ['x', 'xu']
        assert y_dim in ['y', 'yv']
        assert z_dim in ['z', 'zw', 'zsoil']
        assert mode in self.grid_modes_for_lateral_boundaries + ['2d_top', '3d']
        x_palm = self.dimensions[x_dim]["values"]
        y_palm = self.dimensions[y_dim]["values"]
        z_palm = self.dimensions[z_dim]["values"]
        if mode == '2d_left':
            x_palm = self.dimensions[x_dim]["values"][:1]
        if mode == '2d_right':
            x_palm = self.dimensions[x_dim]["values"][-1:]
        if mode == '2d_south':
            y_palm = self.dimensions[y_dim]["values"][:1]
        if mode == '2d_north':
            y_palm = self.dimensions[y_dim]["values"][-1:]
        if mode == '2d_top':
            z_palm = self.dimensions[z_dim]["values"][-1:]
        return PALMGrid3D(
            x_palm=x_palm,
            y_palm=y_palm,
            z_palm=z_palm,
            origin_x=self.global_attributes['origin_x'],
            origin_y=self.global_attributes['origin_y'],
            origin_z=self.global_attributes['origin_z'] if z_dim != 'zsoil' else 0.0,
        )

    def set_palm_terrain(
            self,
            input_terrain_grid: InputGrid2D,
            input_terrain_array: numpy.array,
    ) -> None:
        grid_configs = [
            {"stagger": '', "x_dim": 'x', "y_dim": 'y'},
            {"stagger": 'X', "x_dim": 'xu', "y_dim": 'y'},
            {"stagger": 'Y', "x_dim": 'x', "y_dim": 'yv'},
            {"stagger": 'Z', "x_dim": 'x', "y_dim": 'y'},
        ]
        for grid_config in grid_configs:
            self._palm_terrain_grids[grid_config['stagger']] = self.get_grid_2d(
                x_dim=grid_config['x_dim'],
                y_dim=grid_config['y_dim'],
                mode='2d',
            )
        self._input_terrain_arrays_on_palm_grid[''] = GridTransformer2D.transform(
            input_grid=input_terrain_grid,
            input_array=input_terrain_array,
            output_grid=self.palm_terrain_grids[''],
            method='linear',
        )
        self._palm_terrain_arrays[''] = self.static_driver.get_terrain_array()
        for grid_config in grid_configs[1:]:
            self._input_terrain_arrays_on_palm_grid[grid_config['stagger']] = GridTransformer2D.transform(
                input_grid=input_terrain_grid,
                input_array=input_terrain_array,
                output_grid=self.palm_terrain_grids[grid_config['stagger']],
                method='linear',
            )
            self._palm_terrain_arrays[grid_config['stagger']] = GridTransformer2D.transform(
                input_grid=self.palm_terrain_grids[''],
                input_array=self._palm_terrain_arrays[''].flatten(order='F'),
                output_grid=self.palm_terrain_grids[grid_config['stagger']],
                method='nearest',
            )

    def get_boundary_terrain(
            self,
            mode: str,
            stagger: str = '',
    ):
        assert mode in self.grid_modes_for_lateral_boundaries
        if mode == '2d_left':
            return (
                self.input_terrain_arrays_on_palm_grid[stagger][:, 0],
                self.palm_terrain_arrays[stagger][:, 0],
            )
        if mode == '2d_right':
            return (
                self.input_terrain_arrays_on_palm_grid[stagger][:, -1],
                self.palm_terrain_arrays[stagger][:, -1],
            )
        if mode == '2d_south':
            return (
                self.input_terrain_arrays_on_palm_grid[stagger][0, :],
                self.palm_terrain_arrays[stagger][0, :],
            )
        if mode == '2d_north':
            return (
                self.input_terrain_arrays_on_palm_grid[stagger][-1, :],
                self.palm_terrain_arrays[stagger][-1, :],
            )

    def get_metadata(self) -> dict:
        return dict(
            namelist_metadata=self.namelist.get_metadata(),
            static_driver_metadata=self.static_driver.get_metadata(),
            global_attributes=self.global_attributes,
            terrain_bounds=dict(
                min=self.palm_terrain_arrays[''].min(),
                max=self.palm_terrain_arrays[''].max(),
            )
        )

