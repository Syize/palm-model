import pyproj
import numpy

from promet.api.config import Config
from promet.api.dynamic_driver import DynamicDriver
from promet.api.logging import print_and_log_section, PrometException, print_error
from promet.api.logging import print_and_log_subsection
from promet.api.logging import print_and_log_step
from promet.api.palm_setup import PALMSetup
from promet.api.plugin_input_data import PluginInputData
from promet.api.grid_transformer_3d import GridTransformer3D
from promet.api.terrain_mapper import TerrainMapper
from promet.api.input_grid_3d import InputGrid3D


class Plugin:

    input_data_class = PluginInputData
    filtered_time_slots = dict()
    variable_group_mapping = dict()

    def __init__(
            self,
            config: Config,
            dynamic_driver: DynamicDriver,
            palm_setup: PALMSetup,
    ):
        self.input_data_atmosphere = config['input_data']
        self.input_crs = config['input_crs']
        self.palm_crs = config['palm_crs']
        self.palm_variable_groups = config['palm_variable_groups']
        self.start_time = config['start_time']
        self.end_time = config['end_time']
        self.top_adjustment_index = config['top_adjustment_index']
        self.dynamic_driver = dynamic_driver
        self.palm_setup = palm_setup

        unknown_variable_groups = []
        for variable_group in self.palm_variable_groups:
            if variable_group not in self.dynamic_driver.allowed_variable_groups:
                unknown_variable_groups.append(variable_group)
        if len(unknown_variable_groups) > 0:
            print_error(f'The following variable groups are unknown:')
            for variable_group in unknown_variable_groups:
                print_error(f'   {variable_group}')
            raise PrometException(
                f'Found {len(unknown_variable_groups)} unknown variable group{"s" if len(unknown_variable_groups) > 1 else ""}.',
                id='MV7beQ',
            )

        unusable_variable_groups = []
        for variable_group in self.palm_variable_groups:
            if variable_group not in self.variable_group_mapping:
                unusable_variable_groups.append(variable_group)
        if len(unusable_variable_groups) > 0:
            print_error(f'The following variable groups are not usable with the current plugin:')
            for variable_group in unusable_variable_groups:
                print_error(f'   {variable_group}')
            raise PrometException(
                f'Found {len(unusable_variable_groups)} usable variable group{"s" if len(unusable_variable_groups) > 1 else ""}.',
                id='wnRJPw',
            )

        try:
            self.crs_input = pyproj.CRS.from_string(self.input_crs)
        except pyproj.exceptions.CRSError as e:
            raise PrometException(f'{self.input_crs} unknown', id='GyZ30Q')
        try:
            self.crs_output = pyproj.CRS.from_string(self.palm_crs)
        except pyproj.exceptions.CRSError as e:
            raise PrometException(f'{self.palm_crs} unknown', id='wZHdrg')

        self.transformer_in2out = pyproj.Transformer.from_crs(self.crs_input, self.crs_output, always_xy=True)
        self.transformer_out2in = pyproj.Transformer.from_crs(self.crs_output, self.crs_input, always_xy=True)
        self._input_data = None

    @property
    def input_data(self) -> PluginInputData:
        return self._input_data

    @input_data.setter
    def input_data(self, input_data: PluginInputData):
        self._input_data = input_data

    def initialize(self, verbose: bool = False) -> bool:
        self.input_data = self.input_data_class(
            filepath=self.input_data_atmosphere,
        )
        self.input_data.discover_dataset()

        palm_grid_xu = self.palm_setup.get_grid(x_dim='xu', y_dim='y', z_dim='z')
        palm_grid_yv = self.palm_setup.get_grid(x_dim='x', y_dim='yv', z_dim='z')
        palm_grid_zw = self.palm_setup.get_grid(x_dim='x', y_dim='y', z_dim='zw')
        palm_grid_zsoil = self.palm_setup.get_grid(x_dim='x', y_dim='y', z_dim='zsoil')
        self.input_data.load_grids(
            transformer=self.transformer_in2out,
            bounds_x=(palm_grid_xu.x[0], palm_grid_xu.x[-1]),
            bounds_y=(palm_grid_yv.y[0], palm_grid_yv.y[-1]),
            bounds_z=(0, palm_grid_zw.z[-1]),
            bounds_zsoil=(0, palm_grid_zsoil.z[-1]),
        )

        print_and_log_section('Scanning for eligible input data time_slots...')
        for file_index, date_time in enumerate(self.input_data.available_time_slots):
            if self.start_time <= date_time <= self.end_time:
                self.filtered_time_slots[file_index] = dict(
                    date_time=date_time,
                    start_time_delta=(date_time - self.start_time).total_seconds(),
                )
        if len(self.filtered_time_slots) == 0:
            raise PrometException('No matching input data timeslots were found.', id='dVmArw')
        self.initial_time_slot = min(self.filtered_time_slots.keys())
        if verbose:
            for n, time_slot in self.filtered_time_slots.items():
                print_and_log_step(
                    f'Using time_slot({n}) with date_time: {time_slot["date_time"]}, which is palm_time: {time_slot["start_time_delta"]} seconds'
                )

        self.input_terrain_frame = self.input_data.get_terrain_frame(
            time_index=self.initial_time_slot,
        )
        input_terrain_array = self.input_terrain_frame['array'].flatten(order='F')
        self.palm_setup.set_palm_terrain(
            input_terrain_grid=self.input_terrain_frame['grid'],
            input_terrain_array=input_terrain_array,
        )
        return False

    def process_palm_variables(self, verbose: bool = False) -> bool:
        print_and_log_section('Populating dynamic driver with variables...')

        # first set time values
        self.dynamic_driver.create_time(
            time_slots=self.filtered_time_slots,
        )

        # make sure all variables are created
        for group_name in self.palm_variable_groups:
            group_variables = self.dynamic_driver.get_group_variables(group_name)
            for var_name in group_variables:
                self.dynamic_driver.create_variable(
                    var_name=var_name,
                    verbose=verbose,
                )

        print_and_log_section('Processing initial variables...')
        variable_groups = [v for v in self.palm_variable_groups if v in self.dynamic_driver.get_initial_variable_groups()]
        self.process_timeframe_for_variable_groups(
            variable_groups=variable_groups,
            time_index=self.initial_time_slot,
            verbose=verbose,
        )

        print_and_log_section('Processing forcing variables...')
        for n, time_slot in self.filtered_time_slots.items():
            print_and_log_subsection(
                f'Loading Data for date_time: {time_slot["date_time"]}, which is palm_time: {time_slot["start_time_delta"]} seconds'
            )
            variable_groups = [v for v in self.palm_variable_groups if v in self.dynamic_driver.get_forcing_variable_groups()]
            self.process_timeframe_for_variable_groups(
                variable_groups=variable_groups,
                time_index=n,
                verbose=verbose,
            )

        print_and_log_section('Finished processing all variables!')
        return False

    def process_timeframe_for_variable_groups(
            self,
            variable_groups: list,
            time_index: int,
            verbose: bool = False,
    ):
        for group_name in variable_groups:
            var_frame = self.input_data.get_variable_timeframe(
                input_variable=self.variable_group_mapping[group_name]['input_variable'],
                time_index=time_index,
                terrain_extended=True,
            )
            input_grid = var_frame['grid']
            input_array = var_frame['array'].flatten(order='F')

            stagger = self.dynamic_driver.get_group_stagger(group_name=group_name)

            if isinstance(input_grid, InputGrid3D):
                mapping_mode = 'grid_transformer_3d'
            elif input_grid is None:
                mapping_mode = 'copy'
            else:
                raise NotImplementedError('Unsupported input grid')

            group_type = self.dynamic_driver.get_group_type(group_name=group_name)
            group_variables = self.dynamic_driver.get_group_variables(group_name=group_name)
            group_grid_dimensions = self.dynamic_driver.get_group_grid_dimensions(group_name=group_name)
            for var_name in group_variables:
                grid_mode = self.dynamic_driver.allowed_variables[var_name]['grid_mode']
                if mapping_mode == 'grid_transformer_3d':
                    output_grid = self.palm_setup.get_grid(
                        x_dim=group_grid_dimensions[2],
                        y_dim=group_grid_dimensions[1],
                        z_dim=group_grid_dimensions[0],
                        mode=grid_mode,
                    )
                    output_array = GridTransformer3D.transform(
                        input_grid=input_grid,
                        input_array=input_array,
                        output_grid=output_grid,
                    )

                    if grid_mode in self.palm_setup.grid_modes_for_lateral_boundaries:
                        input_boundary_terrain, palm_boundary_terrain = self.palm_setup.get_boundary_terrain(
                            mode=grid_mode,
                            stagger=stagger,
                        )
                        output_array = TerrainMapper.apply(
                            input_grid=output_grid,
                            input_array=output_array,
                            input_terrain=input_boundary_terrain,
                            output_terrain=palm_boundary_terrain,
                            mode=grid_mode,
                            top_adjustment_index=self.top_adjustment_index,
                        )
                elif mapping_mode == 'copy':
                    output_array = input_array
                else:
                    raise NotImplementedError('Unsupported mapping mode')
                if group_type == 'initial':
                    self.dynamic_driver.set_init_variable(
                        var_name=var_name,
                        array=output_array,
                        verbose=verbose,
                    )
                elif group_type == 'ls_forcing':
                    self.dynamic_driver.set_forcing_variable(
                        var_name=var_name,
                        time_index=time_index - self.initial_time_slot,
                        array=output_array,
                        verbose=verbose,
                    )
                else:
                    raise NotImplementedError(f'Unknown group type {group_type}')

    def get_metadata(self) -> dict:
        return dict(
            input_metadata=self.input_data.get_metadata(),
            filtered_time_slots=self.filtered_time_slots,
        )
