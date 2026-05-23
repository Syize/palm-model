import glob
import netCDF4
import pyproj
from datetime import datetime
from datetime import timedelta
from datetime import timezone
import numpy
import re

from promet.api.logging import PrometException
from promet.api.logging import print_and_log_section
from promet.api.logging import print_and_log_step
from promet.api.plugin_input_data import PluginInputData

from .input_grid_2d import ICONInputGrid2D
from .input_grid_3d import ICONInputGrid3D


def extend_vertically(cells_1d, nz):
    assert len(cells_1d.shape) == 1
    cells_2d_expanded = numpy.expand_dims(cells_1d, axis=0)
    cells_2d_tileed = numpy.tile(cells_2d_expanded, (nz, 1))
    return cells_2d_tileed


def extend_to_surface(cells_2d, surface_1d):
    assert len(cells_2d.shape) == 2
    assert len(surface_1d.shape) == 1
    assert cells_2d.shape[1] == surface_1d.shape[0]
    surface_2d_expanded = numpy.expand_dims(surface_1d, axis=0)
    cells_2d_surfaced = numpy.concatenate((cells_2d, surface_2d_expanded), axis=0)
    return cells_2d_surfaced


class ICONInputData(PluginInputData):

    g = 9.81

    def __init__(
            self,
            filepath: str,
    ):
        print_and_log_section(f'Loading ICON files from glob: {filepath}')
        self.icon_datafiles = list(sorted(glob.glob(filepath)))
        if len(self.icon_datafiles) == 0:
            raise PrometException(f'No ICON data files could be found that match given glob', id='OPOSDw')
        self.grids = dict()
        self.variables = dict()
        self.input_grids_3d = dict()
        self.input_grids_3d_terrain_extended = dict()
        self.crop_mask_2d = dict()
        self.crop_mask_2d_terrain_extended = dict()

    def load_grids(
            self,
            transformer: pyproj.Transformer,
            bounds_x: tuple,
            bounds_y: tuple,
            bounds_z: tuple,
            bounds_zsoil: tuple,
    ):
        with netCDF4.Dataset(self.icon_datafiles[0], mode='r') as ncfile:
            if 'tlon' in ncfile.variables:
                longitude = ncfile.variables['tlon'][0, :]
            elif 'CLON' in ncfile.variables:
                longitude = ncfile.variables['CLON'][0, :]
            else:
                raise PrometException(
                    f'ICON grid longitude not found in input data. One of "tlon" or "CLON" must be provided.',
                    id='Yddrlg',
                )
            if 'tlat' in ncfile.variables:
                latitude = ncfile.variables['tlat'][0, :]
            elif 'CLAT' in ncfile.variables:
                latitude = ncfile.variables['CLAT'][0, :]
            else:
                raise PrometException(
                    f'ICON grid latitude not found in input data. One of "tlat" or "CLAT" must be provided.',
                    id='Yddrlg',
                )
            if 'HSURF' in ncfile.variables:
                hsurf = ncfile.variables['HSURF'][0, :]
            else:
                raise PrometException(f'ICON variable "HSURF" not found in input data.', id='Yddrlg')
            if 'HHL' in ncfile.variables:
                hhl = ncfile.variables['HHL'][0, :, :]
            else:
                raise PrometException(f'ICON variable "HHL" not found in input data.', id='Yddrlg')
            if 'T_SO' in ncfile.variables:
                depth_t_so_name = ncfile.variables['T_SO'].dimensions[1]
                depth_t_so = ncfile.variables[depth_t_so_name][:]
            else:
                raise PrometException(f'ICON variable "T_SO" not found in input data.', id='Yddrlg')
            if 'W_SO' in ncfile.variables:
                depth_w_so_name = ncfile.variables['W_SO'].dimensions[1]
                depth_w_so = ncfile.variables[depth_w_so_name][:]
            else:
                raise PrometException(f'ICON variable "W_SO" not found in input data.', id='Yddrlg')

        x_cells, y_cells = transformer.transform(
            xx=longitude,
            yy=latitude,
        )
        self.x_cells = x_cells
        self.y_cells = y_cells

        self.zw_2d = hhl
        self.z_2d = self.zw_2d[:-1, :] + numpy.diff(self.zw_2d, axis=0) * 0.5

        self.hsurf = hsurf
        self.depth_t_so = depth_t_so
        self.depth_w_so = depth_w_so
        self.d_depth = (self.depth_t_so[1:] - self.depth_w_so) * 2.0
        self.zw_2d_te = extend_to_surface(
            cells_2d=self.zw_2d,
            surface_1d=self.hsurf,
        )
        self.z_2d_te = extend_to_surface(
            cells_2d=self.z_2d,
            surface_1d=self.hsurf,
        )

        x = numpy.asarray(self.x_cells[0:100])
        y = numpy.asarray(self.y_cells[0:100])
        diff_x = x[:, numpy.newaxis] - x[numpy.newaxis, :]
        diff_y = y[:, numpy.newaxis] - y[numpy.newaxis, :]
        dist_squared = diff_x ** 2 + diff_y ** 2
        numpy.fill_diagonal(dist_squared, numpy.inf)
        h_distance = numpy.sqrt(numpy.min(dist_squared))
        bounds_size_h = numpy.sqrt((bounds_x[1] - bounds_x[0]) ** 2 + (bounds_y[1] - bounds_y[0]) ** 2)
        size_factor = h_distance / bounds_size_h

        buffer_h_factor = 2.0 * size_factor
        buffer_v_factor = 0.2
        buffer_x = buffer_h_factor * (bounds_x[1] - bounds_x[0])
        buffer_y = buffer_h_factor * (bounds_y[1] - bounds_y[0])
        buffer_z = buffer_v_factor * (bounds_z[1] - bounds_z[0])
        mask_x = (self.x_cells >= (bounds_x[0] - buffer_x)) & (self.x_cells <= (bounds_x[1] + buffer_x))
        mask_y = (self.y_cells >= (bounds_y[0] - buffer_y)) & (self.y_cells <= (bounds_y[1] + buffer_y))
        mask_zw = (self.zw_2d >= (bounds_z[0] - buffer_z)) & (self.zw_2d <= (bounds_z[1] + buffer_z))
        mask_zw_te = (self.zw_2d_te >= (bounds_z[0] - buffer_z)) & (self.zw_2d_te <= (bounds_z[1] + buffer_z))
        mask_z = (self.z_2d >= (bounds_z[0] - buffer_z)) & (self.z_2d <= (bounds_z[1] + buffer_z))
        mask_z_te = (self.z_2d_te >= (bounds_z[0] - buffer_z)) & (self.z_2d_te <= (bounds_z[1] + buffer_z))

        self.crop_mask_1d = mask_x & mask_y

        self.crop_mask_zw_2d = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.zw_2d.shape[0],
        ) & mask_zw
        self.crop_mask_zw_2d_flat = self.crop_mask_zw_2d.flatten(order='F')

        self.crop_mask_zw_2d_te = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.zw_2d_te.shape[0],
        ) & mask_zw_te
        self.crop_mask_zw_2d_te_flat = self.crop_mask_zw_2d_te.flatten(order='F')

        self.crop_mask_z_2d = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.z_2d.shape[0],
        ) & mask_z
        self.crop_mask_z_2d_flat = self.crop_mask_z_2d.flatten(order='F')

        self.crop_mask_z_2d_te = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.z_2d_te.shape[0],
        ) & mask_z_te
        self.crop_mask_z_2d_te_flat = self.crop_mask_z_2d_te.flatten(order='F')

        self.crop_mask_depth_t_so_2d = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.depth_t_so.shape[0],
        )
        self.crop_mask_depth_t_so_2d_flat = self.crop_mask_depth_t_so_2d.flatten(order='F')

        self.crop_mask_depth_w_so_2d = extend_vertically(
            cells_1d=self.crop_mask_1d,
            nz=self.depth_w_so.shape[0],
        )
        self.crop_mask_depth_w_so_2d_flat = self.crop_mask_depth_w_so_2d.flatten(order='F')

        self.input_grids_3d[''] = self.get_grid_3d(stagger='')
        self.input_grids_3d['Z'] = self.get_grid_3d(stagger='Z')
        self.input_grids_3d_terrain_extended[''] = self.get_grid_3d(stagger='', terrain_extended=True)
        self.input_soil_grid_3d = self.get_soil_grid_3d()
        self.input_soil_grid_3d_terrain_extended = self.get_soil_grid_3d(terrain_extended=True)

        self.input_terrain_grid_3d = self.get_terrain_grid_3d()
        self.input_terrain_frame = self.get_variable_timeframe(
            input_variable='HSURF',
            time_index=0,
        )

    def discover_dataset(self):
        for icon_datafile in self.icon_datafiles:
            print_and_log_step(f'Parsing ICON file {icon_datafile}')
            with netCDF4.Dataset(icon_datafile, mode='r') as ncfile:
                required_global_attributes = [
                    'number_of_grid_used'
                ]
                for required_global_attribute in required_global_attributes:
                    self.attributes[required_global_attribute] = ncfile.getncattr(required_global_attribute)
                time_unit_string = ncfile.variables['time'].units
                match = re.search(r"(\w+) since (\d{4}-\d{1,2}-\d{1,2} \d{2}:\d{2}:\d{2})", time_unit_string)
                if not match:
                    raise PrometException(
                        "Invalid time unit format. Expected format: '{unit} since YYYY-M-D HH:MM:SS'",
                        id='k42Up6',
                    )
                time_unit, reference_time_str = match.groups()
                if time_unit not in ['hours', 'minutes', 'seconds']:
                    raise PrometException(
                        "Invalid time unit format. Expected format: '{unit} since YYYY-M-D HH:MM:SS',"
                        " where unit must be one of ['hours', 'minutes', 'seconds']",
                        id='k42Up6',
                    )
                try:
                    base_time = datetime.strptime(reference_time_str, "%Y-%m-%d %H:%M:%S")
                except ValueError:
                    raise PrometException(
                        f"Invalid time unit format. Expected format: '{time_unit} since YYYY-M-D HH:MM:SS'",
                        id='k42Up6',
                    )
                base_time = base_time.replace(tzinfo=timezone.utc)
                time = base_time + timedelta(**{time_unit: float(ncfile.variables['time'][0])})
                self.available_time_slots.append(time)

        print('')
        for n, time_slot in enumerate(self.available_time_slots):
            print_and_log_step(f'Found time_slot({n}) with date_time: {time_slot}')

    @staticmethod
    def get_netcdf_variable(variable_name: str, ncfile: netCDF4.Dataset, alternative_variable_names: list = None):
        eligible_variable_names = [variable_name, variable_name.lower(), variable_name.upper()]
        eligible_variable_names.extend([f'{v}_2' for v in eligible_variable_names])
        if isinstance(alternative_variable_names, list):
            for alternative_variable_name in alternative_variable_names:
                if isinstance(alternative_variable_name, str):
                    if alternative_variable_name not in eligible_variable_names:
                        eligible_variable_names.extend([alternative_variable_name, alternative_variable_name.lower(), alternative_variable_name.upper()])
        existing_variable_names = list(set(eligible_variable_names) & set(ncfile.variables))
        validated_variable_names = []
        for vn in existing_variable_names:
            if not ncfile.variables[vn].dimensions[1].startswith('plev'):
                validated_variable_names.append(vn)
        if len(validated_variable_names) == 0:
            raise PrometException(f'ICON variable "{variable_name}" not found in input data.', id='Yddrlg')
        if variable_name in validated_variable_names:
            return variable_name
        else:
            print_and_log_step(f'ICON variable "{variable_name}" not found. Using "{validated_variable_names[0]}" instead')
            return validated_variable_names[0]

    def get_variable_timeframe(
            self,
            input_variable,
            time_index,
            terrain_extended=False,
    ):
        with netCDF4.Dataset(self.icon_datafiles[time_index], mode='r') as ncfile:
            if input_variable == 'QV':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile, alternative_variable_names=['q'])
                input_array = ncfile.variables[netcdf_variable_name][0, :, :] / (1 - ncfile.variables[netcdf_variable_name][0, :, :]) # compute mixing-ratio
                input_array = extend_to_surface(
                    cells_2d=input_array,
                    surface_1d=input_array[-1, :],
                )
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_z_2d_te_flat
                timeframe = dict(
                    grid=self.input_grids_3d_terrain_extended[''],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'T':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                netcdf_variable_name_pres = self.get_netcdf_variable('P', ncfile, alternative_variable_names=['pres'])
                input_array = ncfile.variables[netcdf_variable_name][0, :, :] * (100000.0 / ncfile.variables[netcdf_variable_name_pres][0, :, :])**0.286 # compute potential temperature
                input_array = extend_to_surface(
                    cells_2d=input_array,
                    surface_1d=input_array[-1, :],
                )
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_z_2d_te_flat
                timeframe = dict(
                    grid=self.input_grids_3d_terrain_extended[''],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'U':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                input_array = ncfile.variables[netcdf_variable_name][0, :, :]
                input_array = extend_to_surface(
                    cells_2d=input_array,
                    surface_1d=numpy.zeros_like(input_array[-1, :]),
                )
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_z_2d_te_flat
                timeframe = dict(
                    grid=self.input_grids_3d_terrain_extended[''],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'V':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                input_array = ncfile.variables[netcdf_variable_name][0, :, :]
                input_array = extend_to_surface(
                    cells_2d=input_array,
                    surface_1d=numpy.zeros_like(input_array[-1, :]),
                )
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_z_2d_te_flat
                timeframe = dict(
                    grid=self.input_grids_3d_terrain_extended[''],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'W':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile, alternative_variable_names=['wz'])
                input_array = ncfile.variables[netcdf_variable_name][0, :, :]
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_zw_2d_flat
                timeframe = dict(
                    grid=self.input_grids_3d['Z'],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'ART_O3':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                input_array = ncfile.variables[netcdf_variable_name][0, :, :] * 0.6035e6  # ozone mass mixing ratio (Kg/Kg) to ppm
                input_array = extend_to_surface(
                    cells_2d=input_array,
                    surface_1d=input_array[-1, :],
                )
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_z_2d_te_flat
                timeframe = dict(
                    grid=self.input_grids_3d_terrain_extended[''],
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'T_SO':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                input_array = ncfile.variables[netcdf_variable_name][0, :, :]
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_depth_t_so_2d_flat
                timeframe = dict(
                    grid=self.input_soil_grid_3d_terrain_extended,
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'W_SO':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                input_array = ncfile.variables[netcdf_variable_name][0, :, :] * 0.001 / self.d_depth[:, numpy.newaxis]  # Kg/m2 to m3/m3
                input_array_flat = input_array.flatten(order='F')
                crop_mask_flat = self.crop_mask_depth_w_so_2d_flat
                timeframe = dict(
                    grid=self.input_soil_grid_3d,
                    array=input_array_flat[crop_mask_flat],
                    stagger='',
                )
            elif input_variable == 'HSURF':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile)
                timeframe = dict(
                    grid=self.get_grid_2d(),
                    array=ncfile.variables[netcdf_variable_name][0, self.crop_mask_1d],
                    stagger='',
                )
            elif input_variable == 'PS':
                netcdf_variable_name = self.get_netcdf_variable(input_variable, ncfile, alternative_variable_names=['sp'])
                timeframe = dict(
                    grid=self.get_grid_2d(),
                    array=ncfile.variables[netcdf_variable_name][0, self.crop_mask_1d],
                    stagger='',
                )
                timeframe['grid'] = None
                timeframe['array'] = numpy.mean(timeframe['array'])
            else:
                raise PrometException(f'Unknown input variable "{input_variable}".', id='ohHDsw')
        return timeframe

    def get_grid_2d(
            self,
    ) -> ICONInputGrid2D:
        return ICONInputGrid2D.centered(
            x_cells=self.x_cells[self.crop_mask_1d],
            y_cells=self.y_cells[self.crop_mask_1d],
        )

    def get_terrain_array_2d(self) -> ICONInputGrid3D:
        return self.hsurf[self.crop_mask_1d]

    def get_terrain_frame(self, time_index):
        return self.get_variable_timeframe(
            input_variable='HSURF',
            time_index=time_index,
        )

    def get_terrain_grid_3d(self) -> ICONInputGrid3D:
        return ICONInputGrid3D(
            x_mesh_flat=self.x_cells[self.crop_mask_1d],
            y_mesh_flat=self.y_cells[self.crop_mask_1d],
            z_mesh_flat=self.hsurf[self.crop_mask_1d],
        )

    def get_soil_grid_3d(
            self,
            terrain_extended: bool = False,
    ) -> ICONInputGrid3D:
        if terrain_extended:
            return ICONInputGrid3D.centered_zsoil(
                x_cells=self.x_cells[self.crop_mask_1d],
                y_cells=self.y_cells[self.crop_mask_1d],
                zsoil1d=self.depth_t_so,
            )
        else:
            return ICONInputGrid3D.centered_zsoil(
                x_cells=self.x_cells[self.crop_mask_1d],
                y_cells=self.y_cells[self.crop_mask_1d],
                zsoil1d=self.depth_w_so,
            )

    def get_grid_3d(
            self,
            stagger: str = '',
            terrain_extended: bool = False,
    ) -> ICONInputGrid3D:
        assert stagger in ['', 'Z']
        if stagger == '':
            if terrain_extended:
                z_cells = self.z_2d_te
                crop_mask = self.crop_mask_z_2d_te
            else:
                z_cells = self.z_2d
                crop_mask = self.crop_mask_z_2d
        elif stagger == 'Z':
            if terrain_extended:
                z_cells = self.zw_2d_te
                crop_mask = self.crop_mask_zw_2d_te
            else:
                z_cells = self.zw_2d
                crop_mask = self.crop_mask_zw_2d
        else:
            raise PrometException(f'Unknown stagger: {stagger}', id='jdbR1A')
        x_cells_tiled = extend_vertically(cells_1d=self.x_cells, nz=z_cells.shape[0])
        x_flat = x_cells_tiled.flatten(order='F')
        y_cells_tiled = extend_vertically(cells_1d=self.y_cells, nz=z_cells.shape[0])
        y_flat = y_cells_tiled.flatten(order='F')
        z_flat = z_cells.flatten(order='F')
        crop_mask_flat = crop_mask.flatten(order='F')
        x_flat_croped = x_flat[crop_mask_flat]
        y_flat_croped = y_flat[crop_mask_flat]
        z_flat_croped = z_flat[crop_mask_flat]
        return ICONInputGrid3D(
            x_mesh_flat=x_flat_croped,
            y_mesh_flat=y_flat_croped,
            z_mesh_flat=z_flat_croped,
        )

    def get_metadata(self) -> dict:
        return dict(
            attributes=self.attributes,
            available_time_slots=self.available_time_slots,
            input_data_shape=self.zw_2d.shape,
            input_data_cropping=[self.crop_mask_1d.sum(), self.crop_mask_zw_2d.sum()],
        )
