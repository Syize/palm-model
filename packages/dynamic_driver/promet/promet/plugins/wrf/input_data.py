import glob
import netCDF4
import pyproj
from datetime import datetime
from datetime import timezone
import numpy

from promet.api.logging import PrometException
from promet.api.logging import print_and_log_section
from promet.api.logging import print_and_log_step
from promet.api.plugin_input_data import PluginInputData

from .input_grid_2d import WRFInputGrid2D
from .input_grid_3d import WRFInputGrid3D


class WRFInputData(PluginInputData):

    def __init__(
            self,
            filepath: str,
    ):
        print_and_log_section(f'Loading WRF files from glob: {filepath}')
        self.wrf_datafiles = list(sorted(glob.glob(filepath)))
        if len(self.wrf_datafiles) == 0:
            raise PrometException(f'No WRF data files could be found that match given glob', id='DFvnVA')
        self.grids = dict()
        self.variables = dict()
        self.input_grids_3d = dict()
        self.input_grids_3d_terrain_extended = dict()

    def load_grids(
            self,
            transformer: pyproj.Transformer,
            bounds_x: tuple,
            bounds_y: tuple,
            bounds_z: tuple,
            bounds_zsoil: tuple,
    ):
        with netCDF4.Dataset(self.wrf_datafiles[0], mode='r') as ncfile:
            xlong = ncfile.variables['XLONG'][0, :, :]
            xlat = ncfile.variables['XLAT'][0, :, :]
            hgt = ncfile.variables['HGT'][0, :, :]
            ph = ncfile.variables['PH'][0, :, :, :]
            phb = ncfile.variables['PHB'][0, :, :, :]
            zsoil = ncfile.variables['ZS'][0, :]

        x2d_in_t_flat, y2d_in_t_flat = transformer.transform(
            xx=xlong.flatten(order='F'),
            yy=xlat.flatten(order='F'),
        )
        self.x2d_in_t = x2d_in_t_flat.reshape(xlong.shape, order='F')
        self.y2d_in_t = y2d_in_t_flat.reshape(xlat.shape, order='F')
        g = 9.81
        self.zw_3d = (ph + phb) * (1. / g)
        self.z_3d = self.zw_3d[:-1, :, :] + numpy.diff(self.zw_3d, axis=0)
        for k in range(len(self.z_3d)):
            self.z_3d[k, :, :] = (((ph[k, :, :] + phb[k+1, :, :])*0.5) + ((ph[k+1, :, :] + phb[k, :, :])*0.5)) * (1. / g)
        self.hgt = hgt
        self.zsoil = zsoil

        self.i_min = 0
        self.i_max = 0
        for i in range(self.x2d_in_t.shape[1]):
            if numpy.all(self.x2d_in_t[:, i] < bounds_x[0]):
                self.i_min = i
            if numpy.any(self.x2d_in_t[:, i] <= bounds_x[1]):
                self.i_max = i + 2

        self.j_min = 0
        self.j_max = 0
        for j in range(self.y2d_in_t.shape[0]):
            if numpy.all(self.y2d_in_t[j, self.i_min:self.i_max] < bounds_y[0]):
                self.j_min = j
            if numpy.any(self.y2d_in_t[j, self.i_min:self.i_max] <= bounds_y[1]):
                self.j_max = j + 2

        self.i_min = 0
        self.i_max = 0
        for i in range(self.x2d_in_t.shape[1]):
            if numpy.all(self.x2d_in_t[self.j_min:self.j_max, i] < bounds_x[0]):
                self.i_min = i
            if numpy.any(self.x2d_in_t[self.j_min:self.j_max, i] <= bounds_x[1]):
                self.i_max = i + 2

        self.j_min = 0
        self.j_max = 0
        for j in range(self.y2d_in_t.shape[0]):
            if numpy.all(self.y2d_in_t[j, self.i_min:self.i_max] < bounds_y[0]):
                self.j_min = j
            if numpy.any(self.y2d_in_t[j, self.i_min:self.i_max] <= bounds_y[1]):
                self.j_max = j + 2

        self.k_min = 0
        self.k_max = 0
        for k in range(self.zw_3d.shape[0]):
            if numpy.all(self.zw_3d[k, :, :] < bounds_z[0]):
                self.k_min = k
            if numpy.any(self.zw_3d[k, :, :] <= bounds_z[1]):
                self.k_max = k + 1

        self.ks_min = 0
        self.ks_max = 0
        for k in range(self.zsoil.shape[0]):
            if self.zsoil[k] < bounds_zsoil[0]:
                self.ks_min = k
            if self.zsoil[k] <= bounds_zsoil[1]:
                self.ks_max = k + 1
        print_and_log_step(
            f'Cropping atmospheric input data [bottom_top, south_north, west_east] '
            f'\n      with shape={self.zw_3d.shape} '
            f'to [{self.k_min}:{self.k_max}, {self.j_min}:{self.j_max}, {self.i_min}:{self.i_max}]',
        )
        print_and_log_step(
            f'Cropping soil-layer input data: [soil_layers_stag, south_north, west_east] '
            f'\n      with shape={(self.zsoil.shape[0], self.zw_3d.shape[1], self.zw_3d.shape[2])} '
            f'to [{self.ks_min}:{self.ks_max}, {self.j_min}:{self.j_max}, {self.i_min}:{self.i_max}]',
        )

        self.input_grids_3d[''] = self.get_grid_3d(stagger='')
        self.input_grids_3d['X'] = self.get_grid_3d(stagger='X')
        self.input_grids_3d['Y'] = self.get_grid_3d(stagger='Y')
        self.input_grids_3d['Z'] = self.get_grid_3d(stagger='Z')
        self.input_grids_3d_terrain_extended[''] = self.get_grid_3d(stagger='', terrain_extended=True)
        self.input_grids_3d_terrain_extended['X'] = self.get_grid_3d(stagger='X', terrain_extended=True)
        self.input_grids_3d_terrain_extended['Y'] = self.get_grid_3d(stagger='Y', terrain_extended=True)
        self.input_soil_grid_3d = self.get_soil_grid_3d(terrain_extended=False)
        self.input_soil_grid_3d_terrain_extended = self.get_soil_grid_3d(terrain_extended=True)

        self.input_terrain_grid_3d = self.get_terrain_grid_3d()
        self.input_terrain_frame = self.get_variable_timeframe(
            input_variable='HGT',
            time_index=0,
        )

    def discover_dataset(self):
        for wrf_datafile in self.wrf_datafiles:
            print_and_log_step(f'Parsing WRF file {wrf_datafile}')
            with netCDF4.Dataset(wrf_datafile, mode='r') as ncfile:
                required_global_attributes = [
                    'TRUELAT1',
                    'TRUELAT2',
                    'MOAD_CEN_LAT',
                    'STAND_LON',
                    'CEN_LON',
                    'CEN_LAT',
                    'DX',
                    'DY',
                    'WEST-EAST_GRID_DIMENSION',
                    'SOUTH-NORTH_GRID_DIMENSION',
                ]
                for required_global_attribute in required_global_attributes:
                    self.attributes[required_global_attribute] = ncfile.getncattr(required_global_attribute)
                time_string = ncfile.variables['Times'][0].tobytes().decode('utf-8')
                try:
                    time = datetime.strptime(time_string, '%Y-%m-%d_%H:%M:%S')
                except ValueError:
                    raise PrometException(
                        "Invalid time format. Expected format: 'YYYY-MM-DD_HH:MM:SS'",
                        id='hT6S9k',
                    )
                time = time.replace(tzinfo=timezone.utc)
                self.available_time_slots.append(time)

        print('')
        for n, time_slot in enumerate(self.available_time_slots):
            print_and_log_step(f'Found time_slot({n}) with date_time: {time_slot}')

    def get_variable_timeframe(
            self,
            input_variable,
            time_index,
            terrain_extended=False,
    ):
        with netCDF4.Dataset(self.wrf_datafiles[time_index], mode='r') as ncfile:
            if input_variable not in ncfile.variables:
                raise PrometException(f'WRF variable "{input_variable}" not found in input data.', id='dqpKZA')
            if input_variable == 'QVAPOR':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['QVAPOR'].stagger],
                    array=ncfile.variables['QVAPOR'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['QVAPOR'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'T':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['T'].stagger],
                    array=ncfile.variables['T'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max] + 300.0,
                    stagger=ncfile.variables['T'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'U':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['U'].stagger],
                    array=ncfile.variables['U'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max+1],
                    stagger=ncfile.variables['U'].stagger,
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                dirichlet_bc_array = numpy.zeros_like(timeframe['array'][0:1, :, :])
                timeframe['array'] = numpy.concatenate((dirichlet_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'V':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['V'].stagger],
                    array=ncfile.variables['V'][0, self.k_min:self.k_max, self.j_min:self.j_max+1, self.i_min:self.i_max],
                    stagger=ncfile.variables['V'].stagger,
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                dirichlet_bc_array = numpy.zeros_like(timeframe['array'][0:1, :, :])
                timeframe['array'] = numpy.concatenate((dirichlet_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'W':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['W'].stagger],
                    array=ncfile.variables['W'][0, self.k_min:self.k_max+1, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['W'].stagger,
                )
            elif input_variable == 'no':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['no'].stagger],
                    array=ncfile.variables['no'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['no'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'no2':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['no2'].stagger],
                    array=ncfile.variables['no2'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['no2'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'o3':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['o3'].stagger],
                    array=ncfile.variables['o3'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['o3'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'PM10':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['PM10'].stagger],
                    array=ncfile.variables['PM10'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max] * 1.0E-9,
                    stagger=ncfile.variables['PM10'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'PM2_5_DRY':
                timeframe = dict(
                    grid=self.input_grids_3d[ncfile.variables['PM2_5_DRY'].stagger],
                    array=ncfile.variables['PM2_5_DRY'][0, self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max] * 1.0E-9,
                    stagger=ncfile.variables['PM2_5_DRY'].stagger,
                    bc_b='neumann',
                )
                timeframe['grid'] = self.input_grids_3d_terrain_extended[timeframe['stagger']]
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'SMOIS':
                timeframe = dict(
                    grid=self.input_soil_grid_3d,
                    array=ncfile.variables['SMOIS'][0, self.ks_min:self.ks_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['SMOIS'].stagger,
                )
            elif input_variable == 'TSLB':
                timeframe = dict(
                    grid=self.input_soil_grid_3d,
                    array=ncfile.variables['TSLB'][0, 0:self.ks_max, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['TSLB'].stagger,
                )
                timeframe['grid'] = self.input_soil_grid_3d_terrain_extended
                neumann_bc_array = timeframe['array'][0:1, :, :]
                timeframe['array'] = numpy.concatenate((neumann_bc_array, timeframe['array']), axis=0)
            elif input_variable == 'HGT':
                timeframe = dict(
                    grid=self.get_grid_2d(stagger=''),
                    array=ncfile.variables['HGT'][0, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['HGT'].stagger,
                )
            elif input_variable == 'PSFC':
                timeframe = dict(
                    grid=self.get_grid_2d(stagger=''),
                    array=ncfile.variables['PSFC'][0, self.j_min:self.j_max, self.i_min:self.i_max],
                    stagger=ncfile.variables['PSFC'].stagger,
                )
                timeframe['grid'] = None
                timeframe['array'] = numpy.mean(timeframe['array'])
            else:
                raise PrometException(f'Unknown input variable "{input_variable}".', id='6Hb3lQ')
        return timeframe

    def get_grid_2d(
            self,
            stagger: str = '',
    ) -> WRFInputGrid2D:
        assert stagger in ['', 'X', 'Y']
        if stagger == '':
            return WRFInputGrid2D.centered(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
            )
        if stagger == 'X':
            return WRFInputGrid2D.staggered_x(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                dx=self.attributes['DX'],
            )
        if stagger == 'Y':
            return WRFInputGrid2D.staggered_y(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                dy=self.attributes['DY'],
            )

    def get_terrain_array_2d(self) -> WRFInputGrid3D:
        return self.hgt[self.j_min:self.j_max, self.i_min:self.i_max]

    def get_terrain_frame(self, time_index):
        return self.get_variable_timeframe(
            input_variable='HGT',
            time_index=time_index,
        )

    def get_terrain_grid_3d(self) -> WRFInputGrid3D:
        return WRFInputGrid3D.centered_2d(
            x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
            y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
            z2d=self.get_terrain_array_2d(),
        )

    def get_soil_grid_3d(
            self,
            terrain_extended: bool = False,
    ) -> WRFInputGrid3D:
        if terrain_extended:
            zsoil = numpy.concatenate((numpy.array([0.0]), self.zsoil[0:self.ks_max]))
        else:
            zsoil = self.zsoil[0:self.ks_max]
        return WRFInputGrid3D.centered_zsoil(
            x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
            y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
            zsoil1d=zsoil,
        )

    def get_grid_3d(
            self,
            stagger: str = '',
            terrain_extended: bool = False,
    ) -> WRFInputGrid3D:
        assert stagger in ['', 'X', 'Y', 'Z']
        if stagger == '':
            return WRFInputGrid3D.centered(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                z3d=self.z_3d[self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                terrain=self.get_terrain_array_2d() if terrain_extended else None,
            )
        if stagger == 'X':
            return WRFInputGrid3D.staggered_x(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                z3d=self.z_3d[self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                dx=self.attributes['DX'],
                terrain=self.get_terrain_array_2d() if terrain_extended else None,
            )
        if stagger == 'Y':
            return WRFInputGrid3D.staggered_y(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                z3d=self.z_3d[self.k_min:self.k_max, self.j_min:self.j_max, self.i_min:self.i_max],
                dy=self.attributes['DY'],
                terrain=self.get_terrain_array_2d() if terrain_extended else None,
            )
        if stagger == 'Z':
            return WRFInputGrid3D.staggered_z(
                x2d=self.x2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                y2d=self.y2d_in_t[self.j_min:self.j_max, self.i_min:self.i_max],
                zw3d=self.zw_3d[self.k_min:self.k_max+1, self.j_min:self.j_max, self.i_min:self.i_max],
            )

    def get_metadata(self) -> dict:
        return dict(
            attributes=self.attributes,
            available_time_slots=self.available_time_slots,
            input_data_shape=self.zw_3d.shape,
            input_data_cropping=[(self.k_min, self.k_max), (self.j_min, self.j_max), (self.i_min, self.i_max)],
        )
