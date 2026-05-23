import os
from datetime import datetime
import netCDF4
import numpy

from promet.api.logging import print_error
from promet.api.logging import print_and_log_section
from promet.api.logging import print_and_log_step


class DynamicDriver:
    """used for creating, populating and modifying dynamic driver file."""

    default_fill_value_float = -9999.0

    allowed_variables = {
        'init_atmosphere_pt': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the potential temperature',
                'units': 'K',
            },
        },
        'init_atmosphere_qv': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the specific humidity',
                'units': 'kg/kg',
            },
        },
        'init_atmosphere_u': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'xu'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'init_atmosphere_v': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'yv', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'init_atmosphere_w': {
            'datatype': 'f4',
            'file_dimensions': ('zw', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'init_atmosphere_no': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'init_atmosphere_no2': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'init_atmosphere_o3': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of ozone',
                'units': 'ppm',
            },
        },
        'init_atmosphere_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of PM10',
                'units': 'Kg/m3',
            },
        },
        'init_atmosphere_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('z', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of PM2.5',
                'units': 'Kg/m3',
            },
        },
        'init_soil_m': {
            'datatype': 'f4',
            'file_dimensions': ('zsoil', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the soil moisture',
                'units': 'm3/m3',
            },
        },
        'init_soil_t': {
            'datatype': 'f4',
            'file_dimensions': ('zsoil', 'y', 'x'),
            'grid_mode': '3d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'initial volume data of the soil temperature',
                'units': 'K',
            },
        },
        'surface_forcing_surface_pressure': {
            'datatype': 'f4',
            'file_dimensions': ('time',),
            'grid_mode': '0d',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 0,
                'long_name': 'surface pressure forcing',
                'units': 'Pa',
            },
        },
        'ls_forcing_left_pt': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary potential temperature',
                'units': 'K',
            },
        },
        'ls_forcing_right_pt': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary potential temperature',
                'units': 'K',
            },
        },
        'ls_forcing_south_pt': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary potential temperature',
                'units': 'K',
            },
        },
        'ls_forcing_north_pt': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary potential temperature',
                'units': 'K',
            },
        },
        'ls_forcing_top_pt': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary potential temperature',
                'units': 'K',
            },
        },
        'ls_forcing_left_qv': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary specific humidity',
                'units': 'kg/kg',
            },
        },
        'ls_forcing_right_qv': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary specific humidity',
                'units': 'kg/kg',
            },
        },
        'ls_forcing_south_qv': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary specific humidity',
                'units': 'kg/kg',
            },
        },
        'ls_forcing_north_qv': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary specific humidity',
                'units': 'kg/kg',
            },
        },
        'ls_forcing_top_qv': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary specific humidity',
                'units': 'kg/kg',
            },
        },
        'ls_forcing_left_u': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_right_u': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_south_u': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'xu'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_north_u': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'xu'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_top_u': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'xu'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary wind velocity component in x direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_left_v': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'yv'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_right_v': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'yv'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_south_v': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_north_v': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_top_v': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'yv', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary wind velocity component in y direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_left_w': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'zw', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_right_w': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'zw', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_south_w': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'zw', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_north_w': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'zw', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_top_w': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary wind velocity component in z direction',
                'units': 'm/s',
            },
        },
        'ls_forcing_left_no': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_right_no': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_south_no': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_north_no': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_top_no': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary nitrogen monoxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_left_no2': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_right_no2': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_south_no2': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_north_no2': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_top_no2': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary nitrogen dioxide',
                'units': 'ppm',
            },
        },
        'ls_forcing_left_o3': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary ozone',
                'units': 'ppm',
            },
        },
        'ls_forcing_right_o3': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary ozone',
                'units': 'ppm',
            },
        },
        'ls_forcing_south_o3': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary ozone',
                'units': 'ppm',
            },
        },
        'ls_forcing_north_o3': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary ozone',
                'units': 'ppm',
            },
        },
        'ls_forcing_top_o3': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary ozone',
                'units': 'ppm',
            },
        },
        'ls_forcing_left_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary PM10',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_right_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary PM10',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_south_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary PM10',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_north_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary PM10',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_top_pm10': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary PM10',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_left_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_left',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the left model boundary PM2.5',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_right_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'y'),
            'grid_mode': '2d_right',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the right model boundary PM2.5',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_south_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_south',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the south model boundary PM2.5',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_north_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'z', 'x'),
            'grid_mode': '2d_north',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the north model boundary PM2.5',
                'units': 'Kg/m3',
            },
        },
        'ls_forcing_top_pm2.5': {
            'datatype': 'f4',
            'file_dimensions': ('time', 'y', 'x'),
            'grid_mode': '2d_top',
            'fill_value': default_fill_value_float,
            'attributes': {
                'lod': 2,
                'long_name': 'large-scale forcing for the top model boundary PM2.5',
                'units': 'Kg/m3',
            },
        },
    }

    allowed_variable_groups = {
        'init_atmosphere_pt': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_pt',
            ]
        },
        'init_atmosphere_qv': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_qv',
            ]
        },
        'init_atmosphere_u': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'xu'),
            'variable_stagger': 'X',
            'variables': [
                'init_atmosphere_u',
            ]
        },
        'init_atmosphere_v': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'yv', 'x'),
            'variable_stagger': 'Y',
            'variables': [
                'init_atmosphere_v',
            ]
        },
        'init_atmosphere_w': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('zw', 'y', 'x'),
            'variable_stagger': 'Z',
            'variables': [
                'init_atmosphere_w',
            ]
        },
        'init_atmosphere_no': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_no',
            ]
        },
        'init_atmosphere_no2': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_no2',
            ]
        },
        'init_atmosphere_o3': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_o3',
            ]
        },
        'init_atmosphere_pm10': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_pm10',
            ]
        },
        'init_atmosphere_pm2.5': {
            'type': 'initial',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'init_atmosphere_pm2.5',
            ]
        },
        'init_soil_m': {
            'type': 'initial',
            'terrain_following': True,
            'grid_dimensions': ('zsoil', 'y', 'x'),
            'variable_stagger': 'Z',
            'variables': [
                'init_soil_m',
            ]
        },
        'init_soil_t': {
            'type': 'initial',
            'terrain_following': True,
            'grid_dimensions': ('zsoil', 'y', 'x'),
            'variable_stagger': 'Z',
            'variables': [
                'init_soil_t',
            ]
        },
        'surface_forcing_surface_pressure': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': tuple(),
            'variable_stagger': '',
            'variables': [
                'surface_forcing_surface_pressure',
            ]
        },
        'ls_forcing_pt': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_pt',
                'ls_forcing_right_pt',
                'ls_forcing_south_pt',
                'ls_forcing_north_pt',
                'ls_forcing_top_pt',
            ]
        },
        'ls_forcing_qv': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_qv',
                'ls_forcing_right_qv',
                'ls_forcing_south_qv',
                'ls_forcing_north_qv',
                'ls_forcing_top_qv',
            ]
        },
        'ls_forcing_u': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'xu'),
            'variable_stagger': 'X',
            'variables': [
                'ls_forcing_left_u',
                'ls_forcing_right_u',
                'ls_forcing_south_u',
                'ls_forcing_north_u',
                'ls_forcing_top_u',
            ]
        },
        'ls_forcing_v': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'yv', 'x'),
            'variable_stagger': 'Y',
            'variables': [
                'ls_forcing_left_v',
                'ls_forcing_right_v',
                'ls_forcing_south_v',
                'ls_forcing_north_v',
                'ls_forcing_top_v',
            ]
        },
        'ls_forcing_w': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('zw', 'y', 'x'),
            'variable_stagger': 'Z',
            'variables': [
                'ls_forcing_left_w',
                'ls_forcing_right_w',
                'ls_forcing_south_w',
                'ls_forcing_north_w',
                'ls_forcing_top_w',
            ]
        },
        'ls_forcing_no': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_no',
                'ls_forcing_right_no',
                'ls_forcing_south_no',
                'ls_forcing_north_no',
                'ls_forcing_top_no',
            ]
        },
        'ls_forcing_no2': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_no2',
                'ls_forcing_right_no2',
                'ls_forcing_south_no2',
                'ls_forcing_north_no2',
                'ls_forcing_top_no2',
            ]
        },
        'ls_forcing_o3': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_o3',
                'ls_forcing_right_o3',
                'ls_forcing_south_o3',
                'ls_forcing_north_o3',
                'ls_forcing_top_o3',
            ]
        },
        'ls_forcing_pm10': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_pm10',
                'ls_forcing_right_pm10',
                'ls_forcing_south_pm10',
                'ls_forcing_north_pm10',
                'ls_forcing_top_pm10',
            ]
        },
        'ls_forcing_pm2.5': {
            'type': 'ls_forcing',
            'terrain_following': False,
            'grid_dimensions': ('z', 'y', 'x'),
            'variable_stagger': '',
            'variables': [
                'ls_forcing_left_pm2.5',
                'ls_forcing_right_pm2.5',
                'ls_forcing_south_pm2.5',
                'ls_forcing_north_pm2.5',
                'ls_forcing_top_pm2.5',
            ]
        },
    }

    @classmethod
    def get_initial_variables(cls):
        return [var_name for var_name, var_dict in cls.allowed_variables.items() if 'time' not in var_dict['file_dimensions']]

    @classmethod
    def get_forcing_variables(cls):
        return [var_name for var_name, var_dict in cls.allowed_variables.items() if 'time' in var_dict['file_dimensions']]

    @classmethod
    def get_initial_variable_groups(cls):
        return [group_name for group_name, group_dict in cls.allowed_variable_groups.items() if group_dict['type'] == 'initial']

    @classmethod
    def get_forcing_variable_groups(cls):
        return [group_name for group_name, group_dict in cls.allowed_variable_groups.items() if group_dict['type'] == 'ls_forcing']

    @classmethod
    def get_group_type(cls, group_name):
        return cls.allowed_variable_groups[group_name]['type']

    @classmethod
    def is_group_terrain_following(cls, group_name):
        return cls.allowed_variable_groups[group_name]['terrain_following']

    @classmethod
    def get_group_grid_dimensions(cls, group_name):
        return cls.allowed_variable_groups[group_name]['grid_dimensions']

    @classmethod
    def get_group_stagger(cls, group_name):
        return cls.allowed_variable_groups[group_name]['variable_stagger']

    @classmethod
    def get_group_variables(cls, group_name):
        return cls.allowed_variable_groups[group_name]['variables']

    def __init__(
            self,
            filepath: str,
    ):
        self.filepath = os.path.normpath(os.path.expandvars(os.path.expanduser(filepath)))
        self.filemode = 'r'
        self.write_if_exists = False
        self.already_existing_variables = []
        self.already_existing_dimensions = []

    def initialize(
            self,
            append: bool = False,
            overwrite: bool = False,
    ) -> bool:
        exists = os.path.exists(self.filepath)
        if not exists:
            self.filemode = 'x'
            self.write_if_exists = True
        elif exists and append and overwrite:
            self.filemode = 'r+'
            self.write_if_exists = True
        elif exists and not append and overwrite:
            self.filemode = 'w'
            self.write_if_exists = True
        elif exists and append and not overwrite:
            self.filemode = 'r+'
            self.write_if_exists = False
        elif exists and not append and not overwrite:
            print_error(f'File "{self.filepath}" already exists. You can use --append or --overwrite to specify what to do.')
            return True
        if exists:
            with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
                self.already_existing_variables = ncfile.variables
                self.already_existing_dimensions = ncfile.dimensions
        return False

    def set_global_attributes(
            self,
            origin_lat: float,
            origin_lon: float,
            origin_time: datetime,
            origin_x: float,
            origin_y: float,
            origin_z: float,
            rotation_angle: float = 0.0,
            title: str = 'dynamic driver created by promet',
            author: str = 'created with promet (by pecanode.com)',
            verbose: bool = False,
    ) -> bool:
        print_and_log_section('Populating dynamic driver file with global attributes...')
        attributes = dict(
            origin_lat=origin_lat,
            origin_lon=origin_lon,
            origin_time=origin_time,
            origin_x=origin_x,
            origin_y=origin_y,
            origin_z=origin_z,
            rotation_angle=rotation_angle,
            title=title,
            author=author,
            Conventions='CF-1.7',
            creation_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S +00"),
        )
        with netCDF4.Dataset(self.filepath, mode=self.filemode) as ncfile:
            for key, value in attributes.items():
                if self.write_if_exists or key not in ncfile.ncattrs():
                    ncfile.setncattr(key, value)
                    if verbose:
                        print_and_log_step(f'Writing global attribute {key} = {value}')
        return False

    def set_dimensions(
            self,
            x: dict,
            xu: dict,
            y: dict,
            yv: dict,
            z: dict,
            zw: dict,
            zsoil: dict,
            verbose: bool = False,
    ) -> bool:
        print_and_log_section('Populating dynamic driver file with dimensions...')
        dimensions = dict(
            x=x,
            xu=xu,
            y=y,
            yv=yv,
            z=z,
            zw=zw,
            zsoil=zsoil,
        )
        with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
            for key, dimension in dimensions.items():
                if self.write_if_exists or key not in self.already_existing_dimensions:
                    if key not in ncfile.dimensions:
                        ncfile.createDimension(dimname=key, size=dimension["size"])
                    else:
                        assert ncfile.dimensions[key].size == dimension["size"]
                    if key not in ncfile.variables:
                        nc_xu = ncfile.createVariable(key, 'f4', (key,))
                    else:
                        nc_xu = ncfile.variables[key]
                    nc_xu[:] = dimension["values"]
                    nc_xu.long_name = key
                    nc_xu.units = "m"
                    if verbose:
                        print_and_log_step(f'Writing dimension {key} with size={dimension["size"]}')
        return False

    def create_time(self, time_slots):
        with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
            if self.write_if_exists or 'time' not in self.already_existing_dimensions:
                if 'time' not in ncfile.dimensions:
                    ncfile.createDimension(dimname='time', size=len(time_slots))
                else:
                    assert ncfile.dimensions['time'].size == len(time_slots)
                if 'time' not in ncfile.variables:
                    nc_xu = ncfile.createVariable(
                        varname='time',
                        datatype='f4',
                        dimensions=('time',),
                        fill_value=self.default_fill_value_float,
                    )
                else:
                    nc_xu = ncfile.variables['time']
                nc_xu[:] = numpy.array([t['start_time_delta'] for n, t in time_slots.items()])
                nc_xu.long_name = 'time'
                nc_xu.units = "s"

    def create_variable(
            self,
            var_name: str,
            verbose: bool = False,
    ):
        assert var_name in self.allowed_variables
        with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
            if self.write_if_exists or var_name not in self.already_existing_variables:
                if var_name not in ncfile.variables:
                    nc_var = ncfile.createVariable(
                        varname=var_name,
                        datatype=self.allowed_variables[var_name]['datatype'],
                        dimensions=self.allowed_variables[var_name]['file_dimensions'],
                        fill_value=self.allowed_variables[var_name]['fill_value'],
                    )
                    if verbose:
                        print_and_log_step(
                            f'Creating variable {var_name:32} with dimensions={self.allowed_variables[var_name]["file_dimensions"]}'
                        )
                else:
                    nc_var = ncfile.variables[var_name]
                    nc_var.datatype = self.allowed_variables[var_name]['datatype']
                    nc_var.fill_value = self.allowed_variables[var_name]['fill_value']
                    if verbose:
                        print_and_log_step(
                            f'Found existing variable {var_name:32} with dimensions={self.allowed_variables[var_name]["file_dimensions"]}'
                        )
                nc_var.setncatts(self.allowed_variables[var_name]['attributes'])

    def set_init_variable(
            self,
            var_name: str,
            array,
            verbose: bool = False,
    ):
        with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
            assert var_name in ncfile.variables
            if self.write_if_exists or var_name not in self.already_existing_variables:
                if verbose:
                    print_and_log_step(
                        f'Writing variable {var_name}'
                    )
                nc_var = ncfile.variables[var_name]
                nc_var[:] = array

    def set_forcing_variable(
            self,
            var_name: str,
            time_index: int,
            array,
            verbose: bool = False,
    ):
        with netCDF4.Dataset(self.filepath, mode='r+') as ncfile:
            assert var_name in ncfile.variables
            if self.write_if_exists or var_name not in self.already_existing_variables:
                if verbose:
                    print_and_log_step(
                        f'Writing time_index "{time_index}" of variable {var_name}'
                    )
                nc_var = ncfile.variables[var_name]
                nc_var[time_index] = array
