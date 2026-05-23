from promet.api.plugin import Plugin

from .input_data import WRFInputData


class WRFPlugin(Plugin):

    input_data_class = WRFInputData

    variable_group_mapping = {
        'init_atmosphere_pt': {
            'input_variable': 'T',
        },
        'init_atmosphere_qv': {
            'input_variable': 'QVAPOR',
        },
        'init_atmosphere_u': {
            'input_variable': 'U',
        },
        'init_atmosphere_v': {
            'input_variable': 'V',
        },
        'init_atmosphere_w': {
            'input_variable': 'W',
        },
        'init_atmosphere_no': {
            'input_variable': 'no',
        },
        'init_atmosphere_no2': {
            'input_variable': 'no2',
        },
        'init_atmosphere_o3': {
            'input_variable': 'o3',
        },
        'init_atmosphere_pm10': {
            'input_variable': 'PM10',
        },
        'init_atmosphere_pm2.5': {
            'input_variable': 'PM2_5_DRY',
        },
        'init_soil_m': {
            'input_variable': 'SMOIS',
        },
        'init_soil_t': {
            'input_variable': 'TSLB',
        },
        'surface_forcing_surface_pressure': {
            'input_variable': 'PSFC',
        },
        'ls_forcing_pt': {
            'input_variable': 'T',
        },
        'ls_forcing_qv': {
            'input_variable': 'QVAPOR',
        },
        'ls_forcing_u': {
            'input_variable': 'U',
        },
        'ls_forcing_v': {
            'input_variable': 'V',
        },
        'ls_forcing_w': {
            'input_variable': 'W',
        },
        'ls_forcing_no': {
            'input_variable': 'no',
        },
        'ls_forcing_no2': {
            'input_variable': 'no2',
        },
        'ls_forcing_o3': {
            'input_variable': 'o3',
        },
        'ls_forcing_pm10': {
            'input_variable': 'PM10',
        },
        'ls_forcing_pm2.5': {
            'input_variable': 'PM2_5_DRY',
        },
    }
