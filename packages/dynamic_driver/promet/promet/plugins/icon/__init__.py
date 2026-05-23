from promet.api.plugin import Plugin

from .input_data import ICONInputData


class ICONPlugin(Plugin):

    input_data_class = ICONInputData

    variable_group_mapping = {
        'init_atmosphere_pt': {
            'input_variable': 'T',
        },
        'init_atmosphere_qv': {
            'input_variable': 'QV',
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
        'init_atmosphere_o3': {
            'input_variable': 'ART_O3',
        },
        'init_soil_m': {
            'input_variable': 'W_SO',
        },
        'init_soil_t': {
            'input_variable': 'T_SO',
        },
        'surface_forcing_surface_pressure': {
            'input_variable': 'PS',
        },
        'ls_forcing_pt': {
            'input_variable': 'T',
        },
        'ls_forcing_qv': {
            'input_variable': 'QV',
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
        'ls_forcing_o3': {
            'input_variable': 'ART_O3',
        },
    }
