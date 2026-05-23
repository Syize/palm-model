import os
import yaml
from datetime import datetime

from promet.api.logging import PrometException
from promet.api.logging import print_error

class Config:

    default_content = {
        'plugin': ('icon', True, ['icon', 'wrf']),
        'input_data': ('input_data/', True),
        'input_crs': ('EPSG:4326', True),
        'palm_crs': ('EPSG:25832', True),
        'palm_variable_groups': ([
            'init_atmosphere_u',
            'init_atmosphere_v',
            'init_atmosphere_w',
        ], True),
        'start_time': (datetime.strptime('2023-11-02 15:30:00 +00', '%Y-%m-%d %H:%M:%S +00'), True),
        'end_time': (datetime.strptime('2023-11-02 15:30:00 +00', '%Y-%m-%d %H:%M:%S +00'), True),
        'zsoil': ([0.01, 0.25, 0.7, 1.5], False),
        'top_adjustment_index': (-1, False),
    }

    def __init__(
            self,
            filepath: str,
            verbose: bool = False,
    ):
        self.filepath = os.path.normpath(os.path.expandvars(os.path.expanduser(filepath)))
        self.verbose = verbose
        self._content = dict()

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._content})'

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, key):
        return self._content.get(key)

    @property
    def content(self):
        return self._content

    def load(self) -> bool:
        # Open the YAML config file for reading
        try:
            with open(self.filepath, 'r') as config_file:
                try:
                    # Load the YAML data from the file
                    input_dict = yaml.safe_load(config_file)
                except yaml.YAMLError as e:
                    print_error(f'Failed to parse YAML: {e}')
                    e.id = 'TXBIPg'
                    raise e
        except FileNotFoundError as e:
            print_error(f'Unable to find config file: {self.filepath}')
            e.id = 'nRRIHw'
            raise e
        if input_dict is None:
            print_error('Config file is empty.')
            input_dict = dict()
        self._content, content_errors = self.set_defaults_and_check_types(
            input_dict=input_dict,
            default_dict=self.default_content,
        )
        if content_errors:
            print_error(f'Errors in config detected:')
            for content_error in content_errors:
                print_error(f'   {content_error}')
            raise PrometException(f'Found {len(content_errors)} errors in config.', id='yKSsXg')
        return False

    def set_defaults_and_check_types(
            self,
            input_dict: dict,
            default_dict: dict,
            location: str = '',
    ) -> (dict, list):
        result = {}
        errors = []
        for key, default_info in default_dict.items():
            current_location = f"{location}/{key}" if location else key  # Build the full location path
            if isinstance(default_info, tuple):
                if len(default_info) == 2:
                    default_value, mandatory = default_info
                    allowed = []
                elif len(default_info) == 3:
                    default_value, mandatory, allowed = default_info
                else:
                    raise PrometException(
                        f'Config variable "default_content" at location "{current_location}" '
                        f'has unknown tuple length (only 2 or 3 are allowed)', id='yROkXA',
                    )
            else:
                default_value, mandatory, allowed = default_info, False, []
            if key in input_dict:
                if isinstance(input_dict[key], dict):
                    if not isinstance(default_value, dict):
                        errors.append(
                            f'Value at location "{current_location}" is a dictionary, '
                            f'but should be of type "{type(default_value).__name__}".'
                        )
                    else:
                        subresult, suberrors = self.set_defaults_and_check_types(
                            input_dict[key], default_value, current_location
                        )
                        result[key] = subresult
                        errors.extend(suberrors)
                elif isinstance(input_dict[key], type(default_value)):
                    result[key] = input_dict[key]
                else:
                    errors.append(
                        f'Value at location "{current_location}" should be of type "{type(default_value).__name__}", '
                        f'not "{type(input_dict[key]).__name__}".'
                    )
                assert isinstance(allowed, list)
                if len(allowed) > 0:
                    assert all(isinstance(v, type(default_value)) for v in allowed)
                    if input_dict[key] not in allowed:
                        errors.append(
                            f'Value at location "{current_location}" has unknown value. '
                            f'The following values are allowed: {allowed}'
                        )
            elif mandatory:
                errors.append(
                    f"Key at location '{os.path.join(location, key)}' is mandatory but is missing.")
            else:
                result[key] = default_value
        for key, value in input_dict.items():
            if key not in default_dict:
                errors.append(f"Key at location '{os.path.join(location, key)}' is not allowed.")
        return result, errors

