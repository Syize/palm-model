import os
import textwrap
import re
import yaml
import termcolor
import jinja2


project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
docs_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

file_link_url = f'file://{project_dir}'

def add_palm_macros(env):
    env.macro(include_palm_driver_dimensions)
    env.macro(include_palm_driver_global_attributes)
    env.macro(include_palm_driver_variables)
    env.macro(include_palm_namelist)
    env.macro(include_palm_logging_ids)
    env.macro(include_palm_output_quantities)
    env.macro(link_palm_repo_file)


def print_message_to_terminal(message_string, loglevel):
    message_wrapped = textwrap.fill(message_string, width=100, initial_indent=' ' * 9 + '-  ', subsequent_indent=' ' * 12)
    if loglevel == 'info':
        message = termcolor.colored('INFO:', 'green') + message_wrapped[5:]
    elif loglevel == 'warning':
        message = termcolor.colored('WARNING:', 'yellow') + message_wrapped[8:]
    elif loglevel == 'error':
        message = termcolor.colored('ERROR:', 'red') + message_wrapped[6:]
    else:
        raise ValueError(f'Unknown loglevel: "{loglevel}"')
    print(message)


def print_yaml_parameter_warning_to_terminal(path, parameter_type_name, parameter_name, text):
    print_message_to_terminal(
        message_string='reading yaml database "{}" with warning {}. {}'.format(
            termcolor.colored(path, 'magenta'),
            f'at {parameter_type_name} "{termcolor.colored(parameter_name, "cyan")}"',
            text,
        ),
        loglevel='warning'
    )


def print_yaml_parameter_error_to_terminal(path, parameter_name, text):
    print_message_to_terminal(
        message_string='reading yaml database "{}" fails {}. {}'.format(
            termcolor.colored(path, 'magenta'),
            f'at parameter "{termcolor.colored(parameter_name, "cyan")}"',
            text,
        ),
        loglevel='error'
    )


def print_yaml_field_error_to_terminal(
        database_path: str,
        field_identifier: str,
        message: str,
) -> None:
    print_message_to_terminal(
        message_string='reading yaml database "{}" fails {}. {}'.format(
            termcolor.colored(database_path, 'magenta'),
            f'at "{termcolor.colored(field_identifier, "cyan")}"',
            message,
        ),
        loglevel='error'
    )


def validate_yaml_field_generic(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        field_name,
        field_type,
        mandatory: bool = True,
        contains_list: bool = False,
        allowed_values: list = None
) -> int:
    violation_counter = 0
    if not field_name in content_dict and mandatory:
        print_yaml_field_error_to_terminal(
            database_path=database_path,
            field_identifier=field_identifier,
            message=f'The field "{field_name}" is mandatory but missing',
        )
        violation_counter += 1
    if field_name in content_dict:
        if contains_list:
            if not isinstance(content_dict[field_name], list):
                print_yaml_field_error_to_terminal(
                    database_path=database_path,
                    field_identifier=field_identifier,
                    message=f'The field "{field_name}" must have a value of type "list"',
                )
                violation_counter += 1
            else:
                if not all(isinstance(v, field_type) for v in content_dict[field_name]):
                    print_yaml_field_error_to_terminal(
                        database_path=database_path,
                        field_identifier=field_identifier,
                        message=f'The field "{field_name}" must have a value of type "list" which must contain items of type "{str(field_type.__name__)}"',
                    )
                    violation_counter += 1
                else:
                    if allowed_values:
                        local_regex = '|'.join([f'^{v}$' for v in allowed_values])
                        if not all(bool(re.match(local_regex, content_dict[field_name])) for v in content_dict[field_name]):
                            print_yaml_field_error_to_terminal(
                                database_path=database_path,
                                field_identifier=field_identifier,
                                message=f'The field "{field_name}" must have a list of values where each is one of "{allowed_values}"',
                            )
                            violation_counter += 1
        else:
            try:
                isinstance(content_dict[field_name], field_type)
            except TypeError as e:
                print(field_type)
            if not isinstance(content_dict[field_name], field_type):
                print_yaml_field_error_to_terminal(
                    database_path=database_path,
                    field_identifier=field_identifier,
                    message=f'The field "{field_name}" must have a value of type "{str(field_type.__name__)}"',
                )
                violation_counter += 1
            else:
                if allowed_values:
                    local_regex = '|'.join([f'^{v}$' for v in allowed_values])
                    if not bool(re.match(local_regex, content_dict[field_name])):
                        print_yaml_field_error_to_terminal(
                            database_path=database_path,
                            field_identifier=field_identifier,
                            message=f'The field "{field_name}" must have a value that is one of "{allowed_values}"',
                        )
                        violation_counter += 1

    return violation_counter


def validate_yaml_field_mutually_exclusive_generic(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        field_names: list,
        field_types: list,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    present_fields = [key for key in field_names if key in content_dict]
    if len(present_fields) == 0:
        print_yaml_field_error_to_terminal(
            database_path=database_path,
            field_identifier=field_identifier,
            message=f'At least one field from {field_names} is mandatory and must be present.',
        )
        violation_counter += 1
    elif len(present_fields) > 1:
        print_yaml_field_error_to_terminal(
            database_path=database_path,
            field_identifier=field_identifier,
            message=f'Only one field from {field_names} is allowed, but found multiple: {present_fields}',
        )
        violation_counter += 1
    else:
        violation_counter += validate_yaml_field_generic(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=content_dict,
            field_name=present_fields[0],
            field_type=field_types[field_names.index(present_fields[0])],
            mandatory=mandatory,
        )
    return violation_counter


def validate_yaml_field_mandatory(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='mandatory',
        field_type=bool,
        mandatory=mandatory,
    )
    return violation_counter


def validate_yaml_field_type(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
        allowed_values: list = None
) -> int:
    violation_counter = 0
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='type',
        field_type=str,
        mandatory=mandatory,
        allowed_values=allowed_values,
    )
    return violation_counter


def validate_yaml_field_modules(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
        allowed_values: list = None
) -> int:
    violation_counter = 0
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='modules',
        field_type=list,
        mandatory=mandatory,
    )

    if allowed_values:
        local_regex = '|'.join([f'^{v}$' for v in allowed_values])
        for value in content_dict['modules']:
            if not bool(re.match(local_regex, value)):
                print_yaml_field_error_to_terminal(
                    database_path=database_path,
                    field_identifier=field_identifier,
                    message=f'The field "modules" must be a list of values that each are one of "{allowed_values}"',
                )
                violation_counter += 1
    return violation_counter


def validate_yaml_field_size(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    if 'size' in content_dict:
        if not isinstance(content_dict['size'], dict):
            content_dict['size'] = dict(
                value=content_dict['size'],
            )
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='size',
        field_type=dict,
        mandatory=mandatory,
    )
    if 'size' in content_dict:
        violation_counter += validate_yaml_field_mutually_exclusive_generic(
            database_path=database_path,
            field_identifier=f'{field_identifier}.size',
            content_dict=content_dict['size'],
            field_names=['value', 'value_of', 'depends_on'],
            field_types=[int ,str , str],
            mandatory=True,
        )
        if 'value' in content_dict['size']:
            pass
        elif 'value_of' in content_dict['size']:
            content_dict['size']['value'] = 'Value of {}'.format(content_dict['size']['value_of'])
        elif 'depends_on' in content_dict['default']:
            content_dict['size']['value'] = 'Depends on {}'.format(content_dict['size']['depends_on'])
    return violation_counter


def validate_yaml_field_shape(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    if 'shape' in content_dict:
        if isinstance(content_dict['shape'], str):
            if not bool(re.fullmatch(r'\(\[*(?:[a-z_][a-z0-9_]*|\d+)\]*(?:,\[*(?:[a-z_][a-z0-9_]*|\d+)\]*)*\)', content_dict['shape'])):
                print_yaml_field_error_to_terminal(
                    database_path=database_path,
                    field_identifier=field_identifier,
                    message='The field "shape" must have a value of type "tuple" which must contain items of type "int" or "str" separated by a "," and no spaces. '
                            'An example for a shape would be "([zsoil],y,x)" or "(3,y,x)"',
                )
                violation_counter += 1
            else:
                shape_list = []
                for n, v in enumerate(content_dict['shape'][1:-1].split(',')):
                    if bool(re.fullmatch(r'\[*[a-z_][a-z0-9_]*\]*', v)):
                        shape_list.append(v)
                    elif bool(re.fullmatch(r'\d+', v)):
                        shape_list.append(int(v))
                    else:
                        print_yaml_field_error_to_terminal(
                            database_path=database_path,
                            field_identifier=field_identifier,
                            message='The field "shape" must have a value of type "tuple" which must contain items of type "int" or "str" separated by a "," and no spaces. '
                                    f'At index "{n}" the tuple contains the invalid value "{v}"',
                        )
                        violation_counter += 1
                content_dict['shape'] = tuple(shape_list)
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='shape',
        field_type=tuple,
        mandatory=mandatory,
    )
    return violation_counter


def validate_yaml_field_si_unit(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    if 'si-unit' in content_dict:
        if isinstance(content_dict['si-unit'], int) and content_dict['si-unit'] == 1:
            content_dict['si-unit'] = str(content_dict['si-unit'])
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='si-unit',
        field_type=str,
        mandatory=mandatory,
    )
    return violation_counter


def validate_yaml_field_default(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        value_type: type,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    if 'default' in content_dict:
        if not isinstance(content_dict['default'], dict):
            content_dict['default'] = dict(
                value=content_dict['default'],
            )
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='default',
        field_type=dict,
        mandatory=mandatory,
    )
    if 'default' in content_dict:
        if 'value' in content_dict['default']:
            if content_dict['default']['value'] is None:
                content_dict['default']['value'] = 'undefined'
                value_type = str
        violation_counter += validate_yaml_field_mutually_exclusive_generic(
            database_path=database_path,
            field_identifier=f'{field_identifier}.default',
            content_dict=content_dict['default'],
            field_names=['value', 'value_of', 'depends_on'],
            field_types=[value_type ,str , str],
            mandatory=True,
        )
        if 'value' in content_dict['default']:
            pass
        elif 'value_of' in content_dict['default']:
            content_dict['default']['value'] = 'Value of {}'.format(content_dict['default']['value_of'])
        elif 'depends_on' in content_dict['default']:
            content_dict['default']['value'] = 'Depends on {}'.format(content_dict['default']['depends_on'])
    return violation_counter


def validate_yaml_field_description(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='description',
        field_type=dict,
        mandatory=mandatory,
    )
    if 'description' in content_dict:
        violation_counter += validate_yaml_field_generic(
            database_path=database_path,
            field_identifier=f'{field_identifier}.description',
            content_dict=content_dict['description'],
            field_name='short',
            field_type=str,
            mandatory=True,
        )
        violation_counter += validate_yaml_field_generic(
            database_path=database_path,
            field_identifier=f'{field_identifier}.description',
            content_dict=content_dict['description'],
            field_name='long',
            field_type=str,
            mandatory=False,
        )
    return violation_counter


def validate_yaml_field_allowed_values(
        database_path: str,
        field_identifier: str,
        content_dict: dict,
        value_type: type,
        mandatory: bool = True,
) -> int:
    violation_counter = 0
    violation_counter += validate_yaml_field_generic(
        database_path=database_path,
        field_identifier=field_identifier,
        content_dict=content_dict,
        field_name='allowed_values',
        field_type=list,
        mandatory=mandatory,
    )
    if 'allowed_values' in content_dict:
        for n, allowed_value_dict in enumerate(content_dict['allowed_values']):
            violation_counter += validate_yaml_field_mutually_exclusive_generic(
                database_path=database_path,
                field_identifier=f'{field_identifier}.allowed_values[{n}]',
                content_dict=allowed_value_dict,
                field_names=['value', 'value_of', 'depends_on'],
                field_types=[value_type ,str , str],
                mandatory=True,
            )
            if 'value' in allowed_value_dict:
                pass
            elif 'value_of' in allowed_value_dict:
                allowed_value_dict['value'] = 'Value of {}'.format(allowed_value_dict['value_of'])
            elif 'depends_on' in allowed_value_dict:
                allowed_value_dict['value'] = 'Depends on {}'.format(allowed_value_dict['depends_on'])
            violation_counter += validate_yaml_field_generic(
                database_path=database_path,
                field_identifier=f'{field_identifier}.allowed_values[{n}]',
                content_dict=allowed_value_dict,
                field_name='description',
                field_type=str,
                mandatory=True,
            )
    return violation_counter






def validate_yaml_namelist_database(data_path, content_dict):
    """
    <parameter-name>
      category: ''
      type: ''
      shape: ()  # optional
      default:
        value:  # optional
        depends_on:  # optional
        value_of:  # optional
      si-unit:  # optional
      description:
        short: ''
        long: ''  # optional
      allowed_values:  # optional
        - value: 1.0
          description: ''
    """
    valid = True
    for parameter, parameter_dict in content_dict.items():

        if 'category' in parameter_dict and not 'categories' in parameter_dict:
            parameter_dict['categories'] = parameter_dict['category']
            #print_yaml_parameter_warning_to_terminal(
            #    data_path,
            #    parameter,
            #    'Deprecated field "category" is used. Please use field "categories" instead',
            #)
        if not 'categories' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "categories" is mandatory but missing',
            )
        if isinstance(parameter_dict['categories'], str):
            parameter_dict['categories'] = [parameter_dict['categories']]
        if not isinstance(parameter_dict['categories'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "categories" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in parameter_dict['categories']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "categories:" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False

        parameter_type_class = None
        if not 'type' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "type" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(parameter_dict['type'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "type" must have a value of type "str"',
                )
                valid = False
            if not bool(re.fullmatch(r'C(?:\*\d+)?$|^I$|^L$|^R$|^D', parameter_dict['type'])):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "type" must have a value that starts with either C, C*<int>, I, L, R or D',
                )
                valid = False
            if parameter_dict['type'].startswith('C'):
                parameter_type_class = str
            if parameter_dict['type'].startswith('I'):
                parameter_type_class = int
            if parameter_dict['type'].startswith('L'):
                parameter_type_class = bool
            if parameter_dict['type'].startswith('R'):
                parameter_type_class = float

        if 'shape' in parameter_dict:
            if not isinstance(parameter_dict['shape'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "shape" must have a value of type "tuple". '
                    'An example for a shape of a 2D array would be "(10,50)"',
                )
                valid = False
            else:
                if not bool(re.fullmatch(r'\((\d+,)*(\d+)\)', parameter_dict['shape'])):
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The optional field "shape" must have a value of type "tuple" which must contain items '
                        'of type "int". An example for a shape of a 2D array would be "(10,50)"',
                    )
                    valid = False
                else:
                    parameter_dict['shape'] = tuple([int(v) for v in parameter_dict['shape'][1:-1].split(',')])
                    if not all(isinstance(v, int) for v in parameter_dict['shape']):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The optional field "shape" must have a value of type "tuple" '
                            'which must contain items of type "int"',
                        )
                        valid = False

        if not 'default' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "default" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(parameter_dict['default'], dict):
                parameter_dict['default'] = dict(
                    value=parameter_dict['default'],
                )
            if not any(map(lambda k: k in parameter_dict['default'].keys(), ['value', 'value_of', 'depends_on'])):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "default" must contain a key '
                    'that is either "value", "value_of" or "depends_on"',
                )
                valid = False

            if not isinstance(parameter_dict['default'], dict):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "default" must have a value of type "dict"',
                )
                valid = False
            if 'value' in parameter_dict['default']:
                if parameter_dict['default']['value'] is not None:
                    if parameter_type_class == bool:
                        if parameter_dict['default']['value'] in ['.T.', '.TRUE.']:
                            parameter_dict['default']['value'] = True
                        elif parameter_dict['default']['value'] in ['.F.', '.FALSE.']:
                            parameter_dict['default']['value'] = False
                        else:
                            print_yaml_parameter_error_to_terminal(
                                data_path,
                                parameter,
                                'The mandatory field "default.value" must have a value that represents '
                                'a Fortran LOGICAL and can be either .T., .TRUE., .F. or .FALSE.',
                            )
                            valid = False
                    if 'shape' in parameter_dict and isinstance(parameter_dict['default']['value'], list):
                        if not all(isinstance(v, parameter_type_class) for v in parameter_dict['default']['value']):
                            print_yaml_parameter_error_to_terminal(
                                data_path,
                                parameter,
                                'The mandatory field "default.value" must have a value of type "list" '
                                'which must contain items of type "{}"'.format(str(parameter_type_class.__name__)),
                            )
                            valid = False
                    else:
                        if not isinstance(parameter_dict['default']['value'], parameter_type_class):
                            print_yaml_parameter_error_to_terminal(
                                data_path,
                                parameter,
                                'The mandatory field "default.value" must have a value '
                                'of type "{}"'.format(str(parameter_type_class.__name__)),
                            )
                            valid = False
                    if parameter_type_class == bool:
                        parameter_dict['default']['value'] = '.TRUE.' if parameter_dict['default']['value'] else '.FALSE.'
            elif 'value_of' in parameter_dict['default']:
                parameter_dict['default']['value'] = 'Value of {}'.format(parameter_dict['default']['value_of'])
            elif 'depends_on' in parameter_dict['default']:
                parameter_dict['default']['value'] = 'Depends on {}'.format(parameter_dict['default']['depends_on'])

        if 'si-unit' in parameter_dict:
            if not isinstance(parameter_dict['si-unit'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "si-unit" must have a value of type "str"',
                )
                valid = False

        if not 'description' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "description" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(parameter_dict['description'], dict):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "description" must have a value of type "dict"',
                )
                valid = False

            if not 'short' in parameter_dict['description']:
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The field "description.short" is mandatory but missing',
                )
                valid = False
            if not isinstance(parameter_dict['description']['short'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The mandatory field "description.short" must have a value of type "str"',
                )
                valid = False

            if 'long' in parameter_dict['description']:
                if not isinstance(parameter_dict['description']['long'], str):
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The optional field "description.long" must have a value of type "str"',
                    )
                    valid = False

        if 'allowed_values' in parameter_dict:
            if not isinstance(parameter_dict['allowed_values'], list):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "allowed_values" must have a value of type "list"',
                )
                valid = False
            if not all(isinstance(v, dict) for v in parameter_dict['allowed_values']):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "allowed_values" must have a value of type "list" '
                    'which must contain items of type "dict"',
                )
                valid = False

            for allowed_value in parameter_dict['allowed_values']:
                if not 'value' in allowed_value:
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The field "allowed_values[].value" is mandatory but missing',
                    )
                    valid = False
                if 'shape' in parameter_dict and isinstance(allowed_value['value'], list):
                    if not all(isinstance(v, parameter_type_class) for v in allowed_value['value']):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "allowed_values[].value" has a value of type "list" '
                            'which must contain items of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False
                else:
                    if not isinstance(allowed_value['value'], parameter_type_class):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "allowed_values[].value" must have a value '
                            'of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False

                if not 'description' in allowed_value:
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The field "allowed_values[].description" is mandatory but missing',
                    )
                    valid = False
                if not isinstance(allowed_value['description'], str):
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The mandatory field "allowed_values[].description" must have a value of type "str"',
                    )
                    valid = False
    return valid


def render_namelist_to_markdown_as_table(namelist, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Parameter', 'Default','Description')
    output_str += format_str_head.format('-', '-', '-', '-', '-')
    for parameter, parameter_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--{2})'.format(link_path, namelist, parameter) if link_table else parameter,
            parameter_dict['default']['value'] if parameter_dict['default']['value'] is not None else 'undefined',
            parameter_dict['description']['short'],
        )

    output_str += '\n\n'
    output_str += '*[I]: Integer\n'
    output_str += '*[I1]: Integer 8-Bit\n'
    output_str += '*[I2]: Integer 16-Bit\n'
    output_str += '*[I4]: Integer 32-Bit\n'
    output_str += '*[I8]: Integer 64-Bit\n'
    output_str += '*[R]: Real\n'
    output_str += '*[R1]: Real 8-Bit\n'
    output_str += '*[R2]: Real 16-Bit\n'
    output_str += '*[R4]: Real 32-Bit\n'
    output_str += '*[R8]: Real 64-Bit\n'
    output_str += '*[C]: Character\n'
    output_str += '*[L]: Logical\n'
    output_str += '*[D]: Derived Data-Type\n'
    output_str += '*[undefined]: This parameter has no usable default value and probably needs to be set by the user!\n'
    return output_str


def render_namelist_to_markdown(namelist, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = '<br>\n'
    for parameter, parameter_dict in content_dict.items():
        output_str += '#' * heading_level + parameter + ' {#' + '{0}--{1}'.format(namelist, parameter) + '}\n\n'
        output_str += list_format_str.format(
            'Fortran Type', '{} {}'.format(
                parameter_dict['type'],
                '({})'.format(','.join([str(x) for x in parameter_dict['shape']])) if 'shape' in parameter_dict else '',
            )
        )
        output_str += list_format_str.format('Default', '*' + str(parameter_dict['default']['value']) + '*' if
        parameter_dict['default']['value'] is not None else '*undefined*')
        if 'si-unit' in parameter_dict:  # and re.search('^[I,R].*', parameter_dict['type']):
            output_str += list_format_str.format('SI-Unit', parameter_dict['si-unit'])
        output_str += '\n{}\n\n'.format(
            textwrap.indent(parameter_dict['description']['short'], ' ' * 4),
        )
        if 'long' in parameter_dict['description']:
            if parameter_dict['description']['long']:
                output_str += '{}\n\n'.format(
                    textwrap.indent(parameter_dict['description']['long'], ' ' * 4),
                )
        if 'allowed_values' in parameter_dict:
            output_str += '\n{}\n\n'.format(
                textwrap.indent('Currently {} choices are available:'.format(len(parameter_dict['allowed_values'])),
                                ' ' * 4),
            )
            for allowed_value in parameter_dict['allowed_values']:
                output_str += '    - *{}*\n\n{}\n\n'.format(
                    allowed_value['value'],
                    textwrap.indent(allowed_value['description'], ' ' * 8),
                )
        output_str += '\n<br>\n'
    output_str += '\n\n'
    output_str += '*[I]: Integer\n'
    output_str += '*[I1]: Integer 8-Bit\n'
    output_str += '*[I2]: Integer 16-Bit\n'
    output_str += '*[I4]: Integer 32-Bit\n'
    output_str += '*[I8]: Integer 64-Bit\n'
    output_str += '*[R]: Real\n'
    output_str += '*[R1]: Real 8-Bit\n'
    output_str += '*[R2]: Real 16-Bit\n'
    output_str += '*[R4]: Real 32-Bit\n'
    output_str += '*[R8]: Real 64-Bit\n'
    output_str += '*[C]: Character\n'
    output_str += '*[L]: Logical\n'
    output_str += '*[D]: Derived Data-Type\n'
    output_str += '*[undefined]: This parameter has no usable default value and probably needs to be set by the user!\n'
    return output_str


def include_palm_namelist(data_path, categories=['all'], as_table=False, link_table=True, link_path='', heading_level=3):
    call_string = 'include_palm_namelist(\'{}\', categories={}, as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, categories, as_table, link_table, link_path, heading_level,
    )

    namelist = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(categories, str):
        categories = [categories]
    assert isinstance(categories, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/namelists', data_path + '.yml')) as f:
            content_str = f.read()
    except OSError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to open YAML database file. See server terminal output for details.\n'

    try:
        j2env = jinja2.Environment()
        j2env.globals['link_palm_repo_file'] = link_palm_repo_file
        template = j2env.from_string(content_str)
        content_str = template.render()
    except Exception as e:
        print_message_to_terminal(
            message_string='processing yaml database with jinja2 "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to processing YAML database file with jinja2. See server terminal output for details.\n'

    try:
        content_dict = yaml.load(content_str, Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_namelist_database(data_path, content_dict)
        if not valid:
            raise KeyError('namelist database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by category if required
    if 'all' not in categories:
        filtered_content_dict = dict()
        for parameter, parameter_dict in content_dict.items():
            if any(map(lambda c: c in parameter_dict['categories'], categories)):
                filtered_content_dict[parameter] = parameter_dict
        content_dict = filtered_content_dict

    if as_table:
        output_str = render_namelist_to_markdown_as_table(namelist, content_dict, link_table=link_table, link_path=link_path)
    else:
        output_str = render_namelist_to_markdown(namelist, content_dict, heading_level=heading_level)
    return output_str


def validate_yaml_logging_id_database(data_path, content_dict):
    """
    <logging_id>
      loglevel: ''
      message: ''
      description: ''
    """
    valid = True
    for logging_id, logging_id_dict in content_dict.items():

        if not 'loglevel' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "loglevel" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(logging_id_dict['loglevel'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    logging_id,
                    'The mandatory field "loglevel" must have a value of type "str"',
                )
                valid = False
            if not bool(re.match('^INFO$|^WARNING$|^ERROR$', logging_id_dict['loglevel'])):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    logging_id,
                    'The mandatory field "loglevel" must have a value of either INFO, WARNING or ERROR',
                )
                valid = False

        if not 'message' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "message" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(logging_id_dict['message'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    logging_id,
                    'The mandatory field "message" must have a value of type "str"',
                )
                valid = False

        if not 'description' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "description" is mandatory but missing',
            )
            valid = False
        else:
            if not isinstance(logging_id_dict['description'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    logging_id,
                    'The mandatory field "description" must have a value of type "str"',
                )
                valid = False
    return valid


def render_logging_ids_to_markdown_as_table(palm_module, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Logging ID', 'LogLevel', 'Message')
    output_str += format_str_head.format('-', '-', '-')
    for logging_id, logging_id_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--{2})'.format(link_path, palm_module, logging_id) if link_table else logging_id,
            logging_id_dict['loglevel'],
            str(logging_id_dict['message']).rstrip(),
        )

    output_str += '\n'
    return output_str


def render_logging_ids_to_markdown(palm_module, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = ''
    for logging_id, logging_id_dict in content_dict.items():
        output_str += '#' * heading_level + logging_id + ' {#' + '{0}'.format(logging_id) + '}\n\n'
        output_str += list_format_str.format(
            'LogLevel', '*{}*\n\n'.format(logging_id_dict['loglevel'])
        )
        output_str += list_format_str.format(
            'Message', '{}\n\n'.format(
                textwrap.indent(logging_id_dict['message'], ' ' * 4),
            )
        )
        output_str += list_format_str.format(
            'Description', '{}\n\n'.format(
                textwrap.indent(logging_id_dict['description'], ' ' * 4),
            )
        )
    output_str += '\n'
    return output_str


def include_palm_logging_ids(data_path, loglevel=['all'], as_table=False, link_table=True, link_path='', heading_level=4):
    call_string = 'include_palm_logging_ids(\'{}\', loglevel={}, as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, loglevel, as_table, link_table, link_path, heading_level,
    )

    palm_module = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(loglevel, str):
        loglevel = [loglevel]
    assert isinstance(loglevel, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/logging', data_path + '.yml')) as f:
            content_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'
    except yaml.scanner.ScannerError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". This could be caused by using a ":" inside an unquoted string. {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_logging_id_database(data_path, content_dict)
        if not valid:
            raise KeyError('logging_ids database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by loglevel if required
    if 'all' not in loglevel:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['loglevel'], loglevel)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    if as_table:
        output_str = render_logging_ids_to_markdown_as_table(palm_module, content_dict, link_table=link_table, link_path=link_path)
    else:
        output_str = render_logging_ids_to_markdown(palm_module, content_dict, heading_level=heading_level)
    return output_str


def validate_yaml_driver_database(database_path, database_dict):

    valid = True
    violation_counter = 0
    for root_level_field_name in ['global_attributes', 'dimensions', 'variables']:

        violation_counter += validate_yaml_field_generic(
            database_path=database_path,
            field_identifier='root level',
            content_dict=database_dict,
            field_name=root_level_field_name,
            field_type=dict,
            mandatory=True,
        )

    for global_attribute_name, global_attribute_dict in database_dict['global_attributes'].items():
        """
        global_attributes:
            <attribute_name>:
                mandatory: bool
                type: str
                si-unit: str
                description:
                    short: str
                    long: str
        """
        field_identifier = f'global_attributes.{global_attribute_name}'

        violation_counter += validate_yaml_field_mandatory(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=global_attribute_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_type(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=global_attribute_dict,
            mandatory=True,
            allowed_values=['NC_CHAR', 'NC_INT', 'NC_FLOAT'],
        )

        violation_counter += validate_yaml_field_si_unit(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=global_attribute_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_description(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=global_attribute_dict,
            mandatory=True,
        )

    for dimension_name, dimension_dict in database_dict['dimensions'].items():
        """
        dimensions:
            <dimension_name>:
                mandatory: bool
                size:
                    value: 7
                    value_of: str
                    depends_on: str
                description:
                    short: str
                    long: str
        """
        field_identifier = f'dimensions.{dimension_name}'

        violation_counter += validate_yaml_field_mandatory(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=dimension_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_size(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=dimension_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_description(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=dimension_dict,
            mandatory=True,
        )

    for variable_name, variable_dict in database_dict['variables'].items():
        """
        variables:
            <variable_name>:
                modules: list[str]
   #            mandatory: bool
   #            type: str
   #            si-unit: str
   #            shape: tuple
   #            attributes:
   #                <attribute_name>:
   #                    mandatory: bool
   #            default:
   #                value: <type>
   #                value_of: str
   #                depends_on: str
   #            description:
   #                short: str
   #                long: str
   #            allowed_values:
   #                - value: <type>
   #                  value_of: str
   #                  depends_on: str
   #                  description: str
        """
        field_identifier = f'variables.{variable_name}'

        violation_counter += validate_yaml_field_modules(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=True,
            allowed_values=['BIO', 'BCM', 'CCT', 'CHM', 'PAC', 'LPM', 'LSF', 'LSM', 'USM', 'RAD', 'PCM', 'IDM', 'MAS'],
        )

        violation_counter += validate_yaml_field_mandatory(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_type(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=True,
            allowed_values=['NC_CHAR', 'NC_BYTE', 'NC_INT', 'NC_FLOAT'],
        )
        variable_type_class = None
        if 'type' in variable_dict:
            if variable_dict['type'].startswith('NC_CHAR'):
                variable_type_class = str
            if variable_dict['type'].startswith('NC_BYTE'):
                variable_type_class = int
            if variable_dict['type'].startswith('NC_INT'):
                variable_type_class = int
            if variable_dict['type'].startswith('NC_FLOAT'):
                variable_type_class = float
        if variable_type_class is None:
            print(field_identifier)
            raise ValueError('Variable type class unknown')

        violation_counter += validate_yaml_field_shape(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=True,
        )

        violation_counter += validate_yaml_field_si_unit(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=False,
        )

        violation_counter += validate_yaml_field_generic(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            field_name='attributes',
            field_type=dict,
            mandatory=True,
        )
        if 'attributes' in variable_dict and isinstance(variable_dict['attributes'], dict):
            for attribute_name, attribute_dict in variable_dict['attributes'].items():

                violation_counter += validate_yaml_field_mandatory(
                    database_path=database_path,
                    field_identifier=f'{field_identifier}.attributes.{attribute_name}',
                    content_dict=attribute_dict,
                    mandatory=True,
                )

        if variable_type_class is not None:
            violation_counter += validate_yaml_field_default(
                database_path=database_path,
                field_identifier=field_identifier,
                content_dict=variable_dict,
                value_type=variable_type_class,
                mandatory=False,
            )

        violation_counter += validate_yaml_field_description(
            database_path=database_path,
            field_identifier=field_identifier,
            content_dict=variable_dict,
            mandatory=True,
        )

        if variable_type_class is not None:
            violation_counter += validate_yaml_field_allowed_values(
                database_path=database_path,
                field_identifier=field_identifier,
                content_dict=variable_dict,
                value_type=variable_type_class,
                mandatory=False,
            )
    return violation_counter == 0


def render_driver_global_attributes_to_markdown_as_table(driver, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Attribute', 'SI-Unit','Description')
    output_str += format_str_head.format('-', '-', '-', '-', '-')
    for global_attribute_name, global_attribute_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--global_attribute--{2})'.format(link_path, driver, global_attribute_name) if link_table else global_attribute_name,
            global_attribute_dict['si-unit'] if global_attribute_dict['si-unit'] is not None else '',
            global_attribute_dict['description']['short'],
        )

    output_str += '\n'
    return output_str


def render_driver_global_attributes_to_markdown(driver, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = '<br>\n'
    for global_attribute_name, global_attribute_dict in content_dict.items():
        output_str += '#' * heading_level + global_attribute_name + ' {#' + '{0}--global_attribute--{1}'.format(driver, global_attribute_name) + '}\n\n'
        output_str += list_format_str.format('SI-Unit', global_attribute_dict['si-unit'])
        output_str += list_format_str.format('Datatype', global_attribute_dict['type'])
        output_str += list_format_str.format('Mandatory', global_attribute_dict['mandatory'])
        output_str += '\n{}\n\n'.format(
            textwrap.indent(global_attribute_dict['description']['short'], ' ' * 4),
        )
        if 'long' in global_attribute_dict['description']:
            if global_attribute_dict['description']['long']:
                output_str += '{}\n\n'.format(
                    textwrap.indent(global_attribute_dict['description']['long'], ' ' * 4),
                )
        output_str += '\n<br>\n'
    output_str += '\n'
    return output_str


def include_palm_driver_global_attributes(data_path, as_table=False, link_table=True, link_path='', heading_level=3):
    call_string = 'include_palm_driver_global_attributes(\'{}\', as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, as_table, link_table, link_path, heading_level,
    )

    driver = os.path.splitext(os.path.basename(data_path))[0]

    try:
        with open(os.path.join(docs_dir, 'content/data/drivers', data_path + '.yml')) as f:
            content_str = f.read()
    except OSError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to open YAML database file. See server terminal output for details.\n'

    try:
        j2env = jinja2.Environment()
        j2env.globals['link_palm_repo_file'] = link_palm_repo_file
        template = j2env.from_string(content_str)
        content_str = template.render()
    except Exception as e:
        print_message_to_terminal(
            message_string='processing yaml database with jinja2 "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to processing YAML database file with jinja2. See server terminal output for details.\n'

    try:
        content_dict = yaml.load(content_str, Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_driver_database(data_path, content_dict)
        if not valid:
            raise KeyError('driver database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    if as_table:
        output_str = render_driver_global_attributes_to_markdown_as_table(driver, content_dict['global_attributes'], link_table=link_table, link_path=link_path)
    else:
        output_str = render_driver_global_attributes_to_markdown(driver, content_dict['global_attributes'], heading_level=heading_level)
    return output_str


def render_driver_dimensions_to_markdown_as_table(driver, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Dimension', 'Size','Description')
    output_str += format_str_head.format('-', '-', '-', '-', '-')
    for dimension_name, dimension_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--dimension--{2})'.format(link_path, driver, dimension_name) if link_table else dimension_name,
            dimension_dict['size']['value'],
            dimension_dict['description']['short'],
        )

    output_str += '\n'
    return output_str


def render_driver_dimensions_to_markdown(driver, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = '<br>\n'
    for dimension_name, dimension_dict in content_dict.items():
        output_str += '#' * heading_level + dimension_name + ' {#' + '{0}--dimension--{1}'.format(driver, dimension_name) + '}\n\n'
        output_str += list_format_str.format('Size', '*' + str(dimension_dict['size']['value']) + '*')
        output_str += list_format_str.format('Mandatory', dimension_dict['mandatory'])
        output_str += '\n{}\n\n'.format(
            textwrap.indent(dimension_dict['description']['short'], ' ' * 4),
        )
        if 'long' in dimension_dict['description']:
            if dimension_dict['description']['long']:
                output_str += '{}\n\n'.format(
                    textwrap.indent(dimension_dict['description']['long'], ' ' * 4),
                )
        output_str += '\n<br>\n'
    output_str += '\n'
    return output_str


def include_palm_driver_dimensions(data_path, as_table=False, link_table=True, link_path='', heading_level=3):
    call_string = 'include_palm_driver_dimensions(\'{}\', as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, as_table, link_table, link_path, heading_level,
    )

    driver = os.path.splitext(os.path.basename(data_path))[0]

    try:
        with open(os.path.join(docs_dir, 'content/data/drivers', data_path + '.yml')) as f:
            content_str = f.read()
    except OSError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to open YAML database file. See server terminal output for details.\n'

    try:
        j2env = jinja2.Environment()
        j2env.globals['link_palm_repo_file'] = link_palm_repo_file
        template = j2env.from_string(content_str)
        content_str = template.render()
    except Exception as e:
        print_message_to_terminal(
            message_string='processing yaml database with jinja2 "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to processing YAML database file with jinja2. See server terminal output for details.\n'

    try:
        content_dict = yaml.load(content_str, Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_driver_database(data_path, content_dict)
        if not valid:
            raise KeyError('driver database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    if as_table:
        output_str = render_driver_dimensions_to_markdown_as_table(driver, content_dict['dimensions'], link_table=link_table, link_path=link_path)
    else:
        output_str = render_driver_dimensions_to_markdown(driver, content_dict['dimensions'], heading_level=heading_level)
    return output_str


def render_driver_variables_to_markdown_as_table(driver, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Variable', 'Shape','Description')
    output_str += format_str_head.format('-', '-', '-', '-', '-')
    for variable_name, variable_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--variable--{2})'.format(link_path, driver, variable_name) if link_table else variable_name,
            variable_dict['shape'] if variable_dict['shape'] is not None else '',
            variable_dict['description']['short'],
        )

    output_str += '\n'
    return output_str


def render_driver_variables_to_markdown(driver, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = '<br>\n'
    for variable_name, variable_dict in content_dict.items():
        output_str += '#' * heading_level + variable_name + ' {#' + '{0}--variable--{1}'.format(driver, variable_name) + '}\n\n'
        output_str += list_format_str.format('Shape', variable_dict['shape'])
        output_str += list_format_str.format('Datatype', variable_dict['type'])
        output_str += list_format_str.format('Default', '*' + str(variable_dict['default']['value']) + '*')
        if 'si-unit' in variable_dict:  # and re.search('^[I,R].*', parameter_dict['type']):
            output_str += list_format_str.format('SI-Unit', variable_dict['si-unit'])
        output_str += list_format_str.format('Mandatory', variable_dict['mandatory'])
        if len(variable_dict['attributes']) > 0:
            output_str += list_format_str.format('Attributes', '')
            attribute_string_list = []
            for attribute_name, attribute_dict in variable_dict['attributes'].items():
                attribute_string_list.append(
                    '{}{}'.format(
                        attribute_name,
                        ' (mandatory)' if attribute_dict['mandatory'] else '',
                    )
                )
            output_str += '{}\n\n'.format(', '.join(attribute_string_list))
        output_str += '\n{}\n\n'.format(
            textwrap.indent(variable_dict['description']['short'], ' ' * 4),
        )
        if 'long' in variable_dict['description']:
            if variable_dict['description']['long']:
                output_str += '{}\n\n'.format(
                    textwrap.indent(variable_dict['description']['long'], ' ' * 4),
                )
        if 'allowed_values' in variable_dict:
            output_str += '\n{}\n\n'.format(
                textwrap.indent('Currently {} choices are available:'.format(len(variable_dict['allowed_values'])), ' ' * 4),
            )
            for allowed_value in variable_dict['allowed_values']:
                output_str += '    - *{}*\n\n{}\n\n'.format(
                    allowed_value['value'],
                    textwrap.indent(allowed_value['description'], ' ' * 8),
                )

        output_str += '\n<br>\n'
    output_str += '\n'
    return output_str


def include_palm_driver_variables(data_path, modules=['all'], as_table=False, link_table=True, link_path='', heading_level=3):
    call_string = 'include_palm_driver_variables(\'{}\', as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, as_table, link_table, link_path, heading_level,
    )

    driver = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(modules, str):
        modules = [modules]
    assert isinstance(modules, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/drivers', data_path + '.yml')) as f:
            content_str = f.read()
    except OSError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to open YAML database file. See server terminal output for details.\n'

    try:
        j2env = jinja2.Environment()
        j2env.globals['link_palm_repo_file'] = link_palm_repo_file
        template = j2env.from_string(content_str)
        content_str = template.render()
    except Exception as e:
        print_message_to_terminal(
            message_string='processing yaml database with jinja2 "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to processing YAML database file with jinja2. See server terminal output for details.\n'

    try:
        content_dict = yaml.load(content_str, Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_driver_database(data_path, content_dict)
        if not valid:
            raise KeyError('driver database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by module if required
    if 'all' not in modules:
        filtered_content_dict = dict()
        for variable_name, variable_dict in content_dict['variables'].items():
            if any(map(lambda c: c in variable_dict['modules'], modules)):
                filtered_content_dict[variable_name] = variable_dict
        content_dict['variables'] = filtered_content_dict

    if as_table:
        output_str = render_driver_variables_to_markdown_as_table(driver, content_dict['variables'], link_table=link_table, link_path=link_path)
    else:
        output_str = render_driver_variables_to_markdown(driver, content_dict['variables'], heading_level=heading_level)
    return output_str


def validate_yaml_output_quantities_database(data_path, content_dict):
    """
    <output_quantity>
      scope:
        - palm_core
        - land_surface_model
      type:
        - vertical profile
        - 2d-array
        - 3d-array
        - masked array
      si-unit: %
      quantity: coverage of the land surface with bare soil
      remarks: ''
    """
    valid = True
    for output_quantity, output_quantity_dict in content_dict.items():

        if not 'scope' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "scope" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['scope'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "scope" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in output_quantity_dict['scope']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "scope" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False

        if not 'type' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "type" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['type'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in output_quantity_dict['type']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False
        if not all(bool(re.match('^vertical profile$|^2d-array$|^3d-array$|^masked array$', v)) for v in output_quantity_dict['type']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of either "vertical profile", "2d-array", "3d-array" or "masked array"',
            )
            valid = False

        if isinstance(output_quantity_dict['si-unit'], int):
            output_quantity_dict['si-unit'] = str(output_quantity_dict['si-unit'])
        if not 'si-unit' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "si-unit" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['si-unit'], str) and output_quantity_dict['si-unit'] != 1:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "si-unit" must have a value of type "str"',
            )
            valid = False

        if not 'description' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "description" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['description'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "description" must have a value of type "str"',
            )
            valid = False

        if not 'remarks' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "remarks" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['remarks'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "remarks" must have a value of type "str"',
            )
            valid = False
    return valid


def render_output_quantities_to_markdown_as_table(content_dict, show_remarks=True):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Name { .oq-name }', 'SI-Unit { .oq-unit }', 'Description { .oq-description }')
    output_str += format_str_head.format('-', '-', '-')
    for output_quantity, output_quantity_dict in content_dict.items():
        output_quantity_excaped = output_quantity.replace('>', r'\>').replace('_', r'\_').replace('*', r'\*')
        description_excaped = output_quantity_dict['description'].replace('>', r'\>').replace('_', r'\_').replace('*', r'\*')
        remarks_excaped = output_quantity_dict['remarks'].replace('>', r'\>').replace('_', r'\_').replace('*', r'\*')
        output_str += format_str_body.format(
            '{}'.format(output_quantity_excaped),
            output_quantity_dict['si-unit'],
            '{} {}'.format(
                str(description_excaped).rstrip(),
                '<br><br>*'+str(remarks_excaped).rstrip()+'*' if output_quantity_dict['remarks'] and show_remarks else '',
            ),
        )

    output_str += '\n'
    return output_str


def include_palm_output_quantities(data_path, scopes=['all'], types=['all'], show_remarks=True):
    call_string = 'include_palm_logging_ids(\'{}\', scopes={}, types={})'.format(
        data_path, scopes, types
    )

    palm_module = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(scopes, str):
        scopes = [scopes]
    assert isinstance(scopes, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/output_quantities', data_path + '.yml')) as f:
            content_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'
    except yaml.scanner.ScannerError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". This could be caused by using a ":" inside an unquoted string. {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_output_quantities_database(data_path, content_dict)
        if not valid:
            raise KeyError('output_quantities database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by scope if required
    if 'all' not in scopes:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['scope'], scopes)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    # filter by type if required
    if 'all' not in types:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['type'], types)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    output_str = render_output_quantities_to_markdown_as_table(content_dict, show_remarks=show_remarks)
    return output_str


def link_palm_repo_file(link_text, file_path, new_tab=True):
    output_str = f'[{link_text}]({file_link_url}/{file_path})' + '{ target=_blank }' if new_tab else ''
    return output_str
