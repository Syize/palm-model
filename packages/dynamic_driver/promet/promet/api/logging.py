import logging
import hashlib
from pprint import pprint

import numpy
from termcolor import colored
import os
import zlib
import base64
import subprocess
import platform
import json
import http
import traceback
from datetime import datetime
import time

execution_start_time = time.time()

program_root_dir = os.path.dirname(os.path.abspath(os.path.join(__file__, '..', '..')))


class PrometException(Exception):
    def __init__(self, message, id):
        super().__init__(message)
        self.id = id


def get_git_commit_id():
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            capture_output=True,
            cwd=os.path.dirname(os.path.realpath(__file__)),
            text=True,
            check=True,
        ).stdout.strip()
        git_commit_id = f'git_{result}'
    except:
        git_commit_id = 'git_unknown'
    return git_commit_id


try:
    from promet.version_info import git_commit_id
except:
    git_commit_id = get_git_commit_id()


try:
    from promet.version_info import package_version
except:
    try:
        from importlib.metadata import version, metadata

        package_version = version("promet")
    except:
        package_version = "pip_unknown"

def print_and_log_section(s):
    print(f'\n{s}\n')
    logging_database.add_checkpoint(s)
    logging.info(s)

def print_and_log_subsection(s):
    print(f'\n   {s}')
    logging_database.add_checkpoint(s)
    logging.info(f'   {s}')

def print_and_log_step(s):
    print(f'   {s}')
    logging.info(f'   {s}')

def print_info(s, *args, **kwargs):
    if args or kwargs:
        print(colored('INFO: ', 'green') + s.format(*args, **kwargs))
    else:
        print(colored('INFO: ', 'green') + s)


def print_warn(s, *args, **kwargs):
    if args or kwargs:
        print(colored('WARNING: ', 'yellow') + s.format(*args, **kwargs))
        logging_database.add_warning(message=s.format(*args, **kwargs))
    else:
        print(colored('WARNING: ', 'yellow') + s)
        logging_database.add_warning(message=s)


def print_error(s, *args, **kwargs):
    if args or kwargs:
        print(colored('ERROR: ', 'red') + s.format(*args, **kwargs))
        logging_database.add_error(message=s.format(*args, **kwargs))
    else:
        print(colored('ERROR: ', 'red') + s)
        logging_database.add_error(message=s)


def check_version():
    query = '''query promet($current: String!) { check_version(input: {current: $current}) { show message } }'''
    variables = {'current': 'v25.10-rc.1'}
    data = request_data(
        query=query,
        variables=variables,
    )
    if 'check_version' in data:
        show_message = data['check_version']['show']
        message = data['check_version']['message']
        if show_message:
            print(f'{message}')


class CustomEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, datetime):
            #return obj.isoformat()
            return obj.strftime('%Y-%m-%d %H:%M:%S +00')
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return super().default(obj)


def check_for_known_errors() -> bool:
    try:
        json_string = json.dumps(
            logging_database.get_data_package(),
            cls=CustomEncoder,
        )
    except TypeError as e:
        json_string = json.dumps(
            dict(
                elapsed_runtime=time.time()-execution_start_time,
                exceptions=[
                    dict(
                        type='TypeError',
                        id='',
                        message=str(e),
                        traceback=[],
                    )
                ]
            ),
        )
    json_hash = hashlib.sha1(json_string.encode('utf-8')).hexdigest()
    compressed_data = zlib.compress(json_string.encode('utf-8'))
    encoded_data = base64.urlsafe_b64encode(compressed_data).decode('utf-8')
    query = '''query promet($hash: String!, $data: String!) { check_for_known_errors(input: {hash: $hash data: $data}) { catch show message } }'''
    variables = {'hash': json_hash, 'data': encoded_data}
    data = request_data(
        query=query,
        variables=variables,
    )
    catch = False
    if 'check_for_known_errors' in data:
        catch = bool(data['check_for_known_errors']['catch'])
        show_message = bool(data['check_for_known_errors']['show'])
        message = data['check_for_known_errors']['message']
        if show_message:
            print()
            print(f'{message}')
    return catch


class LoggingDatabase:

    def __init__(self):
        self.data = dict(
            namespace=dict(),
            config=dict(),
            checkpoints=[],
            metadata=dict(),
            warnings=[],
            errors=[],
            exceptions=[],
        )

    def get_data_package(self):
        data_package = dict()
        for key, value in self.data.items():
            data_package[key] = value
        data_package['elapsed_runtime'] = time.time()-execution_start_time
        return data_package

    def add_checkpoint(self, message):
        self.data['checkpoints'].append(
            dict(
                elapsed_runtime=time.time()-execution_start_time,
                message=message,
            )
        )

    def add_namespace(self, namespace):
        for key, value in namespace.__dict__.items():
            if key == 'func':
                continue
            self.data['namespace'][key] = value

    def add_config(self, config):
        for key, value in config.content.items():
            if key.startswith('input_data_'):
                if isinstance(value, str):
                    value = os.path.basename(value)
            self.data['config'][key] = value

    def add_metadata(self, key, metadata):
        self.data['metadata'][key] = metadata

    def add_warning(self, message):
        self.data['warnings'].append(message)

    def add_error(self, message):
        self.data['errors'].append(message)

    def add_exception(self, exc: Exception):
        exc_type = type(exc)
        exc_traceback = traceback.format_exception(exc_type, exc, exc.__traceback__)

        for n, line in enumerate(exc_traceback):
            if line.startswith('  File "'):
                abs_path = line.split('"')[1]
                if abs_path.startswith('/'):
                    if abs_path.startswith(program_root_dir):
                        rel_path = os.path.relpath(abs_path, start=program_root_dir)
                        line = line.replace(abs_path, rel_path)
                        exc_traceback[n] = line
                    else:
                        rel_path = os.path.basename(abs_path)
                        line = line.replace(abs_path, rel_path)
                        exc_traceback[n] = line
        if hasattr(exc, 'id'):
            exc_id = exc.id
        else:
            exc_id = 'missing'
        self.data['exceptions'].append(
            dict(
                type=exc_type.__name__,
                id=exc_id,
                message=str(exc),
                traceback=exc_traceback,
            )
        )


def request_data(query, variables):
    try:
        h = dict()
        h['Content-Type'] = 'application/json'
        h['User-Agent'] = f'{platform.platform()} {platform.version()} with Python {platform.python_version()}'
        h['Authorization'] = f'basic cHJvbWV0X3YyNS4xMC1yYy4xOkpVUWVWU0FrTERkTGRocElTWDlFRDBOQg=='
        c = http.client.HTTPConnection('api.palm-model.com', 80, timeout=2)
        c.request('POST', '/v1', body=json.dumps({'query': query, 'variables': variables}), headers=h)
        r = c.getresponse()
        dc = r.read().decode('utf-8')
        dd = json.loads(dc)
        if 'data' in dd and isinstance(dd['data'], dict):
            data = dd['data']
        else:
            data = dict()
    except:
        data = dict()
    return data


logging_database = LoggingDatabase()
