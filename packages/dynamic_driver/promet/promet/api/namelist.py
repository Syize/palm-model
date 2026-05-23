import os
import f90nml
import numpy

from promet.api.logging import PrometException
from promet.api.logging import print_error


class Namelist:

    required_parameters = {
        'initialization_parameters': [
            'nx',
            'ny',
            'nz',
            'dx',
            'dy',
            'dz',
            'dz_stretch_level',
            'dz_stretch_factor',
            'dz_max',
            'initializing_actions',
        ],
        'runtime_parameters': [
            'end_time',
        ]
    }

    parameter_defaults = {
        'dz_stretch_level': 1.0e32,
        'dz_stretch_factor': 1.08,
        'dz_max': 1.0e32,
        'rotation_angle': 0.0,
    }

    def __init__(self, filepath):
        self.filepath = os.path.normpath(os.path.expandvars(os.path.expanduser(filepath)))
        self._content = dict()
        self._parameters = dict()

    def __getitem__(self, key):
        if key in self.parameter_defaults:
            return self._content.get(key, self.parameter_defaults[key])
        else:
            return self._content[key]

    @property
    def parameters(self):
        return self._parameters

    def load(self) -> bool:
        try:
            self._content = f90nml.read(self.filepath)
        except FileNotFoundError as e:
            raise PrometException(
                f'Unable to find namelist file: {self.filepath}',
                id='TEqapQ',
            )
        except UnicodeDecodeError as e:
            raise PrometException(
                f'Unable to open namelist file: {self.filepath}',
                id='VkftLw',
            )
        if self._content is None:
            raise PrometException(
                f'Namelist file is empty.',
                id='kHVBmw',
            )
        for namelist_name in self.required_parameters.keys():
            if self[namelist_name] is None:
                raise PrometException(
                    f'Unable to find required namelist "{namelist_name}"',
                    id='PRBhXg',
                )
        missing_parameters = []
        for namelist_name, namelist_parameters in self.required_parameters.items():
            for namelist_parameter in namelist_parameters:
                try:
                    self._parameters[namelist_parameter] = self[namelist_name][namelist_parameter]
                except KeyError as e:
                    if namelist_parameter in self.parameter_defaults:
                        self._parameters[namelist_parameter] = self.parameter_defaults[namelist_parameter]
                    else:
                        missing_parameters.append(namelist_parameter)
                        print_error(f'Unable to find parameter "{namelist_parameter}" in namelist "{namelist_name}"')
        if len(missing_parameters) > 0:
            raise PrometException(
                f'Missing {len(missing_parameters)} required namelist parameters',
                id='PRBhXg',
            )
        return False

    @staticmethod
    def get_z(dz, nz, dz_stretch_level=1.0e32, dz_stretch_factor=1.08, dz_max=1.0e32):
        z = [dz/2]
        while z[-1] < dz_stretch_level and len(z) < nz:
            z.append(z[-1] + dz)
        while len(z) < nz:
            dz = min(dz * dz_stretch_factor, dz_max)
            z.append(z[-1] + dz)
        assert len(z) == nz
        zw = [0.0]
        for i, dummy in enumerate(z[0:-1]):
            zw.append((z[i] + z[i+1]) * 0.5)
        if z[-1] > dz_stretch_level:
            size_z = z[-1] + ((z[-1] - z[-2]) * 0.5) * dz_stretch_factor
        else:
            size_z = z[-1] + ((z[-1] - z[-2]) * 0.5)
        zw.append(size_z)
        z = [0.0] + z
        k = [k for k in range(nz+1)]
        assert len(z) == len(zw), str(len(z))+' != '+str(len(zw))
        assert len(z) == len(k), str(len(z))+' != '+str(len(k))
        return z, zw, k

    def get_dimension(
            self,
            dimension: str
    ):
        if dimension == 'z':
            z, zw, k = self.get_z(
                dz=self.parameters['dz'],
                nz=self.parameters['nz'],
                dz_stretch_level=self.parameters['dz_stretch_level'],
                dz_stretch_factor=self.parameters['dz_stretch_factor'],
                dz_max=self.parameters['dz_max'],
            )
            z = numpy.array(z[1:])
            result = dict(
                size=z.size,
                long_name='height above origin',
                units='m',
                values=z,
            )
            return result
        else:
            raise ValueError(f'Dimension {dimension} is not supported')

    def get_metadata(self) -> dict:
        return dict(
            nx=self.parameters['nx'],
            ny=self.parameters['ny'],
            nz=self.parameters['nz'],
            dx=self.parameters['dx'],
            dy=self.parameters['dy'],
            dz=self.parameters['dz'],
            dz_stretch_level=self.parameters['dz_stretch_level'],
            dz_stretch_factor=self.parameters['dz_stretch_factor'],
            dz_max=self.parameters['dz_max'],
            end_time=self.parameters['end_time'],
        )
