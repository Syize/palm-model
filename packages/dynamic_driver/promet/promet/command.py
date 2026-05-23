from argparse import Namespace
import logging

from promet.api.config import Config
from promet.api.dynamic_driver import DynamicDriver
from promet.api.logging import logging_database
from promet.api.namelist import Namelist
from promet.api.palm_setup import PALMSetup
from promet.api.static_driver import StaticDriver
from promet.plugins.icon import ICONPlugin
from promet.plugins.wrf import WRFPlugin


class Command:

    """ The Command class.

    This is the base class for any command controller.

    Attributes:
       args (Namespace): The argument namespace for the controller.

    """

    def __init__(self, args):
        self.args = args
        self.failed = False
        self.config = None
        self.plugin = None

    @classmethod
    def create(cls, args):
        """ Create the command object. This is called from the parser. """
        logging.info('{0} called with: {1}'.format(__name__, args))
        command = cls(args)
        return command

    def start(self):
        """ This is called first at program start."""
        logging_database.add_namespace(self.args)

        self.config = Config(
            filepath=self.args.config,
            verbose=self.args.verbose,
        )
        self.failed = self.config.load()
        logging_database.add_config(config=self.config)
        if self.failed:
            return

        self.namelist = Namelist(
            filepath=self.args.namelist,
        )
        self.failed = self.namelist.load()
        if self.failed:
            return

        self.static_driver = StaticDriver(
            filepath=self.args.static_driver,
        )
        self.failed = self.static_driver.load()
        if self.failed:
            return

        self.palm_setup = PALMSetup(
            namelist=self.namelist,
            static_driver=self.static_driver,
        )

        self.dynamic_driver = DynamicDriver(
            filepath=self.args.output_file,
        )

        if self.config['plugin'] == 'icon':
            self.plugin = ICONPlugin(
                config=self.config,
                dynamic_driver=self.dynamic_driver,
                palm_setup=self.palm_setup,
            )
        elif self.config['plugin'] == 'wrf':
            self.plugin = WRFPlugin(
                config=self.config,
                dynamic_driver=self.dynamic_driver,
                palm_setup=self.palm_setup,
            )
        else:
            raise Exception(f'Unknown plugin "{self.config["plugin"]}"')

    def interact(self):
        """ This is called second to handle user interaction."""
        if self.failed:
            return

        self.failed = self.palm_setup.load_global_attributes(verbose=self.args.verbose)
        if self.failed:
            return
        self.failed = self.palm_setup.load_dimensions(
            verbose=self.args.verbose,
            zsoil=self.config['zsoil'],
        )
        if self.failed:
            return
        self.failed = self.dynamic_driver.initialize(
            append=self.args.append,
            overwrite=self.args.overwrite,
        )
        if self.failed:
            return

        self.dynamic_driver.set_global_attributes(
            origin_lat=self.palm_setup.global_attributes['origin_lat'],
            origin_lon=self.palm_setup.global_attributes['origin_lon'],
            origin_time=self.palm_setup.global_attributes['origin_time'],
            origin_x=self.palm_setup.global_attributes['origin_x'],
            origin_y=self.palm_setup.global_attributes['origin_y'],
            origin_z=self.palm_setup.global_attributes['origin_z'],
            rotation_angle=self.palm_setup.global_attributes['rotation_angle'],
            verbose=self.args.verbose,
        )

        self.dynamic_driver.set_dimensions(
            x=self.palm_setup.dimensions['x'],
            xu=self.palm_setup.dimensions['xu'],
            y=self.palm_setup.dimensions['y'],
            yv=self.palm_setup.dimensions['yv'],
            z=self.palm_setup.dimensions['z'],
            zw=self.palm_setup.dimensions['zw'],
            zsoil=self.palm_setup.dimensions['zsoil'],
            verbose=self.args.verbose,
        )

        self.failed = self.plugin.initialize(
            verbose=self.args.verbose,
        )
        logging_database.add_metadata(key='palm_setup', metadata=self.palm_setup.get_metadata())
        logging_database.add_metadata(key='plugin', metadata=self.plugin.get_metadata())

    def execute(self):
        """ This is called third to execute the main task."""
        if self.failed:
            return
        self.failed = self.plugin.process_palm_variables(
            verbose=self.args.verbose,
        )
