import logging
import argcomplete
import sys
import shutil

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from . command import Command
from promet.api.logging import check_for_known_errors
from promet.api.logging import check_version
from promet.api.logging import logging_database


version = '1.0.0'


class PrometArgumentParser(ArgumentParser):
    """ The PrometArgumentParser class.

    An instance of this class represents the argument parser that all entry
    points of this package have in common.

    """
    def __init__(self):

        super().__init__(
            description=f'This is promet {version}\n' +
                        'Developer Support: support@pecanode.com',
            formatter_class=RawTextHelpFormatter,
            epilog=f'Copyright pecanode GmbH (AGPLv3)',
            add_help=True,
            allow_abbrev=False,
        )

        self.set_defaults(func=Command.create)

        self.add_argument(
            '--version',
            action='version',
            version=f'promet {version}',
        )

        self.add_argument(
            '--config',
            action='store',
            dest='config',
            help='path to config file',
            required=True,
            metavar='PATH',
        )

        self.add_argument(
            '--namelist',
            action='store',
            dest='namelist',
            help='path to namelist file of PALM setup',
            required=True,
            metavar='PATH',
        )

        self.add_argument(
            '--static-driver',
            action='store',
            dest='static_driver',
            help='path to static driver file of PALM setup',
            required=True,
            metavar='PATH',
        )

        self.add_argument(
            '--output-file',
            action='store',
            dest='output_file',
            default='dynamic_driver.nc',
            help='desired output filepath of ',
            required=False,
            metavar='PATH',
        )

        self.add_argument(
            '--append',
            action='store_true',
            dest='append',
            help='append to output file if existing',
            required=False,
        )

        self.add_argument(
            '--overwrite',
            action='store_true',
            dest='overwrite',
            help='overwrite output file if existing',
            required=False,
        )

        self.add_argument(
            '--verbose',
            action='store_true',
            dest='verbose',
            help='enable increased verbosity',
            required=False,
        )


class CLI:
    """ Command Line Interface (CLI) class.

    An instance of this class represents the CLI for this package.

    """
    def __init__(self, args=None, namespace=None):
        self.parser = PrometArgumentParser()
        argcomplete.autocomplete(self.parser)  # Autocomplete arguments here
        logging.basicConfig(
            filename='promet.log',
            level=logging.INFO,
            format='%(asctime)-23s %(levelname)-9s: %(message)s',
            #format='%(asctime)-23s %(levelname)-9s %(pathname)s(%(lineno)d): %(message)s',
        )
        self.args = self.parser.parse_args( # Parse arguments here
            args=args,
            namespace=namespace,
        )
        check_version()
        terminal_columns, terminal_lines = shutil.get_terminal_size()
        hline = '#' * min(terminal_columns, 300)
        print(hline)
        print(f'This is promet {version}')
        print(hline)
        logging.info('{0} called with: {1}'.format(__name__, self.args))
        self.command = self.args.func(args=self.args)
        logging.info('%s initialized', __name__)

    def start(self):
        """ Start the program."""
        self.command.start()
        logging.info('%s started', __name__)

    def interact(self):
        """ Interact with the program."""
        logging.info('%s interaction started', __name__)
        self.command.interact()
        logging.info('%s interaction finished', __name__)

    def execute(self):
        """ Execute the program."""
        logging.info('%s execution started', __name__)
        self.command.execute()
        logging.info('%s execution finished', __name__)

    def quit(self):
        """ Quit the program."""
        self.exit_state = int(self.command.failed)
        logging.info('%s quit with exit code %s', __name__, self.exit_state)
        sys.exit(self.exit_state)

    def execute_all(self):
        try:
            self.start()
            self.interact()
            self.execute()
        except Exception as e:
            logging_database.add_exception(e)
            catch = check_for_known_errors()
            if not catch:
                raise e
        except KeyboardInterrupt as e:
            e.id = '6ULuVQ'
            logging_database.add_exception(e)
            catch = check_for_known_errors()
            if not catch:
                raise e
        else:
            check_for_known_errors()
        self.quit()