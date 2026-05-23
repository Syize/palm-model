#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

from promet.cli import CLI

if __name__ == '__main__':
    cli = CLI()
    cli.execute_all()
