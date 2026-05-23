import subprocess
from setuptools.command.install import install
import os

class PrometInstallCommand(install):

    def run(self):
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
        try:
            package_version = self.distribution.get_version()
        except:
            package_version = 'pip_unknown'
        try:
            build_py = self.get_finalized_command('build_py')
            build_dir = build_py.build_lib
            package_name = self.distribution.get_name()
            package_dir = os.path.join(build_dir, package_name)
            git_info_file = os.path.join(package_dir, 'version_info.py')
            with open(git_info_file, 'w') as f:
                f.write(f"git_commit_id = '{git_commit_id}'\n")
                f.write(f"package_version = '{package_version}'\n")
        except:
            pass
        install.run(self)
