# setup.py, config file for distutils
#
# To install this package, execute
#
#   python setup.py install
#
# in this directory.  To upload the latest version to the python
# repository run
#
#   python setup.py sdist --formats gztar,zip upload
#
# The initial version of this file was provided by
# Andrew MacIntyre <Andrew.MacIntyre@acma.gov.au>.

from distutils.core import setup
from distutils.cmd import Command
from sphinx.setup_command import BuildDoc


class TestCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys, subprocess
        raise SystemExit(subprocess.call([sys.executable,
                                          '-m',
                                          'unittest',
                                          '-v',
                                          'test.test_geodesic'
                                          ]))


name = "geographiclib"
version = "1.46"
release = "1.46.0"

setup(name = name,
      version = version,
      description =
        "A translation of the GeographicLib::Geodesic class to Python",
      author = "Charles Karney",
      author_email = "charles@karney.com",
      url = "http://geographiclib.sourceforge.net/",
      packages = ["geographiclib", "geographiclib/test"],
      data_files = [],
      license = "MIT",
      keywords = "gis geographical earth distance geodesic",
      classifiers = [
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Topic :: Scientific/Engineering :: GIS",
          "Topic :: Software Development :: Libraries :: Python Modules",
      ],

      cmdclass={
          'test': TestCommand,
          'build_sphinx': BuildDoc,
          'build_sphinx_latex': BuildDoc,
      },

      command_options={
          # these are optional and override conf.py settings
          'build_sphinx': {
              'project': ('setup.py', name),
              'version': ('setup.py', version),
              'release': ('setup.py', release),
          },
          'build_sphinx_latex': {
              'project': ('setup.py', name),
              'version': ('setup.py', version),
              'release': ('setup.py', release),
          }
      },

      )
