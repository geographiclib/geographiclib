# setup.py, config file for distutils
#
# To build a windows installer for this package, execute
#
#   python setup.py bdist_wininst
#
# in this directory.  The resulting installer will be in the dist
# subdirectory.
#
# $Id$

from distutils.core import setup

setup(name="geographiclib",
      version="0.13",
      description=
        "A translation of the GeographicLib::Geodesic class to Python",
      author="Charles Karney",
      author_email="charles@karney.com",
      url="http://geographiclib.sourceforge.net",
      packages=["geographiclib"],
      data_files=[],
      license="MIT",
      )
