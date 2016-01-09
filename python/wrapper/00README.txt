Here is a small example of exposing one capability of the C++ version
GeographicLib to python.  PyGeographicLib.cpp uses boost-python to wrap
parts of 2 classes, Constants and Geodesic, so that they can be invoked
from python.

To build and install this interface, do

  mkdir BUILD
  cd BUILD
  cmake -D CMAKE_INSTALL_PREFIX=/home/ckarney/.local ..
  make && make install

This assumes that you have installed GeographicLib somewhere that cmake
can find it.  If you want just to use the version of GeographicLib that
you have built in the top-level BUILD directory, include, e.g.,

  -D GeographicLib_DIR=../../BUILD

in the invocation of cmake.

"make install" installs PyGeographicLib.so in

  ~/.local/lib/python2.7/site-packages

which is in the default search path for python 2.7.  To use this in
python, do

  $ python
  >>> from PyGeographicLib import Constants, Geodesic
  >>> geod = Geodesic(Constants.WGS84_a(), Constants.WGS84_f())
  >>> geod.Inverse(40.63972, -73.77889, 1.35917, 103.98944)
  (3.3083172599825987, 177.48589836536408, 15347627.52149483,
   138.05222403153218, 4302458.117647906,
   -0.7373981702672949, -0.7435561002478325,
   123377685741829.42)
  >>> help(Geodesic.Inverse)

Notes:

This precription applies to Linux machines.  Similar steps can be used
on Windows and MacOSX machines.

You will need the packages, boost-python, boost-devel, python, and
python-devel installed.

CMakeLists.txt specifies the version of python to look for (version
2.7).  This must match that used in boost-python.  To check do

  ldd /usr/lib64/libboost_python.so

Acknowledgment:

Thanks to Jonathan Takahashi <jtakahashi@gmail.com> for the sample code
in PyGeographicLib.cpp and the commands needed to compile and link this.
