Installing the package
======================

The native package
------------------

The full `Geographic <http://geographiclib.sourceforge.net>`_ package
can be downloaded from
`sourceforge <http://sourceforge.net/projects/geographiclib/files/distrib>`_.
However the python implementation is available as a stand-alone package.
To install this, run

.. code-block:: sh

  pip install geographiclib

Alternatively downloaded the package directly from
`Python Package Index <http://pypi.python.org/pypi/geographiclib>`_
and install it with

.. code-block:: sh

  tar xpfz geographiclib-1.46.tar.gz
  cd geographiclib-1.46
  python setup.py install

It's a good idea to run the unit tests to verify that the installation
worked OK by running

.. code-block:: sh

  python -m unittest geographiclib.test.test_geodesic

Python wrapper for GeographicLib
--------------------------------

The python package only implements a subset of the C++ library
GeographicLib.  In order to access the full capabilities of
GeographicLib, you can use
`boost-python <http://www.boost.org/doc/libs/release/libs/python/>`_
to wrap the functionality you need into python classes.  You'll need to
install the GeographicLib library and boost and to have the ability to
compile C++ code.  A sample is given in the source distribution for the
C++ library in the file ``wrapper/python/PyGeographicLib.cpp``.  The
file ``00README.txt`` in the same directory provides instructions for
how to compile this.

GeographicLib in various languages
----------------------------------

* C++ (complete library):
  `documentation <../index.html>`__,
  `download <https://sourceforge.net/projects/geographiclib/files/distrib>`_
* C (geodesic routines):
  `documentation <../C/index.html>`__,
  also included with recent versions of
  `proj.4 <https://github.com/OSGeo/proj.4/wiki>`_
* Fortran (geodesic routines):
  `documentation <../Fortran/index.html>`__
* Java (geodesic routines):
  `Maven Central package <http://repo1.maven.org/maven2/net/sf/geographiclib/GeographicLib-Java/>`_,
  `documentation <../java/index.html>`__
* JavaScript (geodesic routines):
  `npm package <https://www.npmjs.com/package/geographiclib>`_,
  `documentation <../js/index.html>`__
* Python (geodesic routines):
  `PyPI package <http://pypi.python.org/pypi/geographiclib>`_,
  `documentation <../python/index.html>`__
* Matlab/Octave (geodesic and some other routines):
  `Matlab Central package <http://www.mathworks.com/matlabcentral/fileexchange/50605>`_,
  `documentation
  <http://www.mathworks.com/matlabcentral/fileexchange/50605/content/Contents.m>`__
* C# (.NET wrapper for complete C++ library):
  `documentation <../NET/index.html>`__
