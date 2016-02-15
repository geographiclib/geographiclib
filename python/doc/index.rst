.. geographiclib documentation master file, created by
   sphinx-quickstart on Sat Oct 24 17:14:00 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

geographiclib
=============

Author: Charles F. F. Karney (charles@karney.com)

Version: |version|

The documentation for other versions is available at
``http://geographiclib.sourceforge.net/m.nn/python/`` for versions numbers
``m.nn`` ≥ 1.46.

Licensed under the MIT/X11 License; see
`LICENSE.txt <http://geographiclib.sourceforge.net/html/LICENSE.txt>`_.

Introduction
============

This is a python implementation of the geodesic routines in
`GeographicLib <http://geographiclib.sourceforge.net>`_.

Although it is maintained in conjunction with the larger C++ library,
this python package can be used independently.

Installation
------------

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

Contents
--------

.. toctree::
   :maxdepth: 2

   geodesics
   interface
   code
   examples

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

Change log
----------

* Version 1.46 (released 2016-02-15)

  * Add Geodesic.DirectLine, Geodesic.ArcDirectLine,
    Geodesic.InverseLine, GeodesicLine.SetDistance, GeodesicLine.SetArc,
    GeodesicLine.s13, GeodesicLine.a13.
  * More accurate inverse solution when longitude difference is close to
    180°.
  * Remove unnecessary functions, CheckPosition, CheckAzimuth,
    CheckDistance.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
