Invoking the compiler from Matlab or Octave (Tested with Linux, Mac OSX,
and Windows):

Start Matlab or Octave and run, e.g.,

  mex -setup
  help geographiclibinterface
  geographiclibinterface
  help geodesicinverse
  geodesicinverse([40.6, -73.8, 51.6, -0.5])
  ans =

     5.1199e+01   1.0782e+02   5.5518e+06

The first command allows you to select the compiler to use (which should
be the same as that used to compile GeographicLib). On Mac OSX and
Matlab R2014b, the setup command is

  mex -setup C++

These routines just offer a simple interface to the corresponding C++
class. Use the help function to get documentation,

  help geodesicdirect

Unfortunately, the help function does not work for compiled functions in
Octave; in this case, just list the .m file, e.g.,

  type geodesicdirect
