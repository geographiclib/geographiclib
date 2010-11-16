/**
 * \file geoidheight.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I../include -L../lib -lGeographic geoidheight.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographicLib geoidheight.cpp

#include "GeographicLib/Geoid.hpp"
#include "mex.h"
#include <string>

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  static char rcsid[]
    = "$Id$";

  if (nrhs < 1)
    mexErrMsgTxt("One input argument required.");
  if (nrhs > 3)
    mexErrMsgTxt("More than three input arguments specified.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument can be specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("latlong coordinates are not of type double.");

  if (mxGetN(prhs[0]) != 2)
    mexErrMsgTxt("latlong coordinates must be M x 2 matrix.");

  int m = mxGetM(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, 3, mxREAL);

  double* lat = mxGetPr(prhs[0]);
  double* lon = lat + m;

  double* h = mxGetPr(plhs[0]);
  double* gradn = h + m;
  double* grade = h + 2*m;

  string geoidname("egm96-5");
  if (nrhs > 1) {
    if (!mxIsChar(prhs[1]))
      mexErrMsgTxt("geoid name must be a string.");
    if (mxGetM(prhs[1]) == 0)
      mexErrMsgTxt("geoid name cannot be empty.");
    if (mxGetM(prhs[1]) != 1)
      mexErrMsgTxt("geoid name cannot be a vector of strings.");
    int n = mxGetN(prhs[1]);
    if (n < 1)
      mexErrMsgTxt("geoid name cannot be empty.");
    mxChar* ptr = mxGetChars(prhs[1]);
    geoidname.resize(n);
    for (int i = 0; i < n; ++i)
      geoidname[i] = ptr[i];
  }
  string geoiddir("");
  if (nrhs > 2) {
    if (!mxIsChar(prhs[2]))
      mexErrMsgTxt("geoid directory must be a string.");
    if (mxGetM(prhs[2]) > 1)
      mexErrMsgTxt("geoid directory cannot be a vector of strings.");
    int n = mxGetN(prhs[2]);
    if (n > 0 && mxGetM(prhs[2]) == 1) {
      mxChar* ptr = mxGetChars(prhs[2]);
      geoiddir.resize(n);
      for (int i = 0; i < n; ++i)
        geoiddir[i] = ptr[i];
    } // else string is empty and do nothing
  }

  try {
    const Geoid g(geoidname, geoiddir);
    for (int i = 0; i < m; ++i) {
      try {
        h[i] = g(lat[i], lon[i], gradn[i], grade[i]);
      }
      catch (const std::exception& e) {
        mexWarnMsgTxt(e.what());
        h[i] = gradn[i] = grade[i] = Math::NaN();
      }
    }
  }
  catch (const std::exception& e) {
    mexErrMsgTxt(e.what());
  }
}
