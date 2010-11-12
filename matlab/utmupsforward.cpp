/**
 * \file utmupsforward.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I../include -L../lib -lGeographic utmupsforward.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographicLib utmupsforward.cpp

#include "GeographicLib/UTMUPS.hpp"
#include "mex.h"

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  static char rcsid[]
    = "$Id$";

  if (nrhs < 1)
    mexErrMsgTxt("One input argument required.");
  if (nrhs > 2)
    mexErrMsgTxt("More than two input arguments specified.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument can be specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("latlong coordinates are not of type double.");

  if (mxGetN(prhs[0]) != 2)
    mexErrMsgTxt("latlong coordinates must be M x 2 matrix.");

  int setzone;
  if (nrhs == 1)
    setzone = UTMUPS::STANDARD;
  else {
    if (!( mxIsDouble(prhs[1]) && mxIsComplex(prhs[1]) &&
           mxGetNumberOfElements(prhs[1]) == 1 ))
      mexErrMsgTxt("setzone is not an integer.");
    double rzone = mxGetScalar(prhs[1]);
    setzone = int(rzone);
    if (double(setzone) != rzone)
      mexErrMsgTxt("setzone is not an integer.");
    if (setzone < UTMUPS::MINPSEUDOZONE || setzone > UTMUPS::MAXZONE)
      mexErrMsgTxt("setzone outside the legal range.");
  }

  int m = mxGetM(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, 6, mxREAL);

  double* lat = mxGetPr(prhs[0]);
  double* lon = lat + m;

  double* x = mxGetPr(plhs[0]);
  double* y = x + m;
  double* zone = x + 2*m;
  double* hemi = x + 3*m;
  double* gamma = x + 4*m;
  double* k = x + 5*m;

  for (int i = 0; i < m; ++i) {
    int ZONE;
    bool HEMI;
    try {
      UTMUPS::Forward(lat[i], lon[i], ZONE, HEMI, x[i], y[i], gamma[i], k[i],
                      setzone);
      zone[i] = double(ZONE);
      hemi[i] = HEMI ? 1 : 0;
    }
    catch (const std::exception& e) {
      mexWarnMsgTxt(e.what());
      x[i] = y[i] = zone[i] = hemi[i] = gamma[i] = k[i] = Math::NaN();
    }
  }
}


