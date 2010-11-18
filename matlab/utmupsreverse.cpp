/**
 * \file utmupsreverse.cpp
 * \brief Matlab mex file for UTM/UPS to geographic conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic utmupsreverse.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographicLib utmupsreverse.cpp

// "$Id$";

#include "GeographicLib/UTMUPS.hpp"
#include "mex.h"

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs != 1)
    mexErrMsgTxt("One input argument required.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument can be specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("utmups coordinates are not of type double.");

  if (mxGetN(prhs[0]) != 4)
    mexErrMsgTxt("utmups coordinates must be M x 4 matrix.");

  int m = mxGetM(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, 4, mxREAL);

  double* x = mxGetPr(prhs[0]);
  double* y = x + m;
  double* zone = x + 2*m;
  double* hemi = x + 3*m;

  double* lat = mxGetPr(plhs[0]);
  double* lon = lat + m;
  double* gamma = lat + 2*m;
  double* k = lat + 3*m;

  for (int i = 0; i < m; ++i) {
    try {
      int ZONE = int(zone[i]);
      if (double(ZONE) != zone[i])
        throw GeographicErr("Zone is not an integer");
      bool HEMI = (hemi[i] != 0);
      if (HEMI && (hemi[i] != 1))
        throw GeographicErr("Hemisphere is not 0 or 1");
      UTMUPS::Reverse(ZONE, HEMI, x[i], y[i], lat[i], lon[i], gamma[i], k[i]);
    }
    catch (const std::exception& e) {
      mexWarnMsgTxt(e.what());
      lat[i] = lon[i] = gamma[i] = k[i] = Math::NaN();
    }
  }
}
