/**
 * \file localcartesianforward.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic localcartesianforward.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographic localcartesianforward.cpp

// "$Id$";

#include "GeographicLib/LocalCartesian.hpp"
#include "mex.h"

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs < 2)
    mexErrMsgTxt("Two input arguments required.");
  else if (nrhs > 4)
    mexErrMsgTxt("More than four input arguments specified.");
  else if (nrhs == 3)
    mexErrMsgTxt("Must specify repicrocal flattening with the major radius.");
  else if (nlhs > 2)
    mexErrMsgTxt("More than two output arguments specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("origin is not of type double.");
  if (!( mxGetM(prhs[0]) == 1 &&
         (mxGetN(prhs[0]) == 2 || mxGetN(prhs[0]) == 3) ))
    mexErrMsgTxt("origin be 1 x 3 or 1 x 2 matrix.");

  double* origin = mxGetPr(prhs[0]);
  double lat0 = origin[0], lon0 = origin[1],
    h0 = mxGetN(prhs[0]) == 3 ? origin[2] : 0;

  if (!( mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) ))
    mexErrMsgTxt("geodetic coordinates are not of type double.");

  if (mxGetN(prhs[1]) != 3)
    mexErrMsgTxt("geodetic coordinates must be M x 3 matrix.");

  double a = Constants::WGS84_a(), r = Constants::WGS84_r();
  if (nrhs == 4) {
    if (!( mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) &&
           mxGetNumberOfElements(prhs[2]) == 1 ))
      mexErrMsgTxt("major radius is not a real scalar.");
    a = mxGetScalar(prhs[2]);
    if (!( mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) &&
           mxGetNumberOfElements(prhs[3]) == 1 ))
      mexErrMsgTxt("reciprocal flattening is not a real scalar.");
    r = mxGetScalar(prhs[3]);
  }

  int m = mxGetM(prhs[1]);

  double* lat = mxGetPr(prhs[1]);
  double* lon = lat + m;
  double* h = lat + 2*m;

  plhs[0] = mxCreateDoubleMatrix(m, 3, mxREAL);
  double* x = mxGetPr(plhs[0]);
  double* y = x + m;
  double* z = x + 2*m;
  double* rot = NULL;
  bool rotp = nlhs == 2;

  if (rotp) {
    plhs[1] = mxCreateDoubleMatrix(m, 9, mxREAL);
    rot = mxGetPr(plhs[1]);
  }

  try {
    std::vector<double> rotv(rotp ? 9 : 0);
    const Geocentric c(a, r);
    if (abs(lat0) > 90)
      throw GeographicErr("Invalid latitude");
    if (lon0 < -180 || lon0 > 360)
      throw GeographicErr("Invalid longitude");
    const LocalCartesian l(lat0, lon0, h0, c);
    for (int i = 0; i < m; ++i) {
      try {
        if (abs(lat[i]) > 90)
          throw GeographicErr("Invalid latitude");
        if (lon[i] < -180 || lon[i] > 360)
          throw GeographicErr("Invalid longitude");
        l.Forward(lat[i], lon[i], h[i], x[i], y[i], z[i], rotv);
        if (rotp) {
          for (int k = 0; k < 9; ++k)
            rot[m * k + i] = rotv[k];
        }
      }
      catch (const std::exception& e) {
        mexWarnMsgTxt(e.what());
        x[i] = y[i] = z[i] = Math::NaN();
        if (rotp) {
          for (int k = 0; k < 9; ++k)
            rot[m * k + i] = Math::NaN();
        }
      }
    }
  }
  catch (const std::exception& e) {
    mexErrMsgTxt(e.what());
  }
}
