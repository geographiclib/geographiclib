/**
 * \file geodesicinverse.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic geodesicinverse.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographic geodesicinverse.cpp

// "$Id$";

#include "GeographicLib/Geodesic.hpp"
#include "mex.h"

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs < 1)
    mexErrMsgTxt("One input argument required.");
  if (nrhs > 3)
    mexErrMsgTxt("More than three input arguments specified.");
  if (nrhs == 2)
    mexErrMsgTxt("Must specify repicrocal flattening with the major radius.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument can be specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("latlong coordinates are not of type double.");

  if (mxGetN(prhs[0]) != 4)
    mexErrMsgTxt("latlong coordinates must be M x 4 matrix.");

  double a = Constants::WGS84_a(), r = Constants::WGS84_r();
  if (nrhs == 3) {
    if (!( mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) &&
           mxGetNumberOfElements(prhs[1]) == 1 ))
      mexErrMsgTxt("major radius is not a real scalar.");
    a = mxGetScalar(prhs[1]);
    if (!( mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) &&
           mxGetNumberOfElements(prhs[2]) == 1 ))
      mexErrMsgTxt("reciprocal flattening is not a real scalar.");
    r = mxGetScalar(prhs[2]);
  }

  int m = mxGetM(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, 7, mxREAL);

  double* lat1 = mxGetPr(prhs[0]);
  double* lon1 = lat1 + m;
  double* lat2 = lat1 + 2*m;
  double* lon2 = lat1 + 3*m;

  double* azi1 = mxGetPr(plhs[0]);
  double* azi2 = azi1 + m;
  double* s12 = azi1 + 2*m;
  double* m12 = azi1 + 3*m;
  double* M12 = azi1 + 4*m;
  double* M21 = azi1 + 5*m;
  double* S12 = azi1 + 6*m;

  try {
    const Geodesic g(a, r);
    for (int i = 0; i < m; ++i) {
      try {
        if (abs(lat1[i]) > 90 || abs(lat2[i]) > 90)
          throw GeographicErr("Invalid latitude");
        if (lon1[i] < -180 || lon1[i] > 360 || lon2[i] < -180 || lon2[i] > 360)
          throw GeographicErr("Invalid longitude");
        g.Inverse(lat1[i], lon1[i], lat2[i], lon2[i],
                  s12[i], azi1[i], azi2[i], m12[i], M12[i], M21[i], S12[i]);
      }
      catch (const std::exception& e) {
        mexWarnMsgTxt(e.what());
        s12[i] = azi1[i] = azi2[i]
          = m12[i] = M12[i] = M21[i] = S12[i] = Math::NaN();
      }
    }
  }
  catch (const std::exception& e) {
    mexErrMsgTxt(e.what());
  }
}
