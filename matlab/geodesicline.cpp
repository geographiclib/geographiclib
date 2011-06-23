/**
 * \file geodesicline.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic geodesicline.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographic geodesicline.cpp

// $Id$

#include <GeographicLib/GeodesicLine.hpp>
#include <mex.h>

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs < 4)
    mexErrMsgTxt("Four input arguments required.");
  else if (nrhs > 6)
    mexErrMsgTxt("More than three input arguments specified.");
  else if (nrhs == 5)
    mexErrMsgTxt("Must specify flattening with the major radius.");
  else if (nlhs > 2)
    mexErrMsgTxt("More than two output arguments specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) &&
         mxGetNumberOfElements(prhs[0]) == 1 ))
    mexErrMsgTxt("lat1 is not a real scalar.");
  double lat1 = mxGetScalar(prhs[0]);

  if (!( mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) &&
         mxGetNumberOfElements(prhs[1]) == 1 ))
    mexErrMsgTxt("lon1 is not a real scalar.");
  double lon1 = mxGetScalar(prhs[1]);

  if (!( mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) &&
         mxGetNumberOfElements(prhs[2]) == 1 ))
    mexErrMsgTxt("azi1 is not a real scalar.");
  double azi1 = mxGetScalar(prhs[2]);

  if (!( mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) ))
    mexErrMsgTxt("distances are not of type double.");

  if (mxGetN(prhs[3]) != 1)
    mexErrMsgTxt("distances must be M x 1 matrix.");

  double a = Constants::WGS84_a(), f = Constants::WGS84_f();
  if (nrhs == 6) {
    if (!( mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4]) &&
           mxGetNumberOfElements(prhs[4]) == 1 ))
      mexErrMsgTxt("major radius is not a real scalar.");
    a = mxGetScalar(prhs[4]);
    if (!( mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5]) &&
           mxGetNumberOfElements(prhs[5]) == 1 ))
      mexErrMsgTxt("flattening is not a real scalar.");
    f = mxGetScalar(prhs[5]);
  }

  int m = mxGetM(prhs[3]);

  double* s12 = mxGetPr(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix(m, 3, mxREAL);
  double* lat2 = mxGetPr(plhs[0]);
  double* lon2 = lat2 + m;
  double* azi2 = lat2 + 2*m;
  double* m12 = NULL;
  double* M12 = NULL;
  double* M21 = NULL;
  double* S12 = NULL;
  bool aux = nlhs == 2;

  if (aux) {
    plhs[1] = mxCreateDoubleMatrix(m, 4, mxREAL);
    m12 = mxGetPr(plhs[1]);
    M12 = m12 + m;
    M21 = m12 + 2*m;
    S12 = m12 + 3*m;
  }

  try {
    const Geodesic g(a, f);
    if (abs(lat1) > 90)
      throw GeographicErr("Invalid latitude");
    if (lon1 < -180 || lon1 > 360)
      throw GeographicErr("Invalid longitude");
    if (azi1 < -180 || azi1 > 360)
      throw GeographicErr("Invalid azimuth");
    const GeodesicLine l(g, lat1, lon1, azi1);
    for (int i = 0; i < m; ++i)
      if (aux)
        l.Position(s12[i], lat2[i], lon2[i], azi2[i],
                   m12[i], M12[i], M21[i], S12[i]);
      else
        l.Position(s12[i], lat2[i], lon2[i], azi2[i]);
  }
  catch (const std::exception& e) {
    mexErrMsgTxt(e.what());
  }
}
