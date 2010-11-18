/**
 * \file mgrsreverse.cpp
 * \brief Matlab mex file for UTM/UPS to MGRS conversions
 *
 * Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic mgrsreverse.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographicLib mgrsreverse.cpp

// "$Id$";

#include "GeographicLib/MGRS.hpp"
#include "mex.h"

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs != 1)
    mexErrMsgTxt("One input argument required.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument can be specified.");

  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("mgrs coordinates should be strings.");

  int n = mxGetN(prhs[0]);
  if (n < 3)
    mexErrMsgTxt("mgrs strings are too short.");

  int m = mxGetM(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, 5, mxREAL);

  mxChar* mgrs = mxGetChars(prhs[0]);

  double* x = mxGetPr(plhs[0]);
  double* y = x + m;
  double* zone = x + 2*m;
  double* hemi = x + 3*m;
  double* prec = x + 4*m;

  string mgrsstr;
  const char* spaces = " \t\n\v\f\r"; // Include comma as a space

  for (int i = 0; i < m; ++i) {
    try {
      mgrsstr.resize(n);
      for (int k = 0; k < n; ++k)
        mgrsstr[k] = mgrs[i + k * m];
      string::size_type pos1 = mgrsstr.find_first_not_of(spaces);
      if (pos1 == string::npos)
        throw GeographicErr("Empty MGRS string");
      string::size_type pos2 = mgrsstr.find_first_of(spaces, pos1);
      mgrsstr = mgrsstr.substr(pos1, pos2);
      int ZONE, PREC;
      bool HEMI;
      MGRS::Reverse(mgrsstr, ZONE, HEMI, x[i], y[i], PREC);
      zone[i] = double(ZONE);
      hemi[i] = HEMI ? 1 : 0;
      prec[i] = double(PREC);
    }
    catch (const std::exception& e) {
      mexWarnMsgTxt(e.what());
      x[i] = y[i] = zone[i] = hemi[i] = prec[i] = Math::NaN();
    }
  }
}
