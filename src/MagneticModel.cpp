/**
 * \file MagneticModel.cpp
 * \brief Implementation for GeographicLib::MagneticModel class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/MagneticModel.hpp>
#include <fstream>
#include <sstream>
#include <GeographicLib/SphericalHarmonic.hpp>
#include <iostream>

#define GEOGRAPHICLIB_MAGNETICMODEL_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_CPP)
RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_HPP)

namespace GeographicLib {

  using namespace std;

  MagneticModel::MagneticModel(const std::string& datafile,
                               const Geocentric& earth)
    : _datafile(datafile)
    , _name("UNKNOWN")
    , _date("UNKNOWN")
    , _t0(Math::NaN<real>())
    , _tmin(Math::NaN<real>())
    , _tmax(Math::NaN<real>())
    , _a(Math::NaN<real>())
    , _minh(-10e3)
    , _maxh(1000e3)
    , _N(-1)
    , _earth(earth)
  {
    ifstream data(_datafile.c_str());
    string line;
    int Nsecular = -1;
    bool header = true;         // Start reading header
    while (getline(data, line)) {
      if (line.empty())
        continue;
      char c = line[0];
      istringstream is(line.substr(1));
      switch (c) {
      case '#':
        continue;
        break;
      case '%':
        {
          // Typical line = %GeoMagRefRad: 6371200
          if (!header)
            throw GeographicErr("Header field found in body");
          std::string key;
          std::string val;
          is >> key;
          if (key == "ModelName:")
            is >> _name;
          else if (key == "ReleaseDate:")
            is >> _date;
          else if (key == "Epoch:")
            is >> _t0;
          else if (key == "ModelStartYear:")
            is >> _tmin;
          else if (key == "ModelEndYear:")
            is >> _tmax;
          else if (key == "IntStaticDeg:")
            is >> _N;
          else if (key == "IntSecVarDeg:")
            is >> Nsecular;
          else if (key == "GeoMagRefRad:")
            is >> _a;
          else if (key == "Normalization:") {
            is >> val;
            if (val != "Schmidt")
              throw GeographicErr("Unknown normalization " + val);
          } else if (key == "SpatBasFunc:") {
            is >> val;
            if (val != "Spherical")
              throw GeographicErr("Unknown basis function " + val);
          }
        }
        break;
      case 'I':
      case 'i':
        if (header) {
          if (_N == -1)
            throw GeographicErr("Model degree not specified");
          if (!(_N > 0))
            throw GeographicErr("Model degree must be positive");
          if (!(Nsecular == -1 || Nsecular == _N))
            throw GeographicErr("Wrong degree for secular terms");
          if (!(_a > 0))
            throw GeographicErr("Reference radius not positive");
          if (Math::isnan(_t0))
            throw GeographicErr("Epoch time not specified");
          int K = (_N + 1) * (_N + 2)/2;
          _G.resize(K, 0); _Gt.resize(K, 0);
          _H.resize(K, 0); _Ht.resize(K, 0);
          header = false;
        }
        {
          // Typical line (m = 0) = I,1,0,-29496.6,,11.6,
          // Typical line (m > 0) = I,1,1,-1586.3,4944.4,16.5,-25.9
          // Fields are  n m Gn,m Hn,m SV-Gn,m SV-Hn,m
          char d;
          int n = -1, m = -1;
          if (!(is >> d >> n >> d >> m))
            throw GeographicErr("Short read on a coefficient line 1");
          if (!(n >= 0 && n <= _N && m >=0 && m <= m))
            throw GeographicErr("Bad order or degree given");
          int k =  m * _N + n - m * (m - 1)/2;
          if (!(is >> d >> _G[k]))
            throw GeographicErr("Short read on a coefficient line 2");
          if (!(is >> d >> _H[k] || m == 0))
            throw GeographicErr("Short read on a coefficient line 3");
          is.clear();           // Clear error from missing _H
          if (!(is >> d >> _Gt[k]))
            throw GeographicErr("Short read on a coefficient line 4");
          if (!(is >> d >> _Ht[k] || m == 0))
            throw GeographicErr("Short read on a coefficient line 5");
        }
        break;
      case '\r':
      case '\n':
        break;
      default:
        throw GeographicErr("Bad initial character in " + _datafile);
      }
    }
    /*
    for (int n = 0; n <= _N; ++n) {
      // Change from Schmidt normalization to fully normalized (P_bar) using
      // P_Schmidt = P_bar / sqrt(2*n + 1).  Also fold in the extra factor of
      // -a in the definition of the potential.
      //      real f = -_a / sqrt(real(2 * n + 1));
      real f = 1;
      for (int m = 0, k = n; m <= n; k += _N - m++) {
        _G[k] *= f;
        _H[k] *= f;
        _Gt[k] *= f;
        _Ht[k] *= f;
      }
      }*/
  }

  void MagneticModel::Field(real lat, real lon, real h, real t, bool diffp,
                            real& Bx, real& By, real& Bz,
                            real& Bxt, real& Byt, real& Bzt) const {
    t -= _t0;
    real x, y, z;
    vector<real> M(9);
    _earth.Forward(lat, lon, h, x, y, z, M);
    real BX, BY, BZ;            // Components in geocentric basis
    SphericalHarmonic::NValue(_N, _G, _H, _N, _Gt, _Ht, t, x, y, z, _a,
                              BX, BY, BZ,
                              SphericalHarmonic::schmidt);
    if (diffp) {
      real BXt, BYt, BZt;
      SphericalHarmonic::NValue
        (_N, _Gt, _Ht, x, y, z, _a, BXt, BYt, BZt, SphericalHarmonic::schmidt);
      Bxt = - _a * (M[0] * BXt + M[3] * BYt + M[6] * BZt);
      Byt = - _a * (M[1] * BXt + M[4] * BYt + M[7] * BZt);
      Bzt = - _a * (M[2] * BXt + M[5] * BYt + M[8] * BZt);
    }
    Bx = - _a * (M[0] * BX + M[3] * BY + M[6] * BZ);
    By = - _a * (M[1] * BX + M[4] * BY + M[7] * BZ);
    Bz = - _a * (M[2] * BX + M[5] * BY + M[8] * BZ);
  }

} // namespace GeographicLib
