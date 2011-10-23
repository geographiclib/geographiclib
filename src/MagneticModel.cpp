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
#include <iostream>
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/MagneticCircle.hpp>
#include <GeographicLib/Utility.hpp>

#define GEOGRAPHICLIB_MAGNETICMODEL_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_CPP)
RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_HPP)

#define MAGNETIC_DEFAULT_PATH "/home/ckarney/geographiclib/magnetic"

namespace GeographicLib {

  using namespace std;

  MagneticModel::MagneticModel(const std::string& name,
                               const Geocentric& earth)
    : _name(name)
    , _description("UNKNOWN")
    , _date("UNKNOWN")
    , _t0(Math::NaN<real>())
    , _tmin(Math::NaN<real>())
    , _tmax(Math::NaN<real>())
    , _a(Math::NaN<real>())
    , _hmin(Math::NaN<real>())
    , _hmax(Math::NaN<real>())
    , _N(-1)
    , _M(-1)
    , _N1(-1)
    , _M1(-1)
    , _earth(earth)
  {
    ReadMetadata(_name);
    int
      K = (_M + 1) * (2*_N - _M + 2) / 2,
      K1 = (_M1 + 1) * (2*_N1 - _M1 + 2) / 2;
    _G.resize(K); _H.resize(K-(_N+1));
    _G1.resize(K1); _H1.resize(K1-(_N1+1));
    _coeff = string(MAGNETIC_DEFAULT_PATH) + "/" + _coeff;
    {
      ifstream coeffstr(_coeff.c_str(), ios::binary);
      if (!coeffstr.good())
        throw GeographicErr("Error opening " + _coeff);
      Utility::readarray<double, real, false>(coeffstr, _G);
      Utility::readarray<double, real, false>(coeffstr, _H);
      Utility::readarray<double, real, false>(coeffstr, _G1);
      Utility::readarray<double, real, false>(coeffstr, _H1);
      int pos = coeffstr.tellg();
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in  " + _coeff);
    }
    _harma = SphericalHarmonic1(_G, _H, _N, _N, _M, _G1, _H1, _N1, _N1, _M1,
                                _a, SphericalHarmonic1::schmidt);
    _harmb = SphericalHarmonic(_G1, _H1, _N1, _N1, _M1,
                               _a, SphericalHarmonic::schmidt);

  }

  void MagneticModel::ReadMetadata(const std::string& name) {
    const char* spaces = " \t\n\v\f\r";
    _meta = string(MAGNETIC_DEFAULT_PATH) + "/" + name + ".meta";
    ifstream metastr(_meta.c_str());
    if (!metastr.good())
      throw GeographicErr("Cannot open " + _meta);
    string line;
    getline(metastr, line);
    if (!(line.size() >= 6 && line.substr(0,5) == "GMMF-"))
      throw GeographicErr(_meta + " does not contain GMMF-n signature");
    string::size_type n = line.find_first_of(spaces, 5);
    if (n != string::npos)
      n -= 5;
    string version = line.substr(5, n);
    if (version != "1")
      throw GeographicErr("Unknown version in " + _meta + ": " + version);
    string key, val;
    int N2 = -1, M2 = -1, N3 = -1, M3 = -1;
    while (getline(metastr, line)) {
      if (!ParseLine(line, key, val))
        continue;
      // Process key words
      if (key == "Name")
        _name = val;
      else if (key == "Description")
        _description = val;
      else if (key == "Radius")
        _a = Utility::readstr<real>(val);
      else if (key == "Epoch")
        _t0 = Utility::readstr<real>(val);
      else if (key == "StartTime")
        _tmin = Utility::readstr<real>(val);
      else if (key == "StopTime")
        _tmax = Utility::readstr<real>(val);
      else if (key == "MinHeight")
        _hmin =  Utility::readstr<real>(val);
      else if (key == "MaxHeight")
        _hmax =  Utility::readstr<real>(val);
      else if (key == "N")
        _N =  Utility::readstr<int>(val);
      else if (key == "M")
        _M =  Utility::readstr<int>(val);
      else if (key == "N1")
        _N1 =  Utility::readstr<int>(val);
      else if (key == "M1")
        _M1 =  Utility::readstr<int>(val);
      else if (key == "N2")
        N2 =  Utility::readstr<int>(val);
      else if (key == "M2")
        M2 =  Utility::readstr<int>(val);
      else if (key == "N3")
        N3 =  Utility::readstr<int>(val);
      else if (key == "M3")
        M3 =  Utility::readstr<int>(val);
      else if (key == "CoefficientFile")
        _coeff = val;
    }
    // Check values
    _M = _M == -1 ? _N : _M;
    _M1 = _M1 == -1 ? _N1 : _M1;
    M2  =  M2 == -1 ?  N2 :  M2;
    M3  =  M3 == -1 ?  N3 :  M3;
    if (!(N2 == -1 && M2 == -1 && N3 == -1 && M3 == -1))
      throw GeographicErr("Can only deal with linear time variation");
    if (!(_M <= _N && _M1 <= _N1 && _N1 <= _N && _M1 <= _M))
      throw GeographicErr("Degree/order mismatch");
    if (!(_a > 0))
      throw GeographicErr("Reference radius must be positive");
    if (!(_t0 > 0))
      throw GeographicErr("Epoch time not defined");
  }

  bool MagneticModel::ParseLine(const std::string& line,
                                std::string& key, std::string& val) {
    const char* spaces = " \t\n\v\f\r";
    string::size_type n0 = line.find_first_not_of(spaces);
    if (n0 == string::npos)
      return false;             // Blank line
    string::size_type n1 = line.find_first_of('#', n0);
    if (n0 == n1)
      return false;             // Only a comment
    val = line.substr(n0, n1 == string::npos ? n1 : n1 - n0);
    n0 = val.find_first_of(spaces);
    key = val.substr(0, n0);
    if (n0 == string::npos) {
      val = "";
      return true;
    }
    n0 = val.find_first_not_of(spaces, n0);
    if (n0 == string::npos) {
      val = "";
      return true;
    }
    n1 = val.find_last_not_of(spaces);
    val = val.substr(n0, n1 + 1 - n0);
    return true;
  }

  void MagneticModel::Field(real lat, real lon, real h, real t, bool diffp,
                            real& Bx, real& By, real& Bz,
                            real& Bxt, real& Byt, real& Bzt) const {
    t -= _t0;
    real x, y, z;
    real M[9];
    _earth.IntForward(lat, lon, h, x, y, z, M);
    real BX, BY, BZ;            // Components in geocentric basis
    _harma(t, x, y, z, BX, BY, BZ);
    if (diffp) {
      real BXt, BYt, BZt;
      _harmb(x, y, z, BXt, BYt, BZt);
      Bxt = - _a * (M[0] * BXt + M[3] * BYt + M[6] * BZt);
      Byt = - _a * (M[1] * BXt + M[4] * BYt + M[7] * BZt);
      Bzt = - _a * (M[2] * BXt + M[5] * BYt + M[8] * BZt);
    }
    Bx = - _a * (M[0] * BX + M[3] * BY + M[6] * BZ);
    By = - _a * (M[1] * BX + M[4] * BY + M[7] * BZ);
    Bz = - _a * (M[2] * BX + M[5] * BY + M[8] * BZ);
  }

  MagneticCircle MagneticModel::Circle(real lat, real h, real t)
    const {
    t -= _t0;
    real x, y, z;
    real M[Geocentric::dim2_];
    _earth.IntForward(lat, 0, h, x, y, z, M);
    // y = 0, cphi = M[7], sphi = M[8];

    return MagneticCircle(_a, M[7], M[8],
                          _harma.Circle(t, x, z, true),
                          _harmb.Circle(x, z, true));
  }

} // namespace GeographicLib
