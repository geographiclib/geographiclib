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
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/MagneticCircle.hpp>
#include <GeographicLib/Utility.hpp>

#define GEOGRAPHICLIB_MAGNETICMODEL_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_CPP)
RCSID_DECL(GEOGRAPHICLIB_MAGNETICMODEL_HPP)

#if !defined(MAGNETIC_DEFAULT_PATH)
#  if defined(GEOGRAPHICLIB_MAGNETIC_PATH)
       // Use cmake supplied value is available
#      define MAGNETIC_DEFAULT_PATH GEOGRAPHICLIB_MAGNETIC_PATH
#  else
#    if defined(_MSC_VER)
#      define MAGNETIC_DEFAULT_PATH \
  "C:/Documents and Settings/All Users/Application Data/GeographicLib/magnetic"
#    else
#      define MAGNETIC_DEFAULT_PATH "/usr/local/share/GeographicLib/magnetic"
#    endif
#  endif
#endif
#if !defined(MAGNETIC_DEFAULT_NAME)
#  define MAGNETIC_DEFAULT_NAME "wmm2010"
#endif

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  MagneticModel::MagneticModel(const std::string& name,const std::string& path,
                               const Geocentric& earth)
    : _name(name)
    , _dir(path)
    , _description("NONE")
    , _date("UNKNOWN")
    , _t0(Math::NaN<real>())
    , _dt0(1)
    , _tmin(Math::NaN<real>())
    , _tmax(Math::NaN<real>())
    , _a(Math::NaN<real>())
    , _hmin(Math::NaN<real>())
    , _hmax(Math::NaN<real>())
    , _Nmodels(1)
    , _norm(SphericalHarmonic::schmidt)
    , _earth(earth)
  {
    if (_dir.empty())
      _dir = DefaultMagneticPath();
    ReadMetadata(_name);
    _G.resize(_Nmodels + 1);
    _H.resize(_Nmodels + 1);
    {
      std::string coeff = _filename + ".cof";
      ifstream coeffstr(coeff.c_str(), ios::binary);
      if (!coeffstr.good())
        throw GeographicErr("Error opening " + coeff);
      char id[idlength_ + 1];
      coeffstr.read(id, idlength_);
      if (!coeffstr.good())
        throw GeographicErr("No header in " + coeff);
      id[idlength_] = '\0';
      if (_id != std::string(id))
        throw GeographicErr("ID mismatch: " + _id + " vs " + id);
      for (int i = 0; i <= _Nmodels; ++i) {
        int nm[2];
        Utility::readarray<int, int, false>(coeffstr, nm, 2);
        int N = nm[0], M = nm[1];
        if (!(N >= M && M >= -1))
          throw GeographicErr("Bad degree and order " +
                              Utility::str(N) + " " + Utility::str(M));
        _G[i].resize(SphericalEngine::coeff::Csize(N, M));
        _H[i].resize(SphericalEngine::coeff::Ssize(N, M));
        Utility::readarray<double, real, false>(coeffstr, _G[i]);
        if (!(_G[i][0] == 0))
          throw GeographicErr("A degree 0 term is not permitted");
        Utility::readarray<double, real, false>(coeffstr, _H[i]);
        _harm.push_back(SphericalHarmonic(_G[i], _H[i], N, N, M, _a, _norm));
      }
      int pos = int(coeffstr.tellg());
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in  " + coeff);
    }
  }

  void MagneticModel::ReadMetadata(const std::string& name) {
    const char* spaces = " \t\n\v\f\r";
    _filename = _dir + "/" + name + ".wmm";
    ifstream metastr(_filename.c_str());
    if (!metastr.good())
      throw GeographicErr("Cannot open " + _filename);
    string line;
    getline(metastr, line);
    if (!(line.size() >= 6 && line.substr(0,5) == "WMMF-"))
      throw GeographicErr(_filename + " does not contain WMMF-n signature");
    string::size_type n = line.find_first_of(spaces, 5);
    if (n != string::npos)
      n -= 5;
    string version = line.substr(5, n);
    if (version != "1")
      throw GeographicErr("Unknown version in " + _filename + ": " + version);
    string key, val;
    while (getline(metastr, line)) {
      if (!Utility::ParseLine(line, key, val))
        continue;
      // Process key words
      if (key == "Name")
        _name = val;
      else if (key == "Description")
        _description = val;
      else if (key == "ReleaseDate")
        _date = val;
      else if (key == "Radius")
        _a = Utility::num<real>(val);
      else if (key == "Type") {
        if (!(val == "Linear" || val == "linear"))
          throw GeographicErr("Only linear models are supported");
      } else if (key == "Epoch")
        _t0 = Utility::num<real>(val);
      else if (key == "DeltaEpoch")
        _dt0 = Utility::num<real>(val);
      else if (key == "NumModels")
        _Nmodels = Utility::num<int>(val);
      else if (key == "MinTime")
        _tmin = Utility::num<real>(val);
      else if (key == "MaxTime")
        _tmax = Utility::num<real>(val);
      else if (key == "MinHeight")
        _hmin =  Utility::num<real>(val);
      else if (key == "MaxHeight")
        _hmax =  Utility::num<real>(val);
      else if (key == "Normalization") {
        if (val == "Full" || val == "full")
          _norm = SphericalHarmonic::full;
        else if (val == "Schmidt" || val == "schmidt")
          _norm = SphericalHarmonic::schmidt;
        else
          throw GeographicErr("Unknown normalization " + val);
      } else if (key == "ByteOrder") {
        if (val == "Big" || val == "big")
          throw GeographicErr("Only little-endian ordering is supported");
        else if (!(val == "Little" || val == "little"))
          throw GeographicErr("Unknown byte ordering " + val);
      } else if (key == "ID")
        _id = val;
      // else unrecognized keywords are skipped
    }
    // Check values
    if (!(_a > 0))
      throw GeographicErr("Reference radius must be positive");
    if (!(_t0 > 0))
      throw GeographicErr("Epoch time not defined");
    if (_tmin >= _tmax)
      throw GeographicErr("Min time exceeds max time");
    if (_hmin >= _hmax)
      throw GeographicErr("Min height exceeds max height");
    if (int(_id.size()) != idlength_)
      throw GeographicErr("Invalid ID");
    if (!(_dt0 > 0)) {
      if (_Nmodels > 1)
        throw GeographicErr("DeltaEpoch must be positive");
      else
        _dt0 = 1;
    }
  }

  void MagneticModel::Field(real t, real lat, real lon, real h, bool diffp,
                            real& Bx, real& By, real& Bz,
                            real& Bxt, real& Byt, real& Bzt) const throw() {
    t -= _t0;
    int n = max(min(int(floor(t / _dt0)), _Nmodels - 1), 0);
    bool interpolate = n + 1 < _Nmodels;
    t -= n * _dt0;
    real x, y, z;
    real M[9];
    _earth.IntForward(lat, lon, h, x, y, z, M);
    real BX0, BY0, BZ0, BX1, BY1, BZ1; // Components in geocentric basis
    _harm[n](x, y, z, BX0, BY0, BZ0);
    _harm[n + 1](x, y, z, BX1, BY1, BZ1);
    if (interpolate) {
      BX1 = (BX1 - BX0) / _dt0;
      BY1 = (BY1 - BY0) / _dt0;
      BZ1 = (BZ1 - BZ0) / _dt0;
    }
    BX0 += t * BX1;
    BY0 += t * BY1;
    BZ0 += t * BZ1;
    if (diffp) {
      Bxt = - _a * (M[0] * BX1 + M[3] * BY1 + M[6] * BZ1);
      Byt = - _a * (M[1] * BX1 + M[4] * BY1 + M[7] * BZ1);
      Bzt = - _a * (M[2] * BX1 + M[5] * BY1 + M[8] * BZ1);
    }
    Bx = - _a * (M[0] * BX0 + M[3] * BY0 + M[6] * BZ0);
    By = - _a * (M[1] * BX0 + M[4] * BY0 + M[7] * BZ0);
    Bz = - _a * (M[2] * BX0 + M[5] * BY0 + M[8] * BZ0);
  }

  MagneticCircle MagneticModel::Circle(real t, real lat, real h) const {
    t -= _t0;
    int n = max(min(int(floor(t / _dt0)), _Nmodels - 1), 0);
    bool interpolate = n + 1 < _Nmodels;
    t -= n * _dt0;
    real x, y, z;
    real M[Geocentric::dim2_];
    _earth.IntForward(lat, 0, h, x, y, z, M);
    // y = 0, cphi = M[7], sphi = M[8];

    return MagneticCircle(_a, M[7], M[8], t, _dt0, interpolate,
                          _harm[n].Circle(x, z, true),
                          _harm[n + 1].Circle(x, z, true));
  }

  void MagneticModel::FieldComponents(real Bx, real By, real Bz,
                                      real Bxt, real Byt, real Bzt,
                                      real& H, real& F, real& D, real& I,
                                      real& Ht, real& Ft,
                                      real& Dt, real& It) throw() {
    H = Math::hypot(Bx, By);
    Ht = H ? (Bx * Bxt + By * Byt) / H : Math::hypot(Bxt, Byt);
    D = (0 - (H ? atan2(-Bx, By) : atan2(-Bxt, Byt))) / Math::degree<real>();
    Dt = (H ? (By * Bxt - Bx * Byt) / Math::sq(H) : 0) / Math::degree<real>();
    F = Math::hypot(H, Bz);
    Ft = F ? (H * Ht + Bz * Bzt) / F : Math::hypot(Ht, Bzt);
    I = (F ? atan2(-Bz, H) : atan2(-Bzt, Ht)) / Math::degree<real>();
    It = (F ? (Bz * Ht - H * Bzt) / Math::sq(F) : 0) / Math::degree<real>();
  }

  std::string MagneticModel::DefaultMagneticPath() {
    string path;
    char* magneticpath = getenv("MAGNETIC_PATH");
    if (magneticpath)
      path = string(magneticpath);
    return path.length() ? path : string(MAGNETIC_DEFAULT_PATH);
  }

  std::string MagneticModel::DefaultMagneticName() {
    string name;
    char* magneticname = getenv("MAGNETIC_NAME");
    if (magneticname)
      name = string(magneticname);
    return name.length() ? name : string(MAGNETIC_DEFAULT_NAME);
  }

} // namespace GeographicLib
