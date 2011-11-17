/**
 * \file GravityModel.cpp
 * \brief Implementation for GeographicLib::GravityModel class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/GravityModel.hpp>
#include <fstream>
#include <GeographicLib/SphericalEngine.hpp>
//#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/Utility.hpp>
#include <iostream>

#define GEOGRAPHICLIB_GRAVITYMODEL_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GRAVITYMODEL_CPP)
RCSID_DECL(GEOGRAPHICLIB_GRAVITYMODEL_HPP)

#if !defined(GRAVITY_DEFAULT_PATH)
#  if defined(GEOGRAPHICLIB_GRAVITY_PATH)
       // Use cmake supplied value is available
#      define GRAVITY_DEFAULT_PATH GEOGRAPHICLIB_GRAVITY_PATH
#  else
#    if defined(_MSC_VER)
#      define GRAVITY_DEFAULT_PATH \
  "C:/Documents and Settings/All Users/Application Data/GeographicLib/gravity"
#    else
#      define GRAVITY_DEFAULT_PATH "/usr/local/share/GeographicLib/gravity"
#    endif
#  endif
#endif
#if !defined(GRAVITY_DEFAULT_NAME)
#  define GRAVITY_DEFAULT_NAME "egm96"
#endif

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  GravityModel::GravityModel(const std::string& name,const std::string& path)
    : _name(name)
    , _dir(path)
    , _description("NONE")
    , _date("UNKNOWN")
    , _amodel(Math::NaN<real>())
    , _GMmodel(Math::NaN<real>())
    , _zeta0(0)
    , _corrmult(1)
    , _norm(SphericalHarmonic::full)
  {
    if (_dir.empty())
      _dir = DefaultGravityPath();
    ReadMetadata(_name);
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
      int N, M;
      SphericalEngine::coeff::readcoeffs(coeffstr, N, M, _C, _S);
      if (!(M < 0 || _C[0] == 0))
        throw GeographicErr("A degree 0 term should be zero");
      _C[0] = 1;                // Include the 1/r term in the sum
      _gravitational = SphericalHarmonic(_C, _S, N, N, M, _amodel, _norm);
      SphericalEngine::coeff::readcoeffs(coeffstr, N, M, _CC, _CS);
      _correction = SphericalHarmonic(_CC, _CS, N, N, M, real(1), _norm);
      int pos = int(coeffstr.tellg());
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in  " + coeff);
    }
    int nmx = _gravitational.Coefficients().nmx();
    // Adjust the normalization of the normal potential to match the model.
    real mult = _earth._GM / _GMmodel;
    real amult = Math::sq(_earth._a / _amodel);
    // The 0th term in _zonal should be is 1 + _dzonal0.  Instead set it to 1
    // to give exact cancelation with the (0,0) term in the model and account
    // for _dzonal0 separately.
    _zonal.resize(0); _zonal.push_back(1);
    _dzonal0 = (_earth.GravitationalConstant() - _GMmodel) / _GMmodel;
    _dzonal0 = 0;               // FOR COMPATIBILITY WITH NGA
    for (int n = 2; n <= nmx; n += 2) {
      mult *= amult;
      real
        r = _C[n],                                         // the model term
        s = - mult * _earth.Jn(n) / sqrt(real(2 * n + 1)), // the normal term
        t = r - s;                                         // the difference
      if (t == r)               // the normal term is negligible
        break;
      _zonal.push_back(0);      // index = n - 1; the odd terms are 0
      _zonal.push_back(s);
    }
    int nmx1 = int(_zonal.size()) - 1;
    _disturbing = SphericalHarmonic1(_C, _S,
                                     _gravitational.Coefficients().N(),
                                     nmx, _gravitational.Coefficients().mmx(),
                                     _zonal,
                                     _zonal, // This is not accessed!
                                     nmx1, nmx1, 0,
                                     _amodel,
                                     SphericalHarmonic1::normalization(_norm));
  }

  void GravityModel::ReadMetadata(const std::string& name) {
    const char* spaces = " \t\n\v\f\r";
    _filename = _dir + "/" + name + ".egm";
    ifstream metastr(_filename.c_str());
    if (!metastr.good())
      throw GeographicErr("Cannot open " + _filename);
    string line;
    getline(metastr, line);
    if (!(line.size() >= 6 && line.substr(0,5) == "EGMF-"))
      throw GeographicErr(_filename + " does not contain EGMF-n signature");
    string::size_type n = line.find_first_of(spaces, 5);
    if (n != string::npos)
      n -= 5;
    string version = line.substr(5, n);
    if (version != "1")
      throw GeographicErr("Unknown version in " + _filename + ": " + version);
    string key, val;
    real a = Math::NaN<real>(), GM = a, omega = a, f = a, J2 = a;
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
      else if (key == "ModelRadius")
        _amodel = Utility::num<real>(val);
      else if (key == "ModelGravity")
        _GMmodel = Utility::num<real>(val);
      else if (key == "AngularVelocity")
        omega = Utility::num<real>(val);
      else if (key == "ReferenceRadius")
        a = Utility::num<real>(val);
      else if (key == "ReferenceGravity")
        GM = Utility::num<real>(val);
      else if (key == "Flattening")
        f = Utility::fract<real>(val);
      else if (key == "DynamicalFormFactor")
        J2 = Utility::fract<real>(val);
      else if (key == "HeightOffset")
        _zeta0 = Utility::fract<real>(val);
      else if (key == "CorrectionMultiplier")
        _corrmult = Utility::fract<real>(val);
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
    if (!(Math::isfinite(_amodel) && _amodel > 0))
      throw GeographicErr("Model radius must be positive");
    if (!(Math::isfinite(_GMmodel) && _GMmodel > 0))
      throw GeographicErr("Model gravitational constant must be positive");
    bool flatp = Math::isfinite(f);
    if (flatp && Math::isfinite(J2))
      throw GeographicErr
        ("Cannot specify both flattening and dynamical form factor");
    if (int(_id.size()) != idlength_)
      throw GeographicErr("Invalid ID");
    _earth = NormalGravity(a, GM, omega, flatp ? f : J2, flatp);
  }

  Math::real GravityModel::InternalT(real X, real Y, real Z,
                                     real& deltaX, real& deltaY, real& deltaZ,
                                     bool gradp) const throw() {
    real
      invR = _dzonal0 ? 1 / Math::hypot(Math::hypot(X, Y),  Z) : 1,
      T = (gradp
           ? _disturbing(-1, X, Y, Z, deltaX, deltaY, deltaZ)
           : _disturbing(-1, X, Y, Z));
    T = (T / _amodel - _dzonal0 * invR) * _GMmodel;
    if (gradp) {
      real f = _GMmodel / _amodel;
      deltaX *= f;
      deltaY *= f;
      deltaZ *= f;
      if (_dzonal0) {
        invR = _GMmodel * _dzonal0 * invR * invR * invR;
        deltaX += X * invR;
        deltaY += Y * invR;
        deltaZ += Z * invR;
      }
    }
    return T;
  }

  Math::real GravityModel::InternalV(real X, real Y, real Z,
                                     real& gX, real& gY, real& gZ,
                                     bool gradp) const throw() {
    real
      V = (gradp
           ? _gravitational(X, Y, Z, gX, gY, gZ)
           : _gravitational(X, Y, Z)),
      f = _GMmodel / _amodel;
    V *= f;
    if (gradp) {
      gX *= f;
      gY *= f;
      gZ *= f;
    }
    return V;
  }

  Math::real GravityModel::Geoid(real lat, real lon) const throw()
  {
    real X, Y, Z;
    _earth.Earth().IntForward(lat, lon, 0, X, Y, Z, NULL);
    real
      gamma = _earth.SurfaceGravity(lat),
      dummy,
      T = InternalT(X, Y, Z, dummy, dummy, dummy, false),
      invR = 1 / Math::hypot(Math::hypot(X, Y),  Z),
      correction = _corrmult * _correction(invR * X, invR * Y, invR * Z);
    return T/gamma + _zeta0 + correction;
  }
                                   
  Math::real GravityModel::Disturbing(real /*lat*/, real /*lon*/, real /*h*/) const throw()
  { return 0; }
  Math::real GravityModel::Disturbing(real /*lat*/, real /*lon*/, real /*h*/,
                                      real& /*gx*/, real& /*gy*/, real& /*gz*/) const throw()
  { return 0; }
  Math::real GravityModel::Gravitational(real /*lat*/, real /*lon*/, real /*h*/,
                                         real& /*gx*/, real& /*gy*/, real& /*gz*/) const throw()
  { return 0; }
  Math::real GravityModel::Normal(real /*lat*/, real /*lon*/, real /*h*/,
                                  real& /*gx*/, real& /*gy*/, real& /*gz*/) const throw()
  { return 0; }
  Math::real GravityModel::Rotational(real /*lat*/, real /*h*/,
                                      real& /*gy*/, real& /*gz*/) const throw()
  { return 0; }
  Math::real GravityModel::Total(real /*lat*/, real /*lon*/, real /*h*/,
                                 real& /*gx*/, real& /*gy*/, real& /*gz*/) const throw()
  { return 0; }

  std::string GravityModel::DefaultGravityPath() {
    string path;
    char* gravitypath = getenv("GRAVITY_PATH");
    if (gravitypath)
      path = string(gravitypath);
    return path.length() ? path : string(GRAVITY_DEFAULT_PATH);
  }

  std::string GravityModel::DefaultGravityName() {
    string name;
    char* gravityname = getenv("GRAVITY_NAME");
    if (gravityname)
      name = string(gravityname);
    return name.length() ? name : string(GRAVITY_DEFAULT_NAME);
  }

} // namespace GeographicLib
