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
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/Utility.hpp>

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
#  define GRAVITY_DEFAULT_NAME "wmm2010"
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
      for (int i = 0; i <= 2; ++i) {
        int nm[2];
        Utility::readarray<int, int, false>(coeffstr, nm, 2);
        int N = nm[0], M = nm[1];
        if (!(N >= M && M >= -1))
          throw GeographicErr("Bad degree and order " +
                              Utility::str(N) + " " + Utility::str(M));
        (i == 0 ? _C : _CC).resize(SphericalEngine::coeff::Csize(N, M));
        (i == 0 ? _S : _SC).resize(SphericalEngine::coeff::Ssize(N, M));
        Utility::readarray<double, real, false>(coeffstr, i == 0 ? _C : _CC);
        if (i == 0 && !(_C[0] == 0))
          throw GeographicErr("A degree 0 term should not be included");
        Utility::readarray<double, real, false>(coeffstr, i == 0 ? _S : _SC);
        if (i == 0)
          _gravitational = SphericalHarmonic(_C, _S, N, N, M, _amodel, _norm);
        else
          _correction = SphericalHarmonic(_CC, _SC, N, N, M, real(1), _norm);
      }
      int pos = int(coeffstr.tellg());
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in  " + coeff);
    }
  }

  void GravityModel::ReadMetadata(const std::string& name) {
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
        _amodel = Utility::num<real>(val);
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
    if (!(_amodel > 0))
      throw GeographicErr("Reference radius must be positive");
    if (int(_id.size()) != idlength_)
      throw GeographicErr("Invalid ID");
  }

  Math::real GravityModel::Geoid(real /*lat*/, real /*lon*/) const throw()
  { return 0; }
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
