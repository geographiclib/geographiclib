/**
 * \file Utility.cpp
 * \brief Implementation for GeographicLib::Utility class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <cstdlib>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#  pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  bool Utility::ParseLine(const string& line, string& key, string& val) {
    const char* spaces = " \t\n\v\f\r";
    key = "";
    string::size_type n = line.find('#');
    val = trim(line.substr(0, n));
    if (val.empty())
      return false;
    n = val.find("=");
    if (n == string::npos) n = val.find_first_of(spaces);
    key = val.substr(0, n);
    if (key.empty()) {
      val = "";
      return false;
    }
    val = trim(val.substr(n + 1));
    return true;
  }

  int Utility::set_digits(int ndigits) {
#if GEOGRAPHICLIB_PRECISION == 5
    if (ndigits <= 0) {
      char* digitenv = getenv("GEOGRAPHICLIB_DIGITS");
      if (digitenv)
        ndigits = strtol(digitenv, NULL, 0);
      if (ndigits <= 0)
        ndigits = 256;
    }
#endif
    return Math::set_digits(ndigits);
  }

} // namespace GeographicLib
