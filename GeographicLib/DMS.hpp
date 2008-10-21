/**
 * \file DMS.hpp
 *
 * Copyright (c) Charles Karney (2008) <charles@karney.com>
 * http://charles.karney.info/geographic
 * and licensed under the LGPL.
 **********************************************************************/

#if !defined(DMS_HPP)
#define DMS_HPP "$Id$"

#include <string>
#include <sstream>

namespace GeographicLib {

  class DMS {
  private:
    static int lookup(const std::string& s, char c) {
      std::string::size_type r = s.find(toupper(c));
      return r == std::string::npos ? -1 : int(r);
    }
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }
    static const std::string hemispheres;
    static const std::string signs;
    static const std::string digits;
    static const std::string dmsindicators;
    static const std::string components[3];

  public:
    enum flag { NONE = 0, LATITUDE = 1, LONGITUDE = 2 };
    enum component { DEGREE = 0, MINUTE = 1, SECOND = 2 };
    static double Decode(const std::string& dms, flag& ind);
    static std::string Encode(double degree,
			      component trailing = DEGREE,
			      unsigned prec = 0,
			      flag ind = NONE);
  };

} // namespace GeographicLib

#endif
