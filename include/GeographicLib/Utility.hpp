/**
 * \file Utility.hpp
 * \brief Header for GeographicLib::Utility class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_UTILITY_HPP)
#define GEOGRAPHICLIB_UTILITY_HPP "$Id$"

#include <GeographicLib/Constants.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

namespace GeographicLib {

  /**
   * \brief Some utility routines for %GeographicLib
   **********************************************************************/
  class GEOGRAPHIC_EXPORT Utility {
  public:
    /**
     * Convert a object of type T to a string.
     *
     * @tparam T the type of the argument.
     * @param x the value to be converted.
     * @return the string representation.
     **********************************************************************/
    template<typename T> static std::string str(T x) {
      std::ostringstream s; s << x; return s.str();
    }

    /**
     * Convert a string to an object of type T.
     *
     * @tparam T the type of the return value.
     * @param s the string to be converted.
     * @return object of type T
     **********************************************************************/
    template<typename T> static T readstr(const std::string& s) {
      T x;
      std::istringstream is(s);
      if (!(is >> x))
        throw GeographicErr("Cannot decode " + s);
      int pos = int(is.tellg()); // Returns -1 at end of string?
      if (!(pos < 0 || pos == int(s.size())))
        throw GeographicErr("Extra text at end of " + s);
      return x;
    }

    /**
     * Lookup up a character in a string
     *
     * @param s the string to be searched.
     * @param c the character to look for.
     * @return the index of the first occurrence character in the string or -1
     *   is the character is not present.
     **********************************************************************/
    static int lookup(const std::string& s, char c) throw() {
      std::string::size_type r = s.find(toupper(c));
      return r == std::string::npos ? -1 : int(r);
    }

    /**
     * Read data of type ExtT from a binary stream to an array of type IntT.
     * The data in the file is in (bigendp ? big : little)-endian format.
     *
     * @tparam ExtT the type of the objects in the binary stream (external).
     * @tparam IntT the type of the objects in the array (internal).
     * @tparam bigendp true if the external storage format is big-endian.
     * @param[in] str the input stream containing the data of type ExtT
     *   (external).
     * @param[out] array the output array of type IntT (internal).
     **********************************************************************/
    template<typename ExtT, typename IntT, bool bigendp>
      static inline void readarray(std::istream& str,
                                   std::vector<IntT>& array) {
      if (sizeof(IntT) == sizeof(ExtT) &&
          std::numeric_limits<IntT>::is_integer ==
          std::numeric_limits<ExtT>::is_integer) {
        // Data is compatible (aside from the issue of endian-ness).
        str.read(reinterpret_cast<char *>(&array[0]),
                 array.size() * sizeof(IntT));
        if (!str.good())
          throw GeographicErr("Failure reading data");
        if (bigendp != Math::bigendian) { // endian mismatch -> swap bytes
          for (int i = array.size(); i--;)
            array[i] = Math::swab<IntT>(array[i]);
        }
      } else {
        const int bufsize = 1024; // read this many values at a time
        ExtT buffer[bufsize];     // temporary buffer
        int k = array.size();     // data values left to read
        int i = 0;                // index into output array
        while (k) {
          int num = (std::min)(k, bufsize);
          str.read(reinterpret_cast<char *>(buffer), num * sizeof(IntT));
          if (!str.good())
            throw GeographicErr("Failure reading data");
          for (int j = 0; j < num; ++j)
            // fix endian-ness and cast to IntT
            array[i++] = IntT(bigendp == Math::bigendian ? buffer[j] :
                              Math::swab<IntT>(buffer[j]));
          k -= num;
        }
      }
      return;
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_UTILITY_HPP
