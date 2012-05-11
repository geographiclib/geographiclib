/**
 * \file Ellipsoid.hpp
 * \brief Header for GeographicLib::Ellipsoid class
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ELLIPSOID_HPP)
#define GEOGRAPHICLIB_ELLIPSOID_HPP \
  "$Id$"

#include <string>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/AlbersEqualArea.hpp>

namespace GeographicLib {

  /**
   * \brief Conversions for ellipsoids
   *
   * Example of use:
   * \include example-Ellipsoid.cpp
   **********************************************************************/

  class GEOGRAPHIC_EXPORT Ellipsoid {
  private:
    typedef Math::real real;
    static const int numit_ = 10;
    real _a, _f, _f1, _f12, _e2, _e12, _n, _b, _stol;
    const TransverseMercatorExact _tm;
    const AlbersEqualArea _au;
    static real tand(real x) throw() {
      return
        std::abs(x) == real(90) ? (x < 0 ?
                                   - TransverseMercatorExact::overflow_
                                   : TransverseMercatorExact::overflow_) :
        std::tan(x * Math::degree<real>());
    }
    static real atand(real x) throw()
    { return std::atan(x) / Math::degree<real>(); }

  public:

    Ellipsoid(real a, real f);
    Math::real ConformalLatitude(real phi) const throw();
    Math::real InverseConformalLatitude(real chi) const throw();
    Math::real ParametricLatitude(real phi) const throw();
    Math::real InverseParametricLatitude(real beta) const throw();
    Math::real GeocentricLatitude(real phi) const throw();
    Math::real InverseGeocentricLatitude(real theta) const throw();
    Math::real AuthalicLatitude(real phi) const throw();
    Math::real InverseAuthalicLatitude(real xi) const throw();
    Math::real RectifyingLatitude(real phi) const throw();
    Math::real InverseRectifyingLatitude(real mu) const throw();
    Math::real IsometricLatitude(real phi) const throw();
    Math::real InverseIsometricLatitude(real psi) const throw();
    Math::real QuarterMeridian() const throw();
    Math::real TransverseRadius(real phi) const throw();
    Math::real CircleRadius(real phi) const throw();
    Math::real MeridionalRadius(real phi) const throw();
    Math::real MeridianDistance(real phi) const throw();
    Math::real Volume() const throw()
    { return (4 * Math::pi<real>()) * Math::sq(_a) * _b / 3; }
    Math::real Area() const throw() {
      return 4 * Math::pi<real>() *
        ((Math::sq(_a) + Math::sq(_b) *
          (_e2 == 0 ? 1 :
           (_e2 > 0 ? Math::atanh(sqrt(_e2)) : atan(sqrt(-_e2))) /
           sqrt(abs(_e2))))/2);
    }
    Math::real EccentricitySq() { return _e2; }
    Math::real SecondEccentricitySq() { return _e12; }
    Math::real ThirdEccentricitySq() { return _e2 / (2 - _e2); }
    Math::real Flattening() { return _f; }
    Math::real SecondFlattening() { return _f / (1 - _f); }
    Math::real ThirdFlattening() { return _n; }

  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ELLIPSOID_HPP
