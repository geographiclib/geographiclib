/**
 * \file GeodesicLine.hpp
 * \brief Header for GeographicLib::GeodesicLine class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESICLINE_HPP)
#define GEOGRAPHICLIB_GEODESICLINE_HPP "$Id: GeodesicLine.hpp 6867 2010-09-11 13:04:26Z karney $"

#include "GeographicLib/Constants.hpp"
#include "GeographicLib/Geodesic.hpp"

namespace GeographicLib {

  /**
   * \brief A geodesic line.
   *
   * GeodesicLine facilitates the determination of a series of points on a
   * single geodesic.  The starting point (\e lat1, \e lon1) and the azimuth \e
   * azi1 are specified in the constructor.  GeodesicLine.Position returns the
   * location of point 2 a distance \e s12 along the geodesic.  Alternatively
   * GeodesicLine.ArcPosition gives the position of point 2 an arc length \e
   * a12 along the geodesic.  An example of use of this class is:
   \code
   #include <iostream>
   #include <iomanip>
   #include <cmath>
   #include "GeographicLib/Geodesic.hpp"
   #include "GeographicLib/GeodesicLine.hpp"

   int main() {
     // Print waypoints between JFK and SIN at
     // approximately 100km intervals.
     double
       lat1 = 40.640, lon1 = -73.779, // JFK
       lat2 = 1.359, lon2 = 177.486;  // SIN
     const GeographicLib::Geodesic&
       g = GeographicLib::Geodesic::WGS84;
     double azi1, azi2,
       a12 = g.Inverse(lat1, lon1, lat2, lon2, azi1, azi2);
     double ds = 100e3;  // Nominal distance between points = 100 km
     int num = std::ceil(a12/ (90 * ds / 10e6)); // 90 deg = 10e6 m
     double da = a12/num;              // Arc length between points
     const GeographicLib::GeodesicLine l(g, lat1, lon1, azi1);
     std::cout << std::fixed << std::setprecision(3);
     for (int i = 0; i <= num; ++i) {
       double lat, lon;
       l.ArcPosition(i * da, lat, lon);
       std::cout << lat << " " << lon << "\n";
     }
   }
   \endcode
   * The default copy constructor and assignment operators work with this
   * class, so that, for example, the previous example could start
   \code
       GeographicLib::GeodesicLine l;
       l = g.Line(lat1, lon1, azi1);
   \endcode
   * Similarly, a vector can be used to hold GeodesicLine objects.
   *
   * The calculations are accurate to better than 12 nm.  (See \ref geoderrors
   * for details.)
   **********************************************************************/

  class GeodesicLine {
  private:
    typedef Math::real real;
    friend class Geodesic;
    static const int nC1 = Geodesic::nC1, nC1p = Geodesic::nC1p,
      nC2 = Geodesic::nC2, nC3 = Geodesic::nC3, nC4 = Geodesic::nC4;

    real _lat1, _lon1, _azi1;
    real _a, _r, _b, _c2, _f1, _salp0, _calp0, _k2,
      _salp1, _calp1, _ssig1, _csig1, _stau1, _ctau1, _somg1, _comg1,
      _A1m1, _A2m1, _A3c, _B11, _B21, _B31, _A4, _B41;
    // index zero elements of _C1a, _C1pa, _C2a, _C3a are unused
    real _C1a[nC1 + 1], _C1pa[nC1p + 1], _C2a[nC2 + 1], _C3a[nC3],
      _C4a[nC4];    // all the elements of _C4a are used
    bool _areap;
    unsigned _caps;

    static inline real sq(real x) throw() { return x * x; }
    enum captype {
      CAP_NONE = Geodesic::CAP_NONE,
      CAP_C1   = Geodesic::CAP_C1,
      CAP_C1p  = Geodesic::CAP_C1p,
      CAP_C2   = Geodesic::CAP_C2,
      CAP_C3   = Geodesic::CAP_C3,
      CAP_C4   = Geodesic::CAP_C4,
      CAP_ALL  = Geodesic::CAP_ALL,
      OUT_ALL  = Geodesic::OUT_ALL,
    };
  public:

    /**
     * Bit masks for what calculations to do.  These masks do double duty.
     * They signify to the GeodesicLine::GeodesicLine constructor and to
     * Geodesic::Line what capabilities should be included in the GeodesicLine
     * object.  They also specify which results to return in the general
     * routines Geodesic::GenDirect and Geodesic::GenInverse routines.  This is
     * merely a duplication of Geodesic::mask.
     **********************************************************************/
    enum mask {
      /**
       * No capabilities, no output.
       * @hideinitializer
       **********************************************************************/
      NONE          = Geodesic::NONE,
      /**
       * Calculate latitude \e lat2.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       * @hideinitializer
       **********************************************************************/
      LATITUDE      = Geodesic::LATITUDE,
      /**
       * Calculate longitude \e lon2.
       * @hideinitializer
       **********************************************************************/
      LONGITUDE     = Geodesic::LONGITUDE,
      /**
       * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
       * include this as a capability to GeodesicLine because this is included
       * by default.)
       * @hideinitializer
       **********************************************************************/
      AZIMUTH       = Geodesic::AZIMUTH,
      /**
       * Calculate distance \e s12.
       * @hideinitializer
       **********************************************************************/
      DISTANCE      = Geodesic::DISTANCE,
      /**
       * Allow distance \e s12 to be used as input in the direct geodesic
       * problem.
       * @hideinitializer
       **********************************************************************/
      DISTANCE_IN   = Geodesic::DISTANCE_IN,
      /**
       * Calculate reduced length \e m12.
       * @hideinitializer
       **********************************************************************/
      REDUCEDLENGTH = Geodesic::REDUCEDLENGTH,
      /**
       * Calculate geodesic scales \e M12 and \e M21.
       * @hideinitializer
       **********************************************************************/
      GEODESICSCALE = Geodesic::GEODESICSCALE,
      /**
       * Calculate area \e S12.
       * @hideinitializer
       **********************************************************************/
      AREA          = Geodesic::AREA,
      /**
       * All capabilities.  Calculate everything.
       * @hideinitializer
       **********************************************************************/
      ALL           = Geodesic::ALL,
    };

    /** \name Constructors
     **********************************************************************/
    ///@{

    /**
     * Constructor for a geodesic line staring at latitude \e lat1, longitude
     * \e lon1, and aziumuth \e azi1 (all in degrees).
     *
     * @param[in] g A Geodesic object used to compute the necessary information
     *   about the GeodesicLine.
     *
     * @param[in] lat1 latitude of point 1 (degrees).
     * @param[in] lon1 longitude of point 1 (degrees).
     * @param[in] azi1 azimuth at point 1 (degrees).
     * @param[in] caps bitor'ed combination of GeodesicLine::mask values
     *   specifying the capabilities the GeodesicLine object should possess,
     *   i.e., which quantities can be returned in calls to
     *   GeodesicLib::Position.
     *
     * The GeodesicLine::mask values are
     * - \e caps |= GeodesicLine::LATITUDE for the latitude \e lat2; this is
     *   added automatically
     * - \e caps |= GeodesicLine::LONGITUDE for the latitude \e lon2
     * - \e caps |= GeodesicLine::AZIMUTH for the latitude \e azi2; this is
     *   added automatically
     * - \e caps |= GeodesicLine::DISTANCE for the distance \e s12
     * - \e caps |= GeodesicLine::REDUCEDLENGTH for the reduced length \e m12
     * - \e caps |= GeodesicLine::GEODESICSCALE for the geodesic scales \e M12
     *   and \e M21
     * - \e caps |= GeodesicLine::AREA for the area \e S12
     * - \e caps |= GeodesicLine::DISTANCE_IN permits the length of the
     *   geodesic to be given in terms of \e s12; without this capability the
     *   length can only be specified in terms of arc length.
     * .
     * The default value of \e caps is GeodesicLine::ALL which turns on all the
     * capabilities.
     *
     * If the point is at a pole, the azimuth is defined by keeping the \e lon1
     * fixed and writing \e lat1 = 90 - \e eps or -90 + \e eps and taking the
     * limit \e eps -> 0 from above.
     **********************************************************************/
    GeodesicLine(const Geodesic& g, real lat1, real lon1, real azi1,
                 unsigned caps = ALL)
      throw();

    /**
     * A default constructor.  If GeodesicLine::Position is called on the
     * resulting object, it returns immediately (without doing any
     * calculations).  The object can be set with a call to Geodesic::Line.
     * Use Init() to test whether object is still in this uninitialized state.
     **********************************************************************/
    GeodesicLine() throw() : _caps(0U) {};
    ///@}

    /** \name Position in terms of distance
     **********************************************************************/
    ///@{

    /**
     * Compute the position of point 2 which is a distance \e s12 (meters)
     * from point 1.
     *
     * @param[in] s12 distance between point 1 and point 2 (meters); it can be
     *   signed.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::LONGITUDE.
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] m12 reduced length of geodesic (meters); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::REDUCEDLENGTH.
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>); requires
     *   that the GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::AREA.
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * The GeodesicLine object \e must have been constructed with \e caps |=
     * GeodesicLine::DISTANCE_IN; otherwise Math::NaN() is returned and no
     * parameters are set.  Requesting a value which the GeodesicLine object is
     * not capable of computing is not an error; the corresponding argument
     * will not be altered.
     *
     * The following functions are overloaded versions of
     * GeodesicLine::Position which omit some of the output parameters.  Note,
     * however, that the arc length is always computed and returned as the
     * function value.
     **********************************************************************/
    Math::real Position(real s12,
                        real& lat2, real& lon2, real& azi2,
                        real& m12, real& M12, real& M21,
                        real& S12) const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE | AZIMUTH |
                         REDUCEDLENGTH | GEODESICSCALE | AREA,
                         lat2, lon2, azi2, t, m12, M12, M21, S12);
    }

    /**
     * See the documentation for GeodesicLine::Position.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2) const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE,
                         lat2, lon2, t, t, t, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::Position.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2) const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE | AZIMUTH,
                         lat2, lon2, azi2, t, t, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::Position.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2, real& m12) const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE |
                         AZIMUTH | REDUCEDLENGTH,
                         lat2, lon2, azi2, t, m12, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::Position.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2, real& M12, real& M21)
      const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE |
                         AZIMUTH | GEODESICSCALE,
                         lat2, lon2, azi2, t, t, M12, M21, t);
    }

    /**
     * See the documentation for GeodesicLine::Position.
     **********************************************************************/
    Math::real Position(real s12,
                        real& lat2, real& lon2, real& azi2,
                        real& m12, real& M12, real& M21)
      const throw() {
      real t;
      return GenPosition(false, s12,
                         LATITUDE | LONGITUDE | AZIMUTH |
                         REDUCEDLENGTH | GEODESICSCALE,
                         lat2, lon2, azi2, t, m12, M12, M21, t);
    }

    ///@}

    /** \name Position in terms of arc length
     **********************************************************************/
    ///@{

    /**
     * Compute the position of point 2 which is an arc length \e a12 (degrees)
     * from point 1.
     *
     * @param[in] a12 arc length between point 1 and point 2 (degrees); it can
     *   be signed.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::LONGITUDE.
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] s12 distance between point 1 and point 2 (meters); requires
     *   that the GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::DISTANCE.
     * @param[out] m12 reduced length of geodesic (meters); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::REDUCEDLENGTH.
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>); requires
     *   that the GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::AREA.
     *
     * Requesting a value which the GeodesicLine object is not capable of
     * computing is not an error; the corresponding argument will not be
     * altered.
     *
     * The following functions are overloaded versions of
     * GeodesicLine::ArcPosition which omit some of the output parameters.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12, real& M12, real& M21,
                     real& S12) const throw() {
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH | DISTANCE |
                  REDUCEDLENGTH | GEODESICSCALE | AREA,
                  lat2, lon2, azi2, s12, m12, M12, M21, S12);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2)
      const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE,
                  lat2, lon2, t, t, t, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12,
                     real& lat2, real& lon2, real& azi2)
      const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH,
                  lat2, lon2, azi2, t, t, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12) const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH | DISTANCE,
                  lat2, lon2, azi2, s12, t, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12) const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH |
                  DISTANCE | REDUCEDLENGTH,
                  lat2, lon2, azi2, s12, m12, t, t, t);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& M12, real& M21)
      const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH |
                  DISTANCE | GEODESICSCALE,
                  lat2, lon2, azi2, s12, t, M12, M21, t);
    }

    /**
     * See the documentation for GeodesicLine::ArcPosition.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12, real& M12, real& M21)
      const throw() {
      real t;
      GenPosition(true, a12,
                  LATITUDE | LONGITUDE | AZIMUTH |
                  DISTANCE | REDUCEDLENGTH | GEODESICSCALE,
                  lat2, lon2, azi2, s12, m12, M12, M21, t);
    }
    ///@}

    /** \name The general position function.
     **********************************************************************/
    ///@{

    /**
     * The general position function.  GeodesicLine::Position and
     * GeodesicLine::ArcPosition are defined in terms of this function.
     *
     * @param[in] arcmode boolean flag determining the meaning of the second
     *   parameter; if arcmode is false, then the GeodesicLine object must have
     *   been constructed with \e caps |= GeodesicLine::DISTANCE_IN.
     * @param[in] s12_a12 if \e arcmode is false, this is the distance between
     *   point 1 and point 2 (meters); otherwise it is the arc length between
     *   point 1 and point 2 (degrees); it can be signed.
     * @param[in] outmask a bitor'ed combination of GeodesicLine::mask values
     *   specifying which of the following parameters should be set.
     * @param[out] lat2 latitude of point 2 (degrees).
     * @param[out] lon2 longitude of point 2 (degrees); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::LONGITUDE.
     * @param[out] azi2 (forward) azimuth at point 2 (degrees).
     * @param[out] s12 distance between point 1 and point 2 (meters); requires
     *   that the GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::DISTANCE.
     * @param[out] m12 reduced length of geodesic (meters); requires that the
     *   GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::REDUCEDLENGTH.
     * @param[out] M12 geodesic scale of point 2 relative to point 1
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] M21 geodesic scale of point 1 relative to point 2
     *   (dimensionless); requires that the GeodesicLine object was constructed
     *   with \e caps |= GeodesicLine::GEODESICSCALE.
     * @param[out] S12 area under the geodesic (meters<sup>2</sup>); requires
     *   that the GeodesicLine object was constructed with \e caps |=
     *   GeodesicLine::AREA.
     * @return \e a12 arc length of between point 1 and point 2 (degrees).
     *
     * The GeodesicLine::mask values possible for \e outmask are
     * - \e outmask |= GeodesicLine::LATITUDE for the latitude \e lat2.
     * - \e outmask |= GeodesicLine::LONGITUDE for the latitude \e lon2.
     * - \e outmask |= GeodesicLine::AZIMUTH for the latitude \e azi2.
     * - \e outmask |= GeodesicLine::DISTANCE for the distance \e s12.
     * - \e outmask |= GeodesicLine::REDUCEDLENGTH for the reduced length \e
     *   m12.
     * - \e outmask |= GeodesicLine::GEODESICSCALE for the geodesic scales \e
     *   M12 and \e M21.
     * - \e outmask |= GeodesicLine::AREA for the area \e S12.
     * .
     * Requesting a value which the GeodesicLine object is not capable of
     * computing is not an error; the corresponding argument will not be
     * altered.  Note, however, that the arc length is always computed and
     * returned as the function value.
     **********************************************************************/
    Math::real GenPosition(bool arcmode, real s12_a12, unsigned outmask,
                           real& lat2, real& lon2, real& azi2,
                           real& s12, real& m12, real& M12, real& M21,
                           real& S12) const throw();

    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{

    /**
     * @return true if the object has been initialized.
     **********************************************************************/
    bool Init() const throw() { return _caps != 0U; }

    /**
     * @return \e lat1 the latitude of point 1 (degrees).
     **********************************************************************/
    Math::real Latitude() const throw() { return Init() ? _lat1 : 0; }

    /**
     * @return \e lon1 the longitude of point 1 (degrees).
     **********************************************************************/
    Math::real Longitude() const throw() { return Init() ? _lon1 : 0; }

    /**
     * @return \e azi1 the azimuth (degrees) of the geodesic line at point 1.
     **********************************************************************/
    Math::real Azimuth() const throw() { return Init() ? _azi1 : 0; }

    /**
     * @return \e azi0 the azimuth (degrees) of the geodesic line as it crosses
     * the equator in a northward direction.
     **********************************************************************/
    Math::real EquatorialAzimuth() const throw() {
      return Init() ? atan2(_salp0, _calp0) / Constants::degree() : 0;
    }

    /**
     * @return \e a1 the arc length (degrees) between the northward equatorial
     * crossing and point 1.
     **********************************************************************/
    Math::real EquatorialArc() const throw() {
      return Init() ? atan2(_ssig1, _csig1) / Constants::degree() : 0;
    }

    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return Init() ? _a : 0; }

    /**
     * @return \e r the inverse flattening of the ellipsoid.  This is the
     *   value inherited from the Geodesic object used in the constructor.  A
     *   value of 0 is returned for a sphere (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return Init() ? _r : 0; }

    /**
     * @return \e caps the computational capabilities that this object was
     *   constructed with.  LATITUDE and AZIMUTH are always included.
     **********************************************************************/
    unsigned Capabilities() const throw() { return _caps; }

    /**
     * @param[in] testcaps a set of bitor'ed GeodesicLine::mask values.
     * @return true if the GeodesicLine object has all these capabilities.
     **********************************************************************/
    bool Capabilities(unsigned testcaps) const throw() {
      testcaps &= OUT_ALL;
      return (_caps & testcaps) == testcaps;
    }
    ///@}

    /** \name Deprecated Functions
     **********************************************************************/
    ///@{

    /**
     * <b>DEPRECATED</b>.  Return the latitude, \e lat2, longitude, \e lon2,
     * and forward azimuth, \e azi2 (degrees) of the point 2 which is a
     * distance, \e s12 (in meters), from point 1.  Also return the reduced
     * length \e m12 (meters).  \e s12 can be signed.  If \e arcmode (default
     * false) is set to true, \e s12 is interpreted as the arc length \e a12
     * (in degrees) on the auxiliary sphere.  Returned value is the arc length
     * \e a12 (degrees) if \e arcmode is false, otherwise it is the distance \e
     * s12 (meters).
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2, real &m12, bool arcmode)
      const throw() {
      if (arcmode) {
        real s12x;
        ArcPosition(s12, lat2, lon2, azi2, s12x, m12);
        return s12x;
      } else
        return Position(s12, lat2, lon2, azi2, m12);
    }

    /**
     * <b>DEPRECATED</b>.  Return the scale of the geodesic line extending an
     * arc length \e a12 (degrees) from point 1 to point 2.  \e M12 (a number)
     * measures the convergence of initially parallel geodesics.  It is defined
     * by the following construction: starting at point 1 proceed at azimuth \e
     * azi1 + 90<sup>o</sup> a small distance \e dt; turn -90<sup>o</sup> and
     * proceed a distance \e s12 (\e not the arc length \e a12); the distance
     * to point 2 is given by \e M12 \e dt.  \e M21 is defined analogously.
     **********************************************************************/
    void Scale(real a12, real& M12, real& M21) const throw() {
      real lat2, lon2, azi2, s12;
      ArcPosition(a12, lat2, lon2, azi2, s12, M12, M21);
    }
    ///@}

  };

} // namespace GeographicLib
#endif
