/**
 * \file Geodesic.hpp
 * \brief Header for GeographicLib::Geodesic class
 *
 * Copyright (c) Charles Karney (2009, 2010) <charles@karney.com>
 * and licensed under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEODESICLINE_HPP)
#define GEOGRAPHICLIB_GEODESICLINE_HPP "$Id$"

#include "GeographicLib/Constants.hpp"
#include "GeographicLib/Geodesic.hpp"

namespace GeographicLib {

  /**
   * \brief A geodesic line.
   *
   * GeodesicLine facilitates the determination of a series of points on a
   * single geodesic.  Geodesic.Line returns a GeodesicLine object with the
   * geodesic defined by by \e lat1, \e lon1, and \e azi1.
   * GeodesicLine.Position returns the \e lat2, \e lon2, \e azi2, and \e m12
   * given \e s12.  An example of use of this class is:
   \verbatim
   // Print positions on a geodesic going through latitude 30,
   // longitude 10 at azimuth 80.  Points at intervals of 10km
   // in the range [-1000km, 1000km] are given.
   GeographicLib::GeodesicLine
   line(GeographicLib::Geodesic::WGS84.Line(30.0, 10.0, 80.0));
   double step = 10e3;
   for (int s = -100; s <= 100; ++s) {
     double lat2, lon2, azi2, m12;
     double s12 = s * step;
     line.Position(s12, lat2, lon2, azi2, m12);
     std::cout << s12 << " " << lat2 << " " << lon2 << " "
          << azi2 << " " << m12 << "\n";
    }
   \endverbatim
   * The default copy constructor and assignment operators work with this
   * class, so that, for example, the previous example could start
   \verbatim
   GeographicLib::GeodesicLine line;
   line = GeographicLib::Geodesic::WGS84.Line(30.0, 10.0, 80.0);
   ...
   \endverbatim
   * Similarly, a vector can be used to hold GeodesicLine objects.
   *
   * This class is also used to give access to two other pieces of information
   * about a geodesic line:
   * - The geodesic scale accessible with GeodesicLine::Scale.
   * - The area under a geodesic accessible with GeodesicLine::Area.
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
      CAP_NONE = 0U,
      CAP_C1 = 1U<<0,
      CAP_C1p= 1U<<1,
      CAP_C2 = 1U<<2,
      CAP_C3 = 1U<<3,
      CAP_C4 = 1U<<4,
      CAP_ALL= 0x1FU,           // Bits 0 thru 4
      OUT_ALL= 0x7F80U,         // Bits 7 thru 14
    };
  public:

    /**
     * Bit masks for what calculations to do.
     * arc to position
     * distance output
     * distance input
     * reduced length
     * scale
     * area

     A1m1f   A2m1f A3f
     C1f C1pf  C2f  C3f C4f
     ARCLENGTH_IN
     DISTANCE_IN             X
     DISTANCE            X
     AZIMUTH
     REDUCEDLENGTH       X         X
     GEODESICSCALE       X         X
     AREA                                   X
     LONGITUDE                         X
    **********************************************************************/
    enum mask {
      NONE          = 0U,
      LATITUDE      = 1U<<7  | CAP_NONE,
      LONGITUDE     = 1U<<8  | CAP_C3,
      AZIMUTH       = 1U<<9  | CAP_NONE,
      DISTANCE      = 1U<<10 | CAP_C1,
      DISTANCE_IN   = 1U<<11 | CAP_C1 | CAP_C1p,
      REDUCEDLENGTH = 1U<<12 | CAP_C1 | CAP_C2,
      GEODESICSCALE = 1U<<13 | CAP_C1 | CAP_C2,
      AREA          = 1U<<14 | CAP_C4,
      ALL           = OUT_ALL| CAP_ALL,
    };

    /** \name Constructors
     **********************************************************************/
    ///@{

    /**
     * Constructor for a geodesic line staring at latitude \e lat1, longitude
     * \e lon1, and aziumuth \e azi1 (all in degrees).  If the point is at a
     * pole, the azimuth is defined by keeping the \e lon1 fixed and writing \e
     * lat1 = 90 - \e eps or -90 + \e eps and taking the limit \e eps -> 0 from
     * above.
     *
     * Return area below the geodesic in
     * meters<sup>2</sup>.  This must be called after GeodesicLine::AreaEnable.
     * The end point of the geodesic is specified as a spherical arc length, \e
     * a12 (in degrees), from the initial point.  The area is computed for the
     * geodesic quadilateral joining the points (0,\e lon1), (\e lat1,\e lon1),
     * (\e lat2,\e lon2), (0,\e lon2), and (0,\e lon1), with a clockwise
     * traversal counted as positive and with the segment from (0,\e lon2) to
     * (0,\e lon1) taken along the equator (even if this isn't the shortest
     * path).  If the geodesic is a meridian passing over a pole, that it is
     * deemed to pass very close to the pole in an easterly direction for the
     * purposes of determining the area.  The following function can be used to
     * compute the area of a geodesic polygon.  \verbatim #include <vector>
     * #include <utility>
     // Return the area of a geodesic polygon with vertices, pts, which is a
     // vector of latitude/longitude pairs.  Clockwise traversal of the polygon
     // counts as positive.  If the polygon encircles a pole one or more times
     // EllipsoidArea()/2 should be added (or subtracted) to the result for
     // each clockwise (or anticlockwise) encirclement.  This is left as an
     // exercise.
     double PolygonArea(const Geodesic& geod,
                        const std::vector<std::pair<double,double> >& pts) {
       int n = pts.size();
       double area = 0;
       for (int i = 0; i < n; ++i) {
         int j = (i + 1) % n;
         double azi1, azi2, s12, m12;
         double a12 = geod.Inverse(pts[i].first, pts[i].second,
                                   pts[j].first, pts[j].second,
                                   s12, azi1, azi2, m12);
         GeodesicLine l(geod.Line(pts[i].first, pts[i].second, azi1));
         l.AreaEnable(geod);
         area += l.Area(a12);
       }
       return area;
     }
     \endverbatim
     *
     * The accuracy of the result is about 0.25 m<sup>2</sup>.
     **********************************************************************/
    GeodesicLine(const Geodesic& g, real lat1, real lon1, real azi1,
                 unsigned caps = ALL)
      throw();

    /**
     * A default constructor.  If GeodesicLine::Position is called on the
     * resulting object, it returns immediately (without doing any
     * calculations).  The object should be set with a call to Geodesic::Line.
     * Use Init() to test whether object is still in this uninitialized state.
     **********************************************************************/
    GeodesicLine() throw() : _caps(0U) {};
    ///@}

    /** \name Position in terms of distance
     **********************************************************************/
    ///@{

    /**
     * The general position routine.
     **********************************************************************/
    Math::real Position(bool arcmode, real a12, unsigned outmask,
                        real& lat2, real& lon2, real& azi2,
                        real& s12, real& m12, real& M12, real& M21,
                        real& S12) const throw();

    /**
     * Return the latitude, \e lat2, and longitude, \e lon2 of the point 2
     * which is a distance, \e s12 (in meters), from point 1.  \e s12 can be
     * signed.  Returned function value is the arc length \e a12 (degrees).
     * The GeodesicLine object should have been constructed with \e caps =
     * CAP_DISTANCE_IN | CAP_LATITUDE | CAP_LONGITUDE.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  \e s12 can be signed.  Returned function value
     * is the arc length \e a12 (degrees).  The GeodesicLine object should have
     * been constructed with \e caps = CAP_DISTANCE_IN | CAP_LATITUDE |
     * CAP_LONGITUDE | CAP_AZIMUTH.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2, real& azi2)
      const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  \e s12 can be signed.  Returned function value
     * is the arc length \e a12 (degrees).  Also return the reduced length \e
     * m12 (meters).  The GeodesicLine object should have been constructed with
     * \e caps = CAP_DISTANCE_IN | CAP_LATITUDE | CAP_LONGITUDE | CAP_AZIMUTH |
     * CAP_REDUCEDLENGTH.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2, real& azi2, real& m12)
      const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  \e s12 can be signed.  Returned function value
     * is the arc length \e a12 (degrees).  Also return the geodesic scales \e
     * M12 and \e M21 (dimensionless).  The GeodesicLine object should have
     * been constructed with \e caps = CAP_DISTANCE_IN | CAP_LATITUDE |
     * CAP_LONGITUDE | CAP_AZIMUTH | CAP_GEODESICSCALE.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2,
                        real& azi2, real& M12, real& M21) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  \e s12 can be signed.  Returned function value
     * is the arc length \e a12 (degrees).  Also return the reduced length \e
     * m12 (meters) and the geodesic scales \e M12 and \e M21 (dimensionless)
     * The GeodesicLine object should have been constructed with \e caps =
     * CAP_DISTANCE_IN | CAP_LATITUDE | CAP_LONGITUDE | CAP_AZIMUTH |
     * CAP_REDUCEDLENGTH | CAP_GEODESICSCALE.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2, real& azi2,
                        real& m12, real& M12, real& M21) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward azimuth,
     * \e azi2 (degrees) of the point 2 which is a distance, \e s12 (in
     * meters), from point 1.  \e s12 can be signed.  Returned function value
     * is the arc length \e a12 (degrees).  Also return the reduced length \e
     * m12 (meters), the geodesic scales \e M12 and \e M21 (dimensionless), and
     * the area \e S12 (meter<sup>2</sup>).  The GeodesicLine object should
     * have been constructed with \e caps = CAP_DISTANCE_IN | CAP_LATITUDE |
     * CAP_LONGITUDE | CAP_AZIMUTH | CAP_REDUCEDLENGTH | CAP_GEODESICSCALE |
     * CAP_AREA.
     **********************************************************************/
    Math::real Position(real s12, real& lat2, real& lon2, real& azi2,
                        real& m12, real& M12, real& M21, real& S12)
      const throw();

    ///@}

    /** \name Position in terms of arc length
     **********************************************************************/
    ///@{

    /**
     * Return the latitude, \e lat2, and longitude, \e lon2 of the point 2
     * which is an arc length, \e a12 (in degree), from point 1.  \e a12, which
     * can be signed, is the arc length on the auxiliary sphere.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2)
      const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2)
      const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.  Also return the distance \e s12 (in
     * meters).
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2, real& s12)
      const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.  Also return the distance \e s12 (in
     * meters) and the reduced length \e m12 (meters).
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.  Also return the distance \e s12 (in
     * meters) and the geodesic scales \e M12 and \e M21 (dimensionless).
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& M12, real& M21) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.  Also return the distance \e s12 (in
     * meters), the reduced length \e m12 (meters), and the geodesic scales \e
     * M12 and \e M21 (dimensionless).
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12, real& M12, real& M21) const throw();

    /**
     * Return the latitude, \e lat2, longitude, \e lon2, and forward
     * azimuth, \e azi2 (degrees) of the point 2 which is an arc length, \e a12
     * (in degree), from point 1.  \e a12, which can be signed, is the arc
     * length on the auxiliary sphere.  Also return the distance \e s12 (in
     * meters), the reduced length \e m12 (meters), the geodesic scales \e M12
     * and \e M21 (dimensionless), and the area \e S12 (meter<sup>2</sup>).
     **********************************************************************/
    void ArcPosition(real a12, real& lat2, real& lon2, real& azi2,
                     real& s12, real& m12, real& M12, real& M21, real& S12)
      const throw();

    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{

    /**
     * Has this object been initialized so that Position can be called?
     **********************************************************************/
    bool Init() const throw() { return _caps != 0U; }


    /**
     * Return total area of ellipsoid in meters<sup>2</sup>.  (Does not require
     * GeodesicLine::AreaEnable to have been called.)  The area of a polygon
     * encircling a pole can be found by adding GeodesicLine::EllipsoidArea()/2
     * to the sum of GeodesicLine::Area for each side of the polygon.
     **********************************************************************/
    Math::real EllipsoidArea() const throw() {
      return 4 * Constants::pi() * _c2;
    }

    /**
     * Return the latitude of point 1 (degrees).
     **********************************************************************/
    Math::real Latitude() const throw() { return Init() ? _lat1 : 0; }

    /**
     * Return the longitude of point 1 (degrees).
     **********************************************************************/
    Math::real Longitude() const throw() { return Init() ? _lon1 : 0; }

    /**
     * Return the azimuth (degrees) of the geodesic line as it passes through
     * point 1.
     **********************************************************************/
    Math::real Azimuth() const throw() { return Init() ? _azi1 : 0; }

    /**
     * Return the azimuth (degrees) of the geodesic line as it crosses the
     * equator in a northward direction.
     **********************************************************************/
    Math::real EquatorialAzimuth() const throw() {
      return Init() ? atan2(_salp0, _calp0) / Constants::degree() : 0;
    }

    /**
     * Return the arc length (degrees) between the northward equatorial
     * crossing and point 1.
     **********************************************************************/
    Math::real EquatorialArc() const throw() {
      return Init() ? atan2(_ssig1, _csig1) / Constants::degree() : 0;
    }

    /**
     * The major radius of the ellipsoid (meters).  This is that value of \e a
     * inherited from the Geodesic object used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return Init() ? _a : 0; }

    /**
     * The inverse flattening of the ellipsoid.  This is that value of \e r
     * inherited from the Geodesic object used in the constructor.  A value of
     * 0 is returned for a sphere (infinite inverse flattening).
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return Init() ? _r : 0; }

    /**
     * The computational capabilities that this object was constructed with.
     * LATITUDE and AZIMUTH are always included.
     **********************************************************************/
    unsigned Capabilities() const throw() { return _caps; }
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
