/**
 * \file Cartesian3.hpp
 * \brief Header for GeographicLib::Triaxial::Cartesian3 class
 *
 * Copyright (c) Charles Karney (2025) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CARTESIAN3_HPP)
#define GEOGRAPHICLIB_CARTESIAN3_HPP 1

#include <utility>
#include <functional>
#include <random>
#include <GeographicLib/Triaxial/Ellipsoid3.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs random
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {
  namespace Triaxial {

  /**
   * \brief Transformations between cartesian and triaxial coordinates
   *
   * The Cartesian3 class supports transformations between cartesian
   * coordinates and various coordinates for a triaxial ellipsoid.  This is
   * covered in Appendices A and B of
   * - C. F. F. Karney,<br>
   *   <a href="https://arxiv.org/abs/2511.01621">
   *   Jacobi's solution for geodesics on a triaxial ellipsoid</a>,<br>
   *   Technical Report, SRI International, Nov. 2025.<br>
   *   <a href="https://arxiv.org/abs/2511.01621">arxiv:2511.01621</a>
   *
   * Besides ellipsoidal coordinates defined in Ellipsoid3, the following
   * coordinates are supported:
   * * geodetic coordinates \f$(\phi, \lambda)\f$ defined by
   *   \f[
   *   \hat{\mathbf U} =
   *   [\cos\phi \cos\lambda, \cos\phi \sin\lambda, \sin\phi]^T,
   *   \f]
   *   where \f$\hat{\mathbf U}\f$ is the normal to the surface of the
   *   ellipsoid.
   * * parametric coordinates \f$(\phi', \lambda')\f$ defined by
   *   \f[
   *   \mathbf R =
   *   [a \cos\phi' \cos\lambda', b \cos\phi' \sin\lambda',
   *   c \sin\phi']^T,
   *   \f]
   * * geocentric coordinates \f$(\phi'', \lambda'')\f$ defined by
   *   \f[
   *   \hat{\mathbf R} =
   *   [\cos\phi'' \cos\lambda'', \cos\phi'' \sin\lambda'', \sin\phi'']^T.
   *   \f]
   * .
   * For each of these 3 coordinates, the "north pole" is at \f$[0, 0, c]^T\f$
   * and the origin for longitudes is \f$[a, 0, 0]^T\f$.  We also define
   * alternate versions (named "geodetic*", etc., where the north pole is
   * placed at \f$[a, 0, 0]^T\f$ and the origin for longitude is \f$[0, 0,
   * -c]\f$.  This latter set of coordinates is appropriate for ellipsoids that
   * are nearly prolate.
   *
   * Directions on the ellipsoid are easily specified in cartesian coordinates
   * as a vector tangent to the surface of the ellipsoid.  This is converted to
   * a heading by defined the angle the vector makes (measured clockwise) from
   * the coordinate-specific north.  This is defined as the direction of a line
   * of constant (coordinate-specific) longitude.  The resulting heading is
   * denoted by \f$\alpha\f$ for ellipsoidal coordinates and by \f$\zeta\f$ for
   * the other coordinates.  The unstarred coordinates all share the same
   * direction for north, and likewise for the starred coordinates.  Note that
   * the lines of constant longitude and latitude are only orthogonal (in
   * general) for ellipsoidal coordinates.
   *
   * Arbitrary points (not necessarily lying on the ellipsoid) an additional
   * "height" is required to specify the position.  For ellipsoidal
   * coordinates, we find the confocal ellipsoid on which the point lies and
   * the height is then defined as \f$H = u - c\f$ where \f$u\f$ is the
   * semiminor axes of the confocal ellipsoid; the ellipsoid latitude and
   * longitude are those for the confocal ellipsoid For the other coordinates
   * systems, we define \f$h\f$ a the height above the closest point on the
   * ellipsoid and the latitude and longitude refer to the closest point.
   *
   * \note The family of confocal ellipsoids has semiaxes \f$[\sqrt{a^2 - c^2 +
   *   u^2}, \sqrt{b^2 - c^2 + u^2}, u]\f$.
   *
   * \note In the function names "any" stands for any of the seven coordinate
   *   systems enumerated by Cartesian3::coord.  "cart2" refers to a point
   *   given in cartesian coordinates that lies on the ellipsoid.  On the other
   *   hand, "cart" refers to an arbitrary point.
   *
   * Example of use:
   * \include example-Cartesian3.cpp
   *
   * <a href="Cart3Convert.1.html">Cart3Convert</a> is a command-line utility
   * providing access to the functionality of Cartestian3.
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT Cartesian3 {
  public:
    /**
     * A type to hold three-dimensional positions and directions in cartesian
     * coordinates.
     **********************************************************************/
    using vec3 = Ellipsoid3::vec3;
  private:
    using real = Math::real;
#if GEOGRAPHICLIB_PRECISION > 3
    // <random> only supports "standard" floating point types
    using random_prec = Math::extended;
#else
    using random_prec = Math::real;
#endif
    using ang = Angle;
    static constexpr int maxit_ = 20;
    static constexpr bool throw_ = true; // exception on convergence failure
    const Ellipsoid3 _t;
    const vec3 _axes, _axes2, _linecc2;
    // mutable because using these objects in a non-const operation
    mutable std::normal_distribution<random_prec> _norm;
    mutable std::uniform_real_distribution<random_prec> _uni;

    static void roty(vec3& R, int n) {
      // require n = -1, 0, 1
      // Prolate convention has major axis in z direction, minor axis in -x
      // direction, median axis is unchanged.
      // If n = 0, do nothing otherwise...
      // With n = +1, multiply by
      //  [ 0  0 -1]
      //  [ 0  1  0]
      //  [ 1  0  0]
      // which transforms original x, y, z to prolate convention.
      // With n = -1, multiply by
      //  [ 0  0  1]
      //  [ 0  1  0]
      //  [-1  0  0]
      // which transforms prolate convention to  original x, y, z.
      if (n != 0) {
        using std::swap;
        R[1+n] = -R[1+n];
        swap(R[0], R[2]);
      }
    }

    template<int n>
    void cart2togeneric(vec3 R, ang& phi, ang& lam, bool alt) const;
    template<int n>
    void generictocart2(ang phi, ang lam, vec3& R, bool alt) const;
    template<int n> ang meridianplane(ang lam, bool alt) const;
    void cardinaldir(vec3 R, ang merid, vec3& N, vec3& E, bool alt) const;
    template<int n>
    void cart2togeneric(vec3 R, vec3 V, ang& phi, ang& lam, ang& zet, bool alt)
      const;
    template<int n>
    void generictocart2(ang phi, ang lam, ang zet, vec3& R, vec3&V, bool alt)
      const;
    real cubic(vec3 R2) const;

    template<int n>
    class funp {
    private:
      // Evaluate
      //   f(p) = sum( (R[0]/(p + l[0]))^n, k = 0..2) - 1
      // and it derivative.
      const real _d;
      const vec3 _r, _l;
    public:
      funp(const vec3& R, const vec3& l)
        : _d(std::numeric_limits<real>::epsilon()/2)
        , _r(R)
        , _l(l)
      {
        static_assert(n >= 1 && n <= 2, "Bad power in funp");
      }
      std::pair<real, real> operator()(real p) const;
    };

    static real cartsolve(const std::function<std::pair<real, real>(real)>& f,
                          real p0, real pscale);
    void carttoellip(vec3 R, Angle& bet, Angle& omg, real& H) const;
    void elliptocart(Angle bet, Angle omg, real H, vec3& R) const;

    // real a() const { return t().a(); } // not needed
    real b() const { return t().b(); }
    real c() const { return t().c(); }
  public:
    /**
     * Enumerator for all the coordinates.
     **********************************************************************/
    enum coord {
      /**
       * Geodetic coordinates, \e phi, \e lam, \e zet \e h;
       * @hideinitializer
       **********************************************************************/
      GEODETIC = 0,
      /**
       * Parametric coordinates, \e phi', \e lam', \e zet, \e h;
       * @hideinitializer
       **********************************************************************/
      PARAMETRIC = 1,
      /**
       * %Geocentric coordinates, \e phi'', \e lam'', \e zet, \e h;
       * @hideinitializer
       **********************************************************************/
      GEOCENTRIC = 2,
      /**
       * Ellipsoidal coordinates, \e beta, \e omg, \e alp, \e H;
       * @hideinitializer
       **********************************************************************/
      ELLIPSOIDAL = 3,
      /**
       * Geodetic coordinates with pole aligned with the major axis.
       * @hideinitializer
       **********************************************************************/
      GEODETIC_X = 4 + GEODETIC,
      /**
       * Parametric coordinates with pole aligned with the major axis.
       * @hideinitializer
       **********************************************************************/
      PARAMETRIC_X = 4 + PARAMETRIC,
      /**
       * %Geocentric coordinates with pole aligned with the major axis.
       * @hideinitializer
       **********************************************************************/
      GEOCENTRIC_X = 4 + GEOCENTRIC,
      /**
       * An alias for GEODETIC;
       * @hideinitializer
       **********************************************************************/
      PLANETODETIC = GEODETIC,
      /**
       * Another alias for GEODETIC;
       * @hideinitializer
       **********************************************************************/
      GEOGRAPHIC = GEODETIC,
      /**
       * An alias for GEOCENTRIC;
       * @hideinitializer
       **********************************************************************/
      PLANETOCENTRIC = GEOCENTRIC,
    };
    /** \name Transformations for points on the ellipsoid.
     **********************************************************************/
    ///@{
    /**
     * Constructor for a triaxial ellipsoid defined by Ellipsoid3 object.
     *
     * @param[in] t the Ellipsoid3 object.
     **********************************************************************/
    Cartesian3(const Ellipsoid3& t);
    /**
     * Constructor for a triaxial ellipsoid with semiaxes.
     *
     * @param[in] a the largest semiaxis.
     * @param[in] b the middle semiaxis.
     * @param[in] c the smallest semiaxis.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * The semiaxes must satisfy \e a &ge; \e b &ge; \e c &gt; 0.
     * If \e a = \e c (a sphere), then the oblate limit is taken.
     **********************************************************************/
    Cartesian3(real a, real b, real c);
    /**
     * Alternate constructor for a triaxial ellipsoid.
     *
     * @param[in] b the middle semiaxis.
     * @param[in] e2 the eccentricity squared \f$e^2 = (a^2 - c^2)/b^2\f$.
     * @param[in] k2 the oblateness parameter squared \f$k^2 = (b^2 - c^2) /
     *  (a^2 - c^2)\f$.
     * @param[in] kp2 the prolateness parameter squared \f$k'^2= (a^2 - b^2) /
     *   (a^2 - c^2)\f$.
     * @exception GeographicErr if the required ordering is semiaxes is
     *   violated.
     *
     * \note The constructor normalizes \e k2 and \e kp2 to ensure then \e k2 +
     * \e kp2 = 1.
     **********************************************************************/
    Cartesian3(real b, real e2, real k2, real kp2);
    ///@}

    /** \name Transformations for points on the ellipsoid.
     **********************************************************************/
    ///@{
    /**
     * Convert latitude and longitude to a point on the surface.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point.
     * @param[in] lon the longitude of the point.
     * @param[out] R the cartesian position on the surface of the ellipsoid.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart2(coord coordin, Angle lat, Angle lon, vec3& R) const;
    /**
     * Convert latitude and longitude in degrees to a point on the surface.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point (in degrees).
     * @param[in] lon the longitude of the point (in degrees).
     * @param[out] R the cartesian position on the surface of the ellipsoid.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart2(coord coordin, real lat, real lon, vec3& R) const {
      anytocart2(coordin, Angle(lat), Angle(lon), R);
    }
    /**
     * Convert a point on the surface to latitude and longitude.
     *
     * @param[in] R the cartesian position on the surface of the ellipsoid.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point.
     * @param[out] lon the longitude of the point.
     * @exception GeographicErr if \e coordout is not recognized.
     **********************************************************************/
    void cart2toany(vec3 R, coord coordout, Angle& lat, Angle& lon) const;
    /**
     * Convert a point on the surface to latitude and longitude in degrees.
     *
     * @param[in] R the cartesian position on the surface of the ellipsoid.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point (in degrees).
     * @param[out] lon the longitude of the point (in degrees).
     * @exception GeographicErr if \e coordout is not recognized.
     **********************************************************************/
    void cart2toany(vec3 R, coord coordout, real& lat, real& lon) const {
      Angle lata, lona; cart2toany(R, coordout, lata, lona);
      lat = real(lata); lon = real(lona);
    }
    /**
     * Convert between latitudes and longitudes.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat1 the \e coordin latitude of the point.
     * @param[in] lon1 the \e coordin longitude of the point.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat2 the \e coordout latitude of the point.
     * @param[out] lon2 the \e coordout longitude of the point.
     * @exception GeographicErr if \e coordin or \e coordout is not recognized.
     **********************************************************************/
    void anytoany(coord coordin, Angle lat1, Angle lon1,
                  coord coordout, Angle& lat2, Angle& lon2) const;
    /**
     * Convert between latitudes and longitudes in degrees.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat1 the \e coordin latitude of the point (in degrees).
     * @param[in] lon1 the \e coordin longitude of the point (in degrees).
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat2 the \e coordout latitude of the point (in degrees).
     * @param[out] lon2 the \e coordout longitude of the point (in degrees).
     * @exception GeographicErr if \e coordin or \e coordout is not recognized.
     **********************************************************************/
    void anytoany(coord coordin, real lat1, real lon1,
                  coord coordout, real& lat2, real& lon2) const {
      Angle lat2a, lon2a;
      anytoany(coordin, Angle(lat1), Angle(lon1), coordout, lat2a, lon2a);
      lat2 = real(lat2a); lon2 = real(lon2a);
    }
    ///@}

    /** \name Transformations for points and directions on the ellipsoid.
     **********************************************************************/
    ///@{
    /**
     * Convert latitude, longitude, and azimuth to cartesian position and
     * direction.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point.
     * @param[in] lon the longitude of the point.
     * @param[in] azi the azimuth of the heading.
     * @param[out] R the cartesian position on the surface of the ellipsoid.
     * @param[out] V the cartesian direction tangent to the ellipsoid.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart2(coord coordin, Angle lat, Angle lon, Angle azi,
                    vec3& R, vec3& V) const;
    /**
     * Convert latitude, longitude, and azimuth in degrees to cartesian
     * position and direction.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point (in degrees).
     * @param[in] lon the longitude of the point (in degrees).
     * @param[in] azi the azimuth of the heading (in degrees).
     * @param[out] R the cartesian position on the surface of the ellipsoid.
     * @param[out] V the cartesian direction tangent to the ellipsoid.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart2(coord coordin, real lat, real lon, real azi,
                    vec3& R, vec3& V) const {
      anytocart2(coordin, Angle(lat), Angle(lon), Angle(azi), R, V);
    }
    /**
     * Convert position and direction on surface to latitude, longitude, and
     * azimuth.
     *
     * @param[in] R the cartesian position on the surface of the ellipsoid.
     * @param[in] V the cartesian direction tangent to the ellipsoid.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point.
     * @param[out] lon the longitude of the point.
     * @param[out] azi the azimuth of the heading.
     * @exception GeographicErr if \e coordout is not recognized.
     **********************************************************************/
    void cart2toany(vec3 R, vec3 V,
                    coord coordout, Angle& lat, Angle& lon, Angle& azi) const;
    /**
     * Convert position and direction on surface to latitude, longitude, and
     * azimuth in degrees.
     *
     * @param[in] R the cartesian position on the surface of the ellipsoid.
     * @param[in] V the cartesian direction tangent to the ellipsoid.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point (in degrees).
     * @param[out] lon the longitude of the point (in degrees).
     * @param[out] azi the azimuth of the heading (in degrees).
     * @exception GeographicErr if \e coordout is not recognized.
     **********************************************************************/
    void cart2toany(vec3 R, vec3 V,
                    coord coordout, real& lat, real& lon, real& azi) const {
      Angle lata, lona, azia; cart2toany(R, V, coordout, lata, lona, azia);
      lat = real(lata); lon = real(lona), azi = real(azia);
    }
    ///@}

    /** \name Transformations for arbitrary points.
     **********************************************************************/
    ///@{
    /**
     * Convert latitude, longitude, and height to a cartesian position.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point.
     * @param[in] lon the longitude of the point.
     * @param[in] h the height (in meters).
     * @param[out] R the cartesian position of the point.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart(coord coordin, Angle lat, Angle lon, real h, vec3& R) const;
    /**
     * Convert latitude, longitude in degrees, and height to a cartesian
     * position.
     *
     * @param[in] coordin one of the coordinate types, Cartesian3::coord.
     * @param[in] lat the latitude of the point (in degrees).
     * @param[in] lon the longitude of the point (in degrees).
     * @param[in] h the height (in meters).
     * @param[out] R the cartesian position of the point.
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void anytocart(coord coordin, real lat, real lon, real h, vec3& R) const {
      anytocart(coordin, Angle(lat), Angle(lon), h, R);
    }
    /**
     * Convert a cartesian position to latitude, longitude, and height.
     *
     * @param[in] R the cartesian position of the point.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point.
     * @param[out] lon the longitude of the point.
     * @param[out] h the height (in meters).
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void carttoany(vec3 R,
                   coord coordout, Angle& lat, Angle& lon, real& h) const;
    /**
     * Convert a cartesian position to latitude, longitude in degrees, and
     * height.
     *
     * @param[in] R the cartesian position of the point.
     * @param[in] coordout one of the coordinate types, Cartesian3::coord.
     * @param[out] lat the latitude of the point (in degrees).
     * @param[out] lon the longitude of the point (in degrees).
     * @param[out] h the height (in meters).
     * @exception GeographicErr if \e coordin is not recognized.
     **********************************************************************/
    void carttoany(vec3 R,
                   coord coordout, real& lat, real& lon, real& h) const {
      Angle lata, lona; carttoany(R, coordout, lata, lona, h);
      lat = real(lata); lon = real(lona);
    }
    ///@}

    /** \name Transferring an arbitrary point onto the ellipsoid.
     **********************************************************************/
    ///@{
    /**
     * Convert a point on the ellipsoid and a height to a cartesian position.
     *
     * @param[in] R2 the cartesian position of the point on the ellipsoid.
     * @param[in] h the height above the ellipsoid (in meters).
     * @param[out] R the cartesian position of the point.
     **********************************************************************/
    void cart2tocart(vec3 R2, real h, vec3& R) const;
    /**
     * Find the closest point on the ellipsoid
     *
     * @param[in] R the cartesian position of the point.
     * @param[out] R2 the cartesian position of the closest point on the
     *   ellipsoid.
     * @param[out] h the height above the ellipsoid (in meters).
     **********************************************************************/
    void carttocart2(vec3 R, vec3& R2, real& h) const;
    ///@}

    /** \name Generating random points on the ellipsoid.
     **********************************************************************/
    ///@{
    /**
     * Generate a random point on the ellipsoid.
     *
     * @tparam G the type of the random generator.
     * @param[inout] g the random generator.
     * @param[out] R a cartesian position uniformly sampled on the surface of
     *   the ellipsoid.
     *
     * See the example listed in the description of this class for an example
     * of using this function.
     *
     * The method of sampling is given by
     * <a href="https://doi.org/10.1007/s11075-023-01628-4"> Marples and
     * Williams (2023)</a> Algorithm 1, based on the general method of
     * <a href="https://doi.org/10.1088/0031-9155/32/10/009"> Williamson
     * (1987)</a>.
     **********************************************************************/
    template <class G> void cart2rand(G& g, vec3& R) const;
    /**
     * Generate a random point and direction on the ellipsoid.
     *
     * @tparam G the type of the random generator.
     * @param[inout] g the random generator.
     * @param[out] R a cartesian position uniformly sampled on the surface of
     *   the ellipsoid.
     * @param[out] V a cartesian direction uniformly sampled tangent to the
     *   ellipsoid.
     **********************************************************************/
    template <class G> void cart2rand(G& g, vec3& R, vec3& V) const;
    ///@}

    /** \name Inspector function
     **********************************************************************/
    ///@{
    /**
     * @return the Ellipsoid3 object for this projection.
     **********************************************************************/
    const Ellipsoid3& t() const { return _t; }
    ///@}
  };

  template<class G> inline void Cartesian3::cart2rand(G& g, vec3& R) const {
    // This uses the simple rejection technique given by Marples and Williams,
    // Num. Alg. (2023), Algorithm 1 based on the general method of Williamson,
    // Phys. Med. Biol. (1987).
    using std::isfinite;
    while (true) {
      while (true) {
        // guaranteed evaluated left to right
        R = {real(_norm(g)), real(_norm(g)), real(_norm(g))};
        Ellipsoid3::normvec(R); // But catch rare cases where |R| = 0
        if (isfinite(R[0])) break;
      }
      R[0] *= _axes[0]; R[1] *= _axes[1]; R[2] *= _axes[2];
      vec3 up{ R[0] / _axes2[0],  R[1] / _axes2[1],  R[2] / _axes2[2] };
      real q = c() * Math::hypot3(up[0], up[1], up[2]);
      if (real(_uni(g)) < q) break;
    }
  }
  template<class G> inline void Cartesian3::cart2rand(G& g, vec3& R, vec3& V)
  const {
    using std::isfinite;
    cart2rand<G>(g, R);
    while (true) {
      // guaranteed evaluated left to right
      V = {real(_norm(g)), real(_norm(g)), real(_norm(g))};
      vec3 up{ R[0] / _axes2[0],  R[1] / _axes2[1],  R[2] / _axes2[2] };
      real u2 = Math::sq(up[0]) + Math::sq(up[1]) + Math::sq(up[2]), // |up|^2
        // (up . V) / |up|^2
        uv = (V[0] * up[0] + V[1] * up[1] + V[2] * up[2])/u2;
      // V - up * (up . V) / |up|^2
      V[0] -= uv * up[0]; V[1] -= uv * up[1]; V[2] -= uv * up[2];
      Ellipsoid3::normvec(V);   // But catch rare cases where |V| = 0
      if (isfinite(V[0])) break;
    }
  }

  } // namespace Triaxial
} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_CARTESIAN3_HPP
