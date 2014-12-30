/**
 * \file JacobiConformal.hpp
 * \brief A class for Jacobi's conformal projection of a triaxial ellipsoid.
 *
 * Copyright (c) Charles Karney (2014) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/EllipticFunction.hpp>

namespace GeographicLib {
  /**
   * \brief Jacobi's conformal projection
   *
   * This is a conformal projection of the ellipsoid to a plane in which the
   * grid lines are straight; see Jacobi, Vorlesungen ueber Dynamik, Sect. 28.
   * The constructor takes the semi-axes of the ellipsoid (which must be in
   * order).  Member functions map the ellipsoidal coordinates &omega; and
   * &beta; separately to \e x and \e y.  Jacobi's coordinates have been
   * multiplied by
   * (<i>a</i><sup>2</sub>&minus;<i>c</i><sup>2</sub>)<sup>1/2</sup> /
   * (2<i>b</i>) so that the customary results are returned in the cases of a
   * sphere or an ellipsoid of revolution.
   **********************************************************************/
  class JacobiConformal {
    Math::real _a, _b, _c, _ab2, _bc2, _ac2;
    EllipticFunction _ex, _ey;
    static void norm(Math::real& x, Math::real& y)
    { Math::real z = Math::hypot(x, y); x /= z; y /= z; }
  public:
    /**
     * Constructor for a trixial ellipsoid with semi-axes
     *
     * @param[in] a
     * @param[in] b
     * @param[in] c
     *
     * The semi-axes must satisfy \e a &ge; \e b &ge; \e c > 0 and \e a > \e c.
     * This form of the constructor cannot be used to specify a sphere (use the
     * next constructor).
     **********************************************************************/
    JacobiConformal(Math::real a, Math::real b, Math::real c)
      : _a(a), _b(b), _c(c)
      , _ab2((_a - _b) * (_a + _b))
      , _bc2((_b - _c) * (_b + _c))
      , _ac2((_a - _c) * (_a + _c))
      , _ex(_ab2 / _ac2 * Math::sq(_c / _b), -_ab2 / Math::sq(_b),
            _bc2 / _ac2 * Math::sq(_a / _b), Math::sq(_a / _b))
      , _ey(_bc2 / _ac2 * Math::sq(_a / _b), +_bc2 / Math::sq(_b),
            _ab2 / _ac2 * Math::sq(_c / _b), Math::sq(_c / _b))
    {
      using std::sqrt;
      if (!(a >= b && b >= c && c > 0 && a > c))
        throw GeographicErr("axes are not in order");
    }
    /**
     * Constructor for a triaxial ellipsoid with semi-axes
     *
     * @param[in] a
     * @param[in] b
     * @param[in] c
     * @param[in] ab the relative magnitude of \e a &minus \e b.
     * @param[in] bc the relative magnitude of \e b &minus \e c.
     *
     * This form can be used to specify a sphere.  The semi-axes must satisfy
     * \e a &ge \e b &ge c > 0.  The ratio \e ab : \e bc must equal
     * (<i>a</i>&minus;<i>b</i>) : (<i>b</i>&minus;<i>c</i>) with \e ab &ge; 0,
     * \e bc &be; 0, and \e ab + \e bc > 0.
     **********************************************************************/
    JacobiConformal(Math::real a, Math::real b, Math::real c,
                    Math::real ab, Math::real bc)
      : _a(a), _b(b), _c(c)
      , _ab2(ab * (_a + _b))
      , _bc2(bc * (_b + _c))
      , _ac2(_ab2 + _bc2)
      , _ex(_ab2 / _ac2 * Math::sq(_c / _b),
            -(_a - _b) * (_a + _b) / Math::sq(_b),
            _bc2 / _ac2 * Math::sq(_a / _b), Math::sq(_a / _b))
      , _ey(_bc2 / _ac2 * Math::sq(_a / _b),
            +(_b - _c) * (_b + _c) / Math::sq(_b),
            _ab2 / _ac2 * Math::sq(_c / _b), Math::sq(_c / _b))
    {
      using std::sqrt;
      if (!(a >= b && b >= c && c > 0 && ab >= 0 &&  bc >= 0 && ab + ab > 0))
        throw GeographicErr("axes are not in order");
    }
    /**
     * @return the quadrant length in the \e x direction
     **********************************************************************/
    Math::real x() const { return Math::sq(_a / _b) * _ex.Pi(); }
    /**
     * The \e x projection
     *
     * @param[in] somg sin(&omega;)
     * @param[in] comg cos(&omega;)
     * @return \e x
     **********************************************************************/
    Math::real x(Math::real somg, Math::real comg) const {
      Math::real somg1 = _b * somg, comg1 = _a * comg; norm(somg1, comg1);
      return Math::sq(_a / _b) * _ex.Pi(somg1, comg1, _ex.Delta(somg1, comg1));
    }
    /**
     * The \e x projection
     *
     * @param[in] omg &omega; (in degrees)
     * @return \e x (in degrees)
     *
     * &omega; must be in (&minus;180&deg;, 180&deg;].
     **********************************************************************/
    Math::real x(Math::real omg) const {
      using std::abs; using std::sin; using std::cos;
      Math::real
        a = omg * Math::degree(),
        somg = abs(omg) == 180 ? 0 : sin(a),
        comg = abs(omg) ==  90 ? 0 : cos(a);
      return x(somg, comg) / Math::degree();
    }
    /**
     * @return the quadrant length in the \e y direction
     **********************************************************************/
    Math::real y() const { return Math::sq(_c / _b) * _ey.Pi(); }
    /**
     * The \e y projection
     *
     * @param[in] sbet sin(&beta;)
     * @param[in] cbet cos(&beta;)
     * @return \e y
     **********************************************************************/
    Math::real y(Math::real sbet, Math::real cbet) const {
      Math::real sbet1 = _b * sbet, cbet1 = _c * cbet; norm(sbet1, cbet1);
      return Math::sq(_c / _b) * _ey.Pi(sbet1, cbet1, _ey.Delta(sbet1, cbet1));
    }
    /**
     * The \e y projection
     *
     * @param[in] bet &beta; (in degrees)
     * @return \e y (in degrees)
     *
     * &beta; must be in (&minus;180&deg;, 180&deg;].
     **********************************************************************/
    Math::real y(Math::real bet) const {
      using std::abs; using std::sin; using std::cos;
      Math::real
        a = bet * Math::degree(),
        sbet = abs(bet) == 180 ? 0 : sin(a),
        cbet = abs(bet) ==  90 ? 0 : cos(a);
      return y(sbet, cbet) / Math::degree();
    }
  };

} // namespace GeographicLib
