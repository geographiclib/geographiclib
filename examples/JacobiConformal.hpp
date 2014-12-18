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
   * The constructor takes the semi-axes of the ellipsoid (which must be
   * scalene).  Member functions map the ellipsoidal coordinates &omega; and
   * &beta; separately to \e x and \e y.
   **********************************************************************/
  class JacobiConformal {
    Math::real _a, _b, _c, _ab2, _bc2, _ac2, _ab, _bc, _ac;
    EllipticFunction _exa, _eya;
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
     * The semi-axes must satisfy \e a > \e b > \e c > 0.
     **********************************************************************/
    JacobiConformal(Math::real a, Math::real b, Math::real c)
      : _a(a), _b(b), _c(c)
      , _ab2((_a - _b) * (_a + _b))
      , _bc2((_b - _c) * (_b + _c))
      , _ac2((_a - _c) * (_a + _c))
      , _exa(+_ab2 / _ac2 * Math::sq(_c / _b), -_ab2 / Math::sq(_b),
             +_bc2 / _ac2 * Math::sq(_a / _b), +Math::sq(_a / _b) )
      , _eya(-_bc2 / _ab2 * Math::sq(_a / _c), -_bc2 / Math::sq(_c),
             +_ac2 / _ab2 * Math::sq(_b / _c), +Math::sq(_b / _c))
    {
      using std::sqrt;
      if (!(a > b && b > c && c > 0))
        throw GeographicErr("axes are not in order");
      _ab = sqrt(_ab2);
      _bc = sqrt(_bc2);
      _ac = sqrt(_ac2);
    }
    /**
     * @return the quadrant length in the \e x direction
     **********************************************************************/
    Math::real x() const {
      return 2 * Math::sq(_a) / (_b * _ac) * _exa.Pi();
    }
    /**
     * The \e x projection
     *
     * @param[in] somg sin(&omega;)
     * @param[in] comg cos(&omega;)
     * @return \e x
     **********************************************************************/
    Math::real x(Math::real somg, Math::real comg) const {
      using std::sqrt; using std::atan;
      Math::real somg1 = _ac * somg, comg1 = _bc * comg;
      norm(somg1, comg1);
      Math::real domg1 = _exa.Delta(somg1, comg1);
      return 2 *
        ( Math::sq(_a) / (_b * _ac) * _exa.Pi(somg1, comg1, domg1) -
          atan(_ab2 * somg1 * comg1 / sqrt(Math::sq(_a * somg1) * _bc2 +
                                           Math::sq(_b * comg1) * _ac2)) );
    }
    /**
     * The \e x projection
     *
     * @param[in] omg &omega; (in degrees)
     * @return \e x
     *
     * &omega; must be in (&minus;180&deg;, 180&deg;].
     **********************************************************************/
    Math::real x(Math::real omg) const {
      using std::abs; using std::sin; using std::cos;
      Math::real
        a = omg * Math::degree(),
        somg = abs(omg) == 180 ? 0 : sin(a),
        comg = abs(omg) ==  90 ? 0 : cos(a);
      return x(somg, comg);
    }
    /**
     * @return the quadrant length in the \e y direction
     **********************************************************************/
    Math::real y() const {
      return  2 * Math::sq(_b) / (_c * _ab) * _eya.Pi();
    }
    /**
     * The \e y projection
     *
     * @param[in] sbet sin(&beta;)
     * @param[in] cbet cos(&beta;)
     * @return \e y
     **********************************************************************/
    Math::real y(Math::real sbet, Math::real cbet) const {
      using std::sqrt;
      Math::real sbet1 = _ab * sbet, cbet1 = _ac * cbet;
      norm(sbet1, cbet1);
      Math::real dbet1 = _eya.Delta(sbet1, cbet1);
      return 2 *
        ( Math::sq(_b) / (_c * _ab) * _eya.Pi(sbet1, cbet1, dbet1) -
          Math::asinh(_bc2 * sbet * cbet / sqrt(Math::sq(_b * sbet) * _ab2 +
                                                Math::sq(_c * cbet) * _ac2 )) );
    }
    /**
     * The \e y projection
     *
     * @param[in] bet &beta; (in degrees)
     * @return \e y
     *
     * &beta; must be in (&minus;180&deg;, 180&deg;].
     **********************************************************************/
    Math::real y(Math::real bet) const {
      using std::abs; using std::sin; using std::cos;
      Math::real
        a = bet * Math::degree(),
        sbet = abs(bet) == 180 ? 0 : sin(a),
        cbet = abs(bet) ==  90 ? 0 : cos(a);
      return y(sbet, cbet);
    }
  };

} // namespace GeographicLib
