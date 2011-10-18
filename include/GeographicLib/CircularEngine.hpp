/**
 * \file CircularEngine.hpp
 * \brief Header for GeographicLib::CircularEngine class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CIRCULARENGINE_HPP)
#define GEOGRAPHICLIB_CIRCULARENGINE_HPP "$Id$"

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Circular Harmonic series
   *
   * Sum a circular harmonic series.
   **********************************************************************/

  class GEOGRAPHIC_EXPORT CircularEngine {
  private:
    typedef Math::real real;
    enum normalization {
      full = SphericalEngine::full,
      schmidt = SphericalEngine::schmidt,
    };
    int _M;
    bool _gradp;
    normalization _norm;
    real _scale, _a, _r, _u, _t;
    std::vector<real> _wc, _ws, _wrc, _wrs, _wtc, _wts;
    real _q, _uq, _uq2;

    Math::real Value(bool gradp, real coslam, real sinlam,
                     real& gradx, real& grady, real& gradz) const;

    static inline void cossin(real x, real& cosx, real& sinx) {
      x = x >= 180 ? x - 360 : (x < -180 ? x + 360 : x);
      real xi = x * Math::degree<real>();
      cosx = std::abs(x) ==   90 ? 0 : cos(xi);
      sinx =          x  == -180 ? 0 : sin(xi);
    }

  public:

    CircularEngine(int M, bool gradp, SphericalEngine::normalization norm,
                   real scale, real a, real r, real u, real t)
      : _M(M)
      , _gradp(gradp)
      , _norm(normalization(norm))
      , _scale(scale)
      , _a(a)
      , _r(r)
      , _u(u)
      , _t(t)
      , _wc(std::vector<real>(_M + 1, 0))
      , _ws(std::vector<real>(_M + 1, 0))
      , _wrc(std::vector<real>(_gradp ? _M + 1 : 0, 0))
      , _wrs(std::vector<real>(_gradp ? _M + 1 : 0, 0))
      , _wtc(std::vector<real>(_gradp ? _M + 1 : 0, 0))
      , _wts(std::vector<real>(_gradp ? _M + 1 : 0, 0))
      {
        _q = _a / _r;
        _uq = _u * _q;
        _uq2 = Math::sq(_uq);
      }
    void SetCoeff(int m, real wc, real ws)
    { _wc[m] = wc; _ws[m] = ws; }
    void SetCoeff(int m, real wc, real ws,
                  real wrc, real wrs, real wtc, real wts) {
      _wc[m] = wc; _ws[m] = ws;
      if (_gradp) {
        _wrc[m] = wrc; _wrs[m] = wrs;
        _wtc[m] = wtc; _wts[m] = wts;
      }
    }
    Math::real operator()(real lam) const {
      real coslam, sinlam;
      cossin(lam, coslam, sinlam);
      return (*this)(coslam, sinlam);
    }
    Math::real operator()(real coslam, real sinlam) const {
      real dummy;
      return Value(false, coslam, sinlam, dummy, dummy, dummy);
    }
    Math::real operator()(real lam,
                          real& gradx, real& grady, real& gradz) const {
      real coslam, sinlam;
      cossin(lam, coslam, sinlam);
      return (*this)(coslam, sinlam, gradx, grady, gradz);
    }
    Math::real operator()(real coslam, real sinlam,
                          real& gradx, real& grady, real& gradz) const {
      return Value(true, coslam, sinlam, gradx, grady, gradz);
    }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_CIRCULARENGINE_HPP
