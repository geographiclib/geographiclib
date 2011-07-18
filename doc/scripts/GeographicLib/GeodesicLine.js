/**
 * GeodesicLine.js
 * Transcription of GeodesicLine.[ch]pp into javascript.
 *
 * See the documentation for the C++ class.  The conversion is a literal
 * conversion from C++.
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * $Id$
 **********************************************************************/

// Load AFTER GeographicLib/Math.js and GeographicLib/Geodesic.js

(function() {
  var g = GeographicLib.Geodesic;
  var l = GeographicLib.GeodesicLine;
  var m = GeographicLib.Math;

  l.GeodesicLine = function(geod, lat1, lon1, azi1, caps) {
    this._a = geod._a;
    this._f = geod._f;
    this._b = geod._b;
    this._c2 = geod._c2;
    this._f1 = geod._f1;
    this._caps = !caps ? g.ALL : (caps | g.LATITUDE | g.AZIMUTH);
    azi1 = g.AngNormalize(azi1);
    // Guard against underflow in salp0
    azi1 = g.AngRound(azi1);
    lon1 = g.AngNormalize(lon1);
    this._lat1 = lat1;
    this._lon1 = lon1;
    this._azi1 = azi1;
    // alp1 is in [0, pi]
    var alp1 = azi1 * m.degree;
    // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
    // problems directly than to skirt them.
    this._salp1 =     azi1  == -180 ? 0 : Math.sin(alp1);
    this._calp1 = Math.abs(azi1) ==   90 ? 0 : Math.cos(alp1);
    var cbet1, sbet1, phi;
    phi = lat1 * m.degree;
    // Ensure cbet1 = +epsilon at poles
    sbet1 = this._f1 * Math.sin(phi);
    cbet1 = Math.abs(lat1) == 90 ? g.tiny_ : Math.cos(phi);
    var t = m.hypot(sbet1, cbet1);
    sbet1 /= t; cbet1 /= t;
    // SinCosNorm(sbet1, cbet1);

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    this._salp0 = this._salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    this._calp0 = m.hypot(this._calp1, this._salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    this._ssig1 = sbet1; this._somg1 = this._salp0 * sbet1;
    this._csig1 = this._comg1 =
      sbet1 != 0 || this._calp1 != 0 ? cbet1 * this._calp1 : 1;
    var t = m.hypot(this._ssig1, this._csig1);
    this._ssig1 /= t; this._csig1 /= t;
    // SinCosNorm(this._ssig1, this._csig1); // sig1 in (-pi, pi]
    var t = m.hypot(this._somg1, this._comg1);
    this._somg1 /= t; this._comg1 /= t;
    // SinCosNorm(this._somg1, this._comg1);

    this._k2 = m.sq(this._calp0) * geod._ep2;
    var eps = this._k2 / (2 * (1 + Math.sqrt(1 + this._k2)) + this._k2);

    if (this._caps & g.CAP_C1) {
      this._A1m1 = g.A1m1f(eps);
      this._C1a = new Array(g.nC1_ + 1);
      g.C1f(eps, this._C1a);
      this._B11 = g.SinCosSeries(true, this._ssig1, this._csig1,
                                 this._C1a, g.nC1_);
      var s = Math.sin(this._B11), c = Math.cos(this._B11);
      // tau1 = sig1 + B11
      this._stau1 = this._ssig1 * c + this._csig1 * s;
      this._ctau1 = this._csig1 * c - this._ssig1 * s;
      // Not necessary because C1pa reverts C1a
      //    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);
    }

    if (this._caps & g.CAP_C1p) {
      this._C1pa = new Array(g.nC1p_ + 1),
      g.C1pf(eps, this._C1pa);
    }

    if (this._caps & g.CAP_C2) {
      this._A2m1 = g.A2m1f(eps);
      this._C2a = new Array(g.nC2_ + 1);
      g.C2f(eps, this._C2a);
      this._B21 = g.SinCosSeries(true, this._ssig1, this._csig1,
                                 this._C2a, g.nC2_);
    }

    if (this._caps & g.CAP_C3) {
      this._C3a = new Array(g.nC3_);
      geod.C3f(eps, this._C3a);
      this._A3c = -this._f * this._salp0 * geod.A3f(eps);
      this._B31 = g.SinCosSeries(true, this._ssig1, this._csig1,
                                 this._C3a, g.nC3_-1);
    }

    if (this._caps & g.CAP_C4) {
      this._C4a = new Array(g.nC4_);    // all the elements of _C4a are used
      geod.C4f(this._k2, this._C4a);
      // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      this._A4 = m.sq(this._a) * this._calp0 * this._salp0 * geod._e2;
      this._B41 = g.SinCosSeries(false, this._ssig1, this._csig1,
                                 this._C4a, g.nC4_);
    }
  }

  // return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
  l.GeodesicLine.prototype.GenPosition = function(arcmode, s12_a12,
                                                  outmask) {
    var vals = {};
    outmask &= this._caps & g.OUT_ALL;
    if (!( arcmode || (this._caps & g.DISTANCE_IN & g.OUT_ALL) )) {
      // Uninitialized or impossible distance calculation requested
      vals.a12 = Number.NaN();
      return vals;
    }

    // Avoid warning about uninitialized B12.
    var sig12, ssig12, csig12, B12 = 0, AB1 = 0;
    if (arcmode) {
      // Interpret s12_a12 as spherical arc length
      sig12 = s12_a12 * m.degree;
      var s12a = Math.abs(s12_a12);
      s12a -= 180 * Math.floor(s12a / 180);
      ssig12 = s12a ==  0 ? 0 : Math.sin(sig12);
      csig12 = s12a == 90 ? 0 : Math.cos(sig12);
    } else {
      // Interpret s12_a12 as distance
      var
      tau12 = s12_a12 / (this._b * (1 + this._A1m1)),
      s = Math.sin(tau12),
      c = Math.cos(tau12);
      // tau2 = tau1 + tau12
      B12 = - g.SinCosSeries(true, this._stau1 * c + this._ctau1 * s,
                                   this._ctau1 * c - this._stau1 * s,
                                   this._C1pa, g.nC1p_);
      sig12 = tau12 - (B12 - this._B11);
      ssig12 = Math.sin(sig12);
      csig12 = Math.cos(sig12);
    }

    var omg12, lam12, lon12;
    var ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = this._ssig1 * csig12 + this._csig1 * ssig12;
    csig2 = this._csig1 * csig12 - this._ssig1 * ssig12;
    if (outmask & (g.DISTANCE | g.REDUCEDLENGTH | g.GEODESICSCALE)) {
      if (arcmode)
        B12 = g.SinCosSeries(true, ssig2, csig2, this._C1a, g.nC1_);
      AB1 = (1 + this._A1m1) * (B12 - this._B11);
    }
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = this._calp0 * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = m.hypot(this._salp0, this._calp0 * csig2);
    if (cbet2 == 0)
      // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
      cbet2 = csig2 = g.tiny_;
    // tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = this._salp0 * ssig2; comg2 = csig2;  // No need to normalize
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = this._salp0; calp2 = this._calp0 * csig2; // No need to normalize
    // omg12 = omg2 - omg1
    omg12 = Math.atan2(somg2 * this._comg1 - comg2 * this._somg1,
                       comg2 * this._comg1 + somg2 * this._somg1);

    if (outmask & g.DISTANCE)
      vals.s12 = arcmode ? this._b * ((1 + this._A1m1) * sig12 + AB1) : s12_a12;

    if (outmask & g.LONGITUDE) {
      lam12 = omg12 + this._A3c *
        ( sig12 + (g.SinCosSeries(true, ssig2, csig2, this._C3a, g.nC3_-1)
                   - this._B31));
      lon12 = lam12 / m.degree;
      // Can't use AngNormalize because longitude might have wrapped multiple
      // times.
      lon12 = lon12 - 360 * Math.floor(lon12/360 + 0.5);
      vals.lon2 = g.AngNormalize(this._lon1 + lon12);
    }

    if (outmask & g.LATITUDE)
      vals.lat2 = Math.atan2(sbet2, this._f1 * cbet2) / m.degree;

    if (outmask & g.AZIMUTH)
      // minus signs give range [-180, 180). 0- converts -0 to +0.
      vals.azi2 = 0 - Math.atan2(-salp2, calp2) / m.degree;

    if (outmask & (g.REDUCEDLENGTH | g.GEODESICSCALE)) {
      var
      ssig1sq = m.sq(this._ssig1),
      ssig2sq = m.sq( ssig2),
      w1 = Math.sqrt(1 + this._k2 * ssig1sq),
      w2 = Math.sqrt(1 + this._k2 * ssig2sq),
      B22 = g.SinCosSeries(true, ssig2, csig2, this._C2a, g.nC2_),
      AB2 = (1 + this._A2m1) * (B22 - this._B21),
      J12 = (this._A1m1 - this._A2m1) * sig12 + (AB1 - AB2);
      if (outmask & g.REDUCEDLENGTH)
        // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
        // accurate cancellation in the case of coincident points.
        vals.m12 = this._b * ((w2 * (this._csig1 * ssig2) -
                               w1 * (this._ssig1 * csig2))
                              - this._csig1 * csig2 * J12);
      if (outmask & g.GEODESICSCALE) {
        vals.M12 = csig12 + (this._k2 * (ssig2sq - ssig1sq) *  ssig2 / (w1 + w2)
                             - csig2 * J12) * this._ssig1 / w1;
        vals.M21 = csig12 - (this._k2 * (ssig2sq - ssig1sq) * this._ssig1 /
                             (w1 + w2)
                             - this._csig1 * J12) * ssig2 / w2;
      }
    }

    if (outmask & g.AREA) {
      var
      B42 = g.SinCosSeries(false, ssig2, csig2, this._C4a, g.nC4_);
      var salp12, calp12;
      if (this._calp0 == 0 || this._salp0 == 0) {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalized
        salp12 = salp2 * this._calp1 - calp2 * this._salp1;
        calp12 = calp2 * this._calp1 + salp2 * this._salp1;
        // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        // salp12 = -0 and alp12 = -180.  However this depends on the sign being
        // attached to 0 correctly.  The following ensures the correct behavior.
        if (salp12 == 0 && calp12 < 0) {
          salp12 = g.tiny_ * this._calp1;
          calp12 = -1;
        }
      } else {
        // tan(alp) = tan(alp0) * sec(sig)
        // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
        // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
        // If csig12 > 0, write
        //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
        // else
        //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
        // No need to normalize
        salp12 = this._calp0 * this._salp0 *
          (csig12 <= 0 ? this._csig1 * (1 - csig12) + ssig12 * this._ssig1 :
           ssig12 * (this._csig1 * ssig12 / (1 + csig12) + this._ssig1));
        calp12 = m.sq(this._salp0) + m.sq(this._calp0) * this._csig1 * csig2;
      }
      vals.S12 = this._c2 * Math.atan2(salp12, calp12) +
        this._A4 * (B42 - this._B41);
    }

    vals.a12 = arcmode ? s12_a12 : sig12 / m.degree;
    return vals;
  }


})();
