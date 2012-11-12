function [s12, azi1, azi2] = geoddistance(lat1, lon1, lat2, lon2, ellipsoid)
%DISTANCE  Distance between points on ellipsoid
%
%   [S12, AZI1, AZI2] = GEODDISTANCE(LAT1, LON1, LAT2, LON2, ELLIPSOID)
%   computes the geodesic distance and azimuth assuming that the points lie
%   on the reference ellipsoid defined by the input ELLIPSOID.  The input
%   latitudes and longitudes, LAT1, LON1, LAT2, LON2, can be scalars or
%   arrays of equal size and must be expressed in degrees.  The ELLIPSOID
%   vector is of the form [a, e], where a is the equatorial radius, e is
%   the eccentricity e = sqrt(a^2 - b^2)/a, and b is the polar semi-axis.
%   The output S12 is expressed in the same distance units as the
%   equatorial radius.  AZI1 and AZI2 are the forward azimuths at the end
%   points in degrees.
%
%   When given a combination of scalar and array inputs, the scalar inputs
%   are automatically expanded to match the size of the arrays.
%
%   This is an implementation of the algorithm given in
%
%     C. F. F. Karney
%     Algorithms for geodesics
%     J. Geodesy (2012)
%     http://dx.doi.org/10.1007/s00190-012-0578-z
%
%   The calculations are carried out as expansions in the eccentricity
%   which are accurate for eccentricities typical of the Earth (i.e.,
%   abs(e) < 0.1).  Note that the algorithms are valid also for slightly
%   prolate ellipsoids (b > a), in which case the eccentricity should be
%   specified as a pure imaginary number.
%
%   This function duplicates some of the functionality of the DISTANCE
%   function in the MATLAB mapping toolbox.  Differences are
%
%     * When the ELLIPSOID argument is omitted, use the WGS84 ellipsoid.
%     * The azimuths at the end points AZI1 and AZI2 are returned.
%     * The solution is accurate to round-off error.
%     * The algorithm converges for all pairs of input points.
%
%   This is the solution of the so-called inverse geodesic problem.  The
%   direct geodesic problem is solved by GEODRECKON.
%
%   The MATLAB implementation is a transcription of the C++ version in
%   GeographicLib http://geographiclib.sf.net.  Note the C++ version has a
%   few additional capabilities (e.g., computing also the reduced length
%   and the ellipsoidal area).  These capabilities are accessible from
%   MATLAB using the wrapper function, GEODESICINVERSE.
%
%   See also GEODESICINVERSE, GEODRECKON.

% Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
% under the MIT/X11 License.  For more information, see
% http://geographiclib.sourceforge.net/
%
% This is a straightforward transcription of the C++ implementation in
% GeographicLib and the C++ source should be consulted for additional
% documentation.  This is a vector implementation and the results returned
% with array arguments are identical to those obtained with multiple calls
% with scalar arguments.  The biggest change was to eliminate the branching
% to allow a vectorized solution.
%
% This file was distributed with GeographicLib 1.27.

  try
    Z = lat1 + lon1 + lat2 + lon2;
    S = size(Z);
    Z = zeros(S);
    lat1 = lat1 + Z; lon1 = lon1 + Z;
    lat2 = lat2 + Z; lon2 = lon2 + Z;
    Z = Z(:);
  catch err
    error('lat1, lon1, s12, azi1 have incompatible sizes');
  end

  degree = pi/180;
  tiny = sqrt(realmin);
  tol0 = eps;
  tol1 = 200 * tol0;
  maxit = 30;
  bisection = -log2(eps) + 1 + 10;

  if nargin < 5,
    a = 6378137;
    f = 1/298.257223563;
    e2 = f * (2 - f);
  else
    a = ellipsoid(1);
    e2 = ellipsoid(2)^2;
    f = e2 / (1 + sqrt(1 - e2));
  end
  f1 = 1 - f;
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  b = a * f1;

  A3x = A3coeff(n);
  C3x = C3coeff(n);

  lon1 = AngNormalize(lon1(:));
  lon12 = AngNormalize(AngNormalize(lon2(:)) - lon1);
  lon12 = AngRound(lon12);
  lonsign = 2 * (lon12 >= 0) - 1;
  lon12 = lonsign .* lon12;
  lonsign(lon12 == 180) = 1;
  lat1 = AngRound(lat1(:));
  lat2 = AngRound(lat2(:));
  swapp = 2 * (abs(lat1) >= abs(lat2)) - 1;
  lonsign(swapp < 0) = - lonsign(swapp < 0);
  [lat1(swapp < 0), lat2(swapp < 0)] = swap(lat1(swapp < 0), lat2(swapp < 0));

  latsign = 2 * (lat1 < 0) - 1;
  lat1 = latsign .* lat1;
  lat2 = latsign .* lat2;

  phi = lat1 * degree;
  sbet1 = f1 * sin(phi); cbet1 = cos(phi); cbet1(lat1 == -90) = tiny;
  [sbet1, cbet1] = SinCosNorm(sbet1, cbet1);

  phi = lat2 * degree;
  sbet2 = f1 * sin(phi); cbet2 = cos(phi); cbet2(abs(lat2) == 90) = tiny;
  [sbet2, cbet2] = SinCosNorm(sbet2, cbet2);

  c = cbet1 < -sbet1 & cbet2 == cbet1;
  sbet2(c) = (2 * (sbet2(c) < 0) - 1) .* sbet1(c);
  c = ~(cbet1 < -sbet1) & abs(sbet2) == - sbet1;
  cbet2(c) = cbet1(c);

  dn1 = sqrt(1 + ep2 * sbet1.^2);
  dn2 = sqrt(1 + ep2 * sbet2.^2);
  lam12 = lon12 * degree;
  slam12 = sin(lam12); slam12(lon12 == 180) = 0; clam12 = cos(lam12);

  sig12 = Z; ssig1 = Z; csig1 = Z; ssig2 = Z; csig2 = Z;
  calp1 = Z; salp1 = Z; calp2 = Z; salp2 = Z;
  s12 = Z; m12 = Z;

  m = lat1 == -90 | slam12 == 0;

  if any(m),
    calp1(m) = clam12(m); salp1(m) = slam12(m);
    calp2(m) = 1; salp2(m) = 0;

    ssig1(m) = sbet1(m); csig1(m) = calp1(m) .* cbet1(m);
    ssig2(m) = sbet2(m); csig2(m) = calp2(m) .* cbet2(m);

    sig12(m) = atan2(max(csig1(m) .* ssig2(m) - ssig1(m) .* csig2(m), 0), ...
                     csig1(m) .* csig2(m) + ssig1(m) .* ssig2(m));

    [s12(m), m12(m), ~] = Lengths(n, sig12(m), ...
                                  ssig1(m), csig1(m), dn1(m), ...
                                  ssig2(m), csig2(m), dn2(m));
    m = m & (sig12 < 1 | m12 >= 0);
    s12(m) = s12(m) * b;
  end

  eq = ~m & sbet1 == 0;
  if f > 0,
    eq = eq & lam12 < pi - f * pi;
  end
  calp1(eq) = 0; calp2(eq) = 0; salp1(eq) = 1; salp2(eq) = 1;
  s12(eq) = a * lam12(eq);

  g = ~eq & ~m;

  [sig12(g), salp1(g), calp1(g), salp2(g), calp2(g)] = ...
      InverseStart(sbet1(g), cbet1(g), dn1(g), sbet2(g), cbet2(g), dn2(g), ...
                   lam12(g), f, A3x);

  s = g & sig12 >= 0;
  s12(s) = b * sig12(s) .* (dn1(s) + dn2(s)) / 2;

  g = g & sig12 < 0;

  ov = Z;
  salp1a = Z + tiny; calp1a = Z + 1;
  salp1b = Z + tiny; calp1b = Z - 1;
  ssig1 = Z; csig1 = Z; ssig2 = Z; csig2 = Z;
  epsi = Z; v = Z; dv = Z;
  trip = Z; numit = Z;
  gsave = g;
  for k = 1 : maxit,
    if ~any(g),
      break;
    end
    numit(g) = k;
    [v(g), dv(g), ...
     salp2(g), calp2(g), sig12(g), ...
     ssig1(g), csig1(g), ssig2(g), csig2(g), epsi(g)] = ...
        Lambda12(sbet1(g), cbet1(g), dn1(g), ...
                 sbet2(g), cbet2(g), dn2(g), ...
                 salp1(g), calp1(g), f, A3x, C3x);
    v = v - lam12;

    c = g & v > 0 & calp1 ./ salp1 > calp1b ./ salp1b;
    salp1b(c) = salp1(c); calp1b(c) = calp1(c);
    c = g & v < 0 & calp1 ./ salp1 < calp1a ./ salp1a;
    salp1a(c) = salp1(c); calp1a(c) = calp1(c);

    c = g & (~(abs(v) > tiny) | ~(trip < 1));
    numit(c & ~(abs(v) <= max(tol1, ov))) = maxit;
    g = g & ~c;

    dalp1 = -v ./ dv;
    sdalp1 = sin(dalp1); cdalp1 = cos(dalp1);
    nsalp1 = salp1 .* cdalp1 + calp1 .* sdalp1;
    calp1 = calp1 .* cdalp1 - salp1 .* sdalp1;
    salp1 = nsalp1;
    c = g & ~(abs(v) >= tol1 & v.^2 >= ov * tol0);
    trip(c) = trip(c) + 1;
    ov = abs(v);

    c = g & ~(dv > 0 & nsalp1 > 0 & abs(dalp1) < pi);
    salp1(c) = (salp1a(c) + salp1b(c))/2;
    calp1(c) = (calp1a(c) + calp1b(c))/2;
    trip(c) = 0;
    ov(c) = 0;

    [salp1(g), calp1(g)] = SinCosNorm(salp1(g), calp1(g));
  end

  g = numit >= maxit;
  for k = maxit + 1 : maxit + bisection,
    if ~any(g),
      break;
    end
    numit(g) = k;
    salp1(g) = (salp1a(g) + salp1b(g))/2;
    calp1(g) = (calp1a(g) + calp1b(g))/2;
    [salp1(g), calp1(g)] = SinCosNorm(salp1(g), calp1(g));
    c = (abs(salp1 - salp1b) < tol0 & calp1 - calp1b < tol0) | ...
        (abs(salp1a - salp1) < tol0 & calp1a - calp1 < tol0);
    g = g & ~c;
    [v(g), dv(g), ...
     salp2(g), calp2(g), sig12(g), ...
     ssig1(g), csig1(g), ssig2(g), csig2(g), epsi(g)] = ...
        Lambda12(sbet1(g), cbet1(g), dn1(g), ...
                 sbet2(g), cbet2(g), dn2(g), ...
                 salp1(g), calp1(g), f, A3x, C3x);
    v - v - lam12;

    c = abs(v) <= 2 * tol0;
    g = g & ~c;
    c = v > 0;
    salp1b(c) = salp1(c); calp1b(c) = calp1(c);
    c = ~c;
    salp1a(c) = salp1(c); calp1a(c) = calp1(c);
  end

  g = gsave;
  [s12(g), ~, ~] = Lengths(epsi(g), sig12(g), ...
                           ssig1(g), csig1(g), dn1(g), ...
                           ssig2(g), csig2(g), dn2(g));

  s12(g) = s12(g) * b;

  [salp1(swapp<0), salp2(swapp<0)] = swap(salp1(swapp<0), salp2(swapp<0));
  [calp1(swapp<0), calp2(swapp<0)] = swap(calp1(swapp<0), calp2(swapp<0));

  salp1 = salp1 .* swapp .* lonsign; calp1 = calp1 .* swapp .* latsign;
  salp2 = salp2 .* swapp .* lonsign; calp2 = calp2 .* swapp .* latsign;

  azi1 = 0 - atan2(-salp1, calp1) / degree;
  azi2 = 0 - atan2(-salp2, calp2) / degree;

  g = numit >= maxit + bisection;
  s12(g) = NaN; azi1(g) = NaN; azi2(g) = NaN;

  s12 = reshape(s12, S); azi1 = reshape(azi1, S); azi2 = reshape(azi2, S);
end

%% UTILITIES

function z = cvmgt(x, y, p)
%CVMGT  Conditional merge of two vectors
%
%   Z = CVMGT(X, Y, P) return a vector Z whose elements are X if P is true
%   and Y otherwise.  P, X, and Y should be the same shape except that X
%   and Y may be scalars.  CVMGT stands for conditional vector merge true
%   (an intrinsic function for the Cray fortran compiler).  It implements
%   the C++ statement
%
%     Z = P ? X : Y;

  z = zeros(size(p));
  if isscalar(x),
    z(p) = x;
  else
    z(p) = x(p);
  end
  if isscalar(y),
    z(~p) = y;
  else
    z(~p) = y(~p);
  end
end

function [a, b] = swap(x, y)
%SWAP  Swap two variables.
%
%   [A, B] = SWAP(X, Y) sets A to Y and B to X.

  a = y;
  b = x;
end

function y = cbrt(x)
%CBRT   The real cube root
%   CBRT(X) is the real cube root of X (assuming X is real).  X
%   can be any shape.

  y = abs(x).^(1/3);
  y(x < 0) = -y(x < 0);
end

function [sinx, cosx] = SinCosNorm(sinx, cosx)
%SINCOSNORM  Normalize sinx and cosx
%
%  [sinx, cosx] = SINCOSNORM(sinx, cosx) normalize sinx and cosx so that
%  sinx^2 + cosx^2 = 1.  sinx and cosx can be any shape.

  r = hypot(sinx, cosx);
  sinx = sinx ./ r;
  cosx = cosx ./ r;
end

function y = SinSeries(sinx, cosx, c)
%SINSERIES  Evaluate a sine series using Clenshaw summation
%
%  Y = SINSERIES(SINX, COSX, C) evaluate
%
%    y = sum(c[i] * sin( 2*i * x), i, 1, n)
%
%  where n is the size of c.  x is given via its sine and cosine in SINX and
%  COSX.  SINX, COSX, and Y are K x 1 arrays.  C is a K x N array.

  if size(sinx, 1) == 0,
    y = [];
    return;
  end
  n = size(c, 2);
  ar = 2 * (cosx - sinx) .* (cosx + sinx); % 2 * cos(2 * x)
  y1 = zeros(size(sinx, 1), 1);
  if mod(n, 2),
    y0 = c(:, n);
    n = n - 1;
  else
    y0 = y1;
  end

  for k = n : -2 : 1,
    y1 = ar .* y0 - y1 + c(:, k);
    y0 = ar .* y1 - y0 + c(:, k-1);
  end
  y = 2 * sinx .* cosx .* y0;
end

function x = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range [-180, 180)
%
%  X = ANGNORMALIZE(X) reduces angles in [-540, 540) to the range
%  [-180, 180).  X can be any shape.

  x(x >= 180) = x(x >= 180) - 360;
  x(x < -180) = x(x < -180) + 360;
end

function y =  AngRound(x)
%ANGROUND  Round tiny values so that tiny values become zero.
%
%  Y = ANGROUND(X) rounds X by adding and subtracting 1/16 to it if it is
%  small.  X and Y can be any shape.

  z = 1/16;
  y = abs(x);
  y(y < z) = z - (z - y(y < z));
  y(x < 0) = -y(x < 0);
end

%% HELPER FUNCTIONS FOR THE INVERSE PROBLEM

function [sig12, salp1, calp1, salp2, calp2] = ...
      InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, f, A3x)
%INVERSESTART  Compute a starting point for Newton's method

  N = size(sbet1, 1);
  f1 = 1 - f;
  e2 = f * (2 - f);
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  tol0 = eps;
  tol1 = 200 * tol0;
  tol2 = sqrt(eps);
  etol2 = 10 * tol2 / max(0.1, sqrt(abs(e2)));
  xthresh = 1000 * tol2;

  sig12 = - ones(N, 1); salp2 = NaN(N, 1); calp2 = NaN(N, 1);
  sbet12 = sbet2 .* cbet1 - cbet2 .* sbet1;
  cbet12 = cbet2 .* cbet1 + sbet2 .* sbet1;
  sbet12a = sbet2 .* cbet1 + cbet2 .* sbet1;
  s = cbet12 >= 0 & sbet12 < 0.5 & lam12 <= pi / 6;
  omg12 = lam12;
  omg12(s) = omg12(s) ./ (f1 * (dn1(s) + dn2(s)) / 2);
  somg12 = sin(omg12); comg12 = cos(omg12);

  salp1 = cbet2 .* somg12;
  t = cbet2 .* sbet1 .* somg12.^2;
  calp1 = cvmgt(sbet12  + t ./ (1 + comg12), ...
                sbet12a - t ./ (1 - comg12), ...
                comg12 >= 0);

  ssig12 = hypot(salp1, calp1);
  csig12 = sbet1 .* sbet2 + cbet1 .* cbet2 .* comg12;

  s = s & ssig12 < etol2;
  salp2(s) = cbet1(s) .* somg12(s);
  calp2(s) = sbet12(s) - cbet1(s) .* sbet2(s) .* somg12(s).^2 ./ ...
      (1 + comg12(s));
  [salp2, calp2] = SinCosNorm(salp2, calp2);
  sig12(s) = atan2(ssig12(s), csig12(s));

  s = ~(s | abs(n) > 0.1 | csig12 >= 0 | ssig12 >= 6 * abs(n) * pi * cbet1.^2);

  if any(s),
    if f >= 0,
      k2 = sbet1(s).^2 * ep2;
      epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
      lamscale = f * cbet1(s) .* A3f(epsi, A3x) * pi;
      betscale = lamscale .* cbet1(s);
      x = (lam12(s) - pi) ./ lamscale;
      y = sbet12a(s) ./ betscale;
    else
      cbet12a = cbet2(s) .* cbet1(s) - sbet2(s) .* sbet1(s);
      bet12a = atan2(sbet12a(s), cbet12a);
      [~, m12b, m0] = Lengths(n, pi + bet12a, ...
                                  sbet1(s), -cbet1(s), dn1(s), ...
                                  sbet2(s), cbet2(s), dn2(s));
      x = -1 + m12b ./ (cbet1(s) .* cbet2(s) .* m0 * pi);
      betscale = cvmgt(sbet12a(s) ./ x, - f * cbet1(s).^2 * pi, x < -0.01);
      lamscale = betscale ./ cbet1(s);
      y = (lam12(s) - pi) ./ lamscale;
    end
    k = Astroid(x, y);
    if f >= 0,
      omg12a = -x .* k ./ (1 + k);
    else
      omg12a = -y .* (1 + k) ./ k;
    end
    omg12a = lamscale .* omg12a;
    somg12 = sin(omg12a); comg12 = -cos(omg12a);
    salp1(s) = cbet2(s) .* somg12;
    calp1(s) = sbet12a(s) - cbet2(s) .* sbet1(s) .* somg12.^2 ./ (1 - comg12);

    str = y > -tol1 & x > -1 - xthresh;
    if any(str),
      salp1s = salp1(s); calp1s = calp1(s);
      if f >= 0,
        salp1s(str) = min(1, -x(str));
        calp1s(str) = -sqrt(1 - salp1s(str).^2);
      else
        calp1s(str) = max(cvmgt(0, -1, x(str) > -tol1), x(str));
        salp1s(str) = sqrt(1 - calp1s(str).^2);
      end
      salp1(s) = salp1s; calp1(s) = calp1s;
    end
  end

  calp1(salp1 <= 0) = 0; salp1(salp1 <= 0) = 1;
  [salp1, calp1] = SinCosNorm(salp1, calp1);
end

function k = Astroid(x, y)
% ASTROID  Solve the astroid equation
%
%   K = ASTROID(X, Y) solves the quartic polynomial Eq. (55)
%
%     K^4 + 2 * K^3 - (X^2 + Y^2 - 1) * K^2 - 2*Y^2 * K - Y^2 = 0
%
%   for the positive root K.  X and Y are column vectors of the same size
%   and the returned value K has the same size.

  k = zeros(size(x, 1), 1);
  p = x.^2;
  q = y.^2;
  r = (p + q - 1) / 6;
  fl1 = ~(q == 0 & r <= 0);
  p = p(fl1);
  q = q(fl1);
  r = r(fl1);
  S = p .* q / 4;
  r2 = r.^2;
  r3 = r .* r2;
  disc = S .* (S + 2 * r3);
  u = r;
  fl2 = disc >= 0;
  T3 = S(fl2) + r3(fl2);
  T3 = T3 + (1 - 2 * (T3 < 0)) .* sqrt(disc(fl2));
  T = cbrt(T3);
  u(fl2) = u(fl2) + T + cvmgt(r2(fl2) ./ T, 0, T ~= 0);
  ang = atan2(sqrt(-disc(~fl2)), -(S(~fl2) + r3(~fl2)));
  u(~fl2) = u(~fl2) + 2 * r(~fl2) .* cos(ang / 3);
  v = sqrt(u.^2 + q);
  uv = u + v;
  fl2 = u < 0;
  uv(fl2) = q(fl2) ./ (v(fl2) - u(fl2));
  w = (uv - q) ./ (2 * v);
  k(fl1) = uv ./ (sqrt(uv + w.^2) + w);
end

function [s12b, m12b, m0] = Lengths(epsi, sig12, ...
                                    ssig1, csig1, dn1, ssig2, csig2, dn2)
%LENGTHS  Compute various lengths associate with a geodesic

  if size(sig12, 1) == 0,
    s12b = [];
    m12b = [];
    m0 = [];
    return;
  end

  C1a = C1f(epsi);
  C2a = C2f(epsi);
  A1m1 = A1m1f(epsi);
  AB1 = (1 + A1m1) .* (SinSeries(ssig2, csig2, C1a) - ...
                       SinSeries(ssig1, csig1, C1a));
  A2m1 = A2m1f(epsi);
  AB2 = (1 + A2m1) .* (SinSeries(ssig2, csig2, C2a) - ...
                       SinSeries(ssig1, csig1, C2a));
  m0 = A1m1 - A2m1;
  J12 = m0 .* sig12 + (AB1 - AB2);
  m12b = dn2 .* (csig1 .* ssig2) - dn1 .* (ssig1 .* csig2) - ...
         csig1 .* csig2 .* J12;
  s12b = (1 + A1m1) .* sig12 + AB1;
end

function [lam12, dlam12, ...
          salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, epsi] = ...
    Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, f, A3x, C3x)
%LAMBDA12  Solve the hybrid problem

  tiny = sqrt(realmin);
  f1 = 1 - f;
  e2 = f * (2 - f);
  ep2 = e2 / (1 - e2);

  calp1(sbet1 == 0 & calp1 == 0) = -tiny;

  salp0 = salp1 .* cbet1;
  calp0 = hypot(calp1, salp1 .* sbet1);

  ssig1 = sbet1; somg1 = salp0 .* sbet1;
  csig1 = calp1 .* cbet1; comg1 = csig1;
  [ssig1, csig1] = SinCosNorm(ssig1, csig1);

  salp2 = cvmgt(salp0 ./ cbet2, salp1, cbet2 ~= cbet1);
  calp2 = cvmgt(sqrt((calp1 .* cbet1).^2 + ...
                     cvmgt((cbet2 - cbet1) .* (cbet1 + cbet2), ...
                           (sbet1 - sbet2) .* (sbet1 + sbet2), ...
                           cbet1 < -sbet1)) ./ cbet2, ...
                abs(calp1), cbet2 ~= cbet1 | abs(sbet2) ~= -sbet1);
  ssig2 = sbet2; somg2 = salp0 .* sbet2;
  csig2 = calp2 .* cbet2;  comg2 = csig2;
  [ssig2, csig2] = SinCosNorm(ssig2, csig2);

  sig12 = atan2(max(csig1 .* ssig2 - ssig1 .* csig2, 0), ...
                csig1 .* csig2 + ssig1 .* ssig2);

  omg12 = atan2(max(comg1 .* somg2 - somg1 .* comg2, 0), ...
                comg1 .* comg2 + somg1 .* somg2);
  k2 = calp0.^2 * ep2;
  epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
  C3a = C3f(epsi, C3x);
  B312 = SinSeries(ssig2, csig2, C3a) - SinSeries(ssig1, csig1, C3a);
  h0 = -f * A3f(epsi, A3x);
  domg12 = salp0 .* h0 .* (sig12 + B312);
  lam12 = omg12 + domg12;

  [~, dlam12, ~] = Lengths(epsi, sig12, ...
                           ssig1, csig1, dn1, ssig2, csig2, dn2);
  dlam12 = dlam12 .* f1 ./ (calp2 .* cbet2);
  z = calp2 == 0;
  dlam12(z) = - 2 * f1 .* dn1(z) ./ sbet1(z);
end

%% SERIES FOR THE GEODESIC PROBLEM

function A1m1 = A1m1f(epsi)
%A1M1F  Evaluate A_1 - 1
%
%  A1M1 = A1M1F(EPSI) evaluates A_1 - 1 using Eq. (17).  EPSI and A1M1 are
%  K x 1 arrays.

  eps2 = epsi.^2;
  t = eps2.*(eps2.*(eps2+4)+64)/256;
  A1m1 = (t + epsi) ./ (1 - epsi);
end

function C1 = C1f(epsi)
%C1F  Evaluate C_{1,k}
%
%  C1 = C1F(EPSI) evaluates C_{1,l} using Eq. (18).  EPSI is an
%  K x 1 array and C1 is a K x 6 array.

  nC1 = 6;
  C1 = zeros(size(epsi, 1), nC1);
  eps2 = epsi.^2;
  d = epsi;
  C1(:,1) = d.*((6-eps2).*eps2-16)/32;
  d = d.*epsi;
  C1(:,2) = d.*((64-9*eps2).*eps2-128)/2048;
  d = d.*epsi;
  C1(:,3) = d.*(9*eps2-16)/768;
  d = d.*epsi;
  C1(:,4) = d.*(3*eps2-5)/512;
  d = d.*epsi;
  C1(:,5) = -7*d/1280;
  d = d.*epsi;
  C1(:,6) = -7*d/2048;
end

function A2m1 =  A2m1f(epsi)
%A2M1F  Evaluate A_2 - 1
%
%  A2M1 = A2M1F(EPSI) evaluates A_2 - 1 using Eq. (42).  EPSI and A2M1 are
%  K x 1 arrays.

  eps2 = epsi.^2;
  t = eps2.*(eps2.*(25*eps2+36)+64)/256;
  A2m1 = t .* (1 - epsi) - epsi;
end

function C2 = C2f(epsi)
%C2F  Evaluate C_{2,k}
%
%  C2 = C2F(EPSI) evaluates C_{2,l} using Eq. (43).  EPSI is an
%  K x 1 array and C2 is a K x 6 array.

  nC2 = 6;
  C2 = zeros(size(epsi, 1), nC2);
  eps2 = epsi.^2;
  d = epsi;
  C2(:,1) = d.*(eps2.*(eps2+2)+16)/32;
  d = d.*epsi;
  C2(:,2) = d.*(eps2.*(35*eps2+64)+384)/2048;
  d = d.*epsi;
  C2(:,3) = d.*(15*eps2+80)/768;
  d = d.*epsi;
  C2(:,4) = d.*(7*eps2+35)/512;
  d = d.*epsi;
  C2(:,5) = 63*d/1280;
  d = d.*epsi;
  C2(:,6) = 77*d/2048;
end

function A3x = A3coeff(n)
%A3COEFF  Evaluate coefficients for A_3
%
%  A3x = A3COEFF(N) evaluates the coefficients of epsilon^l in Eq. (24).  N
%  is a scalar.  A3x is a 1 x 6 array.

  nA3 = 6;
  A3x = zeros(1, nA3);
  A3x(0+1) = 1;
  A3x(1+1) = (n-1)/2;
  A3x(2+1) = (n*(3*n-1)-2)/8;
  A3x(3+1) = ((-n-3)*n-1)/16;
  A3x(4+1) = (-2*n-3)/64;
  A3x(5+1) = -3/128;
end

function A3 = A3f(epsi, A3x)
%A3F  Evaluate A_3
%
%  A3 = A3F(EPSI, A3X) evaluates A_3 using Eq. (24) and the coefficient
%  vector A3X.  EPSI and A3 are K x 1 arrays.  A3X is a 1 x 6 array.

  nA3 = 6;
  A3 = zeros(size(epsi, 1), 1);
  for i = nA3 : -1 : 1,
    A3 = epsi .* A3 + A3x(i);
  end
end

function C3x = C3coeff(n)
%C3COEFF  Evaluate coefficients for C_3
%
%  C3x = C3COEFF(N) evaluates the coefficients of epsilon^l in Eq. (25).  N
%  is a scalar.  A3x is a 1 x 15 array.

  nC3 = 6;
  nC3x = (nC3 * (nC3 - 1)) / 2;
  C3x = zeros(1, nC3x);
  C3x(0+1) = (1-n)/4;
  C3x(1+1) = (1-n*n)/8;
  C3x(2+1) = ((3-n)*n+3)/64;
  C3x(3+1) = (2*n+5)/128;
  C3x(4+1) = 3/128;
  C3x(5+1) = ((n-3)*n+2)/32;
  C3x(6+1) = ((-3*n-2)*n+3)/64;
  C3x(7+1) = (n+3)/128;
  C3x(8+1) = 5/256;
  C3x(9+1) = (n*(5*n-9)+5)/192;
  C3x(10+1) = (9-10*n)/384;
  C3x(11+1) = 7/512;
  C3x(12+1) = (7-14*n)/512;
  C3x(13+1) = 7/512;
  C3x(14+1) = 21/2560;
end

function C3 = C3f(epsi, C3x)
%C3F  Evaluate C_3
%
%  C3 = C3F(EPSI, C3X) evaluates C_{3,l} using Eq. (25) and the coefficient
%  vector C3X.  EPSI is a K x 1 arraya.  C3X is a 1 x 15 array.  C3 is a
%  K x 5 array.

  nC3 = 6;
  nC3x = size(C3x, 2);
  j = nC3x;
  C3 = zeros(size(epsi, 1), nC3 - 1);
  for k = nC3 - 1 : -1 : 1,
    t = C3(:, k);
    for i = nC3 - k : -1 : 1,
      t = epsi .* t + C3x(j);
      j = j - 1;
    end
    C3(:, k) = t;
  end
  mult = ones(size(epsi, 1), 1);
  for k = 1 : nC3 - 1,
    mult = mult .* epsi;
    C3(:, k) = C3(:, k) .* mult;
  end
end
