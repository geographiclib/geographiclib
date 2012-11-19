function [s12, azi1, azi2, S12, m12, M12, M21, a12] = ...
      geoddistance(lat1, lon1, lat2, lon2, ellipsoid)
%GEODDISTANCE  Distance between points on an ellipsoid
%
%   [s12, azi1, azi2] = GEODDISTANCE(lat1, lon1, lat2, lon2)
%   [s12, azi1, azi2, S12, m12, M12, M21, a12] =
%      GEODDISTANCE(lat1, lon1, lat2, lon2, ellipsoid)
%
%   solves the inverse geodesic problem of finding of length and azimuths
%   of the shortest geodesic between points specified by lat1, lon1, lat2,
%   lon2.  The input latitudes and longitudes, lat1, lon1, lat2, lon2, can
%   be scalars or arrays of equal size and must be expressed in degrees.
%   The ellipsoid vector is of the form [a, e], where a is the equatorial
%   radius in meters, e is the eccentricity.  If ellipsoid is omitted, the
%   WGS84 ellipsoid is used.  The output s12 is the distance in meters and
%   azi1 and azi2 are the forward azimuths at the end points in degrees.
%   The other optional outputs, S12, m12, M12, M21, a12 are documented in
%   GEODDOC.  GEODDOC also gives the restrictions on the allowed ranges of
%   the arguments.
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
%   This function duplicates some of the functionality of the DISTANCE
%   function in the MATLAB mapping toolbox.  Differences are
%
%     * When the ellipsoid argument is omitted, use the WGS84 ellipsoid.
%     * The routines work for prolate (as well as oblate) ellipsoids.
%     * The azimuth at the second end point azi2 is returned.
%     * The solution is accurate to round off for abs(e) < 0.2.
%     * The algorithm converges for all pairs of input points.
%     * Additional properties of the geodesic are calcuated.
%
%   See also GEODDOC, GEODRECKON, GEODAREA, GEODESICINVERSE.

% Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
% under the MIT/X11 License.  For more information, see
% http://geographiclib.sourceforge.net/
%
% This file was distributed with GeographicLib 1.27.
%
% This is a straightforward transcription of the C++ implementation in
% GeographicLib and the C++ source should be consulted for additional
% documentation.  This is a vector implementation and the results returned
% with array arguments are identical to those obtained with multiple calls
% with scalar arguments.  The biggest change was to eliminate the branching
% to allow a vectorized solution.

  try
    Z = lat1 + lon1 + lat2 + lon2;
    S = size(Z);
    Z = zeros(S);
    lat1 = lat1 + Z; lon1 = lon1 + Z;
    lat2 = lat2 + Z; lon2 = lon2 + Z;
    Z = Z(:);
  catch err
    error('lat1, lon1, s12, azi1 have incompatible sizes')
  end

  degree = pi/180;
  tiny = sqrt(realmin);
  tol0 = eps;
  tolb = eps * sqrt(eps);
  maxit1 = 20;
  maxit2 = maxit1 + (-log2(eps) + 1) + 10;

  if nargin < 5, ellipsoid = defaultellipsoid; end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  f = e2 / (1 + sqrt(1 - e2));

  f1 = 1 - f;
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  b = a * f1;

  areap = nargout >= 4;
  scalp = nargout >= 6;

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
  s12 = Z; m12 = Z; M12 = Z; M21 = Z; omg12 = Z;

  m = lat1 == -90 | slam12 == 0;

  if any(m)
    calp1(m) = clam12(m); salp1(m) = slam12(m);
    calp2(m) = 1; salp2(m) = 0;

    ssig1(m) = sbet1(m); csig1(m) = calp1(m) .* cbet1(m);
    ssig2(m) = sbet2(m); csig2(m) = calp2(m) .* cbet2(m);

    sig12(m) = atan2(max(csig1(m) .* ssig2(m) - ssig1(m) .* csig2(m), 0), ...
                     csig1(m) .* csig2(m) + ssig1(m) .* ssig2(m));

    [s12(m), m12(m), ~, M12(m), M21(m)] = ...
        Lengths(n, sig12(m), ...
                ssig1(m), csig1(m), dn1(m), ssig2(m), csig2(m), dn2(m), ...
                cbet1(m), cbet2(m), scalp, ep2);
    m = m & (sig12 < 1 | m12 >= 0);
    s12(m) = s12(m) * b;
  end

  eq = ~m & sbet1 == 0;
  if f > 0
    eq = eq & lam12 < pi - f * pi;
  end
  calp1(eq) = 0; calp2(eq) = 0; salp1(eq) = 1; salp2(eq) = 1;
  s12(eq) = a * lam12(eq); sig12(eq) = lam12(eq) / f1; omg12(eq) = sig12(eq);
  m12(eq) = b * sin(omg12(eq)); M12(eq) = cos(omg12(eq)); M21(eq) = M12(eq);

  g = ~eq & ~m;

  [sig12(g), salp1(g), calp1(g), salp2(g), calp2(g)] = ...
      InverseStart(sbet1(g), cbet1(g), dn1(g), sbet2(g), cbet2(g), dn2(g), ...
                   lam12(g), f, A3x);

  s = g & sig12 >= 0;
  dnm = (dn1(s) + dn2(s)) / 2;
  s12(s) = b * sig12(s) .* dnm;
  m12(s) = b * dnm.^2 .* sin(sig12(s) ./ dnm);
  if scalp
    M12(s) = cos(sig12(s) ./ dnm); M21(s) = M12(s);
  end
  omg12(s) = lam12(s) ./ (f1 * dnm);

  g = g & sig12 < 0;

  salp1a = Z + tiny; calp1a = Z + 1;
  salp1b = Z + tiny; calp1b = Z - 1;
  ssig1 = Z; csig1 = Z; ssig2 = Z; csig2 = Z;
  epsi = Z; v = Z; dv = Z;
  numit = Z;
  tripn = Z > 0;
  tripb = tripn;
  gsave = g;
  for k = 0 : maxit2 - 1
    if k == 0 && ~any(g), break, end
    numit(g) = k;
    [v(g), dv(g), ...
     salp2(g), calp2(g), sig12(g), ...
     ssig1(g), csig1(g), ssig2(g), csig2(g), epsi(g), omg12(g)] = ...
        Lambda12(sbet1(g), cbet1(g), dn1(g), ...
                 sbet2(g), cbet2(g), dn2(g), ...
                 salp1(g), calp1(g), f, A3x, C3x);
    v = v - lam12;
    g = g & ~(tripb | abs(v) < ((tripn * 6) + 2) * tol0);
    if ~any(g), break, end

    c = g & v > 0;
    if k <= maxit1
      c = c & calp1 ./ salp1 > calp1b ./ salp1b;
    end
    salp1b(c) = salp1(c); calp1b(c) = calp1(c);

    c = g & v < 0;
    if k <= maxit1
      c = c & calp1 ./ salp1 < calp1a ./ salp1a;
    end
    salp1a(c) = salp1(c); calp1a(c) = calp1(c);

    if k == maxit1, tripn(g) = false; end
    if k < maxit1
      dalp1 = -v ./ dv;
      sdalp1 = sin(dalp1); cdalp1 = cos(dalp1);
      nsalp1 = salp1 .* cdalp1 + calp1 .* sdalp1;
      calp1(g) = calp1(g) .* cdalp1(g) - salp1(g) .* sdalp1(g);
      salp1(g) = nsalp1(g);
      tripn = g & abs(v) <= 16 * tol0;
      c = g & ~(dv > 0 & nsalp1 > 0 & abs(dalp1) < pi);
      tripn(c) = false;
    else
      c = g;
    end

    salp1(c) = (salp1a(c) + salp1b(c))/2;
    calp1(c) = (calp1a(c) + calp1b(c))/2;
    [salp1(g), calp1(g)] = SinCosNorm(salp1(g), calp1(g));
    tripb(c) = (abs(salp1a(c) - salp1(c)) + (calp1a(c) - calp1(c)) < tolb | ...
                abs(salp1(c) - salp1b(c)) + (calp1(c) - calp1b(c)) < tolb);
  end

  g = gsave;
  [s12(g), m12(g), ~, M12(g), M21(g)] = ...
      Lengths(epsi(g), sig12(g), ...
              ssig1(g), csig1(g), dn1(g), ssig2(g), csig2(g), dn2(g), ...
              cbet1(g), cbet2(g), scalp, ep2);

  m12(g) = m12(g) * b;
  s12(g) = s12(g) * b;
  omg12(g) = lam12(g) - omg12(g);

  s12 = 0 + s12;

  if areap
    salp0 = salp1 .* cbet1; calp0 = hypot(calp1, salp1 .* sbet1);
    ssig1 = sbet1; csig1 = calp1 .* cbet1;
    ssig2 = sbet2; csig2 = calp2 .* cbet2;
    k2 = calp0.^2 * ep2;
    epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
    A4 = (a^2 * e2) * calp0 .* salp0;
    [ssig1, csig1] = SinCosNorm(ssig1, csig1);
    [ssig2, csig2] = SinCosNorm(ssig2, csig2);

    C4x = C4coeff(n);
    C4a = C4f(epsi, C4x);
    B41 = SinCosSeries(false, ssig1, csig1, C4a);
    B42 = SinCosSeries(false, ssig2, csig2, C4a);
    S12 = A4 .* (B42 - B41);
    S12(calp0 == 0 | salp0 == 0) = 0;

    l = ~m & omg12 < 0.75 * pi & sbet2 - sbet1 < 1.75;
    alp12 = Z;
    somg12 = sin(omg12(l)); domg12 = 1 + cos(omg12(l));
    dbet1 = 1 + cbet1(l); dbet2 = 1 + cbet2(l);
    alp12(l) = 2 * atan2(somg12 .* (sbet1(l) .* dbet2 + sbet2(l) .* dbet1), ...
                         domg12 .* (sbet1(l) .* sbet2(l) + dbet1 .* dbet2));
    l = ~l;
    salp12 = salp2(l) .* calp1(l) - calp2(l) .* salp1(l);
    calp12 = calp2(l) .* calp1(l) + salp2(l) .* salp1(l);
    s = salp12 == 0 & calp12 < 0;
    salp12(s) = tiny * calp1(s); calp12(s) = -1;
    alp12(l) = atan2(salp12, calp12);
    if e2 == 0
      c2 = a^2;
    elseif e2 > 0
      c2 = (a^2 + b^2 * atanh(sqrt(e2))/sqrt(e2)) / 2;
    else
      c2 = (a^2 + b^2 * atan(sqrt(-e2))/sqrt(-e2)) / 2;
    end
    S12 = 0 + swapp .* lonsign .* latsign .* (S12 + c2 * alp12);
  end

  [salp1(swapp<0), salp2(swapp<0)] = swap(salp1(swapp<0), salp2(swapp<0));
  [calp1(swapp<0), calp2(swapp<0)] = swap(calp1(swapp<0), calp2(swapp<0));
  if scalp
    [M12(swapp<0), M21(swapp<0)] = swap(M12(swapp<0), M21(swapp<0));
  end
  salp1 = salp1 .* swapp .* lonsign; calp1 = calp1 .* swapp .* latsign;
  salp2 = salp2 .* swapp .* lonsign; calp2 = calp2 .* swapp .* latsign;

  azi1 = 0 - atan2(-salp1, calp1) / degree;
  azi2 = 0 - atan2(-salp2, calp2) / degree;
  a12 = sig12 / degree;

  s12 = reshape(s12, S); azi1 = reshape(azi1, S); azi2 = reshape(azi2, S);
  m12 = reshape(m12, S); M12 = reshape(M12, S); M21 = reshape(M21, S);
  a12 = reshape(a12, S);
  if (areap)
    S12 = reshape(S12, S);
  end
end
