function [lat2, lon2, azi2, S12, m12, M12, M21, a12] = ...
      geodreckon(lat1, lon1, s12, azi1, ellipsoid, arcmode)
%GEODRECKON  Point at specified azimuth, range on ellipsoid
%
%   [lat2, lon2, azi2] = GEODRECKON(lat1, lon1, s12, azi1)
%   [lat2, lon2, azi2, S12, m12, M12, M21, a12] =
%     GEODRECKON(lat1, lon1, s12, azi1, ellipsoid, arcmode)
%
%   solves the direct geodesic problem of finding the final point and
%   azimuth given lat1, lon1, s12, and azi1.  The input arguments lat1,
%   lon1, s12, azi1, can be scalars or arrays of equal size.  lat1, lon1,
%   azi1 are given in degrees and s12 in meters.  The ellipsoid vector is
%   of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid is
%   used.  lat2, lon2, and azi2 give the position and forward azimuths at
%   the end point in degrees.  The other outputs, S12, m12, M12, M21, a12
%   are documented in GEODDOC.
%
%   If arcmode is true (default is false), then the input argument s12 is
%   interpreted as the arc length on the auxiliary sphere (in degrees) and
%   the corresponding distance is returned in the final output variable a12
%   (in meters).  The two optional arguments, ellispoid and arcmode, may be
%   given in any order and either or both may be omitted.
%
%   When given a combination of scalar and array inputs, GEODRECKON behaves
%   as though the inputs were expanded to match the size of the arrays.
%   However, in the particular case where LAT1 and AZI1 are the same for
%   all the input points, they should be specified as scalars since this
%   will considerably speed up the calculations.  (In particular a series
%   of points along a single geodesic is efficiently computed by specifying
%   an array for S12 only.)
%
%   This is an implementation of the algorithm given in
%
%     C. F. F. Karney
%     Algorithms for geodesics
%     J. Geodesy (2012)
%     http://dx.doi.org/10.1007/s00190-012-0578-z
%
%   This function duplicates some of the functionality of the RECKON
%   function in the MATLAB mapping toolbox.  Differences are
%
%     * When the ellipsoid argument is omitted, use the WGS84 ellipsoid.
%     * The azimuth at the end point azi2 is returned.
%     * The solution is accurate to round-off error.
%     * The algorithm is non-iterative and thus may be faster.
%     * Redundant calculations are avoided when computing multiple
%       points on a single geodesic.
%     * Additional properties of the geodesic are calcuated.
%
%   See also GEODDOC, GEODDISTANCE, GEODAREA, GEODESICDIRECT, GEODESICLINE.

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
% with scalar arguments.

  try
    S = size(lat1 + lon1 + s12 + azi1);
  catch err
    error('lat1, lon1, s12, azi1 have incompatible sizes');
  end

  degree = pi/180;
  tiny = sqrt(realmin);

  a = 6378137;
  f = 1/298.257223563;
  e2 = f * (2 - f);
  % Default: ellipsoid = [6378137, 0.081819190842621494335];
  wgs84 = [a, sqrt(e2)];

  if nargin <= 4,
    ellipsoid = wgs84; arcmode = false;
  elseif nargin == 5,
    arg5 = ellipsoid(:);
    if length(arg5) == 2,
      ellipsoid = arg5; arcmode = false;
    else
      arcmode = arg5; ellipsoid = wgs84;
    end
  else
    arg5 = ellipsoid(:);
    arg6 = arcmode;
    if length(arg5) == 2,
      ellipsoid = arg5; arcmode = arg6;
    else
      arcmode = arg5; ellipsoid = arg6;
    end
  end
  if length(ellipsoid) ~= 2,
    error('ellipsoid must be a vector of size 2');
  end
  if ~isscalar(arcmode),
    error('arcmode must be true or false');
  end

  f1 = 1 - f;
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  b = a * f1;

  areap = nargout >= 4;
  redlp = nargout >= 5;
  scalp = nargout >= 6;

  A3x = A3coeff(n);
  C3x = C3coeff(n);

  lat1 = lat1(:);
  lon1 = AngNormalize(lon1(:));
  azi1 = AngNormalize(azi1(:));
  s12 = s12(:);

  p = lat1 == 90;
  lon1(p) = lon1(p) + ((lon1(p) < 1) * 2 - 1) * 180;
  lon1(p) = AngNormalize(lon1(p) - azi1(p));
  azi1(p) = -180;

  p = lat1 == -90;
  lon1(p) = AngNormalize(lon1(p) + azi1(p));
  azi1(p) = 0;

  azi1 = AngRound(azi1);

  alp1 = azi1 * degree;
  salp1 = sin(alp1); salp1(azi1 == -180) = 0;
  calp1 = cos(alp1); calp1(abs(azi1) == 90) = 0;
  phi = lat1 * degree;
  sbet1 = f1 * sin(phi);
  cbet1 = cos(phi); cbet1(abs(lat1) == 90) = tiny;
  [sbet1, cbet1] = SinCosNorm(sbet1, cbet1);
  dn1 = sqrt(1 + ep2 * sbet1.^2);

  salp0 = salp1 .* cbet1; calp0 = hypot(calp1, salp1 .* sbet1);
  ssig1 = sbet1; somg1 = salp0 .* sbet1;
  csig1 = cbet1 .* calp1; csig1(sbet1 == 0 & calp1 == 0) = 1; comg1 = csig1;
  [ssig1, csig1] = SinCosNorm(ssig1, csig1);

  k2 = calp0.^2 * ep2;
  epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
  A1m1 = A1m1f(epsi);
  C1a = C1f(epsi);
  B11 = SinCosSeries(true, ssig1, csig1, C1a);
  s = sin(B11); c = cos(B11);
  stau1 = ssig1 .* c + csig1 .* s; ctau1 = csig1 .* c - ssig1 .* s;

  C1pa = C1pf(epsi);
  C3a = C3f(epsi, C3x);
  A3c = -f * salp0 .* A3f(epsi, A3x);
  B31 = SinCosSeries(true, ssig1, csig1, C3a);

  if arcmode,
    sig12 = s12 * degree;
    ssig12 = sin(sig12);
    csig12 = cos(sig12);
    s12a = abs(s12);
    s12a = s12a - 180 * floor(s12a / 180);
    ssig12(s12a == 0) = 0;
    csig12(s12a == 90) = 0;
  else
    tau12 = s12 ./ (b * (1 + A1m1));
    s = sin(tau12); c = cos(tau12);
    B12 = - SinCosSeries(true,  stau1 .* c + ctau1 .* s, ...
                         ctau1 .* c - stau1 .* s, C1pa);
    sig12 = tau12 - (B12 - B11);
    ssig12 = sin(sig12); csig12 = cos(sig12);
  end

  ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
  csig2 = csig1 .* csig12 - ssig1 .* ssig12;
  dn2 = sqrt(1 + k2 .* ssig2.^2);
  if arcmode | redlp | scalp,
    if arcmode,
      B12 = SinCosSeries(true, ssig2, csig2, C1a);
    end
    AB1 = (1 + A1m1) * (B12 - B11);
  end
  sbet2 = calp0 .* ssig2;
  cbet2 = hypot(salp0, calp0 .* csig2);
  cbet2(cbet2 == 0) = tiny;
  somg2 = salp0 .* ssig2; comg2 = csig2;
  salp2 = salp0; calp2 = calp0 .* csig2;
  omg12 = atan2(somg2 .* comg1 - comg2 .* somg1, ...
                comg2 .* comg1 + somg2 .* somg1);
  lam12 = omg12 + ...
          A3c .* ( sig12 + (SinCosSeries(true, ssig2, csig2, C3a) - B31));
  lon12 = lam12 / degree;
  lon12 = AngNormalize2(lon12);
  lon2 = AngNormalize(lon1 + lon12);
  lat2 = atan2(sbet2, f1 * cbet2) / degree;
  azi2 = 0 - atan2(-salp2, calp2) / degree;
  if arcmode,
    a12 = b * ((1 + A1m1) .* sig12 + AB1);
  else
    a12 = sig12 / degree;
  end

  if redlp | scalp,
    A2m1 = A2m1f(epsi);
    C2a = C2f(epsi);
    B21 = SinCosSeries(true, ssig1, csig1, C2a);
    B22 = SinCosSeries(true, ssig2, csig2, C2a);
    AB2 = (1 + A2m1) * (B22 - B21);
    J12 = (A1m1 - A2m1) .* sig12 + (AB1 - AB2);
    if redlp,
      m12 = b * ((dn2 .* (csig1 .* ssig2) - dn1 .* (ssig1 .* csig2)) ...
                 - csig1 .* csig2 .* J12);
      m12 = reshape(m12, S);
    end
    if scalp,
      t = k2 .* (ssig2 - ssig1) .* (ssig2 + ssig1) ./ (dn1 + dn2);
      M12 = csig12 + (t .* ssig2 - csig2 .* J12) .* ssig1 ./ dn1;
      M21 = csig12 - (t .* ssig1 - csig1 .* J12) .* ssig2 ./  dn2;
      M12 = reshape(M12, S); M21 = reshape(M21, S);
    end
  end

  if areap,
    C4x = C4coeff(n);
    C4a = C4f(k2, C4x);
    A4 = (a^2 * e2) * calp0 .* salp0;
    B41 = SinCosSeries(false, ssig1, csig1, C4a);
    B42 = SinCosSeries(false, ssig2, csig2, C4a);
    salp12 = calp0 .* salp0 .* ...
             cvmgt(csig1 .* (1 - csig12) + ssig12 .* ssig1, ...
                   ssig12 .* (csig1 .* ssig12 ./ (1 + csig12) + ssig1), ...
                   csig12 <= 0);
    calp12 = salp0.^2 + calp0.^2 .* csig1 .* csig2;
    s = calp0 == 0 | salp0 == 0;
    salp12(s) = salp2(s) .* calp1(s) - calp2(s) .* salp1(s);
    calp12(s) = calp2(s) .* calp1(s) + salp2(s) .* salp1(s);
    s = s & salp12 == 0 & calp12 < 0;
    salp12(s) = tiny * calp1; cal12(s) = -1;
    if e2 == 0,
      c2 = a^2;
    elseif e2 > 0,
      c2 = (a^2 + b^2 * atanh(sqrt(e2))/sqrt(e2)) / 2;
    else
      c2 = (a^2 + b^2 * atan(sqrt(-e2))/sqrt(-e2)) / 2;
    end
    S12 = c2 * atan2(salp12, calp12) + A4 .* (B42 - B41);
  end

  lat2 = reshape(lat2, S);
  lon2 = reshape(lon2, S);
  azi2 = reshape(azi2, S);

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

function [sinx, cosx] = SinCosNorm(sinx, cosx)
%SINCOSNORM  Normalize sinx and cosx
%
%   [SINX, COSX] = SINCOSNORM(SINX, COSX) normalize SINX and COSX so that
%   SINX^2 + COSX^2 = 1.  SINX and COSX can be any shape.

  r = hypot(sinx, cosx);
  sinx = sinx ./ r;
  cosx = cosx ./ r;
end

function y = SinCosSeries(sinp, sinx, cosx, c)
%SINSCOSERIES  Evaluate a sine or cosine series using Clenshaw summation
%
%   Y = SINCOSSERIES(SINP, SINX, COSX, C) evaluate
%     y = sum(c[i] * sin( 2*i    * x), i, 1, n), if  sinp
%     y = sum(c[i] * cos((2*i-1) * x), i, 1, n), if ~sinp
%
%   where n is the size of C.  x is given via its sine and cosine in SINX
%   and COSX.  SINP is a scalar.  SINX, COSX, and Y are K x 1 arrays.  C is
%   a K x N array.

  if length(sinx) == 0,
    y = [];
    return;
  end
  n = size(c, 2);
  ar = 2 * (cosx - sinx) .* (cosx + sinx);
  y1 = zeros(length(sinx), 1);
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
  if sinp,
    y = 2 * sinx .* cosx .* y0;
  else
    y = cosx .* (y0 - y1);
  end
end

function x = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range [-180, 180)
%
%   X = ANGNORMALIZE(X) reduces angles in [-540, 540) to the range
%   [-180, 180).  X can be any shape.

  x(x >= 180) = x(x >= 180) - 360;
  x(x < -180) = x(x < -180) + 360;
end

function x = AngNormalize2(x)
%ANGNORMALIZE2  Reduce any angle to range [-180, 180)
%
%   X = ANGNORMALIZE(X) reduces arbitrary angles to the range [-180, 180).
%   X can be any shape.

  x = AngNormalize(mod(x, 360));
end

function y = AngRound(x)
%ANGROUND  Round tiny values so that tiny values become zero.
%
%   Y = ANGROUND(X) rounds X by adding and subtracting 1/16 to it if it is
%   small.  X and Y can be any shape.

  z = 1/16;
  y = abs(x);
  y(y < z) = z - (z - y(y < z));
  y(x < 0) = -y(x < 0);
end

%% SERIES FOR THE GEODESIC PROBLEM

function A1m1 = A1m1f(epsi)
%A1M1F  Evaluate A_1 - 1
%
%   A1M1 = A1M1F(EPSI) evaluates A_1 - 1 using Eq. (17).  EPSI and A1M1 are
%   K x 1 arrays.

  eps2 = epsi.^2;
  t = eps2.*(eps2.*(eps2+4)+64)/256;
  A1m1 = (t + epsi) ./ (1 - epsi);
end

function C1 = C1f(epsi)
%C1F  Evaluate C_{1,k}
%
%   C1 = C1F(EPSI) evaluates C_{1,l} using Eq. (18).  EPSI is a K x 1
%   array and C1 is a K x 6 array.

  nC1 = 6;
  C1 = zeros(length(epsi), nC1);
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

function c = C1pf(epsi)
%C1PF  Evaluate C'_{1,k}
%
%   C1P = C1PF(EPSI) evaluates C'_{1,l} using Eq. (21).  EPSI is an
%   K x 1 array and C1 is a K x 6 array.

  nC1p = 6;
  c = zeros(length(epsi), nC1p);
  eps2 = epsi.^2;
  d = epsi;
  c(:,1) = d.*(eps2.*(205*eps2-432)+768)/1536;
  d = d.*epsi;
  c(:,2) = d.*(eps2.*(4005*eps2-4736)+3840)/12288;
  d = d.*epsi;
  c(:,3) = d.*(116-225*eps2)/384;
  d = d.*epsi;
  c(:,4) = d.*(2695-7173*eps2)/7680;
  d = d.*epsi;
  c(:,5) = 3467*d/7680;
  d = d.*epsi;
  c(:,6) = 38081*d/61440;
end

function A2m1 = A2m1f(epsi)
%A2M1F  Evaluate A_2 - 1
%
%   A2M1 = A2M1F(EPSI) evaluates A_2 - 1 using Eq. (42).  EPSI and A2M1 are
%   K x 1 arrays.

  eps2 = epsi.^2;
  t = eps2.*(eps2.*(25*eps2+36)+64)/256;
  A2m1 = t .* (1 - epsi) - epsi;
end

function C2 = C2f(epsi)
%C2F  Evaluate C_{2,k}
%
%   C2 = C2F(EPSI) evaluates C_{2,l} using Eq. (43).  EPSI is an
%   K x 1 array and C2 is a K x 6 array.

  nC2 = 6;
  C2 = zeros(length(epsi), nC2);
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
%   A3x = A3COEFF(N) evaluates the coefficients of epsilon^l in Eq. (24).  N
%   is a scalar.  A3x is a 1 x 6 array.

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
%   A3 = A3F(EPSI, A3X) evaluates A_3 using Eq. (24) and the coefficient
%   vector A3X.  EPSI and A3 are K x 1 arrays.  A3X is a 1 x 6 array.

  nA3 = 6;
  A3 = zeros(length(epsi), 1);
  for i = nA3 : -1 : 1,
    A3 = epsi .* A3 + A3x(i);
  end
end

function C3x = C3coeff(n)
%C3COEFF  Evaluate coefficients for C_3
%
%   C3x = C3COEFF(N) evaluates the coefficients of epsilon^l in Eq. (25).
%   N is a scalar.  C3x is a 1 x 15 array.

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
%   C3 = C3F(EPSI, C3X) evaluates C_{3,l} using Eq. (25) and the
%   coefficient vector C3X.  EPSI is a K x 1 array.  C3X is a 1 x 15 array.
%   C3 is a K x 5 array.

  nC3 = 6;
  nC3x = size(C3x, 2);
  j = nC3x;
  C3 = zeros(length(epsi), nC3 - 1);
  for k = nC3 - 1 : -1 : 1,
    t = C3(:, k);
    for i = nC3 - k : -1 : 1,
      t = epsi .* t + C3x(j);
      j = j - 1;
    end
    C3(:, k) = t;
  end
  mult = ones(length(epsi), 1);
  for k = 1 : nC3 - 1,
    mult = mult .* epsi;
    C3(:, k) = C3(:, k) .* mult;
  end
end

function C4x = C4coeff(n)
%C4COEFF  Evaluate coefficients for C_4
%
%   C4x = C4COEFF(N) evaluates the coefficients of epsilon^l in expansion
%   of the area (Eq. (65) expressed in terms of n and epsi).  N is a
%   scalar.  C4x is a 1 x 21 array.

  nC4 = 6;
  nC4x = (nC4 * (nC4 + 1)) / 2;
  C4x = zeros(1, nC4x);
  C4x(0+1) = (n*(n*(n*(n*(100*n+208)+572)+3432)-12012)+30030)/45045;
  C4x(1+1) = (n*(n*(n*(64*n+624)-4576)+6864)-3003)/15015;
  C4x(2+1) = (n*((14144-10656*n)*n-4576)-858)/45045;
  C4x(3+1) = ((-224*n-4784)*n+1573)/45045;
  C4x(4+1) = (1088*n+156)/45045;
  C4x(5+1) = 97/15015;
  C4x(6+1) = (n*(n*((-64*n-624)*n+4576)-6864)+3003)/135135;
  C4x(7+1) = (n*(n*(5952*n-11648)+9152)-2574)/135135;
  C4x(8+1) = (n*(5792*n+1040)-1287)/135135;
  C4x(9+1) = (468-2944*n)/135135;
  C4x(10+1) = 1/9009;
  C4x(11+1) = (n*((4160-1440*n)*n-4576)+1716)/225225;
  C4x(12+1) = ((4992-8448*n)*n-1144)/225225;
  C4x(13+1) = (1856*n-936)/225225;
  C4x(14+1) = 8/10725;
  C4x(15+1) = (n*(3584*n-3328)+1144)/315315;
  C4x(16+1) = (1024*n-208)/105105;
  C4x(17+1) = -136/63063;
  C4x(18+1) = (832-2560*n)/405405;
  C4x(19+1) = -128/135135;
  C4x(20+1) = 128/99099;
end

function C4 = C4f(k2, C4x)
%C4F  Evaluate C_4
%
%   C4 = C4F(K2, C4X) evaluates C_{4,l} in the expansion for the area
%   (Eq. (65) expressed in terms of n and epsi) using the coefficient
%   vector C4X.  K2 is a K x 1 array.  C4X is a 1 x 15 array.  C4 is a K x
%   5 array.

  nC4 = 6;
  nC4x = size(C4x, 2);
  epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
  j = nC4x;
  C4 = zeros(length(epsi), nC4);
  for k = nC4 : -1 : 1,
    t = C4(:, k);
    for i = nC4 - k : -1 : 0,
      t = epsi .* t + C4x(j);
      j = j - 1;
    end
    C4(:, k) = t;
  end
  mult = ones(length(epsi), 1);
  for k = 2 : nC4,
    mult = mult .* epsi;
    C4(:, k) = C4(:, k) .* mult;
  end
end
