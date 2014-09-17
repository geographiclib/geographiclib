function [lat2, lon2, azi2] = gereckon(lat1, lon1, s12, azi1, ellipsoid)
%GERECKON  Point along great ellipse at specified azimuth and range
%
%   [lat2, lon2, azi2] = GERECKON(lat1, lon1, s12, azi1)
%   [lat2, lon2, azi2] = GERECKON(lat1, lon1, s12, azi1, ellipsoid)
%
%   solves the direct great ellipse problem of finding the final point and
%   azimuth given lat1, lon1, s12, and azi1.  The input arguments lat1,
%   lon1, s12, azi1, can be scalars or arrays of equal size.  lat1, lon1,
%   azi1 are given in degrees and s12 in meters.  The ellipsoid vector is
%   of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by DEFAULTELLIPSOID) is used.  lat2,
%   lon2, and azi2 give the position and forward azimuths at the end point
%   in degrees.
%
%   When given a combination of scalar and array inputs, GERECKON behaves
%   as though the inputs were expanded to match the size of the arrays.
%   However, in the particular case where LAT1 and AZI1 are the same for
%   all the input points, they should be specified as scalars since this
%   will considerably speed up the calculations.  (In particular a series
%   of points along a single geodesic is efficiently computed by specifying
%   an array for S12 only.)
%
%   This is an implementation of the algorithm given in
%
%     https://en.wikipedia.org/wiki/Great_ellipse
%
%   This routine depends on the MATLAB File Exchange package "Geodesics on
%   an ellipsoid of revolution":
%
%     http://www.mathworks.com/matlabcentral/fileexchange/39108
%
%   See also GEDISTANCE, DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2014) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.38.

  if nargin < 4, error('Too few input arguments'), end
  if nargin < 5, ellipsoid = defaultellipsoid; end
  try
    S = size(lat1 + lon1 + s12 + azi1);
  catch err
    error('lat1, lon1, s12, azi1 have incompatible sizes')
  end
  if length(ellipsoid) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  degree = pi/180;
  tiny = sqrt(realmin);

  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  f = e2 / (1 + sqrt(1 - e2));
  f1 = 1 - f;

  lat1 = lat1(:);
  lon1 = AngNormalize(lon1(:));
  azi1 = AngRound(AngNormalize(azi1(:)));
  s12 = s12(:);

  alp1 = azi1 * degree;
  salp1 = sin(alp1); salp1(azi1 == -180) = 0;
  calp1 = cos(alp1); calp1(abs(azi1) == 90) = 0;
  phi = lat1 * degree;
  sbet1 = f1 * sin(phi);
  cbet1 = cos(phi); cbet1(abs(lat1) == 90) = tiny;
  [sbet1, cbet1] = SinCosNorm(sbet1, cbet1);
  [salp1, calp1] = SinCosNorm(salp1 .* sqrt(1 - e2 * cbet1.^2), calp1);
  salp0 = salp1 .* cbet1; calp0 = hypot(calp1, salp1 .* sbet1);
  ssig1 = sbet1; slam1 = salp0 .* sbet1;
  csig1 = cbet1 .* calp1; csig1(sbet1 == 0 & calp1 == 0) = 1; clam1 = csig1;
  [ssig1, csig1] = SinCosNorm(ssig1, csig1);

  k2 = e2 * calp0.^2;
  n = k2 ./ (2 * (1 + sqrt(1 - k2)) - k2);
  A1 = a * (1 + A1m1f(n)) .* (1 - n)./(1 + n);
  C1a = C1f(n);
  B11 = SinCosSeries(true, ssig1, csig1, C1a);
  s = sin(B11); c = cos(B11);
  stau1 = ssig1 .* c + csig1 .* s; ctau1 = csig1 .* c - ssig1 .* s;

  C1pa = C1pf(n);
  tau12 = s12 ./ A1;
  s = sin(tau12); c = cos(tau12);
  B12 = - SinCosSeries(true,  stau1 .* c + ctau1 .* s, ...
                       ctau1 .* c - stau1 .* s, C1pa);
  sig12 = tau12 - (B12 - B11);
  ssig12 = sin(sig12); csig12 = cos(sig12);
  if abs(f) > 0.01
    ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
    csig2 = csig1 .* csig12 - ssig1 .* ssig12;
    B12 =  SinCosSeries(true, ssig2, csig2, C1a);
    serr = (A1/a) .* (sig12 + (B12 - B11)) - s12;
    sig12 = sig12 - serr ./ sqrt(1 - k2 + k2 .* ssig2.^2);
    ssig12 = sin(sig12); csig12 = cos(sig12);
  end

  ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
  csig2 = csig1 .* csig12 - ssig1 .* ssig12;
  sbet2 = calp0 .* ssig2;
  cbet2 = hypot(salp0, calp0 .* csig2);
  cbet2(cbet2 == 0) = tiny;
  slam2 = salp0 .* ssig2; clam2 = csig2;
  salp2 = salp0; calp2 = calp0 .* csig2  .* sqrt(1 - e2 * cbet2.^2);
  lam12 = atan2(slam2 .* clam1 - clam2 .* slam1, ...
                clam2 .* clam1 + slam2 .* slam1);
  lon12 = lam12 / degree;
  lon12 = AngNormalize2(lon12);
  lon2 = AngNormalize(lon1 + lon12);
  lat2 = atan2(sbet2, f1 * cbet2) / degree;
  azi2 = 0 - atan2(-salp2, calp2) / degree;

  lat2 = reshape(lat2, S);
  lon2 = reshape(lon2, S);
  azi2 = reshape(azi2, S);

end
