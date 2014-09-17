function [s12, azi1, azi2] = gedistance(lat1, lon1, lat2, lon2, ellipsoid)
%GEDISTANCE  Great ellipse distance between points on an ellipsoid
%
%   [s12, azi1, azi2] = GEDISTANCE(lat1, lon1, lat2, lon2)
%   [s12, azi1, azi2] = GEDISTANCE(lat1, lon1, lat2, lon2, ellipsoid)
%
%   solves the inverse great ellipse problem of finding of length and
%   azimuths of the great ellipse between points specified by lat1, lon1,
%   lat2, lon2.  The input latitudes and longitudes, lat1, lon1, lat2,
%   lon2, can be scalars or arrays of equal size and must be expressed in
%   degrees.  The ellipsoid vector is of the form [a, e], where a is the
%   equatorial radius in meters, e is the eccentricity.  If ellipsoid is
%   omitted, the WGS84 ellipsoid (more precisely, the value returned by
%   DEFAULTELLIPSOID) is used.  The output s12 is the distance in meters
%   and azi1 and azi2 are the forward azimuths at the end points in
%   degrees.
%
%   When given a combination of scalar and array inputs, the scalar inputs
%   are automatically expanded to match the size of the arrays.
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
%   See also GERECKON, DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2014) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.38.

  if nargin < 4, error('Too few input arguments'), end
  if nargin < 5, ellipsoid = defaultellipsoid; end
  try
    Z = lat1 + lon1 + lat2 + lon2;
    S = size(Z);
    Z = zeros(S);
    lat1 = lat1 + Z; lon1 = lon1 + Z;
    lat2 = lat2 + Z; lon2 = lon2 + Z;
  catch err
    error('lat1, lon1, s12, azi1 have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  degree = pi/180;
  tiny = sqrt(realmin);

  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  f = e2 / (1 + sqrt(1 - e2));

  f1 = 1 - f;

  lon12 = AngDiff(AngNormalize(lon1(:)), AngNormalize(lon2(:)));
  lon12 = AngRound(lon12);

  phi = lat1 * degree;
  sbet1 = f1 * sin(phi); cbet1 = cos(phi); cbet1(lat1 == -90) = tiny;
  [sbet1, cbet1] = SinCosNorm(sbet1, cbet1);

  phi = lat2 * degree;
  sbet2 = f1 * sin(phi); cbet2 = cos(phi); cbet2(abs(lat2) == 90) = tiny;
  [sbet2, cbet2] = SinCosNorm(sbet2, cbet2);

  lam12 = lon12 * degree;
  slam12 = sin(lam12); slam12(lon12 == 180) = 0; clam12 = cos(lam12);

  % Solve great circle
  salp1 = cbet2 .* slam12; calp1 = +cbet1 .* sbet2 - sbet1 .* cbet2 .* clam12;
  salp2 = cbet1 .* slam12; calp2 = -sbet1 .* cbet2 + cbet1 .* sbet2 .* clam12;
  ssig12 = hypot(salp1, calp1);
  csig12 = sbet1 .* sbet2 + cbet1 .* cbet2 .* clam12;
  [salp1, calp1] = SinCosNorm(salp1, calp1);
  [salp2, calp2] = SinCosNorm(salp2, calp2);
  % no need to normalize [ssig12, csig12]

  calp0 = hypot(calp1, salp1 .* sbet1);
  
  ssig1 = sbet1; csig1 = cbet1 .* calp1;
  [ssig1, csig1] = SinCosNorm(ssig1, csig1);
  ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
  csig2 = csig1 .* csig12 - ssig1 .* ssig12;

  k2 = e2 * calp0.^2;
  n = k2 ./ (2 * (1 + sqrt(1 - k2)) - k2);
  C1a = C1f(n);
  A1 = a * (1 + A1m1f(n)) .* (1 - n)./(1 + n);
  s12 = A1 .* (atan2(ssig12, csig12) + ...
               (SinCosSeries(true, ssig2, csig2, C1a) - ...
                SinCosSeries(true, ssig1, csig1, C1a)));
  calp1 = calp1 .* sqrt(1 - e2 * cbet1.^2);
  calp2 = calp2 .* sqrt(1 - e2 * cbet2.^2);
  azi1 = atan2(salp1, calp1) / degree;
  azi2 = atan2(salp2, calp2) / degree;

  s12 = reshape(s12, S); azi1 = reshape(azi1, S); azi2 = reshape(azi2, S);
end
