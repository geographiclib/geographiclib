function [lat, lon, azi, rk] = gnomonicrev(lat0, lon0, x, y, ellipsoid)
%GNOMONICREV  Reverse ellipsoidal gnomonic projection
%
%   [LAT, LON] = GNOMONICREV(LAT0, LON0, X, Y)
%   [LAT, LON, AZI, RK] = GNOMONICREV(LAT0, LON0, X, Y, ELLIPSOID)
%
%   performs the reverse ellipsoidal gnomonic projection of points (X,Y)
%   using (LAT0,LON0) as the center of projection.  These input arguments
%   can be scalars or arrays of equal size.  The ELLIPSOID vector is of the
%   form [a, e], where a is the equatorial radius in meters, e is the
%   eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by DEFAULTELLIPSOID) is used.
%
%   AZI and RK give metric properties of the projection at (LAT,LON); AZI
%   is the azimuth of the geodesic from the center of projection and RK is
%   the reciprocal of the azimuthal scale.  The scale in the radial
%   direction is 1/RK^2.
%
%   In principle, all finite X and Y are allowed.  However, it's possible
%   that the reverse projection fails for very large X and Y (when the
%   geographic position is close to the "horizon").  In that case, NaNs are
%   returned for the corresponding output variables.
%
%   LAT0, LON0, LAT, LON, AZI are in degrees.  The projected coordinates X,
%   Y are in meters (more precisely the units used for the equatorial
%   radius).  RK is dimensionless.
%
%   The ellipsoidal gnomonic projection is an azimuthal projection about a
%   center point.  All geodesics through the center point are projected
%   into straight lines with the correct azimuth relative to the center
%   point.  In addition all geodesics are pass close to the center point
%   are very nearly straight.  The projection is derived in Section 8 of
%
%     C. F. F. Karney,
%     Algorithms for geodesics,
%     J. Geodesy (2012);
%     http://dx.doi.org/10.1007/s00190-012-0578-z
%     Addenda: http://geographiclib.sf.net/geod-addenda.html
%
%   which also includes methods for solving the "intersection" and
%   "interception" problems using the gnomonic projection.
%
%   This routine depends on the MATLAB File Exchange package "Geodesics on
%   an ellipsoid of revolution":
%  
%     http://www.mathworks.com/matlabcentral/fileexchange/39108
%
%   See also GNOMONICFWD, GEODRECKON, DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
% under the MIT/X11 License.  For more information, see
% http://geographiclib.sourceforge.net/
%
% This file was distributed with GeographicLib 1.28.

  try
    Z = lat0 + lon0 + x + y;
    Z = zeros(size(Z));
  catch err
    error('lat0, lon0, x, y have incompatible sizes')
  end
  if nargin < 5, ellipsoid = defaultellipsoid; end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  a = ellipsoid(1);
  numit = 5;
  eps1 = a * 0.01 * sqrt(eps);

  lat0 = lat0 + Z; lon0 = lon0 + Z; x = x + Z; y = y + Z;
  azi0 = atan2(x, y) / (pi/180);
  rho = hypot(x, y);
  s = a * atan(rho / a);
  little = rho <= a;
  rho(~little) = 1 ./ rho(~little);
  g = Z == 0; trip = ~g;
  lat = Z; lon = Z; azi = Z; m = Z; M = Z; ds = Z;
  for k = 1 : numit
    [lat(g), lon(g), azi(g), ~, m(g), M(g)] = ...
        geodreckon(lat0(g), lon0(g), s(g), azi0(g), ellipsoid);
    g = g & ~trip;
    if ~any(g), break, end
    c = little & g;
    ds(c) = (m(c) ./ M(c) - rho(c)) .* M(c).^2;
    c = ~little & g;
    ds(c) = (rho(c) - M(c) ./ m(c)) .* m(c).^2;
    s(g) = s(g) - ds(g);
    trip(g) = abs(ds(g)) < eps1;
  end
  c = ~trip;
  if any(c)
    lat(c) = NaN;
    lon(c) = NaN;
    azi(c) = NaN;
    M(c) = NaN;
  end
  rk = M;
end
