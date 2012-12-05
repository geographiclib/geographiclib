function [x, y, azi, rk] = gnomonicfwd(lat0, lon0, lat, lon, ellipsoid)
%GNOMONICFWD  Forward ellipsoidal gnomonic projection
%
%   [X, Y] = GNOMONICFWD(LAT0, LON0, LAT, LON)
%   [X, Y, AZI, RK] = GNOMONICFWD(LAT0, LON0, LAT, LON, ELLIPSOID)
%
%   performs the forward ellipsoidal gnomonic projection of points
%   (LAT,LON) using (LAT0,LON0) as the center of projection.  These input
%   arguments can be scalars or arrays of equal size.  The ELLIPSOID vector
%   is of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by DEFAULTELLIPSOID) is used.
%
%   AZI and RK give metric properties of the projection at (LAT,LON); AZI
%   is the azimuth of the geodesic from the center of projection and RK is
%   the reciprocal of the azimuthal scale.  The scale in the radial
%   direction is 1/RK^2.
%
%   If the point lies "over the horizon", i.e., if RK <= 0, then NaNs are
%   returned for X and Y (the correct values are returned for AZI and RK).
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
%   See also GNOMONICREV, GEODDISTANCE, DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
% under the MIT/X11 License.  For more information, see
% http://geographiclib.sourceforge.net/
%
% This file was distributed with GeographicLib 1.28.

  try
    Z = lat0 + lon0 + lat + lon;
    Z = zeros(size(Z));
  catch err
    error('lat0, lon0, lat, lon have incompatible sizes')
  end
  if nargin < 5, ellipsoid = defaultellipsoid; end

  [~, azi0, azi, ~, m, M] = geoddistance(lat0, lon0, lat, lon, ellipsoid);
  rho = m ./ M;
  azi0 = azi0 * (pi/180);
  x = rho .* sin(azi0);
  y = rho .* cos(azi0);
  rk = M;
  x(M <= 0) = NaN;
  y(M <= 0) = NaN;
end
