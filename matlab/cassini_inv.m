function [lat, lon, azi, rk] = cassini_inv(lat0, lon0, x, y, ellipsoid)
%CASSINI_INV  Inverse Cassini-Soldner projection
%
%   [LAT, LON] = CASSINI_INV(LAT0, LON0, X, Y)
%   [LAT, LON, AZI, RK] = CASSINI_INV(LAT0, LON0, X, Y, ELLIPSOID)
%
%   performs the inverse Cassini-Soldner projection of points (X,Y) using
%   (LAT0,LON0) as the center of projection.  These input arguments can be
%   scalars or arrays of equal size.  The ELLIPSOID vector is of the form
%   [a, e], where a is the equatorial radius in meters, e is the
%   eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by DEFAULTELLIPSOID) is used.  GEODPROJ
%   gives the restrictions on the allowed ranges of the arguments.  The
%   forward projection is given by CASSINI_FWD.
%
%   AZI and RK give metric properties of the projection at (LAT,LON); AZI
%   is the azimuth of the easting (x) direction and RK is the reciprocal of
%   the northing (y) scale.  The scale in the easing direction is 1.
%
%   LAT0, LON0, LAT, LON, AZI are in degrees.  The projected coordinates X,
%   Y are in meters (more precisely the units used for the equatorial
%   radius).  RK is dimensionless.
%
%   The ellipsoidal azimuthal projection is an azimuthal projection about a
%   center point.  The distance and azimuth from the center point to all
%   other points are correctly represented in the projection.  Section
%   14 of
%
%     C. F. F. Karney,
%     Geodesics on an ellipsoid of revolution (2011),
%     http://arxiv.org/abs/1102.1215
%     Errata: http://geographiclib.sf.net/geod-addenda.html#geod-errata
%
%   describes how to use this projection in the determination of maritime
%   boundaries (finding the median line).
%
%   This routine depends on the MATLAB File Exchange package "Geodesics on
%   an ellipsoid of revolution":
%
%     http://www.mathworks.com/matlabcentral/fileexchange/39108
%
%   See also GEODPROJ, CASSINI_FWD, GEODRECKON, DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2012) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.28.

  if nargin < 4, error('Too few input arguments'), end
  if nargin < 5, ellipsoid = defaultellipsoid; end
  try
    [~] = lat0 + lon0 + x + y;
  catch err
    error('lat0, lon0, x, y have incompatible sizes')
  end

  [lat1, lon1, azi0] = geodreckon(lat0, lon0, y, 0, ellipsoid);
  [lat, lon, azi, ~, ~, rk] = ...
      geodreckon(lat1, lon1, x, azi0 + 90, ellipsoid);
end
