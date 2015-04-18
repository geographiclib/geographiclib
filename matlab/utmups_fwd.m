function [x, y, zone, northp, gam, k] = utmups_fwd(lat, lon, setzone)
%UTM_FWD  Forward UTMUPS projection
%
%   [X, Y, ZONE, NORTHP] = UTMUPS_FWD(LAT, LON)
%   [X, Y, ZONE, NORTHP, GAM, K] = UTMUPS_FWD(LAT, LON, SETZONE)
%
%   See also UTMUPS_INV, UTM_FWD, UPS_FWD

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 2, error('Too few input arguments'), end
  if nargin < 3, setzone = -1; end
  try
    Z = lat + lon + setzone;
    Z = zeros(size(Z));
  catch err
    error('lat, lon, setzone have incompatible sizes')
  end
  lat = lat + Z; lon = lon + Z;
  northp = lat >= 0;
  zone = StandardZone(lat, lon, setzone);
  Z = nan(size(Z));
  x = Z; y = Z; gam = Z; k = Z;
  utm = zone > 0;
  [x(utm), y(utm), gam(utm), k(utm)] = ...
      utm_fwd(zone(utm), northp(utm), lat(utm), lon(utm));
  ups = zone == 0;
  [x(ups), y(ups), gam(ups), k(ups)] = ...
      ups_fwd(northp(ups), lat(ups), lon(ups));
end

function [x, y, gam, k] = utm_fwd(zone, northp, lat, lon)
%UTM_FWD  Forward UTM projection
%
%   [X, Y] = UTM_FWD(ZONE, NORTHP, LAT, LON)
%   [X, Y, GAM, K] = UTM_FWD(ZONE, NORTHP, LAT, LON)
%
%   performs the forward universal transverse Mercator projection of points
%   (LAT,LON) to (X,Y) using ZONE and NORTHP.  LAT and LON can be scalars
%   or arrays of equal size.  ZONE should be an integer in [1,60] and
%   NORTHP is a logical indicating whether the transformation should use
%   the false northing for the northern (NORTHP = true) or southern (NORTHP
%   = false) hemisphere.  The inverse projection is given by UTM_INV.
%
%   GAM and K give metric properties of the projection at (LAT,LON); GAM is
%   the meridian convergence at the point and K is the scale.
%
%   LAT, LON, GAM are in degrees.  The projected coordinates X, Y are in
%   meters.  K is dimensionless.
%
%   This implementation for the UTM projection is based on the series
%   method described in
%
%     C. F. F. Karney, Transverse Mercator with an accuracy of a few
%     nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011);
%     Addenda: http://geographiclib.sf.net/tm-addenda.html
%
%   This extends the series given by Krueger (1912) to sixth order in the
%   flattening.  This is a substantially better series than that used by
%   the MATLAB mapping toolbox.  In particular the errors in the projection
%   are less than 5 nanometers withing 3900 km of the central meridian (and
%   less than 1 mm within 7600 km of the central meridian).  The mapping
%   can be continued accurately over the poles to the opposite meridian.
%
%   This routine depends on the MATLAB File Exchange package "Geodesics on
%   an ellipsoid of revolution":
%
%     http://www.mathworks.com/matlabcentral/fileexchange/39108
%
%   See also GEODPROJ, UTM_INV, TRANMERC_FWD.

% Copyright (c) Charles Karney (2012) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.29.

  if nargin < 4, error('Too few input arguments'), end
  lon0 = -183 + 6 * floor(zone); lat0 = 0;
  fe = 500e3; fn = cvmgt(0,10000e3,logical(northp)); k0 = 0.9996;
  [x, y, gam, k] = tranmerc_fwd(lat0, lon0, lat, lon);
  x = x * k0 + fe; y = y * k0 + fn; k = k * k0;
end

function [x, y, gam, k] = ups_fwd(northp, lat, lon)
%UPS_FWD  Forward UPS projection

  fe = 20e5; fn = 20e5; k0 = 0.994;
  [x, y, gam, k] = polarst_fwd(northp, lat, lon);
  x = x * k0 + fe; y = y * k0 + fn; k = k * k0;
end

function zone = StandardZone(lat, lon, setzone)
  INVALID = -4;
  UTM = -2;
  STANDARD = -1;
  MINZONE = 0;
  MAXZONE = 60;
  UPS = 0;
  if nargin < 3
    setzone = STANDARD;
  end
  zone = floor(setzone) + zeros(size(lat));
  zone(zone > MAXZONE | zone < INVALID) = INVALID;
  g = zone < MINZONE & zone ~= INVALID;
  c = abs(lat) <= 90 & abs(lon) <= 540;
  zone(g & ~c) = INVALID;
  g = g & c;
  c = zone == UTM | (lat >= -80 & lat < 84);
  u = g & c;
  ilon = mod(floor(lon(u)) + 180, 360) - 180;
  z = floor((ilon + 186) / 6);
  % Norway exception
  exception = z == 31 & floor(lat(u) / 8) == 7 & ilon >= 3;
  z(exception) = 32;
  % Svalbard exception
  exception = lat(u) >= 72 & ilon >= 0 & ilon < 42;
  z(exception) = 2 * floor((ilon(exception) + 183)/12) + 1;
  zone(u) = z;
  zone(g & ~c) = UPS;
end
