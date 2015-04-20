function [x, y, zone, isnorth, gam, k] = utmups_fwd(lat, lon, setzone)
%UTMUPS_FWD  Forward UTM/UPS projection
%
%   [x, y, zone, isnorth] = UTMUPS_FWD(lat, lon)
%   [x, y, zone, isnorth, gam, k] = UTMUPS_FWD(lat, lon, setzone)
%
%   performs the forward universal transverse Mercator projection of points
%   (lat,lon) to (x,y) using zone and isnorth.  lat and lon can be scalars
%   or arrays of equal size.  zone should be an integer in [1,60] and
%   isnorth is a logical indicating whether the transformation should use
%   the false northing for the northern (isnorth = true) or southern
%   (isnorth = false) hemisphere.  The inverse projection is given by
%   utmups_inv.
%
%   gam and k give metric properties of the projection at (lat,lon); gam is
%   the meridian convergence at the point and k is the scale.
%
%   lat, lon, gam are in degrees.  The projected coordinates x, y are in
%   meters.  k is dimensionless.
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
  isnorth = lat >= 0;
  zone = StandardZone(lat, lon, setzone);
  Z = nan(size(Z));
  x = Z; y = Z; gam = Z; k = Z;
  utm = zone > 0;
  [x(utm), y(utm), gam(utm), k(utm)] = ...
      utm_fwd(zone(utm), isnorth(utm), lat(utm), lon(utm));
  ups = zone == 0;
  [x(ups), y(ups), gam(ups), k(ups)] = ...
      ups_fwd(isnorth(ups), lat(ups), lon(ups));
  zone(isnan(x)) = -4; isnorth(isnan(x)) = false;
end

function [x, y, gam, k] = utm_fwd(zone, isnorth, lat, lon)
%UTM_FWD  Forward UTM projection
%
%   [x, y] = UTM_FWD(zone, isnorth, lat, lon)
%   [x, y, gam, k] = UTM_FWD(zone, isnorth, lat, lon)

  lon0 = -183 + 6 * floor(zone); lat0 = 0;
  bad = ~(abs(mod(lon - lon0 + 180, 360) - 180) <= 60);
  fe = 5e5; fn = 100e5 * (1-isnorth); k0 = 0.9996;
  [x, y, gam, k] = tranmerc_fwd(lat0, lon0, lat, lon);
  x = x * k0; y = y * k0; k = k * k0;
  bad = bad | ~(abs(x) <= 5e5 & y >= -91e5 & y <= 96e5);
  x = x + fe; y = y + fn;
  x(bad) = nan; y(bad) = nan; gam(bad) = nan; k(bad) = nan;
end

function [x, y, gam, k] = ups_fwd(isnorth, lat, lon)
%UPS_FWD  Forward UPS projection

  fe = 20e5; fn = 20e5; k0 = 0.994;
  [x, y, gam, k] = polarst_fwd(isnorth, lat, lon);
  x = x * k0; y = y * k0; k = k * k0;
  lim = (13 - 5 * isnorth) * 1e5;
  bad = ~(abs(x) <= lim & abs(y) <= lim);
  x = x + fe; y = y + fn;
  x(bad) = nan; y(bad) = nan; gam(bad) = nan; k(bad) = nan;
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
  zone(~(zone >= INVALID & zone <= MAXZONE)) = INVALID;
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
