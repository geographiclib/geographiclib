function [lat, lon, gam, k] = utmups_inv(x, y, zone, isnorth)
%UTMUPS_INV  Forward UTM/UPS projection
%
%   [lat, lon] = UTMUPS_INV(x, y, zone, isnorth)
%   [lat, lon, gam, k] = UTMUPS_INV(x, y, zone, isnorth)
%
%   performs the inverse universal transverse Mercator projection of points
%   (x,y) to (lat,lon) using zone and isnorth.  x and y can be scalars or
%   arrays of equal size.  zone should be an integer in [1,60] and isnorth
%   is a logical indicating whether the transformation should use the false
%   northing for the northern (isnorth = true) or southern (isnorth =
%   false) hemisphere.  The forward projection is given by utmups_fwd.
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
%   See also UTMUPS_FWD, UTM_INV, UPS_INV

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 4, error('Too few input arguments'), end
  try
    Z = x + y + zone + isnorth;
    Z = zeros(size(Z));
  catch err
    error('x, y, zone, isnorth have incompatible sizes')
  end
  x = x + Z; y = y + Z;
  zone = floor(zone) + Z; isnorth = logical(isnorth + Z);
  Z = nan(size(Z));
  lat = Z; lon = Z; gam = Z; k = Z;
  utm = zone > 0 & zone <= 60;
  [lat(utm), lon(utm), gam(utm), k(utm)] = ...
      utm_inv(zone(utm), isnorth(utm), x(utm), y(utm));
  ups = zone == 0;
  [lat(ups), lon(ups), gam(ups), k(ups)] = ...
      ups_inv(isnorth(ups), x(ups), y(ups));
end

function [lat, lon, gam, k] = utm_inv(zone, isnorth, x, y)
%UTM_INV  Forward UTM projection
%
%   [lat, lon] = UTM_INV(zone, isnorth, x, y)
%   [lat, lon, gam, k] = UTM_INV(zone, isnorth, x, y)

  lon0 = -183 + 6 * floor(zone); lat0 = 0;
  fe = 5e5; fn = 100e5 * (1-isnorth); k0 = 0.9996;
  x = x - fe; y = y - fn;
  bad = ~(abs(x) <= 5e5 & y >= -91e5 & y <= 96e5);
  x = x / k0; y = y / k0;
  [lat, lon, gam, k] = tranmerc_inv(lat0, lon0, x, y);
  k = k * k0;
  lat(bad) = nan; lon(bad) = nan; gam(bad) = nan; k(bad) = nan;
end

function [lat, lon, gam, k] = ups_inv(isnorth, x, y)
%UPS_INV  Inverse UPS projection

  fe = 20e5; fn = 20e5; k0 = 0.994;
  x = x - fe; y = y - fn;
  lim = (13 - 5 * isnorth) * 1e5;
  bad = ~(abs(x) <= lim & abs(y) <= lim);
  x = x / k0; y = y / k0;
  [lat, lon, gam, k] = polarst_inv(isnorth, x, y);
  k = k * k0;
  lat(bad) = nan; lon(bad) = nan; gam(bad) = nan; k(bad) = nan;
end
