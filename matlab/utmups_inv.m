function [lat, lon, gam, k] = utmups_inv(x, y, zone, northp)
%UTM_INV  Forward UTMUPS projection
%
%   [LAT, LON] = UTMUPS_INV(X, Y, ZONE, NORTHP)
%   [LAT, LON, GAM, K] = UTMUPS_INV(X, Y, ZONE, NORTHP)
%
%   See also UTMUPS_FWD, UTM_INV, UPS_INV

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 4, error('Too few input arguments'), end
  try
    Z = x + y + zone + northp;
    Z = zeros(size(Z));
  catch err
    error('x, y, zone, northp have incompatible sizes')
  end
  x = x + Z; y = y + Z;
  zone = floor(zone) + Z; northp = logical(northp + Z);
  Z = nan(size(Z));
  lat = Z; lon = Z; gam = Z; k = Z;
  utm = zone > 0 & zone <= 60;
  [lat(utm), lon(utm), gam(utm), k(utm)] = ...
      utm_inv(zone(utm), northp(utm), x(utm), y(utm));
  ups = zone == 0;
  [lat(ups), lon(ups), gam(ups), k(ups)] = ...
      ups_inv(northp(ups), x(ups), y(ups));
end

function [lat, lon, gam, k] = ups_inv(northp, x, y)
%UPS_INV  Inverse UPS projection

  fe = 20e5; fn = 20e5; k0 = 0.994;
  x = (x - fe) / k0; y = (y - fn) / k0;
  [lat, lon, gam, k] = polarst_inv(northp, x, y);
  k = k * k0;
end
