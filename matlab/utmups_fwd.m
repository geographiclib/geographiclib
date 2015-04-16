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

function [x, y, gam, k] = ups_fwd(northp, lat, lon)
%UPS_FWD  Forward UPS projection

  fe = 20e5; fn = 20e5; k0 = 0.994;
  [x, y, gam, k] = polarst_fwd(northp, lat, lon);
  x = x * k0 + fe; y = y * k0 + fn; k = k * k0;
end
