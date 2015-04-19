function [lat, lon, gam, k] = polarst_inv(isnorth, x, y, ellipsoid)
%POLARST_INV  Forward polar stereographic projection
%
%   [lat, lon] = POLARST_INV(isnorth, x, y)
%   [lat, lon, gam, k] = POLARST_INV(isnorth, x, y, ellipsoid)
%
%   See also POLARST_FWD.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 3, error('Too few input arguments'), end
  if nargin < 4, ellipsoid = defaultellipsoid; end
  try
    [~] = isnorth + x + y;
  catch err
    error('isnorth, x, y have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  degree = pi/180;
  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  e2m = 1 - e2;
  c = sqrt(e2m) * exp(eatanhe(1, e2));

  isnorth = 2 * logical(isnorth) - 1;
  rho = hypot(x, y);
  t = rho / (2 * a / c);
  taup = (1 ./ t - t) / 2;
  tau = tauf(taup, e2);
  phi = atan(tau);
  lat =  phi / degree;
  lat(rho == 0) = 90;
  lat = isnorth .* lat;
  lon = 0 - atan2( -x, -isnorth .* y) / degree;
  if nargout > 2
    gam = isnorth .* lon;
    if nargout > 3
      secphi = hypot(1, tau);
      k = (rho / a) .* secphi .* sqrt(e2m + e2 .* secphi.^-2);
      k(rho == 0) = 1;
    end
  end
end
