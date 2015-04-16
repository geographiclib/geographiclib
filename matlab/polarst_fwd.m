function [x, y, gam, k] = polarst_fwd(northp, lat, lon, ellipsoid)
%POLARST_FWD  Forward polar stereographic projection
%
%   [X, Y] = POLARST_FWD(NORTHP, LAT, LON)
%   [X, Y, GAM, K] = POLARST_FWD(NORTHP, LAT, LON, ELLIPSOID)
%
%   See also POLARST_INV.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 3, error('Too few input arguments'), end
  if nargin < 4, ellipsoid = defaultellipsoid; end
  try
    [~] = northp + lat + lon;
  catch err
    error('northp, lat, lon have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  degree = pi/180;
  overflow = 1/eps^2;
  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  e2m = 1 - e2;
  c = sqrt(e2m) * exp(eatanhe(1, e2));

  northp = 2 * logical(northp) - 1;
  lat = lat .* northp;
  phi = lat * degree;
  lam = AngNormalize(lon) .* degree;
  tau  = tan(phi); tau(abs(lat) == 90) = sign(lat(abs(lat) == 90)) * overflow;
  taup = taupf(tau, e2);
  rho = hypot(1, taup) + abs(taup);
  rho(taup >= 0) = cvmgt(1./rho(taup >= 0), 0, lat ~= 90);
  rho = rho * (2 * a / c);
  lon = AngNormalize(lon);
  x = rho .* sin(lam); x(lon == -180) = 0;
  y = -northp .* rho .* cos(lam); x(abs(lon) == 90) = 0;
  if nargout > 2
    gam = northp .* lon;
    if nargout > 3
      secphi = hypot(1, tau);
      k = (rho / a) .* secphi .* sqrt(e2m + e2 .* secphi.^-2);
      k(lat == 90) = 1;
    end
  end
end

