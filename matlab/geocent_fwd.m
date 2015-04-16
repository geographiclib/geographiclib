function [X, Y, Z, M] = geocent_fwd(lat, lon, h, ellipsoid)
%GEOCENT_FWD  Conversion from geographic to geocentric coordinates
%
%   [X, Y, Z] = GEOCENT_FWD(LAT, LON, H)
%   [X, Y, Z, M] = GEOCENT_FWD(LAT, LON, H, ELLIPSOID)
%
%   See also GEOCENT_INV.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 2, error('Too few input arguments'), end
  if nargin < 3, h = 0; end
  if nargin < 4, ellipsoid = defaultellipsoid; end
  try
    z = lat + lon + h;
    z = zeros(size(z));
    lat = lat + z; lon = lon + z; h = h + z;
  catch err
    error('lat, lon, h have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  degree = pi/180;
  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  e2m = 1 - e2;

  lon = AngNormalize(lon);

  phi = lat * degree;
  lam = lon * degree;
  sphi = sin(phi);
  cphi = cos(phi); cphi(abs(lat) == 90) = 0;
  n = a./sqrt(1 - e2 * sphi.^2);
  slam = sin(lam); slam(lon == -180) = 0;
  clam = cos(lam); clam(abs(lon) == 90) = 0;
  Z = (e2m * n + h) .* sphi;
  X = (n + h) .* cphi;
  Y = X .* slam;
  X = X .* clam;

  if nargout > 3
    M = GeoRotation(sphi, cphi, slam, clam);
  end
end
