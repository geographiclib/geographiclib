function [x, y, z, M] = loccart_fwd(lat0, lon0, h0, lat, lon, h, ellipsoid)
  if nargin < 6, h = 0; end
  if nargin < 7, ellipsoid = defaultellipsoid; end
  try
    S = size(lat + lon + h);
    num = prod(S);
    Z = zeros(num, 1);
    lat = lat(:) + Z;
    lon = lon(:) + Z;
    h = h(:) + Z;
  catch err
    error('lat, lon, h have incompatible sizes')
  end
  if ~(isscalar(lat0) && isscalar(lon0) && isscalar(h0))
    error('lat0, lon0, h0 must be scalar')
  end 
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  [X0, Y0, Z0, M0] = geocent_fwd(lat0, lon0, h0, ellipsoid);
  [X , Y , Z , M ] = geocent_fwd(lat , lon , h , ellipsoid);
  r = [X-X0, Y-Y0, Z-Z0] * M0;
  x = reshape(r(:, 1), S); y = reshape(r(:, 2), S); z = reshape(r(:, 3), S);
  if nargout > 3
    for i = 1:num
      M(:,:, i) = M0' * M(:,:, i);
    end
    M = reshape(M, [3, 3, S]);
  end
end
