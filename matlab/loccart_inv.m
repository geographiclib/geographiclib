function [lat, lon, h, M] = loccart_inv(lat0, lon0, h0, x, y, z, ellipsoid)
  if nargin < 7, ellipsoid = defaultellipsoid; end
  try
    S = size(x + y + z);
    num = prod(S);
    Z = zeros(num, 1);
    x = x(:) + Z;
    y = y(:) + Z;
    z = z(:) + Z;
  catch err
    error('x, y, z have incompatible sizes')
  end
  if ~(isscalar(lat0) && isscalar(lon0) && isscalar(h0))
    error('lat0, lon0, h0 must be scalar')
  end 
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  [X0, Y0, Z0, M0] = geocent_fwd(lat0, lon0, h0, ellipsoid);
  r = [x, y, z] * M0';
  X = r(:, 1) + X0; Y = r(:, 2) + Y0; Z = r(:, 3) + Z0;
  [lat , lon , h , M] = geocent_inv(X, Y, Z, ellipsoid);
  lat = reshape(lat, S); lon = reshape(lon, S); h = reshape(h, S);
  if nargout > 3
    for i = 1:num
      M(:,:, i) = M0' * M(:,:, i);
    end
    M = reshape(M, [3, 3, S]);
  end
end
