function ellipsoid = defaultellipsoid
%DEFAULTELLIPSOID  Return the WGS84 ellipsoid
%
%   ELLIPSOID = DEFAULTELLIPSOID
%
%   returns a vector of the equatorial radius and eccentricity for the
%   WGS84 ellipsoid.

  a = 6378137;
  f = 1/298.257223563;
  e = flat2ecc(f);
  ellipsoid = [a, e];
end
