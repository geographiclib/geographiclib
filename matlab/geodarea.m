function [A, P, N] = geodarea(lat, lon, ellipsoid)
%GEODAREA  Surface area of polygon on an ellipsoid
%
%   A = GEODAREA(lat, lon)
%   [A, P, N] = GEODAREA(lat, lon, ellipsoid)
%
%   calculates the surface area A of the geodesic polygon specified by the
%   input vectors lat, lon (in degrees).  The ellipsoid vector is of the
%   form [a, e], where a is the equatorial radius in meters, e is the
%   eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid is used.
%   There is no need to "close" the polygon by repeating the first point.
%   Multiple polygons can be specified by separating the vertices by NaNs
%   in the vectors.  Thus a series of quadrilaterals can be specified as
%   two 5 x K arrays where the 5th row is NaN.  The output, A, is in
%   meters^2.  Counter-clockwise traversal counts as a positive area.  Only
%   simple polygons (which do not intersect themselves) are supported.
%   Also returned are the perimeters of the polygons in P (meters) and the
%   numbers of vertices in N.  GEODDOC gives the restrictions on the
%   allowed ranges of the arguments.
%
%   GEODAREA loosely duplicates the functionality of the AREAINT function
%   in the MATLAB mapping toolbox.  The major difference is that the
%   polygon edges are taken to be geodesics and the area contributed by
%   each edge is computed using a series expansion with is accurate
%   regardless of the length of the edge.  The formulas are derived in
%
%     C. F. F. Karney
%     Algorithms for geodesics
%     J. Geodesy (2012)
%     http://dx.doi.org/10.1007/s00190-012-0578-z
%
%   See also GEODDOC, GEODDISTANCE, GEODRECKON, POLYGONAREA

% Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
% under the MIT/X11 License.  For more information, see
% http://geographiclib.sourceforge.net/
%
% This file was distributed with GeographicLib 1.27.

  if ~isequal(size(lat), size(lon))
    error('lat, lon have incompatible sizes')
  end

  lat1 = lat(:);
  lon1 = lon(:);
  M = length(lat1);
  ind = [0; find(isnan(lat1 + lon1))];
  if length(ind) == 1 || ind(end) ~= M
    ind = [ind; M + 1];
  end
  K = length(ind) - 1;
  A = zeros(K, 1); P = A; N = A; crossings = A;
  if M == 0, return, end

  lat2 = [lat1(2:end, 1); 0];
  lon2 = [lon1(2:end, 1); 0];
  m0 = min(M, ind(1:end-1) + 1);
  m1 = max(1, ind(2:end) - 1);
  lat2(m1) = lat1(m0); lon2(m1) = lon1(m0);

  if nargin < 3, ellipsoid = defaultellipsoid; end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  f = e2 / (1 + sqrt(1 - e2));

  b = (1 - f) * a;
  if e2 == 0
    c2 = a^2;
  elseif e2 > 0
    c2 = (a^2 + b^2 * atanh(sqrt(e2))/sqrt(e2)) / 2;
  else
    c2 = (a^2 + b^2 * atan(sqrt(-e2))/sqrt(-e2)) / 2;
  end
  area0 = 4 * pi * c2;

  [s12, ~, ~, S12] = geoddistance(lat1, lon1, lat2, lon2, ellipsoid);
  cross = transit(lon1, lon2);

  for k = 1 : K
    N(k) = m1(k) - m0(k) + 1;
    P(k) = sum(s12(m0(k):m1(k)));
    A(k) = -sum(S12(m0(k):m1(k)));
    crossings(k) = sum(cross(m0(k):m1(k)));
  end
  p = mod(crossings, 2) ~= 0;
  A(p) = A(p) - ((A(p) < 0) * 2 - 1) * area0 / 2;
  p = A > area0 / 2;
  A(p) = A(p) - area0;
  p = ~p & A <= -area0 / 2;
  A(p) = A(p) + area0;
end

function x = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range [-180, 180)
%
%   X = ANGNORMALIZE(X) reduces angles in [-540, 540) to the range
%   [-180, 180).  X can be any shape.

  x(x >= 180) = x(x >= 180) - 360;
  x(x < -180) = x(x < -180) + 360;
end

function cross = transit(lon1, lon2)
%TRANSIT  Count crossings of prime meridian
%
%   CROSS = TRANSIT(LON1, LON2) return 1 or -1 if crossing prime meridian
%   in east or west direction.  Otherwise return zero.

  lon1 = AngNormalize(lon1);
  lon2 = AngNormalize(lon2);
  lon12 = -AngNormalize(lon1 - lon2);
  cross = zeros(length(lon1), 1);
  cross(lon1 < 0 & lon2 >= 0 & lon12 > 0) = 1;
  cross(lon2 < 0 & lon1 >= 0 & lon12 < 0) = -1;

end
