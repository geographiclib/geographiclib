function d = AngDiff(x, y)
%ANGDIFF  Compute angle difference accurately
%
%   d = ANGDIFF(x, y) computes y - x, reduces the result to (-180,180] and
%   rounds the result.  x and y must be in [-180,180].  x and y can be any
%   compatible shapes.

  [d, t] = sumx(-x, y);
  c = (d - 180) + t > 0;
  d(c) = (d(c) - 360) + t(c);
  c = (d + 180) + t <= 0;
  d(c) = (d(c) + 360) + t(c);
end
