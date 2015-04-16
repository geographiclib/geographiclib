function [x, y] = normalize(x, y)
%NORMALIZE  Normalize sinx and cosx
%
%   [x, y] = NORMALIZE(x, y) normalize X and Y so that X^2 + Y^2 = 1.  X and Y
%   can be any shape.

  r = hypot(x, y);
  x = x ./ r;
  y = y ./ r;
end
