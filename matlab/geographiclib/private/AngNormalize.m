function x = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range [-180, 180)
%
%   x = ANGNORMALIZE(x) reduces angles in [-540, 540) to the range
%   [-180, 180).  x can be any shape.

  x(x >= 180) = x(x >= 180) - 360;
  x(x < -180) = x(x < -180) + 360;
end
