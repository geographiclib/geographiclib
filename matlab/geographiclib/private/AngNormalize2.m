function x = AngNormalize2(x)
%ANGNORMALIZE2  Reduce any angle to range [-180, 180)
%
%   x = ANGNORMALIZE(x) reduces arbitrary angles to the range [-180, 180).
%   x can be any shape.

  x = AngNormalize(rem(x, 360));
end
