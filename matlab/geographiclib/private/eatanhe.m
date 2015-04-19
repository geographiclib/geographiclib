function y = eatanhe(x, e2)
%EATANHE   e*atanh(e*x)
%
%   ATANHEE(X, E2) returns E*atanh(E*X) where E = SQRT(E2)
%   E2 is a scalar; X can be any shape.

  e = sqrt(abs(e2));
  if (e2 >= 0)
    y = e * atanh(e * x);
  else
    y = e * atan(e * x);
  end
end
