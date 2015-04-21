function A3 = A3f(epsi, A3x)
%A3F  Evaluate A_3
%
%   A3 = A3F(epsi, A3x) evaluates A_3 using Eq. (24) and the coefficient
%   vector A3x.  epsi and A3 are K x 1 arrays.  A3x is a 1 x 6 array.

  nA3 = 6;
  A3 = zeros(length(epsi), 1);
  for i = nA3 : -1 : 1
    A3 = epsi .* A3 + A3x(i);
  end
end
