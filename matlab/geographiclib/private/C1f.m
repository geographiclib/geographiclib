function C1 = C1f(epsi)
%C1F  Evaluate C_{1,k}
%
%   C1 = C1F(epsi) evaluates C_{1,l} using Eq. (18).  epsi is a K x 1
%   array and C1 is a K x 6 array.

  coeff = [ ...
      -1, 6, -16, 32,     ... C1[1]/eps^1, polynomial in eps2 of order 2
      -9, 64, -128, 2048, ... C1[2]/eps^2, polynomial in eps2 of order 2
      9, -16, 768,        ... C1[3]/eps^3, polynomial in eps2 of order 1
      3, -5, 512,         ... C1[4]/eps^4, polynomial in eps2 of order 1
      -7, 1280,           ... C1[5]/eps^5, polynomial in eps2 of order 0
      -7, 2048,           ... C1[6]/eps^6, polynomial in eps2 of order 0
          ];
  nC1 = 6;
  C1 = zeros(length(epsi), nC1);
  eps2 = epsi.^2;
  d = epsi;
  o = 1;
  for  l = 1 : nC1
    m = floor((nC1 - l) / 2);
    C1(:,l) = d .* polyval(coeff(o : o + m), eps2) / coeff(o + m + 1);
    o = o + m + 2;
    d = d .* epsi;
  end
end
