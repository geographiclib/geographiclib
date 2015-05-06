function C2 = C2f(epsi)
%C2F  Evaluate C_{2,k}
%
%   C2 = C2F(epsi) evaluates C_{2,l} using Eq. (43).  epsi is an
%   K x 1 array and C2 is a K x 6 array.

  coeff = [ ...
      1, 2, 16, 32,      ... C2[1]/eps^1, polynomial in eps2 of order 2
      35, 64, 384, 2048, ... C2[2]/eps^2, polynomial in eps2 of order 2
      15, 80, 768,       ... C2[3]/eps^3, polynomial in eps2 of order 1
      7, 35, 512,        ... C2[4]/eps^4, polynomial in eps2 of order 1
      63, 1280,          ... C2[5]/eps^5, polynomial in eps2 of order 0
      77, 2048,          ... C2[6]/eps^6, polynomial in eps2 of order 0
          ];
  nC2 = 6;
  C2 = zeros(length(epsi), nC2);
  eps2 = epsi.^2;
  d = epsi;
  o = 1;
  for l = 1 : nC2
    m = floor((nC2 - l) / 2);
    C2(:, l) = d .* polyval(coeff(o : o + m), eps2) / coeff(o + m + 1);
    o = o + m + 2;
    d = d .* epsi;
  end
end
