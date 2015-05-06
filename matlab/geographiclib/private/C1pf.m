function C1p = C1pf(epsi)
%C1PF  Evaluate C'_{1,k}
%
%   C1p = C1PF(epsi) evaluates C'_{1,l} using Eq. (21).  epsi is an
%   K x 1 array and C1 is a K x 6 array.

  coeff = [ ...
      205, -432, 768, 1536,     ... C1p[1]/eps^1, polynomial in eps2 of order 2
      4005, -4736, 3840, 12288, ... C1p[2]/eps^2, polynomial in eps2 of order 2
      -225, 116, 384,           ... C1p[3]/eps^3, polynomial in eps2 of order 1
      -7173, 2695, 7680,        ... C1p[4]/eps^4, polynomial in eps2 of order 1
      3467, 7680,               ... C1p[5]/eps^5, polynomial in eps2 of order 0
      38081, 61440,             ... C1p[6]/eps^6, polynomial in eps2 of order 0
          ];
  nC1p = 6;
  C1p = zeros(length(epsi), nC1p);
  eps2 = epsi.^2;
  d = epsi;
  o = 1;
  for  l = 1 : nC1p
    m = floor((nC1p - l) / 2);
    C1p(:,l) = d .* polyval(coeff(o : o + m), eps2) / coeff(o + m + 1);
    o = o + m + 2;
    d = d .* epsi;
  end
end
