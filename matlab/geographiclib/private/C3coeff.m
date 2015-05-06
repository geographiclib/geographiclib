function C3x = C3coeff(n)
%C3COEFF  Evaluate coefficients for C_3
%
%   C3x = C3COEFF(n) evaluates the coefficients of epsilon^l in Eq. (25).
%   n is a scalar.  C3x is a 1 x 15 array.

  coeff = [ ...
      3, 128,        ... C3[1], coeff of eps^5, polynomial in n of order 0
      2, 5, 128,     ... C3[1], coeff of eps^4, polynomial in n of order 1
      -1, 3, 3, 64,  ... C3[1], coeff of eps^3, polynomial in n of order 2
      -1, 0, 1, 8,   ... C3[1], coeff of eps^2, polynomial in n of order 2
      -1, 1, 4,      ... C3[1], coeff of eps^1, polynomial in n of order 1
      5, 256,        ... C3[2], coeff of eps^5, polynomial in n of order 0
      1, 3, 128,     ... C3[2], coeff of eps^4, polynomial in n of order 1
      -3, -2, 3, 64, ... C3[2], coeff of eps^3, polynomial in n of order 2
      1, -3, 2, 32,  ... C3[2], coeff of eps^2, polynomial in n of order 2
      7, 512,        ... C3[3], coeff of eps^5, polynomial in n of order 0
      -10, 9, 384,   ... C3[3], coeff of eps^4, polynomial in n of order 1
      5, -9, 5, 192, ... C3[3], coeff of eps^3, polynomial in n of order 2
      7, 512,        ... C3[4], coeff of eps^5, polynomial in n of order 0
      -14, 7, 512,   ... C3[4], coeff of eps^4, polynomial in n of order 1
      21, 2560,      ... C3[5], coeff of eps^5, polynomial in n of order 0
          ];
  nC3 = 6;
  nC3x = (nC3 * (nC3 - 1)) / 2;
  C3x = zeros(1, nC3x);
  o = 1;
  k = 1;
  for l = 1 : nC3 - 1
    for j = nC3 - 1 : -1 : l
      m = min(nC3 - j - 1, j);
      C3x(k) = polyval(coeff(o : o + m), n) / coeff(o + m + 1);
      k = k + 1;
      o = o + m + 2;
    end
  end
end
