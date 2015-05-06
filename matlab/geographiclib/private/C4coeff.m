function C4x = C4coeff(n)
%C4COEFF  Evaluate coefficients for C_4
%
%   C4x = C4COEFF(n) evaluates the coefficients of epsilon^l in expansion
%   of the area (Eq. (65) expressed in terms of n and epsi).  n is a
%   scalar.  C4x is a 1 x 21 array.

  coeff = [ ...
      97, 15015,          ... C4[0], coeff of eps^5, polynomial in n of order 0
      1088, 156, 45045,   ... C4[0], coeff of eps^4, polynomial in n of order 1
                          ... C4[0], coeff of eps^3, polynomial in n of order 2
      -224, -4784, 1573, 45045, ...
                          ... C4[0], coeff of eps^2, polynomial in n of order 3
      -10656, 14144, -4576, -858, 45045, ...
                          ... C4[0], coeff of eps^1, polynomial in n of order 4
      64, 624, -4576, 6864, -3003, 15015, ...
                          ... C4[0], coeff of eps^0, polynomial in n of order 5
      100, 208, 572, 3432, -12012, 30030, 45045, ...
      1, 9009,            ... C4[1], coeff of eps^5, polynomial in n of order 0
      -2944, 468, 135135, ... C4[1], coeff of eps^4, polynomial in n of order 1
                          ... C4[1], coeff of eps^3, polynomial in n of order 2
      5792, 1040, -1287, 135135, ...
                          ... C4[1], coeff of eps^2, polynomial in n of order 3
      5952, -11648, 9152, -2574, 135135, ...
                          ... C4[1], coeff of eps^1, polynomial in n of order 4
      -64, -624, 4576, -6864, 3003, 135135, ...
      8, 10725,           ... C4[2], coeff of eps^5, polynomial in n of order 0
      1856, -936, 225225, ... C4[2], coeff of eps^4, polynomial in n of order 1
                          ... C4[2], coeff of eps^3, polynomial in n of order 2
      -8448, 4992, -1144, 225225, ...
                          ... C4[2], coeff of eps^2, polynomial in n of order 3
      -1440, 4160, -4576, 1716, 225225, ...
      -136, 63063,        ... C4[3], coeff of eps^5, polynomial in n of order 0
      1024, -208, 105105, ... C4[3], coeff of eps^4, polynomial in n of order 1
                          ... C4[3], coeff of eps^3, polynomial in n of order 2
      3584, -3328, 1144, 315315, ...
      -128, 135135,       ... C4[4], coeff of eps^5, polynomial in n of order 0
      -2560, 832, 405405, ... C4[4], coeff of eps^4, polynomial in n of order 1
      128, 99099,         ... C4[5], coeff of eps^5, polynomial in n of order 0
          ];
  nC4 = 6;
  nC4x = (nC4 * (nC4 + 1)) / 2;
  C4x = zeros(1, nC4x);
  o = 1;
  k = 1;
  for l = 0 : nC4 - 1
    for j = nC4 - 1 : -1 : l
      m = nC4 - j - 1;
      C4x(k) = polyval(coeff(o : o + m), n) / coeff(o + m + 1);
      k = k + 1;
      o = o + m + 2;
    end
  end
end
