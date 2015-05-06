function G4x = G4coeff(n)
%G4COEFF  Evaluate coefficients for C_4 for great ellipse
%
%   G4x = G4COEFF(n) evaluates the coefficients of epsilon^l in expansion
%   of the greate ellipse area (expressed in terms of n and epsi).  n is a
%   scalar.  G4x is a 1 x 21 array.

  coeff = [ ...
      -13200233, 1537536, ... G4[0], coeff of eps^5, polynomial in n of order 0
                          ... G4[0], coeff of eps^4, polynomial in n of order 1
      138833443, 13938873, 5765760, ...
                          ... G4[0], coeff of eps^3, polynomial in n of order 2
      -135037988, -32774196, -4232371, 5765760, ...
                          ... G4[0], coeff of eps^2, polynomial in n of order 3
      6417449, 3013374, 1012583, 172458, 720720, ...
                          ... G4[0], coeff of eps^1, polynomial in n of order 4
      -117944, -110552, -84227, -41184, -9009, 120120, ...
                          ... G4[0], coeff of eps^0, polynomial in n of order 5
      200, 416, 1144, 6864, 21021, 15015, 90090, ...
      2625577, 1537536,   ... G4[1], coeff of eps^5, polynomial in n of order 0
                          ... G4[1], coeff of eps^4, polynomial in n of order 1
      -39452953, -3753828, 8648640, ...
                          ... G4[1], coeff of eps^3, polynomial in n of order 2
      71379996, 16424252, 1987557, 17297280, ...
                          ... G4[1], coeff of eps^2, polynomial in n of order 3
      -5975241, -2676466, -847847, -136422, 4324320, ...
                          ... G4[1], coeff of eps^1, polynomial in n of order 4
      117944, 110552, 84227, 41184, 9009, 1081080, ...
      -5512967, 15375360, ... G4[2], coeff of eps^5, polynomial in n of order 0
                          ... G4[2], coeff of eps^4, polynomial in n of order 1
      2443153, 208182, 2882880, ...
                          ... G4[2], coeff of eps^3, polynomial in n of order 2
      -3634676, -741988, -76219, 5765760, ...
                          ... G4[2], coeff of eps^2, polynomial in n of order 3
      203633, 80106, 20735, 2574, 1441440, ...
      22397, 439296,      ... G4[3], coeff of eps^5, polynomial in n of order 0
                          ... G4[3], coeff of eps^4, polynomial in n of order 1
      -71477, -5317, 768768, ...
                          ... G4[3], coeff of eps^3, polynomial in n of order 2
      48020, 8372, 715, 1153152, ...
      -5453, 1317888,     ... G4[4], coeff of eps^5, polynomial in n of order 0
      1407, 91, 329472,   ... G4[4], coeff of eps^4, polynomial in n of order 1
      21, 146432,         ... G4[5], coeff of eps^5, polynomial in n of order 0
          ];
  nG4 = 6;
  nG4x = (nG4 * (nG4 + 1)) / 2;
  G4x = zeros(1, nG4x);
  o = 1;
  k = 1;
  for l = 0 : nG4 - 1
    for j = nG4 - 1 : -1 : l
      m = nG4 - j - 1;
      G4x(k) = polyval(coeff(o : o + m), n) / coeff(o + m + 1);
      k = k + 1;
      o = o + m + 2;
    end
  end
end
