function e = flat2ecc(f)
%FLAT2ECC   Convert the flattening of an ellipsoid to its eccentricity
%
%  E = FLAT2ECC(F) returns the eccentricity given the flattening.

  e = sqrt(f .* (2 - f));
end
