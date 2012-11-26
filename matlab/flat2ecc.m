function e = flat2ecc(f)
%FLAT2ECC   Convert the eccentricity of ellipsoid to the flattening
%
%  E = FLAT2ECC(F) returns the eccentricity given the flattening.

  e = sqrt(f .* (2 - f));
end
