function taup = taupf(tau, e2)
%TAUPF   tan(chi)
%
%   TAUPF(TAU, E2) returns tangent of chi in terms of TAU the tangent of phi.
%   E2, the square of the eccentricity, is a scalar; TAUP can be any shape.

  tau1 = hypot(1, tau);
  sig = sinh( eatanhe( tau ./ tau1, e2 ) );
  taup = hypot(1, sig) .* tau - sig .* tau1;
end
