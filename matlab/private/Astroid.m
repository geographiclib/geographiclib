function k = Astroid(x, y)
% ASTROID  Solve the astroid equation
%
%   K = ASTROID(X, Y) solves the quartic polynomial Eq. (55)
%
%     K^4 + 2 * K^3 - (X^2 + Y^2 - 1) * K^2 - 2*Y^2 * K - Y^2 = 0
%
%   for the positive root K.  X and Y are column vectors of the same size
%   and the returned value K has the same size.

  k = zeros(length(x), 1);
  p = x.^2;
  q = y.^2;
  r = (p + q - 1) / 6;
  fl1 = ~(q == 0 & r <= 0);
  p = p(fl1);
  q = q(fl1);
  r = r(fl1);
  S = p .* q / 4;
  r2 = r.^2;
  r3 = r .* r2;
  disc = S .* (S + 2 * r3);
  u = r;
  fl2 = disc >= 0;
  T3 = S(fl2) + r3(fl2);
  T3 = T3 + (1 - 2 * (T3 < 0)) .* sqrt(disc(fl2));
  T = cbrt(T3);
  u(fl2) = u(fl2) + T + cvmgt(r2(fl2) ./ T, 0, T ~= 0);
  ang = atan2(sqrt(-disc(~fl2)), -(S(~fl2) + r3(~fl2)));
  u(~fl2) = u(~fl2) + 2 * r(~fl2) .* cos(ang / 3);
  v = sqrt(u.^2 + q);
  uv = u + v;
  fl2 = u < 0;
  uv(fl2) = q(fl2) ./ (v(fl2) - u(fl2));
  w = (uv - q) ./ (2 * v);
  k(fl1) = uv ./ (sqrt(uv + w.^2) + w);
end
