function [sig12, salp1, calp1, salp2, calp2] = ...
      InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, f, A3x)
%INVERSESTART  Compute a starting point for Newton's method

  N = length(sbet1);
  f1 = 1 - f;
  e2 = f * (2 - f);
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  tol0 = eps;
  tol1 = 200 * tol0;
  tol2 = sqrt(eps);
  etol2 = 0.001 * tol2 / max(0.1, sqrt(abs(e2)));
  xthresh = 1000 * tol2;

  sig12 = - ones(N, 1); salp2 = NaN(N, 1); calp2 = NaN(N, 1);
  sbet12 = sbet2 .* cbet1 - cbet2 .* sbet1;
  cbet12 = cbet2 .* cbet1 + sbet2 .* sbet1;
  sbet12a = sbet2 .* cbet1 + cbet2 .* sbet1;
  s = cbet12 >= 0 & sbet12 < 0.5 & lam12 <= pi / 6;
  omg12 = lam12;
  omg12(s) = omg12(s) ./ (f1 * (dn1(s) + dn2(s)) / 2);
  somg12 = sin(omg12); comg12 = cos(omg12);

  salp1 = cbet2 .* somg12;
  t = cbet2 .* sbet1 .* somg12.^2;
  calp1 = cvmgt(sbet12  + t ./ (1 + comg12), ...
                sbet12a - t ./ (1 - comg12), ...
                comg12 >= 0);

  ssig12 = hypot(salp1, calp1);
  csig12 = sbet1 .* sbet2 + cbet1 .* cbet2 .* comg12;

  s = s & ssig12 < etol2;
  salp2(s) = cbet1(s) .* somg12(s);
  calp2(s) = sbet12(s) - cbet1(s) .* sbet2(s) .* somg12(s).^2 ./ ...
      (1 + comg12(s));
  [salp2, calp2] = SinCosNorm(salp2, calp2);
  sig12(s) = atan2(ssig12(s), csig12(s));

  s = ~(s | abs(n) > 0.1 | csig12 >= 0 | ssig12 >= 6 * abs(n) * pi * cbet1.^2);

  if any(s)
    if f >= 0
      k2 = sbet1(s).^2 * ep2;
      epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
      lamscale = f * cbet1(s) .* A3f(epsi, A3x) * pi;
      betscale = lamscale .* cbet1(s);
      x = (lam12(s) - pi) ./ lamscale;
      y = sbet12a(s) ./ betscale;
    else
      cbet12a = cbet2(s) .* cbet1(s) - sbet2(s) .* sbet1(s);
      bet12a = atan2(sbet12a(s), cbet12a);
      [~, m12b, m0] = ...
          Lengths(n, pi + bet12a, ...
                  sbet1(s), -cbet1(s), dn1(s), sbet2(s), cbet2(s), dn2(s), ...
                  cbet1(s), cbet2(s), false);
      x = -1 + m12b ./ (cbet1(s) .* cbet2(s) .* m0 * pi);
      betscale = cvmgt(sbet12a(s) ./ x, - f * cbet1(s).^2 * pi, x < -0.01);
      lamscale = betscale ./ cbet1(s);
      y = (lam12(s) - pi) ./ lamscale;
    end
    k = Astroid(x, y);
    if f >= 0
      omg12a = -x .* k ./ (1 + k);
    else
      omg12a = -y .* (1 + k) ./ k;
    end
    omg12a = lamscale .* omg12a;
    somg12 = sin(omg12a); comg12 = -cos(omg12a);
    salp1(s) = cbet2(s) .* somg12;
    calp1(s) = sbet12a(s) - cbet2(s) .* sbet1(s) .* somg12.^2 ./ (1 - comg12);

    str = y > -tol1 & x > -1 - xthresh;
    if any(str)
      salp1s = salp1(s); calp1s = calp1(s);
      if f >= 0
        salp1s(str) = min(1, -x(str));
        calp1s(str) = -sqrt(1 - salp1s(str).^2);
      else
        calp1s(str) = max(cvmgt(0, -1, x(str) > -tol1), x(str));
        salp1s(str) = sqrt(1 - calp1s(str).^2);
      end
      salp1(s) = salp1s; calp1(s) = calp1s;
    end
  end

  calp1(salp1 <= 0) = 0; salp1(salp1 <= 0) = 1;
  [salp1, calp1] = SinCosNorm(salp1, calp1);
end
