function [lam12, dlam12, ...
          salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, epsi, domg12] = ...
    Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, f, A3x, C3x)
%LAMBDA12  Solve the hybrid problem

  tiny = sqrt(realmin);
  f1 = 1 - f;
  e2 = f * (2 - f);
  ep2 = e2 / (1 - e2);

  calp1(sbet1 == 0 & calp1 == 0) = -tiny;

  salp0 = salp1 .* cbet1;
  calp0 = hypot(calp1, salp1 .* sbet1);

  ssig1 = sbet1; somg1 = salp0 .* sbet1;
  csig1 = calp1 .* cbet1; comg1 = csig1;
  [ssig1, csig1] = SinCosNorm(ssig1, csig1);

  salp2 = cvmgt(salp0 ./ cbet2, salp1, cbet2 ~= cbet1);
  calp2 = cvmgt(sqrt((calp1 .* cbet1).^2 + ...
                     cvmgt((cbet2 - cbet1) .* (cbet1 + cbet2), ...
                           (sbet1 - sbet2) .* (sbet1 + sbet2), ...
                           cbet1 < -sbet1)) ./ cbet2, ...
                abs(calp1), cbet2 ~= cbet1 | abs(sbet2) ~= -sbet1);
  ssig2 = sbet2; somg2 = salp0 .* sbet2;
  csig2 = calp2 .* cbet2;  comg2 = csig2;
  [ssig2, csig2] = SinCosNorm(ssig2, csig2);

  sig12 = atan2(max(csig1 .* ssig2 - ssig1 .* csig2, 0), ...
                csig1 .* csig2 + ssig1 .* ssig2);

  omg12 = atan2(max(comg1 .* somg2 - somg1 .* comg2, 0), ...
                comg1 .* comg2 + somg1 .* somg2);
  k2 = calp0.^2 * ep2;
  epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
  C3a = C3f(epsi, C3x);
  B312 = SinCosSeries(true, ssig2, csig2, C3a) - ...
         SinCosSeries(true, ssig1, csig1, C3a);
  h0 = -f * A3f(epsi, A3x);
  domg12 = salp0 .* h0 .* (sig12 + B312);
  lam12 = omg12 + domg12;

  [~, dlam12] = ...
      Lengths(epsi, sig12, ...
              ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, false);
  dlam12 = dlam12 .* f1 ./ (calp2 .* cbet2);
  z = calp2 == 0;
  dlam12(z) = - 2 * f1 .* dn1(z) ./ sbet1(z);
end
