function [s12b, m12b, m0, M12, M21] = ...
      Lengths(epsi, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, ...
              cbet1, cbet2, scalp, ep2)
%LENGTHS  Compute various lengths associate with a geodesic

  if isempty(sig12)
    s12b = [];
    m12b = [];
    m0 = [];
    M12 = [];
    M21 = [];
    return
  end

  C1a = C1f(epsi);
  C2a = C2f(epsi);
  A1m1 = A1m1f(epsi);
  AB1 = (1 + A1m1) .* (SinCosSeries(true, ssig2, csig2, C1a) - ...
                       SinCosSeries(true, ssig1, csig1, C1a));
  A2m1 = A2m1f(epsi);
  AB2 = (1 + A2m1) .* (SinCosSeries(true, ssig2, csig2, C2a) - ...
                       SinCosSeries(true, ssig1, csig1, C2a));
  m0 = A1m1 - A2m1;
  J12 = m0 .* sig12 + (AB1 - AB2);
  m12b = dn2 .* (csig1 .* ssig2) - dn1 .* (ssig1 .* csig2) - ...
         csig1 .* csig2 .* J12;
  s12b = (1 + A1m1) .* sig12 + AB1;
  if scalp
    csig12 = csig1 .* csig2 + ssig1 .* ssig2;
    t = ep2 * (cbet1 - cbet2) .* (cbet1 + cbet2) ./ (dn1 + dn2);
    M12 = csig12 + (t .* ssig2 - csig2 .* J12) .* ssig1 ./ dn1;
    M21 = csig12 - (t .* ssig1 - csig1 .* J12) .* ssig2 ./ dn2;
  else
    M12 = sig12 + NaN; M21 = M12;
  end
end
