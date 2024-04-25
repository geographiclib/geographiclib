function gamdist(axes, num)
  t = triaxial(axes);
  [x, v] = t.cart2rand(num);
  [ell, alp, gam] = t.cart2toellip(x, v);
  modgam = gam;
  l = gam > 0;
  modgam(l) = (sqrt(gam(l)/1));
  l = gam < 0;
  modgam(l) = -(sqrt(-gam(l)/1));
  figure(1); hist(gam,100);
  figure(2); hist(modgam,100);
  dat = [sum(gam < 0) / num, t.kp2, sqrt(t.kp2),1-sqrt(t.k2)]
end
