function errs = geod3test(dat)
  t = triaxial(sqrt([2,1,1/2]));
  ell1 = dat(:,1:2);
  alp1 = dat(:,3);
  ell2 = dat(:,4:5);
  alp2 = dat(:,6);
  s12 = dat(:,7);
  ind = dat(:,11);
  [r1, v1] = t.elliptocart2(ell1, alp1);
  [r2, v2] = t.elliptocart2(ell2, alp2);
  s12a = t.distance(r1, r2);
  r2a = r2;
  m = length(s12);
  for k = 1:m
    [r2a(k, :), v2a(k,:)] = t.reckon(r1(k, :), v1(k, :), s12(k));
    [r1a(k, :), v1a(k,:)] = t.reckon(r2(k, :), v2(k, :), -s12(k));
  end
  err1 = abs(s12a - s12);
  err2 = sqrt(sum((r2 - r2a).^2, 2));
  err3 = sqrt(sum((v2 - v2a).^2, 2));
  err4 = sqrt(sum((r1 - r1a).^2, 2));
  err5 = sqrt(sum((v1 - v1a).^2, 2));
  errs = ceil([err1, err2, err3, err4, err5] / (eps/2));
  errs = [ind, errs];
end
