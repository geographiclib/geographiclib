function runtest(n,k)
  d = 1;
  switch d
    case 1
      dat = load('testset.txt');
      t = triaxial(sqrt([2,1,1/2]));
    case 2
      dat = load('testobl.txt');
      t = triaxial([1,1,1/2]);
    case 3
      dat = load('testpro.txt');
      t = triaxial([2,1,1]);
  end
  m = size(dat, 1);
  dat = [dat, [1:m]'];
  dat = dat(k:n:end, :);
  q = size(dat, 1);
  for i = 1:q
    fprintf(1, '%06d %d %d %d %d %d\n', geod3test(t, dat(i,:)));
  end
end
