function runtest(d, n, k)
% Run with, e.g.
%   for ((k = 1; k <= 8; ++k)); do
%     matlab-cli -batch "runtest(5,8,$k)" > errs-$k.txt& sleep 1
%   done
  narginchk(3, 3);
  assert(n >= 1 && k >= 1 && k <= n && d >=1 && d <= 7);
  switch d
    case 1
      dat = load('testset.txt');
      t = triaxial([1,3/2,1/3,2/3]);
    case 2
      dat = load('testobl.txt');
      t = triaxial([1,1,1/2]);
    case 3
      dat = load('testpro.txt');
      t = triaxial([2,1,1]);
    case 4
      dat = load('testspha.txt');
      t = triaxial([1,0,1,0]);
    case 5
      dat = load('testsphb.txt');
      t = triaxial([1,0,2/3,1/3]);
    case 6
      dat = load('testsphc.txt');
      t = triaxial([1,0,1/3,2/3]);
    case 7
      dat = load('testsphd.txt');
      t = triaxial([1,0,0,1]);
  end
  m = size(dat, 1);
  dat = [dat, [1:m]'];
  dat = dat(k:n:end, :);
  q = size(dat, 1);
  for i = 1:q
    fprintf(1, '%06d %d %d %d %d %d\n', geod3test(t, dat(i,:)));
  end
end
