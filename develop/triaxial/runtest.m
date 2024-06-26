function runtest(n,k)
  dat = load('Geod3Test.txt');
  m = size(dat, 1);
  dat = [dat, [1:m]'];
  dat = dat(k:n:end, :);
  q = size(dat, 1);
  for i = 1:q
    fprintf(1, '%06d %d %d %d %d %d\n', geod3test(dat(i,:)));
  end
end
