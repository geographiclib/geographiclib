function datout = prunerand(dat)
  l = dat(:,7) > 0;
  dat = dat(l,:);
  bet1 = abs(dat(:,1)); omg1 = abs(abs(dat(:,2))-90);
  bet2 = abs(dat(:,3)); omg2 = abs(abs(dat(:,4))-90);
  [bet1,bet2] = deal(min(bet1,bet2), max(bet1,bet2));
  [omg1,omg2] = deal(min(omg1,omg2), max(omg1,omg2));
  s12 = dat(:,7);
  key = s12+4*(omg2+100*(bet2+100*(omg1+100*bet1)));
  [skey, k] = sort(key);
  dkey = diff([skey;inf]);
  qstop = find(dkey > 0);
  qstart = [1;qstop(1:end-1)+1];
  n = length(qstop);
  datout = zeros(n, 7);
  for m=1:n
    p = randi([qstart(m),qstop(m)]);
    kk = k(p);
    datout(m,:) = dat(k(p),:);
  end
end
