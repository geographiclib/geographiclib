function [s, newdat, cats] = stats(dat, targetstat)
  if nargin < 2
  targetstat=[
      3    100    100     30     20     30   3000      0
      0    900   1500    100    100    500  10000      0
      0      0   1500    100    100    500  10000      0
      0      0      0    890    100     20  10000      0
      0      0      0      0    800     20  10000      0
      0      0      0      0      0    700  20000      0
      0      0      0      0      0      0 425287      0
      0      0      0      0      0      0      0   3600
         ];
  end
  bet1 = abs(dat(:,1)); omg1 = abs(abs(dat(:,2))-90);
  bet2 = abs(dat(:,4)); omg2 = abs(abs(dat(:,5))-90);
  cat1 = cat(bet1, omg1);
  cat2 = cat(bet2, omg2);
  [cat1, cat2] = deal(min(cat1, cat2), max(cat1, cat2));
  conj = abs(dat(:,1)+dat(:,4)) <= 2 & ...
         abs(abs(dat(:,3))-90) <= 2 & abs(abs(dat(:,6))-90) <= 2 & ...
         cat1 == 7 & cat2 == 7;
  % category 8 if two general point and near conjugate
  cat1(conj) = 8; cat2(conj) = 8;
  s = zeros(8,8);
  newdat = zeros(0,10);
  for q = 1:8
    for p = 1:8
      s(q, p) = sum(cat1 == q & cat2 == p);
      if targetstat(q, p) > 0
        tmp = dat(cat1 == q & cat2 == p, :);
        newdat = [newdat;tmp(randperm(s(q, p), targetstat(q, p)),:)];
      end
    end
  end
  
  cats = [cat1, cat2, [1:length(cat1)]'];
end
function res = cat(bet,omg)
  % 7 categories
  umb = bet == 90 & omg == 90;          % umbilical point
  meda = omg == 90;                     % omg = 0 or 180 medium ellipse
  medb = bet ==90;                      % bet = +/-90 medium ellipse
  equ = bet == 0;                       % equator
  min = omg == 0;                       % minor ellipse
  numb = bet > 85 & omg > 85;           % near umbilical point
  gen = bet > -1;                       % general
  res = 0*bet;
  res(gen) = 7;
  res(numb) = 6;
  res(min) = 5;
  res(equ) = 4;
  res(medb) = 3;
  res(meda) = 2;
  res(umb) = 1;
end
