% Use Graticule.cpp to generate list of points on graticulepoints.m
% Read in graticulepoints.m
% for i in gauss-krueger-graticule thompson-tm-graticule gauss-krueger-convergence-scale gauss-laborde-graticule-a gauss-krueger-graticule-a thompson-tm-graticule-a; do pstopnm -stdout -dpi 600 -ppm $i.eps | pnmcrop | pnmscale 0.125 | pnmflip -cw | pnmpad -white -left 10 -top 10 -bottom 10 -right 10 | pnmquant 256 | pnmtopng -compress 9 -interlace > $i.png;done

lonline{19}=[lonline{19}(1,:); latline{11}(end,:); 0 1];
wlonline{1}=[0 0; ku 0];
wlonline{19}=[ku 0; wlatline{11}(end,:); ku kv];

printp=1;
figure(1);hold off;
for i=1:size(lons,2),
  if mod(lons(i),10) == 0,
    continue;
    plot(lonline{i}(:,1),lonline{i}(:,2),'b-');
  else
    plot(lonline{i}(:,1),lonline{i}(:,2),'g-');
  end
  hold on;
end
for i=1:size(lats,2),
  if mod(lats(i),10) == 0,
    continue;
    plot(latline{i}(:,1),latline{i}(:,2),'b-');
  else
    plot(latline{i}(:,1),latline{i}(:,2),'g-');
  end
  hold on;
end
for i=1:size(lons,2),
  if mod(lons(i),10) == 0,
    plot(lonline{i}(:,1),lonline{i}(:,2),'b-');
  else
    continue;
    plot(lonline{i}(:,1),lonline{i}(:,2),'g-');
  end
  hold on;
end
for i=1:size(lats,2),
  if mod(lats(i),10) == 0,
    plot(latline{i}(:,1),latline{i}(:,2),'b-');
  else
    continue;
    plot(latline{i}(:,1),latline{i}(:,2),'g-');
  end
  hold on;
end
hold off;
title('Gauss-Krueger transverse Mercator graticule (extended)');
xlabel x;
ylabel y;
axis image;axis([0 3.5 -0.5 1]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',-0.5:0.5:1);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf gauss-krueger-graticule.pdf;
  print -depsc gauss-krueger-graticule.eps;
end

figure(2);hold off;
scale=1;
plot(weqat(:,2).*scale,weqat(:,1).*scale,'r-');
hold on;
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    continue;
    plot(wlonline{i}(:,2).*scale,wlonline{i}(:,1).*scale,'b-');
  else
    plot(wlonline{i}(:,2).*scale,wlonline{i}(:,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if mod(wlats(i),10) == 0,
    continue;
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'b-');
  else
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    plot(wlonline{i}(:,2).*scale,wlonline{i}(:,1).*scale,'b-');
  else
    continue;
    plot(wlonline{i}(:,2).*scale,wlonline{i}(:,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if mod(wlats(i),10) == 0,
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'b-');
  else
    continue;
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'g-');
  end
  hold on;
end
hold off;
title('Thompson transverse Mercator graticule (extended)');
xlabel v;
ylabel u;
axis image;axis([0 kv*scale 0 ku*scale]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',0:0.5:1.5);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf thompson-tm-graticule.pdf;
  print -depsc thompson-tm-graticule.eps;
end

figure(3);hold off;
plot(latline{11}(2:end,1),latline{11}(2:end,2),'r-');
hold on;
for i=1:size(ks,2),
  if (mod(ks(i), 5) == 0)
    plot(kline{i}(:,1),kline{i}(:,2),'k-');
  elseif (mod(ks(i), 1) == 0)
    plot(kline{i}(:,1),kline{i}(:,2),'b-');
  else
    plot(kline{i}(:,1),kline{i}(:,2),'g-');
  end
  hold on;
end
for i=1:size(gams,2),
  plot(gamline{i}(:,1),gamline{i}(:,2),'b-');
  hold on;
end
hold off;
title('Gauss-Krueger transverse Mercator convergence and scale');
xlabel x;
ylabel y;
axis image;axis([0 3.5 -0.5 1]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',-0.5:0.5:1);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf gauss-krueger-convergence-scale.pdf;
  print -depsc gauss-krueger-convergence-scale.eps;
end

figure(4);hold off;
scale=2/pi;
ee=0.1;
degree=pi/180;
for lon=[0:10:80 81:1:90],
  if mod(lon,10) == 0,
    continue;
    lat=[0:0.2:10 12:2:80];
    col='b-';
  else
    lat=[0:0.2:10];
    col='g-';
  end
  psi=atanh(sin(lat.*degree)) - ee * atanh(ee * sin(lat.*degree));
  lam=lon*degree;
  x=atan(sinh(psi)./cos(lam));
  y=atanh(sin(lam)./cosh(psi));
  plot(y.*scale,x.*scale,col);
  hold on;
end
for lat=[0:1:10 20:10:80],
  if mod(lat,10) == 0,
    continue;
    lon=[0:2:78 80:0.2:90];
    col='b-';
  else
    lon=[80:0.2:90];
    col='g-';
  end
  psi=atanh(sin(lat.*degree)) - ee * atanh(ee * sin(lat.*degree));
  lam=lon*degree;
  x=atan(sinh(psi)./cos(lam));
  y=atanh(sin(lam)./cosh(psi));
  plot(y.*scale,x.*scale,col);
  hold on;
end
for lon=[0:10:80 81:1:90],
  if mod(lon,10) == 0,
    if lon == 0
      lat=[0 89.999999999];
    elseif lon == 90
      lat=[0.1 89.999999999];
    else
      lat=[0:0.2:10 12:2:80];
    end
    col='b-';
  else
    continue;
    lat=[0:0.2:10];
    col='g-';
  end
  psi=atanh(sin(lat.*degree)) - ee * atanh(ee * sin(lat.*degree));
  lam=lon*degree;
  x=atan(sinh(psi)./cos(lam));
  y=atanh(sin(lam)./cosh(psi));
  plot(y.*scale,x.*scale,col);
  hold on;
end
for lat=[0:1:10 20:10:80],
  if mod(lat,10) == 0,
    lon=[0:2:78 80:0.2:90];
    col='b-';
  else
    continue;
    lon=[80:0.2:90];
    col='g-';
  end
  psi=atanh(sin(lat.*degree)) - ee * atanh(ee * sin(lat.*degree));
  lam=lon*degree;
  x=atan(sinh(psi)./cos(lam));
  y=atanh(sin(lam)./cosh(psi));
  plot(y.*scale,x.*scale,col);
  hold on;
end
hold off;
title('Gauss-Laborde transverse Mercator graticule (scaled)');
xlabel x;
ylabel y;
axis image; axis([0 2.5 0 1]);
set(gca,'XTick',0:0.5:2.5);
set(gca,'YTick',0:0.5:1.0);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf gauss-laborde-graticule-a.pdf;
  print -depsc gauss-laborde-graticule-a.eps;
end

figure(5);hold off;
for i=1:size(lons,2),
  if mod(lons(i),10) == 0,
    continue;
    plot(lonline{i}(1+(lons(i)==90):end,1),lonline{i}(1+(lons(i)==90):end,2),'b-');
  else
    plot(lonline{i}(41:end,1),lonline{i}(41:end,2),'g-');
  end
  hold on;
end
for i=1:size(lats,2),
  if lats(i) < 0,
    continue;
  end
  if mod(lats(i),10) == 0,
    continue;
    plot(latline{i}(:,1),latline{i}(:,2),'b-');
  else
    plot(latline{i}(:,1),latline{i}(:,2),'g-');
  end
  hold on;
end
for i=1:size(lons,2),
  if mod(lons(i),10) == 0,
    plot(lonline{i}(1+(lons(i)==90):end,1),lonline{i}(1+(lons(i)==90):end,2),'b-');
  else
    continue;
    plot(lonline{i}(41:end,1),lonline{i}(41:end,2),'g-');
  end
  hold on;
end
for i=1:size(lats,2),
  if lats(i) < 0,
    continue;
  end
  if mod(lats(i),10) == 0,
    plot(latline{i}(:,1),latline{i}(:,2),'b-');
  else
    continue;
    plot(latline{i}(:,1),latline{i}(:,2),'g-');
  end
  hold on;
end
hold off;
title('Gauss-Krueger transverse Mercator graticule (scaled)');
xlabel x;
ylabel y;
axis image; axis([0 2.5 0 1]);
set(gca,'XTick',0:0.5:2.5);
set(gca,'YTick',0:0.5:1.0);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf gauss-krueger-graticule-a.pdf;
  print -depsc gauss-krueger-graticule-a.eps;
end

figure(6);hold off;
scale=1/ku;
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    continue;
    plot(wlonline{i}(1:end-(wlons(i)==90),2).*scale,wlonline{i}(1:end-(wlons(i)==90),1).*scale,'b-');
  else
    plot(wlonline{i}(41:end,2).*scale,wlonline{i}(41:end,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if wlats(i) < 0,
    continue;
  end
  if mod(wlats(i),10) == 0,
    continue;
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'b-');
  else
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    plot(wlonline{i}(1:end-(wlons(i)==90),2).*scale,wlonline{i}(1:end-(wlons(i)==90),1).*scale,'b-');
  else
    continue;
    plot(wlonline{i}(41:end,2).*scale,wlonline{i}(41:end,1).*scale,'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if wlats(i) < 0,
    continue;
  end
  if mod(wlats(i),10) == 0,
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'b-');
  else
    continue;
    plot(wlatline{i}(:,2).*scale,wlatline{i}(:,1).*scale,'g-');
  end
  hold on;
end
hold off;
title('Thompson transverse Mercator graticule (scaled)');
xlabel x;
ylabel y;
axis image; axis([0 2.5 0 1]);
set(gca,'XTick',0:0.5:2.5);
set(gca,'YTick',0:0.5:1.0);
set(gcf,'position',[300,300,800,700]);
if printp,
  orient landscape;
  print -dpdf thompson-tm-graticule-a.pdf;
  print -depsc thompson-tm-graticule-a.eps;
end
