% Use Graticule.cpp to generate list of points on graticulepoints.m
% Read in graticulepoints.m
% for i in gauss-krueger-graticule thompson-tm-graticule gauss-krueger-convergence-scale; do pstopnm -stdout -dpi 800 -ppm $i.eps | pnmcrop | pnmscale 0.125 | pnmflip -cw | pnmpad -white -left 10 -top 20 -bottom 20 -right 20 | pnmquant 256 | pnmtopng -compress 9 -interlace > $i.png;done

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
title('Gauss-Krueger transverse Mercator graticule');
xlabel x;
ylabel y;
axis image;axis([0 3.5 -0.5 1]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',-0.5:0.5:1);
set(gcf,'position',[300,300,1000,500]);
print -dpdf gauss-krueger-graticule.pdf
print -depsc gauss-krueger-graticule.eps

figure(2);hold off;
plot(weqat(:,2),weqat(:,1),'r-');
hold on;
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    continue;
    plot(wlonline{i}(:,2),wlonline{i}(:,1),'b-');
  else
    plot(wlonline{i}(:,2),wlonline{i}(:,1),'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if mod(wlats(i),10) == 0,
    continue;
    plot(wlatline{i}(:,2),wlatline{i}(:,1),'b-');
  else
    plot(wlatline{i}(:,2),wlatline{i}(:,1),'g-');
  end
  hold on;
end
for i=1:size(wlons,2),
  if mod(wlons(i),10) == 0,
    plot(wlonline{i}(:,2),wlonline{i}(:,1),'b-');
  else
    continue;
    plot(wlonline{i}(:,2),wlonline{i}(:,1),'g-');
  end
  hold on;
end
for i=1:size(wlats,2),
  if mod(wlats(i),10) == 0,
    plot(wlatline{i}(:,2),wlatline{i}(:,1),'b-');
  else
    continue;
    plot(wlatline{i}(:,2),wlatline{i}(:,1),'g-');
  end
  hold on;
end
hold off;
title('Thompson transverse Mercator graticule');
xlabel v;
ylabel u;
axis image;axis([0 kv 0 ku]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',0:0.5:1.5);
set(gcf,'position',[300,300,1000,500]);
print -dpdf thompson-tm-graticule.pdf
print -depsc thompson-tm-graticule.eps

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
set(gcf,'position',[300,300,1000,500]);
print -dpdf gauss-krueger-convergence-scale.pdf
print -depsc gauss-krueger-convergence-scale.eps
