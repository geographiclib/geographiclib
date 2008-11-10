% Use Graticule.cpp to generate list of points on graticule-points.m
% Read in graticule-points.m
% Use print to file, landscape, 0.25in margins all round
% fig 1 -> gauss-krueger-graticule.ps
% fig 2 -> thompson-tm-graticule.ps
% fig 3 -> gauss-krueger-convergence-scale.ps
% for i in gauss-krueger-graticule thompson-tm-graticule gauss-krueger-convergence-scale; do sed -e 's/^\(%%Creator:.*\.\) .*$/\1/' -e 's&^%%Title: .*/&%%Title: &' $i.ps > $i.eps;done
% for i in gauss-krueger-graticule thompson-tm-graticule gauss-krueger-convergence-scale; do ps2pdf $i.eps;done
% for i in gauss-krueger-graticule thompson-tm-graticule gauss-krueger-convergence-scale; do pstopnm -stdout -dpi 600 -pgm $i.eps | pnmcrop | pnmscale 0.15 | pnmflip -cw | pnmpad -white -left 10 -top 20 -bottom 20 -right 20 | pnmtopng -compress 9 > $i.png;done

figure(1);hold off;
for i=1:size(lons,2),
  plot(lonline{i}(:,1),lonline{i}(:,2),'k-');
  hold on;
end
for i=1:size(lats,2),
  plot(latline{i}(:,1),latline{i}(:,2),'k-');
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

figure(2);hold off;
for i=1:size(wlons,2),
  plot(wlonline{i}(:,2),wlonline{i}(:,1),'k-');
  hold on;
end
for i=1:size(wlats,2),
  plot(wlatline{i}(:,2),wlatline{i}(:,1),'k-');
  hold on;
end
plot(weqat(:,2),weqat(:,1),'k-.');
hold off;
title('Thompson transverse Mercator graticule');
xlabel v;
ylabel u;
axis image;axis([0 kv 0 ku]);
set(gca,'XTick',0:0.5:3.5);
set(gca,'YTick',0:0.5:1.5);
set(gcf,'position',[300,300,1000,500]);

figure(3);hold off;
for i=1:size(ks,2),
  plot(kline{i}(:,1),kline{i}(:,2),'k-');
  hold on;
end
for i=1:size(gams,2),
  plot(gamline{i}(:,1),gamline{i}(:,2),'k-');
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

% Save as eps
% for i in 1 2 3; do
%   pstopnm -dpi 2400 fig$i.eps
%   ppmtopgm fig$i.eps001.ppm | pnmcrop | pnmscale 0.1 |
%   pnmflip -cw | pnmtopng -compress 9 > fig${i}b.png
% done &

% Alternatively, print landscape 0.25in margins all round, print to file
% for i in 1 2 3;do ps2pdf fig$i.ps;done

% for i in 1 2 3; do
%   pstopnm -dpi 1200 fig$i.ps
%   ppmtopgm fig${i}001.ppm | pnmcrop | pnmscale 0.1 |
%   pnmflip -cw | pnmtopng -compress 9 > fig${i}.png
% done &
