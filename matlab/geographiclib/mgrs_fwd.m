function mgrs = mgrs_fwd(x, y, zone, northp, prec)
%MGRS_FWD  Convert UTM/UPS coordinates to MGRS
%
%   mgrs = MGRS_FWD(x, y, zone, northp)
%   mgrs = MGRS_FWD(x, y, zone, northp, prec)
%
%   converts from UTM/UPS coordinates to MGRS.  x, y are the easting and
%   northing (in meters); zone is the UTM zone, in [1,60] or 0 for UPS;
%   northp is true (false) for the northern (southern) hemisphere.  prec
%   gives the precision of the grid reference; the default is 5 giving 1m
%   precision.  A value of -1 mean that only the grid zone is returned.
%   The maximum allowed value of prec is 11 (denoting 1um precision).  The
%   MGRS references are returned in a cell array of strings.  x, y, zone,
%   northp, prec can be scalars or arrays of the same size.  Values that
%   can't be converted to MGRS return the "invalid" string "INV".
%
%   See also MGRS_INV.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 5
    prec = 5;
  end
  zone = floor(zone);
  prec = floor(prec);
  try
    s = size(x + y + zone + northp + prec);
  catch err
    error('x, y, zone, northp, prec have incompatible sizes')
  end
  num = prod(s);
  if num == 0, mgrs = cell(0); return, end
  Z = zeros(num, 1);
  x = x(:) + Z; y = y(:) + Z; zone = zone(:) + Z;
  northp = northp(:) + Z; prec = prec(:) + Z;
  prec(~(prec >= -1 & prec <= 11)) = -2;
  mgrs = repmat('INV',num,1);
  if ~any(prec >= -1), mgrs = reshape(cellstr(mgrs), s); return, end
  maxprec = max(prec);
  mgrs = [mgrs, repmat(' ',num, 2 + 2*maxprec)];
  minprec = min(prec(prec >= -1));
  for p = minprec:maxprec
    in = prec == p;
    if ~any(in)
      continue
    end
    t = mgrs_fwd_p(x(in), y(in), zone(in), northp(in), p);
    mgrs(in,1:(5 + 2*p)) = t;
  end
  mgrs = reshape(cellstr(mgrs), s);
end

function mgrs = mgrs_fwd_p(x, y, zone, northp, prec)
  num = size(x, 1);
  delta = 10e-9;
  mgrs = repmat('INV',num,1);
  mgrs = [mgrs, repmat(' ',num, 2 + 2*prec)];
  utm = zone >= 1 & zone <= 60;
  y(utm & ~northp) = y(utm & ~northp) - 100e5;
  northp(utm) = 1;
  utm = utm & x >= 1e5 & x <= 9e5 & y >= -90e5 & y <= 95e5;
  x(utm & x ==  9e5) =  9e5 - delta;
  y(utm & y == 95e5) = 95e5 - delta;
  upsn = zone == 0 &  northp;
  upss = zone == 0 & ~northp;
  upsn = upsn & x >= 13e5 & x <= 27e5 & y >= 13e5 & y <= 27e5;
  x(upsn & x == 27e5) = 27e5 - delta;
  y(upsn & y == 27e5) = 27e5 - delta;
  upss = upss & x >=  8e5 & x <= 32e5 & y >=  8e5 & y <= 32e5;
  x(upss & x == 32e5) = 32e5 - delta;
  y(upss & y == 32e5) = 32e5 - delta;
  t = mgrs_fwd_utm(x(utm), y(utm), zone(utm), prec); mgrs(utm,:) = t;
  t = mgrs_fwd_upsn(x(upsn), y(upsn), prec); mgrs(upsn,1:end-2) = t;
  t = mgrs_fwd_upss(x(upss), y(upss), prec); mgrs(upss,1:end-2) = t;
end

function mgrs = mgrs_fwd_utm(x, y, zone, prec)
  mgrs = char(zeros(length(x), 5 + 2 * prec) + ' ');
  if isempty(x), return, end
  mgrs(:,1) = '0' + floor(zone / 10);
  mgrs(:,2) = '0' + mod(zone, 10);
  ys = y / 1e5;
  latp = 0.901 * ys + ((ys > 0) * 2 - 1) * 0.135;
  late = 0.902 * ys .* (1 - 1.85e-6 * ys .* ys);
  latp(abs(ys) < 1) = 0.9 * ys(abs(ys) < 1);
  late(abs(ys) < 1) = latp(abs(ys) < 1);
  band = LatitudeBand(latp);
  bande =  LatitudeBand(late);
  c = band ~= bande;
  band(c) = LatitudeBand(utmups_inv(x(c), y(c), zone(c), 1));
  latband = 'CDEFGHJKLMNPQRSTUVWX';
  mgrs(:,3) = latband(band + 11);
  if prec < 0, return, end
  xh = floor(x / 1e5); yh = floor(y / 1e5);
  utmcols = ['ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ'];
  utmrow = 'ABCDEFGHJKLMNPQRSTUV';
  mgrs(:,4) = utmcols(mod(zone - 1, 3) * 8 + xh);
  mgrs(:,5) = utmrow(mod(yh + mod(zone - 1, 2) * 5, 20) + 1);
  if prec == 0, return, end
  x = x - 1e5 * xh; y = y - 1e5 * yh;
  xy = formatnum(x, y, prec);
  mgrs(:,5+(1:2*prec)) = xy;
end

function mgrs = mgrs_fwd_upsn(x, y, prec)
  mgrs = char(zeros(length(x), 3 + 2 * prec) + ' ');
  if isempty(x), return, end
  upsband = 'YZ';
  xh = floor(x / 1e5);
  eastp = xh >= 20;
  mgrs(:,1) = upsband(eastp + 1);
  if prec < 0, return, end
  yh = floor(y / 1e5);
  upscols = ['RSTUXYZ', 'ABCFGHJ'];
  upsrow = 'ABCDEFGHJKLMNP';
  mgrs(:,2) = upscols(eastp * 7 + xh - cvmgt(20, 13, eastp) + 1);
  mgrs(:,3) = upsrow(yh - 13 + 1);
  if prec == 0, return, end
  x = x - 1e5 * xh; y = y - 1e5 * yh;
  xy = formatnum(x, y, prec);
  mgrs(:,3+(1:2*prec)) = xy;
end

function mgrs = mgrs_fwd_upss(x, y, prec)
  mgrs = char(zeros(length(x), 3 + 2 * prec) + ' ');
  if isempty(x), return, end
  upsband = 'AB';
  xh = floor(x / 1e5);
  eastp = xh >= 20;
  mgrs(:,1) = upsband(eastp + 1);
  if prec < 0, return, end
  yh = floor(y / 1e5);
  upscols = ['JKLPQRSTUXYZ', 'ABCFGHJKLPQR'];
  upsrow = 'ABCDEFGHJKLMNPQRSTUVWXYZ';
  mgrs(:,2) = upscols(eastp * 12 + xh - cvmgt(20, 8, eastp) + 1);
  mgrs(:,3) = upsrow(yh - 8 + 1);
  if prec == 0, return, end
  x = x - 1e5 * xh; y = y - 1e5 * yh;
  xy = formatnum(x, y, prec);
  mgrs(:,3+(1:2*prec)) = xy;
end

function xy = formatnum(x, y, prec)
  if (prec < 5)
    x = x / 10 ^ (5 - prec); y = y / 10 ^ (5 - prec);
  elseif (prec > 5)
    x = x * 10 ^ (prec - 5); y = y * 10 ^ (prec - 5);
  end
  xy = [num2str(floor(x), ['%0', int2str(prec), 'd']), ...
        num2str(floor(y), ['%0', int2str(prec), 'd'])];
end

function band = LatitudeBand(lat)
  band = max(-10, min(9, floor(lat / 8)));
  band(~(abs(lat) <= 90)) = nan;
end
