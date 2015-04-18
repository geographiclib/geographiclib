function h = geoid_height(lat, lon, geoidname, geoiddir)
  persistent saved_geoid
  if nargin == 1 && isempty(lat)
    saved_geoid = []; h = [];
    return
  end
  if nargin == 3 && isstruct(geoidname)
    h = geoid_height_int(lat, lon, geoidname);
  else
    if nargin < 3
      geoidfile = geoid_file;
    elseif nargin < 4
      geoidfile = geoid_file(geoidname);
    else
      geoidfile = geoid_file(geoidname, geoiddir);
    end
    if ~(isstruct(saved_geoid) && strcmp(saved_geoid.file, geoidfile))
      saved_geoid = geoid_load_file(geoidfile);
    end
    h = geoid_height_int(lat, lon, saved_geoid);
  end
end
function height = geoid_height_int(lat, lon, geoid)
  s = size(lat + lon);
  Z = zeros(prod(s),1);
  lat = lat(:) + Z; lon = lon(:) + Z;
  h = geoid.h; w = geoid.w;
  % lat is in [0, h]
  flat = min(max((90 - lat) * (h - 1) / 180, 0), (h - 1));
  % lon is in [0, w)
  flon = mod(lon * w / 360, w);
  flon(isnan(flon)) = 0;
  ilat = min(floor(flat), h - 2);
  ilon = floor(flon);
  flat = flat - ilat; flon = flon - ilon;
  ind = index(ilon + [0,0,1,1], ilat + [0,1,0,1], w, h);
  hf = double(geoid.im(ind));
  height = (1 - flon) .* ((1 - flat) .* hf(:,1) + flat .* hf(:,2)) + ...
           flon       .* ((1 - flat) .* hf(:,3) + flat .* hf(:,4));
  height = geoid.offset + geoid.scale * height;
  height(~(abs(lat) <= 90 & abs(lon) <= 540)) = nan;
  height = reshape(height, s);
end
function ind = index(ix, iy, w, h)
% return 1-based 1d index to w*h array for 0-based 2d indices (ix,iy)
  c = iy < 0;  iy(c) =           - iy(c); ix(c) = ix(c) + w/2;
  c = iy >= h; iy(c) = 2 * (h-1) - iy(c); ix(c) = ix(c) + w/2;
  ix = mod(ix, w);
  ind = 1 + iy + ix * h;
end
