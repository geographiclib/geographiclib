function band = LatitudeBand(lat)
  band = max(-10, min(9, floor(lat / 8)));
  band(~(abs(lat) <= 90)) = nan;
end

