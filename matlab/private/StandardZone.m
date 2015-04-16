function zone = StandardZone(lat, lon, setzone)

  INVALID = -4;
  UTM = -2;
  STANDARD = -1;
  MINZONE = 0;
  MAXZONE = 60;
  UPS = 0;
  
  if nargin < 3
    setzone = STANDARD;
  end

  zone = floor(setzone) + zeros(size(lat));
  zone(zone > MAXZONE | zone < INVALID) = INVALID;
  g = zone < MINZONE & zone ~= INVALID;
  c = abs(lat) <= 90 & abs(lon) <= 540;
  zone(g & ~c) = INVALID;
  g = g & c;
  c = zone == UTM | (lat >= -80 & lat < 84);
  u = g & c;
  ilon = mod(floor(lon(u)) + 180, 360) - 180;
  z = floor((ilon + 186) / 6);
  band = LatitudeBand(lat(u));
  % Norway exception
  exception = z == 31 & floor(lat(u) / 8) == 7 & ilon >= 3;
  z(exception) = 32;
  % Svalbard exception
  exception = lat(u) >= 72 & ilon >= 0 & ilon < 42; 
  z(exception) = 2 * floor((ilon(exception) + 183)/12) + 1;
  zone(u) = z;
  zone(g & ~c) = UPS;
end
