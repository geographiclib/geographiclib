function height = geoidheight_a(latlong, geoidname, geoiddir)
%geoidheight  Compute geoid height
%
%   height = geoidheight(latlong)
%   height = geoidheight(latlong, geoidname, geoiddir)
%
%   latlong is an M x 2 matrix
%       latitude = latlong(:,1) in degrees
%       longitude = latlong(:,2) in degrees
%   geoidname is the name of the geoid; choices are (default egm96-5)
%       egm84-30  egm84-15
%       egm96-15  egm96-5
%       egm2008-5 egm2008-2_5 egm2008-1
%   geoiddir is the directory containing the geoid models (default empty
%       string meaning system default)
%
%   height is an M x 1 matrix
%       geoidheight = height(:,1) height of geoid in meters
%
  if nargin < 2
    height = geoid_height(latlong(:,1), latlong(:,2));
  elseif nargin < 3
    height = geoid_height(latlong(:,1), latlong(:,2), geoidname);
  else
    height = geoid_height(latlong(:,1), latlong(:,2), geoidname, geoiddir);
  end
end
