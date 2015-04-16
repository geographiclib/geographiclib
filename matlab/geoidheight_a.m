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
    geoidname = getenv('GEOGRAPHICLIB_GEOID_NAME');
    if isempty(geoidname)
      geoidname = 'egm96-5';
    end
  end
  if nargin < 3
    geoiddir = getenv('GEOGRAPHICLIB_GEOID_PATH');
    if isempty(geoiddir)
      geoiddir = getenv('GEOGRAPHICLIB_DATA');
      if isempty(geoiddir)
        if ispc
          geoiddir = '/usr/local/share/GeographicLib';
        else
          geoiddir = 'C:/ProgramData/GeographicLib';
        end
      end
      geoiddir = [geoiddir '/geoids'];
    end
  end
  geoid = geoid_load(geoidname, geoiddir);
  height = geoid_height(latlong(:,1), latlong(:,2), geoid);
end
