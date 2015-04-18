function filename = geoid_file(name, dir)
  if nargin < 1
    name = getenv('GEOGRAPHICLIB_GEOID_NAME');
    if isempty(name)
      name = 'egm96-5';
    end
  end
  if nargin < 2
    dir = getenv('GEOGRAPHICLIB_GEOID_PATH');
    if isempty(dir)
      dir = getenv('GEOGRAPHICLIB_DATA');
      if isempty(dir)
        if ispc
          dir = 'C:/ProgramData/GeographicLib';
        else
          dir = '/usr/local/share/GeographicLib';
        end
      end
      dir = [dir '/geoids'];
    end
  end
  filename = [dir '/' name '.pgm'];
end

