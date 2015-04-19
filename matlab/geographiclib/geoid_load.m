function geoid = geoid_load(name, dir)
%GEOID_LOAD  Data a geoid model
%
%   geoid = GEOID_LOAD
%   geoid = GEOID_LOAD(geoidname)
%   geoid = GEOID_LOAD(geoidname, geoiddir)
%
%   Loads geoid data into the workspace.  The possible geoids are
%
%       egm84-30  egm84-15
%       egm96-15  egm96-5
%       egm2008-5 egm2008-2_5 egm2008-1
%
%   The first part of the name is the geoid model.  The second part gives the
%   resolution of the gridded data (in arc-seconds).
%
%   The geoid can be overridden by specifying geoidname.  If geoidname is not
%   specified, the environment variable GEOGRAPHICLIB_GEOID_NAME is used; if
%   this is not defined then egm96-5 is used.  GEOID_HEIGHT looks in the
%   directory geoiddir for the geoid data; if this is not specified, it uses
%   the environment variable GEOGRAPHICLIB_GEOID_PATH; if this is not defined,
%   it appends "/geoids" to the environment variable GEOGRAPHICLIB_DATA;
%   finally, it tries the default directory names
%   /usr/local/share/GeographicLib/geoids or
%   C:/ProgramData/GeographicLib/geoids.
%
%   The returned geoid can be passed to GEOID_HEIGHT to determine the height
%   of the geoid.
%
%   Information on downloading and installing the data for the supported
%   geoid models is available at
%
%     http://geographiclib.sf.net/html/geoid.html#geoidinst
%
%   See also GEOID_HEIGHT.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  if nargin < 1
    file = geoid_file;
  elseif nargin < 2
    file = geoid_file(name);
  else
    file = geoid_file(name, dir);
  end
  geoid = geoid_load_file(file);
end
