function height = geoidheight(latlong)
%geoidheight  Compute geoid height
%
%   CAUTION: THIS CRASHES MATLAB!!  There is some incompatibility with
%   the Geoid constructor and the Matlab environment.  The crash (with
%   Matlab 2008a under Windows) happens on the second call.
%
%   height = geoidheight(latlong);
%   height = geoidheight(latlong, geoidname);
%   height = geoidheight(latlong, geoidname, geoiddir);
%
%   latlong is an M x 2 matrix
%       latitude = latlong(:,1) in degrees
%       longitude = latlong(:,2) in degrees
%   geoidname is the name of the geoid (default egm96-5)
%   geoiddir is the direcortory containing the geoid models (default empty
%       string meaning system default)
%
%   height is an M x 3 matrix
%       geoidheight = utmups(:,1) in meters
%       gradn = utmups(:,2) gradient of height in northerly direction
%       grade = utmups(:,3) gradient of height in easterly direction
  error('Error: executing .m file instead of compiled routine');
end.
% geoidheight.m
% Matlab .m file for looking up geoid heights
%
% Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
