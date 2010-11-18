function latlong = geodesicdirect(geodesic, a, r)
%geodesicdirect  Solve direct geodesic problem
%
%   latlong = geodesicdirect(geodesic)
%   latlong = geodesicdirect(geodesic, a, r)
%
%   geodesic is an M x 4 matrix
%       latitude of point 1 = latlong(:,1) in degrees
%       longitude of point 1 = latlong(:,2) in degrees
%       azimuth at point 1 = latlong(:,3) in degrees
%       distance = latlong(:,4) in meters
%
%   latlong is an M x 7 matrix
%       latitude of point 2 = geodesic(:,1) in degrees
%       longitude of point 2 = geodesic(:,2) in degrees
%       azimuth at point 2 = geodesic(:,3) in degrees
%       reduced length = geodesic(:,4) in meters
%       geodesic scale 1 to 2 = geodesic(:,5)
%       geodesic scale 2 to 1 = geodesic(:,6)
%       area under geodesic = geodesic(:,7) in meters^2
%
%   a = major radius (meters)
%   r = reciprocal flattening (0 means a sphere)
%
%   This is an interface to the GeographicLib C++ routine
%       Geodesic::Direct
%   See the documentation on this function for more information.
  error('Error: executing .m file instead of compiled routine');
end
% geodesicdirect.m
% Matlab .m file for geographic to UTM/UPS conversions
%
% Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
