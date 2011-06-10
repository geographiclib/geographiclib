function [latlong, aux] = geodesicdirect(geodesic, a, r)
%geodesicdirect  Solve direct geodesic problem
%
%   [latlong, aux] = geodesicdirect(geodesic)
%   [latlong, aux] = geodesicdirect(geodesic, a, r)
%
%   geodesic is an M x 4 matrix
%       latitude of point 1 = latlong(:,1) in degrees
%       longitude of point 1 = latlong(:,2) in degrees
%       azimuth at point 1 = latlong(:,3) in degrees
%       distance to point 2 = latlong(:,4) in meters
%
%   latlong is an M x 3 matrix
%       latitude of point 2 = geodesic(:,1) in degrees
%       longitude of point 2 = geodesic(:,2) in degrees
%       azimuth at point 2 = geodesic(:,3) in degrees
%   aux is an M x 4 matrix
%       reduced length = aux(:,1) in meters
%       geodesic scale 1 to 2 = aux(:,2)
%       geodesic scale 2 to 1 = aux(:,3)
%       area under geodesic = aux(:,4) in meters^2
%
%   a = major radius (meters)
%   r = reciprocal flattening (0 means a sphere)
%   If a and r are omitted, the WGS84 values are used.
%
% This is an interface to the GeographicLib C++ routine
%     Geodesic::Direct
% See the documentation on this function for more information:
% http://geographiclib.sf.net/html/classGeographicLib_1_1Geodesic.html
  error('Error: executing .m file instead of compiled routine');
end
% geodesicdirect.m
% Matlab .m file for solving direct geodesic problem
%
% Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and
% licensed under the LGPL.  For more information, see
% http://geographiclib.sourceforge.net/
%
% $Id$
