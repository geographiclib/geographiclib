function [geocentric, rot] = geocentricforward(geodetic, a, r)
%geocentricforward  Convert geographic coordinates to geocentric
%
%   [geocentric, rot] = geocentricforward(geodetic);
%   [geocentric, rot] = geocentricforward(geodetic, a, r);
%
%   geodetic is an M x 3 matrix
%       lat = geodetic(:,1) in degrees
%       lon = geodetic(:,2) in degrees
%       h = geodetic(:,3) in meters
%
%   geocentric is an M x 3 matrix
%       x = geocentric(:,1) in meters
%       y = geocentric(:,2) in meters
%       z = geocentric(:,3) in meters
%   rot is an M x 9 matrix
%       M = rot(:,1:9) rotation matrix in row major order.  Pre-multiplying
%           a unit vector in local cartesian coordinates (east, north, up)
%           by M transforms the vector to geocentric coordinates.
%
%   a = major radius (meters)
%   r = reciprocal flattening (0 means a sphere)
%   If a and r are omitted, the WGS84 values are used.
%
%   This is an interface to the GeographicLib C++ routine
%       Geocentric::Forward
%   See the documenation on this function for more information.
  error('Error: executing .m file instead of compiled routine');
end
% geocentricforward.m
% Matlab .m file for geographic to geocentric conversions
%
% Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
