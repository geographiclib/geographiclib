function mgrs = mgrsforward(utmups)
%mgrsforward  Convert UTM/UPS coordinates to MGRS
%
%   mgrs = mgrsforward(utmups);
%   mgrs = mgrsforward(utmups, prec);
%
%   utmups is an M x 4 matrix
%       easting = utmups(:,1) in meters
%       northing = utmups(:,2) in meters
%       zone = utmups(:,3)
%       hemi = utmups(:,4)
%
%   zone = 0 for UPS, zone = [1,60] for UTM
%   hemi = 0 for southern hemisphere, hemi = 1 for northern hemisphere
%   prec = half the number of trailing digits in the MGRS string
%          (default 5)
%
%   mgrs is a vector of M strings of length 5 + 2 * prec.
%   For UPS coordinates the string begins with 2 blanks.
  error('Error: executing .m file instead of compiled routine');
end.
% mgrsforward.m
% Matlab .m file for geographic to UTM/UPS conversions
%
% Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
