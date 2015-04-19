function [latlong, scale] = utmupsreverse(utmups)
%UTMUPSREVERSE  Convert UTM/UPS coordinates to geographic
%
%   [latlong, scale] = UTMUPSREVERSE(utmups)
%
%   This is a legacy function to replace a compiled interface function of
%   the same name.  This now calls UTMUPS_INV which is implemented as
%   native Matlab code.
%
%   utmups is an M x 4 matrix
%       easting = utmups(:,1) in meters
%       northing = utmups(:,2) in meters
%       zone = utmups(:,3)
%       hemi = utmups(:,4)
%
%   zone = 0 for UPS, zone = [1,60] for UTM
%   hemi = 0 for southern hemisphere, hemi = 1 for northern hemisphere.
%
%   latlong is an M x 2 matrix
%       latitude = latlong(:,1) in degrees
%       longitude = latlong(:,2) in degrees
%   scale is an M x 2 matrix
%       gamma = scale(:,1) meridian convergence in degrees
%       k = scale(:,2) scale
%
%   See also UTMUPS_INV, UTMUPS_FWD, UTMUPSFORWARD.

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  [lat, lon, gam, k] = ...
      utmups_inv(utmups(:,1), utmups(:,2), utmups(:,3), utmups(:,4));
  latlong = [lat, lon];
  scale = [gam, k];
end
