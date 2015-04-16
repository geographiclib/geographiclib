function [utmups, scale] = utmupsforward_a(latlong, setzone)
%utmupsforward  Convert geographic coordinates to UTM/UPS
%
%   [utmups, scale] = utmupsforward(latlong)
%   [utmups, scale] = utmupsforward(latlong, setzone)
%
%   latlong is an M x 2 matrix
%       latitude = latlong(:,1) in degrees
%       longitude = latlong(:,2) in degrees
%
%   utmups is an M x 4 matrix
%       easting = utmups(:,1) in meters
%       northing = utmups(:,2) in meters
%       zone = utmups(:,3)
%       hemi = utmups(:,4)
%   scale is an M x 2 matrix
%       gamma = scale(:,1) meridian convergence in degrees
%       k = scale(:,2) scale
%
%   zone = 0 for UPS, zone = [1,60] for UTM
%   hemi = 0 for southern hemisphere, hemi = 1 for northern hemisphere
%
%   setzone is an zone override flag with legal values
%        0, use UPS
%        [1,60], use the corresponding UTM zone
%       -1, use the standard assigment (the default)
%       -2, use closest UTM zone
  if nargin < 2
    setzone = -1;
  end
  [x, y, zone, northp, gam, k] = ...
      utmups_fwd(latlong(:,1), latlong(:,2), setzone);
  utmups = [x, y, zone, northp];
  scale = [gam, k];
end
