function mgrs = mgrsforward_a(utmups, prec)
%mgrsforward  Convert UTM/UPS coordinates to MGRS
%
%   mgrs = mgrsforward(utmups)
%   mgrs = mgrsforward(utmups, prec)
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
%   mgrs is a M x 1 cell array of MGRS strings.
  if nargin < 2
    prec = 5;
  end
  mgrs = mgrs_fwd(utmups(:,1), utmups(:,2),utmups(:,3),utmups(:,4),prec);
end
