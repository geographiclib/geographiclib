function gedoc
%GEDOC  Great ellipses on an ellipsoid of revolution
%
%   This package includes two routines GEDISTANCE and GERECKON which the
%   inverse and direct problem for great ellipses on the surface of an
%   ellipsoid of revolution.
%
%   The method involves stretching the ellipse along the axis until it
%   becomes a sphere, solving the corresponding great circle problem on the
%   sphere and mapping the results back to the ellipsoid.  Finding the
%   distance involves computing the arc length of an ellipse and this
%   package uses the series solution employed by MATLAB File Exchange
%   package "Geodesics on an ellipsoid of revolution":
%
%     http://www.mathworks.com/matlabcentral/fileexchange/39108
%
%   Consider two points on the ellipsoid at (lat1, lon1) and (lat2, lon2).
%   The plane containing these points and the center of the ellipsoid
%   intersects the ellipsoid on a great ellipse.  The length of the shorter
%   portion of the great ellipse between the two points is s12 and the
%   great ellipse from point 1 to point 2 has forward azimuths azi1 and
%   azi2 at the two end points.
%
%   Two great ellipse problems can be considered:
%     * the direct problem -- given lat1, lon1, s12, and azi1, determine
%       lat2, lon2, and azi2.  This is solved by GERECKON.
%     * the inverse problem -- given lat1, lon1, lat2, lon2, determine s12,
%       azi1, and azi2.  This is solved by GEDISTANCE.
%
%   The parameters of the ellipsoid are specified by the optional ELLIPSOID
%   argument to the routines.  This is a two-element vector of the form
%   [a,e], where a is the equatorial radius, e is the eccentricity e =
%   sqrt(a^2-b^2)/a, and b is the polar semi-axis.  Typically, a and b are
%   measured in meters and the linear and area quantities returned by the
%   routines are then in meters and meters^2.  However, other units can be
%   employed.  If ELLIPSOID is omitted, then the WGS84 ellipsoid (more
%   precisely, the value returned by DEFAULTELLIPSOID) is assumed [6378137,
%   0.0818191908426215] corresponding to a = 6378137 meters and a
%   flattening f = (a-b)/a = 1/298.257223563.  The flattening and
%   eccentricity are related by
%
%       e = sqrt(f * (2 - f))
%       f = e^2 / (1 + sqrt(1 - e^2))
%
%   (The functions ECC2FLAT and FLAT2ECC implement these conversions.)  For
%   a sphere, set e = 0; for a prolate ellipsoid (b > a), specify e as a
%   pure imaginary number.
%
%   All angles (latitude, longitude, azimuth) are measured in degrees with
%   latitudes increasing northwards, longitudes increasing eastwards, and
%   azimuths measured clockwise from north.  For a point at a pole, the
%   azimuth is defined by keeping the longitude fixed, writing lat =
%   +/-(90-eps), and taking the limit eps -> 0+.
%
%   Restrictions on the inputs:
%     * All latitudes must lie in [-90, 90].
%     * All longitudes and azimuths must lie in [-540, 540).  On output,
%       these quantities lie in [-180, 180).
%     * The distance s12 is unrestricted.  This allows great ellipses to wrap
%       around the ellipsoid.
%     * The equatorial radius, a, must be positive.
%     * The eccentricity, e, should be satisfy abs(e) < 0.2 in order to
%       retain full accuracy (this corresponds to flattenings satisfying
%       abs(f) <= 1/50, approximately).  This condition holds for most
%       applications in geodesy.
%
%    Larger values of e can be used with a corresponding drop in accuracy.
%    The following table gives the approximate maximum error in GEDISTANCE
%    and GERECKON (expressed as a distance) for an ellipsoid with the same
%    major radius as the WGS84 ellipsoid and different values of the
%    flattening.
%
%         |f|     error
%         0.01    25 nm
%         0.02    30 nm
%         0.05    10 um
%         0.1    1.5 mm
%         0.2    300 mm
%
%   In order to compute intermediate points on a great ellipse, proceed as
%   in the following example which plots the track from Sydney to
%   Valparaiso.
%
%       lat1 = -33.83; lon1 = 151.29;
%       lat2 = -33.02; lon2 = -71.64;
%       [s12,  azi1] = gedistance(lat1, lon1, lat2, lon2);
%       [lats, lons] = gereckon(lat1, lon1, s12 * [0:100]/100, azi1);
%       plot(lons+360*(lons<0), lats);
%
%   The restriction on e above arises because the formulation is in terms
%   of series expansions in e^2.  The exact solutions (valid for any e) can
%   be expressed in terms of elliptic integrals.
%
%   See also GEDISTANCE, GERECKON, DEFAULTELLIPSOID, ECC2FLAT, FLAT2ECC.

% Copyright (c) Charles Karney (2014) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.38.

  help gedoc
end
