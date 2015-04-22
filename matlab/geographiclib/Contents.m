% GeographicLib toolbox
% Version 1.42 2015-04-22
%
%   This toolbox provides native MATLAB implementations of a subset of the
%   C++ library, GeographicLib.  (But note that the great ellipse routines
%   are not part of the C++ library.)  All the functions are vectorized and
%   so offer speeds comparable to compiled C++ code when operating on
%   arrays.  For information on any function, use MATLAB's help command.
%   More extensive information is available in the documentation for the
%   GeographicLib, which is available at
%
%       http://geographiclib.sf.net/html
%
%   Some common features of these functions:
%     * Angles (latitude, longitude, azimuth, meridian convergence) are
%       measured in degrees.
%     * Distances are measured in meters, areas in meters^2.
%     * Latitudes must lie in [-90,90] and longitudes and azimuths in
%       [-540,540).  However most routines don't check that this condition
%       holds.  (Exceptions are the grid system and geoid functions.  These
%       return NaNs for invalid inputs.)
%     * The ellipsoid is specified as [a, e], where a = equatorial radius
%       and e = eccentricity.  The eccentricity can be pure imaginary to
%       denote a prolate ellipsoid.
%     * Keep abs(e) < 0.2 (i.e., abs(f) <= 1/50) for full double precision
%       accuracy.
%
%   There is some overlap between this toolbox and MATLAB's Mapping
%   Toolbox.  However, this toolbox offers
%     * Better accuracy.
%     * Treatment of oblate and prolate ellipsoid
%     * Guaranteed convergence for geoddistance.
%     * Calculation of area and differential properties of geodesics.
%     * Ellipsoidal versions of the equidistant azimuthal and gnomonic
%       projections.
%
% Geodesics
%   geoddistance     - Distance between points on an ellipsoid
%   geodreckon       - Point at specified azimuth, range on an ellipsoid
%   geodarea         - Surface area of polygon on an ellipsoid
%
% Projections
%   tranmerc_fwd     - Forward transverse Mercator projection
%   tranmerc_inv     - Inverse transverse Mercator projection
%   polarst_fwd      - Forward polar stereographic projection
%   polarst_inv      - Inverse polar stereographic projection
%   eqdazim_fwd      - Forward azimuthal equidistant projection
%   eqdazim_inv      - Inverse azimuthal equidistant projection
%   cassini_fwd      - Forward Cassini-Soldner projection
%   cassini_inv      - Inverse Cassini-Soldner projection
%   gnomonic_fwd     - Forward ellipsoidal gnomonic projection
%   gnomonic_inv     - Inverse ellipsoidal gnomonic projection
%
% Grid systems
%   utmups_fwd       - Convert to UTM/UPS system
%   utmups_inv       - Convert from UTM/UPS system
%   mgrs_fwd         - Convert UTM/UPS coordinates to MGRS
%   mgrs_inv         - Convert MGRS to UTM/UPS coordinates
%
% Geoid lookup
%   geoid_height     - Compute the height of the geoid
%   geoid_load       - Load a geoid model
%
% Great ellipse
%   gedistance       - Great ellipse distance between points on an ellipsoid
%   gereckon         - Point along great ellipse at specified azimuth and range
%
% Geometric transformations
%   geocent_fwd      - Conversion from geographic to geocentric coordinates
%   geocent_inv      - Conversion from geocentric to geographic coordinates
%   loccart_fwd      - Convert geographic to local cartesian coordinates
%   loccart_inv      - Convert local cartesian to geographic coordinates
%
% Utility
%   defaultellipsoid - Return the WGS84 ellipsoid
%   ecc2flat         - Convert eccentricity to flattening
%   flat2ecc         - Convert flattening to eccentricity
%
% Documentation
%   geoddoc          - Geodesics on an ellipsoid of revolution
%   projinfo         - Projections for an ellipsoid
%   gedoc            - Great ellipses on an ellipsoid of revolution

% Copyright (c) Charles Karney (2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.
