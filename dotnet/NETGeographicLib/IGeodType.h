/**
 * \file IGeodType.hpp
 * \brief Header for GeographicLib::IGeodType interface
 *
 * NETGeographicLib is copyright (c) Scott heiman.
 * GeographicLibb is copyright (c) Charles Karney (2010-2014) 
 * <charles@karney.com> and licensed under the MIT/X11 License.  For 
  * more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

namespace NETGeographicLib
{
    /**
     * \brief The interface class for geoid types.
     *
     * This interface is used by the PolygonAreaT class.  It is inherited by
     * the Geodesoc and GeodesicExact classes.
     *************************************************************************/
    public interface class IGeodType
    {
        /** \name Inspector functions.
         **********************************************************************/
        ///@{

        /**
         * @return \e a the equatorial radius of the ellipsoid (meters).  This is
         *   the value used in the constructor.
         **********************************************************************/
        property double MajorRadius { double get(); }

        /**
         * @return \e f the  flattening of the ellipsoid.  This is the
         *   value used in the constructor.
         **********************************************************************/
        property double Flattening { double get(); }

        /**
         * @return total area of ellipsoid in meters<sup>2</sup>.  The area of a
         *   polygon encircling a pole can be found by adding
         *   GeodesicExact::EllipsoidArea()/2 to the sum of \e S12 for each side of
         *   the polygon.
         **********************************************************************/
        property double EllipsoidArea { double get(); }
        ///@}

        /**
         * The general direct geodesic calculation.  GeodesicExact::Direct and
         * GeodesicExact::ArcDirect are defined in terms of this function.
         *
         * @param[in] lat1 latitude of point 1 (degrees).
         * @param[in] lon1 longitude of point 1 (degrees).
         * @param[in] azi1 azimuth at point 1 (degrees).
         * @param[in] arcmode boolean flag determining the meaning of the second
         *   parameter.
         * @param[in] s12_a12 if \e arcmode is false, this is the distance between
         *   point 1 and point 2 (meters); otherwise it is the arc length between
         *   point 1 and point 2 (degrees); it can be signed.
         * @param[in] outmask a bitor'ed combination of  NETGeographicLib::Mask values
         *   specifying which of the following parameters should be set.
         * @param[out] lat2 latitude of point 2 (degrees).
         * @param[out] lon2 longitude of point 2 (degrees).
         * @param[out] azi2 (forward) azimuth at point 2 (degrees).
         * @param[out] s12 distance between point 1 and point 2 (meters).
         * @param[out] m12 reduced length of geodesic (meters).
         * @param[out] M12 geodesic scale of point 2 relative to point 1
         *   (dimensionless).
         * @param[out] M21 geodesic scale of point 1 relative to point 2
         *   (dimensionless).
         * @param[out] S12 area under the geodesic (meters<sup>2</sup>).
         * @return \e a12 arc length of between point 1 and point 2 (degrees).
         *
         * The  NETGeographicLib::Mask values possible for \e outmask are
         * - \e outmask |= NETGeographicLib::Mask::LATITUDE for the latitude \e lat2;
         * - \e outmask |= NETGeographicLib::Mask::LONGITUDE for the latitude \e lon2;
         * - \e outmask |= NETGeographicLib::Mask::AZIMUTH for the latitude \e azi2;
         * - \e outmask |= NETGeographicLib::Mask::DISTANCE for the distance \e s12;
         * - \e outmask |= NETGeographicLib::Mask::REDUCEDLENGTH for the reduced length \e
         *   m12;
         * - \e outmask |= NETGeographicLib::Mask::GEODESICSCALE for the geodesic scales \e
         *   M12 and \e M21;
         * - \e outmask |= NETGeographicLib::Mask::AREA for the area \e S12;
         * - \e outmask |= NETGeographicLib::Mask::ALL for all of the above.
         * .
         * The function value \e a12 is always computed and returned and this
         * equals \e s12_a12 is \e arcmode is true.  If \e outmask includes
         * GeodesicExact::DISTANCE and \e arcmode is false, then \e s12 = \e
         * s12_a12.  It is not necessary to include  NETGeographicLib::Mask::DISTANCE_IN in
         * \e outmask; this is automatically included is \e arcmode is false.
         **********************************************************************/
        double GenDirect(double lat1, double lon1, double azi1,
                        bool arcmode, double s12_a12, NETGeographicLib::Mask outmask,
                        [System::Runtime::InteropServices::Out] double% lat2,
                        [System::Runtime::InteropServices::Out] double% lon2,
                        [System::Runtime::InteropServices::Out] double% azi2,
                        [System::Runtime::InteropServices::Out] double% s12,
                        [System::Runtime::InteropServices::Out] double% m12,
                        [System::Runtime::InteropServices::Out] double% M12,
                        [System::Runtime::InteropServices::Out] double% M21,
                        [System::Runtime::InteropServices::Out] double% S12);
    };
};