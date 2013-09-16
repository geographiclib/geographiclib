#pragma once
/**
 * \file NETGeographicLib/NETGeographicLib.h
 * \brief Header for NETGeographicLib::NETGeographicLib objects
 *
 * NETGeographicLib is copyright (c) Scott Heiman (2013)
 * GeographicLib is Copyright (c) Charles Karney (2010-2012)
 * <charles@karney.com> and licensed under the MIT/X11 License.
 * For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/
#include <string>

using namespace System;

namespace NETGeographicLib
{
    enum class captype {
      CAP_NONE = 0U,
      CAP_C1   = 1U<<0,
      CAP_C1p  = 1U<<1,
      CAP_C2   = 1U<<2,
      CAP_C3   = 1U<<3,
      CAP_C4   = 1U<<4,
      CAP_ALL  = 0x1FU,
      OUT_ALL  = 0x7F80U,
    };
    /**
    * The version string.
     **********************************************************************/
    public ref class VersionInfo
    {
        VersionInfo() {}
    public:
        static System::String^ GetString();
        static int MajorVersion();
        static int MinorVersion();
        static int Patch();
    };

    /**
     * Bit masks for what calculations to do.  These masks do double duty.
     * They signify to the GeodesicLine::GeodesicLine constructor and to
     * Geodesic::Line what capabilities should be included in the GeodesicLine
     * object.  They also specify which results to return in the general
     * routines Geodesic::GenDirect and Geodesic::GenInverse routines.
     **********************************************************************/
    public enum class Mask {
      /**
       * No capabilities, no output.
       * @hideinitializer
       **********************************************************************/
      NONE          = 0U,
      /**
       * Calculate latitude \e lat2.  (It's not necessary to include this as a
       * capability to GeodesicLine because this is included by default.)
       * @hideinitializer
       **********************************************************************/
      LATITUDE      = 1U<<7  | captype::CAP_NONE,
      /**
       * Calculate longitude \e lon2.
       * @hideinitializer
       **********************************************************************/
      LONGITUDE     = 1U<<8  | captype::CAP_C3,
      /**
       * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
       * include this as a capability to GeodesicLine because this is included
       * by default.)
       * @hideinitializer
       **********************************************************************/
      AZIMUTH       = 1U<<9  | captype::CAP_NONE,
      /**
       * Calculate distance \e s12.
       * @hideinitializer
       **********************************************************************/
      DISTANCE      = 1U<<10 | captype::CAP_C1,
      /**
       * Allow distance \e s12 to be used as input in the direct geodesic
       * problem.
       * @hideinitializer
       **********************************************************************/
      DISTANCE_IN   = 1U<<11 | captype::CAP_C1 | captype::CAP_C1p,
      /**
       * Calculate reduced length \e m12.
       * @hideinitializer
       **********************************************************************/
      REDUCEDLENGTH = 1U<<12 | captype::CAP_C1 | captype::CAP_C2,
      /**
       * Calculate geodesic scales \e M12 and \e M21.
       * @hideinitializer
       **********************************************************************/
      GEODESICSCALE = 1U<<13 | captype::CAP_C1 | captype::CAP_C2,
      /**
       * Calculate area \e S12.
       * @hideinitializer
       **********************************************************************/
      AREA          = 1U<<14 | captype::CAP_C4,
      /**
       * All capabilities, calculate everything.
       * @hideinitializer
       **********************************************************************/
      ALL           = captype::OUT_ALL| captype::CAP_ALL,
    };

    public ref class GeographicErr : public System::Exception
    {
    public:
        GeographicErr( const char* msg ) : System::Exception( gcnew System::String( msg ) ) {}
        GeographicErr( System::String^ msg ) : System::Exception( msg ) {}
    };

    ref class StringConvert
    {
        StringConvert() {}
    public:
        static std::string ManagedToUnmanaged( System::String^ s );
        static System::String^ UnmanagedToManaged( const std::string& s )
        {   return gcnew System::String( s.c_str() ); }
    };
}  // namespace NETGeographicLib
