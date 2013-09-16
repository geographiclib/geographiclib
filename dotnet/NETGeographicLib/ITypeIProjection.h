#pragma once

/*!
THIS INTERFACE IS DEPRECATED
*/
namespace NETGeographicLib
{
    public interface class ITypeIProjection
    {
    public:
        void Forward(double lon0, double lat, double lon,
                     [System::Runtime::InteropServices::Out] double% x,
                     [System::Runtime::InteropServices::Out] double% y,
                     [System::Runtime::InteropServices::Out] double% gamma,
                     [System::Runtime::InteropServices::Out] double% k);

        void Forward(double lon0, double lat, double lon,
                     [System::Runtime::InteropServices::Out] double% x,
                     [System::Runtime::InteropServices::Out] double% y);

        void Reverse(double lon0, double x, double y,
                     [System::Runtime::InteropServices::Out] double% lat,
                     [System::Runtime::InteropServices::Out] double% lon,
                     [System::Runtime::InteropServices::Out] double% gamma,
                     [System::Runtime::InteropServices::Out] double% k);

        void Reverse(double lon0, double x, double y,
                     [System::Runtime::InteropServices::Out] double% lat,
                     [System::Runtime::InteropServices::Out] double% lon);

        property double MajorRadius { double get(); }

        property double Flattening { double get(); }

        property double OriginLatitude { double get(); }

        property double CentralScale { double get(); }
    };
} // namespace NETGeographicLib
