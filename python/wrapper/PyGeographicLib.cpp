#include <boost/python.hpp>
#include <GeographicLib/Geodesic.hpp>

using namespace boost::python;
using namespace GeographicLib;

tuple Geodesic_Inverse(Geodesic& geo,
                       double lat1, double lon1, double lat2, double lon2) {
  double s12, azi1, azi2, m12, M12, M21, S12,
    a12 = geo.Inverse(lat1, lon1, lat2, lon2,
                      s12, azi1, azi2, m12, M12, M21, S12);
  return make_tuple(azi1, azi2, s12, a12, m12, M12, M21, S12);
}

BOOST_PYTHON_MODULE(PyGeographicLib) {

  class_<Constants>("Constants", no_init)
    .def("WGS84_a", &Constants::WGS84_a<double>,
         "The equatorial radius of the WGS84 ellipsoid")
    .staticmethod("WGS84_a")
    .def("WGS84_f", &Constants::WGS84_f<double>,
         "The flattening of the WGS84 ellipsoid")
    .staticmethod("WGS84_f")
    ;

  class_<Geodesic>("Geodesic", init<double, double>())
    .def("Inverse", &Geodesic_Inverse,
         "Solve inverse geodesic problem:\n\
    input: lat1, lon1, lat2, lon2\n\
    output: (azi1, azi2, s12, a12, m12, M12, M21, S12)")
    .def("MajorRadius", &Geodesic::MajorRadius,
         "The equatorial radius of the ellipsoid")
    .def("Flattening", &Geodesic::Flattening,
         "The flattening of the ellipsoid")
    ;

}
