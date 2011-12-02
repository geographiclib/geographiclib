/**
 * \file Gravity.cpp
 * \brief Command line utility for evaluating gravity fields
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Compile and link with
 *   g++ -g -O3 -I../include -I../man -o Gravity \
 *       Gravity.cpp \
 *       ../src/CircularEngine.cpp \
 *       ../src/DMS.cpp \
 *       ../src/Geocentric.cpp \
 *       ../src/GravityCircle.cpp \
 *       ../src/GravityModel.cpp \
 *       ../src/NormalGravity.cpp \
 *       ../src/SphericalEngine.cpp \
 *       ../src/Utility.cpp
 *
 * See the <a href="Gravity.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#include "Gravity.usage"

int main(int argc, char* argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    bool verbose = false;
    std::string dir;
    std::string model = GravityModel::DefaultGravityName();
    std::string istring, ifile, ofile;
    real lat = 0, h = 0;
    bool circle = false;
    int prec = 5;
    enum {
      GRAVITY = 0,
      DISTURBANCE = 1,
      ANOMALY = 2,
      UNDULATION = 3,
    };
    unsigned mode = GRAVITY;
    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        model = argv[m];
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        dir = argv[m];
      } else if (arg == "-G")
        mode = GRAVITY;
      else if (arg == "-D")
        mode = DISTURBANCE;
      else if (arg == "-A")
        mode = ANOMALY;
      else if (arg == "-H")
        mode = UNDULATION;
      else if (arg == "-c") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          DMS::flag ind;
          lat = DMS::Decode(std::string(argv[++m]), ind);
          if (ind == DMS::LONGITUDE)
            throw GeographicErr("Bad hemisphere letter on latitude");
          if (!(std::abs(lat) <= 90))
            throw GeographicErr("Latitude not in [-90d, 90d]");
          h =  Utility::num<real>(std::string(argv[++m]));
          circle = true;
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::num<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-v")
        verbose = true;
      else if (arg == "--input-string") {
        if (++m == argc) return usage(1, true);
        istring = argv[m];
      } else if (arg == "--input-file") {
        if (++m == argc) return usage(1, true);
        ifile = argv[m];
      } else if (arg == "--output-file") {
        if (++m == argc) return usage(1, true);
        ofile = argv[m];
      } else if (arg == "--version") {
        std::cout
          << argv[0]
          << ": $Id$\n"
          << "GeographicLib version " << GEOGRAPHICLIB_VERSION_STRING << "\n";
        return 0;
      } else {
        int retval = usage(!(arg == "-h" || arg == "--help"), arg != "--help");
        if (arg == "-h")
          std::cout<< "\nDefault gravity path = \""
                   << GravityModel::DefaultGravityPath()
                   << "\"\nDefault gravity name = \""
                   << GravityModel::DefaultGravityName()
                   << "\"\n";
        return retval;
      }
    }

    if (!ifile.empty() && !istring.empty()) {
      std::cerr << "Cannot specify --input-string and --input-file together\n";
      return 1;
    }
    if (ifile == "-") ifile.clear();
    std::ifstream infile;
    std::istringstream instring;
    if (!ifile.empty()) {
      infile.open(ifile.c_str());
      if (!infile.is_open()) {
        std::cerr << "Cannot open " << ifile << " for reading\n";
        return 1;
      }
    } else if (!istring.empty()) {
      std::string::size_type m = 0;
      while (true) {
        m = istring.find(';', m);
        if (m == std::string::npos)
          break;
        istring[m] = '\n';
      }
      instring.str(istring);
    }
    std::istream* input = !ifile.empty() ? &infile :
      (!istring.empty() ? &instring : &std::cin);

    std::ofstream outfile;
    if (ofile == "-") ofile.clear();
    if (!ofile.empty()) {
      outfile.open(ofile.c_str());
      if (!outfile.is_open()) {
        std::cerr << "Cannot open " << ofile << " for writing\n";
        return 1;
      }
    }
    std::ostream* output = !ofile.empty() ? &outfile : &std::cout;

    prec = std::min(12, std::max(0, prec));
    int retval = 0;
    try {
      const GravityModel g(model, dir);
      if (circle) {
        if (!Math::isfinite<real>(h))
          throw GeographicErr("Bad height");
        else if (mode == UNDULATION && h != 0)
          throw GeographicErr("Height should be zero for geoid undulations");
      }
      if (verbose) {
        std::cerr << "Gravity file: " << g.GravityFile()      << "\n"
                  << "Name: "         << g.GravityModelName() << "\n"
                  << "Description: "  << g.Description()      << "\n"
                  << "Date & Time: "  << g.DateTime()         << "\n";
      }
      unsigned mask = (mode == GRAVITY ? GravityModel::GRAVITY :
                       (mode == DISTURBANCE ? GravityModel::DISTURBANCE :
                        (mode == ANOMALY ? GravityModel::SPHERICAL_ANOMALY :
                         GravityModel::GEOID_HEIGHT))); // mode == UNDULATION
      const GravityCircle c(circle ? g.Circle(lat, h, mask) : GravityCircle());
      std::string s, stra, strb;
      while (std::getline(*input, s)) {
        try {
          std::istringstream str(s);
          real lon;
          if (circle) {
            if (!(str >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::flag ind;
            lon = DMS::Decode(strb, ind);
            if (ind == DMS::LATITUDE)
              throw GeographicErr("Bad hemisphere letter on " + strb);
            if (lon < -180 || lon > 360)
              throw GeographicErr("Longitude " + strb + "not in [-180d, 360d]");
          } else {
            if (!(str >> stra >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::DecodeLatLon(stra, strb, lat, lon);
            h = 0;
            if (!(str >> h))    // h is optional
              str.clear();
            if (mode == UNDULATION && h != 0)
                throw GeographicErr("Height must be zero for geoid heights");
          }
          if (str >> stra)
            throw GeographicErr("Extra junk in input: " + s);
          switch (mode) {
          case GRAVITY:
            {
              real gx, gy, gz;
              if (circle) {
                c.Gravity(lon, gx, gy, gz);
              } else {
                g.Gravity(lat, lon, h, gx, gy, gz);
              }
              *output << Utility::str<real>(gx, prec) << " "
                      << Utility::str<real>(gy, prec) << " "
                      << Utility::str<real>(gz, prec) << "\n";
            }
            break;
          case DISTURBANCE:
            {
              real deltax, deltay, deltaz;
              if (circle) {
                c.Disturbance(lon, deltax, deltay, deltaz);
              } else {
                g.Disturbance(lat, lon, h, deltax, deltay, deltaz);
              }
              *output << Utility::str<real>(deltax, prec) << " "
                      << Utility::str<real>(deltay, prec) << " "
                      << Utility::str<real>(deltaz, prec) << "\n";
            }
            break;
          case ANOMALY:
            {
              real Dg01, xi, eta;
              if (circle)
                c.SphericalAnomaly(lon, Dg01, xi, eta);
              else
                g.SphericalAnomaly(lat, lon, h, Dg01, xi, eta);
              Dg01 *= 1e5;      // Convert to mGals
              xi *= 3600;       // Convert to arcsecs
              eta *= 3600;
              *output << Utility::str<real>(Dg01, prec) << " "
                      << Utility::str<real>(xi, prec) << " "
                      << Utility::str<real>(eta, prec) << "\n";
            }
            break;
          case UNDULATION:
          default:
            {
              real N = circle ? c.GeoidHeight(lon) : g.GeoidHeight(lat, lon);
              *output << Utility::str<real>(N, prec) << "\n";
            }
            break;
          }
        }
        catch (const std::exception& e) {
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading " << model << ": " << e.what() << "\n";
      retval = 1;
    }
    return retval;
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}
