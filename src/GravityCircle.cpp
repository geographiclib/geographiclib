/**
 * \file GravityCircle.cpp
 * \brief Implementation for GeographicLib::GravityCircle class
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/GravityCircle.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <GeographicLib/Geocentric.hpp>

#define GEOGRAPHICLIB_GRAVITYCIRCLE_CPP "$Id$"

RCSID_DECL(GEOGRAPHICLIB_GRAVITYCIRCLE_CPP)
RCSID_DECL(GEOGRAPHICLIB_GRAVITYCIRCLE_HPP)

#define GRAVITY_DEFAULT_PATH "/home/ckarney/geographiclib/gravity"

namespace GeographicLib {

  using namespace std;

  Math::real GravityCircle::GeoidHeight(real lon) const throw() {
    real clam, slam;
    CircularEngine::cossin(lon, clam, slam);
    real T = _disturbing(clam, slam);
    bool correct = false;
    T = (T / _amodel - (correct ? _dzonal0 : 0) * _invR) * _GMmodel;
    real correction = _corrmult * _correction(clam, slam);
    return T/_gamma0 + correction;
  }

} // namespace GeographicLib
