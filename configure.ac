dnl
dnl Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>

AC_INIT([GeographicLib],[1.42],[charles@karney.com])
AC_CANONICAL_SYSTEM
AC_PREREQ(2.61)
AC_CONFIG_SRCDIR(src/Geodesic.cpp)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

GEOGRAPHICLIB_VERSION_MAJOR=1
GEOGRAPHICLIB_VERSION_MINOR=42
GEOGRAPHICLIB_VERSION_PATCH=0
AC_DEFINE_UNQUOTED([GEOGRAPHICLIB_VERSION_MAJOR],
        [$GEOGRAPHICLIB_VERSION_MAJOR],[major version number])
AC_DEFINE_UNQUOTED([GEOGRAPHICLIB_VERSION_MINOR],
        [$GEOGRAPHICLIB_VERSION_MINOR],[minor version number])
AC_DEFINE_UNQUOTED([GEOGRAPHICLIB_VERSION_PATCH],
        [$GEOGRAPHICLIB_VERSION_PATCH],[patch number])
AC_SUBST(GEOGRAPHICLIB_VERSION_MAJOR)
AC_SUBST(GEOGRAPHICLIB_VERSION_MINOR)
AC_SUBST(GEOGRAPHICLIB_VERSION_PATCH)

dnl
dnl This directive is deprecated by someone, but I prefer to avoid
dnl running autotools if not required explicitly. The reason is
dnl the need to be in sync with autoconf/automake.
dnl
AM_MAINTAINER_MODE

AC_CONFIG_HEADERS(include/GeographicLib/Config-ac.h)

dnl Library code modified:              REVISION++
dnl Interfaces changed/added/removed:   CURRENT++ REVISION=0
dnl Interfaces added:                   AGE++
dnl Interfaces removed:                 AGE=0
LT_CURRENT=14
LT_REVISION=3
LT_AGE=0
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

AC_ARG_PROGRAM
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_LANG_CPLUSPLUS

# Checks for long double
AC_TYPE_LONG_DOUBLE

# Check endianness
AC_C_BIGENDIAN

# Check flags for C++11
AX_CHECK_COMPILE_FLAG([-std=c++11],
        [CXXFLAGS="$CXXFLAGS -std=c++11"],
        [AX_CHECK_COMPILE_FLAG([-std=c++0x],
                [CXXFLAGS="$CXXFLAGS -std=c++0x"])])
# Check for C++11 math functions
AC_TRY_COMPILE([#include <cmath>],
        [return int(std::hypot(3.0, 4.0) + std::expm1(0.5) +
                    std::log1p(2.0) + std::asinh(10.0) + std::atanh(0.8) +
                    std::fma(1.0, 2.0, 3.0) + std::cbrt(8.0)) +
                    std::isfinite(4.0) + std::isnan(0.0);],
        [CXXFLAGS="$CXXFLAGS -DGEOGRAPHICLIB_CXX11_MATH=1"],
        [CXXFLAGS="$CXXFLAGS -DGEOGRAPHICLIB_CXX11_MATH=0"])

# Check for doxygen.
AC_CHECK_PROGS([DOXYGEN], [doxygen])
AM_CONDITIONAL([HAVE_DOXYGEN], [test "$DOXYGEN"])

AC_CHECK_PROGS([POD2MAN], [pod2man])
AC_CHECK_PROGS([POD2HTML], [pod2html])
AC_CHECK_PROGS([COL], [col])
AM_CONDITIONAL([HAVE_PODPROGS], [test "$POD2MAN" -a  "$POD2HTML" -a "$COL"])

dnl
dnl Add here new file to be generated
dnl
AC_CONFIG_FILES([
Makefile
src/Makefile
include/Makefile
tools/Makefile
doc/Makefile
man/Makefile
matlab/Makefile
python/Makefile
cmake/Makefile
examples/Makefile
])
AC_OUTPUT
