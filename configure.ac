dnl
dnl Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
dnl $Id$

AC_INIT([GeographicLib],[1.12],[charles@karney.com])
AC_CANONICAL_SYSTEM
AC_PREREQ(2.61)
AC_CONFIG_SRCDIR(src/Geodesic.cpp)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

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
LT_CURRENT=8
LT_REVISION=0
LT_AGE=3
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

AC_ARG_PROGRAM
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_CXX
AC_PROG_LIBTOOL

# Checks for header files.
AC_CHECK_HEADERS([float.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_LONG_DOUBLE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([strtol])
AC_SEARCH_LIBS([floor],[m])
AC_SEARCH_LIBS([pow],[m])
AC_SEARCH_LIBS([sqrt],[m])

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
])
AC_OUTPUT
