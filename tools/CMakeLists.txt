# $Id$

if (MAINTAINER)
  set (MANDIR ${PROJECT_BINARY_DIR}/man)
else (MAINTAINER)
  set (MANDIR ${PROJECT_SOURCE_DIR}/man)
endif (MAINTAINER)

include_directories (${MANDIR})

add_custom_target (tools ALL)
foreach (TOOL ${TOOLS})

  add_executable (${TOOL} ${TOOL}.cpp)
  if (MAINTAINER)
    add_dependencies (${TOOL} usage)
  endif (MAINTAINER)
  add_dependencies (tools ${TOOL})

  set_source_files_properties (${TOOL}.cpp PROPERTIES
    OBJECT_DEPENDS ${MANDIR}/${TOOL}.usage)

  target_link_libraries (${TOOL} Geographic)

endforeach (TOOL)

if (NOT WIN32)
  set_target_properties (${TOOLS} PROPERTIES
    INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
endif (NOT WIN32)

install (TARGETS ${TOOLS} EXPORT GeographicLibDepends DESTINATION bin)

if (NOT WIN32)
  set (GEOID_SCRIPT geographiclib-get-geoids)
  add_custom_target (get-geoids-script ALL DEPENDS ${GEOID_SCRIPT})
  add_custom_command (OUTPUT ${GEOID_SCRIPT}
    COMMAND
      sed -e "s%@GEOID_DEFAULT_PATH@%${GEOGRAPHICLIB_GEOID_PATH}%"
	${CMAKE_CURRENT_SOURCE_DIR}/${GEOID_SCRIPT}.sh >
	${GEOID_SCRIPT}
    COMMAND chmod +x ${GEOID_SCRIPT}
    DEPENDS ${GEOID_SCRIPT}.sh)
  install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${GEOID_SCRIPT}
    DESTINATION sbin)

  set (GRAVITY_SCRIPT geographiclib-get-gravity)
  add_custom_target (get-gravity-script ALL DEPENDS ${GRAVITY_SCRIPT})
  add_custom_command (OUTPUT ${GRAVITY_SCRIPT}
    COMMAND
      sed -e "s%@GRAVITY_DEFAULT_PATH@%${GEOGRAPHICLIB_GRAVITY_PATH}%"
	${CMAKE_CURRENT_SOURCE_DIR}/${GRAVITY_SCRIPT}.sh >
	${GRAVITY_SCRIPT}
    COMMAND chmod +x ${GRAVITY_SCRIPT}
    DEPENDS ${GRAVITY_SCRIPT}.sh)
  install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${GRAVITY_SCRIPT}
    DESTINATION sbin)

  set (MAGNETIC_SCRIPT geographiclib-get-magnetic)
  add_custom_target (get-magnetic-script ALL DEPENDS ${MAGNETIC_SCRIPT})
  add_custom_command (OUTPUT ${MAGNETIC_SCRIPT}
    COMMAND
      sed -e "s%@MAGNETIC_DEFAULT_PATH@%${GEOGRAPHICLIB_MAGNETIC_PATH}%"
	${CMAKE_CURRENT_SOURCE_DIR}/${MAGNETIC_SCRIPT}.sh >
	${MAGNETIC_SCRIPT}
    COMMAND chmod +x ${MAGNETIC_SCRIPT}
    DEPENDS ${MAGNETIC_SCRIPT}.sh)
  install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${MAGNETIC_SCRIPT}
    DESTINATION sbin)
endif (NOT WIN32)

if (WIN32 AND GEOGRAPHIC_SHARED_LIB)
  get_target_property (GEOGRAPHIC_LIBRARY_FILE Geographic
    "LOCATION_${CMAKE_CFG_INTDIR}")
  add_custom_command (TARGET tools POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
    copy ${GEOGRAPHIC_LIBRARY_FILE} ${CMAKE_CFG_INTDIR}/)
endif (WIN32 AND GEOGRAPHIC_SHARED_LIB)

enable_testing ()

add_test (GeoConvert0 GeoConvert -p -3 -m --input-string "33.3 44.4")
set_tests_properties (GeoConvert0
  PROPERTIES PASS_REGULAR_EXPRESSION "38SMB4484")
add_test (GeoConvert1 GeoConvert -d --input-string "38smb")
set_tests_properties (GeoConvert1
  PROPERTIES PASS_REGULAR_EXPRESSION "32d59'14\\.1\"N 044d27'53\\.4\"E")
add_test (GeoConvert2 GeoConvert -p -2 --input-string "30d30'30\" 30.50833")
set_tests_properties (GeoConvert2
  PROPERTIES PASS_REGULAR_EXPRESSION "30\\.508 30\\.508")
add_test (GeoConvert3 GeoConvert --junk)
set_tests_properties (GeoConvert3 PROPERTIES WILL_FAIL ON)
add_test (GeoConvert4 GeoConvert --input-string garbage)
set_tests_properties (GeoConvert4 PROPERTIES WILL_FAIL ON)
# Check fix for DMS::Decode bug fixed on 2011-03-22
add_test (GeoConvert5 GeoConvert --input-string "5d. 0")
set_tests_properties (GeoConvert5 PROPERTIES WILL_FAIL ON)

add_test (Geod0 Geod -i -p 0 --input-string "40.6 -73.8 49d01'N 2d33'E")
set_tests_properties (Geod0
  PROPERTIES PASS_REGULAR_EXPRESSION "53\\.47022 111\\.59367 5853226")
add_test (Geod1
  Geod -p 0 --input-string "40d38'23\"N 073d46'44\"W 53d30' 5850e3")
set_tests_properties (Geod1
  PROPERTIES PASS_REGULAR_EXPRESSION "49\\.01467 2\\.56106 111\\.62947")
# Check fix for antipodal prolate bug found 2010-09-04
add_test (Geod2
  Geod -i -p 0 -e 6.4e6 -1/150 --input-string "0.07476 0 -0.07476 180")
set_tests_properties (Geod2
  PROPERTIES PASS_REGULAR_EXPRESSION "90\\.00078 90\\.00078 20106193")
# Another check for similar bug
add_test (Geod3
  Geod -i -p 0 -e 6.4e6 -1/150 --input-string "0.1 0 -0.1 180")
set_tests_properties (Geod3
  PROPERTIES PASS_REGULAR_EXPRESSION "90\\.00105 90\\.00105 20106193")
# Check fix for short line bug found 2010-05-21
add_test (Geod4
  Geod -i --input-string "36.493349428792 0 36.49334942879201 .0000008")
set_tests_properties (Geod4 PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 0\\.072")
# Check fix for point2=pole bug found 2010-05-03 (but only with long double)
add_test (Geod5
  Geod -p 0 --input-string "0.01777745589997 30 0 10e6")
set_tests_properties (Geod5
  PROPERTIES PASS_REGULAR_EXPRESSION
  "90\\.00000 -150\\.00000 -180\\.00000;90\\.00000 30\\.00000 0\\.00000")

# Check fix for volatile sbet12a bug found 2011-06-25 (gcc 4.4.4 x86 -O3)
add_test (Geod6 Geod -i --input-string
  "88.202499451857 0 -88.202499451857 179.981022032992859592")
add_test (Geod7 Geod -i --input-string
  "89.262080389218 0 -89.262080389218 179.992207982775375662")
add_test (Geod8 Geod -i --input-string
  "89.333123580033 0 -89.333123580032997687 179.99295812360148422")
set_tests_properties (Geod6
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 20003898.214")
set_tests_properties (Geod7
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 20003925.854")
set_tests_properties (Geod8
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 20003926.881")

# Check fix for volatile x bug found 2011-06-25 (gcc 4.4.4 x86 -O3)
add_test (Geod9 Geod -i --input-string
  "56.320923501171 0 -56.320923501171 179.664747671772880215")
set_tests_properties (Geod9
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 19993558.287")

# Check fix for adjust tol1_ bug found 2011-06-25 (Visual Studio 10 rel + debug)
add_test (Geod10 Geod -i --input-string
  "52.784459512564 0 -52.784459512563990912 179.634407464943777557")
set_tests_properties (Geod10
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 19991596.095")

# Check fix for bet2 = -bet1 bug found 2011-06-25 (Visual Studio 10 rel + debug)
add_test (Geod11 Geod -i --input-string
  "48.522876735459 0 -48.52287673545898293 179.599720456223079643")
set_tests_properties (Geod11
  PROPERTIES PASS_REGULAR_EXPRESSION ".* .* 19989144.774")

# Check fix for pole-encircling bug found 2011-03-16
add_test (Planimeter0 Planimeter --input-string "89 0;89 90;89 180;89 270")
add_test (Planimeter1
  Planimeter -r --input-string "-89 0;-89 90;-89 180;-89 270")
add_test (Planimeter2 Planimeter --input-string "0 -1;-1 0;0 1;1 0")
add_test (Planimeter3 Planimeter --input-string "90 0; 0 0; 0 90")
add_test (Planimeter4 Planimeter -l --input-string "90 0; 0 0; 0 90")
set_tests_properties (Planimeter0
  PROPERTIES PASS_REGULAR_EXPRESSION
  "4 631819\\.8745[0-9]+ 2495230567[78]\\.[0-9]+")
set_tests_properties (Planimeter1
  PROPERTIES PASS_REGULAR_EXPRESSION
  "4 631819\\.8745[0-9]+ 2495230567[78]\\.[0-9]+")
set_tests_properties (Planimeter2
  PROPERTIES PASS_REGULAR_EXPRESSION "4 627598\\.2731[0-9]+ 24619419146.[0-9]+")
set_tests_properties (Planimeter3
  PROPERTIES PASS_REGULAR_EXPRESSION
  "3 30022685\\.[0-9]+ 63758202715511\\.[0-9]+")
set_tests_properties (Planimeter4
  PROPERTIES PASS_REGULAR_EXPRESSION "3 20020719\\.[0-9]+")
# Check fix for Planimeter pole crossing bug found 2011-06-24
add_test (Planimeter5 Planimeter --input-string "89,0.1;89,90.1;89,-179.9")
set_tests_properties (Planimeter5
  PROPERTIES PASS_REGULAR_EXPRESSION
  "3 539297\\.[0-9]+ 1247615283[89]\\.[0-9]+")

# Check fix for AlbersEqualArea::Reverse bug found 2011-05-01
add_test (ConicProj0
  ConicProj -a 40d58 39d56 -l 77d45W -r --input-string "220e3 -52e3")
set_tests_properties (ConicProj0
  PROPERTIES PASS_REGULAR_EXPRESSION
  "39\\.95[0-9]+ -75\\.17[0-9]+ 1\\.67[0-9]+ 0\\.99[0-9]+")

add_test (CartConvert0
  CartConvert -e 6.4e6 1/100 -r --input-string "10e3 0 1e3")
add_test (CartConvert1
  CartConvert -e 6.4e6 -1/100 -r --input-string "1e3 0 10e3")
set_tests_properties (CartConvert0
  PROPERTIES PASS_REGULAR_EXPRESSION
  "85\\.57[0-9]+ 0\\.0[0]+ -6334614\\.[0-9]+")
set_tests_properties (CartConvert1
  PROPERTIES PASS_REGULAR_EXPRESSION
  "4\\.42[0-9]+ 0\\.0[0]+ -6398614\\.[0-9]+")

if (EXISTS ${GEOGRAPHICLIB_GEOID_PATH}/egm96-5.pgm)
  # Check fix for single-cell cache bug found 2010-11-23
  add_test (GeoidEval0 GeoidEval --input-string "0d1 0d1;0d4 0d4")
  set_tests_properties (GeoidEval0
    PROPERTIES PASS_REGULAR_EXPRESSION "^17\\.1602 .*17\\.1511 " )
endif (EXISTS ${GEOGRAPHICLIB_GEOID_PATH}/egm96-5.pgm)
