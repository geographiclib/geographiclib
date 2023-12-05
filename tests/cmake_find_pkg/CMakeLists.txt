cmake_minimum_required (VERSION 3.5.0)
project(geographiclib_cmake_find_pkg)

find_package(GeographicLib CONFIG REQUIRED)

add_executable(main_static main.cpp)
target_link_libraries(main_static PRIVATE GeographicLib::GeographicLib)

add_executable(main_shared main.cpp)
target_link_libraries(main_shared PRIVATE GeographicLib::GeographicLib_SHARED)
