# CMake find_package tests

The purpose of this directory is to test the CONFIG mode CMake export can be found and used successfully. 


## Instructions


On Linux:

```bash
cd /path/to/GeoGraphicLib
cmake -S . -B build
cmake --build build
cmake --install build --prefix install
cd tests/cmake_find_pkg
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
./build/main_static
./build/main_shared
```
