# Builds and runs the self tests on the host:
#   cmake -S host_test -B build-host -GNinja -DCMAKE_BUILD_TYPE=Release && cmake --build build-host && ./build-host/host_test
cmake_minimum_required(VERSION 3.14)
project(robust_weight_transfer_host_test CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(GC_EIGEN_LOCATION ${ROOT}/thirdparty/eigen CACHE STRING "" FORCE)
add_definitions(-DIGL_PARALLEL_FOR_FORCE_SERIAL)

add_subdirectory(${ROOT}/thirdparty/geometry_central geometry_central)

add_executable(host_test main.cpp)
target_include_directories(host_test SYSTEM PRIVATE
  ${ROOT}/thirdparty/eigen
  ${ROOT}/thirdparty/libigl/include
  ${ROOT}/thirdparty/witness-cpp/include)
target_link_libraries(host_test PRIVATE geometry-central)
