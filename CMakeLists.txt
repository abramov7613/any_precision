cmake_minimum_required(VERSION 3.20)

project(z3mp LANGUAGES CXX)

add_library(z3mp STATIC 
  src/util/mpn.h
  src/util/mpn.cpp
  src/util/mpz.h
  src/util/mpz.cpp
  src/util/hash.h
  src/util/hash.cpp
  src/util/bit_util.h
  src/util/bit_util.cpp
)

target_compile_features(z3mp PUBLIC cxx_std_20)

target_compile_definitions(z3mp PUBLIC _MP_INTERNAL)

target_include_directories(z3mp
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

enable_testing()
add_subdirectory(tests)
