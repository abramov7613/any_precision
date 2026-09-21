cmake_minimum_required(VERSION 3.20)

project(any_precision LANGUAGES CXX)

add_library(${PROJECT_NAME} STATIC 
  src/util/mpn.h
  src/util/mpn.cpp
  src/util/mpz.h
  src/util/mpz.cpp
  src/util/hash.h
  src/util/hash.cpp
  src/util/bit_util.h
  src/util/bit_util.cpp
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

target_compile_definitions(${PROJECT_NAME} PUBLIC _MP_INTERNAL)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

enable_testing()
add_subdirectory(tests)
