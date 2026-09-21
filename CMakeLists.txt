cmake_minimum_required(VERSION 3.31)

project(any_precision LANGUAGES CXX)

add_library(${PROJECT_NAME} STATIC
  src/mpn.h
  src/mpn.cpp
  src/mpz.h
  src/mpz.cpp
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

enable_testing()
add_subdirectory(tests)
