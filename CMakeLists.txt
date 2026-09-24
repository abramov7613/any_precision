cmake_minimum_required(VERSION 3.31)

project(any_precision LANGUAGES CXX)

find_package(Boost 1.82 CONFIG QUIET)
if(NOT Boost_FOUND)
  message(STATUS "Boost not found or version too old, fetching Boost::multiprecision")
  include(FetchContent)
  FetchContent_Declare(boost_multiprecision
      URL https://github.com/boostorg/multiprecision/archive/refs/tags/Boost_1_92_0.tar.gz
      URL_HASH SHA256=d53a907779b1ba88c80c15953c24106bf92b1b4731678441a8655e527c9d3df1
  )
  FetchContent_MakeAvailable(boost_multiprecision)
  if(NOT TARGET boost_multiprecision)
      message(FATAL_ERROR "Boost::multiprecision not found")
  endif()
endif()

add_library(${PROJECT_NAME} STATIC
  src/mpn.h
  src/mpz.h
  src/mpz.cpp
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(Boost_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
else()
  target_link_libraries(${PROJECT_NAME} PRIVATE boost_multiprecision)
endif()


enable_testing()
add_subdirectory(tests)
