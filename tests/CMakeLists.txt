find_package(Boost 1.82 CONFIG QUIET)
if(NOT Boost_FOUND)
  message(STATUS "Boost not found or version too old, fetching Boost::multiprecision")
  # Fallback to FetchContent
  include(FetchContent)
  FetchContent_Declare(boost_multiprecision
      URL https://github.com/boostorg/multiprecision/archive/refs/tags/Boost_1_92_0.tar.gz
      URL_HASH SHA256=d53a907779b1ba88c80c15953c24106bf92b1b4731678441a8655e527c9d3df1
  )
  FetchContent_MakeAvailable(boost_multiprecision)
  if(NOT TARGET boost_multiprecision)
      message(FATAL_ERROR "Boost::multiprecision not found")
  endif()
  if(NOT TARGET Boost::multiprecision)
      add_library(Boost::multiprecision ALIAS boost_multiprecision)
  endif()
endif()

find_package(GTest 1.17 QUIET)
if(NOT GTest_FOUND)
  message(STATUS "Google Test not found or version too old, fetching GTest")
  include(FetchContent)
  FetchContent_Declare(googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        v1.17.0
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  if(NOT TARGET GTest::gtest_main)
      message(FATAL_ERROR "Google Test not found")
  endif()
endif()
include(GoogleTest)

add_executable(mpz_tests
    test_mpz.cpp
    test_mpz1.cpp
)

target_compile_features(mpz_tests PRIVATE cxx_std_20)

target_link_libraries(mpz_tests PRIVATE
    z3mp
    GTest::gtest_main
    Boost::multiprecision
)

gtest_discover_tests(mpz_tests)
