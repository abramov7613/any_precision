find_package(Boost 1.82 CONFIG QUIET)
if(NOT Boost_FOUND)
  message(STATUS "Boost not found or version too old, fetching Boost::multiprecision")
  # Fallback to FetchContent
  include(FetchContent)
  FetchContent_Declare(boost_multiprecision
      GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
      GIT_TAG        boost-1.92.0
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
)

target_compile_features(mpz_tests PRIVATE cxx_std_20)

target_link_libraries(mpz_tests PRIVATE
    z3mp
    GTest::gtest_main
    Boost::multiprecision
)

gtest_discover_tests(mpz_tests)
