cmake_minimum_required(VERSION 3.21)

project(ez3i) # extracted_z3_integer

option(BUILD_SHARED_LIBS "Build shared library (ON) or static (OFF)" OFF)
option(EZ3I_BUILD_TESTS "Build unit tests" ${PROJECT_IS_TOP_LEVEL})
option(EZ3I_INSTALL "Generate install targets" ${PROJECT_IS_TOP_LEVEL})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

add_library(${PROJECT_NAME} src/mpz.cpp)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
if(BUILD_SHARED_LIBS)
    target_compile_definitions(${PROJECT_NAME} PRIVATE EZ3I_EXPORTS)
else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC EZ3I_STATIC)
endif()

if(EZ3I_BUILD_TESTS AND PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(tests)
endif()

if(EZ3I_INSTALL)
    include(GNUInstallDirs)

    install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endif()
