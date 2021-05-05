cmake_minimum_required(VERSION 3.2)

project(
  baaraan
  LANGUAGES CXX
  VERSION 0.1.0)
set(${PROJECT_NAME}_VERSION_NAME "Appa")

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_CAP)

set(CMAKE_CXX_STANDARD 17)
set(BAARAAN_DIR ${CMAKE_PROGRAM_PATH})

include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH ${ORIGIN})
endif()

set(ENABLE_TESTS OFF)

file(GLOB CPP_FILES *.cpp)

find_package(
  Boost
  COMPONENTS unit_test_framework
  REQUIRED)
find_package(Armadillo REQUIRED)

include_directories(include)
include_directories(include/dists)
include_directories(include/probs)
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_library(baaraan SHARED ${CPP_FILES})

target_link_libraries(baaraan ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})

if(${ENABLE_TESTS})
  enable_testing()
  add_subdirectory(tests)
endif()

write_basic_package_version_file(
  baaraanConfigVersion.cmake
  VERSION ${CMAKE_PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

# Install scripts
install(
  TARGETS baaraan
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          COMPONENT ${PROJECT_NAME_CAP}_Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT ${PROJECT_NAME_CAP}_Runtime
          NAMELINK_COMPONENT ${PROJECT_NAME_CAP}_Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT ${PROJECT_NAME_CAP}_Development)

# Exposing config variables
set(${PROJECT_NAME_CAP}_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME}: Includes Directory" FORCE)

set(${PROJECT_NAME_CAP}_DEFINITIONS
    ${Boost_DEFINITIONS} ${ARMADILLO_DEFINTIONS}
    CACHE INTERNAL "${PROJECT_NAME}: Definitions" FORCE)
