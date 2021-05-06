cmake_minimum_required(VERSION 3.2)

project(
  baaraan
  LANGUAGES CXX
  VERSION 0.1.0)
set(${PROJECT_NAME}_VERSION_NAME "Appa")

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_CAP)

set(BAARAAN_VERSION ${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR}.
                 ${CMAKE_PROJECT_VERSION_PATCH})
message(STATUS ${BAARAAN_VERSION})

set(CMAKE_CXX_STANDARD 17)
set(BAARAAN_DIR ${CMAKE_PROGRAM_PATH})

include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH ${ORIGIN})
endif()

option(ENABLE_TESTS OFF)

file(GLOB CPP_FILES *.cpp)

find_package(Boost)
find_package(Armadillo REQUIRED)

include_directories(include)
include_directories(include/dists)
include_directories(include/probs)
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_library(baaraan INTERFACE)
add_library(baaraan::baaraan ALIAS baaraan)

target_link_libraries(baaraan INTERFACE ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})

if(${ENABLE_TESTS})
  enable_testing()
  add_subdirectory(tests)
endif()


# define project variables
set(BAARAAN_TARGET_NAME ${PROJECT_NAME})
set(BAARAAN_PACKAGE_NAME ${BAARAAN_TARGET_NAME})
set(BAARAAN_TARGETS_FILENAME "${BAARAAN_PACKAGE_NAME}Targets.cmake")
set(BAARAAN_CONFIG_FILENAME "${BAARAAN_PACKAGE_NAME}Config.cmake")
set(BAARAAN_CONFIGVERSION_FILENAME "${BAARAAN_PACKAGE_NAME}ConfigVersion.cmake")
set(BAARAAN_CONFIG_DESTINATION "cmake")
set(BAARAAN_INCLUDE_DESTINATION "include")

# generate a config and config version file for the package
configure_package_config_file("support/cmake/baaraan-config.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/${BAARAAN_CONFIG_FILENAME}"
                              INSTALL_DESTINATION ${BAARAAN_CONFIG_DESTINATION}
                              PATH_VARS BAARAAN_INCLUDE_DESTINATION)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${BAARAAN_CONFIGVERSION_FILENAME}"
                                 COMPATIBILITY SameMajorVersion)
    
# export the library target and store build directory in package registry
export(TARGETS ${BAARAAN_TARGET_NAME}
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${BAARAAN_TARGETS_FILENAME}")
export(PACKAGE ${BAARAAN_PACKAGE_NAME})

# install library target and config files
install(TARGETS ${BAARAAN_TARGET_NAME}
        EXPORT ${BAARAAN_PACKAGE_NAME})
install(DIRECTORY "include/baaraan"
        DESTINATION ${BAARAAN_INCLUDE_DESTINATION})
install(EXPORT ${BAARAAN_PACKAGE_NAME}
        FILE ${BAARAAN_TARGETS_FILENAME}
        DESTINATION ${BAARAAN_CONFIG_DESTINATION})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${BAARAAN_CONFIG_FILENAME}"
              "${CMAKE_CURRENT_BINARY_DIR}/${BAARAAN_CONFIGVERSION_FILENAME}"
        DESTINATION ${BAARAAN_CONFIG_DESTINATION})

# Exposing config variables
set(${PROJECT_NAME_CAP}_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME}: Includes Directory" FORCE)

set(${PROJECT_NAME_CAP}_DEFINITIONS
    ${Boost_DEFINITIONS} ${ARMADILLO_DEFINTIONS}
    CACHE INTERNAL "${PROJECT_NAME}: Definitions" FORCE)