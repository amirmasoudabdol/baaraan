cmake_minimum_required(VERSION 3.2)

project(mvrandom
        LANGUAGES CXX
        VERSION 0.1.0.0)
set(${PROJECT_NAME}_VERSION_NAME "Appa")

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_CAP)

set(CMAKE_CXX_STANDARD 17)
set(MVRANDOM_DIR ${CMAKE_PROGRAM_PATH})

include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


if(NOT APPLE)
    set(CMAKE_INSTALL_RPATH ${ORIGIN})
endif()


file(GLOB CPP_FILES *.cpp)

find_package(BOOST REQUIRED)
find_package(Armadillo REQUIRED)

include_directories(include)
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIRS})

add_library(mvrandom SHARED ${CPP_FILES})

target_link_libraries(mvrandom ${ARMADILLO_LIBRARIES} ${BOOST_LIBRARIES})

write_basic_package_version_file(
        mvrandomConfigVersion.cmake
        VERSION ${CMAKE_PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

# Install scripts
install(TARGETS mvrandom
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT           ${PROJECT_NAME_CAP}_Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT           ${PROJECT_NAME_CAP}_Runtime
        NAMELINK_COMPONENT  ${PROJECT_NAME_CAP}_Development
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT           ${PROJECT_NAME_CAP}_Development)

# Exposing config variables
set(${PROJECT_NAME_CAP}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include
                                     ${BOOST_INCLUDE_DIRS}
                                     ${ARMADILLO_INCLUDE_DIRS}
        CACHE INTERNAL "${PROJECT_NAME}: Includes Directory" FORCE)

set(${PROJECT_NAME_CAP}_DEFINITIONS ${BOOST_DEFINITIONS}
                                    ${ARMADILLO_DEFINTIONS}
        CACHE INTERNAL "${PROJECT_NAME}: Definitions" FORCE)