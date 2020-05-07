cmake_minimum_required(VERSION 3.2)

project(baaraan_tests CXX)

include(CTest)

file(GLOB TESTSOURCES "*.cpp")

message(STATUS "Configuring tests")

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

#Run through each source
foreach(testSrc ${TESTSOURCES})
    #Extract the filename without an extension (NAME_WE)
    get_filename_component(testName ${testSrc} NAME_WE)

    #Add compile target
    add_executable(${testName} ${testSrc} ${SOURCES})

    #link to Boost libraries AND your targets and dependencies
    target_link_libraries(${testName}
                            ${ARMADILLO_LIBRARIES}
                            ${BOOST_LIBRARIES}
                            Boost::unit_test_framework
                            )

    #I like to move testing binaries into a build/tests directory
    set_target_properties(${testName}
                            PROPERTIES
                            RUNTIME_OUTPUT_DIRECTORY
                            ${CMAKE_CURRENT_SOURCE_DIR}/build/tests)

    #Finally add it to tests execution -
    #Notice the WORKING_DIRECTORY and COMMAND
    add_test(NAME ${testName}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/tests
                COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/tests/${testName} )

endforeach(testSrc)

