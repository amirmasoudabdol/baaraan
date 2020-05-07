cmake_minimum_required(VERSION 3.2)

project(mvrandom_test CXX)

# Setup Tests
if (${ENABLE_TESTS})

    enable_testing()
    file(GLOB TESTSOURCES "${CMAKE_PROGRAM_PATH}/test/*.cpp")

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
                                ${BOOST_UNIT_TEST_FRAMEWORK_LIBRARY}
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

endif()