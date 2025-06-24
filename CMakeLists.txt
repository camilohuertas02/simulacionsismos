cmake_minimum_required(VERSION 3.10)

project(SeismicSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Add executable
add_executable(seismic_sim src/main.cpp)

# Add include directories - FORMA CORREGIDA Y ROBUSTA
target_include_directories(seismic_sim PUBLIC "include")
target_include_directories(seismic_sim PUBLIC "src")


# Add link libraries (if any)
# target_link_libraries(seismic_sim PRIVATE some_library)

# Installation (optional)
# install(TARGETS seismic_sim DESTINATION bin)

# Testing (optional)
# enable_testing()
# add_test(NAME MyTest COMMAND seismic_sim --test_args)
# set_tests_properties(MyTest PROPERTIES PASS_REGULAR_EXPRESSION "All tests passed")

message(STATUS "CMakeLists.txt configured for SeismicSimulation")
