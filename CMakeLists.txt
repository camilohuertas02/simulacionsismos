cmake_minimum_required(VERSION 3.10)

project(SeismicSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Add executable
add_executable(seismic_sim src/main.cpp)

# Add include directories
# This line would be useful if you have header files in src/ that main.cpp includes,
# or if you add other .cpp files in src/ that need to find headers in src/.
# For a single main.cpp that might include headers from a subfolder like src/include,
# you might use: target_include_directories(seismic_sim PUBLIC ${PROJECT_SOURCE_DIR}/src/include)
# For now, assuming main.cpp is self-contained or includes are relative to src/ or standard libraries.
target_include_directories(seismic_sim PUBLIC ${PROJECT_SOURCE_DIR}/src)

# Add link libraries (if any)
# target_link_libraries(seismic_sim PRIVATE some_library)

# Installation (optional)
# install(TARGETS seismic_sim DESTINATION bin)

# Testing (optional)
# enable_testing()
# add_test(NAME MyTest COMMAND seismic_sim --test_args)
# set_tests_properties(MyTest PROPERTIES PASS_REGULAR_EXPRESSION "All tests passed")

message(STATUS "CMakeLists.txt configured for SeismicSimulation")
