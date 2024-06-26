cmake_minimum_required(VERSION 3.22)

project(RBM VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # To export compile_commands.json

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the working directory for the project
set(WORKING_DIR "${CMAKE_BINARY_DIR}" CACHE STRING "The working directory for the project")

# Print some debug information
message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")

# Include testing capabilities
option(BUILD_GMOCK OFF)
option(INSTALL_GTEST OFF)
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929 # release-1.11.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Find libraries
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

# Add source files
add_library(spinlattice STATIC
    ./src/test.cpp
    ./src/graph.cpp
)

add_subdirectory(./lattice)

message(STATUS "Boost include directories: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

# Link libraries
target_link_libraries(spinlattice PUBLIC lattice ${Boost_LIBRARIES})

# PUBLIC is used to propagate the include directories and compile options to the target that links to spinlattice
target_include_directories(spinlattice PUBLIC ./include ${Boost_INCLUDE_DIRS})

add_executable(main main.cpp)

target_link_libraries(main spinlattice)
target_include_directories(main PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) # Add the binary directory to the include path

include(CTest)
enable_testing() # enable testing has to be called in the root CMakeLists.txt

add_subdirectory(./tests)

# Set the path to the test lattices file
set(LATTICES_FILE "${CMAKE_CURRENT_SOURCE_DIR}/lattices.xml") # ${CMAKE_CURRENT_SOURCE_DIR should be project_dir/tests
configure_file(LatticesPath.hpp.in LatticesPath.hpp @ONLY)
