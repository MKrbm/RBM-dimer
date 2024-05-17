cmake_minimum_required(VERSION 3.22)

project(RBM VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) #To export compile_commands.json

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

add_subdirectory(./spinlattice)

add_executable(main main.cpp)

target_link_libraries(main spinlattice)

include(CTest)
enable_testing() #enable testing has to be called in the root CMakeLists.txt

add_subdirectory(./tests)

