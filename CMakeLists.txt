cmake_minimum_required(VERSION 3.22)

project(RBM VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the working directory for the project
set(WORKING_DIR "${CMAKE_BINARY_DIR}" CACHE STRING "The working directory for the project")

# Include testing capabilities
include(CTest)
enable_testing()

# Print some debug information
message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")

# Add the main executable
add_executable(main main.cpp)

# Link libraries (if needed, adjust as necessary)
# target_link_libraries(main worm ${Boost_LIBRARIES})

# Include directories (if needed, adjust as necessary)
# target_include_directories(main PUBLIC "${PROJECT_BINARY_DIR}")

