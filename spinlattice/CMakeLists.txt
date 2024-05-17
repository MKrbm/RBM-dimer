cmake_minimum_required(VERSION 3.22)

project(SpinLattice VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find libraries
# # Find eigen
find_package(Eigen3 REQUIRED)

# Add source files
add_library(spinlattice STATIC
    ./src/test.cpp
)

add_subdirectory(./lattice)

# Link libraries
target_link_libraries(spinlattice lattice)

# Include directories
target_include_directories(spinlattice PUBLIC include)
