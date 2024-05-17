include(GoogleTest)

# Add test executable
add_executable(lattice-gtest test_lattice.cpp)

# Link test executable with GoogleTest and your library
target_link_libraries(lattice-gtest gtest_main spinlattice)

gtest_discover_tests(lattice-gtest)

