include(GoogleTest)

set(TEST_LATTICES_FILE "${CMAKE_CURRENT_SOURCE_DIR}/lattice_examples.xml") # ${CMAKE_CURRENT_SOURCE_DIR should be project_dir/tests
configure_file(TestLatticesPath.hpp.in TestLatticesPath.hpp @ONLY)
# Add test executable
add_executable(lattice-gtest test_lattice.cpp)

# Add the binary directory to the include path
target_include_directories(lattice-gtest PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) # ${CMAKE_CURRENT_BINARY_DIR} should be project_dir/build/tests

# Link test executable with GoogleTest and library
target_link_libraries(lattice-gtest gtest_main spinlattice)

gtest_discover_tests(lattice-gtest)

