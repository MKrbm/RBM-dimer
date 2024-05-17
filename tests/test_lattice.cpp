#include <gtest/gtest.h>
#include <graph.hpp>  // Include the header where create_graph is declared
#include <fstream>
#include "TestLatticesPath.hpp"

TEST(TriangularLatticeTest, BasisAndUnitCell) {
    std::string file = TEST_LATTICES_FILE_PATH;
    std::string lattice_name = "triangular lattice";
    std::string cell_name = "simple2d";

    // Check if the file exists
    std::ifstream test_file(file);
    std::string comment = "File does not exist: " + file;
    comment = comment +  "\nCurrent directory: " + std::filesystem::current_path().string();
    ASSERT_TRUE(test_file.good()) << comment;
    test_file.close();
    try {
        std::vector<size_t> size = {2, 2};  // Assuming a 2D lattice for triangular
        lattice::boundary_t boundary = lattice::boundary_t::periodic;  // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        for (auto& bond : bonds){
            ASSERT_EQ(bond.bond_type, 0);
            ASSERT_EQ(bond.bonds.size(), 2);
        }

        for (const auto& coord_vec : coordinates) {
            ASSERT_EQ(coord_vec.size(), 2);
        }

        std::cout << bonds << std::endl;
        std::cout << coordinates << std::endl;

        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], std::sqrt(3) / 2, 1e-6);

    } catch (const std::exception& e) {
        FAIL() << "Exception during test: " << e.what();
    }
}