#include <gtest/gtest.h>
#include <graph.hpp>  // Include the header where create_graph is declared
#include <fstream>
#include "TestLatticesPath.hpp"

void check_bonds(const std::vector<lattice::bond_t>& bonds, 
                 const std::vector<std::vector<size_t>>& expected_bonds, 
                 const std::vector<int>& bond_types) {
    std::vector<bool> found_bonds(expected_bonds.size(), false);

    ASSERT_EQ(bonds.size(), expected_bonds.size()) << "Number of bonds is not as expected";
    ASSERT_EQ(bond_types.size(), expected_bonds.size()) << "Number of bond types is not as expected";

    for (size_t i = 0; i < bonds.size(); ++i) {
        auto bond = bonds[i].bonds;
        auto it = std::find(expected_bonds.begin(), expected_bonds.end(), bond);
        ASSERT_NE(it, expected_bonds.end()) << "Bond not found in expected bonds: {" << bond[0] << ", " << bond[1] << "}";
        size_t index = std::distance(expected_bonds.begin(), it);
        ASSERT_EQ(bond_types[index], bonds[i].bond_type) << "Bond type not as expected: {" << bond[0] << ", " << bond[1] << "}";
        found_bonds[index] = true;
    }

    for (size_t i = 0; i < found_bonds.size(); ++i) {
        ASSERT_TRUE(found_bonds[i]) << "Bond not found: {" << expected_bonds[i][0] << ", " << expected_bonds[i][1] << "}";
    }
}

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

        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], std::sqrt(3) / 2, 1e-6);

        // order of the sites in the bonds is also tested
        std::vector<std::vector<size_t>> expected_bonds = {{0, 1}, {0, 2}, {1, 0}, {1, 3}, {2, 3}, {2, 0}, {3, 2}, {3, 1}};
        std::vector<int> bond_types = {0, 0, 0, 0, 0, 0, 0, 0};

        check_bonds(bonds, expected_bonds, bond_types);

    } catch (const std::exception& e) {
        FAIL() << "Exception during test: " << e.what();
    }
}

TEST(DimerHexagonalLatticeTest, BasisAndUnitCell) {
    std::string file = TEST_LATTICES_FILE_PATH;
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    // Check if the file exists
    std::ifstream test_file(file);
    std::string comment = "File does not exist: " + file;
    comment = comment +  "\nCurrent directory: " + std::filesystem::current_path().string();
    ASSERT_TRUE(test_file.good()) << comment;
    test_file.close();
    try {
        std::vector<size_t> size = {2, 2};  // Assuming a 2D lattice for dimer hexagonal
        lattice::boundary_t boundary = lattice::boundary_t::periodic;  // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        // 4 sites in the unit cell, 4 unit cells
        ASSERT_EQ(coordinates.size(), 4 * 4);

        // 6 bonds in the unit cell, 4 unit cells
        ASSERT_EQ(bonds.size(), 6 * 4);

        for (auto& bond : bonds){
            ASSERT_EQ(bond.bonds.size(), 2);
        }

        for (const auto& coord_vec : coordinates) {
            ASSERT_EQ(coord_vec.size(), 2);
        }
        // (0,0)
        ASSERT_NEAR(coordinates[0][0], 0.0, 1e-6);
        ASSERT_NEAR(coordinates[0][1], 0.0, 1e-6);

        // (1/3, 0)
        ASSERT_NEAR(coordinates[1][0], 1.0/3, 1e-6);
        ASSERT_NEAR(coordinates[1][1], 0.0, 1e-6);

        // (0.5, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], 0.5 / std::sqrt(3), 1e-6);

        // (0.5 + 1 / 3, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[3][0], 0.5 + 1.0 / 3, 1e-6);
        ASSERT_NEAR(coordinates[3][1], 0.5 / std::sqrt(3), 1e-6);


        std::cout << bonds << std::endl;

        // order of the sites in the bonds is also tested
        std::vector<std::vector<size_t>> expected_bonds = {
            {0, 1}, {1,2}, {2,3}, {3,4}, {3,12}, {2,9},{4,5},{5,6},{6,7},{7,0},{7,8},{6,13},
            {8,9}, {9,10}, {10,11}, {11,12}, {11,4}, {10,1}, {12,13}, {13,14}, {14,15}, {15,8}, {15,0}, {14,5}
        };
        std::vector<int> bond_types = {1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0};

        check_bonds(bonds, expected_bonds, bond_types);

    } catch (const std::exception& e) {
        FAIL() << "Exception during test: " << e.what();
    }
}