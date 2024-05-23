#include <gtest/gtest.h>
#include <graph.hpp> // Include the header where create_graph is declared
#include <fstream>
#include <filesystem>
#include "LatticesPath.hpp"

class LatticeTest : public ::testing::Test {
protected:
    std::string file;
    std::string comment;

    void SetUp() override {
        file = LATTICES_FILE_PATH;
        // Check if the file exists
        std::ifstream test_file(file);
        comment = "File does not exist: " + file;
        comment += "\nCurrent directory: " + std::filesystem::current_path().string();
        ASSERT_TRUE(test_file.good()) << comment;
        test_file.close();
    }
};

void check_bonds(const std::vector<lattice::bond_t> &bonds,
                 const std::vector<std::vector<size_t>> &expected_bonds,
                 const std::vector<int> &bond_types)
{
    std::vector<bool> found_bonds(expected_bonds.size(), false);

    ASSERT_EQ(bonds.size(), expected_bonds.size()) << "Number of bonds is not as expected";
    ASSERT_EQ(bond_types.size(), expected_bonds.size()) << "Number of bond types is not as expected";

    for (size_t i = 0; i < bonds.size(); ++i)
    {
        auto bond = bonds[i].bonds;
        auto it = std::find(expected_bonds.begin(), expected_bonds.end(), bond);
        ASSERT_NE(it, expected_bonds.end()) << "Bond not found in expected bonds: {" << bond[0] << ", " << bond[1] << "}";
        size_t index = std::distance(expected_bonds.begin(), it);
        ASSERT_EQ(bond_types[index], bonds[i].bond_type) << "Bond type not as expected: {" << bond[0] << ", " << bond[1] << "}";
        found_bonds[index] = true;
    }

    for (size_t i = 0; i < found_bonds.size(); ++i)
    {
        ASSERT_TRUE(found_bonds[i]) << "Bond not found: {" << expected_bonds[i][0] << ", " << expected_bonds[i][1] << "}";
    }
}

TEST_F(LatticeTest, TriangularLatticeTest_BasisAndUnitCell)
{
    std::string lattice_name = "triangular lattice";
    std::string cell_name = "simple2d";

    try
    {
        std::vector<size_t> size = {2, 2};                            // Assuming a 2D lattice for triangular
        lattice::boundary_t boundary = lattice::boundary_t::periodic; // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        for (auto &bond : bonds)
        {
            ASSERT_EQ(bond.bond_type, 0);
            ASSERT_EQ(bond.bonds.size(), 2);
        }

        for (const auto &coord_vec : coordinates)
        {
            ASSERT_EQ(coord_vec.size(), 2);
        }

        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], std::sqrt(3) / 2, 1e-6);

        // order of the sites in the bonds is also tested
        std::vector<std::vector<size_t>> expected_bonds = {{0, 1}, {0, 2}, {1, 0}, {1, 3}, {2, 3}, {2, 0}, {3, 2}, {3, 1}};
        std::vector<int> bond_types = {0, 0, 0, 0, 0, 0, 0, 0};

        check_bonds(bonds, expected_bonds, bond_types);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}

TEST_F(LatticeTest, DimerHexagonalLatticeTest_2x2_period)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    try
    {
        std::vector<size_t> size = {2, 2};                            // Assuming a 2D lattice for dimer hexagonal
        lattice::boundary_t boundary = lattice::boundary_t::periodic; // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        // 4 sites in the unit cell, 4 unit cells
        ASSERT_EQ(coordinates.size(), 4 * 4);

        // 6 bonds in the unit cell, 4 unit cells
        ASSERT_EQ(bonds.size(), 6 * 4);

        for (auto &bond : bonds)
        {
            ASSERT_EQ(bond.bonds.size(), 2);
        }

        for (const auto &coord_vec : coordinates)
        {
            ASSERT_EQ(coord_vec.size(), 2);
        }
        // (0,0)
        ASSERT_NEAR(coordinates[0][0], 0.0, 1e-6);
        ASSERT_NEAR(coordinates[0][1], 0.0, 1e-6);

        // (1/3, 0)
        ASSERT_NEAR(coordinates[1][0], 1.0 / 3, 1e-6);
        ASSERT_NEAR(coordinates[1][1], 0.0, 1e-6);

        // (0.5, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], 0.5 / std::sqrt(3), 1e-6);

        // (0.5 + 1 / 3, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[3][0], 0.5 + 1.0 / 3, 1e-6);
        ASSERT_NEAR(coordinates[3][1], 0.5 / std::sqrt(3), 1e-6);

        // order of the sites in the bonds is also tested
        std::vector<std::vector<size_t>> expected_bonds = {
            {0, 1}, {1, 2}, {2, 3}, {3, 4}, {3, 12}, {2, 9}, {4, 5}, {5, 6}, {6, 7}, {7, 0}, {7, 8}, {6, 13}, {8, 9}, {9, 10}, {10, 11}, {11, 12}, {11, 4}, {10, 1}, {12, 13}, {13, 14}, {14, 15}, {15, 8}, {15, 0}, {14, 5}};
        std::vector<int> bond_types = {0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};

        check_bonds(bonds, expected_bonds, bond_types);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}

TEST_F(LatticeTest, DimerHexagonalLatticeTest_1x3_period)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    try
    {
        std::vector<size_t> size = {3, 1};                            // Assuming a 2D lattice for dimer hexagonal
        lattice::boundary_t boundary = lattice::boundary_t::periodic; // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        // 4 sites in the unit cell, 4 unit cells
        ASSERT_EQ(coordinates.size(), 4 * 3);

        // 6 bonds in the unit cell, 4 unit cells
        ASSERT_EQ(bonds.size(), 6 * 3);

        for (auto &bond : bonds)
        {
            ASSERT_EQ(bond.bonds.size(), 2);
        }

        for (const auto &coord_vec : coordinates)
        {
            ASSERT_EQ(coord_vec.size(), 2);
        }
        // (0,0)
        ASSERT_NEAR(coordinates[0][0], 0.0, 1e-6);
        ASSERT_NEAR(coordinates[0][1], 0.0, 1e-6);

        // (1/3, 0)
        ASSERT_NEAR(coordinates[1][0], 1.0 / 3, 1e-6);
        ASSERT_NEAR(coordinates[1][1], 0.0, 1e-6);

        // (0.5, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[2][0], 0.5, 1e-6);
        ASSERT_NEAR(coordinates[2][1], 0.5 / std::sqrt(3), 1e-6);

        // (0.5 + 1 / 3, 0.5 / sqrt(3))
        ASSERT_NEAR(coordinates[3][0], 0.5 + 1.0 / 3, 1e-6);
        ASSERT_NEAR(coordinates[3][1], 0.5 / std::sqrt(3), 1e-6);

        // order of the sites in the bonds is also tested
        std::vector<std::vector<size_t>> expected_bonds = {
            {0, 1}, {1, 2}, {2, 3}, {2, 1}, {3, 4}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {6, 5}, {7, 8}, {7, 8}, {8, 9}, {9, 10}, {10, 11}, {10, 9}, {11, 0}, {11, 0}};

        std::vector<int> bond_types = {
            0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};

        for (size_t i = 0; i < bonds.size(); ++i)
        {
            const auto &bond = bonds[i];
            const auto &expected_bond = expected_bonds[i];
            const auto &bond_type = bond_types[i];
            ASSERT_EQ(bond.bonds[0], expected_bond[0]) << "Bond " << i << " does not have the expected site 1";
            ASSERT_EQ(bond.bonds[1], expected_bond[1]) << "Bond " << i << " does not have the expected site 2";
            ASSERT_EQ(bond.bond_type, bond_type) << "Bond " << i << " does not have the expected bond type";
        }

        // check_bonds(bonds, expected_bonds, bond_types);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}

TEST_F(LatticeTest, DimerHexagonalLatticeTest_2x2_open)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    try
    {
        std::vector<size_t> size = {2, 2};                        // Assuming a 2D lattice for dimer hexagonal
        lattice::boundary_t boundary = lattice::boundary_t::open; // Example boundary condition

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        ASSERT_EQ(coordinates.size(), 4 * 4);
        ASSERT_EQ(bonds.size(), 17); // Counted manually
        for (auto &bond : bonds)
        {
            ASSERT_EQ(bond.bonds.size(), 2);
        }
        for (const auto &coord_vec : coordinates)
        {
            ASSERT_EQ(coord_vec.size(), 2);
        }
        std::vector<std::vector<size_t>> expected_bonds = {
            {0, 1}, {1, 2}, {2, 3}, {3, 4}, {3, 12}, {2, 9}, {4, 5}, {5, 6}, {6, 7}, {7, 0}, {7, 8}, {6, 13}, {8, 9}, {9, 10}, {10, 11}, {11, 12}, {11, 4}, {10, 1}, {12, 13}, {13, 14}, {14, 15}, {15, 8}, {15, 0}, {14, 5}};

        std::vector<int> bond_types = {
            1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};

        // check_bonds(bonds, expected_bonds, bond_types);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}


TEST_F(LatticeTest, DimerHexagonalLatticeTest_2x2_periodic_loops)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    try
    {
        std::vector<size_t> size = {2, 2};                            
        lattice::boundary_t boundary = lattice::boundary_t::periodic; 

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        std::vector<std::vector<size_t>> expected_loops = {
            {0, 1, 2, 9, 8, 7},
            {2, 3, 12, 11, 10, 9},
            {4, 5, 6, 13, 12, 3},
            {6, 7, 8, 15, 14, 13},
            {8, 9, 10, 1, 0, 15},
            {10, 11, 4, 3, 2, 1},
            {12, 13, 14, 5, 4, 11},
            {14, 15, 0, 7, 6, 5}
        };

        ASSERT_EQ(loops.size(), expected_loops.size()) << "Number of loops is not as expected";

        for (size_t i = 0; i < loops.size(); ++i)
        {
            ASSERT_EQ(loops[i].loops.size(), expected_loops[i].size()) << "Loop " << i << " size is not as expected";
            for (size_t j = 0; j < loops[i].loops.size(); ++j)
            {
                ASSERT_EQ(loops[i].loops[j], expected_loops[i][j]) << "Loop " << i << " element " << j << " is not as expected";
            }
        }
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}
TEST_F(LatticeTest, DimerHexagonalLatticeTest_8x9_open_loops)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    try
    {
        std::vector<size_t> size = {8, 9};                            
        lattice::boundary_t boundary = lattice::boundary_t::open; 

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        // 4 sites in the unit cell, 4 unit cells
        // ASSERT_EQ(coordinates.size(), 4 * 4);

        for (auto &loop : loops)
        {
            ASSERT_EQ(loop.loops.size(), 6);

            std::array<double, 2> center = {0.0, 0.0};
            for (const size_t &site : loop.loops)
            {
                center[0] += coordinates[site][0];
                center[1] += coordinates[site][1];
            }
            center[0] /= loop.loops.size();
            center[1] /= loop.loops.size();

            std::cout << "Center: " << center[0] << ", " << center[1] << "\n";

            // check if the distance to all sites is the same
            for (const size_t &site : loop.loops)
            {
                double distance = std::sqrt((coordinates[site][0] - center[0]) * (coordinates[site][0] - center[0]) +
                                            (coordinates[site][1] - center[1]) * (coordinates[site][1] - center[1]));
                ASSERT_NEAR(distance, 1.0 / 3, 1e-6);
            }
        }
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}

TEST_F(LatticeTest, DimerHexagonalLatticeTest_return_error)
{
    std::string lattice_name = "dimer hexagonal lattice";
    std::string cell_name = "dimer hexagonal";

    lattice::boundary_t boundary = lattice::boundary_t::periodic; // Example boundary condition
    try
    {
        std::vector<size_t> size1 = {0, 0}; // size 0 is not allowed
        ASSERT_THROW(create_graph(file, lattice_name, cell_name, size1, boundary), std::invalid_argument);

        std::vector<size_t> size2; // Empty size vector
        ASSERT_THROW(create_graph(file, lattice_name, cell_name, size2, boundary), std::invalid_argument);

        std::vector<size_t> size3 = {1, 3, 3}; // Must be 2D
        ASSERT_THROW(create_graph(file, lattice_name, cell_name, size3, boundary), std::invalid_argument);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Exception during test: " << e.what();
    }
}