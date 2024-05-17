#include <gtest/gtest.h>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <lattice/basis.hpp>
#include <lattice/graph_xml.hpp>

TEST(TriangularLatticeTest, BasisAndUnitCell) {
    std::ifstream is("test/triangular_lattice.xml");
    ASSERT_TRUE(is.is_open()) << "Error opening file triangular_lattice.xml";

    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(is, pt);

    std::string basis_name = "triangular lattice";
    std::string cell_name = "triangular";

    lattice::basis bs;
    read_xml(pt, basis_name, bs);

    lattice::unitcell cell;
    read_xml(pt, cell_name, cell);

    // Example assertions to test the behavior
    EXPECT_EQ(cell.dimension(), 2);
    // EXPECT_EQ(bs.vectors().size(), 2);
    // EXPECT_EQ(bs.vectors()[0], lattice::vector{1, 0});
    // EXPECT_EQ(bs.vectors()[1], lattice::vector{0.5, 0.8660254037844386});
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}