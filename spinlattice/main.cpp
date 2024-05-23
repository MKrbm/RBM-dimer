#include <iostream>
#include <vector>
#include <graph.hpp>
#include <boost/program_options.hpp>
#include <algorithm>
#include <iomanip> // Include this header for std::setprecision
#include "LatticesPath.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    try {
        std::string file;
        std::string lattice_name;
        std::string cell_name;
        std::vector<size_t> size;
        std::string boundary_str;
        lattice::boundary_t boundary;
        std::string lattices_file = LATTICES_FILE_PATH;

        // Define and parse the program options
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("file,f", po::value<std::string>(&file)->default_value(lattices_file), "path to the lattice file")
            ("boundary,b", po::value<std::string>(&boundary_str)->default_value("periodic"), "boundary condition (periodic or open)")
            ("lattice_name,l", po::value<std::string>(&lattice_name)->required(), "name of the lattice")
            ("cell_name,c", po::value<std::string>(&cell_name)->required(), "name of the cell")
            ("L1", po::value<int>()->required(), "size in the first dimension")
            ("L2", po::value<int>()->required(), "size in the second dimension");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }

        po::notify(vm);

        // Replace horizontal bars with spaces in lattice_name and cell_name
        std::replace(lattice_name.begin(), lattice_name.end(), '-', ' ');
        std::replace(cell_name.begin(), cell_name.end(), '-', ' ');

        size.push_back(vm["L1"].as<int>());
        size.push_back(vm["L2"].as<int>());

        if (boundary_str == "periodic") {
            boundary = lattice::boundary_t::periodic;
        } else if (boundary_str == "open") {
            boundary = lattice::boundary_t::open;
        } else {
            throw std::invalid_argument("Invalid boundary condition specified");
        }

        // Create the graph using the function from graph.hpp
        auto [bonds, loops, coordinates] = create_graph(file, lattice_name, cell_name, size, boundary);

        // Print out bond types, bonds, and coordinates
        std::cout << "Bond Types and Bonds:\n";
        for (const auto& bond : bonds) {
            std::cout << "Bond Type: " << bond.bond_type << ", Bonds: {";
            for (size_t i = 0; i < bond.bonds.size(); ++i) {
                std::cout << bond.bonds[i];
                if (i < bond.bonds.size() - 1) std::cout << ", ";
            }
            std::cout << "}\n";
        }

        std::cout << "Loops:\n";
        for (const auto& loop : loops) {
            std::cout << loop << "\n";
        }

        std::cout << "\nCoordinates:\n";
        std::cout << std::fixed << std::setprecision(10); // Set precision for floating-point output
        for (const auto& coord_vec : coordinates) {
            std::cout << "{";
            for (size_t i = 0; i < coord_vec.size(); ++i) {
                std::cout << coord_vec[i];
                if (i < coord_vec.size() - 1) std::cout << ", ";
            }
            std::cout << "}\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}