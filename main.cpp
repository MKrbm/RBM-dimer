#include <iostream>
#include <vector>
#include <graph.hpp>

int main() {
    try {
        std::string file = "../lattices.xml";
        std::string lattice_name = "chain lattice";
        std::string cell_name = "simple1d";
        std::vector<int> size = {10};  // Example size, adjust as needed
        lattice::boundary_t boundary = lattice::boundary_t::periodic;  // Example boundary condition

        // Create the graph using the function from graph.hpp
        // lattice::graph lat = create_graph(file, lattice_name, cell_name, size, boundary);

        // std::cout << "Graph created with dimension: " << lat.dimension() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
