import matplotlib.pyplot as plt
import numpy as np
# from ..lattice import get_lattice_data
import sys
sys.path.append("..")
from spinlattice.lattice import get_lattice_data

def plot_lattices(lattices_data, node_colors, bond_type_colors_list):
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 8))

    for lattice_data, node_color, bond_type_colors in zip(lattices_data, node_colors, bond_type_colors_list):
        coordinates = lattice_data.coordinates
        bonds = lattice_data.bonds
        bond_types = lattice_data.bond_types

        # Plot bonds
        for i, (bond, bond_type) in enumerate(zip(bonds, bond_types)):
            start = coordinates[bond[0]]
            end = coordinates[bond[1]]
            ax.plot([start[0], end[0]], [start[1], end[1]], color=bond_type_colors.get(bond_type, "black"))
            # Add bond label
            mid_point = [(start[0] + end[0]) / 2, (start[1] + end[1]) / 2]
            ax.text(mid_point[0], mid_point[1] + 0.02, str(i), color=bond_type_colors.get(bond_type, "black"), fontsize=8, ha='center', va='bottom')
        # Plot nodes
        for j, coord in enumerate(coordinates):
            ax.plot(coord[0], coord[1], 'o', color=node_color)
            # Add node label
            ax.text(coord[0], coord[1] + 0.05, str(j), color=node_color, fontsize=8, ha='center')

    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Set title
    ax.set_title('Lattice Graph')

    # Set equal scaling
    ax.set_aspect('equal', adjustable='box')

    return fig, ax

# Example usage
unit_cells = ["dimer-triangular", "dimer-hexagonal"]
lattices_data = [
    # get_lattice_data("dimer-hexagonal-lattice", unit_cells[0], [2, 2], "open"),
    get_lattice_data("dimer-hexagonal-lattice", unit_cells[1], [2, 2], "open")
]

# Define node colors and bond type colors for each lattice
node_colors = ["black", "blue"]
bond_type_colors_list = [
    {0: "gray", 1: "red", 2: "green"},
    # {0: "blue", 1: "magenta", 2: "yellow"}
]

output_file = "tests/multiple_lattices_plot.png"
fig, ax = plot_lattices(lattices_data, node_colors, bond_type_colors_list)

# Save plot to file with increased dpi
plt.savefig(output_file, dpi=300)
print(f"Plot saved to {output_file}")

# Show plot interactively if possible
try:
    plt.show()
except Exception as e:
    print(f"Interactive display not supported: {e}")