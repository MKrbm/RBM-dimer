import matplotlib.pyplot as plt
import numpy as np
from lattice import get_lattice_data

def plot_lattice(lattice_data, output_file="pytests/lattice_plot.png"):
    coordinates = lattice_data["coordinates"]
    bonds = lattice_data["bonds"]
    bond_types = lattice_data["bond_types"]

    # Create a figure and axis
    # fig size depends on the number of nodes
    scale = np.sqrt(len(coordinates))
    scale = int(scale)
    print(scale)
    fig, ax = plt.subplots(figsize=(scale, scale))

    # Plot nodes
    for coord in coordinates:
        ax.plot(coord[0], coord[1], 'o', color='black')

    # Define colors for bond types
    bond_type_colors = {
        "0": "blue",
        "1": "red",
        "2": "green",
        # Add more colors if there are more bond types
    }

    # Plot bonds
    for bond, bond_type in zip(bonds, bond_types):
        start = coordinates[bond[0]]
        end = coordinates[bond[1]]
        ax.plot([start[0], end[0]], [start[1], end[1]], color=bond_type_colors.get(bond_type, "black"))

    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Set title
    ax.set_title('Lattice Graph')

    # Set equal scaling
    ax.set_aspect('equal', adjustable='box')

    # Save plot to file with increased dpi
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

    # Show plot interactively if possible
    try:
        plt.show()
    except Exception as e:
        print(f"Interactive display not supported: {e}")

# Example usage
lattice_data = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [5, 10], "open")
plot_lattice(lattice_data)

