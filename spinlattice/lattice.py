import subprocess
from typing import List
from pathlib import Path
import os


def get_lattice_data(lattice_name: str, cell_name: str, size: List[int], boundary: str) -> dict:
    """
    Executes the C++ executable and retrieves the lattice data.

    Parameters:
    - lattice_name: str, name of the lattice.
    - cell_name: str, name of the cell.
    - size: list of int, size in the first and second dimensions.
    - boundary: str, boundary condition (periodic or open).

    Returns:
    - dict: containing 'bonds', 'bond_types', and 'coordinates'.
    """
    executable_path = Path(__file__).parent / "build" / "main"

    if Path(executable_path).is_file() is False:
        raise FileNotFoundError("Executable not found. "
                                "You must compile the C++ code first under the build directory."
                                "Current working directory: " + os.getcwd()
                                )

    if len(size) != 2:
        raise ValueError("Size must be a list of two integers.")
    if not all(isinstance(x, int) for x in size):
        raise ValueError("Size must be a list of two integers.")
    if boundary not in ["periodic", "open"]:
        raise ValueError("Boundary must be 'periodic' or 'open'.")

    args = [
        "--lattice_name", lattice_name,
        "--cell_name", cell_name,
        "--L1", str(size[0]),
        "--L2", str(size[1]),
        "--boundary", boundary
    ]

    try:
        # Run the executable with the provided arguments
        result = subprocess.run([executable_path] + args, capture_output=True, text=True, check=True)
        
        # Parse the output
        output = result.stdout
        bonds = []
        coordinates = []
        bond_types = []

        # Process the output to extract bonds, bond types, and coordinates
        lines = output.splitlines()
        parsing_bonds = False
        parsing_coordinates = False

        for line in lines:
            if "Bond Types and Bonds:" in line:
                parsing_bonds = True
                parsing_coordinates = False
                continue
            elif "Coordinates:" in line:
                parsing_bonds = False
                parsing_coordinates = True
                continue

            if parsing_bonds:
                if "Bond Type:" in line:
                    bond_type = line.split("Bond Type: ")[1].split(", Bonds:")[0]
                    bond_types.append(bond_type)
                    bond_list = line.split("{")[1].split("}")[0].split(", ")
                    bonds.append([int(b) for b in bond_list])
            elif parsing_coordinates:
                if "{" in line and "}" in line:
                    coord_list = line.split("{")[1].split("}")[0].split(", ")
                    coordinates.append([float(c) for c in coord_list])

        return {
            "bonds": bonds,
            "bond_types": bond_types,
            "coordinates": coordinates
        }
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing the C++ program: {e}")
        return None

# Example usage:
# lattice_data = get_lattice_data("lattice_name", "cell_name", [10, 10], "periodic")
# print(lattice_data)