from dataclasses import dataclass, field
from typing import List
import subprocess
from pathlib import Path
import os
import numpy as np

@dataclass
class LatticeData:
    bonds: np.ndarray = field(default_factory=lambda: np.array([]))
    bond_types: np.ndarray = field(default_factory=lambda: np.array([]))
    coordinates: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        if self.bonds.ndim != 2 or self.bonds.shape[1] != 2:
            raise ValueError("bonds must be an N x 2 array")
        if self.coordinates.ndim != 2 or self.coordinates.shape[1] != 2:
            raise ValueError("coordinates must be an N x 2 array")
        if self.bond_types.ndim != 1:
            raise ValueError("bond_types must be an N dimensional vector")
        if self.bonds.shape[0] != self.bond_types.shape[0]:
            raise ValueError("bonds and bond_types must have the same length")

    def __eq__(self, other):
        if not isinstance(other, LatticeData):
            return NotImplemented
        return (np.array_equal(self.bonds, other.bonds) and
                np.array_equal(self.bond_types, other.bond_types) and
                np.array_equal(self.coordinates, other.coordinates))

    def n_sites(self):
        return self.coordinates.shape[0]

def get_lattice_data(lattice_name: str, cell_name: str, size: List[int], boundary: str) -> LatticeData:
    """
    Executes the C++ executable and retrieves the lattice data.

    Parameters:
    - lattice_name: str, name of the lattice.
    - cell_name: str, name of the cell.
    - size: list of int, size in the first and second dimensions.
    - boundary: str, boundary condition (periodic or open).

    Returns:
    - LatticeData: containing 'bonds', 'bond_types', and 'coordinates'.
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
    
    # if lattice_name or cell_name contains space, replace space with _
    lattice_name = lattice_name.replace(" ", "-")
    cell_name = cell_name.replace(" ", "-")

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
                    bond_type = int(line.split("Bond Type: ")[1].split(", Bonds:")[0])
                    bond_types.append(bond_type)
                    bond_list = line.split("{")[1].split("}")[0].split(", ")
                    bonds.append([int(b) for b in bond_list])
            elif parsing_coordinates:
                if "{" in line and "}" in line:
                    coord_list = line.split("{")[1].split("}")[0].split(", ")
                    coordinates.append([float(c) for c in coord_list])

        return LatticeData(
            bonds=np.array(bonds),
            bond_types=np.array(bond_types),
            coordinates=np.array(coordinates)
        )
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.strip() if e.stderr else "An unknown error occurred"
        print("An error occurred while executing the C++ program:")
        # Raise ValueError with the error message
        raise ValueError(error_message)