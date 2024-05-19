import pytest
import subprocess
from unittest.mock import patch
from spinlattice.lattice import get_lattice_data



class TestGetLatticeData:

    @patch("spinlattice.lattice.subprocess.run")
    def test_get_lattice_data_success(self, mock_run):
        # Mock the subprocess.run return value
        mock_run.return_value.stdout = (
            "Bond Types and Bonds:\n"
            "Bond Type: 1, Bonds: {0, 1}\n"
            "Bond Type: 2, Bonds: {2, 3}\n"
            "Coordinates:\n"
            "{0.0, 0.0}\n"
            "{1.0, 1.0}\n"
        )
        mock_run.return_value.returncode = 0

        expected_output = {
            "bonds": [[0, 1], [2, 3]],
            "bond_types": ["1", "2"],
            "coordinates": [[0.0, 0.0], [1.0, 1.0]],
        }

        result = get_lattice_data("lattice_name", "cell_name", [10, 10], "periodic")
        assert result == expected_output

    @patch("spinlattice.lattice.subprocess.run")
    def test_get_lattice_data_success(self, mock_run):
        # Mock the subprocess.run return value
        mock_run.return_value.stdout = (
            "Bond Types and Bonds:\n"
            "Bond Type: 1, Bonds: {0, 1}\n"
            "Bond Type: 0, Bonds: {1, 2}\n"
            "Bond Type: 0, Bonds: {2, 3}\n"
            "Bond Type: 0, Bonds: {2, 9}\n"
            "Bond Type: 0, Bonds: {3, 4}\n"
            "Bond Type: 0, Bonds: {3, 12}\n"
            "Bond Type: 1, Bonds: {4, 5}\n"
            "Bond Type: 0, Bonds: {5, 6}\n"
            "Bond Type: 0, Bonds: {6, 7}\n"
            "Bond Type: 0, Bonds: {6, 13}\n"
            "Bond Type: 0, Bonds: {7, 0}\n"
            "Bond Type: 0, Bonds: {7, 8}\n"
            "Bond Type: 1, Bonds: {8, 9}\n"
            "Bond Type: 0, Bonds: {9, 10}\n"
            "Bond Type: 0, Bonds: {10, 11}\n"
            "Bond Type: 0, Bonds: {10, 1}\n"
            "Bond Type: 0, Bonds: {11, 12}\n"
            "Bond Type: 0, Bonds: {11, 4}\n"
            "Bond Type: 1, Bonds: {12, 13}\n"
            "Bond Type: 0, Bonds: {13, 14}\n"
            "Bond Type: 0, Bonds: {14, 15}\n"
            "Bond Type: 0, Bonds: {14, 5}\n"
            "Bond Type: 0, Bonds: {15, 8}\n"
            "Bond Type: 0, Bonds: {15, 0}\n"
            "Coordinates:\n"
            "{0.0000000000, 0.0000000000}\n"
            "{0.3333333333, 0.0000000000}\n"
            "{0.5000000000, 0.2886751346}\n"
            "{0.8333333333, 0.2886751346}\n"
            "{1.0000000000, 0.0000000000}\n"
            "{1.3333333333, 0.0000000000}\n"
            "{1.5000000000, 0.2886751346}\n"
            "{1.8333333333, 0.2886751346}\n"
            "{0.0000000000, 0.5773502692}\n"
            "{0.3333333333, 0.5773502692}\n"
            "{0.5000000000, 0.8660254038}\n"
            "{0.8333333333, 0.8660254038}\n"
            "{1.0000000000, 0.5773502692}\n"
            "{1.3333333333, 0.5773502692}\n"
            "{1.5000000000, 0.8660254038}\n"
            "{1.8333333333, 0.8660254038}\n"
        )
        mock_run.return_value.returncode = 0

        expected_output = {
            "bonds": [
                [0, 1],
                [1, 2],
                [2, 3],
                [2, 9],
                [3, 4],
                [3, 12],
                [4, 5],
                [5, 6],
                [6, 7],
                [6, 13],
                [7, 0],
                [7, 8],
                [8, 9],
                [9, 10],
                [10, 11],
                [10, 1],
                [11, 12],
                [11, 4],
                [12, 13],
                [13, 14],
                [14, 15],
                [14, 5],
                [15, 8],
                [15, 0],
            ],
            "bond_types": [
                "1",
                "0",
                "0",
                "0",
                "0",
                "0",
                "1",
                "0",
                "0",
                "0",
                "0",
                "0",
                "1",
                "0",
                "0",
                "0",
                "0",
                "0",
                "1",
                "0",
                "0",
                "0",
                "0",
                "0",
            ],
            "coordinates": [
                [0.0000000000, 0.0000000000],
                [0.3333333333, 0.0000000000],
                [0.5000000000, 0.2886751346],
                [0.8333333333, 0.2886751346],
                [1.0000000000, 0.0000000000],
                [1.3333333333, 0.0000000000],
                [1.5000000000, 0.2886751346],
                [1.8333333333, 0.2886751346],
                [0.0000000000, 0.5773502692],
                [0.3333333333, 0.5773502692],
                [0.5000000000, 0.8660254038],
                [0.8333333333, 0.8660254038],
                [1.0000000000, 0.5773502692],
                [1.3333333333, 0.5773502692],
                [1.5000000000, 0.8660254038],
                [1.8333333333, 0.8660254038],
            ],
        }

        result = get_lattice_data("lattice_name", "cell_name", [10, 10], "periodic")
        assert result == expected_output

    @patch("spinlattice.lattice.subprocess.run")
    def test_get_lattice_data_failure(self, mock_run, capfd):
        # Mock the subprocess.run to raise an exception
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "cmd", "An error occurred"
        )

        result = get_lattice_data("lattice_name", "cell_name", [10, 10], "periodic")
        assert result is None

        # Capture the output
        captured = capfd.readouterr()
        assert "An error occurred while executing the C++ program" in captured.out


    def test_1x1_open(self):
        expected_output = {
            "bonds": [
                [0, 1], [1, 2], [2, 3]
            ],
            "bond_types": [
                "1", "0", "0"
            ],
            "coordinates": [
                [0.0000000000, 0.0000000000], [0.3333333333, 0.0000000000], [0.5000000000, 0.2886751346], [0.8333333333, 0.2886751346]
            ]
        }

        result = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [1, 1], "open")
        assert result == expected_output
    
    def test_100x100_periodic(self):
        result = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [100, 100], "periodic")
        assert result is not None
        assert len(result["bonds"]) == 100 * 100 * 6

    def test_improper_boundary_condition(self):
        with pytest.raises(ValueError, match="Boundary must be 'periodic' or 'open'."):
            get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [1, 1], "invalid_boundary")

    def test_incorrect_size_list_length(self):
        with pytest.raises(ValueError, match="Size must be a list of two integers."):
            get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [0], "open")

    def test_size_list_contains_non_integer(self):
        with pytest.raises(ValueError, match="Size must be a list of two integers."):
            get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [1, "a"], "open")

    def test_improper_lattice_name(self):
        # This test will fail because the lattice name is not valid
        result = get_lattice_data("invalid_lattice_name", "dimer-hexagonal", [1, 1], "periodic")
        assert result is None

    def test_improper_cell_name(self):
        # This test will fail because the cell name is not valid
        result = get_lattice_data("dimer-hexagonal-lattice", "invalid_cell_name", [1, 1], "open")
        assert result is None

