import numpy as np
import pytest
import subprocess
from unittest.mock import patch
from spinlattice.lattice import get_lattice_data, LatticeData
from spinlattice.models.dimer import DimerHexagonal


class TestDimerHexagonalTwoByTwoPeriodic:

    def setup_method(self):
        self.lattice = get_lattice_data(
            "dimer-hexagonal-lattice",
            "dimer-hexagonal",
            [
                2,
                2,
            ],
            "periodic",
        )
        self.dimer = DimerHexagonal(self.lattice)

    def test_dimer_exists(self):
        ones = np.ones(self.dimer.n_sites)

        for i, c in enumerate(self.dimer.bond_colors):
            op = self.dimer._dimer_exists(i)
            conn, elem = op.get_conn(ones)
            assert len(conn) == 1, "The dimer_exist operator must be diagonal"
            assert bool(np.round(elem)) == (
                c == -1
            ), f"The operator returns unexpected behavior for bond = {i}"

    def test_potential(self):
        ones = np.ones(self.dimer.n_sites)

        try:
            for i, b in enumerate(self.dimer.bonds):
                op = self.dimer._dimer_potential_local(i)
                conn, elem = op.get_conn(ones)
                assert (
                    len(conn) == 1
                ), "The dimer_potential must connected to only itself"
                assert np.array_equal(
                    conn[0], ones
                ), "The dimer_potential must connected to only itself"
        except Exception as e:
            pytest.fail(f"Unexpected error occurred: {e}")

        # bond 0
        op = self.dimer._dimer_potential_local(0)
        conn, elem = op.get_conn(ones)
        assert np.array_equal(elem, np.array([0.0]))

        # bond 1 should return 0 as well
        op = self.dimer._dimer_potential_local(1)
        conn, elem = op.get_conn(ones)
        assert np.array_equal(elem, np.array([0.0]))

        # bond 3 return 1.0
        op = self.dimer._dimer_potential_local(3)
        conn, elem = op.get_conn(ones)
        assert np.array_equal(elem, np.array([1.0]))

        # bond 9 should return 1.0
        op = self.dimer._dimer_potential_local(9)
        conn, elem = op.get_conn(ones)
        assert len(conn) == 1
        assert np.array_equal(elem, np.array([1.0]))

    def test_dimer_flip(self):
        ones = np.ones(self.dimer.n_sites)
        # Test ones is valid state
        _, elem = self.dimer.projection().get_conn(ones)
        assert np.array_equal(elem, np.array([0.0]))

        op = self.dimer._dimer_flip_local(0)
        conn, elem = op.get_conn(ones)
        assert len(conn) == 1, "The dimer flip operator must be non-diagonal"
        assert not np.array_equal(conn[0], ones)
        assert np.array_equal(elem, np.array([1.0]))
        x_prime = conn[0]

        # Flipping 0th bond reult in non-physical state
        conn, elem = self.dimer.projection().get_conn(x_prime)
        assert not np.array_equal(conn[0], ones)
        assert not np.array_equal(elem, np.array([0.0]))

        # Flipping 3rd bond should result in physical state
        op = self.dimer._dimer_flip_local(3)
        conn, elem = op.get_conn(ones)
        x_prime = conn[0]
        assert not np.array_equal(x_prime, ones)
        conn, elem = self.dimer.projection().get_conn(x_prime)
        assert not np.array_equal(conn[0], ones)
        assert np.array_equal(
            elem, np.array([0.0])
        ), "The dimer flip operator should result in physical state: operator should return 0"
