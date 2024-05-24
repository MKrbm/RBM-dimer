import numpy as np
import pytest
import subprocess
from unittest.mock import patch
from spinlattice.lattice import get_lattice_data, LatticeData
from spinlattice.models.dimer import DimerHexagonal
from scipy import sparse as sp


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
        _, elem = self.dimer.constraint().get_conn(ones)
        assert np.array_equal(elem, np.array([0.0]))

        op = self.dimer._dimer_flip_local(0)
        conn, elem = op.get_conn(ones)
        assert len(conn) == 1, "The dimer flip operator must be non-diagonal"
        assert not np.array_equal(conn[0], ones)
        assert np.array_equal(elem, np.array([1.0]))
        x_prime = conn[0]

        # Flipping 0th bond reult in non-physical state
        conn, elem = self.dimer.constraint().get_conn(x_prime)
        assert not np.array_equal(conn[0], ones)
        assert not np.array_equal(elem, np.array([0.0]))

        # Flipping 3rd bond should result in physical state
        op = self.dimer._dimer_flip_local(3)
        conn, elem = op.get_conn(ones)
        x_prime = conn[0]
        assert not np.array_equal(x_prime, ones)
        conn, elem = self.dimer.constraint().get_conn(x_prime)
        assert not np.array_equal(conn[0], ones)
        assert np.array_equal(
            elem, np.array([0.0])
        ), "The dimer flip operator should result in physical state: operator should return 0"

    def test_effective_projection(self):
        op = self.dimer.effective_projection(100)
        conn, elem = op.get_conn(np.ones(self.dimer.n_sites))
        assert np.array_equal(elem, np.array([0.0])), "The valid dimer configuration should return 0"

        op_effective = self.dimer.effective_projection(1.0)
        eff_array = op_effective.to_sparse()
        op_constraint = self.dimer.constraint()
        constraint_array = op_constraint.to_sparse()

        assert eff_array.shape == constraint_array.shape

        result = (
            - eff_array
            + constraint_array / 2
            - sp.identity(eff_array.shape[0], format="csr")
            * (2 * len(self.dimer.loops) - 4 / 6 * self.dimer.n_bonds)
        )
        assert (result != 0).nnz == 0, "The result should contain only 0's"
        assert result.getnnz() == 0, "The result should contain only 0's"

    def test_effective_hamiltonian(self):

        lattice = get_lattice_data(
            "dimer-hexagonal-lattice",
            "dimer-hexagonal",
            [
                4,
                6,
            ],
            "periodic",
        )
        dimer = DimerHexagonal(lattice)
        ones = np.ones(dimer.n_sites)
        op = dimer.effective_hamiltonian(1.0, 1.0, 100)
        conn, elem = op.get_conn(ones)
        x_primes = conn[1:]
        sections = np.ones(len(x_primes))
        conn_prime, elem_prime = op.get_conn_flattened(x_primes, sections)
        diffs = np.diff(sections)
        assert np.all(diffs == 145)
        assert len(conn_prime) == len(elem_prime)

        non_phys_conf = conn_prime[np.where(elem_prime > 100)[0]]
        # Check dimer.constraint return non_zero elements and phys_conf 
        conn, elem = dimer.constraint().get_conn_flattened(non_phys_conf, np.ones(len(non_phys_conf)))
        # print(elem)
        assert np.all(elem > 0)

        # Check with effective_projection
        op = dimer.effective_projection(100)
        conn, elem = op.get_conn_flattened(conn, np.ones(len(conn)))
        phys_conf = conn[np.where(elem <= 100)[0]]
        _, elem_const = dimer.constraint().get_conn_flattened(phys_conf, np.ones(len(phys_conf)))
        assert np.all(elem_const == 0), "configurations having effective potential <= 100 should be valid dimer configuration"

        # print(phys_conf)
        # conn, elem = dimer.constraint().get_conn_flattened(phys_conf, np.ones(len(phys_conf)))
        # assert np.all(elem == 0)


