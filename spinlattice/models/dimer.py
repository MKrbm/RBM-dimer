from abc import ABC, abstractmethod
import netket as nk
from netket.operator.spin import sigmax, sigmaz, identity
import numpy as np
from typing import List

from ..lattice import LatticeData

LocalOperator = nk.operator._local_operator.LocalOperator


class DimerModel(ABC):
    def __init__(self, lattice_data: LatticeData):
        self.lattice_data = lattice_data
        self.n_sites = lattice_data.n_sites()
        self.hi = nk.hilbert.Spin(s=1 / 2, N=self.n_sites)

    @abstractmethod
    def projection(self) -> LocalOperator:
        """
        Return the localoperator that maps the spin configurations into dimer states
        The definition of projection here is bit twisted i.e. when projection returns 0, 
        the configuration is valid dimer state.
        """
        pass

    def hamiltonian(self, V: float, h: float) -> LocalOperator:
        """
        V : float
            The strength of the dimer-dimer interaction
        h : float
            The strength of the dimer-flip
        """
        pass

    @abstractmethod
    def _dimer_potential_local(self, V: float) -> LocalOperator:
        """
        Return the local operator for the dimer-dimer interaction.

        The definition depends on the lattice
        """
        pass

    def dimer_potential(self, V: float) -> LocalOperator:
        return sum([self._dimer_potential_local(b) for b in range(self.lattice_data.n_bonds())])

    def dimer_flip(self, h: float) -> LocalOperator:
        """
        Return the system operator for the dimer-flip term
        """
        return sum([self._dimer_flip_local(b) for b in range(self.lattice_data.n_bonds())])

    def _dimer_flip_local(self, b: int) -> LocalOperator:
        """
        Flips the dimer at bond b
        """
        bond = self.lattice_data.bonds[b]
        return sigmax(self.hi, bond[0]) * sigmax(self.hi, bond[1])

    @property
    def bonds(self) -> np.ndarray:
        return self.lattice_data.bonds

    @property
    def bond_colors(self) -> np.ndarray:
        return -2 * self.lattice_data.bond_types + 1

    @property
    def loops(self) -> np.ndarray:
        return self.lattice_data.loops


class DimerHexagonal(DimerModel):
    def __init__(self, lattice_data: LatticeData):
        super().__init__(lattice_data)
        self._set_projections()

    def _set_projections(self) -> None:
        loops = self.lattice_data.loops
        bonds = self.lattice_data.bonds
        bond_colors = self.lattice_data.bond_types

        self.projections: List[LocalOperator] = []
        for loop in loops:
            bi = np.where(np.isin(bonds, loop).all(axis=1))[0]  # bond indices
            assert len(bi) == 6
            _local_operator = sum([((-1) ** c) * sigmaz(self.hi, e[0]) * sigmaz(self.hi, e[1]) for e, c in zip(bonds[bi], bond_colors[bi])])
            local_projection = (_local_operator - 4)
            self.projections.append(local_projection)

    def projection(self) -> LocalOperator:
        return sum(self.projections)

    def _dimer_potential_local(self, b: int) -> LocalOperator:
        """
        Return the local operator for the dimer-dimer interaction.
        """
        bonds = self.lattice_data.bonds
        bond_colors = self.lattice_data.bond_types
        bond = bonds[b]
        coords = self.lattice_data.coordinates

        b1 = bond[0]
        b2 = bond[1]

        conn_b1 = np.where(np.isin(bonds, b1).any(axis=1))[0]
        conn_b1 = conn_b1[conn_b1 != b]

        conn_b2 = np.where(np.isin(bonds, b2).any(axis=1))[0]
        conn_b2 = conn_b2[conn_b2 != b]

        loops = self.lattice_data.loops

        cand_zigzag1 = np.concatenate([bonds[conn_b1[0]], bonds[conn_b2[0]]])
        cand_zigzag2 = np.concatenate([bonds[conn_b1[1]], bonds[conn_b2[0]]])

        vis_cnt = np.isin(loops, cand_zigzag1).sum(axis=1)
        if (vis_cnt == 3).sum() >= 2:
            cand_zigzag2 = np.concatenate([bonds[conn_b1[1]], bonds[conn_b2[1]]])
            v1 = self._dimer_exists(conn_b1[0]) * self._dimer_exists(conn_b2[0])
            v2 = self._dimer_exists(conn_b1[1]) * self._dimer_exists(conn_b2[1])
        elif (vis_cnt == 4).sum() == 1:
            cand_zigzag1 = np.concatenate([bonds[conn_b1[0]], bonds[conn_b2[1]]])
            v1 = self._dimer_exists(conn_b1[0]) * self._dimer_exists(conn_b2[1])
            v2 = self._dimer_exists(conn_b1[1]) * self._dimer_exists(conn_b2[1])
        else:
            raise RuntimeError("Algorithm failed to find a valid zigzag")

        ## Assert if cand_zigzag1 and cand_zigzag2 are valid zigzags
        assert (np.isin(loops, cand_zigzag1).sum(axis=1) == 3).sum() >= 2
        assert (np.isin(loops, cand_zigzag2).sum(axis=1) == 3).sum() >= 2
        assert np.unique(np.concatenate([cand_zigzag1, cand_zigzag2])).shape[0] == 6

        return v1 + v2

    def _dimer_exists(self, b: int) -> LocalOperator:
        """
        Check if a dimer exists at bond b
        """
        bonds = self.bonds
        colors = self.bond_colors
        return (1 - colors[b] * sigmaz(self.hi, bonds[b][0]) * sigmaz(self.hi, bonds[b][1])) / 2

