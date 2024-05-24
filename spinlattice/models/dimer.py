from abc import ABC, abstractmethod
import netket as nk
from netket.operator.spin import sigmax,sigmaz, identity
import numpy as np

from ..lattice import LatticeData



class DimerModel(ABC):
    def __init__(self, lattice_data: LatticeData):
        self.lattice_data = lattice_data
        self.n_sites = lattice_data.n_sites()
        self.hi = nk.hilbert.Spin(s=1 / 2, N=self.n_sites)

    @abstractmethod
    def projection(self):
        pass


    @abstractmethod
    def hamiltonian(self, V : float, h : float):
        """
        V : float
            The strength of the dimer-dimer interaction
        h : float
            The strength of the dimer-flip
        """
        pass

    @abstractmethod
    def _dimer_potential_local(self, V : float):
        """
        Return the local operator for the dimer-dimer interaction.

        The definition depends on the lattice
        """
        pass

    def dimer_potential(self, V : float):
        return sum([self._dimer_potential_local(b) for b in range(self.lattice_data.n_bonds())])

    def dimer_flip(self, h : float):
        """
        Return the system operator for the dimer-flip term
        """
        return sum([self._dimer_flip_local(b) for b in range(self.lattice_data.n_bonds())])


    def _dimer_flip_local(self, b : int):
        """
        Flips the dimer at bond b
        """
        bond = self.lattice_data.bonds[b]
        return sigmax(self.hi, bond[0]) * sigmax(self.hi, bond[1]) 

class DimerHexagonal(DimerModel):
    def __init__(self, lattice_data: LatticeData):
        super().__init__(lattice_data)
    
    def projection(self):
        loops = self.lattice_data.loops
        bonds = self.lattice_data.bonds
        bond_colors = self.lattice_data.bond_types

        _projection_bonds = []
        _projection_bond_coords = []

        self.projections = []
        for loop in loops:
            bi = np.where(np.isin(bonds, loop).all(axis=1))[0] # bond indices
            assert len(bi) == 6
            local_projection = sum([((-1) ** c) * sigmaz(self.hi, e[0]) * sigmaz(self.hi, e[1]) for e, c in zip(bonds[bi], bond_colors[bi])])
            local_projection -= 4 * identity(self.hi)
            self.projections.append(local_projection)
        
        return sum(self.projections)
    
    def _dimer_potential_local(self, b : int):
        """
        Return the local operator for the dimer-dimer interaction.
        """
        bonds = self.lattice_data.bonds
        bond = bonds[b]

        # Get all bonds that are connected to each sites
        b1 = bond[0]
        print( np.where(np.isin(bonds, b1).any(axis=1)) )


        
