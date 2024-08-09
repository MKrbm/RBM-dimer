from abc import ABC, abstractmethod, abstractproperty
import netket as nk
from netket.graph._lattice_edge_logic import create_padded_sites
from netket.operator.spin import sigmax, sigmaz, identity
from netket.graph import Lattice, Triangular
import numpy as np
from typing import List

from ..lattice import LatticeData

LocalOperator = nk.operator._local_operator.LocalOperator


class DimerModel(ABC):
    _edges : np.ndarray
    _pos : np.ndarray
    _adj_list : np.ndarray
    _num_plaquettes : int
    lattice : Lattice
    def __init__(self, g: Lattice):
        """
        Currently, we don't consider edge colors, thus, specify edge color as 0
        """
        self.lattice = g
        self._pos = self.lattice.positions
        self._edges = np.array(self.lattice.edges())
        self._adj_list = np.array(self.lattice.adjacency_list())
        self._edges_from_nodeid = [self._create_edges_from_nodeid(node_id) for node_id in range(self.n_sites)]
        self._plaquette_edges = self._get_plaquette_edges()

    def get_node_pos(self, node_id : int | np.ndarray) -> np.ndarray:
        return self._pos[node_id] if isinstance(node_id, np.ndarray) else self._pos[node_id]
        
    def get_nodes_from_edge(self, edge_id : int | np.ndarray) -> np.ndarray:
        """
        return the nodes connected to the ends of edges
        """
        return self._edges[edge_id] if isinstance(edge_id, np.ndarray) else self._edges[edge_id]
    
    def get_connected_nodes(self, node_id : int) -> np.ndarray:
        return self._adj_list[node_id]
    
    def get_edges_from_node(self, node_id : int) -> np.ndarray:
        """
        Return the edges connected to the node
        """
        return self._edges_from_nodeid[node_id]
    
    def get_edge_from_ends(self, source : int, target : int) -> int:
        """
        Return the edge index from the source to the target
        """
        try:
            return self.lattice._igraph.es.select(_source=source, _target=target)[0].index
        except:
            raise ValueError(f"Edge from {source} to {target} not found")

    def get_plaquette_edges(self, plaquette_id : int | np.ndarray) -> np.ndarray:
        return self._plaquette_edges[plaquette_id] if isinstance(plaquette_id, np.ndarray) else self._plaquette_edges[plaquette_id]

    def _create_edges_from_nodeid(self, node_id : int) -> np.ndarray:
        es = []
        for e in self.lattice._igraph.es.select(_source=node_id):
            es.append(e.index)
        return np.array(es)
    
    @abstractmethod
    def is_flippable(self, edge_id : int) -> bool:
        """
        Check if the edge is flippable
        """
        pass

    @abstractmethod
    def is_valid_dimer(self, node_id : int) -> bool:
        """
        Check if the edge is a valid dimer
        """
        pass
    
    def is_valid_configuration(self, x : np.ndarray) -> bool:
        """
        Check if the edge is a valid configuration

        x : np.ndarray
            The spin configuration
        """
        if x.shape != (self.n_dimers,):
            raise ValueError(f"Invalid configuration shape {x.shape}")
        return self._is_valid_configuration(x)
    
    @abstractmethod
    def _is_valid_configuration(self, x : np.ndarray) -> bool:
        """
        Check if the edge is a valid configuration

        x : np.ndarray
            The spin configuration
        """
        pass

    @abstractmethod
    def _get_plaquette_edges(self) -> np.ndarray:
        pass
    
    @property
    def n_sites(self) -> int:
        return self.lattice.n_nodes
    
    @property
    def n_edges(self) -> int:
        return self._edges.shape[0]
    
    @property
    def n_dimers(self) -> int:
        return self.n_edges
    
    @property
    def edges(self) -> np.ndarray:
        return self._edges
    
    @property
    def adj_list(self) -> np.ndarray:
        return self._adj_list
    
    @property
    def positions(self) -> np.ndarray:
        return self._pos
    
    @abstractproperty
    def num_plaquettes(self) -> int:
        """
        Return the number of plaquettes
        """
        return -1



class DimerTrinagular(DimerModel):
    def __init__(self, extent : List[int], pbc : bool = True):
        if not pbc:
            raise NotImplementedError("PBC is not implemented for trinagular lattice")
        g = Triangular(extent, pbc=pbc)
        super().__init__(g)
        self._num_plaquettes = self.n_edges

    @property
    def num_plaquettes(self) -> int:
        return self._num_plaquettes
    
    def is_valid_dimer(self, x : np.ndarray, node_id : int) -> bool:
        """
        Check if the edge is a valid dimer
        """
        edges = self.get_edges_from_node(node_id)
        x_ = x[edges]
        return (x_ == 1).sum() == 1
    
    def _is_valid_configuration(self, x : np.ndarray) -> bool:
        for i in range(self.n_sites):
            if not self.is_valid_dimer(x, i):
                return False
        return True
    
    def is_flippable(self,x:np.ndarray, edge_id : int) -> bool:
        plaquette_edges = self.get_plaquette_edges(edge_id)
        x_ = x[plaquette_edges]
        return (x_ == 1).sum() == 2
    
    def get_flippable_edges(self, x : np.ndarray) -> np.ndarray:
        """
        Return the edges that are flippable
        """
        edges = []
        for i in range(self.n_dimers):
            if self.is_flippable(x, i):
                edges.append(i)
        return np.array(edges)
        

    def _get_plaquette_edges(self) -> np.ndarray:

        pos, labels = create_padded_sites(self.lattice.basis_vectors, self.lattice.extent, self.lattice.site_offsets, self.lattice.pbc, 1)

        vecs = np.zeros((3, 2), dtype = np.float64)
        vecs[:2] = self.lattice.basis_vectors
        vecs[2] = vecs[1] - vecs[0]
        assert np.allclose(np.linalg.norm(vecs, axis=1), 1) 
        edges = np.zeros((self.n_edges, 2, 2), dtype = int)
        for i in range(self.n_sites):
            ci = self._pos[i] 
            for a in range(3):
                cj = ci + vecs[a]
                j = labels[_close_idx_coords(pos, cj)]
                if a != 2:
                    tc1 = ci + vecs[a+1]
                    tc2 = cj - vecs[a+1]
                else:
                    tc1 = ci - vecs[0]
                    tc2 = cj + vecs[0]
                t1 = labels[_close_idx_coords(pos, tc1)]
                t2 = labels[_close_idx_coords(pos, tc2)]
                ei = self.get_edge_from_ends(source=i, target=j)
                try: 
                    edges[ei][0, 0] = self.get_edge_from_ends(source=i, target=t1)
                    edges[ei][0, 1] = self.get_edge_from_ends(source=j, target=t2)
                    edges[ei][1, 0] = self.get_edge_from_ends(source=i, target=t2)
                    edges[ei][1, 1] = self.get_edge_from_ends(source=j, target=t1)
                except:
                    raise RuntimeError(f"Failed to create plaquette edges for edge {i} and {j}")
                    
        return edges
    



def _close_idx_coords(coors : np.ndarray, coor : np.ndarray, dist_atol : float = 1e-5) -> int:
    """
    return the index of the coor in coors that is close enough to coor
    """
    dist = np.linalg.norm(coors - coor, axis=1)
    idx = np.where(dist < dist_atol)[0]
    if len(idx) == 0:
        raise ValueError(f"No coordinate is close enough to {coor}")
    elif len(idx) > 1:
        raise ValueError(f"Multiple coordinates are close enough to {coor}")
    return idx[0]
    
    