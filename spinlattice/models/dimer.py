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
        g = Triangular(extent, pbc=pbc)
        super().__init__(g)
        self._num_plaquettes = self.n_edges

    @property
    def num_plaquettes(self) -> int:
        return self._num_plaquettes
    
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
    
    
        
        
    
    # def _get_edges_from_plaquettes(self, plaquette_id : int) -> np.ndarray:
    #     """
    #     Get the edges connected to the specified plaquette.

    #     Parameters
    #     ----------
    #     plaquette_id : int
    #         The identifier of the plaquette. In this case, it is equivalent to the edge id.

    #     Returns
    #     -------
    #     np.ndarray
    #         An array of edges connected to the specified plaquette.
    #     """

    #     nodes = self.get_nodes_from_edge(plaquette_id)
    #     edges1 = self.get_edges_from_node(nodes[0])
    #     edges2 = self.get_edges_from_node(nodes[1])
    #     edges = set(edges1) | set(edges2)
    #     edges -= set([plaquette_id])
    #     _nodes = self.get_nodes_from_edge(np.array(list(edges)))
    #     vals, cnts = np.unique(_nodes, return_counts=True)
    #     idx = np.where(cnts == 2)[0].tolist()
    #     if len(idx) == 3:
    #         pos = self.get_node_pos(nodes)
    #         # NOTE: Find inappropriate vertex. It should be the one that is parallel to bonds
    #         for i in idx:
    #             pos_v = self.get_node_pos(i)
    #             vec = pos -  pos_v
    #             # Note: Check if the vec is parallel or not
    #             cos_sim = np.dot(vec[0], vec[1]) / (np.linalg.norm(vec[0]) * np.linalg.norm(vec[1]))
    #             if abs(abs(cos_sim) - 1) < 1e-5:
    #                 idx.remove(i)
    #                 break
    #     if len(idx) != 2:
    #         raise RuntimeError(f"Error in get_edges_from_plaquettes: {nodes}")

    #     node1 = int(vals[idx[0]])
    #     node2 = int(vals[idx[1]])
    #     p_edges = np.zeros((2, 2), dtype=int)
    #     p_edges[0, 0] = self.get_edge_from_ends(source=node1, target=nodes[0])
    #     p_edges[0, 1] = self.get_edge_from_ends(source=node2, target=nodes[1])
    #     p_edges[1, 0] = self.get_edge_from_ends(source=node1, target=nodes[1])
    #     p_edges[1, 1] = self.get_edge_from_ends(source=node2, target=nodes[0])
    #     return p_edges
        



        
        



# class DimerModel(ABC):
#     def __init__(self, lattice_data: LatticeData):
#         self.lattice_data = lattice_data
#         self.n_sites = lattice_data.n_sites()
#         self.hi = nk.hilbert.Spin(s=1 / 2, N=self.n_sites)

#     @abstractmethod
#     def constraint(self) -> LocalOperator:
#         """
#         Return the localoperator that maps the spin configurations into dimer states
#         The definition of constraint here is bit twisted i.e. when constraint returns 0, 
#         the configuration is valid dimer state.
#         """
#         pass

#     @abstractmethod
#     def _dimer_potential_local(self, b: int) -> LocalOperator:
#         """
#         Return the local operator for the dimer-dimer interaction.

#         The definition depends on the lattice
#         """
#         pass

#     def effective_projection(self, J: float) -> LocalOperator:
#         """
#         Return the local operator to effectively project the system onto the dimer states
#         """
#         op =  sum([self._effective_projection_local(b, J) for b in range(self.n_bonds)])  - 4 / 6 * self.n_bonds
#         return - J * op
    
#     def _effective_projection_local(self, b: int, J: float) -> LocalOperator:

#         return sigmaz(self.hi, self.bonds[b][0]) * sigmaz(self.hi, self.bonds[b][1]) * self.bond_colors[b]




#     def hamiltonian(self, V: float, h: float) -> LocalOperator:
#         """
#         V : float
#             The strength of the dimer-dimer interaction
#         h : float
#             The strength of the dimer-flip
#         """
#         return self.dimer_potential(V) + self.dimer_flip(h)
    
#     def effective_hamiltonian(self, V: float, h: float, J: float) -> LocalOperator:
#         """
#         V : float
#             The strength of the dimer-dimer interaction
#         h : float
#             The strength of the dimer-flip
#         J : float
#             The strength of the effective projection
#         """
#         return self.effective_projection(J) + self.dimer_potential(V) + self.dimer_flip(h)

#     def dimer_potential(self, V: float) -> LocalOperator:
#         return sum([self._dimer_potential_local(b) for b in range(self.n_bonds)]) * V

#     def dimer_flip(self, h: float) -> LocalOperator:
#         """
#         Return the system operator for the dimer-flip term
#         """
#         return -1 * sum([self._dimer_flip_local(b) for b in range(self.n_bonds)]) * h

#     def _dimer_flip_local(self, b: int) -> LocalOperator:
#         """
#         Flips the dimer at bond b
#         """
#         bond = self.lattice_data.bonds[b]
#         return sigmax(self.hi, bond[0]) * sigmax(self.hi, bond[1])

#     @property
#     def bonds(self) -> np.ndarray:
#         return self.lattice_data.bonds

#     @property
#     def bond_colors(self) -> np.ndarray:
#         return -2 * self.lattice_data.bond_types + 1

#     @property
#     def n_bonds(self) -> int:
#         return len(self.bonds)

#     @property
#     def loops(self) -> np.ndarray:
#         return self.lattice_data.loops


# class DimerHexagonal(DimerModel):
#     def __init__(self, lattice_data: LatticeData):
#         super().__init__(lattice_data)
#         self._set_constraints()

#     def _set_constraints(self) -> None:
#         loops = self.lattice_data.loops
#         bonds = self.lattice_data.bonds
#         bond_colors = self.lattice_data.bond_types

#         self.constraints: List[LocalOperator] = []
#         for loop in loops:
#             bi = np.where(np.isin(bonds, loop).all(axis=1))[0]  # bond indices
#             assert len(bi) == 6
#             _local_operator = sum([((-1) ** c) * sigmaz(self.hi, e[0]) * sigmaz(self.hi, e[1]) for e, c in zip(bonds[bi], bond_colors[bi])])
#             local_constraint = -(_local_operator - 4)
#             self.constraints.append(local_constraint)

#     def constraint(self) -> LocalOperator:
#         return sum(self.constraints)

#     def _dimer_potential_local(self, b: int) -> LocalOperator:
#         """
#         Return the local operator for the dimer-dimer interaction.
#         """
#         bonds = self.lattice_data.bonds
#         bond_colors = self.lattice_data.bond_types
#         bond = bonds[b]
#         coords = self.lattice_data.coordinates

#         b1 = bond[0]
#         b2 = bond[1]

#         conn_b1 = np.where(np.isin(bonds, b1).any(axis=1))[0]
#         conn_b1 = conn_b1[conn_b1 != b]

#         conn_b2 = np.where(np.isin(bonds, b2).any(axis=1))[0]
#         conn_b2 = conn_b2[conn_b2 != b]

#         loops = self.lattice_data.loops

#         cand_zigzag1 = np.concatenate([bonds[conn_b1[0]], bonds[conn_b2[0]]])
#         cand_zigzag2 = np.concatenate([bonds[conn_b1[1]], bonds[conn_b2[0]]])

#         vis_cnt = np.isin(loops, cand_zigzag1).sum(axis=1)
#         if (vis_cnt == 3).sum() >= 2:
#             cand_zigzag2 = np.concatenate([bonds[conn_b1[1]], bonds[conn_b2[1]]])
#             v1 = self._dimer_exists(conn_b1[0]) * self._dimer_exists(conn_b2[0])
#             v2 = self._dimer_exists(conn_b1[1]) * self._dimer_exists(conn_b2[1])
#         elif (vis_cnt == 4).sum() == 1:
#             cand_zigzag1 = np.concatenate([bonds[conn_b1[0]], bonds[conn_b2[1]]])
#             v1 = self._dimer_exists(conn_b1[0]) * self._dimer_exists(conn_b2[1])
#             v2 = self._dimer_exists(conn_b1[1]) * self._dimer_exists(conn_b2[0])
#         else:
#             raise RuntimeError("Algorithm failed to find a valid zigzag")

#         ## Assert if cand_zigzag1 and cand_zigzag2 are valid zigzags
#         assert (np.isin(loops, cand_zigzag1).sum(axis=1) == 3).sum() >= 2
#         assert (np.isin(loops, cand_zigzag2).sum(axis=1) == 3).sum() >= 2
#         assert np.unique(np.concatenate([cand_zigzag1, cand_zigzag2])).shape[0] == 6

#         return v1 + v2

#     def _dimer_exists(self, b: int) -> LocalOperator:
#         """
#         Check if a dimer exists at bond b
#         """
#         bonds = self.bonds
#         colors = self.bond_colors
#         return (1 - colors[b] * sigmaz(self.hi, bonds[b][0]) * sigmaz(self.hi, bonds[b][1])) / 2

