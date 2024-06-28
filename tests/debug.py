import numpy as np
from scipy import sparse as sp
from netket.graph import Triangular, Lattice, Square
import sys
import os
sys.path.append(os.getcwd())
from spinlattice.models.dimer import DimerModel

g = Square(length=3, pbc=True)
dimer = DimerModel(g)

dimer.get_connected_nodes(0)


# class TestDimerAbstractMethods:

#     def setup_method(self):
#         g = Square(length=3, pbc=True)
#         dimer = DimerModel(g)
#         self.dimer = dimer

#     def test_get_nodes_single_node(self):
#         nodes = self.dimer.get_node_pos(0)
#         assert len(nodes) == 1

#     def test_get_nodes_array(self):
#         nodes = self.dimer.get_node_pos(np.array([0, 1]))
#         assert len(nodes) == 2
    
#     def test_nodes_at_edge_single_edge(self):
#         """
#         get_edge return the nodes at the end of the edge for a given edge index
#         """
#         edge = self.dimer.get_nodes_at_edge(0)
#         assert len(edge) == 1
#         assert len(edge[0]) == 2

#     def test_nodes_at_edge_array(self):
#         edge = self.dimer.get_nodes_at_edge(np.array([0, 1]))
#         assert len(edge) == 2
#         assert len(edge[0]) == 2
    
#     def test_n_dimers(self):
#         assert self.dimer.n_dimers == self.dimer.n_bonds
    
#     def test_adj_list_type(self):
#         assert isinstance(self.dimer.adj_list, np.ndarray)
    
#     def test_num_adj_list(self):
#         print(self.dimer.adj_list)
#         assert len(self.dimer.adj_list) == self.dimer.n_sites

#     def test_adj_list(self):
#         """
#         check if get_connected_nodes return the all the nodes connected to node idex
#         """
#         node_idx = 3
#         conn_nodes = self.dimer.get_connected_nodes(node_idx)
#         for e in self.dimer.edges:
#             if node_idx == e[0]:
#                 assert e[1] in conn_nodes
#             elif node_idx == e[1]:
#                 assert e[0] in conn_nodes.tolist()