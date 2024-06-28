import numpy as np
import pytest
from netket.graph import Triangular, Lattice, Square
from spinlattice.models.dimer import DimerTrinagular




class TestDimerAbstractMethods:
    """
    Test abstract methods of DimerModel
    """

    def setup_method(self):
        self.dimer = DimerTrinagular(extent=[3, 3], pbc=True)

    def test_get_nodes_single_node(self):
        pos = self.dimer.get_node_pos(0)
        assert len(pos) == 2

    def test_get_nodes_array(self):
        pos = self.dimer.get_node_pos(np.array([0, 1]))
        assert len(pos) == 2
    
    def test_nodes_at_edge_single_edge(self):
        """
        get_edge return the nodes at the end of the edge for a given edge index
        """
        edge = self.dimer.get_nodes_from_edge(0)
        assert len(edge) == 2

    def test_nodes_at_edge_array(self):
        edges = self.dimer.get_nodes_from_edge(np.array([0, 1, 2]))
        assert edges.shape == (3, 2)
    
    def test_n_dimers(self):
        assert self.dimer.n_dimers == self.dimer.n_edges
    
    def test_adj_list_type(self):
        assert isinstance(self.dimer.adj_list, np.ndarray)
    
    def test_num_adj_list(self):
        print(self.dimer.adj_list)
        assert len(self.dimer.adj_list) == self.dimer.n_sites

    def test_adj_list(self):
        """
        check if get_connected_nodes return the all the nodes connected to node idex
        """
        node_idx = 3
        conn_nodes = self.dimer.get_connected_nodes(node_idx)
        for e in self.dimer.edges:
            if node_idx == e[0]:
                assert e[1] in conn_nodes
            elif node_idx == e[1]:
                assert e[0] in conn_nodes.tolist()
    
    def test_edges_nodeid(self):
        for node_id in range(self.dimer.n_sites):
            edges = self.dimer.get_edges_from_node(node_id)
            nodes = self.dimer.get_connected_nodes(node_id)
            assert len(edges) == len(nodes)
            nodes_edges = self.dimer.get_nodes_from_edge(edges)
            assert np.all(np.isin(nodes, nodes_edges))
    
    
class TestDimerTriangularPeriodic:

    def setup_method(self):
        self.dimer = DimerTrinagular(extent=[3, 3], pbc=True)

    def test_get_nodes_array(self):
        pos = self.dimer.get_node_pos(np.array([0, 1]))
        assert pos.shape == (2, 2)

    def test_n_sites(self):
        assert self.dimer.n_sites == 9
    
    def test_n_bonds(self):
        assert self.dimer.n_edges == 27
    
    def test_n_dimers(self):
        assert self.dimer.n_dimers == 27
    
    def test_plaquette_edges(self):
        plaquette_id = 0
        edges = self.dimer.get_plaquette_edges(plaquette_id)
        assert edges.shape == (2, 2)
        # assert len(edges) == 3
        # assert len(edges[0]) == 2


class TestDimerTriangularPeriodicNetket:

    def setup_method(self):
        self.dimer = DimerTrinagular(extent=[3, 3], pbc=True)
    
    def test_get_edge_from_ends_success(self):
        edge_id = self.dimer.get_edge_from_ends(1, 4)
        assert edge_id == 22
    
    def test_get_edge_from_ends_fail(self):
        with pytest.raises(ValueError):
            self.dimer.get_edge_from_ends(0, 4)
    
    def test_get_edges_from_plaquettes_unorder1(self):
        edges = self.dimer.get_plaquette_edges(0)
        edges_ref = np.array([[6, 1], [22, 10]])
        
        edges = np.sort(edges.reshape(-1))
        edges_ref = np.sort(edges_ref.reshape(-1))
        assert np.all(edges == edges_ref)

        
    def test_get_edges_from_plaquettes_unorder2(self):
        edges = self.dimer.get_plaquette_edges(1)
        edges_ref = np.array([0,10,17,16])
        
        edges = np.sort(edges.reshape(-1))
        edges_ref = np.sort(edges_ref.reshape(-1))
        assert np.all(edges == edges_ref)

    def test_get_edges_from_plaquettes_order(self):
        edges = self.dimer.get_plaquette_edges(1)
        edges_ref = np.array([[10, 17], [16, 0]])
        
        assert np.all(edges == edges_ref)
