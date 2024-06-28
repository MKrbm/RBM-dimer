from netket.graph import Triangular, Lattice
import matplotlib.pyplot as plt
import numpy as np

def plot_graph(graph : Lattice, max_edge_length : float = 1.0):
    """
    Plots the graph with nodes, edges, and labels, skipping edges longer than max_edge_length.

    Parameters:
    - graph: Graph object with positions and edges methods.
    - max_edge_length: float, maximum length of edges to be plotted.
    """
    positions = graph.positions
    edges = graph.edges(0)

    # Creating the plot
    plt.figure(figsize=(12, 12))
    plt.scatter(positions[:, 0], positions[:, 1], color='blue')  # Plot points

    # Plotting and labeling nodes
    for j, (x, y) in enumerate(positions):
        plt.text(x, y + 0.05, str(j), color='green', fontsize=10, ha='center', va='center')

    # Plotting edges with length condition and labels
    for i, edge in enumerate(edges):
        start, end = edge
        # Calculate the Euclidean distance between points
        distance = np.linalg.norm(positions[start] - positions[end])
        if distance <= max_edge_length:
            plt.plot([positions[start][0], positions[end][0]], [positions[start][1], positions[end][1]], 'k-')  # 'k-' for black line
            # Labeling the midpoint of the edge
            mid_point = (positions[start] + positions[end]) / 2
            plt.text(mid_point[0], mid_point[1] + 0.05, str(i), color='red', fontsize=10, ha='center', va='center')


    plt.title('Graph Points with Edges and Labels')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.grid(False)
    plt.axis('equal')
    plt.savefig('./tests/lattice.png')

g = Triangular(extent=[3, 3])
plot_graph(g)

