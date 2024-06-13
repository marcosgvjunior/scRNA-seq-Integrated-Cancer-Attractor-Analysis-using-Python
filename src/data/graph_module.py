# Graph construct, plot, statistics

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class GRN:

    def __init__(self, genes_names, genes_len, matrix_activation_numeric, matrix_inhibition_numeric, debug = False):

        self.genes_names = genes_names
        self.genes_len = genes_len

        # Create directed graphs
        self.graph_activation = nx.DiGraph()
        self.graph_activation.add_nodes_from(self.genes_names)

        # Add edges to graph_activation
        for i in range(self.genes_len):
            for j in range(self.genes_len):
                if matrix_activation_numeric[i, j] != 0:
                    self.graph_activation.add_edge(self.genes_names[i], self.genes_names[j])

        self.graph_inhibition = nx.DiGraph()
        self.graph_inhibition.add_nodes_from(self.genes_names)

        # Add edges to graph_inhibition
        for i in range(self.genes_len):
            for j in range(self.genes_len):
                if matrix_inhibition_numeric[i, j] != 0:
                    self.graph_inhibition.add_edge(self.genes_names[i], self.genes_names[j])

        self.graph_combined = nx.DiGraph()
        self.graph_combined.add_nodes_from(self.genes_names)

        # Add edges to self.graph_combined
        for i in range(self.genes_len):
            for j in range(self.genes_len):
                if matrix_activation_numeric[i, j] + matrix_inhibition_numeric[i, j] != 0:
                    self.graph_combined.add_edge(self.genes_names[i], self.genes_names[j])

        if debug:
            print("Edges in graph_activation:", self.graph_activation.edges())
            print("Edges in graph_inhibition:", self.graph_inhibition.edges())
            print("Edges in graph_combined:", self.graph_combined.edges())

    def graph_print(self):

        color_list = list(mcolors.CSS4_COLORS.values())
        colors = np.random.choice(color_list,  size=self.genes_len, replace=False)

        pos = nx.circular_layout(self.graph_activation)
        options = {
            "node_color": colors,
            "node_size": 300,
            "width": 1,
            "edge_color": "black",  # Set the edge color to black
            "with_labels": True,
            "font_weight": "bold",
            "font_size": 12,
        }

        nx.draw(self.graph_activation, pos, **options)
        plt.show()

        nx.draw(self.graph_inhibition, pos, **options)
        plt.show()

        nx.draw(self.graph_combined, pos, **options)
        plt.show()


        # Set up the plot
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

        # Draw the ad graph
        nx.draw(self.graph_activation, with_labels=True, node_color=colors, ax=axes[0])
        axes[0].set_title("graph_activation")

        # Draw the r graph
        nx.draw(self.graph_inhibition, with_labels=True, node_color=colors, ax=axes[1])
        axes[1].set_title("graph_inhibition")

        # Draw the both graph
        nx.draw(self.graph_combined, with_labels=True, node_color=colors, ax=axes[2])
        axes[2].set_title("graph_combined")

        # Show the plot
        plt.show()