# Data Import and Processing

import sys
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(project_root)

import pandas as pd
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from itertools import combinations
from tqdm import tqdm # for long evaluations
from concurrent.futures import ThreadPoolExecutor # for tsne
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN
from joblib import Parallel, delayed
from string import ascii_uppercase
import scipy.stats as stats # compute mode

class Data:

    # aqui eu teria a opção de inicializar o objeto com cada atributo principal
    # e definir @classmethod def from_file para inicializar com arquivos Data.from_file
    # além disso, para multiplos datasets, faz sentido ver o que seria variável global
    def __init__(self, file_path = "inputs/datapoints_seurat_BT_ALL.xlsx", debug=False):

        np.random.seed(1234)

        # Read the data from files
        gamatrix = pd.read_excel("inputs/act_BT_ALL_seurat.xlsx", header=None).values
        gbmatrix = pd.read_excel("inputs/sup_BT_ALL_seurat.xlsx", header=None).values

        # Extract names and the number of genes
        self.genes_names = gamatrix[0, :]
        self.genes_len = len(self.genes_names)

        if debug == True:
            print([arr.shape for arr in [gamatrix, gbmatrix, self.genes_names]])

        # Define symbolic variables
        a, sa, b, sb = sp.symbols('a sa b sb')

        # Calculate amatrix and bmatrix
        self.matrix_activation_symbolic = (gamatrix[1:, :]).T * a + gamatrix[1:, :] * np.eye(self.genes_len) * (sa - a)
        self.matrix_inhibition_symbolic = (gbmatrix[1:, :]).T * b + gbmatrix[1:, :] * np.eye(self.genes_len) * (sb - b)

        if debug == True:
            print([arr.shape for arr in [self.matrix_activation_symbolic,self.matrix_inhibition_symbolic]])

        # Create a dictionary for variable substitution
        values = {'a': 1, 'sa': 1, 'b': 1, 'sb': 1}

        # # Substitute values in amatrix and bmatrix
        def replace_symbols(array, symbols):
            for symbol, value in symbols.items():
                array = np.where(array == symbol, value, array)
            return array

        self.matrix_activation_numeric = replace_symbols(self.matrix_activation_symbolic, values)
        self.matrix_inhibition_numeric = replace_symbols(self.matrix_inhibition_symbolic, values)

        # Read the data from the file
        datapoints0 = pd.read_excel(file_path)

        # Remove the first row and convert the data to a numpy array
        self.raw_count_matrix = datapoints0.iloc[1:, :].values

        # Remove rows with all zeros, and replace 0.0 with 0 ... remove
        self.count_matrix = np.array([row for row in self.raw_count_matrix if not np.all(row == 0)], dtype=np.float64)
        self.count_matrix[self.count_matrix == 0.0] = 0

        # Calculate the number of points and the mean of self.count_matrix
        self.n_samples = len(self.count_matrix)


    def print_info(self):
        print(f"Gene names:\n {self.genes_names}")
        print(f"Number of genes: {self.genes_len}")
        print(f"Number of samples: {self.n_samples}")


    def generate_marker_genes_combinations(self, markers_list, show_values = False):

        self.markers_indices = [np.where(self.genes_names == name)[0][0] for name in markers_list \
                                                                         if name in self.genes_names]        
        
        # with repetitions
        #combinationsTotal = [(self.markers_indices[j], i) for j in range(len(self.markers_indices)) \
        #                                                             for i in range(1, self.genes_len+1)]
        self.markers_combination_list = list(combinations(self.markers_indices, 2))
        self.markers_combination_list_len = len(self.markers_combination_list)

        # Select the columns corresponding to self.markers_indices
        self.selected_datapoints = self.count_matrix[:,  np.array(self.markers_indices)]

        if show_values:
            print(f"Number of markers combinations: {self.markers_combination_list_len}")
            print(f"Markers combinations numbers: {self.markers_combination_list}")
            print(f"Markers combinations names:\n {self.genes_names[self.markers_combination_list]}")


    def generate_markers_dispersion(self, show_plots = False):

        self.markers_dispersion = [
            [
                (self.count_matrix[i, self.markers_combination_list[c][0]], \
                self.count_matrix[i, self.markers_combination_list[c][1]])
                for i in range(self.n_samples)
            ]
            for c in range(self.markers_combination_list_len)
        ]

        if show_plots:

            # Create plots for each pair of combinations in a grid
            nrows = (self.markers_combination_list_len + 1) // 2 #3  # Number of rows in the grid
            ncols = (self.markers_combination_list_len + 1) // 2 #3  # Number of columns in the grid

            fig, axes = plt.subplots(nrows, ncols, figsize=(2 * ncols, 2 * nrows))

            # Monitor progress with tqdm
            for idxcomb in tqdm(range(self.markers_combination_list_len), desc="Progress"):
                ax = axes[idxcomb // ncols, idxcomb % ncols]
                ax.scatter(
                    *zip(*self.markers_dispersion[idxcomb]),
                    s=10,
                )
                ax.set_xlabel(self.genes_names[self.markers_combination_list[idxcomb][0]])
                ax.set_ylabel(self.genes_names[self.markers_combination_list[idxcomb][1]])

            # To remove empty subplots
            for i in range(idxcomb + 1, nrows * ncols):
                fig.delaxes(axes[i // ncols, i % ncols])

            plt.tight_layout()
            plt.show()


    def generate_histograms_plots(self, nrows = 8, ncols = 5):

        fig, axes = plt.subplots(nrows, ncols, figsize=(2 * ncols, 2 * nrows))

        # Monitor progress with tqdm
        for i in tqdm(range(self.genes_len), desc="Progress"):
            ax = axes[i // ncols, i % ncols]
            ax.hist(self.count_matrix[:, i], bins=np.arange(0, 13, 1), density=True) #,edgecolor='black'
            ax.set_title(self.genes_names[i])

        # To remove empty subplots
        for k in range(i + 1, nrows * ncols):
            fig.delaxes(axes[k // ncols, k % ncols])

        plt.tight_layout()
        plt.show()


    def generate_tsne_plots(self, perplexities = [20, 40, 60, 80, 100], nrows = 2, ncols = 3):

        """ tsne on the 4D dimension """
        
        def tsne_reduce(perplexity):
            tsne = TSNE(n_components=2, perplexity=perplexity, random_state=123)
            data_reduced = tsne.fit_transform(self.selected_datapoints)
            return data_reduced

        # Run t-SNE in parallel for each perplexity value and store the results
        with ThreadPoolExecutor() as executor:
            #results = list(executor.map(tsne_reduce, perplexities))
            results = list(tqdm(executor.map(tsne_reduce, perplexities), total=len(perplexities), desc="Running t-SNE"))

        # Set up subplots for grid display
        fig, axes = plt.subplots(nrows, ncols, figsize=(8, 5))
        fig.tight_layout(pad=3.0)

        # Plot the t-SNE results
        for i, data_reduced in enumerate(results):
            ax = axes[i // ncols, i % ncols]
            ax.scatter(data_reduced[:, 0], data_reduced[:, 1], s=10)
            ax.set_xlabel("tSNE-1")
            ax.set_ylabel("tSNE-2")
            ax.set_title(f"Perplexity: {perplexities[i]}")
            ax.grid(True, linestyle="--", alpha=0.3)

        # To remove empty subplots
        for i in range(len(results), nrows * ncols):
            fig.delaxes(axes[i // ncols, i % ncols])

        plt.show()


    def generate_clustering_reduced_dimension(self, perplexity=60, export_images = False):

        """ clustering on tsne 2D dimension with chosen perplexity """

        def clustering_function(method, number, data_points):
            tsne = TSNE(n_components=2, perplexity=perplexity, random_state=123)
            data_reduced = tsne.fit_transform(data_points)

            if method == "KMeans":
                return (KMeans(n_clusters=number, n_init='auto').fit_predict(data_reduced), data_reduced)
            elif method == "DBSCAN":
                dbscan_clusters = DBSCAN(eps=2.5, min_samples=5).fit_predict(data_reduced)
                # Shift the DBSCAN labels up by one
                dbscan_clusters_shifted = dbscan_clusters #+ 1
                return (dbscan_clusters_shifted, data_reduced)
            else:
                print("Unknown or empty method provided. Please provide a valid method.")
                return None

        clustering_methods = ["KMeans", "DBSCAN"]
        number = 5

        results = Parallel(n_jobs=-1)(delayed(clustering_function)(method, number, self.selected_datapoints) for method in clustering_methods)

        def plot_clusters(clusters, data_reduced, ax):
            unique_clusters = np.unique(clusters)
            for clust in unique_clusters:
                if clust >= 0:
                    label = f'Cluster {ascii_uppercase[clust]}' if clust < len(ascii_uppercase) else f'Cluster {clust}'
                    ax.scatter(data_reduced[clusters == clust, 0], data_reduced[clusters == clust, 1], label=label)
            ax.legend()

        fig, axes = plt.subplots(nrows=1, ncols=len(clustering_methods), figsize=(5 * len(clustering_methods), 4))

        for i, method in enumerate(clustering_methods):
            clusters, data_reduced = results[i]
            plot_clusters(clusters, data_reduced, axes[i])
            axes[i].set_title(f"Method: {method}")

        if export_images:
            # Create the directory if it does not exist
            if not os.path.exists('./pdf_images'):
                os.makedirs('./pdf_images')

            plt.tight_layout()
            # plt.savefig("BeforeSim_tsneClustering.pdf")
            plt.savefig('./pdf_images/BeforeSim_tsneClustering.pdf')
            
        plt.show()


    def generate_clusters_and_statistcs(self, method = 1, number = 5): 

        # Function to calculate dataclust from components and return indices of datapoints in each cluster
        def dataclust_from_components(components, datapoints):
            dataclust = []
            cluster_indices = []  # New list to hold indices of datapoints in each cluster
            unique_components = np.unique(components)
            for comp in unique_components:
                cluster_datapoints = datapoints[components == comp]
                dataclust.append(cluster_datapoints)
                cluster_indices.append(np.where(components == comp))  # Store indices of datapoints in this cluster
            return dataclust, cluster_indices

        # Function to perform common calculations
        def common_calculations(dataclust, cluster_indices, original_datapoints):
            # databasinsmode = [stats.mode(clust)[0] for clust in dataclust]
            databasinsmode = [stats.mode(original_datapoints[indices], keepdims=True)[0] for indices in cluster_indices]
            databasinscenters = [np.mean(original_datapoints[indices], axis=0) for indices in cluster_indices]
            # databasinssigmas = [np.std(original_datapoints[indices], axis=0) for indices in cluster_indices] + 1/100
            databasinssigmas = [np.std(original_datapoints[indices], axis=0) + 1/100 for indices in cluster_indices]  # add 0.01 to each sigma
            dataweights = [len(clust) for clust in dataclust]
            dataweightsperc = [weight / sum(dataweights) for weight in dataweights]
            return databasinsmode, databasinscenters, databasinssigmas, dataweights, dataweightsperc

        if method == 1:
            # Clustering Method 1
            kmeans = KMeans(n_clusters=number, n_init='auto')
            kmeans_components = kmeans.fit_predict(self.selected_datapoints)

            # Get clusters and indices of datapoints in each cluster
            kmeans_dataclust, kmeans_cluster_indices = dataclust_from_components(kmeans_components, self.count_matrix) 

            self.cluster_kmeans_results = {"kmeans_components": kmeans_components, \
                                        "kmeans_dataclust": kmeans_dataclust,
                                        "kmeans_cluster_indices": kmeans_cluster_indices
            }

            # Calculate statistics on all original variables using only datapoints in each cluster
            result = common_calculations(kmeans_dataclust, kmeans_cluster_indices, self.count_matrix) 

            keys = ["clusters_mode", "clusters_centroids", "clusters_stds", \
                                             "clusters_weights", "clusters_proportions"]
            self.cluster_kmeans_statistics = dict(zip(keys, result))

        else:
            self.cluster_DBSCAN_statics = {}
            # Similarly for DBSCAN
            dbscan = DBSCAN(eps=1.2, min_samples=5)
            dbscan_components = dbscan.fit_predict(self.selected_datapoints)
            dbscan_dataclust, dbscan_cluster_indices = dataclust_from_components(dbscan_components, self.count_matrix) 
            results2 = common_calculations(dbscan_dataclust, dbscan_cluster_indices, self.count_matrix)


        # Combine the results in a list
        #results = [results1, results2]









