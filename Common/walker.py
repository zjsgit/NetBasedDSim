

"""
Implementation of tissue-specific graph walk with RWR
"""

import sys
import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize

# convergence criterion - when vector L1 norm drops below 10^(-6)
# (this is the same as the original RWR paper)
CONV_THRESHOLD = 0.000001

class Walker:
    """ Class for multi-graph walk to convergence, using matrix computation.

    Random walk with restart (RWR) algorithm adapted from:

    Kohler S, Bauer S, Horn D, Robinson PN. Walking the interactome for
    prioritization of candidate disease genes. The American Journal of Human
    Genetics. 2008 Apr 11;82(4):949-58.

    Attributes:
    -----------
        og_matrix (np.array) : The column-normalized adjacency matrix
                               representing the original graph LCC, with no
                               nodes removed
        tsg_matrix (np.array): The column-normalized adjacency matrix
                               representing the tissue-specific graph LCC, with
                               unexpressed nodes removed as specified by
                               low_list.
        restart_prob (float) : The probability of restarting from the source
                               node for each step in run_path (i.e. r in the
                               original Kohler paper RWR formulation)
        og_prob (float)      : The probability of walking on the original graph
                               for nodes that are expressed (so, we walk on the
                               TSG with probability 1 - og_prob)
    """

    def __init__(self, network):
        self._build_matrices(network)

    def run_exp(self, source, restart_prob, og_prob):
        """ Run a multi-graph random walk experiment, and print results.

        Parameters:
        -----------
            source (list):        The source node indices (i.e. a list of Entrez
                                  gene IDs)
            restart_prob (float): As above
            og_prob (float):      As above
        """
        self.restart_prob = restart_prob
        # print 'restart_prob is {}'.format(self.restart_prob)
        self.og_prob = og_prob
        # print 'og_prob is {}'.format(self.og_prob)

        # set up the starting probability vector
        p_0 = self._set_up_p0(source)
        diff_norm = 1
        # this needs to be a deep copy, since we're reusing p_0 later
        p_t = np.copy(p_0)

        while (diff_norm > CONV_THRESHOLD):
            # first, calculate p^(t + 1) from p^(t)
            p_t_1 = self._calculate_next_p(p_t, p_0)

            # calculate L1 norm of difference between p^(t + 1) and p^(t),
            # for checking the convergence condition
            diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)

            # then, set p^(t) = p^(t + 1), and loop again if necessary
            # no deep copy necessary here, we're just renaming p
            p_t = p_t_1

        # now, generate and print a rank list from the final prob vector

        nodePercent = {}
        for node, prob in self._generate_rank_list(p_t):
            nodePercent[node]= prob
            # print '{}\t{:.10f}'.format(node, prob)
        return nodePercent


        # # fp = open("./result.txt", 'w+')
        # if node_list:
        #     for node, prob in self._generate_prob_list(p_t, node_list):
        #         print '{}\t{:.10f}'.format(node, prob)
        #     #     fp.write('{}\t{:.10f}\n'.format(node, prob))
        #     # fp.close()
        # else:
        #     for node, prob in self._generate_rank_list(p_t):
        #         print '{}\t{:.10f}'.format(node, prob)
        #     #     fp.write('{}\t{:.10f}\n'.format(node, prob))
        #     # fp.close()


    def _generate_prob_list(self, p_t, node_list):
        gene_probs = dict(zip(self.OG.nodes(), p_t.tolist()))
        for node in node_list:
            yield node, gene_probs[node]

    def _generate_rank_list(self, p_t):
        """ Return a rank list, generated from the final probability vector.

        Gene rank list is ordered from highest to lowest probability.
        """
        gene_probs = zip(self.OG.nodes(), p_t.tolist())
        # sort by probability (from largest to smallest), and generate a
        # sorted list of Entrez IDs
        for s in sorted(gene_probs, key=lambda x: x[1], reverse=True):
            yield s[0], s[1]


    def _calculate_next_p(self, p_t, p_0):
        """ Calculate the next probability vector. """

        epsilon = np.squeeze(np.asarray(np.dot(self.og_matrix, p_t)))
        no_restart = epsilon * (1 - self.restart_prob)
        restart = p_0 * self.restart_prob

        return np.add(no_restart, restart)


    def _set_up_p0(self, source):
        """ Set up and return the 0th probability vector. """
        p_0 = [0] * self.OG.number_of_nodes()

        source = set(source)

        proportion = 1 / float(len(source))

        for index, node in enumerate(self.OG.nodes()):
            if node in source:
                p_0[index] = proportion

        # for source_id in source:
        #     try:
        #         # matrix columns are in the same order as nodes in original nx
        #         # graph, so we can get the index of the source node from the OG
        #         source_index = self.OG.nodes().index(source_id)
        #         p_0[source_index] = 1 / float(len(source))
        #     except ValueError:
        #         sys.exit("Source node {} is not in original graph. Source: {}. Exiting.".format(
        #                   source_id, source))

        return np.array(p_0)


    def _build_matrices(self, network):
        """ Build column-normalized adjacency matrix for each graph.

        NOTE: these are column-normalized adjacency matrices (not nx
              graphs), used to compute each p-vector
        """
        original_graph = self._build_og(network)

        self.OG = original_graph
        og_not_normalized = nx.to_numpy_matrix(original_graph)
        self.og_matrix = self._normalize_cols(og_not_normalized)


    def _build_og(self, network):
        """ Build the original graph, without any nodes removed. """

        G = nx.Graph()
        edge_list = []

        # parse network input
        for line in network:
            if len(line) == 2:
                # assume input graph is a simple edgelist without weights
                edge_list.append((line[0], line[1], float(1)))
            else:
                # assume input graph is a simple edgelist with weights
                edge_list.append((line[0], line[1], float(line[2])))

        G.add_weighted_edges_from(edge_list)

        return G


    def _normalize_cols(self, matrix):
        """ Normalize the columns of the adjacency matrix """
        return normalize(matrix, norm='l1', axis=0)

