# coding: utf-8

# -------------------------------------------------------------------------------
# Name:         node2vec
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/20
# -------------------------------------------------------------------------------

import numpy as np
import networkx as nx
import random
import time
from gensim.models import Word2Vec
import warnings
warnings.filterwarnings(action='ignore')


from Util import FileUtil

class Graph():

    def __init__(self, G, is_directed = False, p = 0.5, q = 0.5):
        self.is_directed = is_directed
        self.p = p
        self.q = q
        self.G = G

    def get_alias_edge(self, src, dst):

        G = self.G
        p = self.p
        q = self.q

        unnormalized_probs = []
        for dst_nbr in sorted(G.neighbors(dst)):
            if dst_nbr == src:
                unnormalized_probs.append(float(G[dst][dst_nbr]['weight']) / p)
            elif G.has_edge(dst_nbr, src):
                unnormalized_probs.append(float(G[dst][dst_nbr]['weight']))
            else:
                unnormalized_probs.append(float(G[dst][dst_nbr]['weight']) / q)
        norm_const = sum(unnormalized_probs)
        normalized_probs = [float(u_prob) / norm_const for u_prob in unnormalized_probs]

        return alias_setup(normalized_probs)

    def preprocess_transition_probs(self):

        G = self.G
        is_directed = self.is_directed

        alias_nodes = {}
        for node in G.nodes():
            unnormalized_probs = [float(G[node][nbr]['weight']) for nbr in sorted(G.neighbors(node))]
            norm_const = sum(unnormalized_probs)
            normalized_probs = [float(u_prob) / norm_const for u_prob in unnormalized_probs]
            alias_nodes[node] = alias_setup(normalized_probs)

        alias_edges = {}
        if is_directed:
            for edge in G.edges():
                alias_edges[edge] = self.get_alias_edge(edge[0], edge[1])
        else:
            time1 = time.clock()
            for index, edge in enumerate(G.edges()):
                if (index+1) % 10000 == 0:
                    time2 = time.clock()
                    print("{}th edge -> {}s".format(index + 1, time2 - time1))
                    time1 = time2

                alias_edges[edge] = self.get_alias_edge(edge[0], edge[1])
                alias_edges[(edge[1], edge[0])] = self.get_alias_edge(edge[1], edge[0])

        self.alias_nodes = alias_nodes
        self.alias_edges = alias_edges

    def node2vec_walk(self, walk_length, start_node):
        '''
        Simulate a random walk starting from start node.
        '''
        G = self.G
        alias_nodes = self.alias_nodes
        alias_edges = self.alias_edges

        walk = [start_node]

        while len(walk) < walk_length:
            cur = walk[-1]
            cur_nbrs = sorted(G.neighbors(cur))
            if len(cur_nbrs) > 0:
                if len(walk) == 1:
                    walk.append(cur_nbrs[alias_draw(alias_nodes[cur][0], alias_nodes[cur][1])])
                else:
                    prev = walk[-2]
                    next = cur_nbrs[alias_draw(alias_edges[(prev, cur)][0],
                                               alias_edges[(prev, cur)][1])]
                    walk.append(next)
            else:
                break

        return walk

    def simulate_walks(self, num_walks, walk_length):

        G = self.G
        walks = []
        nodes = list(G.nodes())
        print("walk iteration:")
        for walk_iter in range(num_walks):
            print("{} -> {}".format(str(walk_iter + 1), str(num_walks)))
            random.shuffle(nodes)
            for node in nodes:
                walks.append(self.node2vec_walk(walk_length=walk_length, start_node=node))

        return walks


def alias_setup(probs):
    K = len(probs)
    q = np.zeros(K)
    J = np.zeros(K, dtype=np.int)

    smaller = []
    larger = []

    for kk, prob in enumerate(probs):
        q[kk] = K * prob
        if q[kk] < 1.0:
            smaller.append(kk)
        else:
            larger.append(kk)

    while len(smaller) > 0 and len(larger):
        small = smaller.pop()
        large = larger.pop()

        J[small] = large
        q[large] = q[large] + q[small] - 1.0
        if q[large] < 1.0:
            smaller.append(large)
        else:
            larger.append(large)

    return J, q


def alias_draw(J, q):
    K = len(J)

    kk = int(np.floor(np.random.rand() * K))
    if np.random.rand() < q[kk]:
        return kk
    else:
        return J[kk]

def read_graph(edges, weighted = False, directed = False):

    G = nx.Graph()

    if weighted:
        G.add_weighted_edges_from(edges)
    else:
        for edge in edges:
            G.add_edge(edge[0], edge[1], weight = 1.0)

    if directed:
        G = G.to_directed()

    return G

def net_embedding(walks, output_file):

    walks = [list(map(str, walk)) for walk in walks]
    model = Word2Vec(walks, size=128, window=10, min_count=0, sg=1, workers=8, iter=1)
    model.wv.save_word2vec_format(output_file)


def random_walk(edges, num_walks = 100, walk_length = 160):

    time1 = time.clock()
    nx_G = read_graph(edges, weighted=True, directed=False)
    print("{} nodes, {} edges.".format(len(nx.nodes(nx_G)), len(nx.edges(nx_G))))

    time2 = time.clock()
    print("It cost {}s to read edges.".format(time2 - time1))

    G = Graph(nx_G, p=1, q=1)
    print("generate transition matrix......")
    G.preprocess_transition_probs()
    time3 = time.clock()
    print("It cost {}s to generate transition matrix.".format(time3 - time2))

    print("begin to random walk......")
    walks = G.simulate_walks(num_walks = num_walks, walk_length = walk_length)
    time4 = time.clock()
    print("It cost {}s to random walk.".format(time4 - time3))

    FileUtil.write2DemList2File(walks, "./node2vec_walks.txt")

    return walks

if __name__ == '__main__':


    pass
