# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         MultiLayerSim
# Description:
# Paper: An Ontology-Independent Representation Learning for Similar Disease Detection Based on Multi-layer Similarity Network
# Author:       JiashuaiZhang
# Date:         2020/5/20
#-------------------------------------------------------------------------------

from gensim.models import Word2Vec
import networkx as nx
import random
import numpy as np
from collections import defaultdict
from fastdtw import fastdtw
import multiprocessing
import math
import tempfile
import time

# --------------------------------------------------------
import warnings
warnings.filterwarnings(action='ignore')

# --------------------------------------------------------
from Util import FileUtil
from Util import NetUtil
from Common import common


def cal_path_sim(disease, edges, save_path_sim = False):

    print("begin calculate similarity based on path...")
    pathway = list(NetUtil.getColNodes(edges, col=1))

    ajaMatrix = np.zeros((len(disease), len(pathway)))
    for line in edges:
        row_index = disease.index(line[0])
        col_index = pathway.index(line[1])
        ajaMatrix[row_index][col_index] = float(1)

    W = np.dot(ajaMatrix, ajaMatrix.T)

    print("construct similarity matrix...")
    pathSim = {}
    sim_matrix = np.zeros((len(disease), len(disease)))
    for i in range(0, len(disease)):
        for j in range(i+1, len(disease)):
            # if W[i][j] != 0:
            pathSim["{}\t{}".format(disease[i], disease[j])] = 2*W[i][j]/(W[i][i]+W[j][j])
            sim_matrix[i][j] = 2 * W[i][j] / (W[i][i] + W[j][j])
            sim_matrix[j][i] = 2 * W[i][j] / (W[i][i] + W[j][j])

    if save_path_sim:
        print("sort the path similarity and save...")
        res = sorted(pathSim.items(), key= lambda x:x[1], reverse=True)
        FileUtil.writeSortedDic2File(res, "./path_sim.txt")

    return sim_matrix

def dtw_distance(s1, s2, eta = 0.5):
    '''
       计算两个度序列之间的距离
       :param a: list，表示一个有序的度序列
       :param b: list，表示一个有序的度序列
       :param eta: 表示一个计算序列间距离的参数
       :return:
       '''
    r, c = len(s1), len(s2)
    D0 = np.zeros((r + 1, c + 1))
    D0[0, 1:] = np.inf
    D0[1:, 0] = np.inf
    D1 = D0[1:, 1:]

    # 生成原始距离矩阵
    for i in range(r):
        for j in range(c):
            D1[i, j] = abs(s1[i] - s2[j])
    # print("原始距离矩阵:\n{}".format(D1))

    M = D1.copy()
    for i in range(r):
        for j in range(c):
            D1[i, j] += min(D0[i, j], D0[i, j + 1], D0[i + 1, j])
    # print("Cost Matrix或者叫累积距离矩阵:\n{}".format(D1))

    # print("the distance of two sequence: {}".format(D1[-1, -1]))
    # ------------------------------------------------------------------------
    dis = 0
    for i in range(r):
        for j in range(c):
            max_value = max([M[i, j], D1[i, j]])
            min_value = min([M[i, j], D1[i, j]])
            dis += (max_value + eta) /(min_value + eta) - 1
            # print(max_value, min_value, (max_value + eta) /(min_value + eta) - 1)

    return dis

def dtw_distance_fast(s1, s2):
    '''
       计算两个度序列之间的距离
       :param a: list，表示一个有序的度序列
       :param b: list，表示一个有序的度序列
       :return:
       '''

    dis, path = fastdtw(s1, s2)

    return dis


def cal_nei_sim(disease, edges, save_nei_sim = False):


    print("begin to calculate similarity based on neighbours...")

    G =  nx.Graph()
    G.add_edges_from(edges) # 将多种生物信息构造成异构矩阵

    print("step 1: epsilon -> 2, calculate first degree sequence and second degree sequence...")
    DegreeSequence1 = []
    DegreeSequence2 = []

    for di in disease:

        neighboursOne = G.neighbors(di) #获取节点的第一层邻居
        degreeOfOne = []
        neghboursTwo = []
        for indexOfNeighbours in neighboursOne:
            degreeOfOne.append(nx.degree(G, indexOfNeighbours)) #保存第一层邻居的degree
            neghboursTwo.extend(G.neighbors(indexOfNeighbours)) #获取第一层邻居节点的邻居
        sortedDegreeOfOne = sorted(degreeOfOne) #对第一层邻居的degree进行排序
        DegreeSequence1.append(sortedDegreeOfOne)

        neghboursTwo = set(neghboursTwo)
        neghboursTwo.remove(di) #去除二层邻居节点的自己

        degreeOfTwo = []
        for indexOfNeighbours in neghboursTwo:
            degreeOfTwo.append(nx.degree(G, indexOfNeighbours)) #保存第二层邻居的degree
        sortedDegreeOfTwo = sorted(degreeOfTwo) #对第一层邻居的degree进行排序
        DegreeSequence2.append(sortedDegreeOfTwo)

    cores = multiprocessing.cpu_count()  # 获取计算机CPU数目
    pool = multiprocessing.Pool(cores)  # 构造一个线程池
    print("step 2: compute neighbour_sim in parallel with {} cpus...".format(cores))

    # 构造一个多线程的任务
    resultsOne = [pool.apply_async(dtw_distance_fast, (DegreeSequence1[i], DegreeSequence1[j])) for i in
                  range(0, len(DegreeSequence1)) for j in range(i + 1, len(DegreeSequence1))]

    # 将成对的第一层degree sequence计算结果存储到数组中
    arrOne = np.zeros((len(DegreeSequence1), len(DegreeSequence1)))
    i = 0
    j = 1
    for r in resultsOne:
        if j == len(DegreeSequence1):
            i += 1
            j = i + 1
        arrOne[i][j] = float(r.get())
        j += 1

    # 构造一个多线程任务
    resultsTwo = [pool.apply_async(dtw_distance_fast, (DegreeSequence2[i], DegreeSequence2[j])) for i in
                  range(0, len(DegreeSequence2)) for j in range(i + 1, len(DegreeSequence2))]

    # 将成对的第二层degree sequence计算结果存储到数组中
    arrTwo = np.zeros((len(DegreeSequence2), len(DegreeSequence2)))
    i = 0
    j = 1
    for r in resultsTwo:
        if j == len(DegreeSequence2):
            i += 1
            j = i + 1
        arrTwo[i][j] = float(r.get())
        j += 1

    # ----------------------------------------------------------------------------
    print("step 3: construct similarity matrix...")
    alpha = 0.5  # a decaying weight factor α in the range between 0 and 1
    NeiSim = {}
    sim_matrix = np.zeros((len(disease), len(disease)))
    for i in range(0, len(disease)):
        for j in range(i + 1, len(disease)):
            distance = math.pow(alpha, 1) * arrOne[i][j] + math.pow(alpha, 2) * arrTwo[i][j]
            NeiSim["{}\t{}".format(disease[i], disease[j])] = math.exp(-distance)
            sim_matrix[i][j] = math.exp(-distance)
            sim_matrix[j][i] = math.exp(-distance)

    if save_nei_sim:
        print("sort the path similarity and save...")
        res = sorted(NeiSim.items(), key=lambda x: x[1], reverse=True)
        FileUtil.writeSortedDic2File(res, "./nei_Sim.txt")

    return sim_matrix

def con_single_layer_net(edges, filter_value = 0.1):

    disease, pathway = NetUtil.getNodes2HeterNet(edges)
    disease = list(disease)

    print("1st col -> {}\t 2nd col -> {}".format(len(disease), len(pathway)))

    nei_sim_matrix = cal_nei_sim(disease, edges, save_nei_sim=True)
    path_sim_matrix = cal_path_sim(disease, edges, save_path_sim=True)

    layer_net = []
    layer_net_result = {}
    for i in range(0, len(disease)):
        for j in range(i + 1, len(disease)):
            fusion_sim = 0.5 * nei_sim_matrix[i][j] + 0.5 * path_sim_matrix[i][j]
            if fusion_sim > filter_value:
                layer_net.append([disease[i], disease[j], fusion_sim])
                layer_net_result["{}\t{}".format(disease[i], disease[j])] = fusion_sim

    print("complete a single layer similarity network...")
    layer_net_result = common.sortDict(layer_net_result)
    FileUtil.writeSortedDic2File(layer_net_result, "./layer_net.txt")

    return layer_net

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


def random_walk(G, start_node, walk_length):
    '''

    :param G: nx.Graph，表示一个网络
    :param start_node: str，表示网络中的一个节点
    :param walk_length: int，表示该节点在网路中走的步数
    :return: list，表示随机游走时的路径
    '''

    walk = [start_node]
    while len(walk) < walk_length:
        cur = walk[-1]
        cur_nbrs = list(G.neighbors(cur)) #获取节点的邻居
        weights = [float(G[cur][cur_nbr]['weight']) for cur_nbr in cur_nbrs] #获取节点邻居对应的权值
        next = random.choices(cur_nbrs, weights, k = 1)[0] #根据权值选出下一步的节点
        walk.append(next)

    return walk

def random_walk_multi_layers(multi_layers_net, walk_iters = 100, walk_legth = 160, save_random_walk = False):

    '''
    :param multi_layers_net: list，表示一个多层的大网络，每一层都是一个小网络，并且每一层的节点都相同
    :param output_file: str，表示输出文件的路径，文件内容为节点及其对应的向量
    :param walk_iters: int，表示单个节点游走的次数
    :param walk_legth: int，表示一个节点在网络中游走的步数
    :param save_random_walk: boolean，表示是否保存随机游走的路径
    :return:
    '''
    multi_layers = defaultdict()
    for i in range(len(multi_layers_net)):
        G = read_graph(multi_layers_net[i], weighted=True)
        print("{} layer -> {} nodes and {} edges.".format(i, len(nx.nodes(G)), len(nx.edges(G))))
        multi_layers[i] = G

    nodes = multi_layers[0].nodes()
    walks = []
    for index, node in enumerate(nodes):
        time1 = time.clock()
        max_weights = list()
        for key, G in multi_layers.items():
            nei_weight = set()
            for nei_node in G.neighbors(node):
                nei_weight.add(float(G[node][nei_node]['weight']))
            max_weights.append(max(nei_weight))
        select_layer = max_weights.index(max(max_weights))
        for walk_iter in range(walk_iters):
            walk = random_walk(multi_layers[select_layer], node, walk_legth)
            walks.append(walk)
        time2 = time.clock()
        print("{} * {} -> {} layer: cost {}s".format(index, node, select_layer + 1, time2 - time1))


    if save_random_walk:
        FileUtil.write2DemList2File(walks, "./random_walk.txt")

    return walks


def calculateDisSim(walks, output_file, save_node_vectors = False):

    print("learn representations...")
    walks = [list(map(str, walk)) for walk in walks]
    model = Word2Vec(walks, size=128, window=10, min_count=0, sg=1, workers=8, iter=1)

    if save_node_vectors:
        temp_walk_fname = "./node_vectors.txt"
    else:
        _, temp_walk_fname = tempfile.mkstemp()

    print(temp_walk_fname)
    model.wv.save_word2vec_format(temp_walk_fname)

    node_vectors = defaultdict(list)
    with open(temp_walk_fname, 'r') as f:
        lines = f.readlines()
        for index in range(1, len(lines)):
            line = lines[index].strip().split(' ')
            vec = line[1:]
            vec = [float(i) for i in vec]
            node_vectors[line[0]] = vec  # key：label_disease的下标，value：对应的vectors

    dis_sim = {}
    disease = list(node_vectors.keys())
    for x in range(0, len(disease)):
        for y in range(x + 1, len(disease)):
            sim = common.cosinValue(node_vectors[disease[x]], node_vectors[disease[y]])
            if sim != 0:
                dis_sim["{}\t{}".format(disease[x], disease[y])] = sim

    FileUtil.writeDic2File(dis_sim, output_file)


if __name__ == '__main__':


    pass