# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         XuanSim
# Description:
# Paper: Prediction of microRNAs associated with human diseases based on weighted k most similar neighbors.
# Inferring the human microRNA functional similarity and functional network based on microRNA-associated diseases.
# Author:       JiashuaiZhang
# Date:         2020/5/15
#-------------------------------------------------------------------------------

from Util import FileUtil
from Util import NetUtil

import time
import networkx as nx
import math
import sys

def getNodeDVByDF(sub_DAG, node, decay_factor = 0.5):
    '''
    根据衰退因子求取子图中节点的DV值
    :param sub_DAG: networkx.DirectGraph，表示一个有向无循环图的子图
    :param node: str，表示该图中的一个节点
    :param decay_factor: float，表示衰退因子。其直在0-1之间，默认为0.5。
    :return:
    '''

    DV = 0
    ancestors = nx.descendants(sub_DAG, node)
    for an in ancestors:
        path_length = nx.shortest_path_length(sub_DAG, node, an)
        DV += math.pow(decay_factor, path_length)

    return DV

def getNodeDVByIT(sub_DAG, node):

    '''
    根据information theory求取子图中节点的DV值
    :param sub_DAG: networkx.DirectGraph，表示一个有向无循环图的子图
    :param node: node: str，表示该图中的一个节点
    :return: float，表示一个节点的DV值
    '''
    DV = 0
    ancestors = nx.descendants(sub_DAG, node)
    for an in ancestors:
        path_length = nx.shortest_path_length(sub_DAG, node, an)
        DV += math.log2(float(path_length + 1))

    return DV



def calculateDisSim(DAG, output_file, selected_diseases = None):

    '''

    :param DAG: 二维list，表示一个有向无循环图
    :param output_file: str，表示结果的存储路径
    :param selected_diseases: set，表示需要计算相似性的疾病集合。默认为None
    :return:
    '''
    begin_time = time.clock()
    diseases = NetUtil.getNodes2HomoNet(DAG)
    print("there are {} diseases in DAG.".format(len(diseases)))

    disease_DAG = nx.DiGraph()
    disease_DAG.add_edges_from(DAG)
    # ----------------------------------------------------------

    DV_diseases_dict = {}
    for di in diseases:
        ancestors = nx.descendants(disease_DAG, di)
        if ancestors:
            DV_disease_dict = {}
            ancestors.add(di)
            sub_graph = disease_DAG.subgraph(list(ancestors))
            sub_graph_nodes = sub_graph.nodes
            for node in sub_graph_nodes:
                DV_disease_dict[node] = getNodeDVByIT(sub_graph, node)
            DV_diseases_dict[di] = DV_disease_dict
        else:
            DV_disease_dict = {}
            DV_disease_dict[di] = 1
            DV_diseases_dict[di] = DV_disease_dict

    # -----------------------------------------------------------------------------------
    if selected_diseases:
        diseases = list(selected_diseases & diseases)
    else:
        diseases = list(diseases)

    print("{} diseases are used to calculate similarity.".format(len(diseases)))
    XuanSim_sim = {}
    for i in range(0, len(diseases)):
        DV_disease_A = DV_diseases_dict[diseases[i]]
        for j in range(i+1, len(diseases)):
            DV_disease_B = DV_diseases_dict[diseases[j]]
            common_diseases = set(DV_disease_A.keys()) & set(DV_disease_B.keys())
            if common_diseases:
                common_DV = 0
                for di in common_diseases:
                    common_DV += DV_disease_A[di] + DV_disease_B[di]
                XuanSim_sim["{}\t{}".format(diseases[i], diseases[j])] = common_DV/(
                    DV_disease_A[diseases[i]] + DV_disease_B[diseases[j]])
        if i % 100 == 0:
            print("{}->{}".format(i, len(diseases)))
    FileUtil.writeDic2File(XuanSim_sim, output_file)
    end_time = time.clock()

    print("XuanSim costs {}s.".format(end_time - begin_time))

    pass


if __name__ == '__main__':
    pass