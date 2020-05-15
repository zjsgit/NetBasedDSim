# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         Resink
# Description:
# Paper: Using information content to evaluate semantic similarity in a taxonomy
# Author:       JiashuaiZhang
# Date:         2020/5/14
#-------------------------------------------------------------------------------

from Util import NetUtil
from Util import FileUtil
from Common import common

from collections import defaultdict
import math
import networkx as nx
import numpy as np
import sys
import time



def getSectionoFromDic(keys, dic):
    '''
    根据keys取出整个dict中的一部分
    :param keys: set
    :param dic: dict，表示整个dict
    :return:
    '''
    newDic = {}
    for key in keys:
        newDic[key] = dic[key]
    return newDic


def getCommonAncesters(DG,node1,node2):
    '''
    获取两个节点在有向无循环图中的共同祖先
    :param DG: networkx.DirectedGraph，表示一个有向无循环图
    :param node1: str，表示一个节点
    :param node2: str，表示另一个节点
    :return:
    '''
    ancestorsOfGoIdA = nx.descendants(DG, node1)  # 获取GoIdA对应的祖先
    ancestorsOfGoIdB = nx.descendants(DG, node2)  # 获取GoIdB对应的祖先
    commonAncestors = ancestorsOfGoIdA & ancestorsOfGoIdB

    return commonAncestors


def calculateDisSim(DAG, output_file, disease_genes = None):
    '''

    :param DAG: 二维list，表示一个有向无循环图
    :param output_file: str，表示结果的存储路径
    :param disease_genes: 二维list，表示disease-gene associations
    :return:
    '''
    begin_time = time.clock()
    diseases = NetUtil.getNodes2HomoNet(DAG)
    disease_DAG = nx.DiGraph()
    disease_DAG.add_edges_from(DAG)

    # ----------------------------------------------------------------------------
    IC = defaultdict()
    if disease_genes:

        diseases_asso, genes = NetUtil.getNodes2HeterNet(disease_genes)
        disease2genes = common.list2DictSet(disease_genes, key= 1, value= 2)

        for di in diseases:
            if di in diseases_asso:
                IC[di] = - math.log2(float(len(disease2genes[di])) / len(genes))
            else:
                IC[di] = 0

        diseases = diseases & diseases_asso
        print("there are {} diseases for similarity based on DAG and associations.".format(len(diseases)))
    else:

        for di in diseases:
            descendants = nx.ancestors(disease_DAG, di)
            if descendants:
                IC[di] = - math.log2(float(len(descendants))/len(diseases))

        print("there are {} diseases for similarity based on DAG.".format(len(diseases)))
    # --------------------------------------------------------------------------------

    print("begin to calculate disease similarity......")

    diseases = list(diseases)
    simi_matrix = np.zeros((len(diseases), len(diseases)))

    for i in range(0, len(diseases)):

        sys.stdout.flush()
        temp_time1 = time.clock()
        di_A = diseases[i]
        for j in range(i + 1, len(diseases)):
            di_B = diseases[j]
            commonAncestors = getCommonAncesters(disease_DAG, di_A, di_B)
            sectionOfDoId2Gene = getSectionoFromDic(commonAncestors, IC)
            newDic = sorted(sectionOfDoId2Gene.items(), key=lambda x: x[1], reverse=True)
            if newDic:
                simi_matrix[i][j] = newDic[0][1]

        temp_time2 = time.clock()
        sys.stdout.write('\r{} -> {}, {}s'.format(i, diseases[i], (temp_time2 - temp_time1)))

    print()
    # ---------------------------------------------------------------------------------------------
    Resnik_simi = {}
    for i in range(0, len(diseases)):
        di_A = diseases[i]
        for j in range(i + 1, len(diseases)):
            di_B = diseases[j]
            if simi_matrix[i][j] > 0:
                Resnik_simi["{}\t{}".format( di_A,  di_B)] = simi_matrix[i][j]

    Resnik_simi = common.normalizeDict(Resnik_simi)
    FileUtil.writeDic2File(Resnik_simi, output_file)

    end_time = time.clock()
    print("ResnikSim costs {}s.".format(end_time - begin_time))

    pass




if __name__ == '__main__':

    FG = nx.DiGraph()
    FG.add_edges_from([(3, 2),
                       (2, 1),
                       (4, 1),
                       (5, 4),
                       (5, 6),
                       (6, 7),
                       (8, 7)
                       ]
                      )
    # print("节点{}的祖先为{}".format(7, nx.descendants(FG, 7)))
    print("节点{}的祖先为{}".format(5, nx.descendants(FG, 5)))
    # print("节点{}的子孙为{}".format(7, nx.ancestors(FG, 7)))
    des = nx.descendants(FG, 5)
    # des.add(5)
    H = FG.subgraph(list(des))

    print(list(H.edges))
    # for node in H.nodes:
    #     print("5 -> {} length is {}.".format(node, nx.shortest_path_length(H, 5, node)))
    # print(H.nodes)
    # print(nx.shortest_path_length(H, 4, 6))
    # print( nx.is_directed_acyclic_graph(H))

    pass