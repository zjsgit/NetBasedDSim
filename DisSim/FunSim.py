# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         FunSim
# Description:
# Paper: SemFunSim: A New Method for Measuring Disease Similarity by Integrating Semantic and Gene Functional Association
# Author:       JiashuaiZhang
# Date:         2020/5/14
#-------------------------------------------------------------------------------

from Util import FileUtil

import networkx as nx
import numpy as np
import time

def calculateDisSim(disease2genes,gene_gene, output_file):
    '''

    :param disease2genes: dict，表示疾病及其对应的基因。key：disease，value：set，genes
    :param gene_gene: 二维lsit，表示一个加权PPI network。必须是加权网络
    :param output_file: str，表示结果写入文件
    :return:
    '''

    time_bigin = time.clock()
    diseases = list(disease2genes.keys()) #表示所有的疾病

    FG = nx.Graph()
    FG.add_weighted_edges_from(gene_gene)

    simi_matrix = np.zeros((len(diseases), len(diseases)))

    print("begin to calculate similarity of {} diseases.".format(len(diseases)))
    for i in range(0, len(diseases)):
        temp_time1 = time.clock()
        DoId2GenesA = disease2genes[diseases[i]]  # 取出i-th GoId对应的genes集合

        for j in range(i + 1, len(diseases)):
            DoId2GenesB = disease2genes[diseases[j]]  # 取出j-th GoId对应的genes集合

            # 计算FG(g)OfA的值
            Sum_valueOfFGA = 0
            for geneOfA in DoId2GenesA:
                if geneOfA in DoId2GenesB:  # 如果gene既在GensA集合中，又在GenesB集合中，则为1
                    Sum_valueOfFGA += 1
                else:
                    FGOfA = []
                    for geneOfB in DoId2GenesB:  # 如果geneOfA既在不GensB集合中，但是在HumanNet中，则取多组边的LLS的最大值
                        if FG.has_edge(geneOfA, geneOfB):
                            FGOfA.append(float(FG[geneOfA][geneOfB]['weight']))
                        else:
                            FGOfA.append(0)
                    Sum_valueOfFGA += float(max(FGOfA))

            # 计算FG(g)OfB的值
            Sum_valueOfFGB = 0
            for geneOfB in DoId2GenesB:
                if geneOfB in DoId2GenesA:  # 如果gene既在GensA集合中，又在GenesB集合中，则为1
                    Sum_valueOfFGB += 1
                else:
                    FGOfB = []
                    for geneOfA in DoId2GenesA:  # 如果geneOfB既在不GensA集合中，但是在HumanNet中，则取多组边的LLS的最大值
                        if FG.has_edge(geneOfA, geneOfB):
                            FGOfB.append(float(FG[geneOfA][geneOfB]['weight']))
                        else:
                            FGOfB.append(0)
                    Sum_valueOfFGB += float(max(FGOfB))

            simi_matrix[i][j] = (Sum_valueOfFGA + Sum_valueOfFGB) / (len(DoId2GenesA) + len(DoId2GenesB))

        temp_time2 = time.clock()
        print("{} -> {}, it costs {}s.".format(i, diseases[i], temp_time2-temp_time1))

    # ----------------------------------------将结果写入文件-----------------------------------------
    print("write similarity result to the file {}".format(output_file))
    # 将FunSim写入到文件中
    FunSimiResult = {}
    for i in range(0, len(diseases)):
        for j in range(i + 1, len(diseases)):
            if simi_matrix[i][j] != 0:
                FunSimiResult["{}\t{}".format(diseases[i], diseases[j])] = simi_matrix[i][j]

    result = sorted(FunSimiResult.items(), key=lambda x: x[1], reverse=True)

    FileUtil.writeSortedDic2File(result, output_file)
    time_end = time.clock()
    print("FunSim costs {}s totally".format(time_end - time_bigin))





if __name__ == '__main__':

    # WG = nx.Graph()
    # edges = [['1',"2", "0.97"], ["4", "5", "0.32"]]
    # WG.add_weighted_edges_from(edges)
    # print(WG.edges())
    #
    # for u, v, weight in WG.edges.data('weight'):
    #     if weight is not None:
    #         print(u, v, weight)


    pass