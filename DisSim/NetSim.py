# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         NetSim
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from Common import walker
import numpy as np
import time
import sys

from Util import FileUtil
from Util import NetUtil

def getCommonNodes(netNode,disease2Genes):
    """
    获取PPI network和disease-gene中共有的gene
    :param netNode: set，表示all nodes of a network
    :param disease2Genes: set，表示all genes related to a disease
    :return: list，all nodes are in network and realated with a disease
    """
    leaveList=[]
    for disease2Gene in disease2Genes:
        if disease2Gene in netNode:
            leaveList.append(disease2Gene)

    return leaveList


def calculateDisSim(seed_list, net, output_path):

    '''

    :param seed_list: dict，表示disease和其对应的genes。key：str，表示disease，value：set，表示genes
    :param net: 二维list，表示一个PPI网络
    :param output_path: str，表示保存结果的路径
    :return:
    '''
    nodes = NetUtil.getNodes2HomoNet(net)

    print("there are {} diseases.".format(len(seed_list)))

    FRValueMatrix = np.zeros((len(seed_list), len(seed_list)))
    rowOfFR = 0
    time1 = time.clock()
    for disease, genesOfDisease in seed_list.items():

        leavaList = getCommonNodes(nodes, genesOfDisease)
        wk = walker.Walker(net)

        if len(leavaList) > 0:
            # run RWR(Random walk and restart),then get the proportion of all nodes
            temp_time1 = time.clock()
            nodesPercent = wk.run_exp(leavaList, 0.7, 1)
            temp_time = time.clock()
            print("{} - {} -> genes = {}, it cost {}s".format(rowOfFR, disease, len(leavaList), temp_time - temp_time1))

            # calculate the FR_GeneSet's value
            colOfFR = 0
            for disease2, genesOfDisease2 in seed_list.items():
                FR = 0
                for gene in genesOfDisease2:
                    if gene in nodes:
                        FR += float(nodesPercent[gene])
                    elif gene in genesOfDisease:
                        FR += 1
                    else:
                        FR += 0
                FRValueMatrix[rowOfFR][colOfFR] = FR
                colOfFR += 1

        rowOfFR += 1

    print("begin to calculate NetSim value")
    # calculate the NetSim value of a pair of disease
    NetSimMatrix = np.zeros((len(seed_list), len(seed_list)))
    NetSimList = list(seed_list.keys())
    rowOfFP = 0
    for disease, genesOfDisease in seed_list.items():

        colOfFP = 0
        diseseGeneNum = len(genesOfDisease)
        for disease2, genesOfDisease2 in seed_list.items():
            if disease is not disease2:
                disease2GeneNum = len(genesOfDisease2)
                NetSimMatrix[rowOfFP][colOfFP] = (FRValueMatrix[rowOfFP][colOfFP] +
                            FRValueMatrix[colOfFP][rowOfFP]) / (diseseGeneNum + disease2GeneNum)
            colOfFP += 1
        rowOfFP += 1

    print("write the 'disease-diesae-value' to a file")
    simiResult = {}
    row, col = np.shape(NetSimMatrix)
    for i in range(0, row):
        for j in range(i + 1, col):
            if NetSimMatrix[i][j] > 0:
                simiResult['{}\t{}'.format(NetSimList[i], NetSimList[j])] = NetSimMatrix[i][j]

    sortedSimiResult = sorted(simiResult.items(), key=lambda x: x[1], reverse=True)
    FileUtil.writeSortedDic2File(sortedSimiResult,  output_path)
    print("end")
    time2 = time.clock()
    print("NetSim total cost {}s.".format(time2-time1))

    pass


if __name__ == '__main__':


    pass