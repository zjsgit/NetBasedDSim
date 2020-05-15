# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         MicrobeSim
# Description:
# Paper: An analysis of human microbe–disease associations
# Author:       JiashuaiZhang
# Date:         2020/5/14
#-------------------------------------------------------------------------------

from Common import common
from Util import FileUtil

from collections import defaultdict
import math
import numpy as np
import time

def calculateDisSim(disease_microbe, output_file):
    '''

    :param disease_microbe: 二维list，表示disease-microbe network
    :param output_file: str，表示结果保存的路径
    :return:
    '''

    begin_time = time.clock()
    microbe2disease = defaultdict(set)
    disease2microbe = defaultdict(set)

    for line in disease_microbe:
        microbe2disease[line[0]].add(line[1])
        disease2microbe[line[1]].add(line[0])

    print("there are {} diseases and {} microbes in disease-microbe.".format(
        len(disease2microbe.keys()), len(microbe2disease.keys())))

    diseases = list(disease2microbe.keys())
    microbes = list(microbe2disease.keys())
    weight = np.zeros((len(disease2microbe.keys()), len(microbe2disease)))
    E = np.ones((len(disease2microbe.keys()), len(microbe2disease)))

    for line in disease_microbe:
        indexRow = diseases.index(line[1])
        indexCol = microbes.index(line[0])

        weight[indexRow][indexCol] += 1
        if line[3] == "increase":
            E[indexRow][indexCol] = 1
        elif line[3] == "decrease":
            E[indexRow][indexCol] = -1

    for indexRow in range(0, len(diseases)):
        for indexCol in range(0, len(microbes)):
            # print math.log(diseaseNum / len(n.get(microbeList[indexCol], 2)))
            weight[indexRow][indexCol] *= E[indexRow][indexCol] * math.log2(float(len(diseases))/len(microbe2disease[microbes[indexCol]]))

    # ------------------------------------------------------------------
    MicrobeSim = {}
    for i in range(0, len(diseases)):
        for j in range(i + 1, len(diseases)):
            cosine_value = common.cosinValue(weight[i], weight[j])
            if cosine_value != 0:
                MicrobeSim["{}\t{}".format(diseases[i], diseases[j])] = cosine_value

    MicrobeSim = common.sortDict(MicrobeSim)
    FileUtil.writeSortedDic2File(MicrobeSim, output_file)
    end_time = time.clock()

    print("MicrobeSim costs {}s".format(end_time - begin_time))

    pass


if __name__ == '__main__':
    pass