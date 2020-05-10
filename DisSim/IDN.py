# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         IDN
# Description: 根据论文The integrated disease network实现该方法
#
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

import math
import numpy as np

from Common import common

def calculateDisSim(lines):

    disease2item = common.list2DictSet(lines, key = 1, value = 2)
    item2disease = common.list2DictSet(lines, key = 2, value = 1)

    disease = list(disease2item.keys())
    item = list(item2disease.keys())
    print("there are {} disease, {} item.".format(len(disease), len(item)))

    # -----------------------------------construct disease vectors----------------------------------------
    print("step 1: construct disease vectors...")
    vectors = np.zeros((len(disease), len(item)))

    for line in lines:
        row = disease.index(line[0])
        col = item.index(line[1])
        idf = math.log2(float(len(disease)) / len(item2disease[line[1]]))
        vectors[row][col] = idf

    # -----------------------------calculate disease similarity-----------------------------
    print("step 2: calculate disease similarity...")
    dis_sim = {}
    for di1 in range(0, len(disease)):
        for di2 in range(di1 + 1, len(disease)):
            cos_value = common.cosinValue(vectors[di1], vectors[di2])
            if cos_value > 0:
                dis_sim["{}\t{}".format(disease[di1], disease[di2])] = cos_value

    return dis_sim



if __name__ == '__main__':
    pass