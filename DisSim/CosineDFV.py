# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         CosineDFV
# Description:
# Paper:
# Author:       JiashuaiZhang
# Date:         2020/5/15
#-------------------------------------------------------------------------------

from Util import FileUtil

from Common import common
import numpy as np
import math
import time

def calculateDisSim(disease_symptom, output_file):

    begin_time = time.clock()
    disease2symptom = common.list2DictSet(disease_symptom, key=2, value=1)
    symptom2disease = common.list2DictSet(disease_symptom, key=1, value=2)

    diseases = list(disease2symptom.keys())
    symptoms = list(symptom2disease.keys())
    print("there are {} diseases and {} symptoms in disease-symptom".format(len(diseases), len(symptoms)))

    # --------------------------------------------------------------------------------------

    weight = np.zeros((len(symptoms), len(diseases)))
    for line in disease_symptom:
        row_index = symptoms.index(line[0])
        col_index = diseases.index(line[1])

        weight[row_index][col_index] = float(line[2]) * math.log2(
            float(len(diseases))/ len(symptom2disease[line[0]]))

    weight = weight.transpose()
    # ----------------------------------------------------------------------

    CosineDFV_sim = {}
    for i in range(0, len(diseases)):
        temp_time1 = time.clock()
        for j in range(i+1,  len(diseases)):
            cosine_value = common.cosinValue(weight[i], weight[j])
            if cosine_value != 0:
                CosineDFV_sim["{}\t{}".format(diseases[i], diseases[j])] = cosine_value

        temp_time2 = time.clock()
        print("{} -> {} costs {}s".format(i, diseases[i], temp_time2 - temp_time1))


    FileUtil.writeDic2File(CosineDFV_sim, output_file)
    end_time = time.clock()

    print("CosineDFV costs {}s".format(end_time - begin_time))



    pass

if __name__ == '__main__':
    pass