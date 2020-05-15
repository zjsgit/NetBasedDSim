# coding: utf-8

#-------------------------------------------------------------------------------
# Name:
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/4/6
#-------------------------------------------------------------------------------

'''
进行ROC评估
'''


from sklearn import metrics
import numpy as np
import pandas as pd
import random
from Util import FileUtil
import matplotlib.pyplot as plt
import os

'''
绘制ROC曲线
'''
def drawROC(FPRs, TPRs, labels):
    lw = 2
    plt.figure(figsize=(8, 5))
    colors = ['darkorange', 'gray', 'red', 'tan', 'lime','black','yellow','blue','darkred','salmon','lightblue']
    for i in range(0, len(FPRs)):

        roc_auc = metrics.auc(FPRs[i],TPRs[i])
        plt.plot(FPRs[i],TPRs[i], color=colors[i], lw=lw, label="{}'s area is {:.2f}".format(labels[i], roc_auc))  ###假正率为横坐标，真正率为纵坐标做曲线

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')  #标准线
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()


'''
将计算出来的SimilarityResult分为两部分，一部分是benchmark data，另一部分是随机抽取的部分
获取ROC曲线的横坐标和纵坐标
'''
def evaluate_by_benchmark(benchmark, lines, times = 50):
    selectResult = []  # 保存挑选出来的标准数据
    noSelect = []  # 保存除了标准数据的其他数据
    selectNum = len(benchmark)
    index = 0
    for simiEdge in lines:
        if index < selectNum:
            j = 0
            for benchEdge in benchmark:
                if benchEdge[0] in simiEdge and benchEdge[1] in simiEdge:
                    selectResult.append(simiEdge)
                    index += 1
                    break
                else:
                    j += 1
            if j == selectNum:
                noSelect.append(simiEdge)
        else:
            noSelect.append(simiEdge)

    print("similarity result总包含了benchmark{}行中的{}行".format(len(benchmark) , len(selectResult)))
    # =======================将随机数据集和benchmark数据集合并在一起 =================================
    trueNum = len(selectResult)
    not_benchmark_number = 700
    print("从非标准数据集中随机选出{}疾病对，和benchmark中{}疾病对，共{}。混合后排序，计算ROC".
          format(not_benchmark_number, trueNum, not_benchmark_number + trueNum))

    ROC = []

    # --------------------------------重复随机选择--------------------------------------------
    randomdata = range(1, len(noSelect))  # 给出随机的范围
    for i in range(0, times):
        randomlist = random.sample(randomdata, not_benchmark_number)
        selectRandomResult = []
        for number in randomlist:
            selectRandomResult.append(noSelect[number])
        selectRandomResult.extend(selectResult)
        # print("从非标准数据集中随机选出{}疾病对，和benchmark中{}疾病对，共{}。混合后排序，计算ROC".
        #       format(not_benchmark_number, trueNum, len(selectRandomResult)))
        # ----------------------------对合并组数据按照相似值进行排序 -----------------
        sortedResult = sorted(selectRandomResult, key=lambda x: float(x[2]), reverse=True)

        TPR = []
        FPR = []
        allTopNum = [10, 30, 70, 100, 150, 200, 300, 400, 500, 600, 700, len(sortedResult)]
        for topNum in allTopNum:

            # 计算四个标准值
            TP = 0
            index = 0
            for line in sortedResult:
                index += 1

                for benchEdge in benchmark:
                    if benchEdge[0] in line and benchEdge[1] in line:
                        TP += 1
                        break

                if index == topNum:
                    break

            FN = trueNum - TP
            FP = topNum - TP
            TN = not_benchmark_number - FP

            TPR.append(float(TP) / (TP + FN))
            FPR.append(float(FP) / (FP + TN))

            # print "TP--{}\tTN--{}\tFP--{}\tFN--{}".format(TP, TN, FP, FN)
            # print "TPR is {}".format(float(TP) / (TP + FN))
            # print "FPR is {}".format(float(FP) / (FP + TN))

            # print (TPR)
            # print (FPR)
        try:
            roc_auc = metrics.auc(FPR, TPR)
            ROC.append(roc_auc)
            print("{}->AUC={}".format(i, roc_auc))
        except:
            pass

    averageROC = sum(ROC) / len(ROC)
    print("循环{}次后，the average AUC is {}.-------------------".format(times, averageROC))

    # return TPR, FPR


def drawAverageAUC(averageAUCs, names):

    plt.figure(figsize=(8, 5))
    colors = ['darkorange', 'gray', 'red', 'tan', 'lime', 'black', 'yellow', 'blue', 'darkred', 'salmon', 'lightblue']
    selectColors = colors[0:len(names)]
    plt.bar(range(len(averageAUCs)), averageAUCs, color = selectColors, tick_label = names)
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    plt.xlabel('Methods')
    plt.ylabel('Average AUC')
    plt.title('the average AUC of different methods')
    # plt.legend(loc="lower right")
    plt.show()

def get_basic_info(file_path_benchmark, file_path_simi):
    '''
    比较相似性结果和benchmark dataset中疾病信息
    :param file_path_benchmark:
    :param file_path_simi:
    :return:
    '''

    benchmark_dataframe = pd.read_csv(file_path_benchmark, sep="\t", names=['d1','d2'])
    disease1 = set(benchmark_dataframe['d1']) and set(benchmark_dataframe['d2'])
    print("benchmark中共有{}个疾病，{}个相似疾病对".format(len(disease1), len(benchmark_dataframe)))


    simi_dataframe = pd.read_csv( file_path_simi, sep="\t", names=['d1','d2', 'score'])
    disease2 = set(simi_dataframe['d1']) and set(simi_dataframe['d2'])
    print("simi中共有{}个疾病，{}个相似疾病对".format(len(disease2), len(simi_dataframe)))
    # print(disease1)
    # print(disease2)
    print("benchmark dataset独有的疾病对为{}对".format(len(disease1 - disease2)))


def get_top_number_match(file_path_benchmark, file_path_simi, top_number = 500):

    '''
    选取similarity result的top number，看能否和benchmark dataset比较
    :param file_path_benchmark:
    :param file_path_simi:
    :param top_number:
    :return:
    '''

    benchmark_dataframe = pd.read_csv(file_path_benchmark, sep="\t", names=['d1', 'd2'])
    disease_pairs1 = benchmark_dataframe.values.tolist()

    simi_dataframe = pd.read_csv(file_path_simi, sep="\t", names=['d1', 'd2', 'score'])
    disease_pairs2 = simi_dataframe.values.tolist()

    top_match_number = 0
    if len(simi_dataframe) < top_number:
        print("top number is higher than the lines of similarity result.")
        os._exit(0)
    else:
        for i in range(0, top_number):
            simi_line = disease_pairs2[i]
            # print(simi_line)
            for j in range(0, len(disease_pairs1)):
                if simi_line[0] in disease_pairs1[j] and simi_line[1] in disease_pairs1[j]:
                    top_match_number += 1
                    break
            if i % 10000 == 0:
                print("top {} match {}.".format(i, top_match_number))


if __name__ == "__main__":
    file_path1 = "./benchmark_MeSH_RADAR.txt"
    file_path3 = "./benchmark_DOID.txt"

    file_path2 = "../Result/ModuleSim_sim.txt"
    file_path4 = "../Result/standard_Resnik_result.txt"
    file_path5 = "../Result/mpDisNet_sim.txt"
    file_path6 = "../Result/FunSim_sim.txt"


    # # get_basic_info( file_path1,  file_path3)
    # get_top_number_match(file_path1,  file_path2, top_number= 110000)

    # --------------------------------------读取以MeSH id为标签的标准集---------------------------------------
    BenChmark_MeSH1 = FileUtil.readFile2List(file_path3)

    simi_Result = FileUtil.readFile2List(file_path4 )
    evaluate_by_benchmark(BenChmark_MeSH1, simi_Result, times=10)

