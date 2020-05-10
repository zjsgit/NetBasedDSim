# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         FileUtil
# Description:  对文件进行读取和写入的操作
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------
import sys
import os
from collections import defaultdict

def readFile2DictSet(input_file, sep = "\t", header = False):
    '''

    :param input_file: 输入文件的路径， 文件中默认为disease\tbio-item格式
    :param sep:一行中的隔离形式，默认为\t
    :param header: 每一列的标题，默认为无
    :return: 返回一个字典，key为disease，value为bio-item的set
    '''
    print("the input file is {}.".format(input_file))
    try:
        with open(input_file, encoding= "utf-8", mode="r") as f:
            lines = f.readlines()
    except:
        print("the file path is wrong!!!")
        sys.exit(0)

    dict_set = defaultdict(set)
    for line in lines:
        if header:
            next(lines)
        else:
            line_info = line.strip().split(sep)
            dict_set[line_info[0]].add(line_info[1])

    return dict_set


def readFile2List(input_file, header = False, sep = "\t"):
    '''

    :param input_file: 输入文件的路径
    :param header: 每一列的标题，默认为无
    :param sep: 行信息切割符
    :return: 一个list
    '''
    print("the input file is {}.".format(input_file))
    try:
        with open(input_file, encoding = "utf-8", mode="r") as f:
            lines = f.readlines()
    except:
        print("the file path is wrong!!!")
        sys.exit(0)

    data_list = []
    for line in lines:
        if header:
            next(lines)
        else:
            data_list.append(line.strip().split(sep))

    return data_list


# -------------------------write section---------------------------------------
def writeList2File(dataset, output_file):
    '''

    :param dataset:
    :param output_file: the file path to save data
    :return:
    '''
    if os.path.isfile(output_file):
        print("output file is not a file!!!")
        sys.exit(0)
    else:
        file = os.path.split(output_file)
        if not os.path.exists(file[0]):
            print("mkdirs->{}".format(file[0]))
            os.makedirs(file[0])

        # 在文件夹没有问题的情况下，将数据写入文件
        f = open(output_file, mode="w+", encoding="utf-8")
        for line in dataset:
            f.write(str(line) + "\n")
        f.close()
        print("{} has been written.".format(output_file))


def writeSortedDic2File(sortedDic, output_file):
    '''

    :param sortedDic:
    :param output_file:
    :return:
    '''
    if os.path.isfile(output_file):
        print("output file is not a file!!!")
        sys.exit(0)
    else:
        file = os.path.split(output_file)
        if not os.path.exists(file[0]):
            print("mkdirs->{}".format(file[0]))
            os.makedirs(file[0])

    # -------------------------------------------------------------
    f = open(output_file, mode = "w+", encoding="utf-8")
    for line in sortedDic:
        f.write("{}\t{}\n".format(line[0], line[1]))
    f.close()
    print("{} has been written.".format(output_file))


def writeDic2File(dic, output_file):
    '''

    :param dic:
    :param output_file:
    :return:
    '''
    if os.path.isfile(output_file):
        print("output file is not a file!!!")
        sys.exit(0)
    else:
        file = os.path.split(output_file)
        if not os.path.exists(file[0]):
            print("mkdirs->{}".format(file[0]))
            os.makedirs(file[0])

    # ----------------------------------------------------------
    f = open(output_file, mode = "w+", encoding = "utf-8")
    for key, value in dic.items():
        f.write("{}\t{}\n".format(key, value))
    f.close()
    print("{} has been written.".format(output_file))


if __name__ == '__main__':
    pass