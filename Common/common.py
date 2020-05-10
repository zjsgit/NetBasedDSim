# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         common
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from collections import defaultdict

def list2DictSet(lines, key = 1, value = 2):
    '''

    :param lines: 表示多行数据，每行数据有两列
    :param sep: 行数据切割符号
    :param key: 表示字典中key的列
    :param value: 表示字典中value中的列
    :return:
    '''

    dict_set = defaultdict(set)

    for line in lines:
        dict_set[line[key -1]].add(line[value-1])

    return dict_set

def cosinValue(vec1, vec2):
    '''
    计算两个疾病向量的cosine value
    :param vec1:
    :param vec2:
    :return:
    '''
    if len(vec1)!= len(vec2):
        print("error input, x and y is not in the same space")
        return

    result1 = 0
    result2 = 0
    result3 = 0
    for v in range(0, len(vec1)):
        result1 += vec1[v] * vec2[v]
        result2 += vec1[v] ** 2
        result3 += vec2[v] ** 2

    return result1/((result2 * result3) ** 0.5)

def sort_dict(dict_sim):
    '''
    对疾病对根据疾病相似性值进行排序
    :param dict_sim:
    :return:
    '''
    sorted_dis_sim = sorted(dict_sim.items(), key=lambda x: x[1], reverse=True)

    return sorted_dis_sim


if __name__ == '__main__':
    pass