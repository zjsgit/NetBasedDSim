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

def sortDict(dict_sim):
    '''
    对疾病对根据疾病相似性值进行排序
    :param dict_sim: dict，表示一个需要排序的字典，key：str，value：float
    :return:
    '''
    sorted_dis_sim = sorted(dict_sim.items(), key=lambda x: x[1], reverse=True)

    return sorted_dis_sim

def normalizeList(data_list, col = 1):
    '''
    对list按照某一列进行归一化处理
    :param data_list: 二维list
    :param col: int，表示要标准化的列。所选列必须为数值
    :return:
    '''
    nor_col = [float(x[col]) for x in data_list]
    min_value = min(nor_col)
    max_value = max(nor_col)
    new_nor_col = [(x-min_value)/(max_value - min_value) for x in nor_col]
    for i in range(0, len(new_nor_col)):
        data_list[i][col] = new_nor_col[i]

    return data_list

def normalizeDict(data_dict, over_value = 0):
    '''
    对字典按照value进行归一化处理
    :param data_dict: dict，表示要进行归一化的dict。value必须为数值
    :param  over_value: float，表示阈值。超过此值则保存，否则不保存
    :return:
    '''
    new_data_dict = {}
    values = list(data_dict.values())
    min_value = min(values)
    max_value = max(values)
    for key, value in data_dict.items():
        new_value = (float(value) - min_value) / (max_value - min_value)
        if new_value > over_value:
            new_data_dict [key] = new_value

    return new_data_dict


if __name__ == '__main__':
    pass