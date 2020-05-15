# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         NetUtil
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/11
#-------------------------------------------------------------------------------

from collections import defaultdict

def labelNode(node, label):
    '''
    为一类节点加上标签，便于在异构网络中随机游走
    :param node: list，表示需要加标签的节点
    :param label: str，表示标签，通常是一个字母，如i，f，r等
    :return: dict，key：节点名字，value：节点标签+节点在list中的下标
    '''
    node_name = {} #表示加标签后的节点
    for index, value in enumerate(node):
        node_name[value.upper()] = label + str(index + 1)

    return node_name

def label2HeterogeneousNet(edges, lable1_node, label2_node):

    '''
    获取异构网络中加标签后的节点对应的邻居
    :param edges: 二维list，表示network中的边
    :param lable1_node: dict，表示加标签的一类节点。key：节点名字，value：加标签后节点名字
    :param label2_node: dict，表示加标签的另外一类节点。key：节点名字，value：加标签后节点名字
    :return:
    '''

    label1_node_neighbours = defaultdict(list) #表示网络中节点1的所有邻居
    label2_node_neighbours = defaultdict(list)  # 表示网络中节点2的所有邻居
    for edge in edges:
        try:
            label1_node_neighbours[lable1_node[edge[0].strip().upper()]].append(label2_node[edge[1].upper()])
            label2_node_neighbours[label2_node[edge[1].upper()]].append(lable1_node[edge[0].strip().upper()])
        except:
            pass

    return label1_node_neighbours, label2_node_neighbours


def label2HomogeneousNet(edges, lable_node):

    '''
    获取同构网络中加标签后的节点对应的邻居
    :param edge: 二维list，表示network中的边
    :param lable_node: dict，表示加标签的一类节点。key：节点名字，value：加标签后节点名字
    :return:
    '''
    lebel_node_neighbours = defaultdict(list) #表示网络中节点的所有邻居
    for line in edges:
        try:
            lebel_node_neighbours[lable_node[line[0].strip().upper()]].append(lable_node[line[1].upper()])
            lebel_node_neighbours[lable_node[line[1].upper()]].append(lable_node[line[0].strip().upper()])
        except:
            pass

    del_node = []
    for label_node, neighbours in lebel_node_neighbours.items():
        if len(neighbours) <= 2:
            del_node.append(label_node)

    for node in del_node:
        del lebel_node_neighbours[node]

    return lebel_node_neighbours

def getNodes2HeterNet(edges):
    '''
    获取异构网络中的两类节点
    :param edges: 二维list，表示一个网络的所有边
    :return:
    '''
    nodes1_set = set()
    nodes2_set = set()

    for edge in edges:
        nodes1_set.add(edge[0].upper())
        nodes2_set.add(edge[1].upper())

    return nodes1_set, nodes2_set


def getNodes2HomoNet(edges):

    '''
    获取同构网络中的所有节点
    :param edges: 二维list，表示一个网络的所有边
    :return:
    '''
    nodes_set = set() #表示网络中的节点

    for edge in edges:
        nodes_set.add(edge[0].upper())
        nodes_set.add(edge[1].upper())

    return nodes_set



if __name__ == '__main__':
    pass