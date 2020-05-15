# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         mpDisNet2
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/12
#-------------------------------------------------------------------------------

import pandas as pd
import random
import tempfile
import os
from collections import defaultdict

from Common import common
from Util import FileUtil
from Util import NetUtil


def random_walk_disease(d1, dm, mg, gg, gm, md, type =1):
    '''
    按照特定的路径进行随机游走
    :param d1: 路径的起始节点
    :param dm:
    :param mg:
    :param gg:
    :param gm:
    :param md:
    :param type: int，表示使用的随机游走的方式。默认为1，即d-m-g-g-m-d；2，即d-m-g-m-d
    :return:
    '''
    if type == 1:
        m1 = random.choice(dm[d1])
        g1 = random.choice(mg[m1])
        g2 = random.choice(gg[g1])
        m2 = random.choice(gm[g2])
        d2 = random.choice(md[m1])
        return [d1, m1, g1, g2, m2, d2]
    else:
        m1 = random.choice(dm[d1])
        g = random.choice(mg[m1])

        m2 = random.choice(gm[g])
        d2 = random.choice(md[m2])

        return [d1, m1, g, m2, d2]


def calculateDisSim(dm_edge, mg_edge, gg_edge, output_file,
                    walk_length = 1000, path_type = 1, save_vectors = False):
    '''

    :param dm_edge: 二维list，表示disease-miRNA network
    :param mg_edge: 二维list，表示miRNA-gene network
    :param gg_edge: 二维list，表示gene-gene network
    :param output_file: str，表示结果保存的路径
    :param walk_length: int，表示随机游走中单个节点走的步数。默认为1000
    :param path_type: int，表示随机游走的特定路径。默认为1
    :param save_vectors: boolean，表示是否保存disease vectors。默认为False，表示不保存疾病特征
    :return:
    '''

    dm_d, dm_m = NetUtil.getNodes2HeterNet(dm_edge)
    d_name_id = NetUtil.labelNode(dm_d, "i")
    m_name_id = NetUtil.labelNode(dm_m, "f")

    print("有{}个disease加标签".format(len(d_name_id.keys())))
    print("有{}个miRNA加标签".format(len(m_name_id.keys())))
    # ------------------------------------------------------------


    gg_g = NetUtil.getNodes2HomoNet(gg_edge)
    mg_m, mg_g = NetUtil.getNodes2HeterNet(mg_edge)
    g_name = mg_g & gg_g
    g_name = ['a' + str(i) for i in g_name]

    print("有{}个gene加标签".format(len(g_name)))

    # -----------------------------------------------------------------------------------
    #
    dm = defaultdict(list)
    md = defaultdict(list)

    for line in dm_edge:
        if line[1].upper() in dm_m:
            dm[d_name_id[line[0].upper().strip()]].append(m_name_id[line[1].upper()])
            md[m_name_id[line[1].upper()]].append(d_name_id[line[0].upper().strip()])

    # ----------------------------------------------------------------------------------------

    gg = defaultdict(list)
    for line in gg_edge:
        g1 = 'a' + str(line[0])
        g2 = 'a' + str(line[1])
        if (g1 in g_name) & (g2 in g_name):
            gg[g1].append(g2)
            gg[g2].append(g1)

    gene_del = []
    for gene in gg.keys():
        if len(gg[gene]) == 2:
            gene_del.append(gene)
    for gene in gene_del:
        del gg[gene]

    g_name = list(gg.keys())
    # ---------------------------------------------------------------
    test_mg_m = set()
    mg = defaultdict(list)
    gm = defaultdict(list)

    for line in mg_edge:
        g = 'a' + str(line[1])
        if (line[0].upper() in dm_m) & (g in g_name):
            mg[m_name_id[line[0].upper().strip()]].append(g)
            gm[g].append(m_name_id[line[0].upper().strip()])

    print("mg中有{}个miRNA, 标签后{}个".format(len(test_mg_m), len(mg.keys())))

    # -------------------------------------------------------------------------------------

    print("mg中miRNA标签后{}个".format(len(mg.keys())))
    disease_list = []
    for d in dm.keys():
        ms = dm[d]
        for m in ms:
            if set(mg[m]) & set(g_name):
                disease_list.append(d)
    disease_list = list(set(disease_list))

    print("mg中miRNA标签后{}个".format(len(mg.keys())))
    # print("there are {} diseases, {} genes in random walk.".format(len(disease_list ), len(g_name )))
    print("there {} miRNA in mg, {} miRNA in dm".format(len(mg.keys()), len(md.keys())))

    total_walk = []
    ii = 0
    for disease in disease_list:
        jj = 0
        for i in range(walk_length):
            # print(str(ii) + ' ' + str(i))
            temp = [disease]
            for k in range(50):
                j = 1
                while ((j == 1) & (jj != walk_length)):
                    try:
                        temp2 = random_walk_disease(disease, dm, mg, gg, gm, md, path_type)
                        temp.extend(temp2[1:])
                        disease = temp2[-1]
                    except:
                        j = 1
                        jj += 1
                    else:
                        j = 0
                if jj == walk_length:
                    break
            total_walk.append(temp)
        ii += 1

    # random_walk_path = "./Result/inputomim2.txt"
    # with open(random_walk_path, 'w') as f:
    #     for lines in total_walk:
    #         for line in lines:
    #             f.write(line + ' ')
    #         f.write('\n')

    print("learn representations...")

    _, temp_walk_fname = tempfile.mkstemp()

    print(temp_walk_fname)
    with open(temp_walk_fname, 'w') as f:
        for walk in total_walk:
            for line in walk:
                f.write(line + ' ')
            f.write('\n')

    _, temp_node_vec_fname = tempfile.mkstemp()

    statement = "Common/metapath2vec++ -train {} -output {} -pp 1 -size 128 -window 7 -negative 5 -threads 32".format(
        temp_walk_fname, temp_node_vec_fname)

    print(statement)
    os.system(statement)

    print("\ncalculate disease similarity...")
    node_vectors_path = temp_node_vec_fname + ".txt"
    node_vectors = defaultdict(list)
    with open(node_vectors_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split(' ')
            if line[0].startswith('i'):
                vec = line[1:]
                vec = [float(i) for i in vec]
                node_vectors[line[0]] = vec  # key：label_disease的下标，value：对应的vectors

    if save_vectors:
        FileUtil.writeDicSet2File(node_vectors, "./Result/mpDisNet_node_vectors.txt")
        print()

    new_label_disease = {value: key for key, value in d_name_id.items()}
    dis_sim = {}
    for x in range(0, len(disease_list)):
        for y in range(x + 1, len(disease_list)):
            sim = common.cosinValue(node_vectors[disease_list[x]], node_vectors[disease_list[y]])
            if sim != 0:
                dis_sim["{}\t{}".format(new_label_disease[disease_list[x]],
                                        new_label_disease[disease_list[y]])] = sim

    FileUtil.writeDic2File(dis_sim, output_file)


if __name__ == '__main__':
    pass