# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         MimMinner
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/25
#-------------------------------------------------------------------------------

import xml.sax
import xml.dom.minidom
from collections import defaultdict
import sys
import nltk
from collections import OrderedDict
import math
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

from Common import connect_mongodb
from Util import FileUtil

def process_mesh(input_file):
    '''
    对MeSH desc_数字系列的文件进行处理，获取A C类对应的ID、同义词等信息
    :param input_file: str，表示xml文件的路径
    :return: tree2ynonyms，list，表示A C类对应的层次信息和同义词信息
    '''

    mesh_tree = defaultdict(list)
    mesh_synonyms = defaultdict(list)

    print("begin to read xml...")
    DOMTree = xml.dom.minidom.parse(input_file)
    collection = DOMTree.documentElement
    # if collection.hasAttribute("LanguageCode"):
    #     print("Root element: {}".format(collection.getAttribute("LanguageCode")))

    DescriptorRecords = collection.getElementsByTagName("DescriptorRecord")

    for index, descriptor in enumerate(DescriptorRecords):
        # if descriptor.hasAttribute("DescriptorClass"):
        #     print("DescriptorClass: {}".format(descriptor.getAttribute("DescriptorClass")))

        descriptorUI = descriptor.getElementsByTagName("DescriptorUI")[0]
        descriptorName = descriptor.getElementsByTagName("DescriptorName")[0].getElementsByTagName("String")[0]

        sys.stdout.write("\r{}-descriptorUI->{}, descriptorName->{}".format( index,
            descriptorUI.childNodes[0].data, descriptorName.childNodes[0].data))

        try:
            treeNumberList = descriptor.getElementsByTagName("TreeNumberList")[0].getElementsByTagName("TreeNumber")
            for treeNumber in treeNumberList:
                # print("treeNumber->{}".format(treeNumber.childNodes[0].data))
                lable = treeNumber.childNodes[0].data

                # 提取出标签为A或C的MeSH id，及其对应的同义词
                if str(lable).startswith("A") or str(lable).startswith("C"):
                    mesh_tree[descriptorUI.childNodes[0].data].append(treeNumber.childNodes[0].data)
                    # ----------------------------------------------------------------------------------------------
                    try:
                        conceptList = descriptor.getElementsByTagName("ConceptList")[0].getElementsByTagName("Concept")
                        for concept in conceptList:
                            termList = concept.getElementsByTagName("TermList")[0].getElementsByTagName("Term")
                            for term in termList:
                                # print("term->{}".format(term.getElementsByTagName("String")[0].childNodes[0].data))
                                mesh_synonyms[descriptorUI.childNodes[0].data].append(
                                    term.getElementsByTagName("String")[0].childNodes[0].data)
                    except:
                        pass
                    #-----------------------------------------------------------------------------------------------
        except:
            pass

    print()
    print("there are {} mesh_synonyms, {} mesh_tree.".format(len(mesh_synonyms), len(mesh_tree)))
    # FileUtil.writeDicList2File(mesh_tree, "./Result/mesh_tree.txt")
    # FileUtil.writeDicList2File(mesh_synonyms, "./Result/mesh_synonyms.txt")

    tree_id2synonyms = dict()
    for key, synonyms in mesh_synonyms.items():
        # 对所有Mesh term：分词，提取词干，处理大小写
        tokenizer = nltk.RegexpTokenizer(r"\w+(?:[-]\w+)*|'|[-.(]+|\S\w*")
        del_list = ['syndrom', 'disease', 'of', 'the', 'and', 'in', 'at', 'to', 'by', '.', ',', '(', ')', 's', "'", ';',
                    ':']
        filtered_synonyms = set()
        for synonym in synonyms:
            record = tokenizer.tokenize(synonym)
            record = sterm(record)
            record = [word.lower() for word in record]
            new_synonym = ""
            for word in record:
                if word not in del_list:
                    new_synonym += word + " "
            filtered_synonyms.add(new_synonym.strip())

        for treeNumber in mesh_tree.get(key):
            tree_id2synonyms[treeNumber] = list(filtered_synonyms)

    tree_id2synonyms = OrderedDict(sorted(tree_id2synonyms.items(), key = lambda x:x[0]))

    with open("./Result/tree_id2synonyms.txt", "w", encoding="utf-8") as file:
        for key, values in tree_id2synonyms.items():
            synonyms = "\t"
            synonyms = synonyms.join(values)
            file.write("{}\t{}\n".format(key, synonyms))

    return tree_id2synonyms

def get_phenotype_description(omim_id, db, col_name):
    '''
    根据OMIM id从数据库中获取对应的文本描述
    :param omim_id: 需要查询的OMIM id
    :param db: 连接Mongodb的数据库
    :param col_name: 集合的名字
    :return:
    '''
    collection = db[col_name]
    col = collection.find({'_id': omim_id})

    phenotype_description = ""
    for item in col: # 根据查询结果进行判断，看数据库中是否有该表型信息
        if item != "":
            if item.get('TX'):
                for tx, description in item['TX'].items():
                    for des in description:
                        phenotype_description += des + " "


            if item.get('CS'):
                for cs, description in item['CS'].items():
                    phenotype_description += cs + " "
                    for des in description:
                        phenotype_description += des + " "

            return phenotype_description.lstrip()

        else:
            print("no {} in mongodb".find(omim_id))
            return ""

def process_omim(omim_ids):
    '''
    获取要计算的omim ids的文本描述，并进行过滤
    :param omim_ids: 表示要计算相似性的OMOM id, list
    :return: dict，key：omim id；value：对omim的文本描述
    '''
    db = connect_mongodb.get_mongodb_connection("MimMinner")
    phenotypes_info = {} # 表示omim id到文本描述的映射

    tokenizer = nltk.RegexpTokenizer(r"\w+(?:[-]\w+)*|'|[-.(.[./]+|\S\w*")
    delete_list = ['syndrom', 'disease', 'of', 'the', 'and', 'in', 'at', 'to', 'by', '.', ',', '(', ')', 's', "'", ';',
                   ':', '[', ']', '/']

    for omim_id in omim_ids:
        phenotype_info = get_phenotype_description(omim_id[0], db, 'omim')
        if phenotype_info is not None:
            phenotype_info = tokenizer.tokenize(phenotype_info)
            phenotype_info = sterm(phenotype_info)
            phenotype_info = [word.lower() for word in phenotype_info if word not in delete_list]
            phenotypes_info[omim_id[0]] = phenotype_info

    print("there are {} omim ids in mongodb".format(len(phenotypes_info)))

    with open("./Result/phenotypes_info.txt", "w", encoding="utf-8") as file:
        for key, values in phenotypes_info.items():
            description = " "
            description = description.join(values)
            file.write("{}\t{}\n".format(key, description))

    return phenotypes_info

def sterm(record):
    '''
    提取词干，删除指定字符
    :param record:
    :return:
    '''
    delete_list = ['of', 'the', 'and', 'in', 'at', 'to', 'by', '.', ',', '(', ')', 's', "'", ';',':','/','[',']','-','%']

    recordlist = []
    for word in record:
        if word not in delete_list:
            recordlist.append(word)
    recordlist = [nltk.PorterStemmer().stem(word) for word in recordlist]

    return recordlist

def calculate_counter(treeid,treeID,mesh_count,isCalculate):
    treeidlist = treeid.split('.')
    treeidlen = len(treeidlist)
    sta_treeid = []
    childNodes = []
    actual_count = float(mesh_count[treeid])
    childContri = 0.0

    for stid in treeID:
        if stid.startswith(treeid) :
            sta_treeid.append(stid)

    for staid in sta_treeid:
        stalist = staid.split('.')
        if len(stalist) == (treeidlen+1):
            childNodes.append(staid)
    if len(childNodes) > 0:
        for childNode in childNodes:
            childContri += calculate_counter(childNode,sta_treeid,mesh_count,isCalculate)
        result = actual_count + childContri/len(childNodes)
    else:
        result = actual_count

    mesh_count[treeid] = result
    isCalculate[treeid] = True

    return result

def calculateDisSim(phenotypes_info, tree_id2synonyms, output_file):
    '''
    根据mesh和omim中的数据计算表型的相似性
    :param phenotypes_info: 表示omim数据库中的表型信息，dict；key：omim id，value：对该表型的描述
    :param tree_id2synonyms: 表示mesh tree的节点及对应的表型同义词，dict；key：tree id，value：表型同义词
    :param output_file: 表型相似性的存储路径
    :return:
    '''

    # -------------------------------计算actual acount----------------------------------
    print("-----actual acount-------")
    omim_ids = list(phenotypes_info.keys())
    tree_ids = list(tree_id2synonyms.keys())

    actual_count = np.zeros((len(omim_ids), len(tree_ids)))
    for index_omim, omim_id in enumerate(omim_ids):
        description = phenotypes_info[omim_id]
        for index_tree, tree_id in enumerate(tree_ids):
            synonyms = tree_id2synonyms[tree_id]
            for synonym in synonyms:
                actual_count[index_omim][index_tree] += description.count(synonym)

        sys.stdout.write("\r{}->{}".format(index_omim, omim_id))

    np.savetxt("./Result/actual_count.txt", actual_count, delimiter="\t", fmt="%d")

    # ---------------------根据hirarchy_count-----------
    print("\n-----hirarchy_count-------")
    hiera_count = actual_count
    for index_omim, omim_id in enumerate(omim_ids):
        is_calculate = OrderedDict()
        tree_id_count = OrderedDict()
        for index_tree, tree_id in enumerate(tree_ids):
            tree_id_count[tree_id] = hiera_count[index_omim][index_tree]
            is_calculate[tree_id] = False

        for tree_id in tree_ids:
            if is_calculate[tree_id] == False:
                calculate_counter(tree_id, tree_ids, tree_id_count, is_calculate)

        for index_tree, tree_id in enumerate(tree_id_count.keys()):
            if is_calculate[tree_id] == True:
                hiera_count[index_omim][index_tree] = tree_id_count[tree_id]

        sys.stdout.write("\r{}->{}".format(index_omim, omim_id))

    np.savetxt("./Result/hierachy_count.txt", hiera_count, delimiter="\t", fmt="%f")

    # ----------------------------------计算global weight------------------------------------
    print("\n-----global weight-------")
    gwc_global = OrderedDict()
    for tree_id in tree_ids:
        gwc_global[tree_id] = 0

    mostCount = []
    for index_omim, omim_id in enumerate(omim_ids):
        most_occur = 0
        for index_tree, tree_id in enumerate(tree_ids):
            meshNum = float(actual_count[index_omim][index_tree])
            if meshNum > most_occur:
                most_occur = meshNum
            if meshNum > 0:
                gwc_global[tree_id] += 1

        mostCount.append(most_occur)

    for key in gwc_global.keys():
        recordNum = gwc_global[key]
        if recordNum > 0:
            gwc_global[key] = math.log2(len(omim_ids) / recordNum)
        else:
            gwc_global[key] = 0.0

    # ------------------------------计算local weight------------------------------
    print("-----local weight-------")
    gwc = list(gwc_global.values())
    weight_count = np.zeros((len(omim_ids), len(tree_ids)))
    for index_omim, omim_id in enumerate(omim_ids):
        mf = mostCount[index_omim]
        cal_list = []
        for index_tree, tree_id in enumerate(tree_ids):
            cal_list.append(float(hiera_count[index_omim][index_tree]))
        gwc_cal = np.array(gwc) * np.array(cal_list)
        gwc_list = gwc_cal.tolist()
        cal_result = []
        for score in gwc_list:
            if score > 0:
                cal_result.append(0.5 + 0.5 * (score / mf))
            else:
                cal_result.append(score)
        for index_tree, tree_id in enumerate(tree_ids):
            weight_count[index_omim][index_tree] = cal_result[index_tree]

    np.savetxt("./Result/weight_count.txt", weight_count, delimiter="\t", fmt="%f")

    # ---------------------------------计算phetypes similarity---------------------------------
    print("-----phetypes similarity-------")
    # 根据cosine计算疾病的相似性
    similarity_socre = cosine_similarity(weight_count)
    similarity_result = {}
    for i in range(len(omim_ids)):
        for j in range(i+1, len(omim_ids)):
            if similarity_socre[i][j] != 0:
                similarity_result["{}\t{}".format(omim_ids[i], omim_ids[j])] = similarity_socre[i][j]

    # 对疾病相似性排序
    similarity_result = dict(sorted(similarity_result.items(),key = lambda x:x[1], reverse=True))
    FileUtil.writeDic2File(similarity_result, output_file)


if __name__ == "__main__":

    pass
