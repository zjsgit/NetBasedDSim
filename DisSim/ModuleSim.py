# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         ModuleSim
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from igraph import Graph
import math
import time

def read_interactome(interactomefile, weights, directed):
    """
    read a file into an igraph object using igraph read_Ncol method
    :param interactomefile: a file path of a "ncol" file contains interactome
    :param weights: same as "weights" in igraph object construction method:
    True, False, "auto", "if_present"
    :param directed: True or False
    :return: an igraph object
    """
    with open(interactomefile, mode='r') as f:
        g = Graph.Read_Ncol(f, names=True, weights=weights, directed=directed)
        return g


def similarity_cal_spavgn(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """
    use average shortest path and normalize to cal
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gfilter: True or False, use graph to filter disease gene assos or not
    :param gncutoff: gene number cut off, when filter is True, only diseases whose
    number of associated genes in graph is no less than gncutoff will be calculated
    :param transformdistance: True or False, use transform distance or divide
    :return: a dict, (key-value: string-dict<string-float>)
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")
    dselfavg = {}
    for d in diseases:
        avgavg = 0.0
        for g in dgassos_new[d]:
            avgavg += sim_geneset2gene_avg(g, dgassos_new[d], sim_gene2gene)
        dselfavg[d] = avgavg/len(dgassos_new[d])

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[i]], sim_gene2gene)
            osim = (simsum / (len(dgassos_new[diseases[i]]) + len(dgassos_new[diseases[j]])))
            navg = (dselfavg[diseases[i]] + dselfavg[diseases[j]]) / 2

            result[diseases[i]][diseases[j]] = osim/navg
        print("---------------------------------------cost time:", str(time.time()-now))
    return result

def sim_geneset2gene_avg(g, gset, sim_gene2gene):
    sims = []
    for i in gset:
        sims.append(sim_gene2gene[g][i])
    if len(sims) > 0:
        return sum(sims)/len(sims)
    else:
        return 0.0

def sim_gene2gene_shortestpath(dgs, graph, transformdistance=True):
    nodenames = graph.vs['name']
    sps = graph.shortest_paths(source=nodenames, target=nodenames, weights=None, mode=3)
    result = {}
    for n in nodenames:
        result[n] = {}
    if transformdistance is True:
        for i in range(0, len(nodenames)):
            for j in range(i, len(nodenames)):
                result[nodenames[i]][nodenames[j]] = transformed_distance(sps[i][j])
                result[nodenames[j]][nodenames[i]] = result[nodenames[i]][nodenames[j]]
    else:
        for i in range(0, len(nodenames)):
            for j in range(i, len(nodenames)):
                result[nodenames[i]][nodenames[j]] = 1 / 2**(sps[i][j])
                result[nodenames[j]][nodenames[i]] = result[nodenames[i]][nodenames[j]]
    # supply
    dglefts = set(dgs).difference(set(nodenames))
    for dg in dglefts:
        result[dg] = {}
        for eg in dgs:
            result[dg][eg] = 0.0
        result[dg][dg] = 1.0
    dginters = set(dgs).intersection(set(nodenames))
    for dg in dginters:
        for eg in dglefts:
            result[dg][eg] = 0.0
    return result

def transformed_distance(shortestpathdis=0, a=1, b=1):
    return a*math.exp(-b*shortestpathdis)

if __name__ == '__main__':


    pass