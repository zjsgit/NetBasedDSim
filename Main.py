# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         Main
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from DisSim import IDN
from Util import FileUtil
from DisSim import ModuleSim
from DisSim import NetSim
from DisSim import FunSim
from DisSim import mpDisNet
from DisSim import ResinkSim
from DisSim import MicrobeSim
from DisSim import CosineDFV
from DisSim import XuanSim
from Evaluation import benchmark_evaluation

def DisSim_by_IDN():

    input_file = "./Dataset/disease-gene.txt"
    output_file = "./Result/IDN_gene.txt"
    lines = FileUtil.readFile2List(input_file)
    dis_sim = IDN.calculateDisSim(lines)
    FileUtil.writeDic2File(dis_sim, output_file)

def DisSim_by_ModuleSim():

    input_file1 = "./Dataset/disease-gene.txt"
    input_file2 = "./Dataset/PPI.txt"
    output_file = "./Result/ModuleSim_sim.txt"

    dis2gene = FileUtil.readFile2DictSet(input_file1)
    gene_net = ModuleSim.read_interactome(input_file2, False, False)
    print("number of vertices: {}, number of edges: {}".format(gene_net.vcount(), gene_net.ecount()))

    sims = ModuleSim.similarity_cal_spavgn(dis2gene, gene_net)
    FileUtil.write_sims(sims, output_file)

def DisSim_by_NetSim():
    input_file1 = "./Dataset/disease-gene_SIDD.txt"
    input_file2 = "./Dataset/HumanNet_symbol.txt"
    output_file = "./Result/NetSim_sim_DO.txt"

    disease2gene = FileUtil.readFile2DictSet(input_file1, header=True)
    gene_gene = FileUtil.readFile2List(input_file2)
    NetSim.calculateDisSim(disease2gene, gene_gene, output_file)


def DisSim_by_FunSim():

    input_file1 = "./Dataset/disease-gene_SIDD.txt"
    input_file2 = "./Dataset/HumanNet_symbol_weighted.txt"
    output_file = "./Result/FunSim_sim.txt"

    disease2genes = FileUtil.readFile2DictSet(input_file1)
    weighted_PPI = FileUtil.readFile2List(input_file2)

    FunSim.calculateDisSim(disease2genes, weighted_PPI, output_file)

    # ------------------------ evaluation ---------------------------------------
    file_path1 = "./Evaluation/benchmark_DOID.txt"
    BenChmark_DO = FileUtil.readFile2List(file_path1)

    simi_Result = FileUtil.readFile2List(output_file)
    benchmark_evaluation.evaluate_by_benchmark(BenChmark_DO, simi_Result, times=10)

def DisSim_by_mpDisNet():

    input_file1 = "./Dataset/mpDisNet/omim_dm.txt"
    input_file2 = "./Dataset/mpDisNet/omim_mg.txt"
    input_file3 = "./Dataset/mpDisNet/PPI.txt"


    disease_miRNA = FileUtil.readFile2List(input_file1)
    miRNA_gene = FileUtil.readFile2List(input_file2)
    gene_gene = FileUtil.readFile2List(input_file3)

    output_file = "./Result/mpDisNet_sim.txt"
    mpDisNet.calculateDisSim(disease_miRNA, miRNA_gene, gene_gene, output_file, path_type=2)

    # ------------------------ evaluation ---------------------------------------
    file_path1 = "./Evaluation/benchmark_MeSH_RADAR.txt"
    BenChmark_MeSH1 = FileUtil.readFile2List(file_path1)

    simi_Result = FileUtil.readFile2List(output_file)
    benchmark_evaluation.evaluate_by_benchmark(BenChmark_MeSH1, simi_Result, times=10)

def DisSim_by_Resink():

    input_file1 = "./Dataset/DO_DAG.txt"
    input_file3 = "./Dataset/disease-gene_SIDD.txt"
    output_file = "./Result/Resink_sim.txt"
    DO_DAG = FileUtil.readFile2List(input_file1, header= True)
    disease_genes = FileUtil.readFile2List(input_file3)

    ResinkSim.calculateDisSim(DO_DAG, output_file)

    # -------------------------------------------------------------------
    file_path1 = "./Evaluation/benchmark_DOID.txt"
    BenChmark_DO = FileUtil.readFile2List(file_path1)

    simi_Result = FileUtil.readFile2List(output_file)
    benchmark_evaluation.evaluate_by_benchmark(BenChmark_DO, simi_Result, times=10)

def DisSim_by_MicrobeSim():

    input_file1 = "./Dataset/disease-microbe.txt"
    output_file = "./Result/MicrobeSim_sim.txt"

    disease_microbe = FileUtil.readFile2List(input_file1, header= True)
    MicrobeSim.calculateDisSim(disease_microbe, output_file)


def DisSim_by_CosineDFV():

    input_file1 = "./Dataset/disease-symptom.txt"
    output_file = "./Result/CosineDFV_sim.txt"

    disease_symptom = FileUtil.readFile2List(input_file1, header= True)
    CosineDFV.calculateDisSim(disease_symptom, output_file)

def DisSim_by_XuanSim():

    input_file1 = "./Dataset/DO_DAG.txt"
    input_file2 = "./Dataset/disease-gene_SIDD.txt"
    output_file = "./Result/XuanSim_sim.txt"

    DO_DAG = FileUtil.readFile2List(input_file1, header=True)
    diseases2genes = FileUtil.readFile2DictSet(input_file2)
    XuanSim.calculateDisSim(DO_DAG, output_file, selected_diseases = set(diseases2genes.keys()))

    file_path1 = "./Evaluation/benchmark_DOID.txt"
    BenChmark_DO = FileUtil.readFile2List(file_path1)

    simi_Result = FileUtil.readFile2List(output_file)
    benchmark_evaluation.evaluate_by_benchmark(BenChmark_DO, simi_Result, times=10)


if __name__ == '__main__':

    DisSim_by_XuanSim()

    pass