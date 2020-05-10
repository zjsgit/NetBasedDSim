# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         Main
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from DisSim import IDN
from Util import FileUtil

def DisSim_by_IDN():

    input_file = "./Dataset/disease-gene.txt"
    output_file = "./Result/IDN_gene.txt"
    lines = FileUtil.readFile2List(input_file)
    dis_sim = IDN.calculateDisSim(lines)
    FileUtil.writeDic2File(dis_sim, output_file)




if __name__ == '__main__':

    DisSim_by_IDN()

    pass