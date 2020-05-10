# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         data_process
# Description:  对数据进行预处理，便于方法的集成
# 预处理的行为有： 统一格式，数据过滤等
# Author:       JiashuaiZhang
# Date:         2020/5/10
#-------------------------------------------------------------------------------

from Util import FileUtil

def line2upper(input_file, output_file):

    lines = FileUtil.readFile2List(input_file)
    new_lines = [line.upper() for line in lines]
    FileUtil.writeList2File(new_lines, output_file)



if __name__ == '__main__':

    # 将disease-miRNA中的数据都改写为大写格式
    # input_file1 = "./disease-miRNA.txt"
    # output_file1 = "./disease-miRNA2.txt"
    # line2upper(input_file1, output_file1)

    # 将disease-pathway中的数据都改写为大写格式
    # input_file1 = "./disease-pathway.txt"
    # output_file1 = "./disease-pathway2.txt"
    # line2upper(input_file1, output_file1)

    # 将disease-lncRNA中的数据都改写为大写格式
    input_file1 = "./disease-lncRNA.txt"
    output_file1 = "./disease-lncRNA2.txt"
    line2upper(input_file1, output_file1)


    pass