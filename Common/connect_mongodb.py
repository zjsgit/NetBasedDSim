# coding: utf-8

#-------------------------------------------------------------------------------
# Name:         connect_mongodb
# Description:  
# Author:       JiashuaiZhang
# Date:         2020/10/14
#-------------------------------------------------------------------------------

from pymongo import MongoClient

def get_mongodb_connection(database):
    '''
    根据数据库名返回该数据库的连接
    :param database: 要使用的数据库名称
    :return:
    '''
    client = MongoClient("localhost", 27017)
    db = client.get_database(database)
    return db


if __name__ == '__main__':

    pass