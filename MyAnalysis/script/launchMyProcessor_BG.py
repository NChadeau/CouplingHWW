#!/usr/bin/env python

import os
import sys

import MyProcessor_BG

if __name__ == "__main__":

    os.environ["MARLIN_DLL"] = '/home/chadeau/Documents/MyAnalysis/lib/libMyAnalysis.so'

    processID = 402010
    inputPath = "/home/chadeau/Documents/background/higgsstrahlung/eR.pL/n23n23h"

    fileList = [
        f'{inputPath}/{file}' for file in os.listdir(inputPath)
        if f'{processID}' in file
    ]
    fileList.sort()

    print('File list : ')
    print(fileList)

    params = MyProcessor_BG.Params()
    params.outputFileName = f'{processID}.root'
    params.maxRecordNumber = 0
    #params.skip = 0
    MyProcessor_BG.launch(params, fileList)
