#!/usr/bin/env python

import os
import sys

import MyProcessor

if __name__ == "__main__":

    os.environ["MARLIN_DLL"] = '/home/chadeau/Documents/MyAnalysis/lib/libMyAnalysis.so'

    processID = 402008
    inputPath = "/home/chadeau/Documents/n1n1h/eR.pL"

    fileList = [
        f'{inputPath}/{file}' for file in os.listdir(inputPath)
        if f'{processID}' in file
    ]
    fileList.sort()

    print('File list : ')
    print(fileList)

    params = MyProcessor.Params()
    params.outputFileName = f'{processID}.root'
    params.maxRecordNumber = 0
    #params.skip = 0
    MyProcessor.launch(params, fileList)
