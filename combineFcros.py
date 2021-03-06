# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:28:19 2019
Script to add the output from an fcross script to an existing dictionary
@author: luket
"""
import pickle
import pandas as pd
changesDict = pickle.load( open( "temp\\changesDict.p", "rb" ) )
#this method takes in an fcros file and returns dictionaries with the min and max significance value
def minMaxFcros(fCrosOutputFile, up, key):
    #grab the input file
    upDownRegData = pd.read_csv(fCrosOutputFile)
    #for each row in the input file
    counter = 0
    for index, row in upDownRegData.iterrows():
        #parse the gene symbol
        symbol = (row[0].split(" ")[0])
        symbol = symbol.strip('"')
        if symbol == '':
            break
        #parse the fcros value
        fVal = float(row[0].split(" ")[1])
        #store max/min
        minKey = 'min fcross significance (' + key+ ')'
        maxKey = 'max fcross significance (' + key+ ')'
        if minKey in changesDict[symbol].keys():
            if float(fVal) < changesDict[symbol][minKey][0]:
                changesDict[symbol][minKey] = ([fVal])
            if float(fVal) > changesDict[symbol][maxKey][0]:
                changesDict[symbol][maxKey] = ([fVal])
        else:
            counter+=1
            changesDict[symbol][minKey] = ([fVal])
            changesDict[symbol][maxKey] = ([fVal])
    print(counter)
locations = ["temp\\upregulatedECMaster.csv","temp\\upregulatedNnMaster.csv",
             "temp\\downregulatedECMaster.csv","temp\\downregulatedNnMaster.csv"]
upDown = [True, True, False, False]
keys = ['EC','Nn','EC','Nn']
for i in range(4):
    minMaxFcros(locations[i],upDown[i],keys[i])
pickle.dump( changesDict, open( "temp\\changesDict.p", "wb" ) )