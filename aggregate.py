# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 18:54:50 2019
Gene Analysis Script
@author: luket
"""
#import section
import pandas as pd
import numpy as np
import os
import pickle
fileName = os.listdir("input")[0]
#read in the data file
df = pd.read_csv("input\\" + fileName)
#define the column names to compare expression
identifier = df.columns[6]
measureStem = df.columns[2]
measureEC = df.columns[3]
measureNn = df.columns[4]
biologicalProcessGO = df.columns[8]
cellularComponentGO = df.columns[9]
molecularFunctionGO = df.columns[10]
alignment = df.columns[11]

#initialize a dictionary to store both types of gene expression differences
changesDict = {}
#for each row in the dataframe
for index, row in df.iterrows():
    myID = str(row[identifier]).replace(" ","")
    #if the gene code is in the dictionary
    if(myID in changesDict.keys()):
        #append the difference of endothelial cell expression as a value
        changesDict[myID]["EC-iPs"].append(row[measureEC]-row[measureStem])
        #append the difference of neuronal cell expression as a value
        changesDict[myID]["Nn-iPs"].append(row[measureNn]-row[measureStem])
    #if the gene code is not in the dictionary create the key/value pairs for EC and Nn
    else:
        changesDict[myID] = {"EC-iPs": [row[measureEC]-row[measureStem]],
                   "Nn-iPs": [row[measureNn]-row[measureStem]],"biologicalProcessGO": [row[biologicalProcessGO]],
                   "cellularComponentGO": [row[cellularComponentGO]],"molecularFunctionGO": [row[molecularFunctionGO]],
                   "alignment": [row[alignment]]}
#store the max and min fold changes and calculate the mean value of gene expression
counter = 0
for key in changesDict:
    counter +=1
    for condition in ["EC-iPs","Nn-iPs"]:
        values = changesDict[key][condition]
        changesDict[key]["minFold" + condition] = 2**min(values)
        changesDict[key]["maxFold"+condition] = 2**max(values)
        if not np.isnan(np.mean(values)):
            changesDict[key]["mean"+condition] = np.mean(values)   
print(counter)
#output the result from the file
pickle.dump( changesDict, open( "temp\\changesDict.p", "wb" ) )