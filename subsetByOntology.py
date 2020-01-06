# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 22:41:01 2020
This script outputs group genes based on their ontology
#API call is used to generate list of descendants of epigenetic genes
@author: luket
"""
import requests, sys
import json
import pickle
#this method generates a list of ontologies for use in subsetting from an API
#here we use a url to receive a list of descendants
def generateSubsetOntology(url):
    requestURL = url
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    responseBody = r.text
    responseDict = json.loads(responseBody)
    ontologyList = list()
    for i in range(0,(responseDict['numberOfHits'])):
        ontologyDict = responseDict['results'][i]['descendants']
        for ont in ontologyDict:
            ontologyList.append(ont[3:])
    return ontologyList
#urls generates for use in generate method
epigeneticRequestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A0040029/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
#fileNames for outputData
outNames = ['outputdata\\SignalingSubset.csv','outputdata\\TranscriptionSubset.csv','outputdata\\EpigeneticSubset.csv','outputdata\\MetabolicSubset.csv',
               'outputdata\\CellAdhesionSubset.csv' ,'outputdata\\ExtracelullarMatrixProteinSubset.csv']
#generated using online documentation
signalGO = ["0007165","0023033"]
transcriptionGO = ["0003700", "0000130", "0001071", "0001130", "0001131", "0001151", "0001199", "0001204"]
#replace explicit assignment with API call
#epigeneticGO = ["0040029"]
epigeneticGO = generateSubsetOntology(epigeneticRequestURL)
metabolismGO = ["0008152","0044236","0044710"]
adhesionGO = ["0007155","0098602"]
extracellularGO = ["0031012"]
ontologySubsets = [signalGO, transcriptionGO, epigeneticGO, metabolismGO, adhesionGO, extracellularGO]
pickle.dump( ontologySubsets, open( "temp\\ontology.p", "wb" ) )