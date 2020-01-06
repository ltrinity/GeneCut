# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 23:04:03 2020
This script generates the output csv files by subset or master as well as handles plotting of the data
@author: luket
"""
#this method takes in three parameters and outputs process csv and figures
#this first parameter is a boolean indicating if we are subsetting
#the second parameter is a string the name of the output file
#the third parameter is a list of all the gene ontology terms to include
def outputMasterOrSubset(includingByCategory, outFileName, geneOntologiesIncluded, APICall):
    #initialize the figure and draw the axis
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))
    ax1.axhline(y=0, xmin=-15, xmax=15,color='black',linestyle='-',zorder = 1)
    ax1.axvline(x=0, ymin=-15, ymax=15,color='black',linestyle='-',zorder = 1)
    #place circles
    for r in np.linspace(0,15,15):
        circ = plt.Circle((0, 0), radius=r, color='g',fill=False,zorder=1)
        ax1.add_artist(circ)
    ax1.set_xlabel('EC - iPS',fontsize=25)
    ax1.set_ylabel('Nn - iPS',fontsize=25)
    ax1.tick_params(labelsize=22)
    #parse title from outfile name
    plt.suptitle(outFileName[11:len(outFileName)-10] + ' gene expression changes by quadrant',fontsize=30)
    #generate the output dataframe
    outputDF = pd.DataFrame(columns=['Gene Symbol', 
                                     'GO Biological Process', 
                                     'GO Cellular Component',
                                     'GO Molecular Function','Quadrant', 
                                     'mean(EC-iPS)',
                                     'fcross FC (EC)',
                                     'min fcross significance (EC)', 
                                     'max fcross significance (EC)',
                                     'mean(Nn-iPS)',
                                     'fcross FC (Nn)',
                                     'min fcross significance (Nn)', 
                                     'max fcross significance (Nn)'])
    #for each key from the input file
    for key in meanEC:
        #python fold changes
        ECMaxfloatFoldChange = (ECmaxfoldDict[key])
        ECMinfloatFoldChange = (ECminfoldDict[key])
        NnMaxfloatFoldChange = (NnmaxFoldDict[key])
        NnMinfloatFoldChange = (NnminFoldDict[key])
        #initialize max min fcros values
        ECMaxfcros = float((ECMaxDict[key.replace(" ","")])[0])
        ECMinfcros = float((ECMinDict[key.replace(" ","")])[0])
        NnMaxfcros = float((NnMaxDict[key.replace(" ","")])[0])
        NnMinfcros = float((NnMinDict[key.replace(" ","")])[0])
        #baded on our criteria (see readme)
        if ECMaxfcros > 0.95 or ECMinfcros < 0.05 or NnMaxfcros > 0.95 or NnMinfcros < 0.05 or \
         float(NnMaxfloatFoldChange) >= 1.1479 or float(NnMinfloatFoldChange) <= 0.8521 \
         or float(ECMaxfloatFoldChange) >= 1.1479 or float(ECMinfloatFoldChange) <= 0.8521: 
                 #get the row of this gene
                rowIndex = (df.index[df['Gene Symbol'] == key].tolist())
                #store the ontology annotations
                bioOntology = (df.iloc[rowIndex[0]]['Gene Ontology Biological Process'])
                cellularOntology = (df.iloc[rowIndex[0]]['Gene Ontology Cellular Component'])
                molecOntology = (df.iloc[rowIndex[0]]['Gene Ontology Molecular Function'])
                #parse the annotations for consistency
                bioOntArray = str(bioOntology).replace('///','//').split(' // ')
                molecOntArray = str(molecOntology).replace('///','//').split(' // ')
                cellularOntArray = str(cellularOntology).replace('///','//').split(' // ')
                #assume we include (for the master file)
                include = True
                #track all gene ontology terms
                #if we are subsetting assume we do not include
                if includingByCategory:
                    include = False
                    #test for inclusion in subset across ontologies
                    if APICall:
                        for i  in range(0, len(bioOntArray)-1, 3):
                            if(bioOntArray[i] in geneOntologiesIncluded):  
                                include = True
                        for i  in range(0, len(molecOntArray)-1, 3):
                            if(molecOntArray[i] in geneOntologiesIncluded):
                                include = True
                        for i  in range(0, len(cellularOntArray)-1, 3):
                            if(cellularOntArray[i] in geneOntologiesIncluded):
                                include = True
                    else:
                        for i  in range(0, len(bioOntArray)-1, 3):
                            for string in geneOntologiesIncluded:
                                if(string in bioOntArray[i]):  
                                    include = True
                        for i  in range(0, len(molecOntArray)-1, 3):
                            for string in geneOntologiesIncluded:
                                if(string in molecOntArray[i]):
                                    include = True
                        for i  in range(0, len(cellularOntArray)-1, 3):
                            for string in geneOntologiesIncluded:
                                if(string in cellularOntArray[i]):
                                    include = True
                #if we are including the gene
                if include:
                    #store the symbol, scatter the point
                    gene = df.iloc[rowIndex[0]]['Gene Symbol']
                    ax1.scatter(meanEC[gene],meanNn[gene],zorder=2)
                    #initialize and determine quadrant
                    quadrant = -1
                    if meanEC[gene] > 0 and meanNn[gene] > 0:
                        quadrant = 1
                    if meanEC[gene] < 0 and meanNn[gene] > 0:
                        quadrant = 2
                    if meanEC[gene] < 0 and meanNn[gene] < 0:
                        quadrant = 3
                    if meanEC[gene] > 0 and meanNn[gene] < 0:
                        quadrant = 4
                    #form output row
                    rowDF = pd.DataFrame({'Gene Symbol':[gene],
                                          'GO Biological Process':[bioOntology],
                                          'GO Cellular Component':[cellularOntology],
                                          'GO Molecular Function':[molecOntology],
                                          'Quadrant':[quadrant], 
                                          'mean(EC-iPS)':[meanEC[gene]],
                                          'fcross FC (EC)':['[' + str(np.round(ECMinfloatFoldChange,3)) + ', '
                                                     + str(np.round(ECMaxfloatFoldChange,3)) + ']'],
                                          'min fcross significance (EC)':[ECMinfcros],
                                          'max fcross significance (EC)':[ECMaxfcros],
                                          'mean(Nn-iPS)':[meanNn[gene]],
                                          'fcross FC (Nn)': ['[' + str(np.round(NnMinfloatFoldChange,3)) + ', '
                                                     + str(np.round(NnMaxfloatFoldChange,3)) + ']'],
                                          'min fcross significance (Nn)':[NnMinfcros],
                                          'max fcross significance (Nn)':[NnMaxfcros],
                                          })
                    outputDF=outputDF.append(rowDF)
    #if we are subsetting output the figure
    if includingByCategory:
        quadrantTextXLocations = [1,-5,-5,1]
        quadrantTextYLocations = [2,2,-1.2,-1.2]
        countArray = [0,0,0,0]
        for i in range(1,5):
            quadDF = outputDF.where(outputDF['Quadrant'] == i)
            quadDF = quadDF.dropna()
            quadDF.to_csv(outFileName[:len(outFileName)-4] + "Quadrant" + str(i) + outFileName[len(outFileName)-4:],index=False)
            print(outFileName[:len(outFileName)-4] + "Quadrant" + str(i) + outFileName[len(outFileName)-4:] + " - " + str(quadDF.count()['Gene Symbol']))
            ax1.annotate("Q" + str(i) + ": " + str(quadDF.count()['Gene Symbol']),(quadrantTextXLocations[i-1],quadrantTextYLocations[i-1]),fontsize=20,bbox=dict(edgecolor = 'black',facecolor='white', alpha=1),
                         zorder=3)
            countArray[i-1] = quadDF.count()['Gene Symbol']
        # Data to plot
        labels = ["Q1","Q2","Q3","Q4"]
        sizes = countArray
        colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
        # Plot
        ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=350)
        ax2.set_xlabel('Proportion',fontsize=25)
        plt.axis('equal')
        plt.savefig(outFileName[0:11] + "Figures\\" + outFileName[11:len(outFileName)-4] + "Scatterplot.png",transparent=False,dpi=100,bbox_inches = "tight")
        plt.show()
    else:
        outputDF.to_csv(outFileName,index=False)
        print(outFileName  + " - " + str(outputDF.count()['Gene Symbol']))
#outputMasterOrSubset(True,outNames[0], signalGO, False)
#outputMasterOrSubset(True,outNames[1], transcriptionGO, False)
#outputMasterOrSubset(True,outNames[2], generateSubsetOntologies(urls[2]), True)
outputMasterOrSubset(True,outNames[3], metabolismGO , False)
outputMasterOrSubset(True,outNames[4], adhesionGO, False)
#outputMasterOrSubset(True,outNames[5], extracellularGO, False)
#outputMasterOrSubset(False,'outputdata\\ExpressionAnalysisMaster.csv', {}, False)
