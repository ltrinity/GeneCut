#GeneCut
Gene expression analysis using Python and R Scripts
Contact the author with any questions: ltrinity@uvm.edu

1. Place a .csv file meeting the required specifications in the input folder

2. Run the R script fcrosAnalysis.R -> Check the temp folder for the four output files*
*R is set up by default to output the top N down/ug regulated genes
*This has been modified to output min and max significance values for all genes

3. Run the aggregate.py script -> Check the temp folder for the file ChangesDict.p

4. Run the combineFcros.py script to incorporate the data from the R output into the ChangesDict.p dataframe


5. Run the subsetByOntology.py script to generate the ontologies by which to subset genes
-> Check the temp folder for the file ontology.p

Filtering Requirements:
Fcros significance level < 0.05 or > 0.95 or
2**([Nn-Diff](normalized) - [HUV-iPS](normalized) >= 1.1479 or <= 0.8521 or
2**([EC-Diff](normalized) - [HUV-iPS](normalized) >= 1.1479 or <= 0.8521

Files (Gene Count):
ExpressionAnalysisMaster (Total: 17,336)
Signaling Subsets (Q1:680, Q2:255, Q3:352, Q4:370) -> GO:0007165, GO:0023033 
Transcription Subsets (Q1:363, Q2:193, Q3:226, Q4:117) -> GO:0003700, GO:0000130, GO:0001071, GO:0001130, GO:0001131, GO:0001151, GO:0001199, GO:0001204 
Epigenetic Subsets (Q1:27, Q2:16, Q3:40, Q4:9) ->  Descendants of GO:0040029
Metabolism Subsets (Q1:218, Q2:109, Q3:220, Q4:100) -> GO:0008152, GO:0044236, GO:0044710 
Cell Adhesion Subsets (Q1:198, Q2:93, Q3:110, Q4:121) -> GO:0007155, GO:0098602 
Extracellular Matrix Protein Subsets (Q1:68, Q2:34, Q3:44, Q4:45) -> GO:0031012

