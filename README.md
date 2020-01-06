#GeneCut
Gene expression analysis using Python and R Scripts
Contact the author with any questions: ltrinity@uvm.edu

1. Place a .csv file meeting the required specifications in the input folder

2. Run the R script fcrosAnalysis.R -> Check the temp folder for the four output files*
*R is set up by default to output the top N down/ug regulated genes
*This has been modified to output min and max significance values for all genes

3. Run the aggregate.py script -> Check the temp folder for the output file ChangesDict.p

4. Run the combineFcros.py script to incorporate the data from the R output into the ChangesDict.p dataframe

