** This is a Kmeans implementation for the third Scientific Data Management exercise **


1. Preprocess the csv file to drop non-numeric columns:

  * using bash:                 
    cut -f1-10 -d"," LSH-nmi.csv > LSH-nmi-corrected.csv

  * using excel:
    open csv with excel, drop last column, save as csv


2. Run the program 

  java -jar Kmeans.jar -i LSH-nmi-corrected.csv

  The default parameters are set correctly.
  If you want to play around, you can use -h for help with the program usage.

  The program writes the elapsed seconds to console output.
  
  (i) The maximum number of iterations is limited to 10. 
      A message if it converged or not is written to stderr.

(i) Default parameters:
  k = 15, iterations = 10, LLOYD algorithm, random partitioning of datapoints, dimensions = 10
