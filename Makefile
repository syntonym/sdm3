compile:
	javac -classpath lib/commons-cli-1.4.jar KMeans.java

run:
	java -classpath .:lib/commons-cli-1.4.jar KMeans -testdata LSH-nmi-corrected.csv
