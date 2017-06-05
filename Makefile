compile:
	javac -classpath lib/commons-cli-1.4.jar KMeans.java

compile_debug:
	javac -g -classpath lib/commons-cli-1.4.jar KMeans.java

run: compile
	java -classpath .:lib/commons-cli-1.4.jar KMeans -testdata LSH-nmi-adapted.csv -p 5 -r true

debug: compile_debug
	jdb -classpath .:lib/commons-cli-1.4.jar KMeans -testdata LSH-nmi-adapted.csv -p 5 -r true
