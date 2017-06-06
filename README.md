All figures are created from the data anylysis from Figures.ipynb. Data gathering was done in Analysis.ipynb.
The files are jupyter notebooks with code  written in python.

# Running the jar

You can run the jar via

`java -jar KMeans.jar`

The code we used to gather the data uses the un-jared version.
If you want to use our data gather scripts you must be able to compile, see the points below.

If you want to change some parameters simply type in `java -jar KMeans.jar -help` to see how that works.
example:

`java -jar KMeans.jar -p 7 -r true -width_auto true`

# Compiling the Java code

To compile the java code you need a JDK. This was tested against java1.8. If you have `java` in your PATH, and `make` is installed
you can compile by running `make compile` and run by `make run`. Otherwise see the Makefile for commandline parameters.

# Data Analysis Installation 

You need python > 3.3 and pip and the venv module. The code is tested with 3.6 but earlier versions should work.
To install python and pip use your OS software packaging e.g. `apt` or follow the instructions on python.org.

To run the code (instead of just viewing it) you need additional dependencies listed in `requirements.txt`.
For easy install execute `python3 setup.py` which will install all needed dependencies for you in a virtual environment.

To gather data you also need to have the java code compiled, see the bullet point above.
