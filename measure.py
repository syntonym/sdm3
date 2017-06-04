import subprocess

def call(p):
    return ("java", "-classpath", ".:lib/commons-cli-1.4.jar", "KMeans", "-testdata", "LSH-nmi-adapted.csv", "-p", str(p))

for p in range(0, 1):
    process = subprocess.run(call(p), stdout=subprocess.PIPE)
    with open("output.raw", mode="ab") as f:
        f.write(process.stdout)
