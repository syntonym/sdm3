package at.ac.univie.PRG3;

import java.io.File;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Random;

public class Main {

    String path = "LSH data.txt";
    double[] data;
    double[] centroids, tmpCentroids;
    int[] clusterIDs, numClusterPoints;
    int points;
    final int clusters = 15;
    final int dimension = 10;
    int dataHash[];
    int[] clusterHashes;

    double bucketSize = 30;

    int funcGroups = 4;
    int funcs = funcGroups * funcGroups;
    double[] fHash;


    public void loadData(String path) {
        try {
            String raw = new String(Files.readAllBytes(Paths.get(path)));
            String [] sp = raw.split("(,|\n)"); //split in one big continuous array
            data = new double[sp.length];
            centroids = new double[clusters * dimension]; //centroid array
            tmpCentroids = new double[clusters * dimension];
            points = sp.length /dimension; //number of points
            clusterIDs = new int[points];
            numClusterPoints = new int[clusters];
            fHash = new double[funcs * dimension]; //hash functions
            dataHash = new int[funcs*points]; //hash values of data
            clusterHashes = new int[funcs * clusters]; //hash value of clusters



            double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;

            //parse data to double and find min/max
            for (int i = 0; i< sp.length; ++i) {
                data[i] = Double.parseDouble(sp[i]);
                if(data[i] < min)
                    min = data[i];
                if(data[i] > max)
                    max = data[i];
            }


            Random rng = new Random();
            //starting centroids
            for(int i = 0; i < clusters * dimension; ++i) {
                centroids[i] = rng.nextDouble() * (max - min) + min;
            }

            //random hash functions
            for(int i = 0; i < funcs * dimension; ++i) {
                fHash[i] = (rng.nextDouble()-0.5) *2;
            }


        } catch (Exception e) {
            System.err.println("can't read file.");
            System.exit(1);
        }
    }

    private void saveData(String path) {
        File file = new File(path);
        FileOutputStream fos = null;
        try {
            if (!file.exists())
                file.createNewFile();

            fos = new FileOutputStream(file);
            for(int i = 0; i < points; ++i) {
                String s = "";
                for(int j = 0; j < dimension; ++j) {
                    s += data[i * dimension + j] + ",";
                }
                s += clusterIDs[i] + "\n";
                fos.write(s.getBytes());
            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if(fos != null) {
                try {
                    fos.flush();
                    fos.close();
                } catch (Exception ignored) {}
            }
        }
    }

    private double distance(int dataPoint, int centroid) {
        double tmp, dist = 0;
        for(int i = 0; i < dimension; ++i) {
            tmp = data[dataPoint + i] - centroids[centroid + i];
            dist += tmp*tmp;
        }
        return dist;
    }

    private void kmeans() {
        double minDist, tmp;
        int clusterID;
        boolean converged = false;

        for(int i = 0; i < tmpCentroids.length; ++i) {
            tmpCentroids[i] = 0;
        }

        while (!converged) {
            converged = true;

            //init tmp
            for(int i = 0; i < numClusterPoints.length; ++i)
                numClusterPoints[i] = 0;

            for (int i = 0; i < points; ++i) { //for each point
                //calc min distance
                minDist = distance(i * dimension, 0);
                clusterID = 0;
                for (int j = 1; j < clusters; ++j) {
                    tmp = distance(i * dimension, j * dimension);
                    if (tmp < minDist) {
                        minDist = tmp;
                        clusterID = j;
                    }
                }
                //if cluster changed
                if(clusterIDs[i] != clusterID) {
                    converged = false;
                    clusterIDs[i] = clusterID;
                }

                for(int k = 0; k < dimension; ++k) {
                    tmpCentroids[clusterID * dimension + k] += data[i * dimension + k];
                }
                ++numClusterPoints[clusterID];
            }


            //update centroids
            for(int i = 0; i < tmpCentroids.length; ++i) {
                if(numClusterPoints[i/dimension] != 0) {
                    centroids[i] = tmpCentroids[i] / numClusterPoints[i / dimension];
                    tmpCentroids[i] = 0;
                }
            }
        }
    }

    private int hash(int point, boolean isData, int func) {
        double hashed = 0;
        if(isData) { //differentiate for array
            for(int i = 0; i < dimension; ++i) {
                hashed += data[point + i] * fHash[func * dimension + i]; //calc projection
            }
        } else {
            for(int i = 0; i < dimension; ++i) {
                hashed += centroids[point + i] * fHash[func * dimension + i]; //calc projection
            }
        }
        hashed = hashed / bucketSize; //normalize bucket size

        return (int)Math.floor(hashed); //return bucket
    }

    private void kmeansLsh() {
        double minDist, tmp;
        int clusterID = -1, changed;
        boolean converged = false;
        ArrayList<Integer> sameClusters = new ArrayList<>();


        for(int i = 0; i < tmpCentroids.length; ++i) {
            tmpCentroids[i] = 0;
        }

        //calc hash for all points
        for(int i = 0; i < points; ++i) {
            for(int j = 0; j < funcs; ++j) {
                dataHash[i * funcs + j] = hash(i * dimension, true, j);
            }
        }

        while (!converged) {
            converged = true;
            changed = 0;

            //init tmp variables
            for(int i = 0; i < numClusterPoints.length; ++i)
                numClusterPoints[i] = 0;

            //calc cluster hashes
            for(int i = 0; i < clusters; ++i) {
                for (int j = 0; j < funcs; ++j) {
                    clusterHashes[i*funcs + j] = hash(i * dimension, false, j);
                }
            }



            //for all points
            for (int i = 0; i < points; ++i) {
                boolean same = true;
                sameClusters.clear(); //point is in no cluster
                for(int k = 0; k < clusters; ++k) { //check every cluster point
                    for (int j = 0; j < funcs; ++j) { //check all functions
                        if(dataHash[i * funcs + j] != clusterHashes[k* funcs + j]) { //if in different buckets
                            j = ((j/funcGroups) + 1) * funcGroups; //skip to next function group
                            if(j >= funcs) { //if all func groups tested
                                same = false;
                                break;
                            }
                        }
                        if((j + 1) % funcGroups == 0) { //if full func group is passed
                            same = true;
                            break; //abort calc
                        }
                    }
                    if(same) { //add cluster to point
                        sameClusters.add(k);
                    }
                }

                if(sameClusters.size() != 1) { //if more or less than one cluster in bucket
                    if(sameClusters.size() == 0) { //if no center found, calc minDist over all clusters
                        minDist = distance(i * dimension, 0);
                        clusterID = 0;
                        for (int j = 1; j < clusters; ++j) {
                            tmp = distance(i * dimension, j * dimension);
                            if (tmp < minDist) {
                                minDist = tmp;
                                clusterID = j;
                            }
                        }
                    } else { //find nearest cluster centroid in bucket
                        minDist = -Double.MAX_VALUE;
                        clusterID = 0;
                        for (int j: sameClusters) {
                            tmp = distance(i * dimension, j * dimension);
                            if (tmp < minDist) {
                                minDist = tmp;
                                clusterID = j;
                            }
                        }
                    }
                } else { //only one cluster in bucket
                    clusterID = sameClusters.get(0);
                }


                    if (clusterIDs[i] != clusterID) {
                        converged = false;
                        clusterIDs[i] = clusterID;
                        ++changed;
                    }


                //if(clusterID != -1) {
                    for (int k = 0; k < dimension; ++k) {
                        tmpCentroids[clusterID * dimension + k] += data[i * dimension + k];
                    }
                    ++numClusterPoints[clusterID];
                //}
            }


            //update centroids
            for(int i = 0; i < tmpCentroids.length; ++i) {
                if(numClusterPoints[i/dimension] != 0) {
                    centroids[i] = tmpCentroids[i] / numClusterPoints[i / dimension];
                    tmpCentroids[i] = 0;
                }
            }
            //System.out.println("changed: " + changed);
        }
    }

    public static void main(String[] args) {
        Main prg = new Main();
        prg.run(args);
    }

    //print for debug
    private void printClusters() {
        for (int i = 0; i< clusters; ++i) {
            System.out.println("cluster " + i + ": ");
            for(int j = 0; j < dimension; ++j) {
                System.out.println("   " + centroids[i* dimension + j]);
            }

        }

    }

    public void run(String[] args) {
        loadData(path);
        long start, end;

        try {
            double bSize = Double.parseDouble(args[0]);
            bucketSize = bSize;
            int fGroup = Integer.parseInt(args[1]);
            funcGroups = fGroup;
            funcs = funcGroups * funcGroups;
        } catch (Exception ignore) {}

        System.out.println("calculation without lsh started (n = " + points + ")");
        start = System.nanoTime();
        kmeans();
        end = System.nanoTime();

        long duration = end - start;
        System.out.println("calculation without lsh finished, duration = " + (duration / 1.0e9) +  "s");
        saveData(path + ".out");

        loadData(path); //reload data, new random centroids
        System.out.println("calculation with lsh started (n = " + points + ", bucketsize = " + bucketSize + ")");
        start = System.nanoTime();
        kmeansLsh();
        end = System.nanoTime();

        duration = end - start;
        System.out.println("calculation with lsh finished, duration = " + (duration / 1.0e9) +  "s");
        saveData(path + ".lsh.out");



    }
}
