import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.IllegalArgumentException;
import java.lang.NullPointerException;

import java.lang.*;
import java.util.*;
import org.apache.commons.cli.*;

public class KMeans {

    private String dataPath = "LSH-nmi-adapted.csv";
    private boolean randomizedHashes = false;
    private Integer tries = 6;
    private int run = 0;
    private Integer p = null;
    private boolean bucketWidthAuto;
    private double bucketWidth = -1.0;
    private double[] bucketWidths = new double[10];
    private int amountHashFuncs = 10;
    private double startInitialisationTime = 0.0;
    private double endInitialisationTime = 0.0;
    private double hashFuncs[][] =
    {
        {1., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 1., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 1., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 1., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 1., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 1., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 1.},
    };

    private int amountBuckets;
    private Integer cnt_dist;
    private Integer cnt_LSH;
    private Integer cnt_centroidNaive;

    // save the bucket borders for each function, ignoring the first border (0)
    private double buckets [][] = new double [amountHashFuncs][amountBuckets];
    
    public Integer pointsClusterMap[];
    public Integer correctPointsClusterMap[];

    /**
     * Main access point.
     *
     * to work with unstatic data this only creates a new KMeans object and calls it "real" main method.
     */
    public static void main(String[] args) {
        KMeans m = new KMeans();
        m.unstatic_main(args);
    }

    /**
     * Main Method.
     *
     * This is the "real" main method to access unstatic members of KMeans.
     * There should always only be a single object of KMeans. Running unstatic_main() multiple
     * times is not supported.
     */
    public void unstatic_main(String[] args) {

        // parse command line parameters and set various parameters
        parseCommandLine(args);

        // calculate random hash functions
        if (randomizedHashes) {
			Random randNumber = new Random();
            for (int i = 0; i<amountHashFuncs; i++) {
                for (int j = 0; j<amountHashFuncs; j++) {
                    hashFuncs[i][j] = randNumber.nextDouble();
                }
            }
        }

        double[][] data;

        // read maindata
        try {
            data = readFile(dataPath);
        } catch (NullPointerException e) {
            System.err.println("Couldn't find the testdata. Please provide the testdata via -testdata <path>.\nAborting...");
            System.exit(1);
            return;
        }

		//calculate bucket widths automatically
		if (bucketWidthAuto) {
			bucketWidth = -1.0;
			double[] sd = standardDeviation(data);
			for (int i = 0; i<10; i++) bucketWidths[i] = sd[i];
		} else {
        	for (int i = 0; i<amountHashFuncs; i++) bucketWidths[i]=bucketWidth;
		}

        pointsClusterMap = new Integer[data.length];
        
        double startTime;
        double hashTime;
        double endTime;
        double timeKMeans;
        double timeHashing;

        // run multiple times for better measurements 
        for (int j=0; j<tries; j++) {

            cnt_dist = 0;
            cnt_LSH = 0;
            cnt_centroidNaive = 0;

            // Record time for initial hashing and algorithm

            startTime = System.currentTimeMillis();
            ArrayList<HashMap<Integer, Set<Integer>>> buckets = hash(data);
            hashTime = System.currentTimeMillis();
            // Run the main algorithm
            pointsClusterMap = algorithm(data, buckets, p);
            endTime = System.currentTimeMillis();

            timeHashing = hashTime - startTime;
            timeKMeans = endTime - startTime;

            // Calculate NMI 
            
            ArrayList<Integer> processedPointsClusterMap = new ArrayList<Integer>();
            ArrayList<Integer> processedCorrectPointsClusterMap = new ArrayList<Integer>();

            for (int i = 0; i < pointsClusterMap.length; ++i) {
                processedPointsClusterMap.add(pointsClusterMap[i]);
            }
                    
            for (int i = 0; i < correctPointsClusterMap.length; ++i) {
                processedCorrectPointsClusterMap.add(correctPointsClusterMap[i]);
            }

            double nmiValue = NMI (processedPointsClusterMap, processedCorrectPointsClusterMap);

            // Output only after the first run to warm up JVM
            // Output in JSON for further analysis in python

            if ( j > 0 || tries == 1) {
                System.out.print("{\"p\": " + p + ",\n");
                System.out.print("\"bucket width\": " + bucketWidth + ",\n");
                System.out.print("\"runs\": " + run + ",\n");
                System.out.print("\"NMI\": " + nmiValue + ",\n");
                System.out.print("\"bucketWidthsAuto\": " + bucketWidthAuto + ",\n");
                System.out.print("\"randomizedHashes\": " + randomizedHashes + ",\n");
                System.out.print("\"time_initialisation\": " + (endInitialisationTime - startInitialisationTime) + ",\n");
                System.out.print("\"time_hashing\": " + timeHashing + ",\n");
                System.out.print("\"LSH_count\": " + cnt_LSH + ",\n");
                System.out.print("\"centroidNaive_count\": " + cnt_centroidNaive + ",\n");
                System.out.print("\"distance_cnt\": " + cnt_dist + ",\n");
                System.out.print("\"time\": " + timeKMeans + "}\n");
            }
        }

    }

    private void parseCommandLine(String[] args) {

        // command line parsing

        CommandLine config = readArgs(args);
        p = null;
        randomizedHashes = false;
		bucketWidthAuto = false;

		//command line inputs
        try {
		    p = (Integer.valueOf(config.getOptionValue("p")));
        } catch (Exception e) {
            System.out.println("Error commandline parsing: p\nAssuming default parameter p=6");
        } finally {
            if (p == null) {
                p = 6;
                bucketWidthAuto = false;
                bucketWidth = 33;
                randomizedHashes  = true;
            }
        }

		try {
		    bucketWidth = (Double.valueOf(config.getOptionValue("width")));
        } catch (Exception e) {
            System.out.println("Error commandline parsing: bucketWidth\nAssuming default paramter bucketWidth=30");
        } finally {
			if (bucketWidth == -1.0) bucketWidth = 30; //default
        }

		try {
		    randomizedHashes = (Boolean.valueOf(config.getOptionValue("r")));
        } catch (Exception e) {
            System.out.println("Error commandline parsing: r\nAssuming default parameter false");
        }

		try {
			bucketWidthAuto = (Boolean.valueOf(config.getOptionValue("width_auto")));
        } catch (Exception e) {
            System.out.println("Error commandline parsing: widths_auto\nAssuming default parameter false");
        }

		try {
			String tmp = config.getOptionValue("testdata");
            if (tmp == null) {
                throw new NullPointerException();
            } else {
                dataPath = tmp;
            }
        } catch (Exception e) {
            System.out.println("Error commandline parsing: testdata\nDefaulting to ./LSH-nmi-corrected.csv");
        }
    }

    /**
     * Command line parsing
     */
    private CommandLine readArgs(String[] args) {
        Options options = new Options();

        options.addOption("testdata", true, "Path to the data to readin");
        options.addOption("help", false, "Shows this help");
        options.addOption("p", true, "How many buckets are needed for shortcut");
        options.addOption("r", true, "true -> random hash matrix, otherwise identity matrix");
        options.addOption("width", true, "bucket width");
		options.addOption("width_auto", true, "true -> bucket widths calculated automatically");
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException ex) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("KMeans", options);
            System.exit(1);
        }

        if (cmd.hasOption("help")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("KMeans", options);
            System.exit(0);
        }


        return cmd;
    }

    /**
     * readFile
     *
     * this method reads in the csv-file, parse it
     * and returns the corresponding double[][] values
     * @param path String - path to file
     **/
    private double[][] readFile(String path) {

        // Einlesen des Files und spliten
        FileReader myFile = null;
        BufferedReader buff = null;
        final List<String> lines = new ArrayList<String>();

        try {
            myFile = new FileReader(path);
            buff = new BufferedReader(myFile);
            String line;

            while ((line = buff.readLine()) != null) {
                lines.add(line);
            }
        } catch (IOException|NullPointerException e) {
            System.err.println("Error when reading the file. Is the path to testdata valid?" + 
                    "\nCurrent path: " + dataPath + "\nAborting...");
            System.exit(1);
        } finally {
            try {
                buff.close();
                myFile.close();
            } catch (IOException|NullPointerException e) {
                System.err.println("Error when closing the files.\nAborting...");
                System.exit(1);
            }
        }

        final String[][] valuesArray = new String[lines.size()][];
        int cnt = 0;
        for (final String line : lines) {
            valuesArray[cnt++] = line.split(",");
        }

        correctPointsClusterMap = new Integer[lines.size()];
        // the -1 is necessary since the last column in the source csv file is the corresponding centroid
        double[][] valuesDouble = new double[lines.size()][valuesArray[0].length-1]; 

        for (int i=0; i<lines.size(); ++i) {
            for (int j=0; j<valuesArray[0].length; ++j) {
                if (j < 10) 
                    valuesDouble[i][j] = Double.parseDouble(valuesArray[i][j]);
                else 
                    correctPointsClusterMap[i] = Integer.parseInt(valuesArray[i][j]);
            }
        }

        return valuesDouble;
    }

    /**
     * Hash all points once for fast bucket access.
     *
     * returns  represent the hashed data points in a different data format.
     *          For each hashfunction buckets has a seperate hashmap.
     *          These hashmap map the bucketID to the set of point indexes which are in this bucket.
     *          So hash(...).get(hash_index).get(bucket_index) gets the points which belong to bucket bucket_index
     *          for the hashfunction hash_index.
     *
     */
    private ArrayList<HashMap<Integer, Set<Integer>>> hash(double[][] points) {
        ArrayList<HashMap<Integer, Set<Integer>>> buckets = new ArrayList<HashMap<Integer, Set<Integer>>>();

        for (int funci=0; funci<hashFuncs.length; funci++) {
            buckets.add(funci, new HashMap<>());
        }

        double[] min = new double[points[0].length];
        for (int i=0; i<min.length; i++) {
            min[i] = Double.POSITIVE_INFINITY;
        }
        double[] max = new double[points[0].length];
        for (int i=0; i<max.length; i++) {
            max[i] = Double.NEGATIVE_INFINITY;
        }

        for (int i=0; i<points.length; i++) {
            for (int j=0; j<points[i].length; j++) {
                if (min[j] > points[i][j]) {
                    min[j] = points[i][j];
                }
                if (max[j] < points[i][j]) {
                    max[j] = points[i][j];
                }
            }

            for (int funci=0; funci<hashFuncs.length; funci++) {
                Integer bucketIndex = getBucket(points[i], funci);
                Set<Integer> set = buckets.get(funci).get(bucketIndex);
                if (set == null) {
                    set = new HashSet<>();
                    buckets.get(funci).put(bucketIndex, set);
                }
                set.add(i);
            }

        }

        return buckets;
    }

	/**
     * getBucket
     *
     * this method calculates the hash value of a point and a hash function as
     * a vector-vector-product and returns the resulting bucket number
     * @param point   double[] - data point
     * @param int    double[] - index of hash function
     **/
    private Integer getBucket (double point[], int func) {
        if (point.length != hashFuncs[func].length) {
            throw new IllegalArgumentException("vector dimensions have to match!");
        }

        double sum = 0.0;
        for (int i = 0; i < point.length; ++i) {
            sum += (point[i] * hashFuncs[func][i]);
        }

        Integer r = (int) (sum / bucketWidths[func]);
        if (r == null) {
            throw new NullPointerException();
        }
        return r;
    }
        
    
    
    

    /**
     * Run the algorithm a single time.
     * 
     * points   are the data points as returned by readFile(). 
     * buckets  represent the hashed data points in a different data format.
     *          For each hashfunction buckets has a seperate hashmap.
     *          These hashmap map the bucketID to the set of point indexes which are in this bucket.
     *          So buckets.get(hash_index).get(bucket_index) gets the points which belong to bucket bucket_index
     *          for the hashfunction hash_index.
     * p        parameter how many buckets need to agree to form a field as discussed in the PDF.
     *
     * returns an array that maps pointindex to cluster index. So algorithm(...)[10] is the clusterID of the tenth point.
     */
    private Integer[] algorithm (double[][] points, ArrayList<HashMap<Integer, Set<Integer>>> buckets, Integer p) {

        /*
         * centroid initialisation and hashing
         */

        int clusters = 15;
        int dimension = 10;
        int max = points.length;

        double centroids[][] = new double[clusters][dimension];

        // Initialisation of Centroids

        startInitialisationTime = System.currentTimeMillis();
        Random rand = new Random();
        int randomNum = rand.nextInt(max +1);
        centroids[0] = points[randomNum].clone();

        for (int cluster_index=1; cluster_index<clusters; cluster_index++) {

            double total_min = Double.POSITIVE_INFINITY;
            double total_max = Double.NEGATIVE_INFINITY;

            double distance_total = 0.0;
            double distance [] = new double[points.length];

            for (int i=0; i<points.length; i++) {
                double min = Double.POSITIVE_INFINITY;
                for (int j=0; j<cluster_index; j++) {
                    double d = distance2(centroids[j], points[i]);
                    if (d < min) {
                        min = d;
                    }
                }

                distance[i] = min;
                distance_total += min;

                if (min < total_min) {
                    total_min = min;
                }
                if (min > total_max) {
                    total_max = min;
                }
            }

            for (int i=0; i<points.length; i++) {
                double randomNumber = distance_total * rand.nextDouble();
                if (randomNumber < distance[i]) {
                    centroids[cluster_index] = points[i].clone();
                    break;
                } else {
                    distance_total -= distance[i];
                }
            }
        }
        endInitialisationTime = System.currentTimeMillis();

        Integer centroidBuckets[][] = new Integer[clusters][amountHashFuncs];

        // calculate the hash value for every centroid for every hash function
        for (int i = 0; i < clusters; ++i) {
            for (int j = 0; j < hashFuncs.length; ++j) {
                centroidBuckets[i][j] = getBucket (centroids[i], j);
            }
        }


        /*
         * let's get serious now...
         * Main Algorithm
         */

        // on each index we save the centroid index for the corresponding point
        Integer pointsClusterMap [] = new Integer[points.length];

        // dirty is true if any point changed clusters
        boolean dirty = true;
        run = 0;

        while (dirty) {
            run++;
            dirty = false;

            // isOnlyCentroid is true on index i if centroid i is the only centroid in its field.
            Boolean isOnlyCentroid [] = new Boolean[clusters];
            Integer fieldID [] = new Integer[clusters];

            for (int i = 0; i < clusters; ++i) {
                isOnlyCentroid[i] = true;
                fieldID[i] = i;
            }

            // remember already processed points 
            Set<Integer> processed = new HashSet<>(points.length/100);

            // check which centroids are alone in its field
            for (int i = 0; i < clusters; ++i) {
                for (int j = i + 1; j < clusters; ++j) {
                    int sim = 0;
                    for (int k = 0; k < amountHashFuncs; ++k) {
                        if (centroidBuckets[i][k] == centroidBuckets[j][k]) {
                            sim += 1;

                            if (sim >= p) {
                                fieldID[j] = fieldID[i];
                                isOnlyCentroid[i] = false;
                                isOnlyCentroid[j] = false;
								break;
                            }
                        }
                    }
                }
            }

            // assign all corresponding points to all lonely centroids
            for (int i = 0; i <  clusters; ++i) {
                if (isOnlyCentroid[i]) {
                    // we get amountHashFuncs (default 10) sets of Integers, which
                    // represent the index of all points in the corresponding bucket
                    ArrayList<Set<Integer>> allBucketPoints = new ArrayList<Set<Integer>>(amountHashFuncs);
                    for (int j = 0; j < amountHashFuncs; ++j) {
                        allBucketPoints.add(j, buckets.get(j).get(centroidBuckets[i][j]));
                    }

                    // we calculate the intersection of all sets to obtain
                    // all points which are in the exact same field as the centroid
                    Set<Integer> bucketPoints = allBucketPoints.get(0);
                    if (bucketPoints == null) {
                        bucketPoints = new HashSet<>();
                    }
                    for (int j = 1; j < amountHashFuncs; ++j) {
                        Set<Integer> tmpPoints = allBucketPoints.get(j);
                        if (tmpPoints == null) {
                            bucketPoints = new HashSet<>();
                            break;
                        }
                        bucketPoints.retainAll (tmpPoints);
                    }

                    // finally we assign the clusterIds to the points
                    for (Integer bi : bucketPoints) {
                        if (pointsClusterMap[bi] == null || pointsClusterMap[bi] != i) {
                            pointsClusterMap[bi] = i;
                            dirty = true;
                            cnt_LSH += clusters;
                        }
                    }
                    processed.addAll(bucketPoints);
                }
            }


            // assign points naively to centroids in field with more than one cluster
            for (int i = 0; i < clusters; ++i) {
                if (!isOnlyCentroid[i]) {
                    
                    // we get amountHashFuncs (default 10) sets of Integers, which
                    // represent the index of all points in the corresponding bucket
                    ArrayList<Set<Integer>> allBucketPoints = new ArrayList<Set<Integer>>(amountHashFuncs);
                    for (int j = 0; j < amountHashFuncs; ++j) {
                        allBucketPoints.add(j, buckets.get(j).get(centroidBuckets[i][j]));
                    }

                    // we calculate the intersection of all sets to obtain
                    // all points which are in the exact same field as the centroids are
                    Set<Integer> bucketPoints = allBucketPoints.get(0);
                    for (int j = 1; j < amountHashFuncs; ++j) {
                        bucketPoints.retainAll (allBucketPoints.get(j));
                    }

                    // calculate the minimum to clusters which are in the same bucket
                    for (Integer point_index : bucketPoints) {
                        double min_distance = Double.POSITIVE_INFINITY;
                        int min_centroid_index = -1;

                        for (int centroid_index = 0; centroid_index < clusters; centroid_index++) {
                            if (fieldID[centroid_index] != fieldID[i]) {
                                // count how many distance calculations we save
                                cnt_centroidNaive++;
                                continue;
                            }

                            double d = distance(points[point_index], centroids[centroid_index]);
				            cnt_dist++;
                            if (d < min_distance) {
                                min_distance = d;
                                min_centroid_index = centroid_index;
                            }
                        }

                        if (pointsClusterMap[point_index] == null || pointsClusterMap[point_index] != min_centroid_index) {
                            pointsClusterMap[point_index] = min_centroid_index;
                            dirty = true;
                        }
                        processed.add(point_index);
                    }

                }
            }

            // Calculate centroid for any point that is not already processed
            
            for (int point_index=0; point_index < points.length; point_index++ ) {
                // skip already processed points
                if (processed.contains(point_index)) {
                    continue;
                }
                double min_distance = Double.POSITIVE_INFINITY;
                int min_centroid_index = -1;

                for (int centroid_index=0; centroid_index < clusters; centroid_index++) {
                    double d = distance(points[point_index], centroids[centroid_index]);
                    cnt_dist++;
                    if (d < min_distance) {
                        min_distance = d;
                        min_centroid_index = centroid_index;
                    }
                }

                if (pointsClusterMap[point_index] == null || min_centroid_index != pointsClusterMap[point_index] ) {
                    pointsClusterMap[point_index] = min_centroid_index;
                    dirty = true;
                }
            }

            // reset centroids

            for (int i=0; i<centroids.length; i++) {
                for (int j=0; j<dimension; j++) {
                    centroids[i][j] = 0.0;
                }
            }

            // --- Recalculate Centroids ---
            
            // count how many points are assigned to a centroid
            Integer weight [] = new Integer[points.length];

            for (int i=0; i<weight.length; i++) {
                weight[i] = 0;
            }

            // sum all points assigned to the centroid
            for (int point_index=0; point_index < points.length; point_index++) {
                Integer centroid_index = pointsClusterMap[point_index];
                weight[centroid_index]++;

                for (int i=0; i<dimension; i++) {
                    centroids[centroid_index][i] += points[point_index][i];
                }
            }

            // Divide by weight to calculate the average

            for (int i=0; i<centroids.length; i++) {
                for (int j=0; j<dimension; j++) {
                    if (weight[i] != 0) {
                        centroids[i][j] /= weight[i];
                    }
                }
            }

            // recalculate the hash value for every centroid for every hash function
            for (int i = 0; i < clusters; ++i) {
                for (int j = 0; j < hashFuncs.length; ++j) {
                    centroidBuckets[i][j] = getBucket (centroids[i], j);
                }
            }
        }
 
        return pointsClusterMap;
    }

    /**
     * distance between to datapoints
     *
     * this method returns the euclidian distance between two datapoints
     * @param a    double[]
     * @param b    double[]
     **/
    private double distance(double[] a, double[] b) {

        double distance = 0;

        for (int i=0; i<a.length; ++i) {
            distance += Math.pow(a[i]-b[i],2.0);
        }

        return Math.sqrt(distance);
    }

	 /**
     * squared distance between to datapoints
     *
     * this method returns the squared distance between two datapoints
     * @param a    double[]
     * @param b    double[]
     **/
     private double distance2(double[] a, double[] b) {

        double distance = 0;

        for (int i=0; i<a.length; ++i) {
            distance += Math.pow(a[i]-b[i],2.0);
        }

        return distance;
    }


	/**
     * standard deviations
     *
     * this method calculates the standard deviation of a double[][],
	 * where the first field stores all datapoints and the second field
	 * their corresponding dimensions. standard deviations are calculated for each dimension
	 * and stored in a double[]
	 *
     * @param a    double[][]
     **/
	private double[] standardDeviation(double[][] a) {

		double mean;
		double variance;
		double[] sd = new double[a[0].length];

		for (int dim = 0; dim < a[0].length; dim++) {

			mean = 0;
			variance = 0;

			for (int i = 0; i<a.length; i++) {
				mean += a[i][dim];
			}

			mean /= a.length;

			for (int j = 0; j<a.length; j++) {
				variance += Math.pow(a[j][dim]-mean,2.0);
			}

			variance /= a.length;

			sd[dim] = Math.sqrt(variance);

		}

		return sd;

	}

    /* copied and pasted, original function by martin perdacher */
    public static double NMI(ArrayList<Integer> one, ArrayList<Integer> two){
		if(one.size()!=two.size()){
			throw new IllegalArgumentException("Sizes don't match!");
		}
		int maxone = Collections.max(one);
		int maxtwo = Collections.max(two);

		double[][] count = new double[maxone+1][maxtwo+1];
		for(int i=0;i<one.size();i++){
			count[one.get(i)][two.get(i)]++;
		}
		//i<maxone=R
		//j<maxtwo=C
		double[] bj = new double[maxtwo+1];
		double[] ai = new double[maxone+1];

		for(int m=0;m<(maxtwo+1);m++){
			for(int l=0;l<(maxone+1);l++){
				bj[m]=bj[m]+count[l][m];
			}
		}
		for(int m=0;m<(maxone+1);m++){
			for(int l=0;l<(maxtwo+1);l++){
				ai[m]=ai[m]+count[m][l];
			}
		}

		double N=0;
		for(int i=0;i<ai.length;i++){
			N=N+ai[i];
		}
		double HU = 0;
		for(int l=0;l<ai.length;l++){
			double c=0;
			c=(ai[l]/N);
			if(c>0){
				HU=HU-c*Math.log(c);
			}
		}

		double HV = 0;
		for(int l=0;l<bj.length;l++){
			double c=0;
			c=(bj[l]/N);
			if(c>0){
				HV=HV-c*Math.log(c);
			}
		}
		double HUstrichV=0;
		for(int i=0;i<(maxone+1);i++){
			for(int j=0;j<(maxtwo+1);j++){
				if(count[i][j]>0){
					HUstrichV=HUstrichV-count[i][j]/N*Math.log(((count[i][j])/(bj[j])));
				}
			}
		}
		double IUV = HU-HUstrichV;
		double reto = IUV/(Math.max(HU, HV));

		return reto;
	}
}
