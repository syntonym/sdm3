import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.IllegalArgumentException;
import java.lang.NullPointerException;

public class KMeans {

    private Double[] bucketWidths = {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.};
    private int amountHashFuncs = 10;
    private Double hashFuncs[][] =
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

    // save the bucket borders for each function, ignoring the first border (0)
    private Double buckets [][] = new Double [amountHashFuncs][amountBuckets];

    public static void main(String[] args) {

        KMeans m = new KMeans();
        CommandLine config = m.readArgs(args);

        try {
            Double[][] data = m.readFile(config.getOptionValue("testdata"));
        } catch (NullPointerException e) {
            System.err.println("Konnte das File nicht finden\n" + e);
            System.exit(1);
        }

        double startTime;
    	double endTime;
		double timeKMeans;

		startTime = System.currentTimeMillis();
		    //place your function here
		endTime = System.currentTimeMillis();

		timeKMeans = endTime - startTime;

		System.out.print("time: " + timeKMeans + "\n");

    }

    private CommandLine readArgs(String[] args) {
        Options options = new Options();

        options.addOption("testdata", true, "Path to the data to readin");
        options.addOption("help", false, "Shows this help");
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException ex) {
            System.out.println(ex);
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
     * and returns the corresponding Double[][] values
	 * @param path String - path to file
	 **/
    private Double[][] readFile(String path) {

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
            System.err.println("Konnte das File nicht finden\n" + e);
            System.exit(1);
        } finally {
            try {
                buff.close();
                myFile.close();
            } catch (IOException|NullPointerException e) {
                System.err.println("Error1 :" + e);
                System.exit(1);
            }
        }

        final String[][] valuesArray = new String[lines.size()][];
        int cnt = 0;
        for (final String line : lines) {
            valuesArray[cnt++] = line.split(",");
        }

        Double[][] valuesDouble = new Double[lines.size()][valuesArray[0].length];

        for (int i=0; i<lines.size(); ++i) {
            for (int j=0; j<valuesArray[0].length; ++j) {
                valuesDouble[i][j] = Double.parseDouble(valuesArray[i][j]);
            }
        }

        return valuesDouble;
    }

    /**
     * Erstes Argument Hashfunktion
     * zweites Argument Bucket
     * TODO:
     *  - calculate the max and the min hash values
     *    to obtain bucket borders, which should be saved
     *    in the field buckets
     */
    private ArrayList<HashMap<Integer, Set<Integer>>> hash(Double[][] points) {
        ArrayList<HashMap<Integer, Set<Integer>>> buckets = new ArrayList<HashMap<Integer, Set<Integer>>>();

        for (int funci=0; funci<hashFuncs.length; funci++) {
            buckets.add(funci, new HashMap<>());
        }

        Double[] min = new Double[points[0].length];
        for (int i=0; i<min.length; i++) {
            min[i] = Double.POSITIVE_INFINITY;
        }
        Double[] max = new Double[points[0].length];
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
                set.add(bucketIndex);
            }

        }

        return buckets;
    }

    /**
     *  calculates the hash value of a point and a hash function as
     *  a vector-vector-product
     */
    private Integer getBucket (Double point[], int func) {
        if (point.length != hashFuncs[func].length) {
            throw new IllegalArgumentException("vector dimensions have to match!");
        }

        Double sum = 0.0;
        for (int i = 0; i < point.length; ++i) {
            sum += (point[i] * hashFuncs[func][i]);
        }

        return new Double(sum / (bucketWidths[func])).intValue();
    }


    private void algorithmus(Double[][] points, List<Integer>[][] buckets) {

        /*
         * centroid initialisation and hashing
         */

        int clusters = 15;
        int dimension = 10;
        int max = points.length;

        Double centroids[][] = new Double[clusters][dimension];

        Random rand = new Random();
        int randomNum;

        for (int i = 0; i < clusters; ++i) {
            // take a point at a random index
            // position from points and use it as centroid
            randomNum = rand.nextInt(max + 1);
            centroids[i] = points[randomNum];
        }

        Integer centroidBuckets[][] = new Integer[clusters][amountHashFuncs];

        // calculate the hash value for every centroid for every hash function
        for (int i = 0; i < clusters; ++i) {
            for (int j = 0; j < hashFuncs.length; ++j) {
                centroidBuckets[i][j] = getBucket (centroids[i], j);
            }
        }


        /*
         * let's get serious now...
         */

        // on each index we save the centroid index for the corresponding point
        Integer pointsClusterMap [] = new Integer[points.length];

        // isOnlyCentroid is true on index i if centroid i is the only centroid in its field.
        Boolean isOnlyCentroid [] = new Boolean[clusters];

        for (int i = 0; i < clusters; ++i) {
            isOnlyCentroid[i] = true;
        }

        // check which centroids are alone in its field
        for (int i = 0; i < clusters; ++i) {
            for (int j = i + 1; j < clusters; ++j) {
                int sim = 0;
                for (int k = 0; k < amountHashFuncs; ++k) {
                    if (centroidBuckets[i][k] == centroidBuckets[j][k]) {
                        sim += 1;
                        if (sim == amountHashFuncs) {
                            isOnlyCentroid[i] = false;
                            isOnlyCentroid[j] = false;
                        }
                    }
                }
            }
        }

        // assign all corresponding points to all lonely centroids
        for (int i = 0; i <  clusters; ++i) {
            if (isOnlyCentroid[i]) {
                for (int j = 0; j < amountHashFuncs; ++j) {
                    List<Integer> pointsBuf = buckets[j][centroidBuckets[i][j]];
                    for (int k = 0; k < pointsBuf.size(); ++k) {
                        pointsClusterMap[pointsBuf.get(k)] = i;
                    }
                }
            }
        }
    }

    /**
	 * distance between to datapoints
	 *
	 * this method returns the euclidian distance between two datapoints
	 * @param a    Double[]
	 * @param b    Double[]
	 **/
    private Double distance(Double[] a, Double[] b) {

        double distance = 0;

        for (int i=0; i<a.length; ++i) {
            distance += Math.pow(a[i]-b[i],2);
        }

        return Math.sqrt(distance);
    }


}
