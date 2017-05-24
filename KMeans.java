import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.IllegalArgumentException;

public class KMeans {

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

		Double[][] data = m.readFile(config.getOptionValue("testdata"));

	}

	private CommandLine readArgs(String[] args) {
		Options options = new Options();

		options.addOption("testdata", true, "Path to the data to readin");
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		return cmd;
	}

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
				System.out.println(line); // kontrolle was eingelesen

		        lines.add(line);
		    }
		} catch (IOException e) {
			System.err.println("Error :" + e);
		} finally {
			try {
				buff.close();
				myFile.close();
			} catch (IOException e) {
				System.err.println("Error :" + e);
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
	private List<Integer>[][] hash(Double[][] points) {
		return null;
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
	    
	    Integer bucket = 0;
	    while (sum > buckets[func][bucket]) {
	        ++bucket;
	    }
	    
	    return bucket;
	}

	private void algorithmus(Double[][] points, List<Integer>[][] buckets) {
	
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
        
        Integer centroidBuckets[][] = new Integer[clusters][dimension];
        
        // calculate the hash value for every centroid for every hash function
        for (int i = 0; i < clusters; ++i) {
            for (int j = 0; j < hashFuncs.length; ++j) {
                centroidBuckets[i][j] = getBucket (centroids[i], j);
            }
        }
	}

	private Double distance(Double[] a, Double[] b) {

		double distance = 0;

		for (int i=0; i<a.length; ++i) {
			distance += Math.pow(a[i]-b[i],2);
		}

		return Math.sqrt(distance);
	}


}
