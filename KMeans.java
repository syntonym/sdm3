import java.util.*;
import java.lang.*;
import org.apache.commons.cli.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.IllegalArgumentException;

public class KMeans {

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
 
		final Double[][] valuesArray = new Double[lines.size()][];
		int cnt = 0;
		for (final String line : lines) valuesArray[cnt++] = Double.parseDouble(line.split(","));
 
		return valuesArray;
	}

	/**
	 * Erstes Argument Hashfunktion
	 * zweites Argument Bucket
	 */
	private List<Integer>[][] hash(Double[][] points) {
		return null;
	}
	
	/**
	 *  calculates the hash value of a point and a hash function as
	 *  a vector-vector-product
	 */
	private Integer calcHash (Double point[], Double func[]) {
	    if (point.length != func.length) {
	        throw new IllegalArgumentException("vector dimensions have to match!");
        }
        
	    Integer sum;
	    for (int i = 0; i < point.length; ++i) {
	        sum += (point[i] * func[i]);
	    }
	    
	    return sum;
	}

	private void algorithmus(Double[][] points, List<Integer>[][] buckets) {
	
	    int clusters = 15;
	    int dimension = 10;
	    int max = points.length;
	    
	    Double centroids[][] = new Double[clusters][dimension];
	    
	    int randomNum;
	    
	    for (int i = 0; i < clusters; ++i) {
	        // take a point at a random index 
	        // position from points and use it as centroid
            randomNum = rand.nextInt(max + 1);
            centroids[i] = points[randomNum];
        }
        
        Double cHashValues[][] = new Double[clusters][dimension];
        
        // calculate the hash value for every centroid for every hash function
        for (int i = 0; i < clusters; ++i) {
            for (int j = 0; j < hashFuncs.length; ++j) {
                cHashValues[i][j] = calcHash (centroids[i], hashFuncs[j]);
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
