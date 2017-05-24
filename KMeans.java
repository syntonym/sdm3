
public class KMeans {

    private Double hashFuncs[][] = 
        {
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}
            {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}
            {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}
            {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}
        };

	public static void main(String[] args) {

		HashMap<String, Object> config = readArgs(args);

		Double[][] data = readFile(config.get("path"));

	}

	private HashMap<String, Object> readArgs(String[] args) {
		return null;
	}

	private Double[][] readFile(String path) {
		return null;
	}

	/**
	 * Erstes Argument Hashfunktion
	 * zweites Argument Bucket
	 */
	private List<Integer>[][] hash(Double[][] points) {
	}
	
	/**
	 *  calculates the hash value of a point and a hash function as
	 *  a vector-vector-product
	 */
	private Integer calcHash (Double point[], Double func[]) {
	    if (point.length != func.length) {
	        throw IllegalArgumentException("vector dimensions have to match!");
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
	    
	    Double centroids[clusters][dimension];
	    
	    int randomNum;
	    
	    for (int i = 0; i < clusters; ++i) {
	        // take a point at a random index 
	        // position from points and use it as centroid
            randomNum = rand.nextInt(max + 1);
            centroids[i] = points[randomNum];
        }
        
        Double cHashValues[clusters][dimension];
        
        // calculate the hash value for every centroid for every hash function
        for (int i = 0; i < clusters; ++i) {
            for (int j = 0; j < hashFuncs.length; ++j) {
                cHashValues[i][j] = calcHash (centroids[i], hashFuncs[j]);
            }
        }
	}

	private Double distance(Double[] a, Double[] b) {
		return null;
	}
}
