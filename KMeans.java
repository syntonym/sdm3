
public class KMeans {

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

	private void algorithmus(Double[][] points, List<Integer>[][] buckets) {
	}

	private Double distance(Double[] a, Double[] b) {

		double distance = 0;

		for (int i=0; i<a.length; ++i) {
			distance += Math.pow(a[i]-b[i],2);
		}

		return Math.sqrt(distance);
	}


}
