import java.lang.System;
import java.util.Random;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {

	public enum Strategy {
		LLOYD, MACQUEEN
	}

	public enum Initialisation {
		RANDOM_PARTITION,
		RANDOM_CLUSTER_CENTERS
	}

	public static void main(String[] args) {
		int n;
		int dimensions;
		int k;
		Strategy strat;
		Initialisation init;
		ArrayList<ArrayList<Integer>> points_it = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Point>> centroid_it = new ArrayList<ArrayList<Point>>();
		
	    
	    if (args.length == 5) {
		    n = Integer.parseInt(args[0]);
		    dimensions = Integer.parseInt(args[1]);
		    k = (int) Integer.parseInt(args[2]);
		    strat= Strategy.valueOf(args[3]);
		    init = Initialisation.valueOf(args[4]);
	    } else {
	        n = 4000;
	        dimensions = 2;
	        k = 4;
	        strat = Strategy.LLOYD;
	        init = Initialisation.RANDOM_CLUSTER_CENTERS;
	    }
		
		Point[] points = generateData(k, n, dimensions);
		visualize(points, init, strat, n, k, points_it, centroid_it);

		Scanner reader = new Scanner(System.in);  // Reading from System.in
		System.out.println("Enter a number (k <- number of clusters): ");
		k = reader.nextInt(); // Scans the next token of the input as an int.
		Visualizer.updateK(new Integer(k).toString());

		algoKMeans(points, init, strat, k, n, points_it, centroid_it);
	}

	/**
	 * return generated point data (exercise part 1)
	 *
	 * @param clusters      number of clusters generated (K on the assignment sheet)
	 * @param points        total number of points (N on the assignment sheet)
	 * @param dimensions    number of dimensions in which a single points lives (D on the assignment sheet)
	 * @param deviation_max maximum of deviation. Should be chosen in such a way that clusters overlap.
	 * @param mean_max      maxiumum of mean. Should be chosen in such a way that clusters overlap.
	 **/
	public static Point[] generateData(int clusters, int points, int dimensions, double deviation_max, double mean_max) {
		Point[] r_points = new Point[points];
		Random r = new Random();
		int cluster_id;
		Point p;

		if (deviation_max == 0.0) {
			deviation_max = 5;
		}
		if (mean_max == 0.0) {
			mean_max = 100;
		}

		double[][][] cluster = new double[clusters][2][dimensions];
		double deviation;
		double mean;

		for (int i=0; i<cluster.length; i++) {

			deviation = r.nextDouble() * deviation_max;
			cluster[i][0] = new double[1];
			cluster[i][0][0] = deviation;
			cluster[i][1] = new double[dimensions];

			for (int k=0; k<dimensions; k++) {
				mean = r.nextDouble() * mean_max;
				cluster[i][1][k] = mean;
			}
		}

		for (int i=0; i<r_points.length; i++) {
			cluster_id = r.nextInt(clusters);
			p = new Point(dimensions);
			p.original_category = cluster_id;
			for (int d=0; d<dimensions; d++) {
				p.values[d] = r.nextGaussian() * cluster[cluster_id][0][0] + cluster[cluster_id][1][d];
			}
			r_points[i] = p;
		}

		return r_points;
	}

	public static Point[] generateData(int clusters, int points, int dimensions) {
		return generateData(clusters, points, dimensions, 0.0, 0.0);
	}

	/**
	 * Applies the k-means algorithm on the input data points.
	 * 
	 * The cluster should be indicated by the point.category attribute.
	 * A -1 represents "not-clustered". Points with the same integer in
	 * the point.category attribute should belong to the same cluster.
	 *
	 * Internally you can ofcourse use whatever datastructure you like.
	 *
	 * Return nothing (mutate the array).
	 * 
	 * @param points Array of points to cluster
	 *
	 **/
	public static void algoKMeans(Point[] points, Initialisation initialisation, Strategy strategy, int k, int n, ArrayList<ArrayList<Integer>> points_it, ArrayList<ArrayList<Point>> centroid_it) {

		int count = 0;

		if (strategy == Strategy.LLOYD && initialisation == Initialisation.RANDOM_CLUSTER_CENTERS) {

			Point[] centroids = new Point[k];
			for (int z = 0; z<k; z++) centroids[z] = new Point(points[0].dim);
			int randomIndex;
			int i = 0;

			//find k random start centroids
			while (i<k) {
				randomIndex = (int) ((Math.random()*n));

				if (points[randomIndex].category == -1) {

					for (int l=0; l<points[randomIndex].dim; l++) {
						centroids[i].values[l] = points[randomIndex].values[l];
					}
					centroids[i].category = i;
					points[randomIndex].category = i++;
				}

			}

			Boolean change = true;

			//iterate until there is no change in category
			while (change) {

				//centroid_it.add(centroids);
				centroid_it.add(new ArrayList<Point>());
				for (int j = 0; j < centroids.length; j++) {
					centroid_it.get(count).add(new Point(centroids[j]));
				}

				change = false;

				//Zuordnung der Datenpunkte
				for (int j = 0; j<n; j++) {

					double distance = Double.MAX_VALUE;
					double distance_old = distance;
					double[] vector = new double[points[0].dim];
					points[j].category_change = points[j].category;

					for (int l = 0; l < k; l++) {

						for (int b=0; b<points[0].dim; b++) {	
							vector[b] = centroids[l].values[b] - points[j].values[b];
						}

						//distance between point and centroid
						distance = 0;
						for (int c=0; c<points[0].dim; c++) {
							distance += Math.sqrt(Math.pow(vector[c],2));
						}
			
						if (distance < distance_old) {
							points[j].category = centroids[l].category;
							distance_old = distance;
						}
					}

					if (points[j].category != points[j].category_change) {
							change = true;
					}
				}

				points_it.add(new ArrayList<Integer>());
				for (int j = 0; j < points.length; j++) {
					points_it.get(count).add(new Integer(points[j].category));
				}

				//Neu Berechnung der Centroide
				for (int z=0; z<points[0].dim; z++) {
					for (int y = 0; y<k; y++) {
						double mean = 0;
						double elementCount = 0;
	
						for (int j = 0; j<n; j++) {
							if (points[j].category == y) {
								mean += points[j].values[z];
								elementCount++;
							}
						}
						centroids[y].values[z] = mean / elementCount;
	
					}	
				}
				count++;
				System.out.println("Iterate " + count + "\n");	
	
			}

		} else {
			if (strategy == Strategy.LLOYD && initialisation == Initialisation.RANDOM_PARTITION) {
				Point[] centroids = new Point[k];
				for (int z = 0; z<k; z++) {
					centroids[z] = new Point(points[0].dim);
					centroids[z].category = z;
				}

				int randomIndex;
				int i = 0;

				//random partition to clusters
				for (int a=0; a<n; a++) {
					points[a].category = (int) ((Math.random()*k));
				}

				//assure that every cluster has at least one member
				for (int a=0; a<n; a++) {
					randomIndex = (int) ((Math.random()*n));
					points[randomIndex].category = a;
				}
				
				Boolean change = true;

				//iterate until there is no change in category
				while (change) {

					change = false;

					//Neu Berechnung der Centroide
					for (int z=0; z<points[0].dim; z++) {
						for (int y = 0; y<k; y++) {
							double mean = 0;
							double elementCount = 0;
	
							for (int j = 0; j<n; j++) {
								if (points[j].category == y) {
									mean += points[j].values[z];
									elementCount++;
								}
							}	
							centroids[y].values[z] = mean / elementCount;
	
						}
	
					}

					centroid_it.add(new ArrayList<Point>());
					for (int j = 0; j < centroids.length; j++) {
						centroid_it.get(count).add(new Point(centroids[j]));
					}

					points_it.add(new ArrayList<Integer>());
					for (int j = 0; j < points.length; j++) {
						points_it.get(count).add(new Integer(points[j].category));
					}

					//Zuordnung der Datenpunkte
					for (int j = 0; j<n; j++) {

						double distance = Double.MAX_VALUE;
						double distance_old = distance;
						double[] vector = new double[points[0].dim];
						points[j].category_change = points[j].category;

						for (int l = 0; l < k; l++) {

							for (int b=0; b<points[0].dim; b++) {
								vector[b] = centroids[l].values[b] - points[j].values[b];
							}

							//distance between point and centroid
							distance = 0;
							for (int c=0; c<points[0].dim; c++) {
								distance += Math.sqrt(Math.pow(vector[c],2));
							}
			
							if (distance < distance_old) {
								points[j].category = centroids[l].category;
								distance_old = distance;
							}
						}
	
						if (points[j].category != points[j].category_change) {
								change = true;
						}
					}

					count++;
					System.out.println("Iterate " + count + "\n");
				}
			} else {
				if (strategy == Strategy.MACQUEEN && initialisation == Initialisation.RANDOM_CLUSTER_CENTERS) {
					Point[] centroids = new Point[k];
					for (int z = 0; z<k; z++) centroids[z] = new Point(points[0].dim);
					int randomIndex;
					int i = 0;

					//find k random start centroids
					while (i<k) {
						randomIndex = (int) ((Math.random()*n));

						if (points[randomIndex].category == -1) {
							for (int l=0; l<points[0].dim; l++) {
								centroids[i].values[l] = points[randomIndex].values[l];
							}
							centroids[i].category = i;
							points[randomIndex].category = i++;
						}

					}

					Boolean change = true;

					//iterate until there is no change in category
					while (change) {

						change = false;

						centroid_it.add(new ArrayList<Point>());
						for (int j = 0; j < centroids.length; j++) {
							centroid_it.get(count).add(new Point(centroids[j]));
						}

						//Zuordnung der Datenpunkte
						for (int j = 0; j<n; j++) {

							double distance = Double.MAX_VALUE;
							double distance_old = distance;;
							double[] vector = new double[points[0].dim];
							points[j].category_change = points[j].category;

							for (int l = 0; l < k; l++) {

								for (int b=0; b<points[0].dim; b++) {
									vector[b] = centroids[l].values[b] - points[j].values[b];
								}

								//distance between point and centroid
								distance = 0;
								for (int c=0; c<points[0].dim; c++) {
									distance += Math.sqrt(Math.pow(vector[c],2));
								}
			
								if (distance < distance_old) {
									points[j].category = centroids[l].category;
									distance_old = distance;
								}
							}

							if (points[j].category != points[j].category_change) {
								change = true;
				
								//if there's a change -> calculate centroids
								for (int z=0; z<points[0].dim; z++) { 
									for (int y = 0; y<k; y++) {
										double mean = 0;
										double elementCount = 0;
	
										for (int x = 0; x<n; x++) {
											if (points[x].category == y) {
												mean += points[x].values[z];
												elementCount++;
											}
										}
										centroids[y].values[z] = mean / elementCount;
	
									}	
								}

							}
						}

						points_it.add(new ArrayList<Integer>());
						for (int j = 0; j < points.length; j++) {
							points_it.get(count).add(new Integer(points[j].category));
						}

						count++;
						System.out.println("Iterate " + count + "\n");	
					}
				} else {
					if (strategy == Strategy.MACQUEEN && initialisation == Initialisation.RANDOM_PARTITION) {
						
						Point[] centroids = new Point[k];
						for (int z = 0; z<k; z++) {
							centroids[z] = new Point(points[0].dim);
							centroids[z].category = z;
						}

						int randomIndex;
						int i = 0;

						//random partition to clusters
						for (int a=0; a<n; a++) {
							points[a].category = (int) ((Math.random()*k));
						}

						//calculate centroids
						for (int z=0; z<points[0].dim; z++) {
							for (int y = 0; y<k; y++) {
								double mean = 0;
								double elementCount = 0;
	
								for (int j = 0; j<n; j++) {
									if (points[j].category == y) {
										mean += points[j].values[z];
										elementCount++;
									}
								}	
								centroids[y].values[z] = mean / elementCount;
	
							}	
						}
				
						Boolean change = true;
	
						//iterate until there is no change in category
						while (change) {

							centroid_it.add(new ArrayList<Point>());
							for (int j = 0; j < centroids.length; j++) {
								centroid_it.get(count).add(new Point(centroids[j]));
							}

							points_it.add(new ArrayList<Integer>());
							for (int j = 0; j < points.length; j++) {
								points_it.get(count).add(new Integer(points[j].category));
							}

							change = false;

							
							//Zuordnung der Datenpunkte
							for (int j = 0; j<n; j++) {

								double distance = Double.MAX_VALUE;
								double distance_old = distance;
								double[] vector = new double[points[0].dim];
								points[j].category_change = points[j].category;

								for (int l = 0; l < k; l++) {

									for (int b=0; b<points[0].dim; b++) {
										vector[b] = centroids[l].values[b] - points[j].values[b];
									}

									//distance between point and centroid
									distance = 0;
									for (int c=0; c<points[0].dim; c++) {
										distance += Math.sqrt(Math.pow(vector[c],2));
									}
			
									if (distance < distance_old) {
										points[j].category = centroids[l].category;
										distance_old = distance;
									}
								}
	
								if (points[j].category != points[j].category_change) {
										change = true;

										//if there's a change -> calculate centroids
										for (int z=0; z<points[0].dim; z++) {
											for (int y = 0; y<k; y++) {
												double mean = 0;
												double elementCount = 0;
	
												for (int x = 0; x<n; x++) {
													if (points[x].category == y) {
														mean += points[x].values[z];
														elementCount++;
													}
												}
												centroids[y].values[z] = mean / elementCount;
	
											}
	
										}
								}
							}
			
							count++;
							System.out.println("Iterate " + count + "\n");	
	
						}
					}	
				}

			}
		}
	}

	/**
	 * Visualize the clusters.
	 *
	 * Show all points. Each cluster should have a seperate colour.
	 * If you want to use a non-java visualisation solution output the
	 * data you need to a file and in a seperate program do your
	 * visualisation.
	 *
	 * The visualisiation only needs to be 2D. You can use whatever
	 * projection you like.
	 *
	 * If the backend supports some kind of title include the initialisation
	 * strategy and the update strategy in some kind.
	 *
	 * @param points         Array of points to visualise.
	 * @param initialisation strategy to initialisation k-means algo
	 * @param strategy       strategy to update k-means algo
	 **/
	public static void visualize(Point[] points, Initialisation initialisation, Strategy strategy, int n, int k, ArrayList<ArrayList<Integer>> points_it, ArrayList<ArrayList<Point>> centroid_it) {

	Visualizer visualizer = new Visualizer();
	    
	Visualizer.setInformation(
		strategy.name(), 
		initialisation.name(), 
		new Integer(n).toString(), 
		new Integer(k).toString(),
		points,
		points_it,
		centroid_it		
	);
        
        
	Visualizer.main("Visualizer");
	}
}

