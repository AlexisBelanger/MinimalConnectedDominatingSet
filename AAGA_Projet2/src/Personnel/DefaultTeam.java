package Personnel;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import Papier.DefaultTeam.Edge;
import Papier.DefaultTeam.TaggedPoint;

public class DefaultTeam {
	//Sort function : 

	private ArrayList<Edge> sort(ArrayList<Edge> edges) {
		if (edges.size()==1) return edges;

		ArrayList<Edge> left = new ArrayList<Edge>();
		ArrayList<Edge> right = new ArrayList<Edge>();
		int n=edges.size();
		for (int i=0;i<n/2;i++) { left.add(edges.remove(0)); }
		while (edges.size()!=0) { right.add(edges.remove(0)); }
		left = sort(left);
		right = sort(right);

		ArrayList<Edge> result = new ArrayList<Edge>();
		while (left.size()!=0 || right.size()!=0) {
			if (left.size()==0) { result.add(right.remove(0)); continue; }
			if (right.size()==0) { result.add(left.remove(0)); continue; }
			if (left.get(0).distance() < right.get(0).distance()) result.add(left.remove(0));
			else result.add(right.remove(0));
		}
		return result;
	}

	private boolean isLeaf(ArrayList<Edge>list, Point p) {
		ArrayList<Edge> cpt = new ArrayList<Edge>();

		for(Edge tmp : list) {
			if(tmp.p == p || tmp.q == p) {
				cpt.add(tmp);
			}
			if(cpt.size()==2) {
				return false;
			}
		}
		return true;

	}

	public int distance(Point p, Point q) {
		return  (int)(Math.sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y)));

	}



	public class Edge{
		private Point p;
		private Point q;

		Edge(Point p, Point q){
			this.p = p;
			this.q = q;
		}

		public Point getP() {
			return p;
		}

		public Point getQ() {
			return q;
		}

		public double distance() {
			return  (Math.sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y)));

		}
	}


	public class TaggedPoint{
		Point p;
		int tag;

		public TaggedPoint(Point p, int tag) {
			this.p = p;
			this.tag = tag;
		}

		public int getTag() {
			return tag;
		}

		public void setTag(int tag) {
			this.tag = tag;
		}

		public Point getP() {
			return p;
		}
	}

	public boolean isSolution(Edge e, ArrayList<TaggedPoint> l) {
		Point p1 = e.p;
		Point p2 = e.q;
		int x=Integer.MIN_VALUE;
		int y=Integer.MIN_VALUE;
		int cpt = 0;
		for (TaggedPoint tp : l) {
			if(tp.p.equals(p1)) {
				x = tp.tag;
				cpt++;
			}
			if(tp.p.equals(p2)) {
				y = tp.tag;
				cpt++;
			}
			if(cpt == 2) {
				break;
			}
		}
		if(x == y) {
			return false;
		}
		for (TaggedPoint tp : l) {
			if(tp.tag == x) {
				tp.tag = y;
			}
		}
		return true;

	}


	public ArrayList<Tree2D> edgesToTree(ArrayList<Edge> l, Point root) {

		ArrayList<Tree2D> subtree = new ArrayList<Tree2D>();
		ArrayList<Edge> lbis = new ArrayList<Edge>(l);
		//lbis.addAll(l);
		for(Edge e : l) {
			if(e.p == root) {
				lbis.remove(e);
				if(!lbis.isEmpty()) {
					subtree.add(new Tree2D(e.q, edgesToTree(lbis, e.q)));	  
				}else{
					subtree.add(new Tree2D(e.q, new ArrayList<Tree2D>()));
				}
			}
			if(e.q == root) {
				lbis.remove(e);
				if(!lbis.isEmpty()) {
					subtree.add(new Tree2D(e.p, edgesToTree(lbis, e.p)));

				}else {
					subtree.add(new Tree2D(e.p, new ArrayList<Tree2D>()));

				}
			}
		}
		return subtree;

	}


	// Calcul du chemin minimum d'un point a un autre

	public double dist(ArrayList<Point> points,int[][] paths, int i, int j) {
		if(paths[i][j]==i) {
			return 0;
		}else {
			double d = 1;
			return d + dist(points,paths, paths[i][j], j );
		}

	}

	/**
	 * 
	 * @param points set of all points
	 * @param paths
	 * @param i first point index
	 * @param j second point index
	 * @return list of indexes of points traveled to go from i to j
	 */



	public ArrayList<Integer> pathFinder(ArrayList<Point> points,int[][] paths, int i, int j) {
		ArrayList<Integer> liste = new ArrayList<Integer>();
		while(paths[i][j]!=j) {
			liste.add(paths[i][j]);
			int tmp = paths[i][j];
			i=j;
			j=tmp;
		}
		return liste;

	}

	public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		int[][] paths=new int[points.size()][points.size()];
		double[][] distances=new double[points.size()][points.size()];
		for(int i=0;i<distances.length;i++) {
			for(int j=0;j<distances.length;j++) {
				double d = (Point.distance(points.get(i).x, points.get(i).y, points.get(j).x, points.get(j).y));
				if(d>= edgeThreshold) {
					distances[i][j] = Double.MAX_VALUE;
					paths[i][j] = j;
				}else{
					distances[i][j] = d;
					paths[i][j] = j;
				}
			}
		}
		for (int k=0;k<paths.length;k++) {
			for (int i=0;i<paths.length;i++) {
				for(int j=0;j<paths.length;j++) {
					double d = distances[i][k]+distances[k][j];
					if(d < distances[i][j] && d > 0) {
						distances[i][j] = d;
						paths[i][j] = paths[i][k];
					}
				}
			}
		}


		return paths;
	}
	// __________________________________________ //


	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[][] paths = calculShortestPaths(points, edgeThreshold);


		int[] indexHitPoint = new int[hitPoints.size()];

		for (int i = 0; i<hitPoints.size(); i++) {
			indexHitPoint[i]= points.indexOf(hitPoints.get(i));
		}

		//First Kruskal

		ArrayList<Point> points2 = new ArrayList<Point>(hitPoints);
		ArrayList<Edge> edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));

		edges = sort(edges);

		ArrayList<Edge> solution = new ArrayList<Edge>();
		ArrayList<TaggedPoint> TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}

		for(Edge a :edges) {
			if(points2.isEmpty()) {
				break;
			}
			if( isSolution(a, TaggedP) ) {
				points2.remove(a.p);
				solution.add(a);
			}
		}

		// End

		//Replacement of the edge by a shorter path


		for(Edge e : solution) {
			Point tmp1 = e.p;
			Point tmp2 = e.q;
			int i1 = points.indexOf(tmp1);
			int i2 = points.indexOf(tmp2);
			ArrayList<Integer> l = pathFinder(points, paths, i1, i2); 
			for(int i=0; i<l.size(); i++) {
				if(!(hitPoints.contains(points.get(l.get(i))))){
					hitPoints.add(points.get(l.get(i)));
				}
			}  
		}

		points2 = new ArrayList<Point>(hitPoints);
		edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		// Second Kruskal




		edges = sort(edges);

		solution = new ArrayList<Edge>();
		TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}
		int total = 0;
		for(Edge e :edges) {
			if(points2.isEmpty() ) {
				break;
			}
			if( isSolution(e, TaggedP) ) {
				total+= dist(points, paths, points.indexOf(e.p), points.indexOf(e.q));
				points2.remove(e.p);
				solution.add(e);
			}
		}

		// Solution
		double cost = 0;
		for(Edge e: solution) {
			cost+= e.distance();
		}
		System.out.println(cost);
		return new Tree2D(solution.get(0).p, edgesToTree(solution, solution.get(0).p));


	}

	//__________________________Steiner width budget______________________________//

	public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[][] paths = calculShortestPaths(points, edgeThreshold);
		ArrayList<Point> savedHP = new ArrayList<Point>(hitPoints);

		int[] indexHitPoint = new int[hitPoints.size()];

		for (int i = 0; i<hitPoints.size(); i++) {
			indexHitPoint[i]= points.indexOf(hitPoints.get(i));
		}

		//First Kruskal


		ArrayList<Point> points2 = new ArrayList<Point>(hitPoints);
		ArrayList<Edge> edges = new ArrayList<Edge>();	
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		edges = sort(edges);


		ArrayList<Edge> solution = new ArrayList<Edge>();
		ArrayList<TaggedPoint> TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}

		for(Edge a :edges) {
			if(points2.isEmpty()) {
				break;
			}
			if( isSolution(a, TaggedP) ) {
				points2.remove(a.p);
				solution.add(a);
			}
		}



		// End

		//Replacement of the edge by a shorter path


		for(Edge e : solution) {
			Point tmp1 = e.p;
			Point tmp2 = e.q;
			int i1 = points.indexOf(tmp1);
			int i2 = points.indexOf(tmp2);
			ArrayList<Integer> l = pathFinder(points, paths, i1, i2); 
			for(int i=0; i<l.size(); i++) {
				if(!(hitPoints.contains(points.get(l.get(i))))){
					hitPoints.add(points.get(l.get(i)));
				}
			}  
		}

		points2 = new ArrayList<Point>(hitPoints);
		edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		// Second Kruskal




		edges = sort(edges);

		solution = new ArrayList<Edge>();
		TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}
		for(Edge e :edges) {
			if(points2.isEmpty() ) {
				break;
			}
			if( isSolution(e, TaggedP) ) {
				//cost+= dist(points, paths, points.indexOf(e.p), points.indexOf(e.q));
				points2.remove(e.p);
				solution.add(e);
			}
		}
		//End kruskal_______________________________________________//


		double cost = 0;
		for(Edge e: solution) {
			cost+= e.distance();
		}

		Edge root = solution.get(0);
		ArrayList<Edge> removedE = new ArrayList<Edge>();
		ArrayList<Edge> copySol = new ArrayList<Edge>(solution);
		System.out.println(cost);

		while(cost > 1664) {
			Point MostFar = null;
			double distance = 0;

			for (Point p : hitPoints) {
				if(distance(root.p,p)>distance) {
					distance = dist(points, paths, points.indexOf(root.p), points.indexOf(p));
					MostFar = p;
				}
				if(distance(root.q,p)>distance) {
					distance = dist(points, paths, points.indexOf(root.q), points.indexOf(p));
					MostFar = p;
				}
			}
			Edge toRemove = null;

			for(Edge e : solution) {
				if((isLeaf(solution,e.p) && e.p == MostFar) || isLeaf(solution,e.q) &&  e.q == MostFar) {
					cost -= e.distance();
					toRemove = e;
				}
			}
			solution.remove(toRemove);
			copySol.remove(toRemove);
			hitPoints.remove(MostFar);
			boolean aSuppr = true;

			while(aSuppr) { //delete useless edges
				boolean test = false;
				Edge tmpe = null;
				for(Edge e : copySol) {
					if(isLeaf(copySol,e.p)|| isLeaf(copySol,e.q)) {
						if(!savedHP.contains(e.p) || !savedHP.contains(e.q)  ) {
							test = true;
							tmpe = e;
							break;
						}
					}
				}
				if(test) {
					removedE.add(tmpe);
					copySol.remove(tmpe);
					cost -= tmpe.distance();
				}else {
					break;
				}
			}

		}	

		double costsol2 = 0;
		for(Edge e: copySol) {
			costsol2+=e.distance();
		}

		System.out.println("\n"+costsol2);
		
		return new Tree2D(copySol.get(0).p, edgesToTree(copySol, copySol.get(0).p));

	}
	
	
	private static int nbfile = 0;

	public boolean isValid(ArrayList<Point> ens, ArrayList<Point> sol, int edgeThreshold) {
		ArrayList<Point> points = (ArrayList<Point>) ens.clone();
		for(Point p: sol) {
			for(Point q: neighbor(p,points,edgeThreshold)) {
				points.remove(q);
			}
			points.remove(p);
		}
		return points.isEmpty();
	}

	private ArrayList<Point> neighbor(Point p, ArrayList<Point> vertices, int edgeThreshold){
		ArrayList<Point> result = new ArrayList<Point>();

		for (Point point:vertices) if (point.distance(p)<edgeThreshold && !point.equals(p)) result.add((Point)point.clone());

		return result;
	}

	public ArrayList<Point> improve(ArrayList<Point> points, ArrayList<Point> sol, int edgeThreshold){
		ArrayList<Point> reste = ((ArrayList<Point>) points.clone());
		Collections.shuffle(reste);
		Collections.shuffle(sol);
		ArrayList<Point> clone_sol = (ArrayList<Point>)sol.clone();


		for(int i = 0; i<sol.size(); i++) {
			Point p = sol.get(i);

			for(int j=i+1; j< sol.size(); j++) {
				Point q = sol.get(j);

				if(p.distance(q)>2.5*edgeThreshold) continue;

				for(Point r : reste) {
					if(r.distance(q)>2*edgeThreshold || r.distance(p)>2*edgeThreshold) continue;

					clone_sol.remove(p);
					clone_sol.remove(q);
					clone_sol.add(r);

					if(isValid(points,clone_sol, edgeThreshold)) {

						return clone_sol;

					}else {

						clone_sol.add(p);
						clone_sol.add(q);
						clone_sol.remove(r);
					}
				}
			}
		}
		return clone_sol;
	}




	public ArrayList<Point> calculDominatingSet(ArrayList<Point> points, int edgeThreshold) {

		int pass = 0;
		/*for(int i = 0; i < pass; i++)nbfile++;
		String filename1 = "result" + nbfile + ".points";
		ArrayList<Point> old_res = readFromFile(filename1);
		if(!old_res.isEmpty()) {
			nbfile++;
			return old_res;
		}*/

		System.out.println("begin");
		ArrayList <Point> final_sol = new ArrayList<Point>();
		int nb_points = 1000;

		int essai = 0;
		while(true) {
			ArrayList<Point> sol = new ArrayList<Point>();
			ArrayList<Point> tmp = new ArrayList<Point>();
			ArrayList<Point> res = new ArrayList<Point>();

			tmp = (ArrayList<Point>) points.clone();

			while(!isValid(points,sol, edgeThreshold)) {

				ArrayList<Point> clone = (ArrayList<Point>)tmp.clone();
				
				for(Point p: clone) {
					if(neighbor(p,tmp,edgeThreshold).size()==0) {
						sol.add(p);
						tmp.remove(p);
						res.add(p);
					}
				}

				Collections.shuffle(tmp);
				Point s = new Point();
				double nbminN = 1000;
				for(Point p : tmp) {
					int nrs = neighbor(p,tmp,edgeThreshold).size();

					if(nrs<nbminN) {
						nbminN = nrs;
						s = p;
					}
				}
				//System.out.println("etape 2");
				Point s2 = new Point();
				int neimax = -1;
				//System.out.println(neighbor(s,tmp,edgeThreshold).size());
				for(Point q : neighbor(s,tmp,edgeThreshold)) {
					int n = neighbor(q,tmp,edgeThreshold).size();
					if(n>neimax) {
						s2 = q;
						neimax = n;
						//	System.out.println("ah");
					}
				}



				//s = tmp.get(0);
				sol.add(s2);

				for(Point p_tmp: neighbor(s2,tmp,edgeThreshold)) {
					if(tmp.contains(p_tmp)) {
						tmp.remove(p_tmp);
						res.add(p_tmp);
					}
				}
				tmp.remove(s2);

			}
			//Iterations de improve
			int i = 0;
			while(true) {
				ArrayList<Point> resu = improve(points, sol, edgeThreshold);

				if(resu.size()<sol.size()) {
					sol.clear();
					sol.addAll(resu);
				}else {
					break;
				}
				i++;
			}
			if(sol.size()<nb_points) {
				final_sol = sol;
				nb_points = sol.size();
				essai = 0;
			}
			if(essai == 1000 || nb_points<69) {
				break;
			}else {
				essai++;
			}
			System.out.println("essai : "+ essai + "	best score : "+ nb_points);
		}
		//Fin algo
		System.out.println("best score: "+nb_points);
		String filename = "result"+nbfile;
		saveToFile(filename,final_sol);
		nbfile++;
		if(nbfile == 100) nbfile = 0;
		return final_sol;
	}

	//FILE PRINTER
	private void saveToFile(String filename,ArrayList<Point> result){
		int index=0;
		try {
			while(true){
				BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(filename+".points")));
				try {
					input.close();
				} catch (IOException e) {
					System.err.println("I/O exception: unable to close "+filename+Integer.toString(index)+".points");
				}
				index++;
			}
		} catch (FileNotFoundException e) {
			printToFile(filename+".points",result);
		}
	}
	private void printToFile(String filename,ArrayList<Point> points){
		try {
			PrintStream output = new PrintStream(new FileOutputStream(filename));
			int x,y;
			for (Point p:points) output.println(Integer.toString((int)p.getX())+" "+Integer.toString((int)p.getY()));
			output.close();
		} catch (FileNotFoundException e) {
			System.err.println("I/O exception: unable to create "+filename);
		}
	}

	//FILE LOADER
	private ArrayList<Point> readFromFile(String filename) {
		String line;
		String[] coordinates;
		ArrayList<Point> points=new ArrayList<Point>();
		try {
			BufferedReader input = new BufferedReader(
					new InputStreamReader(new FileInputStream(filename))
					);
			try {
				while ((line=input.readLine())!=null) {
					coordinates=line.split("\\s+");
					points.add(new Point(Integer.parseInt(coordinates[0]),
							Integer.parseInt(coordinates[1])));
				}
			} catch (IOException e) {
				System.err.println("Exception: interrupted I/O.");
			} finally {
				try {
					input.close();
				} catch (IOException e) {
					System.err.println("I/O exception: unable to close "+filename);
				}
			}
		} catch (FileNotFoundException e) {
			System.err.println("Input file not found.");
		}
		return points;
	}
	
	
	
	
	
	


}





