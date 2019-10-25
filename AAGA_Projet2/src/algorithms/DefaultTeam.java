package algorithms;

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
import java.util.Scanner;
import java.util.Set;

import Papier.Colour;
import Papier.MarkedPoint;

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
	//---------------------------------------------------------------------------------------------------------------//
	//----------------------------------------Calcul de l'arbre de Steiner-------------------------------------------//
	//---------------------------------------------------------------------------------------------------------------//
	
	
	public ArrayList<Point> calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		System.out.println("Debut Steiner");
		
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

		
		System.out.println("FIN PREMIER KRUSKAL");
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
		System.out.println("FIN DEUXIEME KRUSKAL");

		// Solution
		double cost = 0;
		for(Edge e: solution) {
			cost+= e.distance();
		}
		
		ArrayList<Point> sol = new ArrayList<Point>();
		for(Edge e : solution) {
			if(!sol.contains(e.p)) sol.add(e.p);
			if(!sol.contains(e.q)) sol.add(e.q);
		}
		return sol;
		//return new Tree2D(solution.get(0).p, edgesToTree(solution, solution.get(0).p));


	}

	//---------------------------------------------------------------------------------------------------------------//
	//----------------------------------------Calcul de l'ensemble dominant------------------------------------------//
	//---------------------------------------------------------------------------------------------------------------//
	
	
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
	
	
	private ArrayList<MarkedPoint> neighborMarked(MarkedPoint p, ArrayList<MarkedPoint> vertices, int edgeThreshold){
		ArrayList<MarkedPoint> result = new ArrayList<MarkedPoint>();

		for (MarkedPoint point:vertices) if (point.distance(p)<edgeThreshold && !point.equals(p)) result.add((MarkedPoint)point.clone());

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
	
	
	
	public boolean isMISValid(ArrayList<Point> blackpoints, int edgeThreshold) {
		
		for(Point p : blackpoints) {
			if(neighbor(p, blackpoints, edgeThreshold).size()>0) {
				return false;
			}
		}
		return true;
		
	}
	
	
	public void treeVisit(ArrayList<MarkedPoint> points, MarkedPoint p, int edgeThreshold) {
		p.visit();
		for(MarkedPoint q : neighborMarked(p, points, edgeThreshold)) {
			if(!q.isVisited()) {
				treeVisit(points, points.get(points.indexOf(q)), edgeThreshold);
			}
		}
	}
	
	
	

//--------------------------------------------------------------- Algos -----------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//------------------ MIS Basique + Steiner (papier) ----------------------//
	//------------------------------------------------------------------------//
	
//	public ArrayList<Point> calculConnectedDominatingSet(ArrayList<Point> points, int edgeThreshold) {
//
//		
//		
//		
//		System.out.println("MARQUAGE....");
//
//		Collections.shuffle(points);
//
//		
//	    
//	    //Creation du MIS
//	    ArrayList<Point> blackPoints = new ArrayList<Point>();
//	    ArrayList<Point> greyPoints = new ArrayList<Point>();
//	    ArrayList<Point> reste = (ArrayList<Point>) points.clone();
//	    ArrayList<ArrayList<MarkedPoint>> composants = new ArrayList<ArrayList<MarkedPoint>>();
//	    
//	    
//	    blackPoints.add(reste.get(0));
//	    for(Point p : neighbor(points.get(0), points, edgeThreshold)) {
//	    	reste.remove(p);
//	    	greyPoints.add(p);
//
//	    }
//	    reste.remove(points.get(0));
//	    
//	    
//	    while(!reste.isEmpty()) {
//	    	
//	    	ArrayList<Point> clone_reste = (ArrayList<Point>) reste.clone();
//	    	Collections.shuffle(clone_reste);
//	    	for(int n = 0; n<clone_reste.size(); n++) {
//	    		Point p = clone_reste.get(n);
//	    		ArrayList<Point> voisins = neighbor(p, points, edgeThreshold);
//	    		boolean isGreyNeigh = false;
//	    		boolean noBlackNeigh = true;
//	    		
//	    		for(Point q : voisins) {
//	    			if(greyPoints.contains(q)) {
//	    				isGreyNeigh=true;
//	    			}
//	    			if(blackPoints.contains(q)) {
//	    				noBlackNeigh=false;
//	    			}
//	    		}
//	    		if(isGreyNeigh && noBlackNeigh) {
//	    			blackPoints.add(p);
//	    			for(Point q : voisins) {
//	    				
//	    				if(reste.contains(q)) {
//	    					greyPoints.add(q);
//	        				reste.remove(q);
//	    				}
//	    				
//	    			}
//	    			reste.remove(p);
//	    			
//	    		}
//	    	}
//	    }
//
//	    int comp = 0;
//	    ArrayList<MarkedPoint> blacks = new ArrayList<MarkedPoint>();
//	    
//	    for(Point p : blackPoints) {
//	    	MarkedPoint point = new MarkedPoint(p.x, p.y, Colour.BLACK, comp);
//	    	composants.add(new ArrayList<MarkedPoint>());
//	    	composants.get(comp).add(point);
//	    	comp++;
//	    	blacks.add(point);
//	    }
//	    
//	    
//	    
//	    ArrayList<MarkedPoint> greys = new ArrayList<MarkedPoint>();
//	    for(Point p : greyPoints) {
//	    	greys.add(new MarkedPoint(p.x, p.y, Colour.GREY));
//	    }
//	    
//	    
//	    
//	    
//	    Set<MarkedPoint> blues = new HashSet<MarkedPoint>();
//	  
//	    System.out.println("MARQUAGE TERMINE");
//	    ArrayList<Point> sol = new ArrayList<Point>();
//	    for(MarkedPoint p : blacks) {
//	    		sol.add(p);
//	    }
//	    System.out.println("CREATION LISTE TERMINE");	    
//	    System.out.println("CREATION DES NOEUDS BLEUS.....");
//	    
//	    
//	    for(int i =5; i>1; i--) {
//	    	System.out.println("i = "+i);
//	    	
//    		for(MarkedPoint p : greys) {
//    			if(p.getColour()==Colour.GREY) {// je check si le noeud est gris, histoire de pas risque de retraiter un bleu
//    											// mais normalement pas besoin
//	    			
//    				//---------------------------- j'utilise ca pour trouver le nombre de black-blue component differents
//    				Set<Integer> diff = new HashSet<Integer>(); 
//	    			ArrayList<MarkedPoint> voisins = neighborMarked(p, blacks, edgeThreshold);
//	    			for(MarkedPoint q : voisins) {
//	    				diff.add(q.getComp());
//	    			}
//	    			//----------------------------
//	    			
//	    			if(diff.size()==i || diff.size()>5) {
//	    				p.setColour(Colour.BLUE); // si il a assez de black blue comp, je l'ajoute
//	    				blues.add(p);
//	    				for(int j =1; j< voisins.size(); j++) { // je parcours ses voisins noirs
//	    					
//	    					if(voisins.get(j).getComp()!=voisins.get(0).getComp()) { // pour eviter qu'un voisin du meme composant se reajoute
//		    					for(MarkedPoint q : composants.get(voisins.get(j).getComp())) {
//		    						q.setComponent(voisins.get(0).getComp()); //je MaJ le numero de composant du point
//		    					}
//		    					// et ici j'ajoute tous les points de l'arraylist correspondant au composant du point que je regarde
//		    					// a l'arraylist du composant du premier point noir des voisins
//		    					composants.get(voisins.get(0).getComp()).addAll(composants.get((voisins.get(j).getComp())));
//		    					//et je les retire de leur ancienne liste
//		    					composants.get((voisins.get(j).getComp())).clear(); 
//	    					}   					
//	    				}
//	    			}
//    			}
//    		}
//		}
//	    System.out.println("CREATION DES NOEUDS BLEUS TERMINEE");
//	    
//	   
//	    ArrayList<Point> dom_sol = new ArrayList<Point>();
//	    dom_sol.addAll(blues);
//	    dom_sol.addAll(blacks);
//
//	    return dom_sol;
//
//	}
	
			//----------------------------------------------------------------//
			//--------------- MIS Ameliore + Steiner (papier) ----------------//
			//----------------------------------------------------------------//
	
	public ArrayList<Point> calculConnectedDominatingSet(ArrayList<Point> points, int edgeThreshold) {
		System.out.println("MARQUAGE....");
		
		Collections.shuffle(points);
		

	    //Creation du MIS
	    ArrayList<Point> blackPoints = new ArrayList<Point>();
	    ArrayList<Point> greyPoints = new ArrayList<Point>();
	    ArrayList<Point> reste = (ArrayList<Point>) points.clone();
	    ArrayList<ArrayList<MarkedPoint>> composants = new ArrayList<ArrayList<MarkedPoint>>();
	    
	    int ind_max = 0;
	    int n_max = 0;
	    
	    for(int i = 0; i< points.size(); i++) {
	    	ArrayList<Point> v = neighbor(points.get(i), points, edgeThreshold);
	    	if(v.size()>n_max) {
	    		ind_max = i;
	    		n_max = v.size();
	    	}
	    }
	    blackPoints.add(reste.get(ind_max));
	    for(Point p : neighbor(points.get(ind_max), points, edgeThreshold)) {
	    	reste.remove(p);
	    	greyPoints.add(p);
	
	    }
	    reste.remove(points.get(ind_max));
	    
	    
	    while(!reste.isEmpty()) {
	    	
	    	ArrayList<Point> clone_reste = (ArrayList<Point>) reste.clone();
	    	Collections.shuffle(clone_reste);
	    	int indice = 0;
	    	int nb_voisins_blancs = 0;
	    	for(int n = 0; n<clone_reste.size(); n++) {
	    		Point p = clone_reste.get(n);
	    		ArrayList<Point> voisins = neighbor(p, points, edgeThreshold);
	    		boolean isGreyNeigh = false;
	    		boolean noBlackNeigh = true;
	    		
	    		for(Point q : voisins) {
	    			if(greyPoints.contains(q)) {
	    				isGreyNeigh=true;
	    			}
	    			if(blackPoints.contains(q)) {
	    				noBlackNeigh=false;
	    			}
	    		}
	    		if(isGreyNeigh && noBlackNeigh) {
	    			int nb_voisins = neighbor(p, clone_reste, edgeThreshold).size();
	    			if(nb_voisins > nb_voisins_blancs) {
	    				indice = n;
	    				nb_voisins_blancs = nb_voisins;
	    			}
	    		}
	    	}
	    	Point p = clone_reste.get(indice);
	    	ArrayList<Point> voisins = neighbor(p, points, edgeThreshold);
	    	blackPoints.add(p);
			reste.remove(p);
			for(Point q : voisins) {
				
				if(reste.contains(q)) {
					greyPoints.add(q);
					reste.remove(q);
				}
				
			}
	    	
	    }
	    
	    // Verification du MIS
	    if(isMISValid(blackPoints, edgeThreshold)) {
	    	int comp = 0;
		    ArrayList<MarkedPoint> blacks = new ArrayList<MarkedPoint>();
		    
		    for(Point p : blackPoints) {
		    	MarkedPoint point = new MarkedPoint(p.x, p.y, Colour.BLACK, comp);
		    	composants.add(new ArrayList<MarkedPoint>());
		    	composants.get(comp).add(point);
		    	comp++;
		    	blacks.add(point);
		    }
		    
		    
		    
		    ArrayList<MarkedPoint> greys = new ArrayList<MarkedPoint>();
		    for(Point p : greyPoints) {
		    	greys.add(new MarkedPoint(p.x, p.y, Colour.GREY));
		    }
		    
		    
		    
		    
		    Set<MarkedPoint> blues = new HashSet<MarkedPoint>();
		  
		    System.out.println("MARQUAGE TERMINE");
		    ArrayList<Point> sol = new ArrayList<Point>();
		    for(MarkedPoint p : blacks) {
		    		sol.add(p);
		    }
		    System.out.println("CREATION LISTE TERMINE");	    
		    System.out.println("CREATION DES NOEUDS BLEUS.....");
		    
		    
		    for(int i =5; i>1; i--) {
		    	System.out.println("i = "+i);
		    	
	    		for(MarkedPoint p : greys) {
	    			if(p.getColour()==Colour.GREY) {// je check si le noeud est gris, histoire de pas risque de retraiter un bleu
	    											// mais normalement pas besoin
		    			
	    				//---------------------------- j'utilise ca pour trouver le nombre de black-blue component differents
	    				Set<Integer> diff = new HashSet<Integer>(); 
		    			ArrayList<MarkedPoint> voisins = neighborMarked(p, blacks, edgeThreshold);
		    			for(MarkedPoint q : voisins) {
		    				diff.add(q.getComp());
		    			}
		    			//----------------------------
		    			
		    			if(diff.size()==i || diff.size()>5) {
		    				p.setColour(Colour.BLUE); // si il a assez de black blue comp, je l'ajoute
		    				blues.add(p);
		    				for(int j =1; j< voisins.size(); j++) { // je parcours ses voisins noirs
		    					
		    					if(voisins.get(j).getComp()!=voisins.get(0).getComp()) { // pour eviter qu'un voisin du meme composant se reajoute
			    					for(MarkedPoint q : composants.get(voisins.get(j).getComp())) {
			    						q.setComponent(voisins.get(0).getComp()); //je MaJ le numero de composant du point
			    					}
			    					// et ici j'ajoute tous les points de l'arraylist correspondant au composant du point que je regarde
			    					// a l'arraylist du composant du premier point noir des voisins
			    					composants.get(voisins.get(0).getComp()).addAll(composants.get((voisins.get(j).getComp())));
			    					//et je les retire de leur ancienne liste
			    					composants.get((voisins.get(j).getComp())).clear(); 
		    					}   					
		    				}
		    			}
	    			}
	    		}
			}
		    System.out.println("CREATION DES NOEUDS BLEUS TERMINEE");

		    ArrayList<Point> dom_sol = new ArrayList<Point>();
		    dom_sol.addAll(blues);
		    dom_sol.addAll(blacks);
		    
			ArrayList<MarkedPoint> visit = new ArrayList<MarkedPoint>();
		    for(Point p : dom_sol) {
		    	visit.add(new MarkedPoint(p.x, p.y));
		    }
		    treeVisit(visit, visit.get(0), edgeThreshold);
		    for(MarkedPoint p : visit) {
		    	if(!p.isVisited()) {
		    		System.out.println("SOLUTION INVALIDE");
		    		return null;
		    	}
		    }
		    return dom_sol;
	    }else{
	    	System.out.println("MIS invalide");
	    	return null;//la
	    }		 
	}	
}





