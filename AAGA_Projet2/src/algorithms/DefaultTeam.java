package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import utils.Colour;
import utils.MarkedPoint;

public class DefaultTeam {


	
	//---------------------------------------------------------------------------------------------------------------//
	//----------------------------------------Calcul de l'ensemble dominant------------------------------------------//
	//---------------------------------------------------------------------------------------------------------------//
	
	

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
	
	/**
	 * Fonction de voisinnage prenant des markedPoints
	 */
	private ArrayList<MarkedPoint> neighborMarked(MarkedPoint p, ArrayList<MarkedPoint> vertices, int edgeThreshold){
		ArrayList<MarkedPoint> result = new ArrayList<MarkedPoint>();

		for (MarkedPoint point:vertices) if (point.distance(p)<edgeThreshold && !point.equals(p)) result.add((MarkedPoint)point.clone());

		return result;
	}

	/**
	 * verification de la validite du MIS (propriete two hops)
	 */
	public boolean isMISValid(ArrayList<Point> blackpoints, int edgeThreshold) {
		
		for(Point p : blackpoints) {
			if(neighbor(p, blackpoints, edgeThreshold).size()>0) {
				return false;
			}
		}
		return true;
		
	}
	
	/**
	 * parcours en profondeur du graphe, utile pour tester la validité de la solution 
	 */
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
	
	public ArrayList<Point> calculConnectedDominatingSet(ArrayList<Point> points, int edgeThreshold) {

		
		
		
		System.out.println("MARQUAGE....");

		Collections.shuffle(points);

		
	    
	    //Creation du MIS
	    ArrayList<Point> blackPoints = new ArrayList<Point>();
	    ArrayList<Point> greyPoints = new ArrayList<Point>();
	    ArrayList<Point> reste = (ArrayList<Point>) points.clone();
	    ArrayList<ArrayList<MarkedPoint>> composants = new ArrayList<ArrayList<MarkedPoint>>();
	    
	    
	    blackPoints.add(reste.get(0));
	    for(Point p : neighbor(points.get(0), points, edgeThreshold)) {
	    	reste.remove(p);
	    	greyPoints.add(p);

	    }
	    reste.remove(points.get(0));
	    
	    
	    while(!reste.isEmpty()) {
	    	System.out.println(reste.size());
	    	ArrayList<Point> clone_reste = (ArrayList<Point>) reste.clone();
	    	Collections.shuffle(clone_reste);
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
	    			blackPoints.add(p);
	    			for(Point q : voisins) {
	    				
	    				if(reste.contains(q)) {
	    					greyPoints.add(q);
	        				reste.remove(q);
	    				}
	    				
	    			}
	    			reste.remove(p);
	    			
	    		}
	    	}
	    }

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
    			if(p.getColour()==Colour.GREY) {
    											
	    			
    				//identification du nombre de black-blue components différents
    				Set<Integer> diff = new HashSet<Integer>(); 
	    			ArrayList<MarkedPoint> voisins = neighborMarked(p, blacks, edgeThreshold);
	    			for(MarkedPoint q : voisins) {
	    				diff.add(q.getComp());
	    			}
	    			
	    			
	    			if(diff.size()==i || diff.size()>5) {
	    				p.setColour(Colour.BLUE);
	    				blues.add(p);
	    				for(int j =1; j< voisins.size(); j++) {
	    					
	    					if(voisins.get(j).getComp()!=voisins.get(0).getComp()) {
		    					for(MarkedPoint q : composants.get(voisins.get(j).getComp())) {
		    						q.setComponent(voisins.get(0).getComp());
		    					}
		    					
		    					composants.get(voisins.get(0).getComp()).addAll(composants.get((voisins.get(j).getComp())));
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

	    return dom_sol;

	}
	
			//----------------------------------------------------------------//
			//---------------- MIS Papier + Steiner (papier) -----------------//
			//----------------------------------------------------------------//
	
//	public ArrayList<Point> calculConnectedDominatingSet(ArrayList<Point> points, int edgeThreshold) {
//		System.out.println("MARQUAGE....");
//		
//		Collections.shuffle(points);
//		
//
//	    //Creation du MIS
//	    ArrayList<Point> blackPoints = new ArrayList<Point>();
//	    ArrayList<Point> greyPoints = new ArrayList<Point>();
//	    ArrayList<Point> reste = (ArrayList<Point>) points.clone();
//	    ArrayList<ArrayList<MarkedPoint>> composants = new ArrayList<ArrayList<MarkedPoint>>();
//	    
//	    int ind_max = 0;
//	    int n_max = 0;
//	    
//	    for(int i = 0; i< points.size(); i++) {
//	    	ArrayList<Point> v = neighbor(points.get(i), points, edgeThreshold);
//	    	if(v.size()>n_max) {
//	    		ind_max = i;
//	    		n_max = v.size();
//	    	}
//	    }
//	    blackPoints.add(reste.get(ind_max));
//	    for(Point p : neighbor(points.get(ind_max), points, edgeThreshold)) {
//	    	reste.remove(p);
//	    	greyPoints.add(p);
//	
//	    }
//	    reste.remove(points.get(ind_max));
//	    
//	    
//	    while(!reste.isEmpty()) {
//	    	System.out.println(reste.size());
//	    	ArrayList<Point> clone_reste = (ArrayList<Point>) reste.clone();
//	    	Collections.shuffle(clone_reste);
//	    	int indice = 0;
//	    	int nb_voisins_blancs = 0;
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
//	    			int nb_voisins = neighbor(p, clone_reste, edgeThreshold).size();
//	    			if(nb_voisins > nb_voisins_blancs) {
//	    				indice = n;
//	    				nb_voisins_blancs = nb_voisins;
//	    			}
//	    		}
//	    	}
//	    	Point p = clone_reste.get(indice);
//	    	ArrayList<Point> voisins = neighbor(p, points, edgeThreshold);
//	    	blackPoints.add(p);
//			reste.remove(p);
//			for(Point q : voisins) {
//				
//				if(reste.contains(q)) {
//					greyPoints.add(q);
//					reste.remove(q);
//				}
//				
//			}
//	    	
//	    }
//	    
//	    // Verification du MIS
//	    if(isMISValid(blackPoints, edgeThreshold)) {
//	    	int comp = 0;
//		    ArrayList<MarkedPoint> blacks = new ArrayList<MarkedPoint>();
//		    
//		    for(Point p : blackPoints) {
//		    	MarkedPoint point = new MarkedPoint(p.x, p.y, Colour.BLACK, comp);
//		    	composants.add(new ArrayList<MarkedPoint>());
//		    	composants.get(comp).add(point);
//		    	comp++;
//		    	blacks.add(point);
//		    }
//		    
//		    
//		    
//		    ArrayList<MarkedPoint> greys = new ArrayList<MarkedPoint>();
//		    for(Point p : greyPoints) {
//		    	greys.add(new MarkedPoint(p.x, p.y, Colour.GREY));
//		    }
//		    
//		    
//		    
//		    
//		    Set<MarkedPoint> blues = new HashSet<MarkedPoint>();
//		  
//		    System.out.println("MARQUAGE TERMINE");
//		    ArrayList<Point> sol = new ArrayList<Point>();
//		    for(MarkedPoint p : blacks) {
//		    		sol.add(p);
//		    }
//		    System.out.println("CREATION LISTE TERMINE");	    
//		    System.out.println("CREATION DES NOEUDS BLEUS.....");
//		    
//		    
//		    for(int i =5; i>1; i--) {
//		    	System.out.println("i = "+i);
//		    	
//	    		for(MarkedPoint p : greys) {
//	    			if(p.getColour()==Colour.GREY) {// je check si le noeud est gris, histoire de pas risque de retraiter un bleu
//	    											// mais normalement pas besoin
//		    			
//	    				//---------------------------- j'utilise ca pour trouver le nombre de black-blue component differents
//	    				Set<Integer> diff = new HashSet<Integer>(); 
//		    			ArrayList<MarkedPoint> voisins = neighborMarked(p, blacks, edgeThreshold);
//		    			for(MarkedPoint q : voisins) {
//		    				diff.add(q.getComp());
//		    			}
//		    			//----------------------------
//		    			
//		    			if(diff.size()==i || diff.size()>5) {
//		    				p.setColour(Colour.BLUE); // si il a assez de black blue comp, je l'ajoute
//		    				blues.add(p);
//		    				for(int j =1; j< voisins.size(); j++) { // je parcours ses voisins noirs
//		    					
//		    					if(voisins.get(j).getComp()!=voisins.get(0).getComp()) { // pour eviter qu'un voisin du meme composant se reajoute
//			    					for(MarkedPoint q : composants.get(voisins.get(j).getComp())) {
//			    						q.setComponent(voisins.get(0).getComp()); //je MaJ le numero de composant du point
//			    					}
//			    					// et ici j'ajoute tous les points de l'arraylist correspondant au composant du point que je regarde
//			    					// a l'arraylist du composant du premier point noir des voisins
//			    					composants.get(voisins.get(0).getComp()).addAll(composants.get((voisins.get(j).getComp())));
//			    					//et je les retire de leur ancienne liste
//			    					composants.get((voisins.get(j).getComp())).clear(); 
//		    					}   					
//		    				}
//		    			}
//	    			}
//	    		}
//			}
//		    System.out.println("CREATION DES NOEUDS BLEUS TERMINEE");
//
//		    ArrayList<Point> dom_sol = new ArrayList<Point>();
//		    dom_sol.addAll(blues);
//		    dom_sol.addAll(blacks);
//		    
//			ArrayList<MarkedPoint> visit = new ArrayList<MarkedPoint>();
//		    for(Point p : dom_sol) {
//		    	visit.add(new MarkedPoint(p.x, p.y));
//		    }
//		    treeVisit(visit, visit.get(0), edgeThreshold);
//		    for(MarkedPoint p : visit) {
//		    	if(!p.isVisited()) {
//		    		System.out.println("SOLUTION INVALIDE");
//		    		return null;
//		    	}
//		    }
//		    return dom_sol;
//	    }else{
//	    	System.out.println("MIS invalide");
//	    	return null;//la
//	    }		 
//	}	





}





