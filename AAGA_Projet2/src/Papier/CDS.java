package Papier;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;

import javax.lang.model.element.NestingKind;

public class CDS {
	
	public static ArrayList<MarkedPoint> neighbor(MarkedPoint p, ArrayList<MarkedPoint> vertices, int edgeThreshold){
		ArrayList<MarkedPoint> result = new ArrayList<MarkedPoint>();

		for (Point point:vertices) if (point.distance(p)<edgeThreshold && !point.equals(p)) result.add((MarkedPoint)point.clone());

		return result;
	}

	
	
	
	
	
	public static void main(String[] args) {
		
		
		int edgeThreshold = 100;
		String test = "test_1";
		
		ArrayList<MarkedPoint> points = new ArrayList<MarkedPoint>();
		
		//Lecture du fichier de coordonnées et ajout des points à la liste
		File f = new File("ressources/"+test+".txt");
	    try {
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String coords;
			while ((coords = reader.readLine()) != null) {
			    String[] c = coords.split(" ");
				points.add(new MarkedPoint(Integer.parseInt(c[0]), Integer.parseInt(c[1]), Colour.WHITE));
			} 
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	    
	    
	    
	    //Creation du MIS
	    ArrayList<MarkedPoint> dominatingPoints = new ArrayList<MarkedPoint>();
	    ArrayList<MarkedPoint> dominatedPoints = new ArrayList<MarkedPoint>();
	    
	    for(MarkedPoint p : points) {
	    	if(p.getColour()== Colour.WHITE) {
	    		p.setColour(Colour.BLACK);
	    		dominatingPoints.add(p);
	    		for(MarkedPoint q : neighbor(p, points, edgeThreshold)) {
	    			points.get(points.indexOf(q)).setColour(Colour.GREY);	
	    			dominatedPoints.add(q);
	    		}
	    	}
	    }
	    
	    
	    for(int i = 5; i>1; i++) {
	    	
	    }
	    
	    
	    
	    
	   
		
		
		
		
		
	}

}
