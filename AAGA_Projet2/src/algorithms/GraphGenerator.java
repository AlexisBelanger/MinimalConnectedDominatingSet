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

public class GraphGenerator {

		

		private static void saveToFile(String filename,ArrayList<Point> result){
			int index=0;
			try {
				while(true){
					BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(filename+Integer.toString(index)+".points")));
					try {
						input.close();
					} catch (IOException e) {
						System.err.println("I/O exception: unable to close "+filename+Integer.toString(index)+".points");
					}
					index++;
				}
			} catch (FileNotFoundException e) {
				printToFile(filename+Integer.toString(index)+".points",result);
			}
		}
		private static void printToFile(String filename,ArrayList<Point> points){
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
	
	
		
		public static void main(String[] args) {
			
			
			int[] taille = {500,1000,2500,5000};
			
			for(int i = 0; i<taille.length; i++) {
				ArrayList<Point> ensemble = new ArrayList<Point>();
				for(int j = 0; j<taille[i]; j++) {
					
					
					Point p = new Point(new Point((int)(Math.random()*800)+300,(int)(Math.random()*700)+100));
					while(ensemble.contains(p)) {
						p = new Point(new Point((int)(Math.random()*800)+300,(int)(Math.random()*700)+100));
					}
					ensemble.add(p);
				}
				String filename = "input_nb_points_"+taille[i];
				saveToFile(filename,ensemble);
				
				System.out.println(i+" done");
			}
			
			
			
			
			
		}
	

}
