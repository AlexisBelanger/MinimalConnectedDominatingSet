package Papier;

import java.awt.Point;

public class MarkedPoint extends Point {
	private Colour colour;
	private int component;
	private boolean visited;
	
	
	public MarkedPoint(int posX, int posY, Colour c) {
		super(posX, posY);
		this.colour=c;
		component = -1; 
		this.visited = false;
	}
	
	public MarkedPoint(int posX, int posY, Colour c, int comp) {
		super(posX, posY);
		this.colour=c;
		component = comp; 
		this.visited = false;
	}
	
	
	public MarkedPoint(int posX, int posY) {
		super(posX, posY);
		this.colour= colour.WHITE;
		component = -1; 
		this.visited = false;
	}
	
	
	public void setColour(Colour c) {
		this.colour=c;
	}
	
	public Colour getColour() {
		return colour;
	}
	
	public int getComp() {
		return component;
	}

	public void setComponent(int c) {
		component = c;
	}
	
	
	public String toString() {
		return "[x: "+this.getX() +" y: "+this.getY()+" couleur = "+this.getColour()+" composant = "+ this.getComp()+ " visited = "+visited+ "]";
	}
	
	public boolean isVisited() {
		return visited;
	}
	
	public void visit() {
		this.visited = true;
	}
	
}
