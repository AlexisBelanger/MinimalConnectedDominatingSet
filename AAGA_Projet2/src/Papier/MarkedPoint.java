package Papier;

import java.awt.Point;

public class MarkedPoint extends Point {
	private Colour colour;
	private int component;
	
	
	public MarkedPoint(int posX, int posY, Colour c) {
		super(posX, posY);
		this.colour=c;
		component = -1; 
	}
	
	public MarkedPoint(int posX, int posY, Colour c, int comp) {
		super(posX, posY);
		this.colour=c;
		component = comp; 
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
		return "[x: "+this.getX() +" y: "+this.getY()+" couleur = "+this.getColour()+" composant = "+ this.getComp()+"]";
	}
	
}
