package comp557.a3;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import com.jogamp.opengl.GL2;

/**
 * Simple face class, which also contains a few extra members for geometry processing.
 */
public class Face {    

    // Brendon Keirle, 260685377
	
    /** Face normal */
    public Vector3d n = new Vector3d();

    /** Some half edge on the face */
    HalfEdge he;
    /** Child vertex in the middle of the face */
    
    /** Area of the face */
    double area;
    
    /** Gradient of the quantities u values (heat) at the vertices */
    Vector3d gradu = new Vector3d();

    /** Center of the face */
    Point3d c = new Point3d();

    /** 
     * Constructs a face from a half edge, and computes the normal, center, and area.
     * @param he
     */
    public Face( HalfEdge he ) {
        this.he = he;
        HalfEdge loop = he;
        do {
            loop.leftFace = this;
            loop = loop.next;
        } while ( loop != he );
        computeNormal();
        computeCenter();
    }
    
    /** 
     * Computes the center of the face
     */
    public void computeCenter() {
    	c.set(0,0,0);
    	HalfEdge loop = he;
    	do {
    		c.add( loop.head.p );
    		loop = loop.next;
    	} while ( loop != he );
    	c.scale( 1.0/3.0 );
    }
    
    /**
     * Computes the non-normalized normal vector
     */
    public void computeNormal() {
    	Point3d p0 = he.head.p;
        Point3d p1 = he.next.head.p;
        Point3d p2 = he.next.next.head.p;
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();
        v1.sub(p1,p0);
        v2.sub(p2,p1);
        n.cross( v1,v2 );
        area = 0.5* n.length();
        n.normalize();
    }
    
    /**
     *  Draw the gradient stored for this face
     *  Draws a line offset from the surface by h times the normal, with the segment
     *  having one end at center of the triangle and the other at a displacement of gradu
     *  scaled by s 
     */
    public void drawGradu( GL2 gl, double h, double s ) {
        // TODO: 7 Write code to draw the line as described in the javadoc above.
    	// use immediate mode (i.e., glBegin and glEnd with GL_LINES, and a call to glVertex
    	// to specfiy each end point.  Note that the line colour and light disabling is done
    	// by the calling function prior to calling this method for all faces.


        Vector3d start = new Vector3d(this.c);
        Vector3d new_normal = new Vector3d(this.n);

        new_normal.scale(h);
        start.add(new_normal);
        if (gradu.length() > 0)
			   gradu.normalize();


        gl.glBegin(GL2.GL_LINES);
        gl.glVertex3d(start.x, start.y, start.z);

        double x = start.x - s*gradu.x;
        double y = start.y - s*gradu.y;
        double z = start.z - s*gradu.z;
        gl.glVertex3d(x, y, z);

        gl.glEnd();

    }
    
}
