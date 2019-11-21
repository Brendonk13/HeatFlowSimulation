package comp557.a3;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

/**
 * Half edge data structure.
 * Maintains a list of faces (i.e., one half edge of each) to allow
 * for easy display of geometry.
 */
public class HEDS {
    
    // Brendon Keirle, 260685377

	// arraylist of faces for each vertex
    Map<Integer,List<Integer>> faces_per_vertex = new TreeMap<Integer,List<Integer>>();

    /** List of faces */
    ArrayList<Face> faces = new ArrayList<Face>();

    /** List of vertices */
    ArrayList<Vertex> vertices;

    /** Convenience member for keeping track of half edges you make or need */
    Map<String,HalfEdge> halfEdges = new TreeMap<String,HalfEdge>();

    /**
     * Builds a half edge data structure from the polygon soup   
     * @param soup
     */
    public HEDS( PolygonSoup soup ) {
        vertices = soup.vertexList;
        ArrayList<int[]> used_face_list;
        
		used_face_list = new ArrayList<int[]>();
		triangulate(used_face_list, soup.faceList);

        for (int[] face: used_face_list) {

        	// here I set the attributes: (head, next) during each iteration.
        	HalfEdge first = new HalfEdge();
        	first.head = vertices.get(face[0]);
        	first.head.he = first;
        	HalfEdge prev = first;
        	HalfEdge he = null;
            // i = 1 since we already did first^
        	for (int i=1; i<face.length; i++) {
        		Vertex v = vertices.get(face[i]);
        		he = new HalfEdge();
        		he.head = v;
                v.he = he;
                v.index = face[i]; 
        		prev.next = he;
        		if (i == 1) {
                    halfEdges.put(face[face.length-1] + "," + face[0], prev);
                } else {
                    halfEdges.put(face[i-2] + "," + face[i-1], prev);
                }

                // for next iteration
        		prev = he;
        	}   	     
            // finish linking initial HE: first
        	he.next = first;
        	halfEdges.put(face[face.length-2] + "," + face[face.length-1], he);
        }


        // Find the twin of all half edges and link them
        // I store the twin of "i,j" at key: "j,i"
        for (int[] face: used_face_list) {

        	for (int i=0; i<face.length-1; i++) {
                HalfEdge he = halfEdges.get(face[i] + "," + face[i+1]);
        		he.twin = halfEdges.get(face[i+1] + "," + face[i]);
        		halfEdges.put(face[i] + "," + face[i+1], he);
        	}
    		HalfEdge he = halfEdges.get(face[face.length - 1] + "," + face[0]);
    		he.twin = halfEdges.get(face[0] + "," + face[face.length - 1]);
    		halfEdges.put(face[face.length - 1] + "," + face[0], he);

    		// this sets he.leftFace for all H.E's in the face
			faces.add(new Face(he));
        }

        // TODO: 3 Compute vertex normals
        computeVertexNormals();
    }


    private void triangulate(ArrayList<int[]> new_list, ArrayList<int[]> old_list){
        int ct;
        int idx = 0;
		int new_vert_idx;
        for (int[] og_face : old_list){
			// if length = 3, this is a triangle, don't need to add a new face
        	if (og_face.length == 3) {
                new_list.add(idx, new int[]{og_face[0], og_face[1], og_face[2]} );
                idx += 1;
                continue;
            }
           
        	vertices.add(newVertInMiddleFace(og_face));
        	
        	// index of vertex we just appended to vertices list.
            new_vert_idx = vertices.size() - 1;
            ct = 0;
            while(og_face.length - ct > 0){
                new_list.add(idx, new int[]{new_vert_idx, og_face[ct], og_face[(ct+1) % og_face.length] });
                ct++;
                idx++;
            }
        }
    }

    // this computes the center of the face and then inserts new vertex here.
    private Vertex newVertInMiddleFace(int[] og_face) {
        Vertex new_vert = new Vertex();
        new_vert.p.set(0,0,0);
        for (int vert_ind : og_face) {
            // new_vert.p = sum(points in face)
            new_vert.p.add(vertices.get(vert_ind).p);
        }
        // center is the sum of all points, scaled by the number of points in the polygon
        new_vert.p.scale( 1.0/og_face.length );
        // index is size since we have yet to add it to vertices list
        new_vert.index = vertices.size(); 
        return new_vert;
    }

    
    private void computeVertexNormals() {
    	Vector3d avg_normal;
		HalfEdge curr;
    	for (Vertex curr_v : vertices) {
    		avg_normal = new Vector3d(0,0,0);
    		curr = curr_v.he;
    		do {
				avg_normal.add(curr.leftFace.n);
				curr = curr.next.twin;
    		} while(curr != curr_v.he);
    		avg_normal.normalize();
    		curr_v.n = avg_normal;
    	}
    }


    /**
     * Helper function for creating a half edge, and pairing it up with its twin if the
     * twin half edge has already been created.
     * @param soup 
     * @param i tail vertex index
     * @param j head vertex index
     * @return the half edge, paired with its twin if the twin was created.
     * @throws Exception
     */
    private HalfEdge createHalfEdge( PolygonSoup soup, int i, int j ) throws Exception {
        String p = i+","+j;
        if ( halfEdges.containsKey( p ) ){
            throw new Exception("non orientable manifold");
        }
        String twin = j+","+i;
        HalfEdge he = new HalfEdge();
        he.head = soup.vertexList.get(j);
        he.head.he = he; // make sure the vertex has at least one half edge that points to it.
        he.twin = halfEdges.get( twin );
        if ( he.twin != null ) he.twin.twin = he;
        halfEdges.put( p, he );        
        return he;        
    }    

    /** 
     * Reset the solutions for heat and distance
     */
    public void resetHeatAndDistanceSolution() {
    	for ( Vertex v : vertices ) {
    		v.u0 = v.constrained? 1 : 0;
    		v.ut = v.u0;
    		v.phi = 0;
    	}
    }

    /** 
     * Perform a specified number of projected Gauss-Seidel steps of the heat diffusion equation.
     * The current ut values stored in the vertices will be refined by this call.
     * @param GSSteps number of steps to take
     * @param t solution time
     */
    public void solveHeatFlowStep( int GSSteps, double t ) {    	
    	// Solve (A - t LC) u_t = u_0 with constrained vertices holding their ut value fixed
    	// Note that this is a backward Euler step of the heat diffusion.
        // double[] left_m = new double[vertices.size()];
        // double[] areas;
    	for ( int i = 0; i < GSSteps; i++ ) {
    		for ( Vertex v : vertices ) {
    			if ( v.constrained ) continue;  // do nothing for the constrained vertex!
    			// TODO: 7 write inner loop code for the PGS heat solve

                double denom = v.area - t*v.LCii;
                double RHS = v.area*v.u0;
                // double RHS = 0;
                HalfEdge curr = v.he.next;
                int j = 0;

                do {
                    RHS += (t*v.LCij[j]) * curr.head.ut;
                    j++;
                    curr = curr.twin.next;

                } while(curr != v.he.next);
                v.ut = RHS / denom;
    		}	
    	}
    }

    /**
     * Compute the gradient of heat at each face
     */
    public void updateGradu() {
    	// TODO: 8 update the gradient of u from the heat values, i.e., f.gradu for each Face f

        Vertex opp_point;
        Vector3d tmp = new Vector3d();
        for (Face f : faces){
            HalfEdge curr = f.he;
            f.gradu.set(0,0,0);

            for(int i=0; i<3; i++){
                opp_point = curr.head;
                // tmp = second_he - fst_he;
                tmp = (Vector3d) curr.prev().getVecToHead();
                tmp.negate();
                
                tmp.cross(f.n, tmp);
                // tmp = ut * (tmp x normal_face)
                tmp.scale(opp_point.ut);
                // sum += tmp*opp_point.ut
                f.gradu.add(tmp);

                curr = curr.next;
            }
            
        }
    }


    /** 
     * Compute the divergence of normalized gradients at the vertices
     */
    public void updateDivx() {
    	// TODO: 9 Update the divergence of the normalized gradients, ie., v.divX for each Vertex v
    	Vector3d gradu;
        for (Vertex curr_v : vertices){

            HalfEdge curr = curr_v.he;
            curr_v.divX = 0.0;
            do {
                Vector3d end1 = (Vector3d) curr.getVecToHead();
                end1.negate();
                double theta2 = angleWithNext(curr.prev());

                curr = curr.next;
                double theta1 = angleWithNext(curr);

                Vector3d end2 = (Vector3d) curr.getVecToHead();
                gradu = curr.leftFace.gradu;
                if (gradu.length() > 0)
                   gradu.normalize(); 

                 curr_v.divX += faceDivergence(gradu, end1, end2, theta1, theta2);
                curr = curr.twin;
            } while(curr != curr_v.he);

            curr_v.divX = curr_v.divX / 2;
        }
    }

        private double faceDivergence(Vector3d gradu, Vector3d end1, Vector3d end2, double theta1, double theta2){
            double cot_theta1 = 1/Math.tan(theta1);
            double term1 = cot_theta1 * end1.dot(gradu);

            double cot_theta2 = 1/Math.tan(theta2);
            double term2 = cot_theta2 * end2.dot(gradu);

            return term1 + term2;
        }


    /** Keep track of the maximum distance for debugging and colour map selection */
    double maxphi = 0 ;

    /**
     * Solves the distances
     * Uses Poisson equation, Laplacian of distance equal to divergence of normalized heat gradients.
     * This is step III in Algorithm 1 of the Geodesics in Heat paper, but here is done iteratively 
     * with a Gauss-Seidel solve of some number of steps to refine the solution whenever this method 
     * is called.
     * @param GSSteps number of Gauss-Seidel steps to take
     */
    public void solveDistanceStep( int GSSteps ) {		
    	for ( int i = 0; i < GSSteps; i++ ) {
    		for ( Vertex v : vertices ) {
    			// TODO: 10 Implement the inner loop of the Gauss-Seidel solve to compute the distances to each vertex, phi

                double RHS = v.divX;
                HalfEdge curr = v.he.next;
                int j = 0;
                double sum = 0;
                do {
                    //RHS += (t*v.LCij[j]) * curr.head.ut;
                    if (Double.isNaN(curr.head.phi)){
                        curr.head.phi = 0;
                    }
                    RHS -= v.LCij[j] * curr.head.phi;
                    j++;
                    curr = curr.twin.next;
                } while(curr != v.he.next);

                v.phi = RHS / v.LCii;
    		}    		
    	}

    	// Note that the solution to step III is unique only up to an additive constant,
    	// final values simply need to be shifted such that the smallest distance is zero. 
    	// We also identify the max phi value here to identify the maximum geodesic and to 
    	// use adjusting the colour map for rendering
    	double minphi = Double.MAX_VALUE;
    	maxphi = Double.MIN_VALUE;
		for ( Vertex v : vertices ) {
			if ( v.phi < minphi ) minphi = v.phi;
			if ( v.phi > maxphi ) maxphi = v.phi;
		}	
		maxphi -= minphi;
		for ( Vertex v : vertices ) {
			v.phi -= minphi;
		}
    }


    /**
     * Computes the cotangent Laplacian weights at each vertex.
	 * You can assume no boundaries and a triangular mesh! 
	 * You should store these weights in an array at each vertex,
	 * and likewise store the associated "vertex area", i.e., 1/3 of
	 * the surrounding triangles and NOT scale your Laplacian weights
	 * by the vertex area (see heat solve objective requirements).
     */
    public void computeLaplacian() {
    	for ( Vertex v : vertices ) {
    		// TODO: 6 Compute the Laplacian and store as vertex weights, and cotan operator diagonal LCii and off diagonal LCij terms.
    		v.area = 0;
    		v.LCii = 0;
    		v.LCij = new double[ v.valence() ];
    		double cotangent_sum = 0;
            double prev_angle = 0;
            double next_angle = 0;
            HalfEdge curr_he = v.he;
            int j_ind = 0;
            do {

                v.area += curr_he.leftFace.area;

                HalfEdge prev = curr_he.twin;

                // prev_angle = alpha
                prev_angle = angleWithNext(prev.next);
                // next_angle = beta
                next_angle = angleWithNext(curr_he.next);
                // compute area of face!! -- store area(i) = 1/3*sum(areas)
                
                curr_he = curr_he.next.twin;
                
                // compute sum of cotangents
				cotangent_sum = 1/Math.tan(prev_angle) + 1/Math.tan(next_angle);
                v.LCij[j_ind] = cotangent_sum / 2;
                v.LCii += v.LCij[j_ind];

                j_ind++;

            } while(curr_he != v.he);

            v.area = v.area / 3;
            v.LCii = - v.LCii;

    	}
        
    }

    /** 
     * Computes the angle between the provided half edge and the next half edge
     * @param he specify which half edge
     * @return the angle in radians
     */
    private double angleWithNext( HalfEdge he ) {
    	// TODO: 6 Implement this function to compute the angle with next edge... you'll want to use this in a few places
        
        if (he.angle_with_next != -360){
            return he.angle_with_next; 
        }

    	// v1 = from this head to next head
    	Vector3d v1 = (Vector3d) he.next.getVecToHead();
        Vector3d v2 = (Vector3d) he.getVecToHead();
    	// v2 = from this head to prev head
        v2.negate();
    	
        he.angle_with_next = v1.angle(v2);
    	return he.angle_with_next;
    }

    /**
     * Legacy drawing code for the half edge data structure by drawing each of its faces.
     * Legacy in that this code uses v1immediate mode OpenGL.  Per vertex normals are used
     * to draw the smooth surface if they are set in the vertices. 
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        for ( Face face : faces ) {
            HalfEdge he = face.he;
            if ( he.head.n == null ) { // don't have per vertex normals? use the face
                gl.glBegin( GL2.GL_POLYGON );
                Vector3d n = he.leftFace.n;
                gl.glNormal3d( n.x, n.y, n.z );
                HalfEdge e = he;
                do {
                	Point3d p = e.head.p;
                    gl.glVertex3d( p.x, p.y, p.z );
                    e = e.next;
                } while ( e != he );
                gl.glEnd();
            } else {
                gl.glBegin( GL2.GL_POLYGON );                
                HalfEdge e = he;
                do {
                	Point3d p = e.head.p;
                    Vector3d n = e.head.n;
                    gl.glNormal3d( n.x, n.y, n.z );
                    gl.glVertex3d( p.x, p.y, p.z );
                    e = e.next;
                } while ( e != he );
                gl.glEnd();
            }
        }
    }

}


