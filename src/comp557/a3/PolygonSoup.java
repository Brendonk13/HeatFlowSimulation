package comp557.a3;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import java.lang.Math;

/**
 * Simple implementation of a loader for a polygon soup
 */
public class PolygonSoup {


    // Brendon Keirle, 260685377

    /** List of vertex objects used in the mesh. */
    public ArrayList<Vertex> vertexList = new ArrayList<Vertex>();

    /** List of faces, where each face is a list of integer indices into the vertex list. */
    public ArrayList<int[]> faceList = new ArrayList<int[]>();

    /** Map for keeping track of how many n-gons we have for each n */
    // ie: map = {'number-of-4-vertex-faces' : 20, '# 3-vertex-faces' : 2}
    private TreeMap<Integer,Integer> faceSidesHistogram = new TreeMap<Integer,Integer>();

    /** A string summarizing the face sides histogram */
    public String soupStatistics;

    /**
     * Creates a polygon soup by loading an OBJ file
     * @param file
     */
    public PolygonSoup(String file) {
        try {
            FileInputStream fis = new FileInputStream( file );
            InputStreamReader isr = new InputStreamReader( fis );
            BufferedReader reader = new BufferedReader( isr );
            String line;
            while ((line = reader.readLine()) != null) {
                if ( line.startsWith("v ") ) {
                    parseVertex(line);
                } else if ( line.startsWith("f ") ) {
                    parseFace(line);
                } 
            }
            reader.close();
            isr.close();
            fis.close();

            soupStatistics = file + "\n" + "faces = " +faceList.size() + "\nverts = " + vertexList.size() + "\n";
            for ( Map.Entry<Integer,Integer> e : faceSidesHistogram.entrySet() ) {
                soupStatistics += e.getValue() + " ";
                if ( e.getKey() == 3 ) {
                    soupStatistics += "triangles\n";
                } else if ( e.getKey() == 4 ) {
                    soupStatistics += "quadrilaterals\n";
                } else {
                    soupStatistics += e.getKey() + "-gons\n";
                }
            }
            System.out.println( soupStatistics );

            // TODO: 1 compute a bounding box and scale and center the geometry
            makeBoundingBox();

         } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void makeBoundingBox(){

    	// dim_lengths is array of maximum displacements between points
        double[] dim_lengths = centerVertices();
        double max_len = -10;
        // int max_dim = -1;
        for (int i=0; i<3; i++){
            if (dim_lengths[i] > max_len) {
                max_len = dim_lengths[i];
                //max_dim = i;
            }
        }

        double scale = 10 / max_len;
        for (Vertex v : vertexList){
            v.p.x = v.p.x * scale;
            v.p.y = v.p.y * scale;
            v.p.z = v.p.z * scale;
        }
    }

    private double[] centerVertices(){

            double[] x_bound = new double[2];
            double[] y_bound = new double[2];
            double[] z_bound = new double[2];

            for (Vertex v : vertexList){
                checkNewMax(v.p.x, x_bound);
                checkNewMax(v.p.y, y_bound);
                checkNewMax(v.p.z, z_bound);
            }
            double[] middle_bounds = getMiddlePoint(x_bound, y_bound, z_bound);
            // transform vertices according 
            for (Vertex v : vertexList){
                v.p.x -= middle_bounds[0];
                v.p.y -= middle_bounds[1];
                v.p.z -= middle_bounds[2];
            }

            double[] dim_lengths = new double[] {
                x_bound[1] - x_bound[0],
                y_bound[1] - y_bound[0],
                z_bound[1] - z_bound[0]
            };
            return dim_lengths;
    }

    
    private void checkNewMax(double coord, double[] bound){
        // update bounds if new lower bound found
        if (coord < bound[0]){
            bound[0] = coord;
        }
        // update bounds if new upper bound found 
        if (coord > bound[1]){
            bound[1] = coord;
        }
    }

    private double[] getMiddlePoint(double[] x_bound, double[] y_bound, double[] z_bound){
        return new double[] {
            x_bound[0] + (x_bound[1] - x_bound[0]) / 2,
            y_bound[0] + (y_bound[1] - y_bound[0]) / 2,
            z_bound[0] + (z_bound[1] - z_bound[0]) / 2
        };
    }


    /**
     * Parses a vertex definition from a line in an obj file, and 
     * directly inserts it into the vertex list.
     * Assumes that there are three components.
     * @param newline
     * @return a new vertex object
     */
    private Vertex parseVertex(String newline) {        
        // Remove the tag "v "
        newline = newline.substring(2, newline.length());
        StringTokenizer st = new StringTokenizer(newline, " ");
        Vertex v = new Vertex();
        v.p.x = Double.parseDouble(st.nextToken());
        v.p.y = Double.parseDouble(st.nextToken());
        v.p.z = Double.parseDouble(st.nextToken());
        v.index = vertexList.size();
        vertexList.add( v );
        return v;
    }

    /**
     * Gets the list of indices for a face from a string in an obj file.
     * Simply ignores texture and normal information for simplicity
     * @param newline
     * @return list of indices
     */
    private int[] parseFace(String newline) {
        // Remove the tag "f "
        newline = newline.substring(2, newline.length());
        // vertex/texture/normal tuples are separated by spaces.
        StringTokenizer st = new StringTokenizer(newline, " ");
        int count = st.countTokens();
        int v[] = new int[count];
        for (int i = 0; i < count; i++) {
            // first token is vertex index... we'll ignore the rest (if it exists)
            StringTokenizer st2 = new StringTokenizer(st.nextToken(),"/");
            v[i] = Integer.parseInt(st2.nextToken()) - 1; // want zero indexed vertices!            
        }
        Integer n = faceSidesHistogram.get( count );
        if ( n == null ) {
            faceSidesHistogram.put( count, 1 );
        } else {
            faceSidesHistogram.put( count, n + 1 );
        }
        faceList.add( v );
        return v;
    }    

    /**
     * Draw the polygon soup using legacy immediate mode OpenGL
     * @param drawable
     */
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        // assume triangular faces!
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();
        Vector3d n = new Vector3d();
        for ( int[] faceVertex : faceList ) {
            Point3d p0 = vertexList.get( faceVertex[0] ).p;
            Point3d p1 = vertexList.get( faceVertex[1] ).p;
            Point3d p2 = vertexList.get( faceVertex[2] ).p;
            v1.sub( p1,p0 );
            v2.sub( p2,p1 );
            n.cross( v1, v2 );
            gl.glBegin( GL2.GL_POLYGON );
            gl.glNormal3d( n.x, n.y, n.z );
            for ( int i = 0; i < faceVertex.length; i++ ) {
                Point3d p = vertexList.get( faceVertex[i] ).p;
                gl.glVertex3d( p.x, p.y, p.z );
            }
            gl.glEnd();
        }        
    }
}
