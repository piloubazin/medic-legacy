package edu.jhmi.rad.medic.libraries;

import java.io.*;
import java.util.*;

import edu.jhmi.rad.medic.structures.BinaryTree;
import edu.jhmi.rad.medic.structures.CriticalPointLUT;
import edu.jhmi.rad.medic.utilities.Numerics;

/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectProcessing {
	
	// no data: used as a library of functions
	private final static int MaxObject = 1000000;
	
	public final static int SUPERIOR =  10;
	public final static int SUPEQUAL =  11;
	public final static int INFERIOR =  12;
	public final static int INFEQUAL =  13;
	public final static int EQUAL =  	14;
	public final static int UNEQUAL =  	15;
	public final static int NONE =      25;
	public final static int AND =       26;
	public final static int OR =        27;
	public final static int XOR =       28;
	
	public final static float SQR2 = (float)Math.sqrt(2.0f);
	public final static float SQR3 = (float)Math.sqrt(3.0f);
   // object manipulation
    
    public static final boolean[][][] objectFromImage(float[][][] img, int nx, int ny, int nz, float level1, float level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromImage(float[][][] img, int nx, int ny, int nz, float level, int type) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type==INFERIOR) && (img[x][y][z] >=level) ) 	obj[x][y][z] = false;
            if ( (type==INFEQUAL) && (img[x][y][z] > level) ) 	obj[x][y][z] = false;
            if ( (type==SUPERIOR) && (img[x][y][z] <=level) ) 	obj[x][y][z] = false;
            if ( (type==SUPEQUAL) && (img[x][y][z] < level) ) 	obj[x][y][z] = false;
			if ( (type==EQUAL) 	  && (img[x][y][z]!=level) ) 	obj[x][y][z] = false;
            if ( (type==UNEQUAL)  && (img[x][y][z]==level) ) 	obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[] objectFromImage(float[] img, int nx, int ny, int nz, float level, int type) {
        boolean[] obj = new boolean[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = true;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = false;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = false;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = false;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = false;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = false;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = false;
        }
        return obj;
    }
    public static final float[] floatObjectFromImage(float[] img, int nx, int ny, int nz, float level, int type) {
        float[] obj = new float[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = 1.0f;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = 0.0f;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = 0.0f;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = 0.0f;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = 0.0f;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = 0.0f;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = 0.0f;
        }
        return obj;
    }
    public static enum Comparator {SUP,SUP_EQ,INF,INF_EQ};
    public static final boolean[][][] objectFromImage(float[][][] img,float thresh,Comparator comp) {
        int x,y,z;
        int nx=img.length;
        int ny=img[0].length;
        int nz=img[0][0].length;
        boolean[][][] obj = new boolean[nx][ny][nz];
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	switch(comp){
        		case SUP:obj[x][y][z]=(img[x][y][z]>thresh);break;
        		case SUP_EQ:obj[x][y][z]=(img[x][y][z]>=thresh);break;
        		case INF:obj[x][y][z]=(img[x][y][z]<thresh);break;
        		case INF_EQ:obj[x][y][z]=(img[x][y][z]<=thresh);break;
        	}
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(int[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(byte[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(byte[][][] img, int nx, int ny, int nz, int level, int type) {
		return objectFromLabelImage(img, nx,ny,nz, level, level, type, type);
	}
    public static final boolean[][][] objectFromLabelImage(int[][][] img, int nx, int ny, int nz, int level, int type) {
		return objectFromLabelImage(img, nx,ny,nz, level, level, type, type);
	}
	public static final boolean[][] objectFromLabelImage(int[][] img, int nx, int ny, int level1, int level2, int type1, int type2) {
        int x,y;
        boolean[][] obj = new boolean[nx][ny];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            obj[x][y] = true;
            if ( (type1==INFERIOR) && (img[x][y] >=level1) ) obj[x][y] = false;
            if ( (type1==INFEQUAL) && (img[x][y] > level1) ) obj[x][y] = false;
            if ( (type1==SUPERIOR) && (img[x][y] <=level1) ) obj[x][y] = false;
            if ( (type1==SUPEQUAL) && (img[x][y] < level1) ) obj[x][y] = false;
            if ( (type1==EQUAL)    && (img[x][y]!= level1) ) obj[x][y] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y] <=level2) ) obj[x][y] = false;
            if ( (type2==SUPEQUAL) && (img[x][y] < level2) ) obj[x][y] = false;
            if ( (type2==INFERIOR) && (img[x][y] >=level2) ) obj[x][y] = false;
            if ( (type2==INFEQUAL) && (img[x][y] > level2) ) obj[x][y] = false;
            if ( (type2==EQUAL)    && (img[x][y]!= level2) ) obj[x][y] = false;
        }
        return obj;
    }
	public static final byte[][][] objectFromLabelImageToByte(byte[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = 0;
        }
        return obj;
    }
	public static final byte[][][] objectFromLabelImageToByte(byte[][][] img, int nx, int ny, int nz, int[] maskinds) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = 1;
            if (maskinds[img[x][y][z]]<0) obj[x][y][z] = 0;
        }
        return obj;
    }
	public static final byte[][][] objectFromLabelImageToByteZwise(byte[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (z=0;z<nz;z++) for (y=0;y<ny;y++) for (x=0;x<nx;x++) {
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = 0;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        int xyz =0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[xyz] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[xyz] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[xyz] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[xyz] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[xyz]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[xyz] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[xyz] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[xyz] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[xyz] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[xyz]!= level2) ) obj[x][y][z] = 0;
            
            xyz++;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[] img, int nx, int ny, int nz, int[] maskinds) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        
        int xyz =0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	
            if (img[xyz]>=0 && maskinds[img[xyz]]>=0) obj[x][y][z] = 1;
            
            xyz++;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByteZwise(byte[] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        int xyz =0;
        for (z=0;z<nz;z++) for(y=0;y<ny;y++) for (x=0;x<nx;x++)  {
        	
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[xyz] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[xyz] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[xyz] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[xyz] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[xyz]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[xyz] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[xyz] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[xyz] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[xyz] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[xyz]!= level2) ) obj[x][y][z] = 0;
            
            xyz++;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByte(img,nx,ny,nz,level,level,type,type);
    }
    public static final byte[][][] objectFromLabelImageToByteZwise(byte[] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByteZwise(img,nx,ny,nz,level,level,type,type);
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[][][] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByte(img,nx,ny,nz,level,level,type,type);
    }
    public static final byte[][][] objectFromLabelImageToByteZwise(byte[][][] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByteZwise(img,nx,ny,nz,level,level,type,type);
    }
    public static final boolean[][][] objectFromLabelImage(int[][][] img, int nx, int ny, int nz, int[] groups) {
		int x,y,z;
		boolean[][][] obj = new boolean[nx][ny][nz];
		Arrays.sort(groups);
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			if(Arrays.binarySearch(groups, img[x][y][z])>=0){
				obj[x][y][z]=true;
			}
		}
		return obj;
	}
    public static final float[][][] labelFromImage(float[][][] img, int nx, int ny, int nz, float level1, float level2, int type1, int type2) {
		int x,y,z;
		float[][][] label = new float[nx][ny][nz];
		boolean[][][] obj = objectFromImage(img, nx, ny, nz, level1, level2, type1, type2);
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++)
			if (obj[x][y][z]) label[x][y][z] = 1.0f;
			else label[x][y][z] = 0.0f;
		obj = null;
		return label;
	}
    public static final boolean[][] objectFromImageXSlice(float[][][] img, int nx, int ny, int nz, int x, float level1, float level2, int type1, int type2) {
        int y,z;
        boolean[][] obj = new boolean[ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[y][z] = false;
        }
        return obj;
    }
    
    public static final boolean[][] objectFromImageYSlice(float[][][] img, int nx, int ny, int nz, int y, float level1, float level2, int type1, int type2) {
        int x,z;
        boolean[][] obj = new boolean[nx][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
            obj[x][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][z] = false;
        }
        return obj;
    }
    
    public static final boolean[][] objectFromImageZSlice(float[][][] img, int nx, int ny, int nz, int z, float level1, float level2, int type1, int andor, int type2) {
        int x,y;
        boolean[][] obj = new boolean[nx][ny];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            if (andor==AND) {
                // and: if either of the properties is false, then no object
                obj[x][y] = true;
                if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y] = false;
                if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y] = false;
                if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y] = false;
                if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y] = false;

                if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y] = false;
                if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y] = false;
                if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y] = false;
                if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y] = false;
            } else {
                // or: if either of the properties is true, then object
                obj[x][y] = false;
                if ( (type1==INFERIOR) && (img[x][y][z] < level1) ) obj[x][y] = true;
                if ( (type1==INFEQUAL) && (img[x][y][z] <=level1) ) obj[x][y] = true;
                if ( (type1==SUPERIOR) && (img[x][y][z] > level1) ) obj[x][y] = true;
                if ( (type1==SUPEQUAL) && (img[x][y][z] >=level1) ) obj[x][y] = true;

                if ( (type2==SUPERIOR) && (img[x][y][z] > level2) ) obj[x][y] = true;
                if ( (type2==SUPEQUAL) && (img[x][y][z] >=level2) ) obj[x][y] = true;
                if ( (type2==INFERIOR) && (img[x][y][z] < level2) ) obj[x][y] = true;
                if ( (type2==INFEQUAL) && (img[x][y][z] <=level2) ) obj[x][y] = true;
            }
        }
        return obj;
    }
	public static final boolean[] objectFromLabelImage(byte[] img, int nx, int ny, int nz, byte level, int type) {
        boolean[] obj = new boolean[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = true;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = false;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = false;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = false;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = false;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = false;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = false;
        }
        return obj;
    }

    public static final int countLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[x][y][z]);
                Nlb++;
            }
        }
        return Nlb;
    }
    public static final int countLabels(int[] label, int nx, int ny) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        Nlb = 0;
        for (xyz=0;xyz<nx*ny;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    public static final int countLabels(int[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(float[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
			int xyz = x+nx*y+nx*ny*z;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(int[][] label, int nx, int ny) {
        int x,y,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
		Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int[] listLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    
    public static final int[] listLabels(int[][] label, int nx, int ny) {
        int x,y,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }

    public static final float[] listLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }

    public static final float[] listLabels(float[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final int[] listLabels(int[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        int ind = -1;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	ind = x + y*nx + z*nx*ny;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[ind]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[ind]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final int[] listLabels(byte[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        int ind = -1;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	ind = x + y*nx + z*nx*ny;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[ind]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,(int)label[ind]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final float[] listOrderedLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final int[] listOrderedLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] listOrderedLabels(int[][] label, int nx, int ny) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++)  {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final int[] listOrderedLabels(int[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] listOrderedLabels(int[] label, int nx, int ny) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final short[] listOrderedLabels(short[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Short> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(short)0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((short)Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		short[] list = new short[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final short[] listOrderedLabels(short[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Short> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(short)0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((short)Nlb,label[xyz]);
                Nlb++;
            }
        }
		short[] list = new short[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final byte[] listOrderedLabels(byte[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Byte> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(byte)0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((byte)Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		byte[] list = new byte[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final boolean[][][] largestObjectFromLabel(int[][][] label, int nlb, int nx, int ny, int nz) {
        int x,y,z,n;
        boolean[][][] obj = new boolean[nx][ny][nz];
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (n=0;n<nlb;n++) Nobj[n]=0;
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            Nobj[ label[x][y][z] ]++;
        }
        size=0;best=0;
        for (n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
				obj[x][y][z] = (label[x][y][z]==best);
			}
		else
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
				obj[x][y][z] = false;
			}
        return obj;
    }
    
   public static final boolean[] largestObjectFromLabel(int[] label, int nlb, int nx, int ny, int nz) {
        boolean[] obj = new boolean[nx*ny*nz];
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (int n=0;n<nlb;n++) Nobj[n]=0;
        
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
            Nobj[ label[xyz] ]++;
        }
        size=0;best=0;
        for (int n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				obj[xyz] = (label[xyz]==best);
			}
		else
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				obj[xyz] = false;
			}
        return obj;
    }
    
    public static final boolean[][][] largestObjectFromLabel(int[][][] label, int nx, int ny, int nz) {
    	int nlabels = countLabels(label,nx,ny,nz);
    	int[] lablist = listLabels(label,nx,ny,nz);
    	
    	if(nlabels>0){
    		return largestObjectFromLabel(label, nlabels, nx,ny,nz);
    	}else{
    		return null;
    	}
    }
    
    public static final boolean[][] largestObjectFromLabel(int[][] label, int nlb, int nx, int ny) {
        int x,y,n;
        boolean[][] obj = new boolean[nx][ny];
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (n=0;n<nlb;n++) Nobj[n]=0;
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            Nobj[ label[x][y] ]++;
        }
        size=0;best=0;
        for (n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
				obj[x][y] = (label[x][y]==best);
			}
		else
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
				obj[x][y] = false;
			}
        return obj;
    }
    
    public static final boolean[][] largestObjectFromLabel(int[][] label, int nx, int ny) {
		return largestObjectFromLabel(label, countLabels(label,nx,ny), nx,ny);
    }
    
    
    public static final boolean[][][] connectedObjectAt(boolean[][][] object, int x0, int y0, int z0, int nx, int ny, int nz) {
        boolean[][][] obj = new boolean[nx][ny][nz];
        
    	int[][][] label = ObjectProcessing.connected6Object3D(object, nx, ny, nz);
		 
        int lb = label[x0][y0][z0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			obj[x][y][z] = (label[x][y][z]==lb);
		}
		return obj;
    }
    
    /*
     * @brief Fill holes that are not 'conn' connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     * @param conn Connectivity to use.
     */
    public static final boolean[][][] removeHoles(boolean[][][] object, int nx, int ny, int nz, int conn) {		
		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					object[x][y][z] = !object[x][y][z];
				}
		
    	int[][][] lb;
		
		if (conn == 6) {
			lb = ObjectProcessing.connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = ObjectProcessing.connected18Object3D(object, nx, ny, nz);
		} else if (conn == 26) {
			lb = ObjectProcessing.connected26Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		
		object = ObjectProcessing.largestObjectFromLabel(lb, nx, ny, nz);
		
		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					object[x][y][z] = !object[x][y][z];
				}
		
		return object;
    }
    
    
    /*
     * @brief Fill holes that are not 6 connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] removeHoles6(boolean[][][] object, int nx, int ny, int nz) {
    	return removeHoles(object, nx, ny, nz, 6);
    }
    
    
    /*
     * @brief Fill holes that are not 18 connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] removeHoles18(boolean[][][] object, int nx, int ny, int nz) {
    	return removeHoles(object, nx, ny, nz, 18);
    }
    
    
    /*
     * @brief Returns the largest 'conn' connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     * @param conn Connectivity to use.
     */
    public static final boolean[][][] largestObject(boolean[][][] object, int nx, int ny, int nz, int conn) {		
    	int[][][] lb;
		
		if (conn == 6) {
			lb = ObjectProcessing.connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = ObjectProcessing.connected18Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		
		object = ObjectProcessing.largestObjectFromLabel(lb, nx, ny, nz);
		
		return object;
    }
    
    
    /*
     * @brief Returns the largest 6 connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] largest6Object(boolean[][][] object, int nx, int ny, int nz) {		
    	return ObjectProcessing.largestObject(object, nx, ny, nz, 6);
    }
    
    
    /*
     * @brief Returns the largest 18 connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] largest18Object(boolean[][][] object, int nx, int ny, int nz) {		
    	return ObjectProcessing.largestObject(object, nx, ny, nz, 18);
    }
    
    
    
	// complete image connectivity computations
	
	/** 
	 *	Connected components of an object.
     *  2D images: 4-connectivity
	 */
	public static final int[][] connected4Object2D(boolean img[][], int nx, int ny) {
		int Nlabel;
		int[][] label = new int[nx][ny];
		ArrayList<Integer>   lb = new ArrayList();
		int lbMin;
		int x,y,c,i,j,k;
		int Nlb;
		int[]   connect = new int[4];
		int Nconnect;
		int AddLabel=0;
	
		// the input is a 3x3 binary image (0 out, 1 in)
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = 0;
			}
		}
		
		lb.add(0,0);
		Nlabel = 1;
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				if (img[x][y]) {
					// object point: neighbors ?
					Nconnect = 0;
                    for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) {
                        if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) ) {
                            if (i*i+j*j < 2) {
                                if (label[x+i][y+j] > 0) {
                                    connect[Nconnect] = lb.get( label[x+i][y+j] );
                                    Nconnect++;
                                }
                            }
                        }
                    }
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						lbMin = lb.get(connect[0]);
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb.get(connect[k]) );
						for (k=0;k<Nconnect;k++) lb.set(connect[k],lbMin);
						label[x][y] = lbMin;
					} else {
						// new, unconnected region
						label[x][y] = Nlabel;
						lb.add(Nlabel,Nlabel);
						Nlabel++;
						/* not needed with ArrayLists
						// check if the number of labels is above the threshold
						if (Nlabel>=lb.length-1) {
							int[] tmp = new int[2*lb.length];
							for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
							lb = tmp;
						}
						*/
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb.get(c)!=c) c = lb.get(c);
			lb.set(k, c);
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb.get(k)==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label image
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = lb2[ lb.get( label[x][y] ) ];
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	
	/** 
	 *	Connected components of an object.
     *  2D images: 8-neighborhood 
	 */
	public static final int[][] connected8Object2D(boolean img[][], int nx, int ny) {
		int Nlabel;
		int[][] label = new int[nx][ny];
		int[]   lb = new int[MaxObject];
		int lbMin;
		int x,y,c,i,j,k;
		int Nlb;
		int[]   connect = new int[4];
		int Nconnect;
		int AddLabel=0;
	
		// the input is a 3x3 binary image (0 out, 1 in)
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = 0;
			}
		}
		
		lb[0] = 0;
		Nlabel = 1;
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				if (img[x][y]) {
					// object point: neighbors ?
					// object point: neighbors ?
					Nconnect = 0;
                    for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) {
                        if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) ) {
                            if (i*i+j*j < 3) {
                                if (label[x+i][y+j] > 0) {
                                    connect[Nconnect] = lb[ label[x+i][y+j] ];
                                    Nconnect++;
                                }
                            }
                        }
                    }
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						lbMin = lb[connect[0]];
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb[connect[k]]);
						for (k=0;k<Nconnect;k++) lb[connect[k]] = lbMin;
						label[x][y] = lbMin;
					} else {
						// new, unconnected region
						label[x][y] = Nlabel;
						lb[Nlabel] = Nlabel;
						Nlabel++;
						// check if the number of labels is above the threshold
						if (Nlabel>=lb.length-1) {
							int[] tmp = new int[2*lb.length];
							for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
							lb = tmp;
						}
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label image
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = lb2[ lb[ label[x][y] ] ];
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}

	/** 
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[][][] connected6Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                    if (label[x+i][y+j][z+k] > 0) {
                                        connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[x][y][z] = lbMin;
						} else {
							// new, unconnected region
							label[x][y][z] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
											lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	/** 
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[] connected6Object3D(boolean img[], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[]   label = new int[nx*ny*nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (img[xyz]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                    if (label[xyz+i+nx*j+nx*ny*k] > 0) {
                                        connect[Nconnect] = lb[ label[xyz+i+nx*j+nx*ny*k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[xyz] = lbMin;
						} else {
							// new, unconnected region
							label[xyz] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (lb2[ lb[ label[xyz] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>0) {
										if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>lb2[ lb[ label[xyz] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ];
											lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ] = lb2[ lb[ label[xyz] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ label[xyz] ] ];
        }    
		
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	/**
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
     */
    public static final int[][][] connected6Object3D(boolean img[][][]) {
    	return connected6Object3D(img, img.length,img[0].length,img[0][0].length);
    }
    /**
	 *	Connected components of an object.
     *  3D images: 26-neighborhood
     */
    public static final int[][][] connected26Object3D(boolean img[][][]) {
    	return connected26Object3D(img, img.length,img[0].length,img[0][0].length);
    }
	/**
	 *	Connected components of an object.
     *  3D images: 18-neighborhood
     */
    public static final int[][][] connected18Object3D(boolean img[][][]) {
    	return connected18Object3D(img, img.length,img[0].length,img[0][0].length);
    }
	
    /**
	 *	Connected components of an object.
     *  3D images: 18-neighborhood
     */
    public static final int[][][] connected18Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[18];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 3) {
                                    if (label[x+i][y+j][z+k] > 0) {
                                        connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[x][y][z] = lbMin;
                        } else {
                            // new, unconnected region
                            label[x][y][z] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 3) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
											lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }

    /**
     *	Connected components of an object.
     *  3D images: 26-neighborhood
     */
    public static final int[][][] connected26Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (label[x+i][y+j][z+k] > 0) {
                                    connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                    Nconnect++;
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[x][y][z] = lbMin;
                        } else {
                            // new, unconnected region
                            label[x][y][z] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
										//System.out.print("labelling problem!");
										int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
										lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
										// if there are higher labels, decrease them
										for (int n=1;n<Nlabel;n++) {
											if (lb2[n] > badLb) lb2[n]--;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }   
	
	/**
     *  Counts the different boundaries
	 *	and returns the Euler characteristic
     */
    public static final int eulerCharacteristic(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int Nv = 0;
		int Ne = 0;
		int Nf = 0;
		
		// hyp: the object is away from the image boundary
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// count all boundaries : 6-C
			if (img[x][y][z]!=img[x-1][y][z]) Nf++;
			if (img[x][y][z]!=img[x][y-1][z]) Nf++;
			if (img[x][y][z]!=img[x][y][z-1]) Nf++;
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// edges : XY
			if (img[x][y][z]!=img[x-1][y-1][z]) Ne++;
			else if (img[x-1][y][z]!=img[x][y-1][z]) Ne++;
			else if ( (img[x][y][z]==img[x-1][y-1][z]) && (img[x][y][z]!=img[x][y-1][z]) && (img[x][y][z]!=img[x-1][y][z]) ) {
				Ne+=2;
			} else if ( (img[x-1][y][z]==img[x][y-1][z]) && (img[x-1][y][z]!=img[x][y][z]) && (img[x-1][y][z]!=img[x-1][y-1][z]) ) {
				Ne+=2;
			}
			// edges : YZ
			if (img[x][y][z]!=img[x][y-1][z-1]) Ne++;
			else if (img[x][y-1][z]!=img[x][y][z-1]) Ne++;
			else if ( (img[x][y][z]==img[x][y-1][z-1]) && (img[x][y][z]!=img[x][y-1][z]) && (img[x][y][z]!=img[x][y][z-1]) ) {
				Ne+=2;
			} else if ( (img[x][y-1][z]==img[x][y][z-1]) && (img[x][y-1][z]!=img[x][y][z]) && (img[x][y-1][z]!=img[x][y-1][z-1]) ) {
				Ne+=2;
			}
			// edges : ZX
			if (img[x][y][z]!=img[x-1][y][z-1]) Ne++;
			else if (img[x][y][z-1]!=img[x-1][y][z]) Ne++;
			else if ( (img[x][y][z]==img[x-1][y][z-1]) && (img[x][y][z]!=img[x][y][z-1]) && (img[x][y][z]!=img[x-1][y][z]) ) {
				Ne+=2;
			} else if ( (img[x][y][z-1]==img[x-1][y][z]) && (img[x][y][z-1]!=img[x][y][z]) && (img[x][y][z-1]!=img[x-1][y][z-1]) ) {
				Ne+=2;
			}
		}
		int Nobj,Nbg,N6,N18,N26;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// classify by cases : look into 2x2 neighborhood, count elements and connections
			Nobj=0;Nbg=0;N6=0;N18=0;N26=0;
			if (img[x][y][z]) Nobj++; else Nbg++;
			if (img[x-1][y][z]) Nobj++; else Nbg++;
			if (img[x][y-1][z]) Nobj++; else Nbg++;
			if (img[x][y][z-1]) Nobj++; else Nbg++;
			if (img[x][y-1][z-1]) Nobj++; else Nbg++;
			if (img[x-1][y][z-1]) Nobj++; else Nbg++;
			if (img[x-1][y-1][z]) Nobj++; else Nbg++;
			if (img[x-1][y-1][z-1]) Nobj++; else Nbg++;
			
			// use the smallest one
			if (Nobj<=Nbg) {
				if (Nobj==0) { 
					// do nothing 
				} else if (Nobj==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (img[x][y][z]) && (img[x-1][y][z]) ) N6++;
					if ( (img[x][y][z]) && (img[x][y-1][z]) ) N6++;
					if ( (img[x][y][z]) && (img[x][y][z-1]) ) N6++;
					if ( (img[x-1][y][z]) && (img[x-1][y-1][z]) ) N6++;
					if ( (img[x-1][y][z]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x][y-1][z]) && (img[x-1][y-1][z]) ) N6++;
					if ( (img[x][y-1][z]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x][y][z-1]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x][y][z-1]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x-1][y-1][z]) ) N6++;
					//System.out.print("case: obj "+Nobj+", ("+N6+" |");
						
					if ( (Nobj==2) && (N6==1) ) Nv++;
					else if ( (Nobj==3) && (N6==2) ) Nv++;
					else if ( (Nobj==4) && (N6==4) ) Nv++;
					else if ( (Nobj==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (img[x][y][z]) && (img[x][y-1][z-1]) ) N18++;
						if ( (img[x][y][z]) && (img[x-1][y][z-1]) ) N18++;
						if ( (img[x][y][z]) && (img[x-1][y-1][z]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x][y-1][z]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x][y][z-1]) ) N18++;
						if ( (img[x][y-1][z]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x][y-1][z]) && (img[x][y][z-1]) ) N18++;
						if ( (img[x][y][z-1]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x][y-1][z-1]) && (img[x-1][y][z-1]) ) N18++;
						if ( (img[x-1][y][z-1]) && (img[x-1][y-1][z]) ) N18++;
						if ( (img[x-1][y-1][z]) && (img[x][y-1][z-1]) ) N18++;
						//System.out.println(" "+N18+")");
						
						if ( (Nobj==2) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==3) ) {
							if (cObj==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nobj==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nobj==2) {
								if (cObj==26) Nv+=0;
								else Nv+=2;
							} else {
								System.out.println("!:"+Nobj+", "+N6+", "+N18);
							}
						}
					}
				}
			} else {
				if (Nbg==0) { 
					// do nothing 
				} else if (Nbg==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (!img[x][y][z]) && (!img[x-1][y][z]) ) N6++;
					if ( (!img[x][y][z]) && (!img[x][y-1][z]) ) N6++;
					if ( (!img[x][y][z]) && (!img[x][y][z-1]) ) N6++;
					if ( (!img[x-1][y][z]) && (!img[x-1][y-1][z]) ) N6++;
					if ( (!img[x-1][y][z]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x][y-1][z]) && (!img[x-1][y-1][z]) ) N6++;
					if ( (!img[x][y-1][z]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x][y][z-1]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x][y][z-1]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x-1][y-1][z]) ) N6++;
					
					if ( (Nbg==2) && (N6==1) ) Nv++;
					else if ( (Nbg==3) && (N6==2) ) Nv++;
					else if ( (Nbg==4) && (N6==4) ) Nv++;
					else if ( (Nbg==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (!img[x][y][z]) && (!img[x][y-1][z-1]) ) N18++;
						if ( (!img[x][y][z]) && (!img[x-1][y][z-1]) ) N18++;
						if ( (!img[x][y][z]) && (!img[x-1][y-1][z]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x][y-1][z]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x][y][z-1]) ) N18++;
						if ( (!img[x][y-1][z]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x][y-1][z]) && (!img[x][y][z-1]) ) N18++;
						if ( (!img[x][y][z-1]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x][y-1][z-1]) && (!img[x-1][y][z-1]) ) N18++;
						if ( (!img[x-1][y][z-1]) && (!img[x-1][y-1][z]) ) N18++;
						if ( (!img[x-1][y-1][z]) && (!img[x][y-1][z-1]) ) N18++;
						
						if ( (Nbg==2) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==3) ) {
							if (cBg==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nbg==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nbg==2) {
								if (cBg==26) Nv+=0;
								else Nv+=2;
							} else {
								System.out.println("!:"+Nbg+", "+N6+", "+N18);
							}
						}
					}
				}
			}

		}
		
//		System.out.println("Surface: "+Nv+" vertices, "+Ne+" edges, "+Nf+" faces");
			
		return Nv-Ne+Nf;
	}

	/**
     *  compute the thickness of the object at each point
	 *	(not properly functional yet)
     */
    public static final float[][][] thicknessMap(boolean img[][][], int nx, int ny, int nz) {
		float [][][] map = new float[nx][ny][nz];
		float thickness;
		int n,m;
		float min;
		
		// hyp: the object is away from the image boundary
		min = nx+ny+nz;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]) {
				// check for pointwise and edgewise connections : maybe not necessary
				// (at least, it picks-up too many cases)
				/*
				// corners: 8
				if (   (img[x][y][z]==img[x+1][y+1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x+1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x+1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y+1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x-1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x-1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y-1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x+1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x+1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y+1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x+1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x+1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else  if (   (img[x][y][z]==img[x-1][y-1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x-1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x-1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y-1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x+1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x+1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y+1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x-1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x-1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y-1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x-1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x-1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
				else
				// edges: 12
				// X,Y
				if (   (img[x][y][z]==img[x+1][y+1][z]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x+1][y-1][z]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y+1][z]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y-1][z]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else
	
				// Y,Z
				if (   (img[x][y][z]==img[x][y+1][z+1]) 
					&& (img[x][y][z]!=img[x][y][z+1])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y-1][z+1]) 
					&& (img[x][y][z]!=img[x][y][z+1])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y+1][z-1]) 
					&& (img[x][y][z]!=img[x][y][z-1])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y-1][z-1]) 
					&& (img[x][y][z]!=img[x][y][z-1])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else
					
				// Z,X	
				if (   (img[x][y][z]==img[x+1][y][z+1]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y][z+1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x+1][y][z-1]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y][z-1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y][z+1]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y][z+1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y][z-1]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y][z-1]) ) map[x][y][z] = 0.0f;
				else {
					*/
					// no 0 thick junction: find depth along all directions
					thickness = 0.0f; map[x][y][z] = nx+ny+nz;
					// compute thickness along the x,y,z directions
					n=1;m=1;
					while ( (x+n<nx) && (img[x][y][z]==img[x+n][y][z]) ) n++;
					while ( (x-m>=0) && (img[x][y][z]==img[x-m][y][z]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (y+n<ny) && (img[x][y][z]==img[x][y+n][z]) ) n++;
					while ( (y-m>=0) && (img[x][y][z]==img[x][y-m][z]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (img[x][y][z]==img[x][y][z+n]) ) n++;
					while ( (z-m>=0) && (img[x][y][z]==img[x][y][z-m]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					
					// compute along the edges
					/* not needed ? 
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (img[x][y][z]==img[x+n][y+n][z]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (img[x][y][z]==img[x-m][y-m][z]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y-n>=0) && (img[x][y][z]==img[x+n][y-n][z]) ) n++;
					while ( (x-m>=0) && (y+m<ny) && (img[x][y][z]==img[x-m][y+m][z]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (z+n<nz) && (img[x][y][z]==img[x+n][y][z+n]) ) n++;
					while ( (x-m>=0) && (z-m>=0) && (img[x][y][z]==img[x-m][y][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (z-n>=0) && (img[x][y][z]==img[x+n][y][z-n]) ) n++;
					while ( (x-m>=0) && (z+m<nz) && (img[x][y][z]==img[x-m][y][z+m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (y+n<ny) && (img[x][y][z]==img[x][y+n][z+n]) ) n++;
					while ( (z-m>=0) && (y-m>=0) && (img[x][y][z]==img[x][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (y-n>=0) && (img[x][y][z]==img[x][y-n][z+n]) ) n++;
					while ( (z-m>=0) && (y+m<ny) && (img[x][y][z]==img[x][y+m][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					
					// along corners
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (z+n<nz) && (img[x][y][z]==img[x+n][y+n][z+n]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (z-m>=0) && (img[x][y][z]==img[x-m][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (z-n>=0) && (img[x][y][z]==img[x+n][y+n][z-n]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (z+m<nz) && (img[x][y][z]==img[x-m][y-m][z+m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y-n>=0) && (z+n<nz) && (img[x][y][z]==img[x+n][y-n][z+n]) ) n++;
					while ( (x-m>=0) && (y+m<ny) && (z-m>=0) && (img[x][y][z]==img[x-m][y+m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x-n>=0) && (y+n<ny) && (z+n<nz) && (img[x][y][z]==img[x-n][y+n][z+n]) ) n++;
					while ( (x+m<nx) && (y-m>=0) && (z-m>=0) && (img[x][y][z]==img[x+m][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					*/
					
				//}
				if (map[x][y][z] < min) min = map[x][y][z];
			} else {
				map[x][y][z] = 0.0f;
			}
		}
		System.out.println("minimum thickness: "+min);
		return map;
	}
	
	/**
     *  checks that the object has all points being 6 connected to each neighbors
	 *	(not properly functional yet)
     */
    public static final float[][][] fullConnection(boolean img[][][], int nx, int ny, int nz) {
		float [][][] res = new float[nx][ny][nz];
		int[] count = new int[8];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			res[x][y][z] = 0;
			if (img[x][y][z]) {
				res[x][y][z]++;
				for (int n=0;n<8;n++) count[n] = 0;
				// count the number of 8-voxel cubes
				for (int i=0;i<=1;i++) for (int j=0;j<=1;j++) for (int l=0;l<=1;l++) {
					if (img[x+i][y+j][z+l]) count[0]++;
					if (img[x-i][y+j][z+l]) count[1]++;
					if (img[x+i][y-j][z+l]) count[2]++;
					if (img[x+i][y+j][z-l]) count[3]++;
					if (img[x+i][y-j][z-l]) count[4]++;
					if (img[x-i][y+j][z-l]) count[5]++;
					if (img[x-i][y-j][z+l]) count[6]++;
					if (img[x-i][y-j][z-l]) count[7]++;
				}
				for (int n=0;n<8;n++) if (count[n]==8) res[x][y][z]++;
			}
		}
		return res;
	}
	
	/**
     *  compute the highest connectivity junction for each pixel.
	 *	<p>
	 *	The result is 6 (only 6-C neighbors), 18 (only 6 and 18-C neighbors)
	 *	or 26 (6, 18 and 26-C neighbors).
     */
    public static final float[][][] connectivity(boolean img[][][], int nx, int ny, int nz) {
		float [][][] obj = new float[nx][ny][nz];
		int[] Nb = new int[8];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]) {
				obj[x][y][z] = 6;
				if ( img[x+1][y+1][z] && !(img[x+1][y][z] || img[x][y+1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y+1][z+1] && !(img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 18;
				if ( img[x+1][y][z+1] && !(img[x][y][z+1] || img[x+1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x-1][y+1][z] && !(img[x-1][y][z] || img[x][y+1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y-1][z+1] && !(img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 18;
				if ( img[x+1][y][z-1] && !(img[x][y][z-1] || img[x+1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x+1][y-1][z] && !(img[x+1][y][z] || img[x][y-1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y+1][z-1] && !(img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 18;
				if ( img[x-1][y][z+1] && !(img[x][y][z+1] || img[x-1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x-1][y-1][z] && !(img[x-1][y][z] || img[x][y-1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y-1][z-1] && !(img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 18;
				if ( img[x-1][y][z-1] && !(img[x][y][z-1] || img[x-1][y][z]) ) obj[x][y][z] = 18;
				
				// note: for 26-C, test only the case with no neighbors 
				// (otherwise, one of the neighbors show up as 18-C)
				if ( img[x+1][y+1][z+1] && !(img[x+1][y][z] || img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y+1][z+1] && !(img[x-1][y][z] || img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y-1][z+1] && !(img[x+1][y][z] || img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y+1][z-1] && !(img[x+1][y][z] || img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y-1][z-1] && !(img[x+1][y][z] || img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y+1][z-1] && !(img[x-1][y][z] || img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y-1][z+1] && !(img[x-1][y][z] || img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y-1][z-1] && !(img[x-1][y][z] || img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;

			} else {
				obj[x][y][z] = 0;
			}
		}
		return obj;
	}
	
	
	/**
	 *  find well-composed and not well-composed regions: groups imgects with relations
	 */
    public static final float[][][] wellComposed(boolean[][][] img, int nx, int ny, int nz) {
		float[][][] obj = new float[nx][ny][nz];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) obj[x][y][z] = 0.0f;
			else if (isWellComposed(img,x,y,z)) obj[x][y][z] = 1.0f;
			else obj[x][y][z] = 2.0f;
		}
		return obj;
	}
				
    public static final boolean isWellComposed(boolean[][] obj2d, int xi, int yi) {
		boolean[][][] obj = new boolean[3][3][3];
		
		for (int x=0;x<3;x++) for (int y=0;y<3;y++) for (int z=0;z<3;z++) {
			obj[x][y][z] = obj2d[x][y];	
		}
				
		return ObjectProcessing.isWellComposed(obj,xi,yi,1);
	}
				
	public static final boolean isWellComposed(boolean[][][] img, int x, int y, int z) {
		if (img[x][y][z]) {
			// 18-C
			if (img[x-1][y-1][z] && !(img[x-1][y][z] || img[x][y-1][z]) ) return false;
			if (img[x][y-1][z-1] && !(img[x][y-1][z] || img[x][y][z-1]) ) return false;
			if (img[x-1][y][z-1] && !(img[x][y][z-1] || img[x-1][y][z]) ) return false;
		
			if (img[x-1][y+1][z] && !(img[x-1][y][z] || img[x][y+1][z]) ) return false;
			if (img[x][y-1][z+1] && !(img[x][y-1][z] || img[x][y][z+1]) ) return false;
			if (img[x+1][y][z-1] && !(img[x][y][z-1] || img[x+1][y][z]) ) return false;
		
			if (img[x+1][y-1][z] && !(img[x+1][y][z] || img[x][y-1][z]) ) return false;
			if (img[x][y+1][z-1] && !(img[x][y+1][z] || img[x][y][z-1]) ) return false;
			if (img[x-1][y][z+1] && !(img[x][y][z+1] || img[x-1][y][z]) ) return false;
				
			if (img[x+1][y+1][z] && !(img[x+1][y][z] || img[x][y+1][z]) ) return false;
			if (img[x][y+1][z+1] && !(img[x][y+1][z] || img[x][y][z+1]) ) return false;
			if (img[x+1][y][z+1] && !(img[x][y][z+1] || img[x+1][y][z]) ) return false;
		
			// 26-C
			if (img[x-1][y-1][z-1] && !( (img[x-1][y][z] && img[x-1][y-1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z-1]) 
									  || (img[x][y-1][z] && img[x][y-1][z-1]) 
									  || (img[x][y-1][z] && img[x-1][y-1][z]) 
									  || (img[x][y][z-1] && img[x-1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y-1][z-1]) ) ) return false;
		
			if (img[x-1][y-1][z+1] && !( (img[x-1][y][z] && img[x-1][y-1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z+1]) 
									  || (img[x][y-1][z] && img[x][y-1][z+1]) 
									  || (img[x][y-1][z] && img[x-1][y-1][z]) 
									  || (img[x][y][z+1] && img[x-1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y-1][z+1]) ) ) return false;
		
			if (img[x-1][y+1][z-1] && !( (img[x-1][y][z] && img[x-1][y+1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z-1]) 
									  || (img[x][y+1][z] && img[x][y+1][z-1]) 
									  || (img[x][y+1][z] && img[x-1][y+1][z]) 
									  || (img[x][y][z-1] && img[x-1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y+1][z-1]) ) ) return false;
		
			if (img[x+1][y-1][z-1] && !( (img[x+1][y][z] && img[x+1][y-1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z-1]) 
									  || (img[x][y-1][z] && img[x][y-1][z-1]) 
									  || (img[x][y-1][z] && img[x+1][y-1][z]) 
									  || (img[x][y][z-1] && img[x+1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y-1][z-1]) ) ) return false;
		
			if (img[x-1][y+1][z+1] && !( (img[x-1][y][z] && img[x-1][y+1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z+1]) 
									  || (img[x][y+1][z] && img[x][y+1][z+1]) 
									  || (img[x][y+1][z] && img[x-1][y+1][z]) 
									  || (img[x][y][z+1] && img[x-1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y+1][z+1]) ) ) return false;
		
			if (img[x+1][y-1][z+1] && !( (img[x+1][y][z] && img[x+1][y-1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z+1]) 
									  || (img[x][y-1][z] && img[x][y-1][z+1]) 
									  || (img[x][y-1][z] && img[x+1][y-1][z]) 
									  || (img[x][y][z+1] && img[x+1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y-1][z+1]) ) ) return false;
		
			if (img[x+1][y+1][z-1] && !( (img[x+1][y][z] && img[x+1][y+1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z-1]) 
									  || (img[x][y+1][z] && img[x][y+1][z-1]) 
									  || (img[x][y+1][z] && img[x+1][y+1][z]) 
									  || (img[x][y][z-1] && img[x+1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y+1][z-1]) ) ) return false;
		
			if (img[x+1][y+1][z+1] && !( (img[x+1][y][z] && img[x+1][y+1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z+1]) 
									  || (img[x][y+1][z] && img[x][y+1][z+1]) 
									  || (img[x][y+1][z] && img[x+1][y+1][z]) 
									  || (img[x][y][z+1] && img[x+1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y+1][z+1]) ) ) return false;
		}
		
		return true;
	}
	
	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectBoundary(boolean[][][] obj, int nx, int ny, int nz) {
		return objectInsideBoundary(obj,nx,ny,nz,6);
	}
	
	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectInsideBoundary(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		boolean[][][]		boundary = new boolean[nx][ny][nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boundary[x][y][z] = false;
			if (obj[x][y][z]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (!obj[x+i][y+j][z+l]) boundary[x][y][z] = true;
						}
					}
				}
			}
		}
		return boundary;
	}

	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectOutsideBoundary(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		boolean[][][]		boundary = new boolean[nx][ny][nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boundary[x][y][z] = false;
			if (!obj[x][y][z]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (obj[x+i][y+j][z+l]) boundary[x][y][z] = true;
						}
					}
				}
			}
		}
		return boundary;
	}
	
	public static boolean[][][] createObjectMask(float[][][] image, float val, int nx, int ny, int nz) {
		int 		x,y,z;
		boolean[][][]  	objMask;

		// uses only values over the threshold, if mask used
		objMask = new boolean[nx][ny][nz];
		for (x=0;x<nx;x++)for (y=0;y<ny;y++)for (z=0;z<nz;z++) {
			if (image[x][y][z] < val ){
				objMask[x][y][z] = false;
			}else{
				System.out.println("x: " + x);
				System.out.println("y: " + y);
				System.out.println("z: " + z);
				objMask[x][y][z] = true;
			}
		}
		// remove the boundary from the computations
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++) {
				objMask[x][y][0] = false;
				objMask[x][y][nz-1] = false;
			}
		for (y=0;y<ny;y++)
			for (z=0;z<nz;z++) {
				objMask[0][y][z] = false;
				objMask[nx-1][y][z] = false;
			}
		for (z=0;z<nz;z++)
			for (x=0;x<nx;x++) {
				objMask[x][0][z] = false;
				objMask[x][ny-1][z] = false;
			}

		return objMask;
	} // createObjectMask

	/**
	 *  compute a simplistic mean curvature function on the object boundary
	 */
	public static final float[][][] objectCurvature(boolean[][][] obj, int nx, int ny, int nz) {
		float   	num, den;
		float		u,xp,xm,yp,ym,zp,zm;
		float		xpyp,xpym,xpzp,xpzm,xmyp,xmym,xmzp,xmzm,ypzp,ypzm,ymzp,ymzm;
		float		xpypzp,xpypzm,xpymzp,xpymzm,xmypzp,xmypzm,xmymzp,xmymzm;
		float   	ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		float		val;
		float[][][]		curv = new float[nx][ny][nz];
		float[][][]		level = new float[nx][ny][nz];
		int 		boundary;
		
		// create a 'levelset' map
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (obj[x][y][z]) {
				// check for boundary
				boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
					else if ( (i*i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
					else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
				}
				// attribute a levelset for the different boundaries
				if (boundary==0) level[x][y][z] = -1.0f;
				else level[x][y][z] = 0.0f;
				/*
				else if (boundary==6) level[x][y][z] = 0.0f;
				else if (boundary==18) level[x][y][z] = -SQR2+1.0f;
				else if (boundary==26) level[x][y][z] = -SQR3+1.0f;
				*/
			} else {
				level[x][y][z] = 1.0f;
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (level[x][y][z]==0) {
				// set the level set values
				u = level[x][y][z];
				
				// 6 - neighbors		
				xp = level[x+1][y][z];
				xm = level[x-1][y][z];
				yp = level[x][y+1][z];
				ym = level[x][y+1][z];
				zp = level[x][y][z+1];
				zm = level[x][y][z-1];
				
				// 18-neighbors		
				xpyp = level[x+1][y+1][z];
				xmyp = level[x-1][y+1][z];
				xpym = level[x+1][y-1][z];
				xmym = level[x-1][y-1][z];
				xpzp = level[x+1][y][z+1];
				xmzp = level[x-1][y][z+1];
				xpzm = level[x+1][y][z-1];
				xmzm = level[x-1][y][z-1];
				ypzp = level[x][y+1][z+1];
				ymzp = level[x][y-1][z+1];
				ypzm = level[x][y+1][z-1];
				ymzm = level[x][y-1][z-1];
				
				// 26-neighbors
				xpypzm = level[x+1][y+1][z-1];
				xpymzm = level[x+1][y-1][z-1];
				xmypzm = level[x-1][y+1][z-1];
				xmymzm = level[x-1][y-1][z-1];
				xpypzp = level[x+1][y+1][z+1];
				xpymzp = level[x+1][y-1][z+1];
				xmypzp = level[x-1][y+1][z+1];
				xmymzp = level[x-1][y-1][z+1];
			
				// central differences ?
				ux = 0.25f*( xp + xpyp + xpzp + xpypzp 
							-xm - xmyp - xmzp - xmypzp );
				
				uy = 0.25f*( yp + xpyp + ypzp + xpypzp 
							-ym - xpym - ymzp - xpymzp );
				
				uz = 0.25f*( zp + xpzp + ypzp + xpypzp 
							-zm - xpzm - ypzm - xpypzm );
		
				uxx = 0.0625f*( 4*xp + 2*xpyp + 2*xpzp + 2*xpym + 2*xpzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*xm + 2*xmyp + 2*xmzp + 2*xmym + 2*xmzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*yp + 2*zp + 2*ym + 2*zm + ypzp + ymzp + ypzm + ymzm);
							 
				uyy = 0.0625f*( 4*yp + 2*xpyp + 2*ypzp + 2*xmyp + 2*ypzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*ym + 2*xpym + 2*ymzp + 2*xmym + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*zp + 2*xm + 2*zm + xpzp + xmzp + xpzm + xmzm);
							 
				uzz = 0.0625f*( 4*zp + 2*xpzp + 2*ypzp + 2*xmzp + 2*ymzp + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*zm + 2*xpzm + 2*ypzm + 2*xmzm + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*yp + 2*xm + 2*ym + xpyp + xmyp + xpym + xmym);
							 
		
				uxy = 0.0625f*( 2*xpyp + 2*xmym + xpypzp + xpypzm + xmymzp + xmymzm
							  - 2*xpym - 2*xmyp - xpymzp - xpymzm - xmypzp - xmypzm );
							 
				uyz = 0.0625f*( 2*ypzp + 2*ymzm + xpypzp + xmypzp + xpymzm + xmymzm
							  - 2*ymzp - 2*ypzm - xpymzp - xmymzp - xpypzm - xmypzm );
							 
				uzx = 0.0625f*( 2*xpzp + 2*xmzm + xpypzp + xpymzp + xmypzm + xmymzm
							  - 2*xmzp - 2*xpzm - xpypzm - xpymzm - xmypzp - xmymzp );
							 
				// 3D mean curvature
				num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
				den = (float)Math.sqrt(ux*ux + uy*uy + uz*uz);
				den = 2.0f*den*den*den;
				
				if (den>0) {
					curv[x][y][z] = num/den;
				} else {
					curv[x][y][z] = 0.0f;
				}
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}

	/* not meaningful
	public static final float[][][] objectNeighborRatio(boolean[][][] obj, int nx, int ny, int nz) {
		float[][][]		factor = new float[nx][ny][nz];
		int 		boundary;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int Nb=0;
			if (obj[x][y][z]) {
				// check for boundary
				boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
					else if ( (i*i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
					else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
				}
				
				if (boundary>0) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
						if (obj[x+i][y+j][z+l]) {
							Nb++;
						}
					}
				}
			}
			factor[x][y][z] = (Nb-26.0f)/25.0f;
		}
		return factor;
	}
	*/
	/* not such a good idea..
	public static final float[][][] objectDiscreteCurvature(boolean[][][] obj, int nx, int ny, int nz) {
		float   	vx,vy,vz, ngb;
		int			den;
		float		Pmax,proj;
		float[][][]		curv = new float[nx][ny][nz];
		int 		boundary;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			Pmax = 0.0f;
			if (obj[x][y][z]) {
				// check for boundary
				boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
					else if ( (i*i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
					else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
				}
				
				if (boundary>0) {
					// mean direction vector
					vx = 0.0f; vy = 0.0f; vz = 0.0f;
					ngb = 0.0f;
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
						if (obj[x+i][y+j][z+l]) {
							den = i*i+j*j+l*l;
							if (den==1) { vx -= i; vy -= j; vz -= l; ngb ++; }
							else if (den==2) { vx -= i/SQR2; vy -= j/SQR2; vz -= l/SQR2; ngb ++; }
							else if (den==3) { vx -= i/SQR3; vy -= j/SQR3; vz -= l/SQR3; ngb ++; }
						}
					}
					if (ngb > 0.0f) {
						vx = vx/ngb; vy = vy/ngb; vz = vz/ngb;
					} else {
						System.out.print("!!");
					}
				
					// minimum projection
					Pmax = -1.0f;
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
						if (obj[x+i][y+j][z+l]) {
							den = i*i+j*j+l*l;
							proj = -1.0f;
							if (den==1) { proj = i*vx + j*vy + l*vz; }
							else if (den==2) { proj = (i*vx + j*vy + l*vz)/SQR2; }
							else if (den==3) { proj = (i*vx + j*vy + l*vz)/SQR3; }
							if (proj>Pmax) Pmax = proj;
						}
					}
				}
			}
			curv[x][y][z] = Pmax;	
		}
	
		return curv;
	}
	*/
	
	/**
     *  compute the object volume
     */
    public static final int volume(boolean img[][][], int nx, int ny, int nz) {
		int count=0;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) count++;
		}

		return count;
	}
	
	/**
     *  compute the object volume
     */
    public static final int volume(int img[][][], int lb, int nx, int ny, int nz) {
		int count=0;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]==lb) count++;
		}

		return count;
	}
	
	
	/**
     *  compute the number of separate object parts
     */
    public static final int countParts(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int[][][] lb = null;
		
			 if (cObj== 6) lb = connected6Object3D(img,nx,ny,nz);
		else if (cObj==18) lb = connected18Object3D(img,nx,ny,nz);
		else if (cObj==26) lb = connected26Object3D(img,nx,ny,nz);
		
		return (countLabels(lb,nx,ny,nz)-1);
	}
	
	/**
     *  compute the number of holes in the object
     */
    public static final int countHoles(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int[][][] lb = null;
		boolean[][][] obj = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			obj[x][y][z] = !img[x][y][z];
		}
			 if (cBg== 6) lb = connected6Object3D(obj,nx,ny,nz);
		else if (cBg==18) lb = connected18Object3D(obj,nx,ny,nz);
		else if (cBg==26) lb = connected26Object3D(obj,nx,ny,nz);

		return (countLabels(lb,nx,ny,nz)-2);
	}
	
	/**
     *  compute a binary operation (and, or, xor) on two boolean images.
     */
    public static final boolean[][][] binaryOperation(boolean img1[][][], boolean img2[][][], int operator, int nx, int ny, int nz) {
		boolean[][][] res = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			res[x][y][z] = false;
			if ( (operator==AND) && (img1[x][y][z] && img2[x][y][z]) ) res[x][y][z]=true;
			else if ( (operator==OR) && (img1[x][y][z] || img2[x][y][z]) ) res[x][y][z]=true;
			else if ( (operator==XOR) && ( (img1[x][y][z] && !img2[x][y][z]) || (!img1[x][y][z] && img2[x][y][z]) ) ) res[x][y][z] = true;
		}

		return res;
	}
	
	/**
     *  compute a binary operation (and, or, xor) on two boolean images.
     */
    public static final boolean[][][] inverse(boolean img[][][], int nx, int ny, int nz) {
		boolean[][][] res = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			res[x][y][z] = !img[x][y][z];
		}

		return res;
	}
    
    /**
     * Returns the value of the "centile" for the obj and intensity image
     * the number of histogram bins can be specified using Nbins
     */
    public static final float objectCentile(boolean[][][] obj, float[][][] intensity, float centile, int Nbins, int nx, int ny, int nz){
		// get min max
		float Imin = intensity[0][0][0];
		float Imax = intensity[0][0][0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
			if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
		}
		
		float[] hist = new float[Nbins+1];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (obj[x][y][z]) {
				int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
				if (bin<0) bin = 0;
				if (bin>=Nbins) bin = Nbins-1;
				hist[bin]++;
				hist[Nbins]++;	// total number
			}
		}
		float count = 0.0f;
		int bin = 0;
		while (count<centile*hist[Nbins]) {
			count += hist[bin];
			bin++;
		}
		float val = Imin + (bin-1)/(float)Nbins*(Imax-Imin);

		return val;
	}
    
    /**
     * computes the confusion matrix / contingency table of two 3d discrete valued images.
     * Many overlap metrics (Dice, Jaccard, TP, FP , etc etc, can be computed from it.
     * 
     * The i,j element of the output array is the number of voxels for which ref=i and seg=j. 
     * 
     * @return
     */
    public static final int[][] confusionMatrix(int[][][] seg, int[][][] ref, int nx, int ny, int nz){
    	
    	int[] seglabels = listOrderedLabels(seg,nx,ny,nz);
    	int[] reflabels = listOrderedLabels(ref,nx,ny,nz);
    	ArrayList<Integer> alllabels = new ArrayList<Integer>();
    	for(int i=0; i<seglabels.length; i++){
    		alllabels.add(seglabels[i]);
    	}
    	// add any labels from ref not in seg
    	for(int i=0; i<reflabels.length; i++){
    		if(!alllabels.contains(reflabels[i])){
    			alllabels.add(reflabels[i]);
    		}
    	}
    	// sort so binary search will work 
    	Collections.sort(alllabels);
    	int nlabels = alllabels.size();
    	
    	System.out.println("labels:"+alllabels);
    	
    	int[][] confmtx = new int[nlabels][nlabels];
    	
    	int i=-1;
    	int j=-1;
    	for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    		i=Collections.binarySearch(alllabels, ref[x][y][z]); // idx of ref label
    		j=Collections.binarySearch(alllabels, seg[x][y][z]); // idx of seg label
    		
    		if(i<0 || j<0){
    		System.out.println("ref[x][y][z]: " + ref[x][y][z]);
    		System.out.println("seg[x][y][z]: " + seg[x][y][z]);
    		System.out.println("i: " + i);
    		System.out.println("j: " + j);
    		}
    		confmtx[i][j]++;
    	}
    	return confmtx;
    }
					
	/**
	*  compute a skeleton using a fast marching method
	*/
	public static final boolean[][][] simpleSkeleton(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		float			val=0.0f, ngb=0.0f;
		int				x,y,z;
		int[]			vec = new int[3];
		boolean[][][]	skeleton, critical;
		float[][][]		dist;
		int max;
		BinaryTree		tree;
		CriticalPointLUT lut;
		float epsilon = 1e-3f;
		float mindist = 2.0f;
		
		/*
		dist = fastMarchingDistance(obj,nx,ny,nz,connectivity);
		float maxdist = ImageFunctions.maximum(dist,nx,ny,nz);
		*/
		
		tree = new BinaryTree(nx*ny*nz, 3, BinaryTree.MINTREE, BinaryTree.ADAPTATIVE);
		
		if (connectivity==6) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
		else if (connectivity==18) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
		else lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			return null;
		}
		
		/*
		if (connectivity==6) max=2;
		else if (connectivity==18) max=3; 
		else max=4;
		*/
		max=2;
		
		// no processing outside the object
		skeleton = new boolean[nx][ny][nz];
		critical = new boolean[nx][ny][nz];
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			skeleton[x][y][z] = obj[x][y][z];
			critical[x][y][z] = false;
		}
		
		// pre-compute the distance function: must use same algorithm
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if (obj[x][y][z]) {
				// check for inside boundary
				int boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<max) {
							 if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
						else if ( (i* i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
						else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
					}
				}
				// set the distances
					 if (boundary==6)  val = 0.5f;
				else if (boundary==18) val = 0.5f*SQR2;
				else if (boundary==26) val = 0.5f*SQR3;
				// record in the sorting tree
				if (boundary>0) {
					vec[0] = x; vec[1] = y; vec[2] = z;
					tree.addValue(val,vec);
				}
			}
		}
		
		// propagate
		while (tree.isNotEmpty()) {
			// get the next value
			x = tree.getFirstIndex(0);
			y = tree.getFirstIndex(1);
			z = tree.getFirstIndex(2);
			val = tree.getFirst();
			tree.removeFirst();
			
			if (skeleton[x][y][z]) {
				
				// check the topology
				if (lut.get(lut.keyFromPattern(skeleton,x,y,z))) {
				
					// check for endpoint
					int nb = -1;
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
						if (skeleton[x+i][y+j][z+l]) nb++;
					}
					if (nb>1 || val<mindist) {
						
						// set the distance, update label
						//dist[x][y][z] = val;
						skeleton[x][y][z] = false;
						critical[x][y][z] = false;
		
						// find the neighbors
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
							if (i*i+j*j+l*l<max) {
								if (critical[x+i][y+j][z+l]) {
									ngb = val + epsilon;
									//ngb = dist[x+i][y+j][z+l] + epsilon;
									vec[0] = x+i; vec[1] = y+j; vec[2] = z+l;
									tree.addValue(ngb,vec);
								} else if (skeleton[x+i][y+j][z+l]) {
										 if (i*i+j*j+l*l==1) ngb = val + 1.0f;
									else if (i*i+j*j+l*l==2) ngb = val + SQR2;
									else if (i*i+j*j+l*l==3) ngb = val + SQR3;
									vec[0] = x+i; vec[1] = y+j; vec[2] = z+l;
									tree.addValue(ngb,vec);
								}
							}
						}
					}
				} else {
					//dist[x][y][z] = val;
					critical[x][y][z] = true;
				}
			}
		}
		// clean up
		tree.finalize(); tree = null;
		
		return skeleton;
	}//simpleSkeleton
	
	/**
	*  compute a feature transform based on Voronoi diagrams
	*	from: Maurer, Qi and Rghavan, PAMI 25:2, 2003
	*/
	public static final short[][][][] voronoiFeatureTransform(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = new short[nx][ny][nz][3];
		short[][] 		gx = new short[nx][3];
		short[][] 		gy = new short[ny][3];
		short[][] 		gz = new short[nz][3];
		
		// initialize the ft
		for (short z=0;z<nz;z++) {
			for (short y=0;y<ny;y++) {
				for (short x=0;x<nx;x++) {
					if (obj[x][y][z]) {
						ft[x][y][z][0] = x;
						ft[x][y][z][1] = y;
						ft[x][y][z][2] = z;
					} else {
						ft[x][y][z][0] = -1;
						ft[x][y][z][1] = -1;
						ft[x][y][z][2] = -1;
					}
				}
				// compute Voronoi 1D
				computePartialVoronoiDiagramX(ft, gx, y, z, nx, ny, nz);
			}
			for (short x=0;x<nx;x++) {
				// compute Voronoi 2D
				computePartialVoronoiDiagramY(ft, gy, x, z, nx, ny, nz);		
			}
		}
		for (short x=0;x<nx;x++) {
			for (short y=0;y<ny;y++) {
				// compute Voronoi 3D
				computePartialVoronoiDiagramZ(ft, gz, x, y, nx, ny, nz);
			}
		}
		gx = null;
		gy = null;
		gz = null;
		return ft;
	}
	
	private static final void computePartialVoronoiDiagramX(short[][][][] ft, short[][] g, short y, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			
			if (ft[i][y][z][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[i][y][z][0];
					g[l][1] = ft[i][y][z][1];
					g[l][2] = ft[i][y][z][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[i][y][z],xi,0)) l--;
					l++;
					g[l][0] = ft[i][y][z][0];
					g[l][1] = ft[i][y][z][1];
					g[l][2] = ft[i][y][z][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[i][y][z][0] = g[l][0];
			ft[i][y][z][1] = g[l][1];
			ft[i][y][z][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramY(short[][][][] ft, short[][] g, short x, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			
			if (ft[x][i][z][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[x][i][z][0];
					g[l][1] = ft[x][i][z][1];
					g[l][2] = ft[x][i][z][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[x][i][z],xi,1)) l--;
					l++;
					g[l][0] = ft[x][i][z][0];
					g[l][1] = ft[x][i][z][1];
					g[l][2] = ft[x][i][z][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[x][i][z][0] = g[l][0];
			ft[x][i][z][1] = g[l][1];
			ft[x][i][z][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramZ(short[][][][] ft, short[][] g, short x, short y, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			
			if (ft[x][y][i][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[x][y][i][0];
					g[l][1] = ft[x][y][i][1];
					g[l][2] = ft[x][y][i][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[x][y][i],xi,2)) l--;
					l++;
					g[l][0] = ft[x][y][i][0];
					g[l][1] = ft[x][y][i][1];
					g[l][2] = ft[x][y][i][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[x][y][i][0] = g[l][0];
			ft[x][y][i][1] = g[l][1];
			ft[x][y][i][2] = g[l][2];
		}
		return;
	}
	
	/**
	*  compute a feature transform based on Voronoi diagrams
	*	from: Maurer, Qi and Rghavan, PAMI 25:2, 2003
	*/
	public static final short[][] voronoiFeatureTransform(boolean[] obj, int nx, int ny, int nz) {
		short[][]		ft = new short[nx*ny*nz][3];
		short[][] 		gx = new short[nx][3];
		short[][] 		gy = new short[ny][3];
		short[][] 		gz = new short[nz][3];
		
		// initialize the ft
		for (short z=0;z<nz;z++) {
			for (short y=0;y<ny;y++) {
				for (short x=0;x<nx;x++) {
					int xyz = x+nx*y+nx*ny*z;
					if (obj[xyz]) {
						ft[xyz][0] = x;
						ft[xyz][1] = y;
						ft[xyz][2] = z;
					} else {
						ft[xyz][0] = -1;
						ft[xyz][1] = -1;
						ft[xyz][2] = -1;
					}
				}
				// compute Voronoi 1D
				computePartialVoronoiDiagramX(ft, gx, y, z, nx, ny, nz);
			}
			for (short x=0;x<nx;x++) {
				// compute Voronoi 2D
				computePartialVoronoiDiagramY(ft, gy, x, z, nx, ny, nz);		
			}
		}
		for (short x=0;x<nx;x++) {
			for (short y=0;y<ny;y++) {
				// compute Voronoi 3D
				computePartialVoronoiDiagramZ(ft, gz, x, y, nx, ny, nz);
			}
		}
		gx = null;
		gy = null;
		gz = null;
		return ft;
	}
	
	private static final void computePartialVoronoiDiagramX(short[][] ft, short[][] g, short y, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			int iyz = i+nx*y+nx*ny*z;
			
			if (ft[iyz][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[iyz][0];
					g[l][1] = ft[iyz][1];
					g[l][2] = ft[iyz][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[iyz],xi,0)) l--;
					l++;
					g[l][0] = ft[iyz][0];
					g[l][1] = ft[iyz][1];
					g[l][2] = ft[iyz][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			int iyz = i+nx*y+nx*ny*z;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[iyz][0] = g[l][0];
			ft[iyz][1] = g[l][1];
			ft[iyz][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramY(short[][] ft, short[][] g, short x, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			int xiz = x+nx*i+nx*ny*z;
			
			if (ft[xiz][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[xiz][0];
					g[l][1] = ft[xiz][1];
					g[l][2] = ft[xiz][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[xiz],xi,1)) l--;
					l++;
					g[l][0] = ft[xiz][0];
					g[l][1] = ft[xiz][1];
					g[l][2] = ft[xiz][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			int xiz = x+nx*i+nx*ny*z;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[xiz][0] = g[l][0];
			ft[xiz][1] = g[l][1];
			ft[xiz][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramZ(short[][] ft, short[][] g, short x, short y, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			int xyi = x+nx*y+nx*ny*i;
			
			if (ft[xyi][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[xyi][0];
					g[l][1] = ft[xyi][1];
					g[l][2] = ft[xyi][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[xyi],xi,2)) l--;
					l++;
					g[l][0] = ft[xyi][0];
					g[l][1] = ft[xyi][1];
					g[l][2] = ft[xyi][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			int xyi = x+nx*y+nx*ny*i;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[xyi][0] = g[l][0];
			ft[xyi][1] = g[l][1];
			ft[xyi][2] = g[l][2];
		}
		return;
	}
	
	public static final float[][][] voronoiFeatureSquaredDistance(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[][][] dist = new float[nx][ny][nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = ftDistance(pt, ft[x][y][z]);
			} else {
				dist[x][y][z] = -1;
			}
		}
		ft = null;
		
		return dist;
	}
	
	public static final float[][][] voronoiFeatureDistance(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[][][] dist = new float[nx][ny][nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = (float)Math.sqrt(ftDistance(pt, ft[x][y][z]));
			} else {
				dist[x][y][z] = -1;
			}
		}
		ft = null;
		
		return dist;
	}
	
	private static final boolean removeVoronoiFeature(short[] u, short[] v, short[] w, short[] Rd, int d) {
		
		float duR = 0.0f; 
		float dvR = 0.0f; 
		float dwR = 0.0f; 
		for (int i=0;i<3;i++) if (i!=d) {
			duR += (u[i]-Rd[i])*(u[i]-Rd[i]);
			dvR += (v[i]-Rd[i])*(v[i]-Rd[i]);
			dwR += (w[i]-Rd[i])*(w[i]-Rd[i]);
		}
		return ( (w[d]-u[d])*dvR - (w[d]-v[d])*duR - (v[d]-u[d])*dwR - (w[d]-u[d])*(w[d]-v[d])*(v[d]-u[d]) > 0 );
	}
	
	private static final float ftDistance(short[] u, short[] v) {
		
		return (u[0]-v[0])*(u[0]-v[0])
			  +(u[1]-v[1])*(u[1]-v[1])
			  +(u[2]-v[2])*(u[2]-v[2]);
	}
	
	public static final float[][][] signedDistanceFunction(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[][][] dist = new float[nx][ny][nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (!obj[x][y][z] && ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = (float)Math.sqrt(ftDistance(pt, ft[x][y][z]))-0.5f;
			} else {
				dist[x][y][z] = 0;
			}
		}
		// same on the other part
		boolean[][][] bg = new boolean[nx][ny][nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			bg[x][y][z] = !obj[x][y][z];
		}
		ft = voronoiFeatureTransform(bg,nx,ny,nz);
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (!bg[x][y][z] && ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = -(float)Math.sqrt(ftDistance(pt, ft[x][y][z]))+0.5f;
			} else if (!bg[x][y][z]) {
				dist[x][y][z] = 0;
			}
		}
		
		ft = null;
		
		return dist;
	}
	
	public static final float[] signedDistanceFunction(boolean[] obj, int nx, int ny, int nz) {
		short[][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[] dist = new float[nx*ny*nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (obj[xyz] && ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = (float)Math.sqrt(ftDistance(pt, ft[xyz]));
			} else {
				dist[xyz] = 0;
			}
		}
		// same on the other part
		boolean[] bg = new boolean[nx*ny*nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			bg[xyz] = !obj[xyz];
		}
		ft = voronoiFeatureTransform(bg,nx,ny,nz);
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (bg[xyz] && ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = -(float)Math.sqrt(ftDistance(pt, ft[xyz]));
			} else {
				dist[xyz] = 0;
			}
		}
		
		ft = null;
		
		return dist;
	}
	
	/**
	 *	check for convex points / concave points:
	 *  count the number of added boundaries (26-C).
	 */
	private static final float convexityScore(boolean[][][] obj, int x, int y, int z) {
		float boundaries=0;
		
		// 6-C
		if (obj[x-1][y][z]) boundaries++;
		if (obj[x+1][y][z]) boundaries++;
		if (obj[x][y-1][z]) boundaries++;
		if (obj[x][y+1][z]) boundaries++;
		if (obj[x][y][z-1]) boundaries++;
		if (obj[x][y][z+1]) boundaries++;
		
		// 18-C + implied 18-C
		if (obj[x-1][y-1][z]) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) ) boundaries++;
		if (obj[x-1][y+1][z]) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) ) boundaries++;
		if (obj[x+1][y-1][z]) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) ) boundaries++;
		if (obj[x+1][y+1][z]) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) ) boundaries++;
		
		if (obj[x][y-1][z-1]) boundaries++;
		else if ( (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;
		if (obj[x][y-1][z+1]) boundaries++;
		else if ( (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;
		if (obj[x][y+1][z-1]) boundaries++;
		else if ( (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;
		if (obj[x][y+1][z+1]) boundaries++;
		else if ( (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;
		
		if (obj[x-1][y][z-1]) boundaries++;
		else if ( (obj[x][y][z-1]) && (obj[x-1][y][z]) ) boundaries++;
		if (obj[x+1][y][z-1]) boundaries++;
		else if ( (obj[x][y][z-1]) && (obj[x+1][y][z]) ) boundaries++;
		if (obj[x-1][y][z+1]) boundaries++;
		else if ( (obj[x][y][z+1]) && (obj[x-1][y][z]) ) boundaries++;
		if (obj[x+1][y][z+1]) boundaries++;
		else if ( (obj[x][y][z+1]) && (obj[x+1][y][z]) ) boundaries++;
		
		// 26-C and implied 26-C
		if (obj[x-1][y-1][z-1]) boundaries++;
		else if ( (obj[x-1][y-1][z]) && (obj[x][y-1][z-1]) ) boundaries++;
		else if ( (obj[x][y-1][z-1]) && (obj[x-1][y][z-1]) ) boundaries++;
		else if ( (obj[x-1][y][z-1]) && (obj[x-1][y-1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x+1][y-1][z-1]) boundaries++;
		else if ( (obj[x+1][y-1][z]) && (obj[x][y-1][z-1]) ) boundaries++;
		else if ( (obj[x][y-1][z-1]) && (obj[x+1][y][z-1]) ) boundaries++;
		else if ( (obj[x+1][y][z-1]) && (obj[x+1][y-1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x-1][y+1][z-1]) boundaries++;
		else if ( (obj[x-1][y+1][z]) && (obj[x][y+1][z-1]) ) boundaries++;
		else if ( (obj[x][y+1][z-1]) && (obj[x-1][y][z-1]) ) boundaries++;
		else if ( (obj[x-1][y][z-1]) && (obj[x-1][y+1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x-1][y-1][z+1]) boundaries++;
		else if ( (obj[x-1][y-1][z]) && (obj[x][y-1][z+1]) ) boundaries++;
		else if ( (obj[x][y-1][z+1]) && (obj[x-1][y][z+1]) ) boundaries++;
		else if ( (obj[x-1][y][z+1]) && (obj[x-1][y-1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x-1][y+1][z+1]) boundaries++;
		else if ( (obj[x-1][y+1][z]) && (obj[x][y+1][z+1]) ) boundaries++;
		else if ( (obj[x][y+1][z+1]) && (obj[x-1][y][z+1]) ) boundaries++;
		else if ( (obj[x-1][y][z+1]) && (obj[x-1][y+1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x+1][y-1][z+1]) boundaries++;
		else if ( (obj[x+1][y-1][z]) && (obj[x][y-1][z+1]) ) boundaries++;
		else if ( (obj[x][y-1][z+1]) && (obj[x+1][y][z+1]) ) boundaries++;
		else if ( (obj[x+1][y][z+1]) && (obj[x+1][y-1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x+1][y+1][z-1]) boundaries++;
		else if ( (obj[x+1][y+1][z]) && (obj[x][y+1][z-1]) ) boundaries++;
		else if ( (obj[x][y+1][z-1]) && (obj[x+1][y][z-1]) ) boundaries++;
		else if ( (obj[x+1][y][z-1]) && (obj[x+1][y+1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x+1][y+1][z+1]) boundaries++;
		else if ( (obj[x+1][y+1][z]) && (obj[x][y+1][z+1]) ) boundaries++;
		else if ( (obj[x][y+1][z+1]) && (obj[x+1][y][z+1]) ) boundaries++;
		else if ( (obj[x+1][y][z+1]) && (obj[x+1][y+1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;

		return boundaries/26.0f;
	}

	/**
     *  compute the object convexity map.
     */
    public static final float[][][] convexityMap(boolean img[][][], int nx, int ny, int nz) {
		float[][][] score = new float[nx][ny][nz];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			score[x][y][z] = convexityScore(img,x,y,z);
		}

		return score;
	}
	
	/**
     *  compute the object mean convexity
     */
    public static final float meanConvexity(boolean img[][][], int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += convexityScore(img,x,y,z);
					count++;
				}
			}
		}

		return score/(float)count;
	}
	
	/**
     *  compute the object mean absolute convexity
     */
    public static final float absConvexity(boolean img[][][], int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += Numerics.abs(convexityScore(img,x,y,z)-0.5f);
					count++;
				}
			}
		}

		return score/(float)count;
	}
	
	/**
     *  compute the object convexity variance
     */
    public static final float stdConvexity(boolean img[][][], float mean, int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += (convexityScore(img,x,y,z)-mean)*(convexityScore(img,x,y,z)-mean);
					count++;
				}
			}
		}

		return (float)Math.sqrt(score/(float)(count-1));
	}
	
	/**
     *  compute the object center
     */
    public static final float[] center(boolean img[][][], int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(int img[][][], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    /**
     * Computes the center of mass of the signed distance function given by img.
     * Negative values are inside of the object.
     */
    public static final float[] center(float img[][][], int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<=0) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(int img[], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		int ind = -1;
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			ind = x + y*nx + z*nx*ny;
			if (img[ind]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(byte img[], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		int ind = -1;
		for (int n=0;n<3;n++) center[n] = 0.0f;
		ind = 0 ;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[ind]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
			ind++;
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
	
	/**
     *  compute the object deviation from the center
     */
    public static final float[] deviation(boolean img[][][], float[] center, int nx, int ny, int nz) {
		float[] std = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) std[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				std[0] += (x-center[0])*(x-center[0]);
				std[1] += (y-center[1])*(y-center[1]);
				std[2] += (z-center[2])*(z-center[2]);
				count++;
			}
		}
		if (count>1) {
			std[0] = (float)Math.sqrt( std[0]/(count-1) );
			std[1] = (float)Math.sqrt( std[1]/(count-1) );
			std[2] = (float)Math.sqrt( std[2]/(count-1) );
		}
		return std;
	}

	/**
     *  compute the mean object distance to a point
     */
    public static final float meanDistance(boolean img[][][], float x0, float y0, float z0, int nx, int ny, int nz, float rx, float ry, float rz) {
		float dist = 0.0f;
		float count = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist += Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				count++;
			}
		}
		if (count>0) dist = dist/count;
		return dist;
	}

	/**
     *  compute the variance of object distance to a point
     */
    public static final float stdDistance(boolean img[][][], float x0, float y0, float z0, float mean, int nx, int ny, int nz, float rx, float ry, float rz) {
		float var = 0.0f;
		float count = 0.0f;
		float dist;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				var += (dist-mean)*(dist-mean);
				count++;
			}
		}
		if (count>1) var = var/(count-1);
		return (float)Math.sqrt(var);
	}

	/**
     *  compute the mean object distance to a point
     */
    public static final float[] meanDirection(boolean img[][][], float x0, float y0, float z0, int nx, int ny, int nz, float rx, float ry, float rz) {
		float[] dir = new float[3];
		for (int i=0;i<3;i++) dir[i] = 0.0f;
		float count = 0.0f;
		float dist = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				dir[0] += (x-x0)*rx/dist;
				dir[1] += (y-y0)*ry/dist;
				dir[2] += (z-z0)*rz/dist;
				count++;
			}
		}
		if (count>0) {
			dist = (float)Math.sqrt( dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2] );
			for (int i=0;i<3;i++) dir[i] = dir[i]/dist;
		}
		return dir;
	}

	/**
     *  compute the variance of object distance to a point
     */
    public static final float[] stdDirection(boolean img[][][], float x0, float y0, float z0, float[] mean, int nx, int ny, int nz, float rx, float ry, float rz) {
		float[] var = new float[3];
		for (int i=0;i<3;i++) var[i] = 0.0f;
		float count = 0.0f, dist;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				var[0] += ( (x-x0)*rx/dist - mean[0] )*( (x-x0)*rx/dist - mean[0] );
				var[1] += ( (y-y0)*ry/dist - mean[1] )*( (y-y0)*ry/dist - mean[1] );
				var[2] += ( (z-z0)*rz/dist - mean[2] )*( (z-z0)*rz/dist - mean[2] );
				count++;
			}
		}
		if (count>1) {
			for (int i=0;i<3;i++) var[i] = (float)Math.sqrt(var[i]/(count-1));
		}
		return var;
	}

	/**
     *  compute the area on the boundary between two objects.
	 */
    public static final float sharedBoundaryArea(float img[][][], float id1, float id2, int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]==id1) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if (img[x+i][y+j][z+l]==id2) {
						// shared boundary point: add the corresponding surface patch
						if (i*i==1) {
							area += ry*rz;
						} else if (j*j==1) {
							area += rz*rx;
						} else if (l*l==1) {
							area += rx*ry;
						}
					}
				}
			}
		}
		return area;
	}

	/**
     *  compute the area on the boundary between two objects.
	 */
    public static final float sharedBoundaryArea(byte img[][][], byte id1, byte id2, int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]==id1) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if (img[x+i][y+j][z+l]==id2) {
						// shared boundary point: add the corresponding surface patch
						if (i*i==1) {
							area += ry*rz;
						} else if (j*j==1) {
							area += rz*rx;
						} else if (l*l==1) {
							area += rx*ry;
						}
					}
				}
			}
		}
		return area;
	}

	/**
     *  compute the area on the object boundary
	 */
    public static final float boundaryArea(boolean img[][][], int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
						if (!img[x+i][y+j][z+l]) {
							// boundary point: add the corresponding surface patch
							if (i*i==1) {
								area += ry*rz;
							} else if (j*j==1) {
								area += rz*rx;
							} else if (l*l==1) {
								area += rx*ry;
							}
						}
					}
				}
			}
		}
		return area;
	}


}//ObjectProcessing class

