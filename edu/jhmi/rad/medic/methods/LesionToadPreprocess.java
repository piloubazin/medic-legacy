package edu.jhmi.rad.medic.methods;

import java.io.*;
import java.util.*;

import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.structures.*;
import edu.jhmi.rad.medic.utilities.*;

/**
 *
 *  This algorithm performs a simple alignment of the image with a topology template
 *	<p>
 *	The algorithm handles all the little things needed for image cropping, conversion
 * 	to an image buffer, padding, etc.
 *
 *	@version    March 2007
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class LesionToadPreprocess {
		
	// data buffers
	private 	float[][][][]		images;  			// original images
	private 	byte[][][]			topology;  			// topology template
	private		float[]			transform; 				// rigid transform from atlas space to image space
	private 	int				nix,niy,niz;   			// image dimensions
	private 	int				mix,miy,miz;   			// cropped image dimensions
	private 	float			rix,riy,riz;   			// image resolutions
	private 	int				x0i,y0i,z0i;   			// cropped image origin
	private 	int				xNi,yNi,zNi;   			// cropped image final point
	private 	float			xCi,yCi,zCi;   			// coordinates of the center
	private 	float			xCt,yCt,zCt;   			// coordinates of the center
	private 	int				ntx,nty,ntz;   			// template dimensions
	private 	float			rtx,rty,rtz;   			// template resolutions
	private   	float			orx,ory,orz;			// image axis orientation
	private		int				nc;						// number of image channels
	private		boolean			cropped;				// check whether the images are the originals or the cropped ones
	private		boolean			transformed;			// check whether the topology is the originals or the cropped ones
	private		boolean			normalized;				// check whether the images are normalized in [0,1]
	
	// structure parameters
	private 	int 		nobj;    	// number of structures
	private		byte		bgLabel;    // label of the background (usually 1)
	private		float		bgRatio;	// ratio of image intensity used as background
	private		int			bs;	// the amount of extra borders
	private		float[]		Imin,Imax;	// image min,max
	private		float[]		Ilow,Ihigh;	// image robust normalised min,max
	//private		float[][]	sliceLow, sliceHigh;
	// for debug and display
	static final boolean		debug				=	false;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
	
	public static final int ORI_L2R = 1;
	public static final int ORI_R2L = 2;
	public static final int ORI_A2P = 3;
	public static final int ORI_P2A = 4;
	public static final int ORI_I2S = 5;
	public static final int ORI_S2I = 6;
	
	public static final int AXIAL = 10;
	public static final int CORONAL = 20;
	public static final int SAGITTAL = 30;
	
	//private static final     boolean  robustMaxMin = false;
	private static final     boolean  robustMaxMin = true;
	
	/**
	 *  constructor
	 */
	public LesionToadPreprocess(float[][] img_, int nc_,
									int nix_, int niy_, int niz_,
									float rix_, float riy_, float riz_,
									int orient_, int orx_, int ory_, int orz_,
									byte[][][] topo_, 
									int ntx_, int nty_, int ntz_,
									float rtx_, float rty_, float rtz_,
									int nobjs_, byte bgLb_, float bgRatio_,
									int bdSz_) {
		if (debug)
			if (robustMaxMin)
				System.out.println("Maximum Minimum Mode : ROBUST" );
			else
				System.out.println("Maximum Minimum Mode : NORMAL" );
										
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		xCi = nix/2.0f;
		yCi = niy/2.0f;
		zCi = niz/2.0f;
		cropped = false;
		normalized = false;
		
		// create the image
		images = new float[nc][nix][niy][niz];
		for (int c=0;c<nc;c++) for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			images[c][x][y][z] = img_[c][x + nix*y + nix*niy*z];
		}
		
		topology = topo_;
		ntx = ntx_;
		nty = nty_;
		ntz = ntz_;
		rtx = rtx_;
		rty = rty_;
		rtz = rtz_;
		xCt = ntx/2.0f;
		yCt = nty/2.0f;
		zCt = ntz/2.0f;
		transformed = false;
		
		nobj = nobjs_;
		bgLabel = bgLb_;
		bgRatio = bgRatio_;
		bs = bdSz_;
		
		// orientation: initialize a rotation
		transform = new float[6];
		// note : assumes a rotation around the image center
		if (orient_==AXIAL) {
			transform[0] = 0.0f;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
			
			if (orx_==ORI_L2R) orx = -1.0f;
			else orx = 1.0f;
			if (ory_==ORI_P2A) ory = -1.0f;
			else ory = 1.0f;
			if (orz_==ORI_S2I) orz = -1.0f;
			else orz = 1.0f;
		} else if (orient_==CORONAL) {
			transform[0] = -ISQRT2;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
			
			if (orx_==ORI_L2R) orx = -1.0f;
			else orx = 1.0f;
			if (ory_==ORI_I2S) ory = -1.0f;
			else ory = 1.0f;
			if (orz_==ORI_P2A) orz = -1.0f;
			else orz = 1.0f;
		} else if (orient_==SAGITTAL) {
			transform[0] = -0.5f;
			transform[1] = -0.5f;
			transform[2] = 0.5f;
			
			if (orx_==ORI_P2A) orx = -1.0f;
			else orx = 1.0f;
			if (ory_==ORI_I2S) ory = -1.0f;
			else ory = 1.0f;
			if (orz_==ORI_R2L) orz = -1.0f;
			else orz = 1.0f;
		} else {
			// default is axial
			transform[0] = 0.0f;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
			
			orx = 1.0f;
			ory = 1.0f;
			orz = 1.0f;
		}
		transform[3] = 0.0f;
		transform[4] = 0.0f;
		transform[5] = 0.0f;
		
		// init the image cropping					
		x0i = 0;
		y0i = 0;
		z0i = 0;
		xNi = nix-1;
		yNi = niy-1;
		zNi = niz-1;
		
		// compute min,max amd robust min,max
		Imin = new float[nc];
		Imax = new float[nc];
		Ilow = new float[nc];
		Ihigh = new float[nc];
		/*sliceLow = new float[nc][niz];
		sliceHigh = new float[nc][niz];*/
				for (int c=0;c<nc;c++) {
			Imin[c] = images[c][0][0][0];
			Imax[c] = images[c][0][0][0];
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				//map any negative value in the image to zero
				if (images[c][x][y][z] < 0.0f){
					images[c][x][y][z] = 0.0f;
					Imin[c]  = 0.0f;
				}else{
					if (images[c][x][y][z] > Imax[c]) Imax[c] = images[c][x][y][z];
					if (images[c][x][y][z] < Imin[c]) Imin[c] = images[c][x][y][z];
				}
			}
			if (robustMaxMin){
				Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, nix, niy, niz );
				Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, nix, niy, niz );
				float[][][] shuflled_image  = new float[niz][nix][niy];
				for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
					shuflled_image[z][x][y] = images[c][x][y][z];
				}
				/*for (int z=0;z<niz;z++){
					sliceLow[c][z] = ImageFunctions.robustMinimum(shuflled_image[z], 0.01f, 2, nix, niy );
					sliceHigh[c][z] = ImageFunctions.robustMaximum(shuflled_image[z], 0.01f, 2, nix, niy );
				}*/
				
				
				
			}else{
				Ilow[c] = 0.0f;
				Ihigh[c] = 1.0f;	
			}
				
					
			 
			
			if (debug) System.out.println("image "+(c+1)+" range: ["+Imin[c]+", "+Imax[c]+"]\n");
			if (robustMaxMin && debug) System.out.println("image "+(c+1)+" eff. range: ["+Ilow[c]+", "+Ihigh[c]+"]\n");
		}
		
		if (debug) System.out.println("initialisation: done\n");
	}

	final public void finalize() {
		images = null;
		transform = null;
		System.gc();
	}
	
	public final float[][][][] getImages() { return images; }
    
	public final void setTransform(float[] trans) {
        if (trans.length==6) transform = trans;
		else System.err.println("wrong transform type");
    }
    
    public final float[] getTransform() {
        return transform;
    }

    public final float[] getIntensityMax() {
    	if (robustMaxMin)
    		return Ihigh;
    	else
    		return Imax;
    }
    
    public float[] getIntensityScale() {
		float[] scale = new float[nc];
		
		for (int c=0;c<nc;c++) 
			if (robustMaxMin)
				scale[c] = (Ihigh[c]-Ilow[c]);
			else
				scale[c] = (Imax[c]-Imin[c]);
		
		return scale;
	}

    public final float[] exportTransform() {
		float[] trans = new float[6];
		for (int i=0;i<6;i++)
			trans[i] = transform[i];
        return trans;
    }
	
	public final String displayTransform() {
		String info;
		
		info = "transform: ("+transform[0]+", "+transform[1]+", "+transform[2]+", "
						     +transform[3]+", "+transform[4]+", "+transform[5]+")\n";
		
		return info;
	}
	
	/** current dimensions */
	public final int[] getCurrentImageDimensions() {
		int[] dim = new int[3];
        if (cropped) {
			dim[0] = mix; dim[1] = miy; dim[2] = miz;
		} else {
			dim[0] = nix; dim[1] = niy; dim[2] = niz;
		}
		if (debug) System.out.print("current image dimensions: "+dim[0]+" x "+dim[1]+" x "+dim[2]+"\n");
        
		return dim;
    }
	
	public final int[] getOriginalImageArrayDimensions(int size) {
		int[] dim = new int[4];
        dim[0] = nix; dim[1] = niy; dim[2] = niz; dim[3] = size;
		
		return dim;
    }
	
	public final int getOriginalImageSize() {
		return nix*niy*niz;
    }
	
	public final int[] getOriginalImageDimensions() {
		int[] dim = new int[3];
        dim[0] = nix; dim[1] = niy; dim[2] = niz;
		
		return dim;
    }
	
	public final int[] getCroppedImageDimensions() {
		int[] dim = new int[3];
        dim[0] = mix; dim[1] = miy; dim[2] = miz;
		
		return dim;
    }
	
	/** current resolutions, with the sign corresponding to the axis orientation */
	public final float[] getSignedImageResolutions() {
		float[] res = new float[3];
        res[0] = rix/orx; res[1] = riy/ory; res[2] = riz/orz;
		
		if (debug) System.out.print("current image resolutions: "+res[0]+" x "+res[1]+" x "+res[2]+"\n");
        
		return res;
    }
	
	public final float[] getImageResolutions() {
		float[] res = new float[3];
        res[0] = rix; res[1] = riy; res[2] = riz;
		
		if (debug) System.out.print("current image resolutions: "+res[0]+" x "+res[1]+" x "+res[2]+"\n");
        
		return res;
    }
	
	public final float[] getImageOrientations() {
		float[] ori = new float[3];
        ori[0] = orx; ori[1] = ory; ori[2] = orz;
		
		if (debug) System.out.print("current image orientations: "+ori[0]+" x "+ori[1]+" x "+ori[2]+"\n");
        
		return ori;
    }
	
	/** current dimensions */
	public final void updateTransformedTemplate(byte[][][] tpl) {	
		topology = tpl;
		transformed = true;
	}	
	
    /** sets the bounding box for cropping from a vector image */
	public void findCroppingBoundaries() {
		x0i = nix-1;
        xNi = 0;
        y0i = niy-1;
        yNi = 0;
        z0i = niz-1;
        zNi = 0;

		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		// check all images 
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			for (int c=0;c<nc;c++) {
				//if ( ( (normalized) && (images[c][x][y][z]> bgRatio) ) 
					//|| (images[c][x][y][z]> Imin[c] + bgRatio*(Imax[c]-Imin[c]) ) ) {
				if ( ( (normalized) && (images[c][x][y][z]> bgRatio) ) 
						|| (robustMaxMin && images[c][x][y][z]> bgRatio*(Ihigh[c]-Ilow[c]) ) 
						|| ( !robustMaxMin && images[c][x][y][z]> Imin[c] + bgRatio*(Imax[c]-Imin[c]) ) ) {
					if (x < x0i) x0i = x;
					if (x > xNi) xNi = x;
					if (y < y0i) y0i = y;
					if (y > yNi) yNi = y;
					if (z < z0i) z0i = z;
					if (z > zNi) zNi = z;
				}
			}
		}
		
		// check for the topology: from topology to image space
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			if (topology[x][y][z]>bgLabel) {
				// get image coordinates
				float[] XI = computeImageCoordinates(x,y,z,transform,rotation.getMatrix());
				
				if (XI[0] < x0i) x0i = Numerics.floor(XI[0]);
				if (XI[0] > xNi) xNi = Numerics.ceil(XI[0]);
				if (XI[1] < y0i) y0i = Numerics.floor(XI[1]);
				if (XI[1] > yNi) yNi = Numerics.ceil(XI[1]);
				if (XI[2] < z0i) z0i = Numerics.floor(XI[2]);
				if (XI[2] > zNi) zNi = Numerics.ceil(XI[2]);
			}
		}
        // debug
        if (debug) System.out.print("boundaries: ["+x0i+","+xNi+"] ["+y0i+","+yNi+"] ["+z0i+","+zNi+"]\n");
        
        return;
    }
	
	/** crop the images */
	public void cropImages() {
		if (cropped) return;
		
		mix = xNi-x0i+1+2*bs;
		miy = yNi-y0i+1+2*bs;
		miz = zNi-z0i+1+2*bs;
		
		// update the transform parameters
		xCi = mix/2.0f;
		yCi = miy/2.0f;
		zCi = miz/2.0f;
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		float[][] R = rotation.getMatrix();
		transform[3] += R[0][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[0][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[0][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		transform[4] += R[1][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[1][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[1][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		transform[5] += R[2][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[2][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[2][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		
		float[][][][] smaller = new float[nc][mix][miy][miz];

		for (int c=0;c<nc;c++) {
			for (int x=0;x<mix;x++) for (int y=0;y<miy;y++) for (int z=0;z<miz;z++) {
				//if (normalized) smaller[c][x][y][z] = 0.0f;
				//else smaller[c][x][y][z] = Imin[c];
				smaller[c][x][y][z] = 0.0f;
			}
			for (int x=x0i;x<=xNi;x++) {
				for (int y=y0i;y<=yNi;y++) {
					for (int z=z0i;z<=zNi;z++) {
						if ( (x<0) || (x>=nix) || (y<0) || (y>=niy) || (z<0) || (z>=niz) )
							smaller[c][x-x0i+bs][y-y0i+bs][z-z0i+bs] = 0.0f;
						else
							smaller[c][x-x0i+bs][y-y0i+bs][z-z0i+bs] = images[c][x][y][z];
					}
				}
			}
		}
		// replace the original images
		images = smaller;
		
		cropped = true;
		
		return;
	}
	
	/** uncrop the images */
	public void uncropImages() {
		if (!cropped) return;
		
		// update the transform parameters
		xCi = nix/2.0f;
		yCi = niy/2.0f;
		zCi = niz/2.0f;
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		float[][] R = rotation.getMatrix();
		transform[3] -= R[0][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[0][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[0][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		transform[4] -= R[1][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[1][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[1][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		transform[5] -= R[2][0]*(mix/2.0f+x0i-bs-nix/2.0f)*rix/orx+R[2][1]*(miy/2.0f+y0i-bs-niy/2.0f)*riy/ory+R[2][2]*(miz/2.0f+z0i-bs-niz/2.0f)*riz/orz;
		
		float[][][][] larger = new float[nc][nix][niy][niz];

		for (int c=0;c<nc;c++) {

			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				//if (normalized) larger[c][x][y][z] = 0.0f;
				//else larger[c][x][y][z] = Imin[c];
				larger[c][x][y][z] = 0.0f;
			}		
			for (int x=Numerics.max(x0i-bs,0);x<=Numerics.min(xNi+bs,nix-1);x++) {
				for (int y=Numerics.max(y0i-bs,0);y<=Numerics.min(yNi+bs,niy-1);y++) {
					for (int z=Numerics.max(z0i-bs,0);z<=Numerics.min(zNi+bs,niz-1);z++) {
						larger[c][x][y][z] = images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs];
					}
				}
			}
		}
		images = larger;
		cropped = false;
		return;
	}

	/** uncrop the image, send it to a 1D buffer*/
	public float[] uncropAndBuffer(float[][][] src, float bgVal) {
		float[] larger = new float[nix*niy*niz];

		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			larger[x + nix*y + nix*niy*z] = bgVal;
		}		
		for (int x=Numerics.max(x0i-bs,0);x<=Numerics.min(xNi+bs,nix-1);x++) {
			for (int y=Numerics.max(y0i-bs,0);y<=Numerics.min(yNi+bs,niy-1);y++) {
				for (int z=Numerics.max(z0i-bs,0);z<=Numerics.min(zNi+bs,niz-1);z++) {
					larger[x + nix*y + nix*niy*z] = src[x-x0i+bs][y-y0i+bs][z-z0i+bs];
				}
			}
		}
		
		return larger;
	}

	/** uncrop the image, send it to a 1D buffer*/
	public byte[] uncropAndBuffer(byte[][][] src, byte bgVal) {
		byte[] larger = new byte[nix*niy*niz];

		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			larger[x + nix*y + nix*niy*z] = bgVal;
		}		
		for (int x=Numerics.max(x0i-bs,0);x<=Numerics.min(xNi+bs,nix-1);x++) {
			for (int y=Numerics.max(y0i-bs,0);y<=Numerics.min(yNi+bs,niy-1);y++) {
				for (int z=Numerics.max(z0i-bs,0);z<=Numerics.min(zNi+bs,niz-1);z++) {
					larger[x + nix*y + nix*niy*z] = src[x-x0i+bs][y-y0i+bs][z-z0i+bs];
				}
			}
		}
		
		return larger;
	}
	
	public short[] uncropAndBuffer(short[][][] src, short bgVal) {
		short[] larger = new short[nix*niy*niz];

		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			larger[x + nix*y + nix*niy*z] = bgVal;
		}		
		for (int x=Numerics.max(x0i-bs,0);x<=Numerics.min(xNi+bs,nix-1);x++) {
			for (int y=Numerics.max(y0i-bs,0);y<=Numerics.min(yNi+bs,niy-1);y++) {
				for (int z=Numerics.max(z0i-bs,0);z<=Numerics.min(zNi+bs,niz-1);z++) {
					larger[x + nix*y + nix*niy*z] = src[x-x0i+bs][y-y0i+bs][z-z0i+bs];
				}
			}
		}
		
		return larger;
	}

	/** crop the images */
	public float[][][][] generateCroppedImages() {
		if (cropped) return images;
		
		int mx = xNi-x0i+1+2*bs;
		int my = yNi-y0i+1+2*bs;
		int mz = zNi-z0i+1+2*bs;
		
		float[][][][] smaller = new float[nc][mx][my][mz];

		for (int c=0;c<nc;c++) {
			for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
				//if (normalized) smaller[c][x][y][z] = 0.0f;
				//else smaller[c][x][y][z] = Imin[c];
				smaller[c][x][y][z] = 0.0f;
			}
			for (int x=x0i;x<=xNi;x++) {
				for (int y=y0i;y<=yNi;y++) {
					for (int z=z0i;z<=zNi;z++) {
						smaller[c][x-x0i+bs][y-y0i+bs][z-z0i+bs] = images[c][x][y][z];
					}
				}
			}
		}
		
		return smaller;
	}
	
    /** 
	 *	returns the transformed topology template.
	 *	the template is cropped and a border is added.
	 *	note: this transformation should only be attempted once!!
	 */
	public final void transformTopology() {
		if (transformed) return;
		
		// new dimensions
		int mx = xNi-x0i+1+2*bs;
		int my = yNi-y0i+1+2*bs;
		int mz = zNi-z0i+1+2*bs;
		
		byte[][][] smaller = new byte[mx][my][mz];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = bgLabel;
		}
		for (int x=x0i;x<=xNi;x++) for (int y=y0i;y<=yNi;y++) for (int z=z0i;z<=zNi;z++) {
			// coordinates in topology space
			float[] XT = computeTemplateCoordinates(x,y,z,transform,rotation.getMatrix());
			
			smaller[x-x0i+bs][y-y0i+bs][z-z0i+bs] = ImageFunctions.nearestNeighborInterpolation(topology,bgLabel,XT[0],XT[1],XT[2],ntx,nty,ntz);
		}
		topology = smaller;
		transformed = true;
		
		return;
	} // transformTopology

    /** 
	 *	returns the transformed topology template.
	 *	the template is cropped and a border is added.
	 *	note: this transformation should only be attempted once!!
	 */
	public final byte[][][] generateTransformedTopology() {
		if (transformed) return topology;
		
		// new dimensions
		int mx = xNi-x0i+1+2*bs;
		int my = yNi-y0i+1+2*bs;
		int mz = zNi-z0i+1+2*bs;
		
		byte[][][] smaller = new byte[mx][my][mz];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = bgLabel;
		}
		for (int x=x0i;x<=xNi;x++) for (int y=y0i;y<=yNi;y++) for (int z=z0i;z<=zNi;z++) {
			// coordinates in topology space
			float[] XT = computeTemplateCoordinates(x,y,z,transform,rotation.getMatrix());
			
			smaller[x-x0i+bs][y-y0i+bs][z-z0i+bs] = ImageFunctions.nearestNeighborInterpolation(topology,bgLabel,XT[0],XT[1],XT[2],ntx,nty,ntz);
		}
		return smaller;
	} // transformTopology

    /** 
	 *  transform an image aligned with the topology
	 *	the template is cropped and a border is added.
	 */
	public final float[][][] transformNewImage(float[][][] img, float bgVal) {
		
		// new dimensions
		int mx = xNi-x0i+1+2*bs;
		int my = yNi-y0i+1+2*bs;
		int mz = zNi-z0i+1+2*bs;
		
		float[][][] smaller = new float[mx][my][mz];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = bgVal;
		}
		for (int x=x0i;x<=xNi;x++) for (int y=y0i;y<=yNi;y++) for (int z=z0i;z<=zNi;z++) {
			// coordinates in topology space
			float[] XT = computeTemplateCoordinates(x,y,z,transform,rotation.getMatrix());
			
			smaller[x-x0i+bs][y-y0i+bs][z-z0i+bs] = ImageFunctions.linearInterpolation(img,bgVal,XT[0],XT[1],XT[2],ntx,nty,ntz);
		}
		return smaller;
	} // transformTopology

	/** 
	 *	computes the transformed coordinates from template to image space
	 */
	private final float[] computeImageCoordinates(int x,int y,int z,float[] trans, float[][] rot) {
		float[] XI = new float[3];
		XI[0] = (rot[0][0]*((x-xCt)*rtx-trans[3])+rot[1][0]*((y-yCt)*rty-trans[4])+rot[2][0]*((z-zCt)*rtz-trans[5]))*orx/rix + xCi;
		XI[1] = (rot[0][1]*((x-xCt)*rtx-trans[3])+rot[1][1]*((y-yCt)*rty-trans[4])+rot[2][1]*((z-zCt)*rtz-trans[5]))*ory/riy + yCi;
		XI[2] = (rot[0][2]*((x-xCt)*rtx-trans[3])+rot[1][2]*((y-yCt)*rty-trans[4])+rot[2][2]*((z-zCt)*rtz-trans[5]))*orz/riz + zCi;
		
		return XI;
	}
	
	/** 
	 *	computes the transformed coordinates from image to template space
	 */
	private final float[] computeTemplateCoordinates(int x,int y,int z,float[] trans, float[][] rot) {
		float[] XT = new float[3];
		
		XT[0] = (rot[0][0]*(x-xCi)*rix/orx+rot[0][1]*(y-yCi)*riy/ory+rot[0][2]*(z-zCi)*riz/orz+trans[3])/rtx + xCt;
		XT[1] = (rot[1][0]*(x-xCi)*rix/orx+rot[1][1]*(y-yCi)*riy/ory+rot[1][2]*(z-zCi)*riz/orz+trans[4])/rty + yCt;
		XT[2] = (rot[2][0]*(x-xCi)*rix/orx+rot[2][1]*(y-yCi)*riy/ory+rot[2][2]*(z-zCi)*riz/orz+trans[5])/rtz + zCt;
		
		return XT;
	}
	
	/** normalize the images */
	
	public void normalizeImages() {
		if (normalized) return;
		
		if (robustMaxMin)
			/* normalize over the robust min, max  but keeping 0 as the min value 
		   	(equivalent span for all images, but different means if needed, and masking is preserved) */
			for (int c=0;c<nc;c++) for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				images[c][x][y][z] = (images[c][x][y][z])/(Ihigh[c]-Ilow[c]);
				//images[c][x][y][z] = (images[c][x][y][z])/(sliceHigh[c][z]-sliceLow[c][z]);
			}
		else
			// Normalize image to [0,1]
			for (int c=0;c<nc;c++) for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				images[c][x][y][z] = (images[c][x][y][z]-Imin[c])/(Imax[c]-Imin[c]);
			}
	
		normalized = true;
		
		return;
	}
	
	/*public void normalizeImages() {
		if (normalized) return;
		

		for (int c=0;c<nc;c++) for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {

			images[c][x][y][z] = (images[c][x][y][z]-Imin[c])/(Imax[c]-Imin[c]);
			
		}
		normalized = true;
		
		return;
	}*/

	/** align the topology and the intensity center of the images */
	public void alignImagesAndTopolgyCenters(byte[] label, float[][] centroid, int classes) {
		float pix,piy,piz;
		float ptx,pty,ptz;
		float den, w;
		pix = 0; piy = 0; piz = 0; den = 0;
		for (int c=0;c<nc;c++) for (int x=x0i;x<xNi;x++) for (int y=y0i;y<yNi;y++) for (int z=z0i;z<zNi;z++) {
			if (normalized) {
				w = images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs];
			} else {
				
				if (robustMaxMin)
					w = (images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs])/(Ihigh[c]-Ilow[c]);
				else
					w = (images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs]-Imin[c])/(Imax[c]-Imin[c]);
			}
			pix += w*x;
			piy += w*y;
			piz += w*z;
			den += w;
		}
		pix /= den;
		piy /= den;
		piz /= den;
		
		ptx = 0; pty = 0; ptz = 0; den = 0;
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			for (int k=0;k<classes;k++) if (topology[x][y][z]==label[k]) {
				for (int c=0;c<nc;c++) {
					ptx += centroid[c][k]*x;
					pty += centroid[c][k]*y;
					ptz += centroid[c][k]*z;
					den += centroid[c][k];
				}
			}
		}
		ptx /= den; 
		pty /= den;
		ptz /= den;
		
		// transform: translation
		transform[3] = (ptx-xCt)*rtx - (pix-xCi)*rix; 
		transform[4] = (pty-yCt)*rty - (piy-yCi)*riy; 
		transform[5] = (ptz-zCt)*rtz - (piz-zCi)*riz; 
		
		return;
	}

	/** 
	 *  compute the FCM membership functions for initial tissues
	 */
    final public void initialMemberships(float[][][][] mems,float[][] centroid, int classes) {
        float dist, dist0;
        float den,num;
        float mean,count;
		boolean	used;
		
		float[] grouping = new float[classes];
		
        /*if (cropped) {
			mems = new float[mix][miy][miz][classes];
		} else {
			mems = new float[nix][niy][niz][classes];
		}*/	
		for (int x=x0i;x<xNi;x++) for (int y=y0i;y<yNi;y++) for (int z=z0i;z<zNi;z++) {
			den = 0;
			for (int k=0;k<classes;k++) {
				num = 0;
				for (int c=0;c<nc;c++) {
					num += (images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs]-centroid[c][k])
						  *(images[c][x-x0i+bs][y-y0i+bs][z-z0i+bs]-centroid[c][k]);
				}
				// invert the result
				if (num>ZERO) num = 1.0f/num;
				else num = INF;

				mems[x-x0i+bs][y-y0i+bs][z-z0i+bs][k] = num;
				den += num;
				
				// grouping: count the number of classes with same centroids
				grouping[k] = 0.0f;
				for (int m=0;m<classes;m++)
					if (centroid[0][m]==centroid[0][k])
						grouping[k]++;
			}
			
			// normalization
			if (den>0.0f) {
				for (int k=0;k<classes;k++) {
					mems[x-x0i+bs][y-y0i+bs][z-z0i+bs][k] = 
						grouping[k]*mems[x-x0i+bs][y-y0i+bs][z-z0i+bs][k]/den;
				}
			}
			//normTime += System.currentTimeMillis()-start;				
		}
       // return mems;
    }// initialMemberships

	/**
	 *	compute the cluster centroids given the template
	 *	using template and image intensity 
	 */
    public final float[][] initialCentroids(float[][] prior, int classes) {
        float num,den,imgT;
		float Nobj;
		float[][] centroid = new float[nc][classes];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		// find image centroids
		
		if (robustMaxMin){
			for (int c=0; c<nc; c++)
				for (int k =0; k<classes; k++)
					centroid[c][k] = Ilow[c]/(Ihigh[c]-Ilow[c]) + prior[c][k];
		}else{
			for (int c=0;c<nc;c++) {
				if (cropped) {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, mix, miy, miz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, mix, miy, miz );
				} else {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, nix, niy, niz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, nix, niy, niz );
				}
				if (debug) System.out.println("image centroid range: ["+(Imin[c]+Ilow[c]*(Imax[c]-Imin[c]) )+", "
																   +(Imin[c]+Ihigh[c]*(Imax[c]-Imin[c]) )+"]\n");
				// adjust the first centroids
				for (int k=0;k<classes;k++) {
				//centroid[c][k] = Ilow + ratio[k]*(Ihigh-Ilow);
					centroid[c][k] = Ilow[c] + prior[c][k]*(Ihigh[c]-Ilow[c]);
				
				}
			}
		}
		if (debug) {
			String info = "";;
			for (int c=0;c<nc;c++) {
				if (robustMaxMin)
					info += "image "+c+" centroids: ("+(centroid[c][0]*(Ihigh[c]-Ilow[c]) );
				else
					info += "image "+c+" centroids: ("+(Imin[c]+centroid[c][0]*(Imax[c]-Imin[c]) );
					
				for (int k=1;k<classes;k++) 
					if (robustMaxMin)
						info += ", "+(centroid[c][k]*(Ihigh[c]-Ilow[c]) );
					else
						info += ", "+(Imin[c]+centroid[c][k]*(Imax[c]-Imin[c]) );
						
				info += ")\n";
			}
			
			System.out.println(info);
		}
		return centroid;
    } // initCentroids
    
	/**
	 *	compute the lesion centroid given the regular centroids
	 */
    public final float[] initialLesionCentroid(float[] prior) {
        float num,den,imgT;
		float Nobj;
		float[] lesion = new float[nc];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		if (robustMaxMin){
			for (int c=0;c<nc;c++) {
				//lesion[c] = Ilow[c] + prior[c]*(Ihigh[c]-Ilow[c]);
				lesion[c] = Ilow[c]/(Ihigh[c]-Ilow[c]) + prior[c];
			}
		}else{
			for (int c=0;c<nc;c++) {
				if (cropped) {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, mix, miy, miz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, mix, miy, miz );
				} else {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, nix, niy, niz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, nix, niy, niz );
				}
			
				if (debug) System.out.println("image centroid range: ["+(Imin[c]+Ilow[c]*(Imax[c]-Imin[c]) )+", "
																   +(Imin[c]+Ihigh[c]*(Imax[c]-Imin[c]) )+"]\n");
			}	
			for (int k=0;k<nc;k++) {
				lesion[k] = Ilow[k] + prior[k]*(Ihigh[k]-Ilow[k]);
			}
		}
		
		if (debug) {
			String info = "";;
			for (int c=0;c<nc;c++) {
				if (robustMaxMin)
					info += "image "+c+" lesion: ("+(lesion[c]*(Ihigh[c]-Ilow[c]) );
				else
					info += "image "+c+" lesion: ("+(Imin[c]+lesion[c]*(Imax[c]-Imin[c]) );
				info += ")\n";
			}
			System.out.println(info);
		}
		
		return lesion;
    } // initLesionCentroids
    
    public final float[] initialBlackHoleCentroid(float[] prior) {
        float num,den,imgT;
		float Nobj;
		float[] blackHole = new float[nc];
		
		RotationMatrix rotation = new RotationMatrix();
		rotation.setParameters(transform[0],transform[1],transform[2]);
		
		if (robustMaxMin){
			for (int c=0;c<nc;c++) {
				//lesion[c] = Ilow[c] + prior[c]*(Ihigh[c]-Ilow[c]);
				blackHole[c] = Ilow[c]/(Ihigh[c]-Ilow[c]) + prior[c];
			}
		}else{
			for (int c=0;c<nc;c++) {
				if (cropped) {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, mix, miy, miz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, mix, miy, miz );
				} else {
					Ilow[c] = ImageFunctions.robustMinimum(images[c], 0.01f, 2, nix, niy, niz );
					Ihigh[c] = ImageFunctions.robustMaximum(images[c], 0.01f, 2, nix, niy, niz );
				}
			
				if (debug) System.out.println("image centroid range: ["+(Imin[c]+Ilow[c]*(Imax[c]-Imin[c]) )+", "
																   +(Imin[c]+Ihigh[c]*(Imax[c]-Imin[c]) )+"]\n");
			}	
			for (int k=0;k<nc;k++) {
				blackHole[k] = Ilow[k] + prior[k]*(Ihigh[k]-Ilow[k]);
			}
		}
		
		if (debug) {
			String info = "";;
			for (int c=0;c<nc;c++) {
				if (robustMaxMin)
					info += "image "+c+" lesion: ("+(blackHole[c]*(Ihigh[c]-Ilow[c]) );
				else
					info += "image "+c+" lesion: ("+(Imin[c]+blackHole[c]*(Imax[c]-Imin[c]) );
				info += ")\n";
			}
			System.out.println(info);
		}
		
		return blackHole;
    } // initBlackHoleCentroids
    

}
