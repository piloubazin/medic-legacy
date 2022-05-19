package edu.jhmi.rad.medic.methods;

import java.io.*;
import java.util.*;

import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.structures.*;
import edu.jhmi.rad.medic.utilities.*;
 
/**
 *
 *  This algorithm handles registration algorithms
 *	of images with the Demons algorithm
 *	and a multiscale / multigrid technique
 *	(original Demons method, not robust)
 *
 *
 *	@version    November 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class DemonsToadsAtlasWarping {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// convenience tags
	private static final	int   X=0;
	private static final	int   Y=1;
	private static final	int   Z=2;
	private static final	int   T=3;
	
	// flag labels
	public static final	int   COMPOSITIVE 		= 101;
	public static final	int   DIFFEOMORPHIC 	= 102;
	
	public static final	int   FIXED 			= 201;
	public static final	int   MOVING 			= 202;
	public static final	int   SYMMETRIC 		= 203;
	
	public static final	int   GAUSS_FLUID 		= 301;
	public static final	int   GAUSS_DIFFUSION 	= 302;
	public static final	int   GAUSS_MIXED	 	= 303;
	public static final	int   DIV_FLUID 		= 304;
	public static final	int   DIV_DIFFUSION 	= 305;
	public static final	int   GAUSS_DIV_FLUID 		= 306;
	public static final	int   GAUSS_DIV_DIFFUSION 	= 307;
	public static final	int   ALL		 		= 401;
	public static final	int   BEST_MEM		 	= 402;
	public static final	int   BEST_ATLAS		= 403;
	public static final	int   SEGMENTED		 	= 404;
	public static final	int   SUM			 	= 405;
	public static final	int   PRODUCT		 	= 406;
	public static final	int   BEST_PRODUCT	 	= 407;
	
	// data buffers
	private 	float[][][][]		atlas;  		// source image
	private		int				na;				// number of atlas images
	private 	float[][][][]	mems;			// target image: best memberships
	private 	byte[][][][]	lb;			// target image: best labels
	private		byte[][][]		seg;			// hard segmentation (with topology)
	private		boolean[]		registered;		// objects used for registration
	private		int				nlb;				// number of best targets
	private		int				ncl;				// number of best targets
	private		int[]			group;				// groupings of classes
	private		int				ng;
	private		float[][][][]		c;				// added transform from target space to original space (multiscale)
	private		float[][][][]		s;				// added transform from target space to original space (multiscale)
	private		float[][][][]		u;				// added transform update from target space to original space (multiscale)
	private static	int		nax,nay,naz;   	// image dimensions (pyramid)
	private static	int		ntx,nty,ntz;   	// target dimensions (pyramid)
	private static	int		nsx,nsy,nsz;   	// target dimensions (pyramid)
	private static	float 	scale;
	private static	float	rax,ray,raz;   	// image resolutions (no pyramid)
	private static	float	rtx,rty,rtz;   	// target resolutions (no pyramid)
	private 		float[][]	transform;		// prior transform matrax (for instance from prior alignment)		
	private 	float[][][]		atlassum;  		// source image: sum
	private		float			atlasScale;
	private		float[]			atlasmax;
	
	// parameters
	private		float		smoothingKernel; 	// smoothing kernel size
	private		float		spatialScale; 		// scale of transformation
	private		boolean		useMask;
	private		float  		maskingThreshold;			// background threshold
	private		int		regType;
	private		int		forceType;
	private		int		fieldType;
	private		int		distType = ALL;
	
	// computation variables
	private		float			sigma2; 		// scale of transformation
	

	private		float[][]		gaussKernel;
	private		int				gx,gy,gz;
    
	private		float[][]		divKernelX, divKernelY, divKernelZ;

	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=false;
	static final boolean		verbose=false;
    
	/**
	 *  constructor
     *  note: the membership function mems_ must have sum = 0 on the masked areas
	 */
	public DemonsToadsAtlasWarping(float[][][][] img_atlas_, int na_,
									float[][][][] mems_, byte[][][][] lb_, int nlb_, int ncl_,
									byte[][][] seg_, boolean[] reged_,
									int nax_, int nay_, int naz_,
									float rax_, float ray_, float raz_,
									int ntx_, int nty_, int ntz_,
									float rtx_, float rty_, float rtz_,
									float smoothing_,
									float spScale_,
									boolean useMask_,
									float maskVal_,
									float scale_, int Ni_, int Nt_,
									int reg_, int force_, int field_,
									float atlas_scale_,
									float[][] trans_) {
	
		atlas = img_atlas_;
		na = na_;
		
		mems = mems_;
		lb = lb_;
		nlb = nlb_;
		ncl = ncl_;
		seg = seg_;
		registered = reged_;	
		scale = scale_;
		
		smoothingKernel = smoothing_;
		spatialScale = spScale_;
		sigma2 = spatialScale*spatialScale;
		
		useMask = useMask_;
		maskingThreshold = maskVal_;

		atlasScale = atlas_scale_;
		
		transform = trans_;
		
		nax = nax_;
		nay = nay_;
		naz = naz_;
		
		rax = rax_;
		ray = ray_;
		raz = raz_;
		
		ntx = ntx_;
		nty = nty_;
		ntz = ntz_;
		
		rtx = rtx_;
		rty = rty_;
		rtz = rtz_;
			
		nsx = Numerics.ceil(ntx/scale);
		nsy = Numerics.ceil(nty/scale);
		nsz = Numerics.ceil(ntz/scale);
		
		regType = reg_;
		forceType = force_;
		fieldType = field_;

		// compute atlas sum
		atlassum = new float[nax][nay][naz];
		atlasmax = new float[na];
		for (int n=0;n<na;n++) atlasmax[n] = 0.0f;
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			atlassum[x][y][z] = 0.0f;
			for (int n=0;n<na;n++) {
				atlassum[x][y][z] += atlas[n][x][y][z];
				atlasmax[n] = Numerics.max(atlasmax[n], atlas[n][x][y][z]);
			}
		}
		
		
		//initialize the smoothing kernel
		gaussKernel = ImageFunctions.separableGaussianKernel(smoothingKernel/rtx,smoothingKernel/rty,smoothingKernel/rtz);
		gx = (gaussKernel[X].length-1)/2;
		gy = (gaussKernel[Y].length-1)/2;
		gz = (gaussKernel[Z].length-1)/2;
		
		divKernelX = new float[3][];
		divKernelX[X] = new float[3];
		divKernelX[Y] = new float[1];
		divKernelX[Z] = new float[1];
		divKernelX[X][0] = 1.0f/4.0f;
		divKernelX[X][1] = 1.0f/2.0f;
		divKernelX[X][2] = 1.0f/4.0f;
		divKernelX[Y][0] = 1.0f;
		divKernelX[Z][0] = 1.0f;
				
		divKernelY = new float[3][];
		divKernelY[X] = new float[1];
		divKernelY[Y] = new float[3];
		divKernelY[Z] = new float[1];
		divKernelY[X][0] = 1.0f;
		divKernelY[Y][0] = 1.0f/4.0f;
		divKernelY[Y][1] = 1.0f/2.0f;
		divKernelY[Y][2] = 1.0f/4.0f;
		divKernelY[Z][0] = 1.0f;
				
		divKernelZ = new float[3][];
		divKernelZ[X] = new float[1];
		divKernelZ[Y] = new float[1];
		divKernelZ[Z] = new float[3];
		divKernelZ[X][0] = 1.0f;
		divKernelZ[Y][0] = 1.0f;
		divKernelZ[Z][0] = 1.0f/4.0f;
		divKernelZ[Z][1] = 1.0f/2.0f;
		divKernelZ[Z][2] = 1.0f/4.0f;

		
		isWorking = true;
		if (debug) System.out.println("Demons:initialisation ("+distType+")\n");
	}

	final public void finalize() {
		s = null; u = null; c = null;
		System.gc();
	}
    
   public final float[][][][] getCurrentTransform() {
	   return s;
   }
   public final float[][][][] getCurrentUpdate() {
		return u;
	}
    
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
	
    /** initialize the transform from previous estimate, if exists */
    public final void initializeTransform() {
				
		// initialization: start from prior transform or zero
		s = new float[3][nsx][nsy][nsz];
		if (transform==null) {
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				s[X][x][y][z] = x*scale*rtx/rax;
				s[Y][x][y][z] = y*scale*rty/ray;
				s[Z][x][y][z] = z*scale*rtz/raz;
			}
		} else {
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				s[X][x][y][z] = transform[X][X]*x*scale + transform[X][Y]*y*scale + transform[X][Z]*z*scale + transform[X][T];
				s[Y][x][y][z] = transform[Y][X]*x*scale + transform[Y][Y]*y*scale + transform[Y][Z]*z*scale + transform[Y][T];
				s[Z][x][y][z] = transform[Z][X]*x*scale + transform[Z][Y]*y*scale + transform[Z][Z]*z*scale + transform[Z][T];
			}
		}
				
		// update always zero
		u = new float[3][nsx][nsy][nsz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			u[X][x][y][z] = 0.0f;
			u[Y][x][y][z] = 0.0f;
			u[Z][x][y][z] = 0.0f;
		}
		
		// composite update always zero
		c = new float[3][nsx][nsy][nsz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			c[X][x][y][z] = 0.0f;
			c[Y][x][y][z] = 0.0f;
			c[Z][x][y][z] = 0.0f;
		}
		
		return;
    }
    
	
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public void registerImageToTarget() {
		
		if (debug) System.out.println("update: ");
		
		float meanDiff = 0.0f;
		for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
			
			// compute the update field
			float xs = s[X][x][y][z];
			float ys = s[Y][x][y][z];
			float zs = s[Z][x][y][z];
		
			float xsmx = s[X][x-1][y][z];
			float ysmx = s[Y][x-1][y][z];
			float zsmx = s[Z][x-1][y][z];
	
			float xspx = s[X][x+1][y][z];
			float yspx = s[Y][x+1][y][z];
			float zspx = s[Z][x+1][y][z];
	
			float xsmy = s[X][x][y-1][z];
			float ysmy = s[Y][x][y-1][z];
			float zsmy = s[Z][x][y-1][z];
	
			float xspy = s[X][x][y+1][z];
			float yspy = s[Y][x][y+1][z];
			float zspy = s[Z][x][y+1][z];
	
			float xsmz = s[X][x][y][z-1];
			float ysmz = s[Y][x][y][z-1];
			float zsmz = s[Z][x][y][z-1];
	
			float xspz = s[X][x][y][z+1];
			float yspz = s[Y][x][y][z+1];
			float zspz = s[Z][x][y][z+1];
			
			float xt = x*scale;
			float yt = y*scale;
			float zt = z*scale;
			
			u[X][x][y][z] = 0.0f;
			u[Y][x][y][z] = 0.0f;
			u[Z][x][y][z] = 0.0f;
			
			float den = 0.0f;
			if (lb!=null) {
				for (int n=0;n<nlb;n++) {
					int best = ImageFunctions.nearestNeighborInterpolation(lb,(byte)0,xt,yt,zt,n,ntx,nty,ntz,nlb);
					float diff = (ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb)
								- ImageFunctions.linearClosestInterpolation(atlas[best],xs,ys,zs,nax,nay,naz));
				
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[best],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[best],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[best],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[best],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[best],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[best],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += diff*diff/sigma2 + J2;
					
					u[X][x][y][z] += diff*Jx;
					u[Y][x][y][z] += diff*Jy;
					u[Z][x][y][z] += diff*Jz;
					
					meanDiff += Numerics.abs(diff);
				}
				u[X][x][y][z] /= Numerics.max(ZERO,den);
				u[Y][x][y][z] /= Numerics.max(ZERO,den);
				u[Z][x][y][z] /= Numerics.max(ZERO,den);
			} else if (distType==ALL) {
				for (int n=0;n<ncl;n++) if (registered[n]) {
					float diff = (ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb)
								- mappedAtlasInterpolation(n,xs,ys,zs) );
				
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(mappedAtlasInterpolation(n,xspx,yspx,zspx)-mappedAtlasInterpolation(n,xsmx,ysmx,zsmx));
					Jy = 0.5f/ray*(mappedAtlasInterpolation(n,xspy,yspy,zspy)-mappedAtlasInterpolation(n,xsmy,ysmy,zsmy));
					Jz = 0.5f/raz*(mappedAtlasInterpolation(n,xspz,yspz,zspz)-mappedAtlasInterpolation(n,xsmz,ysmz,zsmz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += diff*diff/sigma2 + J2;
					
					u[X][x][y][z] += diff*Jx;
					u[Y][x][y][z] += diff*Jy;
					u[Z][x][y][z] += diff*Jz;
					
					meanDiff += Numerics.abs(diff);
				}
				u[X][x][y][z] /= Numerics.max(ZERO,den);
				u[Y][x][y][z] /= Numerics.max(ZERO,den);
				u[Z][x][y][z] /= Numerics.max(ZERO,den);
			} else if (distType==PRODUCT) {
				for (int n=0;n<ncl;n++) if (registered[n]) {
					float mem = ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb);
					float atl = ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz);
				
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[n],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[n],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[n],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[n],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[n],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[n],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += mem*mem/sigma2 + J2;
					
					u[X][x][y][z] += mem*Jx;
					u[Y][x][y][z] += mem*Jy;
					u[Z][x][y][z] += mem*Jz;
					
					meanDiff += Numerics.abs(atl*mem);
				}
				u[X][x][y][z] /= Numerics.max(ZERO,den);
				u[Y][x][y][z] /= Numerics.max(ZERO,den);
				u[Z][x][y][z] /= Numerics.max(ZERO,den);
			} else if (distType==BEST_PRODUCT) {
				float mem = 0.0f;
				int nbest = -1;
				for (int n=0;n<ncl;n++) if (registered[n]) {
					float val = ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb);
					if (val>mem) {
						mem = val;
						nbest = n;
					}
				}
				if (nbest>-1) {
					
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += mem*mem/sigma2 + J2;
					
					u[X][x][y][z] += mem*Jx;
					u[Y][x][y][z] += mem*Jy;
					u[Z][x][y][z] += mem*Jz;
					
					meanDiff += Numerics.abs(mem);
					
					u[X][x][y][z] /= Numerics.max(ZERO,den);
					u[Y][x][y][z] /= Numerics.max(ZERO,den);
					u[Z][x][y][z] /= Numerics.max(ZERO,den);
				}
			} else if (distType==BEST_MEM) {
				float mem = 0.0f;
				int nbest = -1;
				for (int n=0;n<ncl;n++) if (registered[n]) {
					float val = ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb);
					if (val>mem) {
						mem = val;
						nbest = n;
					}
				}
				if (nbest>-1) {
					
					float diff = (mem - ImageFunctions.linearClosestInterpolation(atlas[nbest],xs,ys,zs,nax,nay,naz));
				
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += diff*diff/sigma2 + J2;
					
					u[X][x][y][z] += diff*Jx;
					u[Y][x][y][z] += diff*Jy;
					u[Z][x][y][z] += diff*Jz;
					
					meanDiff += Numerics.abs(diff);
					
					u[X][x][y][z] /= Numerics.max(ZERO,den);
					u[Y][x][y][z] /= Numerics.max(ZERO,den);
					u[Z][x][y][z] /= Numerics.max(ZERO,den);
				}
			} else if (distType==BEST_ATLAS) {
				float atl = 0.0f;
				int nbest = -1;
				for (int n=0;n<ncl;n++) if (registered[n]) {
					float val = ImageFunctions.linearClosestInterpolation(atlas[n],xt,yt,zt,ntx,nty,ntz);
					if (val>atl) {
						atl = val;
						nbest = n;
					}
				}
					
				if (nbest>-1) {
					float diff = (ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,nbest,ntx,nty,ntz,nlb) - atl);
					
					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += diff*diff/sigma2 + J2;
					
					u[X][x][y][z] += diff*Jx;
					u[Y][x][y][z] += diff*Jy;
					u[Z][x][y][z] += diff*Jz;
					
					meanDiff += Numerics.abs(diff);
					
					u[X][x][y][z] /= Numerics.max(ZERO,den);
					u[Y][x][y][z] /= Numerics.max(ZERO,den);
					u[Z][x][y][z] /= Numerics.max(ZERO,den);
				}
			} else if (distType==SEGMENTED) {
				byte nbest = ImageFunctions.nearestNeighborInterpolation(seg,(byte)0,xt,yt,zt,ntx,nty,ntz);

				if (registered[nbest]) {
						
					float diff = (ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,nbest,ntx,nty,ntz,nlb)
								- ImageFunctions.linearClosestInterpolation(atlas[nbest],xs,ys,zs,nax,nay,naz));

					float Jx, Jy, Jz;
					// not symmetric: too complicated
					Jx = 0.5f/rax*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspx,yspx,zspx,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmx,ysmx,zsmx,nax,nay,naz));
					Jy = 0.5f/ray*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspy,yspy,zspy,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmy,ysmy,zsmy,nax,nay,naz));
					Jz = 0.5f/raz*(ImageFunctions.linearClosestInterpolation(atlas[nbest],xspz,yspz,zspz,nax,nay,naz)-ImageFunctions.linearClosestInterpolation(atlas[nbest],xsmz,ysmz,zsmz,nax,nay,naz));
				
					float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
					den += diff*diff/sigma2 + J2;
					
					u[X][x][y][z] += diff*Jx;
					u[Y][x][y][z] += diff*Jy;
					u[Z][x][y][z] += diff*Jz;
					
					meanDiff += Numerics.abs(diff);
					
					u[X][x][y][z] /= Numerics.max(ZERO,den);
					u[Y][x][y][z] /= Numerics.max(ZERO,den);
					u[Z][x][y][z] /= Numerics.max(ZERO,den);
				}
			} else if (distType==SUM) {
				float mem = 0.0f, atl = 0.0f;
				float Jpx = 0.0f, Jpy = 0.0f, Jpz = 0.0f;
				float Jmx = 0.0f, Jmy = 0.0f, Jmz = 0.0f;
				float sAtl = 0.0f;
				float sJpx = 0.0f, sJpy = 0.0f, sJpz = 0.0f;
				float sJmx = 0.0f, sJmy = 0.0f, sJmz = 0.0f;
				
				for (int n=0;n<ncl;n++) {
					sAtl += ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz);
					
					sJpx += ImageFunctions.linearClosestInterpolation(atlas[n],xspx,yspx,zspx,nax,nay,naz);
					sJpy += ImageFunctions.linearClosestInterpolation(atlas[n],xspy,yspy,zspy,nax,nay,naz);
					sJpz += ImageFunctions.linearClosestInterpolation(atlas[n],xspz,yspz,zspz,nax,nay,naz);
				
					sJmx += ImageFunctions.linearClosestInterpolation(atlas[n],xsmx,ysmx,zsmx,nax,nay,naz);
					sJmy += ImageFunctions.linearClosestInterpolation(atlas[n],xsmy,ysmy,zsmy,nax,nay,naz);
					sJmz += ImageFunctions.linearClosestInterpolation(atlas[n],xsmz,ysmz,zsmz,nax,nay,naz);
					
					if (registered[n]) {
						mem += ImageFunctions.linearClosestInterpolation(mems,xt,yt,zt,n,ntx,nty,ntz,nlb);
						atl += ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz);
				
						Jpx += ImageFunctions.linearClosestInterpolation(atlas[n],xspx,yspx,zspx,nax,nay,naz);
						Jpy += ImageFunctions.linearClosestInterpolation(atlas[n],xspy,yspy,zspy,nax,nay,naz);
						Jpz += ImageFunctions.linearClosestInterpolation(atlas[n],xspz,yspz,zspz,nax,nay,naz);
					
						Jmx += ImageFunctions.linearClosestInterpolation(atlas[n],xsmx,ysmx,zsmx,nax,nay,naz);
						Jmy += ImageFunctions.linearClosestInterpolation(atlas[n],xsmy,ysmy,zsmy,nax,nay,naz);
						Jmz += ImageFunctions.linearClosestInterpolation(atlas[n],xsmz,ysmz,zsmz,nax,nay,naz);
					}
				}
				// rescale the atlas
				atl /= Numerics.max(ZERO,sAtl);
				Jpx /= Numerics.max(ZERO,sJpx);
				Jpy /= Numerics.max(ZERO,sJpy);
				Jpz /= Numerics.max(ZERO,sJpz);
				Jmx /= Numerics.max(ZERO,sJmx);
				Jmy /= Numerics.max(ZERO,sJmy);
				Jmz /= Numerics.max(ZERO,sJmz);
				
				float diff = mem-atl;
				
				float Jx = 0.5f/rax*(Jpx-Jmx);
				float Jy = 0.5f/ray*(Jpy-Jmy);
				float Jz = 0.5f/raz*(Jpz-Jmz);
				
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				den += diff*diff/sigma2 + J2;
					
				u[X][x][y][z] += diff*Jx;
				u[Y][x][y][z] += diff*Jy;
				u[Z][x][y][z] += diff*Jz;
					
				meanDiff += Numerics.abs(diff);
				
				u[X][x][y][z] /= Numerics.max(ZERO,den);
				u[Y][x][y][z] /= Numerics.max(ZERO,den);
				u[Z][x][y][z] /= Numerics.max(ZERO,den);
			}
		}
		meanDiff /= (nsx*nsy*nsz);
		
		if (regType==GAUSS_FLUID || regType==GAUSS_MIXED || regType==GAUSS_DIV_FLUID) {
			if (debug) System.out.println("GAUSS_FLUID regularization \n");
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFunctions.separableConvolution(u[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Y] = ImageFunctions.separableConvolution(u[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Z] = ImageFunctions.separableConvolution(u[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		} else if (regType==DIV_FLUID || regType==GAUSS_DIV_FLUID) {
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFunctions.separableConvolution(u[X],nsx,nsy,nsz,divKernelX,1,0,0);
			u[Y] = ImageFunctions.separableConvolution(u[Y],nsx,nsy,nsz,divKernelY,0,1,0);
			u[Z] = ImageFunctions.separableConvolution(u[Z],nsx,nsy,nsz,divKernelZ,0,0,1);
		}
		
		// compose the transformations
		if (fieldType==COMPOSITIVE) {
			if (debug) System.out.println("compose with current transform \n");
							
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				float xu = x+u[X][x][y][z];
				float yu = y+u[Y][x][y][z];
				float zu = z+u[Z][x][y][z];
				
				// note: if outside, extrapolate as X+u
				c[X][x][y][z] = ImageFunctions.linearClosestInterpolation(s[X],xu,yu,zu,nsx,nsy,nsz) - x*scale;
				c[Y][x][y][z] = ImageFunctions.linearClosestInterpolation(s[Y],xu,yu,zu,nsx,nsy,nsz) - y*scale;
				c[Z][x][y][z] = ImageFunctions.linearClosestInterpolation(s[Z],xu,yu,zu,nsx,nsy,nsz) - z*scale;
			}
		}
		
		if (regType==GAUSS_DIFFUSION || regType==GAUSS_MIXED || regType==GAUSS_DIV_DIFFUSION) {
			if (debug) System.out.println("GAUSS_DIFFUSION regularization \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFunctions.separableConvolution(c[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Y] = ImageFunctions.separableConvolution(c[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Z] = ImageFunctions.separableConvolution(c[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		} else if (regType==DIV_DIFFUSION || regType==GAUSS_DIV_DIFFUSION) {
			if (debug) System.out.println("DIV_DIFFUSION regularization \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFunctions.separableConvolution(c[X],nsx,nsy,nsz,divKernelX,1,0,0);
			c[Y] = ImageFunctions.separableConvolution(c[Y],nsx,nsy,nsz,divKernelY,0,1,0);
			c[Z] = ImageFunctions.separableConvolution(c[Z],nsx,nsy,nsz,divKernelZ,0,0,1);
		}
					
		float meanC = 0.0f;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			s[X][x][y][z] = x*scale + c[X][x][y][z];
			s[Y][x][y][z] = y*scale + c[Y][x][y][z];
			s[Z][x][y][z] = z*scale + c[Z][x][y][z];
			
			meanC += c[X][x][y][z]*c[X][x][y][z]+c[Y][x][y][z]*c[Y][x][y][z]+c[Z][x][y][z]*c[Z][x][y][z];
		}
		meanC /= (nsx*nsy*nsz);
		
		if (debug) System.out.println("convergence "+meanC+" -> "+meanDiff+"\n");

        return;
    } // 
    
    /** 
	 *	returns the transformed coordinates
	 */
	public final float[] getCurrentMapping(int x, int y, int z) {
		float xs = x/scale;
		float ys = y/scale;
		float zs = z/scale;
		
		float[] Xs = new float[3];
		Xs[X] = ImageFunctions.linearClosestInterpolation(s[X],xs,ys,zs,nsx,nsy,nsz);
		Xs[Y] = ImageFunctions.linearClosestInterpolation(s[Y],xs,ys,zs,nsx,nsy,nsz);
		Xs[Z] = ImageFunctions.linearClosestInterpolation(s[Z],xs,ys,zs,nsx,nsy,nsz);
		return Xs;
	}// getCurrentMapping

    /** 
	 *	returns the transformed image
	 */
	public final float[][][][] exportTransformedAtlas() {
		float 	xs,ys,zs;
		float[][][][]	atl;
		
		atl = new float[na][ntx][nty][ntz];
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			xs = ImageFunctions.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageFunctions.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageFunctions.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			// compute interpolated values
			for (int n=0;n<na;n++) {
				atl[n][x][y][z] = mappedAtlasInterpolation(n,xs,ys,zs);
				//atl[n][x][y][z] = ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz);
			}
		}
		return atl;
	} // exportTransformedImage

    /** 
	 *	returns the transform field v defined as X' = X+v (=s)
	 */
	public final float[][][][] exportTransformField() {
		float 	xs,ys,zs;
		float[][][][]	vec = new float[3][ntx][nty][ntz];
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			xs = ImageFunctions.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageFunctions.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageFunctions.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			vec[X][x][y][z] = xs-x;
			vec[Y][x][y][z] = ys-y;
			vec[Z][x][y][z] = zs-z;			
 		}
		return vec;
	} // exportTransformedImage

	/**
	 *	linear interpolation, with value outside the image
	 */
	/* 
	private float linearGroupedAtlasInterpolation(float x, float y, float z, int g) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nax-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,nay-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,naz-2);
		
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = -1e30f;
		for (int n=0;n<ncl;n++) if (group[n]==g) {
			val = Numerics.max(val,
				  nalpha*nbeta*ngamma*atlas[n][x0][y0][z0]
				+ alpha*nbeta*ngamma*atlas[n][x0+1][y0][z0]
				+ nalpha*beta*ngamma*atlas[n][x0][y0+1][z0]
				+ nalpha*nbeta*gamma*atlas[n][x0][y0][z0+1] 
				+ alpha*beta*ngamma*atlas[n][x0+1][y0+1][z0] 
				+ nalpha*beta*gamma*atlas[n][x0][y0+1][z0+1] 
				+ alpha*nbeta*gamma*atlas[n][x0+1][y0][z0+1] 
				+ alpha*beta*gamma*atlas[n][x0+1][y0+1][z0+1] );
		}
		return val;
	}
	*/

	/*	
	private final float mappedAtlas(float val) {
		return Numerics.bounded(5.0f*(val-0.5f)+0.5f,0.0f,1.0f);
	}
	*/
	
	private final float mappedAtlasInterpolation(int n, float xs, float ys, float zs) {
		/*
		return ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz)
				/Numerics.max(ZERO,ImageFunctions.linearClosestInterpolation(atlassum,xs,ys,zs,nax,nay,naz));
		*/
		/*
		return Numerics.bounded(atlasScale*(ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f)
				/Numerics.bounded(atlasScale*(atlasmax[n]-0.5f)+0.5f,0.0f,1.0f);
		*/
		//return Numerics.bounded(atlasScale*(ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz)/atlasmax[n]-0.5f)+0.5f,0.0f,1.0f);
		return Numerics.bounded(atlasScale*(ImageFunctions.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
	}
}
