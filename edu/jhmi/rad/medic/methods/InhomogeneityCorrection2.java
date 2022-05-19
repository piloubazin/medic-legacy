package edu.jhmi.rad.medic.methods;

import java.io.*;
import java.util.*;

import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.utilities.*;
 
import Jama.Matrix;

/**
 *
 *  This algorithm handles the inhomogeneity correction
 *	for segmentation algorithms like FANTASM and TOADS.
 *	<p>
 *	The method is based on estimating the parameters of a low degree polynomial.
 *
 *
 *	@version    December 2004
 *	@author     Pierre-Louis Bazin
 *	@see		SegmentationFCM
 *	@see		FuzzyToads
 *	@see		LightMarchingSegmentation
 *
 */
 
public class InhomogeneityCorrection2 {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// data buffers
	private 	float[][][]			image;  			// original image
	private 	float[][][][]		mems;				// membership function
	private 	byte[][][][]		best;				// membership function labels
	private 	float[]				centroids;			// cluster centroids
	private 	float[][][]			field;  			// inhomogeneity field
	private 	float[][][][]		fields;  			// separate inhomogeneity fields
	private 	float[][]			transform = null;	// transformation matrix
	private 	boolean[][][]		mask;   			// image mask: true for data points
	private 	byte[][][]			segmentation;   	// image segmentation: for discarding masked objects, if used
	private		static	int			nx,ny,nz;   		// image dimensions
	private		static	float		rx,ry,rz;   		// image resolutions
	private		static	int			dimensions;			// image dimensions (2D or 3D)
	
	// parameters
	private 	int 		clusters;   	// number of clusters
	private 	int 		classes;    	// number of classes in original membership: > clusters if outliers
	private 	int 		approx;    		// number of classes in approximation
    private		int			degree;			// polynomial function degree
	private		float		scale;		// adding a smoothing term so that local variations are small

	private		byte[]		templateLabel;	// the label for each class, if used
	private		String[]	objType;		// the type of object for each class, if used
	
	// variable metrics
	private		float	pcsf = 1.0f;	
	private		float	pgm = 1.0f;	
	private		float	pwm = 1.0f;	
	private		PowerTable	powercsf,powergm,powerwm;
				
	// computation variables
	private float[]			pol;		// the array for the polynomial basis
	private	int				Np;			// the number of polynomial coefficients
	
	private	int				subsample = 3;	// the sub-sampling to speed up the computations
	
	// type of correction: image field, centroid field or separate centroid fields
	private	int					factorType;	
	public static final	int		IMAGE = 1;
	public static final	int   	CENTROIDS = 2;
	private static final	int   	OPTIMIZED_IMAGE = 3;
	private static final	int   	OPTIMIZED_CENTROIDS = 4;
	
	private int					backgroundType;
	public static final int		OBJECT_BASED = 4;
	public static final int		CLASS_BASED = 5;
	private int					background;			// the index of the background region, if any
	
	private int					q;
	public static final int		Q = -1;
	private float				pq;
	private PowerTable			powerq;
	
	private int					normType;
	public static final int		EQUAL_VARIANCE=1;
	public static final int		SEPARATE_VARIANCES=2;
	private float[]				var;
	
	private int					fieldType;
	public static final int		GLOBAL = 1;
	public static final int		SEPARATE = 2;
	
	private int					algorithmType;
	public static final int		CHEBYSHEV=1;
	public static final int		SPLINES=2;
	
	private int					coordinateType;
	public static final int		FIXED=1;
	public static final int		TRANSFORMED=2;
	
	private int					memType;
	public static final int		EXACT=1;
	public static final int		APPROX=2;
		
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=false;
	static final boolean		verbose=true;
	
	/**
	 *  constructor for Fantasm
	 *	note: all images passed to the algorithm are just linked, not copied
	 */
	public InhomogeneityCorrection2(float[][][] image_, float[][][][] mems_, float[] cent_,
									int classes_, int clusters_,
									boolean[][][] mask_,
									int deg_, float scale_,
									int fact_, int field_, int algo_, int coord_, int norm_,
									float q_, float[] opt_, float[] var_,
									int nx_, int ny_, int nz_,
									float rx_, float ry_, float rz_) {
		
		image = image_;
		mask = mask_;
		mems = mems_;
		centroids = cent_;
		
		degree = deg_;
		if (degree==1) Np = 4;
		else if (degree==2) Np = 10;
		else if (degree==3) Np = 20;
		else if (degree==4) Np = 35;
		else {
			isWorking = false;
			return;
		}		
		
		scale = scale_;
		
		factorType = fact_;
		fieldType = field_;
		algorithmType = algo_;
		coordinateType = coord_;
		normType = norm_;
		var = var_;
		
		powercsf = null;
		powergm = null;
		powerwm = null;
		
		pq = q_;
		if (pq!=0 && pq!=1 && pq!=2 && pq!=3) {
			powerq = new PowerTable(0.0f , 1.0f , 0.000001f , pq );
			q = Q;
		} else {
			powerq = null;
			q = Numerics.round(pq);
		}
		
		backgroundType = CLASS_BASED;
		background = -1;
		templateLabel = null;
		objType = null;
		segmentation = null;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		classes = classes_;
		clusters = clusters_;
		approx = 0;
		memType = EXACT;
		best = null;
		
		// init all the new arrays
		try {
			if (fieldType==SEPARATE) {
				fields = new float[nx][ny][nz][clusters];
				field = null;
			} else {
				field = new float[nx][ny][nz];
				fields = null;
			}
			if (algorithmType==CHEBYSHEV) pol = new float[Np];
			else pol = null;
			
			if (coordinateType==TRANSFORMED) transform = new float[3][4];
			else transform = null;
			
		} catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;

		// init values
		if (fieldType==SEPARATE) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int k=0;k<clusters;k++) {
				fields[x][y][z][k] = 1.0f;
			}
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				field[x][y][z] = 1.0f;
			}
		}	
		
		if (coordinateType==TRANSFORMED) {
			for (int i=0;i<3;i++) {
				for (int j=0;j<4;j++)
					transform[i][j] = 0.0f;
				transform[i][i] = 1.0f;
			}
		}
		if (debug) System.out.println("IC:initialisation\n");
	}

	/**
	 *  constructor for Toads
	 *	note: all images passed to the algorithm are just linked, not copied
	 */
	public InhomogeneityCorrection2(float[][][] image_, float[][][][] mems_, float[] cent_,
									int clusters_,
									byte[][][] classif_, String[] type_,
									int deg_, float scale_,
									int fact_, int field_, int algo_, int coord_, int norm_,
									float q_, float[] opt_, float[] var_,
									int nx_, int ny_, int nz_,
									float rx_, float ry_, float rz_) {
		
		image = image_;
		mask = null;
		mems = mems_;
		centroids = cent_;
		best = null;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		classes = clusters_;
		clusters = clusters_;
		approx = 0;
		memType = EXACT;
		
		degree = deg_;
		if (degree==1) Np = 4;
		else if (degree==2) Np = 10;
		else if (degree==3) Np = 20;
		else if (degree==4) Np = 35;
		else {
			isWorking = false;
			return;
		}		
		
		scale = scale_;
		
		factorType = fact_;
		fieldType = field_;
		algorithmType = algo_;
		coordinateType = coord_;
		normType = norm_;
		var = var_;
		
		pq = q_;
		if (pq!=0 && pq!=1 && pq!=2 && pq!=3) {
			powerq = new PowerTable(0.0f , 1.0f , 0.000001f , pq );
			q = Q;
		} else {
			powerq = null;
			q = Numerics.round(pq);
		}
		
		backgroundType = OBJECT_BASED;
		objType = type_;
		segmentation = classif_;
		background = -1;
		for (int k=0;k<classes;k++) if (objType[k].contains("mask")) {
			background = k;
		}
		if (debug) System.out.println("Background label: "+background);
		
		// init all the new arrays
		try {
			if (fieldType==SEPARATE) {
				fields = new float[nx][ny][nz][clusters];
				field = null;
			} else {
				field = new float[nx][ny][nz];
				fields = null;
			}
			if (algorithmType==CHEBYSHEV) pol = new float[Np];
			else pol = null;
			
			mask = new boolean[nx][ny][nz];
			
			if (coordinateType==TRANSFORMED) transform = new float[3][4];
			else transform = null;

			for (int k=0;k<clusters;k++) if (objType[k].startsWith("optimized_")) {
				if (factorType==IMAGE) factorType = OPTIMIZED_IMAGE;
				if (factorType==CENTROIDS) factorType = OPTIMIZED_CENTROIDS;
				
				if (opt_!=null && opt_.length==3) {
					pcsf = opt_[0];
					pgm = opt_[1];
					pwm = opt_[2];
				} else {
					pcsf = 1.0f;
					pgm = 1.0f;
					pwm = 1.0f;
				}
				powercsf = new PowerTable(0.0f , 2.0f , 0.000001f , pcsf );
				powergm = new PowerTable(0.0f , 2.0f , 0.000001f , pgm );
				powerwm = new PowerTable(0.0f , 2.0f , 0.000001f , pwm );
			}
			
		} catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;

		// init values
		if (fieldType==SEPARATE) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int k=0;k<clusters;k++) {
				fields[x][y][z][k] = 1.0f;
			}
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				field[x][y][z] = 1.0f;
			}
		}	
		
		if (coordinateType==TRANSFORMED) {
			for (int i=0;i<3;i++) {
				for (int j=0;j<4;j++)
					transform[i][j] = 0.0f;
				transform[i][i] = 1.0f;
			}
		}
		
		if (debug) System.out.println("IC:initialisation\n");
	}

	/**
	 *  constructor for Toads with approximations
	 *	note: all images passed to the algorithm are just linked, not copied
	 */
	public InhomogeneityCorrection2(float[][][] image_, float[][][][] mems_, byte[][][][] best_, float[] cent_,
									int clusters_, int approx_,
									byte[][][] classif_, String[] type_,
									int deg_, float scale_,
									int fact_, int field_, int algo_, int coord_, int norm_,
									float q_, float[] opt_, float[] var_,
									int nx_, int ny_, int nz_,
									float rx_, float ry_, float rz_) {
		
		image = image_;
		mask = null;
		mems = mems_;
		best = best_;
		centroids = cent_;
		
		memType = APPROX;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		classes = clusters_;
		clusters = clusters_;
		approx = approx_;
		
		degree = deg_;
		if (degree==1) Np = 4;
		else if (degree==2) Np = 10;
		else if (degree==3) Np = 20;
		else if (degree==4) Np = 35;
		else {
			isWorking = false;
			return;
		}		
		
		scale = scale_;
		
		factorType = fact_;
		fieldType = field_;
		algorithmType = algo_;
		coordinateType = coord_;
		normType = norm_;
		var = var_;
		
		pq = q_;
		if (pq!=0 && pq!=1 && pq!=2 && pq!=3) {
			powerq = new PowerTable(0.0f , 1.0f , 0.000001f , pq );
			q = Q;
		} else {
			powerq = null;
			q = Numerics.round(pq);
		}
		
		backgroundType = OBJECT_BASED;
		objType = type_;
		segmentation = classif_;
		background = -1;
		for (int k=0;k<classes;k++) if (objType[k].contains("mask") ) {
			background = k;
		}
		if (debug) System.out.println("background label = "+background);
		
		// init all the new arrays
		try {
			if (fieldType==SEPARATE) {
				fields = new float[nx][ny][nz][clusters];
				field = null;
			} else {
				field = new float[nx][ny][nz];
				fields = null;
			}
			if (algorithmType==CHEBYSHEV) pol = new float[Np];
			else pol = null;
			
			mask = new boolean[nx][ny][nz];
			
			if (coordinateType==TRANSFORMED) transform = new float[3][4];
			else transform = null;

			for (int k=0;k<clusters;k++) if (objType[k].startsWith("optimized_")) {
				if (factorType==IMAGE) factorType = OPTIMIZED_IMAGE;
				if (factorType==CENTROIDS) factorType = OPTIMIZED_CENTROIDS;
				
				if (opt_!=null && opt_.length==3) {
					pcsf = opt_[0];
					pgm = opt_[1];
					pwm = opt_[2];
				} else {
					pcsf = 1.0f;
					pgm = 1.0f;
					pwm = 1.0f;
				}
				powercsf = new PowerTable(0.0f , 2.0f , 0.000001f , pcsf );
				powergm = new PowerTable(0.0f , 2.0f , 0.000001f , pgm );
				powerwm = new PowerTable(0.0f , 2.0f , 0.000001f , pwm );
			}
			
		} catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;

		// init values
		if (fieldType==SEPARATE) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int k=0;k<clusters;k++) {
				fields[x][y][z][k] = 1.0f;
			}
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				field[x][y][z] = 1.0f;
			}
		}	
		
		if (coordinateType==TRANSFORMED) {
			for (int i=0;i<3;i++) {
				for (int j=0;j<4;j++)
					transform[i][j] = 0.0f;
				transform[i][i] = 1.0f;
			}
		}
		
		if (debug) System.out.println("IC:initialisation\n");
	}

	public InhomogeneityCorrection2(float[][][] image_, float[][][][] mems_, float[] cent_,
									int clusters_,
									byte[][][] classif_, byte[] tempLabel_, String[] type_,
									int deg_, 
									int nx_, int ny_, int nz_,
									float rx_, float ry_, float rz_) {
		this(image_,mems_,cent_,clusters_,classif_,type_,
									deg_,0.0f,IMAGE,GLOBAL,CHEBYSHEV,FIXED,EQUAL_VARIANCE,
									2.0f,null,null,
									nx_,ny_,nz_,rx_,ry_,rz_);
	}
		
	public InhomogeneityCorrection2(float[][][] image_, float[][][][] mems_, float[] cent_,
									int classes_, int clusters_,
									boolean[][][] mask_,
									int deg_, 
									int nx_, int ny_, int nz_,
									float rx_, float ry_, float rz_) {
		this(image_,mems_,cent_,classes_,clusters_,mask_,
				deg_,0.0f,IMAGE,GLOBAL,CHEBYSHEV,FIXED,EQUAL_VARIANCE,
				2.0f,null,null,
				nx_,ny_,nz_,rx_,ry_,rz_);
	}
		
		
	/** clean-up: destroy membership and centroid arrays */
	public final void finalize() {
		field = null;
		pol = null;
		transform = null;
		System.gc();
	}
	
    /** accessor for computed data */ 
    public final float[][][] getField() { return field; }
	public final float[][][][] getFields() { return fields; }
	/** accessor for computed data */
	public final void importTransform(float[][] trans_) { 
		for (int i=0;i<3;i++) for (int j=0;j<4;j++) transform[i][j] = trans_[i][j]; 
	}
	/** accessor for computed data */
	public final void setMemberships(float[][][][] mem) { mems = mem; }
	/** accessor for computed data */
	public final void setMask(boolean[][][] msk) { mask = msk; }
	/** accessor for computed data */
	public final void setImage(float[][][] img) { image = img; }
    /** computation flags */
	public final boolean isWorking() { return isWorking; }
	/** computation flags */
	public final boolean isCompleted() { return isCompleted; }
	
	/** 
	 *	export field 
	 */
	public final float[][][] exportField() {
		int 	x,y,z;
		float[][][]	Field = new float[nx][ny][nz];
		
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			Field[x][y][z] = field[x][y][z];
		}
		return Field;
	} // exportField

	/** 
	 *	export fields 
	 */
	public final float[][][][] exportFields() {
		int 	x,y,z;
		float[][][][]	Fields = new float[nx][ny][nz][clusters];
		
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (int k=0;k<clusters;k++) {
			Fields[x][y][z][k] = fields[x][y][z][k];
		}
		return Fields;
	} // exportFields

	/** 
	 *	export field for a single class 
	 */
	public final float[][][] exportFields(int k) {
		int 	x,y,z;
		float[][][]	Fields = new float[nx][ny][nz];
		
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			Fields[x][y][z] = fields[x][y][z][k];
		}
		return Fields;
	} // exportFields

    /** 
	 *  compute the inhomogeneity field using a Chebyshev polynomial basis
	 *	(main function to call)
	 */
	 final public void computeCorrectionField(int deg) {
		 switch (algorithmType) {
			 case CHEBYSHEV: 	compute3DPolynomialField(deg);	break;
			 case SPLINES:	 	compute3DSplineField( (degree+1)/(deg+1)*scale); break;
			 default:			return;
		 }		 
	 }
    
     /** 
	 *  compute the inhomogeneity field using a Chebyshev polynomial basis
	 *	(main function to call)
	 */
	 final public void computeCorrectionFields(int deg, float scale) {
		 switch (algorithmType) {
			 case CHEBYSHEV: 	compute3DSeparatePolynomialFields(deg);	break;
			 //case SPLINES:	 	compute3DSeparateSplineFields(deg*scale/degree);break;
			 default:			return;
		 }		 
	 }
    
    /** 
	 *  compute the inhomogeneity field using a smoothing spline basis
	 *  with adequate normalization
	 */
    final private void compute3DSplineField(float scale) {
		int progress = 0;
        int mod = nx*ny*nz/100; // mod is 1 percent of length

        long inner_loop_time;
		if (verbose) inner_loop_time = System.currentTimeMillis();

		// compute the corresponding scales
		if (scale==0) return;
		float lambdaX = 0.5f*(scale/rx)*(scale/rx)*(scale/rx)*(scale/rx);
		float lambdaY = 0.5f*(scale/ry)*(scale/ry)*(scale/ry)*(scale/ry);
		float lambdaZ = 0.5f*(scale/rz)*(scale/rz)*(scale/rz)*(scale/rz);

		// center on zero
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			field[x][y][z] += -1.0f;
		}
		
		// update mask, if needed
		computeBackgroundMask();
		
		// compute new samples: only on a subset of data
		for (int x=0;x<nx;x+=subsample) for (int y=0;y<ny;y+=subsample) for (int z=0;z<nz;z+=subsample) {
			if (mask[x][y][z]) {
				// valid image point
				float num = 0.0f;
				float den = 0.0f;
				for (int k=0;k<clusters;k++) if (k!=background) {
					num += splineVectorTerm(x,y,z,k);
					den += splineMatrixTerm(x,y,z,k);
				}
				field[x][y][z] = num/Numerics.max(ZERO,den);
			} else {
				field[x][y][z] = 0.0f;
			}
		}
		// spline smoothing
		SplineProcessing algorithm = new SplineProcessing();
		algorithm.samplesToSmoothingCoefficients(field,nx,ny,nz,lambdaX,lambdaY,lambdaZ);
		
		// final field
		float avg = algorithm.splineBasedImage(field,nx,ny,nz,3);
								
		// re-normalize
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			field[x][y][z] += 1.0f - avg;
		}
		
        if (debug) System.out.print("inner loop time: (milliseconds): " + (System.currentTimeMillis()-inner_loop_time) +"\n"); 

        return;
    } // computeNormalizedPolynomialField
    
    final private void compute3DPolynomialField(int deg) {
		switch (memType) {
			case EXACT:		computeExact3DPolynomialField(deg); break;
			case APPROX:	computeApprox3DPolynomialField(deg); break;
			default:		computeExact3DPolynomialField(deg); break;
		}
	}
	
	/** 
	 *  compute the inhomogeneity field using a Chebyshev polynomial basis
	 *  with adequate normalization and options (generic)
	 */
    final private void computeExact3DPolynomialField(int deg) {
        int progress = 0;
        int mod = nx*ny*nz/100; // mod is 1 percent of length

        long inner_loop_time;
		if (verbose) inner_loop_time = System.currentTimeMillis();

		// compute the degree and size
		int np = 0;
		if (deg>degree) {
			deg = degree;
			np = Np;
		} else 
		if (deg==0) return;   else 
		if (deg==1) np = 4;   else 
		if (deg==2) np = 10;  else 
		if (deg==3) np = 20;  else 
		if (deg>=4) np = 35;
			
		// factor = pol*mems*centroids*pol; img = pol*image
		Matrix img = new Matrix(np,1,0.0f);
        Matrix factor = new Matrix(np,np,0.0f);
	
		progress = 0;
        mod = nx*ny*nz/10; // mod is 1 percent of length

		computeBackgroundMask();
		
		for (int x=0;x<nx;x+=subsample) for (int y=0;y<ny;y+=subsample) for (int z=0;z<nz;z+=subsample) {
			if (mask[x][y][z]) {
				// valid image point
				compute3DChebyshevCoefficients(x,y,z,deg);		
				for (int k=0;k<clusters;k++) if (k!=background) {
					float vect = polynomialVectorTerm(x,y,z,k,0);
					float matr = polynomialMatrixTerm(x,y,z,k,0);
					
					for (int n=0;n<np;n++) {
						img.set(n,0, img.get(n,0) + vect*pol[n] );
						for (int m=0;m<np;m++) {
							factor.set(n,m, factor.get(n,m) + matr*pol[n]*pol[m] );
						}
					}
				}	
			}
		}
		for (int n=0;n<np;n++) {
			img.set(n,0, img.get(n,0)/(float)(nx*ny*nz) );
			for (int m=0;m<np;m++)
				factor.set(n,m, factor.get(n,m)/(float)(nx*ny*nz) );
		}
		
		if (debug) {
			//System.out.println("factor: "+factor.matrixToString(5,2)+"\n");
			//System.out.println("img: "+img.matrixToString(5,2)+"\n");
		}
		
		// coeff = factor^-1 img
		Matrix coeff = factor.solve(img);
        
		//if (debug) System.out.println("coeff: "+coeff.matrixToString(5,2)+"\n");
		
		// transfer to the field
		float min = INF;
		float max = -INF;
		float diff = 0.0f;
		
		float mean = 0.0f;
		float den = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (mask[x][y][z]) {
				float val = 0.0f;
			
				compute3DChebyshevCoefficients(x,y,z,deg);		
				for (int n=0;n<np;n++) {
					val += coeff.get(n,0)*pol[n];
				}
				mean += val;
				den++;
				
				if (val < min) min = val;
				else if (val > max) max = val;
				diff = Numerics.abs(field[x][y][z]-val);	
				
				field[x][y][z] = val;	
			}
		}
		mean = mean/den;
		// normalize and fill the masked area with 1 rather than 0
		// in case of moving masks
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (mask[x][y][z]) {
				field[x][y][z] += 1.0f - mean;
				//field[x][y][z] /= mean;
			}
		}
		
		// output
		if (verbose) {
			System.out.println("IC: degree "+deg+", <min|avg|max> : <"+min+"|"+mean+"|"+max+">\n");
		}
		
		// clean-up
		factor = null;
		img = null;
		coeff = null;
		
        if (debug) System.out.print("inner loop time: (milliseconds): " + (System.currentTimeMillis()-inner_loop_time) +"\n"); 

        return;
    } // computeExact3DPolynomialField
    
	/** 
	 *  compute the inhomogeneity field using a Chebyshev polynomial basis
	 *  with adequate normalization and options (generic).
	 *  This uses approximated memberships.
	 */
    final private void computeApprox3DPolynomialField(int deg) {
        int progress = 0;
        int mod = nx*ny*nz/100; // mod is 1 percent of length

        long inner_loop_time;
		if (verbose) inner_loop_time = System.currentTimeMillis();

		// compute the degree and size
		int np = 0;
		if (deg>degree) {
			deg = degree;
			np = Np;
		} else 
		if (deg==0) return;   else 
		if (deg==1) np = 4;   else 
		if (deg==2) np = 10;  else 
		if (deg==3) np = 20;  else 
		if (deg>=4) np = 35;
			
		// factor = pol*mems*centroids*pol; img = pol*image
		Matrix img = new Matrix(np,1,0.0f);
        Matrix factor = new Matrix(np,np,0.0f);
	
		progress = 0;
        mod = nx*ny*nz/10; // mod is 1 percent of length

		computeBackgroundMask();
		
		for (int x=0;x<nx;x+=subsample) for (int y=0;y<ny;y+=subsample) for (int z=0;z<nz;z+=subsample) {
			if (mask[x][y][z]) {
				// valid image point
				compute3DChebyshevCoefficients(x,y,z,deg);		
				for (int a=0;a<approx;a++) {
					int k = best[x][y][z][a];
					if (k!=background) {
						float vect = polynomialVectorTerm(x,y,z,k,a);
						float matr = polynomialMatrixTerm(x,y,z,k,a);
						
						for (int n=0;n<np;n++) {
							img.set(n,0, img.get(n,0) + vect*pol[n] );
							for (int m=0;m<np;m++) {
								factor.set(n,m, factor.get(n,m) + matr*pol[n]*pol[m] );
							}
						}
					}
				}	
			}
		}
		for (int n=0;n<np;n++) {
			img.set(n,0, img.get(n,0)/(float)(nx*ny*nz) );
			for (int m=0;m<np;m++)
				factor.set(n,m, factor.get(n,m)/(float)(nx*ny*nz) );
		}
		
		if (debug) {
			//System.out.println("factor: "+factor.matrixToString(5,2)+"\n");
			//System.out.println("img: "+img.matrixToString(5,2)+"\n");
		}
		
		// coeff = factor^-1 img
		Matrix coeff = factor.solve(img);
        
		//if (debug) System.out.println("coeff: "+coeff.matrixToString(5,2)+"\n");
		
		// transfer to the field
		float min = INF;
		float max = -INF;
		float diff = 0.0f;
		
		float mean = 0.0f;
		float den = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (mask[x][y][z]) {
				float val = 0.0f;
			
				compute3DChebyshevCoefficients(x,y,z,deg);		
				for (int n=0;n<np;n++) {
					val += coeff.get(n,0)*pol[n];
				}
				mean += val;
				den++;
				
				if (val < min) min = val;
				else if (val > max) max = val;
				diff = Numerics.abs(field[x][y][z]-val);	
				
				field[x][y][z] = val;	
			}
		}
		mean = mean/den;
		// normalize and fill the masked area with 1 rather than 0
		// in case of moving masks
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (mask[x][y][z]) {
				field[x][y][z] += 1.0f - mean;
				//field[x][y][z] /= mean;
			}
		}
		
		// output
		if (verbose) {
			System.out.println("IC: degree "+deg+", <min|avg|max> : <"+min+"|"+mean+"|"+max+">\n");
		}
		
		// clean-up
		factor = null;
		img = null;
		coeff = null;
		
        if (debug) System.out.print("inner loop time: (milliseconds): " + (System.currentTimeMillis()-inner_loop_time) +"\n"); 

        return;
    } // compute3DApproxPolynomialField
    
	/** 
	 *  compute the inhomogeneity field using a Chebyshev polynomial basis
	 *  with adequate normalization (one field per class)
	 */
    final private void compute3DSeparatePolynomialFields(int deg) {
		int progress = 0;
        int mod = nx*ny*nz/100; // mod is 1 percent of length

        long inner_loop_time;
		if (verbose) inner_loop_time = System.currentTimeMillis();
		
		// compute the degree and size
		int np = 0;
		if (deg>degree) {
			deg = degree;
			np = Np;
		} else 
		if (deg==0) return;   else 
		if (deg==1) np = 4;   else 
		if (deg==2) np = 10;  else 
		if (deg==3) np = 20;  else 
		if (deg>=4) np = 35;

		computeBackgroundMask();
		
		for (int k=0;k<clusters;k++) if (k!=background) {
			// factor = pol*mems*centroids*pol; img = pol*image
			Matrix img = new Matrix(np,1,0.0f);
			Matrix factor = new Matrix(np,np,0.0f);
	
			for (int x=0;x<nx;x+=subsample) for (int y=0;y<ny;y+=subsample) for (int z=0;z<nz;z+=subsample) {
				if (mask[x][y][z]) {
					// valid image point
					compute3DChebyshevCoefficients(x,y,z,deg);		
					
					float vect = polynomialVectorTerm(x,y,z,k,0);
					float matr = polynomialMatrixTerm(x,y,z,k,0);
					
					for (int n=0;n<np;n++) {
						img.set(n,0, img.get(n,0) + vect*pol[n] );
						for (int m=0;m<np;m++) {
									factor.set(n,m, factor.get(n,m) + matr*pol[n]*pol[m] );
						}
					}	
				}
			}
			for (int n=0;n<np;n++) {
				img.set(n,0, img.get(n,0)/(float)(nx*ny*nz) );
				for (int m=0;m<np;m++)
					factor.set(n,m, factor.get(n,m)/(float)(nx*ny*nz) );
			}
			if (debug) {
				//System.out.println("factor: "+factor.matrixToString(5,2)+"\n");
				//System.out.println("img: "+img.matrixToString(5,2)+"\n");
			}
			// coeff = factor^-1 img
			Matrix coeff = factor.solve(img);
			
			//if (debug) System.out.println("coeff: "+coeff.matrixToString(5,2)+"\n");
			
			// transfer to the field
			float mean = 0.0f;
			float den = 0.0f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				float val = 0.0f;
				if (mask[x][y][z]) {
					compute3DChebyshevCoefficients(x,y,z,deg);		
					for (int n=0;n<np;n++) {
						val += coeff.get(n,0)*pol[n];
					}
					mean += val;
					den++;
					fields[x][y][z][k] = val;
				}
			}
			mean = mean/den;
			// normalize and fill the masked area with 1 rather than 0
			// in case of moving masks
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
				fields[x][y][z][k] += 1.0f - mean;
			}
			
			// clean-up
			factor = null;
			img = null;
			coeff = null;
		}		
		
        if (debug) System.out.print("inner loop time: (milliseconds): " + (System.currentTimeMillis()-inner_loop_time) +"\n"); 

        return;
    } // computeNormalizedPolynomialField

	private final float optimizedApproximationWeight(int x, int y, int z, int k) {
		float weight = 1.0f;
		
		if (objType[k].equals("optimized_csf")) {
			float dist = (field[x][y][z]*image(x,y,z)-centroids[k])*(field[x][y][z]*image(x,y,z)-centroids[k]);
			if (dist>ZERO) weight = (float)(powercsf.lookup(dist, pcsf)/dist);
		} else
		if (objType[k].equals("optimized_gm")) {
			float dist = (field[x][y][z]*image(x,y,z)-centroids[k])*(field[x][y][z]*image(x,y,z)-centroids[k]);
			if (dist>ZERO) weight = (float)(powergm.lookup(dist, pgm)/dist);
		} else
		if (objType[k].equals("optimized_wm")) {
			float dist = (field[x][y][z]*image(x,y,z)-centroids[k])*(field[x][y][z]*image(x,y,z)-centroids[k]);
			if (dist>ZERO) weight = (float)(powerwm.lookup(dist, pwm)/dist);
		}
		return weight;
	}
	
	private final float polynomialVectorTerm(int x, int y, int z, int k, int a) {
		switch (factorType) {
			case IMAGE: 				return membership(x,y,z,k,a)*image(x,y,z)*centroids[k]; 
			case CENTROIDS: 			return membership(x,y,z,k,a)*image(x,y,z)*centroids[k]; 
			case OPTIMIZED_IMAGE: 		return membership(x,y,z,k,a)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*centroids[k]; 
			case OPTIMIZED_CENTROIDS: 	return membership(x,y,z,k,a)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*centroids[k]; 
			default:					return 0.0f;
		}
	}
	
	private final float polynomialMatrixTerm(int x, int y, int z, int k, int a) {
		switch (factorType) {
			case IMAGE: 				return membership(x,y,z,k,a)*image(x,y,z)*image(x,y,z); 
			case CENTROIDS: 			return membership(x,y,z,k,a)*centroids[k]*centroids[k]; 
			case OPTIMIZED_IMAGE: 		return membership(x,y,z,k,a)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*image(x,y,z); 
			case OPTIMIZED_CENTROIDS: 	return membership(x,y,z,k,a)*optimizedApproximationWeight(x,y,z,k)*centroids[k]*centroids[k]; 
			default:					return 0.0f;
		}
	}
	
	private final float splineVectorTerm(int x, int y, int z, int k) {
		switch (factorType) {
			case IMAGE: 				return membership(x,y,z,k,0)*image(x,y,z)*(centroids[k]-image(x,y,z)); 
			case CENTROIDS: 			return membership(x,y,z,k,0)*image(x,y,z)*(centroids[k]-image(x,y,z)); 
			case OPTIMIZED_IMAGE: 		return membership(x,y,z,k,0)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*(centroids[k]-image(x,y,z)); 
			case OPTIMIZED_CENTROIDS: 	return membership(x,y,z,k,0)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*(centroids[k]-image(x,y,z)); 
			default:					return 0.0f;
		}
	}
	
	private final float splineMatrixTerm(int x, int y, int z, int k) {
		switch (factorType) {
			case IMAGE: 				return membership(x,y,z,k,0)*image(x,y,z)*image(x,y,z); 
			case CENTROIDS: 			return membership(x,y,z,k,0)*centroids[k]*centroids[k]; 
			case OPTIMIZED_IMAGE: 		return membership(x,y,z,k,0)*optimizedApproximationWeight(x,y,z,k)*image(x,y,z)*image(x,y,z); 
			case OPTIMIZED_CENTROIDS: 	return membership(x,y,z,k,0)*optimizedApproximationWeight(x,y,z,k)*centroids[k]*centroids[k]; 
			default:					return 0.0f;
		}
	}
	
	private final float membership(int x, int y, int z, int k, int a) {
		switch (memType) {
			case EXACT:		return exactMembership(x,y,z,k,a);
			case APPROX:	return approxMembership(x,y,z,k,a);
			default:		return exactMembership(x,y,z,k,a);
		}
	}
		
	private final float exactMembership(int x, int y, int z, int k, int a) {
		float mem;
		switch (q) {
			case 2: 	mem = mems[x][y][z][k]*mems[x][y][z][k]; break;
			case 1: 	mem = mems[x][y][z][k]; break;
			case 3: 	mem = mems[x][y][z][k]*mems[x][y][z][k]*mems[x][y][z][k]; break;
			case Q: 	mem = (float)powerq.lookup(mems[x][y][z][k]); break;
			default: 	mem = 1.0f; break;
		}
		switch (normType) {
			case EQUAL_VARIANCE:		return mem;
			case SEPARATE_VARIANCES:	return mem/var[k];
			default:					return mem;
		}
	}
	
	private final float approxMembership(int x, int y, int z, int k, int a) {
		float mem;
		switch (q) {
			case 2: 	mem = mems[x][y][z][a]*mems[x][y][z][a]; break;
			case 1: 	mem = mems[x][y][z][a]; break;
			case 3: 	mem = mems[x][y][z][a]*mems[x][y][z][a]*mems[x][y][z][a]; break;
			case Q: 	mem = (float)powerq.lookup(mems[x][y][z][a]); break;
			default: 	mem = 1.0f; break;
		}
		switch (normType) {
			case EQUAL_VARIANCE:		return mem;
			case SEPARATE_VARIANCES:	return mem/var[k];
			default:					return mem;
		}
	}
	
	private final float image(int x, int y, int z) {
		float img;
		if (coordinateType==TRANSFORMED) {
			float xT = (transform[0][0]*x*rx + transform[0][1]*y*ry + transform[0][2]*z*rz + transform[0][3])/rx;
			float yT = (transform[1][0]*x*rx + transform[1][1]*y*ry + transform[1][2]*z*rz + transform[1][3])/ry;
			float zT = (transform[2][0]*x*rx + transform[2][1]*y*ry + transform[2][2]*z*rz + transform[2][3])/rz;
			
			return ImageFunctions.linearInterpolation(image,xT,yT,zT,nx,ny,nz);
		} else {
			return image[x][y][z];
		}
	}
	
	/** compute the 3D Chebyshev polynomial basis for one point */
	private final void compute3DChebyshevCoefficients(int x, int y, int z, int deg) {
		float rx,ry,rz;
		
		for (int n=0;n<Np;n++) {
			pol[n] = 0.0f;
		}
		
		rx = x/(float)(nx-1);
		ry = y/(float)(ny-1);
		rz = z/(float)(nz-1);

		if (deg==1) {
			pol[0] = 1.0f;
			// x
			pol[1] = rx;
			pol[2] = ry;
			pol[3] = rz;
		} else if (deg==2) {
			pol[0] = 1.0f;
			// x
			pol[1] = rx;
			pol[2] = ry;
			pol[3] = rz;
			// x^2
			pol[4] = 2.0f*rx*rx - 1.0f;
			pol[5] = 2.0f*ry*ry - 1.0f;
			pol[6] = 2.0f*rz*rz - 1.0f;
			// xy
			pol[7] = pol[1]*pol[2];
			pol[8] = pol[2]*pol[3];
			pol[9] = pol[3]*pol[1];
		} else if (deg==3) {
			pol[0] = 1.0f;
			// x
			pol[1] = rx;
			pol[2] = ry;
			pol[3] = rz;
			// x^2
			pol[4] = 2.0f*rx*rx - 1.0f;
			pol[5] = 2.0f*ry*ry - 1.0f;
			pol[6] = 2.0f*rz*rz - 1.0f;
			// xy
			pol[7] = pol[1]*pol[2];
			pol[8] = pol[2]*pol[3];
			pol[9] = pol[3]*pol[1];
			// x^3
			pol[10] = 4.0f*rx*rx*rx - 3.0f*rx;
			pol[11] = 4.0f*ry*ry*ry - 3.0f*ry;
			pol[12] = 4.0f*rz*rz*rz - 3.0f*rz;
			// xy^2
			pol[13] = pol[4]*pol[2];
			pol[14] = pol[1]*pol[5];
			pol[15] = pol[5]*pol[3];
			pol[16] = pol[2]*pol[6];
			pol[17] = pol[6]*pol[1];
			pol[18] = pol[3]*pol[4];
			// xyz
			pol[19] = pol[1]*pol[2]*pol[3];
		} else if (deg>=4) {
			pol[0] = 1.0f;
			// x
			pol[1] = rx;
			pol[2] = ry;
			pol[3] = rz;
			// x^2
			pol[4] = 2.0f*rx*rx - 1.0f;
			pol[5] = 2.0f*ry*ry - 1.0f;
			pol[6] = 2.0f*rz*rz - 1.0f;
			// xy 
			pol[7] = pol[1]*pol[2];
			pol[8] = pol[2]*pol[3];
			pol[9] = pol[3]*pol[1];
			// x^3
			pol[10] = 4.0f*rx*rx*rx - 3.0f*rx;
			pol[11] = 4.0f*ry*ry*ry - 3.0f*ry;
			pol[12] = 4.0f*rz*rz*rz - 3.0f*rz;
			// x^2y 
			pol[13] = pol[4]*pol[2];
			pol[14] = pol[1]*pol[5];
			pol[15] = pol[5]*pol[3];
			pol[16] = pol[2]*pol[6];
			pol[17] = pol[6]*pol[1];
			pol[18] = pol[3]*pol[4];
			// xyz 
			pol[19] = pol[1]*pol[2]*pol[3];
			// x^4
			pol[20] = 8.0f*rx*rx*rx*rx - 8.0f*rx*rx + 1.0f;
			pol[21] = 8.0f*ry*ry*ry*ry - 8.0f*ry*ry + 1.0f;
			pol[22] = 8.0f*rz*rz*rz*rz - 8.0f*rz*rz + 1.0f;
			// x^3y 
			pol[23] = pol[10]*pol[2];
			pol[24] = pol[10]*pol[3];
			pol[25] = pol[11]*pol[3];
			pol[26] = pol[11]*pol[1];
			pol[27] = pol[12]*pol[1];
			pol[28] = pol[12]*pol[2];
			// x^2y^2 
			pol[29] = pol[4]*pol[5];
			pol[30] = pol[5]*pol[6];
			pol[31] = pol[6]*pol[4];
			// x^2yz 
			pol[32] = pol[4]*pol[2]*pol[3];
			pol[33] = pol[1]*pol[5]*pol[3];
			pol[34] = pol[1]*pol[2]*pol[6];
		} 
		return;
	}

	final private void computeBackgroundMask() {
        if (backgroundType==CLASS_BASED) return;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (segmentation[x][y][z]==background) mask[x][y][z] = false;
			else mask[x][y][z] = true;
		}
		return;
	}

}
