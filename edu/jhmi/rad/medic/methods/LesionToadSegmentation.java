package edu.jhmi.rad.medic.methods;

import edu.jhmi.rad.medic.structures.*;
import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.utilities.*;


/**
 *
 *  This is modified version of TOADS adapated for handling MS lesions in WM
 *	
 *
 *	@version    Dec 2011
 *	@author     Navid Shiee
 *		
 *
 */
 
public class LesionToadSegmentation {
	
	// object types
	public static final byte	OBJECT		=	1;
	public static final byte	MASK		=	2;
	public static final byte	OUTLIER		=	3;
	public static final byte	LESION		=	4;
	public static final byte	OUTMASK		=	5;
	public static final byte    BLACKHOLE   =   6;
	public static final byte	OPTIMIZED_CSF		=	7;
	public static final byte	OPTIMIZED_GM		=	8;
	public static final byte	OPTIMIZED_WM		=	9;
	public static final byte	RAYLEIGH		=10;
	public static final byte	OPTIMIZED_LESION		=	11;
	public static final byte	OPTIMIZED_BLACKHOLE		=	12;
			
	
	private  static final float minLesionVolume = 10;
	private boolean prone_lesions = true;
	// numerical quantities
	private static final	float   INF=1e10f;
	private static final	float   ZERO=1e-10f;
	private	static final 	float	INVSQ2 = 1.0f/(float)Math.sqrt(2.0f); 
	private	static final 	float	INVSQ3 = 1.0f/(float)Math.sqrt(3.0f);
	private	static final 	float	SQR2 = (float)Math.sqrt(2.0f);
	private	static final 	float	SQR3 = (float)Math.sqrt(3.0f);
	private static final	float	SUM26 = 6.0f+12.0f*INVSQ2+8.0f*INVSQ3;
	private static final	float 	LOGPI2MPI2 	= 	(float)Math.log(Math.PI/2.0f) - (float)Math.PI/2.0f;
	
	
	// invalid label flag, must be <0
	private	static	final	byte	EMPTY = -1;
			
	// data and membership buffers
	private 	float[][][][] 	images;  			// original images
	private 	float[][][][]   mems;   			// membership functions for each class
	private 	float[][][][]   mems_lesion_corrected;   			// membership functions for each class
	private 	byte[][][]		classification;  	// hard classification
	private     short[][][]     lesionClassification; // hard classifcation for lesion class based on topology
	private     short[][][]     lesionClassification_mem; // hard classifcation for lesion class based on membership
	private 	float[][][]   	ordering;   		// distance functions computed during the propagation
	private 	float[][] 		centroid;   		// mean intensity for each class
	private 	float[][] 		stddev;   			// stddev of intensity for each class
	private     float[][][]     covariance;         // Full Covarinace Matrix for each class   
	private     float[]         channelVar;			// Variance on each channel
	private     float[]         channelWeight;		// Weight of each channel in computing distances
	private     float[]         lesionWeight;
	private 	BinaryHeap4D  	tree;		   		// the binary tree used for the fast marching
	private static	int 		nx,ny,nz,nc;   		// images dimensions
	private static	int 	  	classes;			// total number of classes (i.e. except outliers)
	private		byte[]			templateLabel;		// intensity labels for the template images
	private		byte[]			objType;			// The type of object for each class
	private		float[]			volume;				// Count the volume of each object
	private		float[]			Ihigh;				// The min max of the centroid range
	private     float[]         Iscale;             // Scale by which images has been normalized
	private		float[]			lesionCentroid;		// the centroids for the lesion class
	private     float[]         blackHoleCentroid;  // the centroids for the T1 Black Hole class
	private		float[][][]		lesionMembership;	// the membership for the lesion class
	private		float[][][]     blackHoleMem ;
	private     float[][][]     lesionRelation;		// Relationship function for lesion
	private 	float[] 		lesionStDev;	   	// stddev of intensity for the lesion class
	private     float[][]       lesionCov;          //Covariance Matrix for lesions  
	private     float[]         blackHoleStdDev;    // stddev of intensity for the T1 Black Hole class
	private     short           lesionClass; 		// id of the class corresponds to lesion
	private		float[][][]		outlierMembership;	// the membership for the outlier class
	private 	float[][] 		priorCentroid;   		// prior intensity for each class
	private 	float[][]		varweight;
	private		float[]			smoothingFactor;
		
	
	private		CriticalPointLUT	lut;				// the LUT for critical points
		
	// parameters
	private 	float   		smoothing;   	// RFCM smoothing parameter
	private 	float			limit;			// The boundary for the skeleton
	private		float			volumeRatio;	// The limit for the volume loss in thinning
	private		float			maxDistance;	// The limit between regular and forced values, in the approx case
	private		float			outlier;		// The distance at which a value is considered an outlier
	
	
	private 	short[][][]		relationMap;   	// labels for the relations between neighboring objects (intensity only, simple)
	private		boolean[][]		relationList;		// the relations ordered by sets
	private 	byte[][][]		segmentation;   	// current classification for relations
	private		byte[][][]      max_mem_segmentation;

	private		int				Nrelations;
	private		byte[][][] 	skeleton;
	
//	 variable metrics
	private		float[]	pcsf;	
	private		float[]	pgm;	
	private		float[]	pwm;	
	private 	float[] plesion;
	private		PowerTable[]	powercsf,powergm,powerwm;
	private 	PowerTable[]	powerlesion;
	
	// new stuff
	boolean[][][][]				available;		// to check wether the point can be added to the tree
	
	public static final int		RC = 0;
	public static final int		RX = 1;
	public static final int		RY = 2;
	public static final int		RZ = 3;
	
	public 				int		evolutionMode;
	public static final int		NONE = 0;
	public static final int		NORMAL_DISTANCE = 1;
	public static final int		CONSTRAINED_DISTANCE = 2;
	
	private		byte[][][][]	boundaryDist;		// the outside distance to a boundary, for better relation weighting

	// atlas prior
	private		float			priorCoefficient = 0.05f;	// the amount of the energy spent on priors 
	private		float			priorScale = 0.1f;			// the scale for the centroid difference in priors

	private		DemonToadDeformableAtlas atlas;
	private		String[]		modality;
	int t1,t2,flair ; 
	private		float[][] 		similarity;
	private     float[]         lesionSimilarity;
	private     float[]			blackHoleSimilarity;
	
    //	 centroid constraints
	private		String			centroidMode;
	private		float			centroidSmoothness;
	private		int[][]			centroidLink;
	
	//Lesion Modelling
	private		byte[]			structureDistances;	
	private     byte[][][]     VentDist;			//Distance from Ventricles hard segmentation
	private     byte[][][]     GMDist;			    //Distance from Gray Matter hard segmentation
	private     byte[][][]     CSFDist;			//Distance from CSF hard segmentation
	private     byte[][][]     BSTEMDist;			//Distance from Brain Stem hard segmentation
	private     float[][][]     ModifedVENTDist;
	private     boolean[][][]   VENT;               //Hard segmentation of Ventircles based on membership
	private     boolean[][][]   GM;                 //Hard segmentation of Gray Matter based on membership
	private     boolean[][][]   CSF;                //Hard segmentation of CSF based on membership
	private     boolean[][][]   BSTEM;               //Hard segmentation ofBSTEM based on membership
	private     boolean  GMflag = false;
	private     boolean  useLesionMax = false;  
	private     byte    MaxVentDist = 9;
	private     byte   MaxGMDist = 5;
	private     byte   MaxCSFDist = 5;
	private     byte   MaxBstemDist = 5;
	private     byte   MaxModifiedVentDist = 10;
	private   	boolean   useLesionWeight = false; //ratio of previous value of lesionCentroids that affect the new value
	private		float	spread = 0.5f;    // Conventional Distance inside structure    
	private		float	distanceFactor=0;
	
	// Norm used for computing distances
	
	private     String normType = "Diagonal Covariance";
	
	// option flags
	private 	boolean 		useField;
	private 	float[][][][]	field;  			// inhomogeneity field
	
	private 	boolean         hasLesion;
	private		boolean			hasMPRAGE = false;
	
	
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=false;
	static final boolean		verbose=true;
	 
	
    
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public LesionToadSegmentation(float[][][][] images_, DemonToadDeformableAtlas atlas_,
						int nx_, int ny_, int nz_, int nc_,
						float smooth_, float out_,
						float firstLimit_, float lastLimit_, 
						int steps_, float spread_, 
						short maxGMDist_, short maxBstemDist_, short maxVentDist_,
						float firstPriorCoeff_, float priorScale_,
						String[] imgModal_,
						String algoMode, String normType_, String connectivityType,
						String centMode, float centSmo_) {
		images = images_;
		atlas = atlas_;
		classification = atlas.getTemplate();
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nc = nc_;
				
		
		// object information
		classes = atlas.getNumber();
		objType = new byte[classes];
		hasLesion = false;
		smoothingFactor = new float[classes];
		for (int k=0;k<classes;k++) {
			smoothingFactor[k] =1.0f;
			if (atlas.getTopology()[k].equals("obj"))   objType[k] = OBJECT;
			else if (atlas.getTopology()[k].equals("mask"))  objType[k] = MASK;
			else if (atlas.getTopology()[k].equals("out"))   objType[k] = OUTLIER;
			else if (atlas.getTopology()[k].equals("lesion")){
				objType[k] = LESION;
				hasLesion = true;
			}else if (atlas.getTopology()[k].equals("outmask")) objType[k] = OUTMASK;
			else if (atlas.getTopology()[k].equals("blackhole")){
				objType[k] = BLACKHOLE;
				hasLesion = true;
			}else if (atlas.getTopology()[k].equals("optimized_csf"))   objType[k] = OPTIMIZED_CSF;
			else if (atlas.getTopology()[k].equals("optimized_gm"))   objType[k] = OPTIMIZED_GM;
			else if (atlas.getTopology()[k].equals("optimized_wm"))   objType[k] = OPTIMIZED_WM;
			else if (atlas.getTopology()[k].equals("rayleigh"))   objType[k] = RAYLEIGH;
			else if (atlas.getTopology()[k].equals("optimized_lesion")){
				objType[k] = OPTIMIZED_LESION;
				hasLesion = true;
			}
			else if (atlas.getTopology()[k].equals("optimized_blackHole")){
				objType[k] = OPTIMIZED_BLACKHOLE;
				hasLesion = true;
			}
			else objType[k] = OBJECT;
		}
		templateLabel = atlas.getLabels();
		
		modality = imgModal_;
		
		t1 = t2 = flair= -1;
		for (int w=0; w<nc; w++)
			if ( (modality[w].equals("T1")) || (modality[w].equals("T1_SPGR")) || (modality[w].equals("T1_MPRAGE")) || (modality[w].equals("MPRAGE")) )
				t1 = w;
			else if ( modality[w].equals("T2"))
				t2 = w;
			else if (modality[w].equals("FLAIR"))
				flair =w;
		if (debug) System.out.println("T1= " + t1 + " T2= "+ t2 + " FLAIR= "+ flair);
		
		centroidMode = centMode;
		centroidSmoothness = centSmo_;
		centroidLink = atlas.getIntensityGroups(modality,nc);
		
		normType = normType_;
		
		// coefficient (assumes normalized imagess)
		if (hasLesion)
			smoothing = smooth_/6.0f/(float)(classes+1);
		else
			smoothing = smooth_/6.0f/(float)(classes);
		
		outlier = out_;
		
		limit = firstLimit_;
		volumeRatio = lastLimit_;
		
				
		priorCoefficient = firstPriorCoeff_/(float)classes;
		priorScale = priorScale_;
		
		spread = spread_;
		MaxGMDist = (byte)maxGMDist_;
		MaxBstemDist = (byte)maxBstemDist_;
		MaxVentDist = (byte)maxVentDist_;
		lesionClass = -1;
		
		
		if (hasLesion) for (short k=0;k<classes;k++){
			if (objType[k]==LESION || objType[k] == BLACKHOLE || objType[k] == OPTIMIZED_LESION){
				lesionClass = k;
			}
		}
		
		if (algoMode.equals("no distance"))
			evolutionMode = NONE;
		else if (algoMode.equals("normal distance"))
			evolutionMode = NORMAL_DISTANCE;
		else
			evolutionMode = CONSTRAINED_DISTANCE;
		//evolutionMode = NORMAL_DISTANCE;
		
		if (debug){
			for (int i=0;i<objType.length;i++) System.out.println(":"+objType[i]+":");
			System.out.print("\n");
		}
		
			
		useField = false;
		field = new float[nc][][][];
		for (int c=0;c<nc;c++) field[c] = null;
		
			
		// init all the arrays
		try {
			mems = new float[nx][ny][nz][classes];
			mems_lesion_corrected = new float[nx][ny][nz][classes];
			lesionClassification = new short [nx][ny][nz];
			lesionClassification_mem = new short [nx][ny][nz];
			centroid = new float[nc][classes];
			stddev = new float[nc][classes];
			covariance = new float[classes][nc][nc];
			Ihigh = new float[nc];
			Iscale = new float[nc];
			lesionCentroid = new float[nc];
			blackHoleCentroid = new float[nc];
			lesionMembership = new float[nx][ny][nz];
			blackHoleMem = new float[nx][ny][nz];
			lesionRelation = new float[nx][ny][nz];
			lesionStDev = new float[nc];
			lesionCov = new float[nc][nc];
			blackHoleStdDev = new float[nc];
			channelVar = new float[nc];
			channelWeight = new float[nc];
			lesionWeight = new float[nc];
			outlierMembership = new float[nx][ny][nz];
			structureDistances = new byte[classes*nx*ny*nz];
			VentDist = new byte[nx][ny][nz];
			GMDist = new byte[nx][ny][nz];
			CSFDist = new byte[nx][ny][nz];
			BSTEMDist = new byte[nx][ny][nz];
			VENT = new boolean[nx][ny][nz];
			GM = new boolean[nx][ny][nz];
			BSTEM = new boolean[nx][ny][nz];
			ModifedVENTDist = new float[nx][ny][nz];
			CSF= new boolean[nx][ny][nz];
			boundaryDist = new byte[nx][ny][nz][classes];
			ordering = new float[nx][ny][nz];
			tree = new BinaryHeap4D(nx*ny+ny*nz+nz*nx, BinaryHeap4D.MAXTREE);
			if (connectivityType.equals("26/6")) lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
			else if (connectivityType.equals("6/26")) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
			else if (connectivityType.equals("18/6")) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
			else if (connectivityType.equals("6/18")) lut = new CriticalPointLUT("critical618LUT.raw.gz",200);
			else if (connectivityType.equals("6/6")) lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
			else lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
			lut.loadCompressedPattern();
			relationMap = new short[nx][ny][nz];
			segmentation = new byte[nx][ny][nz]; 
			max_mem_segmentation = new byte[nx][ny][nz];
			skeleton = new byte[nx][ny][nz];
			available = new boolean[nx][ny][nz][classes];
			volume = new float[classes];
			similarity = new float[classes][classes];
			lesionSimilarity = new float[classes];
			blackHoleSimilarity = new float[classes];
			pcsf = new float[nc];
			pgm = new float[nc];
			pwm = new float[nc];
			plesion = new float[nc];
			priorCentroid = new float[nc][classes];

			
			varweight = atlas.getIntensityVariancePriors(modality,nc);
			//Optimzed q stuff
			boolean optimize = false;
			for (int k=0;k<classes;k++) if (atlas.getTopology()[k].startsWith("optimized_")) {
				optimize =true;
				powercsf = new PowerTable[nc];
				powergm = new PowerTable[nc];
				powerwm = new PowerTable[nc];
				break;
			}
			if (optimize)for (int c=0; c<nc; c++) {
				if (debug) System.out.println(modality[c] + " : " + pcsf[c] + " "+pgm[c]+" "+pwm[c] );
				powercsf[c] = new PowerTable(0.0f , 4.0f , 0.000001f , pcsf[c] );
				powergm[c] = new PowerTable(0.0f , 4.0f , 0.000001f , pgm[c] );
				powerwm[c] = new PowerTable(0.0f , 4.0f , 0.000001f , pwm[c] );
			}
			
			//Special case for mprage only
			for (int c=0; c<nc; c++) if (modality[c].equals("T1_MPRAGE")){
				hasMPRAGE = true;
			}
			if (hasMPRAGE){
				powercsf = new PowerTable[nc];
				for (int n=0;n<nc;n++){ 
					pcsf[n] = atlas.getOptimizedFactor(modality[n])[0];
				}
				for (int c=0; c<nc; c++)
					powercsf[c] = new PowerTable(0.0f , 4.0f , 0.000001f , pcsf[c] );
			}
			powerlesion = new PowerTable[nc];
			for (int n=0;n<nc;n++){ 
				plesion[n] = atlas.getOptimizedFactor(modality[n])[3];
				if (debug) System.out.print(modality[n] + " " + plesion[n]+ " | ");
			}
			if (debug) System.out.println();
			for (int c=0; c<nc; c++)
				powerlesion[c] = new PowerTable(0.0f , 4.0f , 0.000001f , plesion[c] );
			
      } catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;
		
		// init values
		for (int k=0;k<classes;k++) volume[k] = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			ordering[x][y][z] = 0.0f;
			lesionClassification[x][y][z] =0;
			lesionClassification_mem[x][y][z] =0;
			segmentation[x][y][z] = EMPTY;
			skeleton[x][y][z] = EMPTY;
			relationMap[x][y][z] = EMPTY;	
			VentDist[x][y][z] = (byte)(MaxVentDist + 1);
			GMDist[x][y][z] = (byte)(MaxGMDist + 1) ;
			BSTEMDist[x][y][z] = (byte)(MaxBstemDist + 1);
			ModifedVENTDist[x][y][z] = 0;
			CSFDist[x][y][z] = 0;
			GM[x][y][z] = false;
			BSTEM[x][y][z] = false;
			CSF[x][y][z]=false;
			for (byte k=0;k<classes;k++) {
				mems[x][y][z][k] = 0.0f;
				mems_lesion_corrected[x][y][z][k] = 0.0f;
				if (classification[x][y][z] == templateLabel[k]) {
					volume[k]++;	
					segmentation[x][y][z] = k;
					skeleton[x][y][z] = k;
					relationMap[x][y][z] = k;
					for (int c=0;c<nc;c++) {
						centroid[c][k] += images[c][x][y][z];
					}
				}
				available[x][y][z][k] = false;
				boundaryDist[x][y][z][k] = 0;
			}
			if (segmentation[x][y][z]==EMPTY) {
				System.out.println("error: empty regions in original template - aborting");
				isWorking = false;
				return;
			}
			lesionMembership[x][y][z] = 0.0f;
			lesionRelation[x][y][z] = 1.0f;
			outlierMembership[x][y][z] = 0.5f;
			blackHoleMem[x][y][z] = 0.0f;
		}
		for (int c=0;c<nc;c++) {
			for (int k=0;k<classes;k++) {
				centroid[c][k] = centroid[c][k]/volume[k];
				stddev[c][k] = 1.0f;
			}
			Ihigh[c] = 1.0f;
			Iscale[c] = 1.0f;
			lesionCentroid[c] = 0.0f;
			lesionStDev[c] = 1.0f;
			blackHoleCentroid[c] = 0.0f;
			blackHoleStdDev[c] = 1.0f;
			channelVar[c] = 1.0f;
			channelWeight[c] = (float)(1.0f/nc);
			lesionWeight[c] = (float)(1.0f/nc);
		}
		for (int k=0;k<classes;k++) {
			lesionSimilarity[k] = 1.0f;
			blackHoleSimilarity[k] =1.0f;
			for (int m=0;m<classes;m++) {
				if (k==m) similarity[k][m] = 0.0f;
				else similarity[k][m] = 1.0f;
			}
			for (int c1=0; c1<nc; c1++)
				for (int c2=0; c2<nc; c2++) if (c1==c2)
					covariance[k][c1][c2] = 1.0f;
				else
					covariance[k][c1][c2] = 0.0f;
					
		}
		if (hasLesion){ 
			lesionSimilarity[lesionClass] = 0.0f;
			blackHoleSimilarity[lesionClass] = 0.0f;
		}
		
		for (int c1=0; c1<nc; c1++)
			for (int c2=0; c2<nc; c2++) if (c1==c2)
				lesionCov[c1][c2] = 1.0f;
			else
				lesionCov[c1][c2] = 0.0f;
		if (debug) System.out.println("initialization\n");
	}
		
	public void finalize() {
		images = null;
		mems = null;
		mems_lesion_corrected=null;
		centroid = null;
		tree = null;
		lut = null;
		boundaryDist = null;
		VentDist = null;
		ModifedVENTDist = null;
		GMDist =null;
		System.gc();
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
        images = null;
		tree.finalize();
		tree = null;
		lut.finalize();
		lut = null;
		System.gc();
	}
    
	public final void setGMflag() { GMflag = true;}
	public final void reSetGMflag() { GMflag = false;}
	public final void setDistanceFactor (float factor_) { distanceFactor = factor_;}
	public final float[] getCentroids(int c) { return centroid[c]; }
	public final float[][] getCentroids() { return centroid; }
	//public final void setCentroids(float[][] cent) { centroid = cent; }
	public final void setLesionWeight() { useLesionWeight = true; }
	public final void reSetLesionWeight() { useLesionWeight = false; }
	final public void setOptimizedFactors(String[] modality, int nc) {
		for (int n=0;n<nc;n++) {
			pcsf[n] = atlas.getOptimizedFactor(modality[n])[0];
			pgm[n]  = atlas.getOptimizedFactor(modality[n])[1];
			pwm[n]  = atlas.getOptimizedFactor(modality[n])[2];
			plesion[n] = atlas.getOptimizedFactor(modality[n])[3];
			powercsf[n] = new PowerTable(0.0f , 4.0f , 0.000001f , pcsf[n] );
			powergm[n] = new PowerTable(0.0f , 4.0f , 0.000001f , pgm[n] );
			powerwm[n] = new PowerTable(0.0f , 4.0f , 0.000001f , pwm[n] );
			powerlesion[n] = new PowerTable(0.0f , 4.0f , 0.000001f , plesion[n] );
		}
		return;
	}
	public final float[][][][] getMemberships() { return mems; }
	public final float[][][][] getCorrectdMemberships() { return mems_lesion_corrected;}
	public final float[][][] getMembership(int k) {
		float[][][] mem = new float[nx][ny][nz];
		for (int x =0; x<nx; x++)for (int y=0; y<ny; y++)for (int z=0; z<nz; z++)
			mem[x][y][z] = mems[x][y][z][k];
		return mem; 
	}
	public final void setMemberships(float[][][][] m_) { mems = m_; }
	
	public final float[][][] getLesionMemberships() { return lesionMembership; }
	public final void setLesionMemberships(float[][][] m_) { lesionMembership = m_; }

	public final float[][][] getBlackHoleMem() { return blackHoleMem; }
	public final void setBlackHoleMem(float[][][] m_) { blackHoleMem = m_; }
	
	public final void setOrdering(float[][][] o_) { ordering = o_; }
	public final byte[][][] getClassification() { return classification; }
	public final byte[][][] getSegmentation() {return segmentation;}
	public final short[][][] getLesionClassification() { return lesionClassification; }
	public final float[][][] getLesionRelation() { return lesionRelation; }

	public final void setIntensityMax(float[] I) { Ihigh = I; }
	public final void setIntensityScale(float[] I) { Iscale = I; }
	
	public final void setGMDistance (byte d) { MaxGMDist = d; }
	
	public final float[] getLesionCentroid() { return lesionCentroid ; }
	public final void setLesionCentroid(float[] cent) { lesionCentroid = cent; }
	
	public final float[] getBlackHoleCentroid() { return blackHoleCentroid ; }
	public final void setBlackHoleCentroid(float[] cent) { blackHoleCentroid = cent; }
	
	public final void setPriorCoefficients(float pC, float pS) { 
		priorCoefficient = pC;
		priorScale = pS;
	}
	
	public void setModifiedVentDistance (byte dist_) { MaxVentDist = dist_;}
	public final void resetClassification() { 
		classification = atlas.getTemplate();
	}

	public final void initCentroids(float[][] cent) { 
		for (int  c=0;c<nc;c++) {
			for (int k=0;k<classes;k++) {
				centroid[c][k] = cent[c][k];
				priorCentroid[c][k] = cent[c][k];
			}
		}
	}
	
	/** add inhomogeneity correction */
    public final void addInhomogeneityCorrection(float[][][] field_, int c) {
        field[c] = field_;
        useField = true;
    }
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
	
	private final boolean isNotEmpty() {
		return tree.isNotEmpty();
	}
	
	/**
	 *  initialize the labels for a new iteration
	 */
	public final void initGrowingLabels() {
		float   order;
		float	mem;
		
		// fix for homeomorphic version
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (skeleton[x][y][z]==EMPTY) {
				for (int n=0;n<classes;n++) {
					available[x][y][z][n] = true;
				}
			} else {
				for (int n=0;n<classes;n++) {
					available[x][y][z][n] = false;
				}
			}
		}

		// init the labels from the skeleton
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (byte k=0;k<classes;k++) {
				if (available[x][y][z][k]) {
					if ( (skeleton[x-1][y][z]==k) || (skeleton[x+1][y][z]==k)
					  || (skeleton[x][y-1][z]==k) || (skeleton[x][y+1][z]==k)
					  || (skeleton[x][y][z-1]==k) || (skeleton[x][y][z+1]==k) ) {
					  
					  	available[x][y][z][k]=false;
						order = ordering[x][y][z] + maxMembershipRatio(x,y,z,k,bestMembership(x,y,z));
						tree.addValue(order,x,y,z,k);
					}
				}
			}
		}
	}//initGlobalGrowingLabels
		
	/**
	 *  initialize the labels for a new iteration
	 */
	public final void initCompetitionLabels() {
		float   order;
		
		
		// at first, everything can change
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (int n=0;n<classes;n++) {
				available[x][y][z][n] = true;
			}
			skeleton[x][y][z] = segmentation[x][y][z];
		}

		// init the labels from the segmentation
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			for (byte k=0;k<classes;k++) {
				if (segmentation[x][y][z]==k) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						// are there candidates ?
						if (segmentation[x+i][y+j][z+l]!=k) {
							int m = segmentation[x+i][y+j][z+l];
							
							//#### LESIONTOADS #######
							if (m==lesionClass)
								order = (mems[x+i][y+j][z+l][m]+lesionMembership[x+i][y+j][z+l]+blackHoleMem[x+i][y+j][z+l])/(mems[x+i][y+j][z+l][k]+ZERO);
							else if (k==lesionClass)
								order = mems[x+i][y+j][z+l][m]/(mems[x+i][y+j][z+l][k]+lesionMembership[x+i][y+j][z+l]+blackHoleMem[x+i][y+j][z+l]+ZERO);
							else
								order = mems[x+i][y+j][z+l][m]/(mems[x+i][y+j][z+l][k]+ZERO);
							
							if (order > 0) {
								// m is a better label for X; schedule a change
								tree.addValue(order,x,y,z,k+classes*m);
							}
						}
					}
				}
			}
		}
	}//initCompetitionLabels
		
	/**
	 *  initialize the labels for a new iteration
	 */
	public final void initThinningLabels() {
		boolean isBoundary = false;
		float   order=0.0f;
		
		//reset the volume here before first (k==0)
		for (int k=0;k<classes;k++) 
				volume[k] = 0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			skeleton[x][y][z] = segmentation[x][y][z];
				
			for (int k=0;k<classes;k++) {
				if (segmentation[x][y][z]==k) {
					//#####LESIONTOADS######
					
					if (k==lesionClass)
						if (useLesionMax)
							volume[k] += Numerics.max(mems[x][y][z][k], lesionMembership[x][y][z], blackHoleMem[x][y][z]);
						else
							volume[k]+=mems[x][y][z][k]+lesionMembership[x][y][z] + blackHoleMem[x][y][z];
					else
						volume[k] += mems[x][y][z][k];
				}
				// available only inside k
				if (segmentation[x][y][z]==k) {
					available[x][y][z][k] = true;
				} else {
					available[x][y][z][k] = false;
				}
			}
		}
						
		
	
	
		int[] Nbound = new int[classes];
		for (int k=0;k<classes;k++) Nbound[k] = 0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			for (int k=0;k<classes;k++) {
				if (segmentation[x][y][z]==k) {
					// neighbors outside the object ?
					isBoundary = false;
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						if (segmentation[x+i][y+j][z+l]!=k) {
							isBoundary = true;
						}
					}
					
					if (isBoundary) {
						// no search for best change yet
						available[x][y][z][k] = false;
						
						// best to use the membership of the object to be removed
						order = - 1.0f/maxMembershipRatio(x,y,z,k,bestMembership(x,y,z));
						tree.addValue(order,x,y,z,k);
						Nbound[k]++;
					}
				}
			}
		}
		
	}//initGlobalThinningLabels
		
/**
 * Methods for comouting distances from boundary of different structures
 * These distances used to form relationship function for lesions
 */
	
	private void init1Ddistance(){
		long t = System.currentTimeMillis();
		byte[] labels = atlas.getLabels();
		int id;
		for (int k=0; k<classes; k++) for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			id = k*nx*ny*nz+x*ny*nz+y*nz+z;
			if (classification[x][y][z]==labels[k])
				structureDistances[id] = 0; //distance inside object set to zero
			else if (structureDistances[id] !=1)
				structureDistances[id] = -1; /// initial distance -1
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
				if ( classification[x+i][y+j][z+l]==labels[k]){
					if ( structureDistances[id]== -1) //check if it is inside
						structureDistances[id] = 1;
				}
			}
		}
		if (debug) System.out.println("1D initilaiztion time: " + (System.currentTimeMillis() - t));
		labels = null;
			
	}
	public final void compute1Ddistances() {
		
		
		
		
	// diffuse the coefficients till get to the desired point
	int id1,id2;	
	long t = System.currentTimeMillis();
	for (int k=0;k<classes;k++) {
		if (debug) System.out.println("Computing distances from "+atlas.getNames()[k]+" ...");
		boolean stop =false;
		byte count = 1;
		while (!stop) {
			stop = true;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				id1 = k*nx*ny*nz+x*ny*nz+y*nz+z;
				if (structureDistances[id1]==count) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						id2 = k*nx*ny*nz+(x+i)*ny*nz+(y+j)*nz+z+l;
						if ( ( structureDistances[id2]==-1 )) {
							structureDistances[id2] = (byte)(count+1);//*(count+1.0f);
							stop = false;
						}
					}
				}
			}
			count = (byte)(count+1);
			if (count >= 10)
				stop = true;
				
		}
		
	}
	if (debug) System.out.println("1D distance computation: "+ (System.currentTimeMillis()-t));
	return;
			
}
	
	
	final public float[][][] fastDistance (boolean[][][] obj, float maxDist, int nx, int ny, int nz){
		float[][][] dist = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (obj[x][y][z])
				dist[x][y][z] = 0.0f; //distance inside object set to zero
			else
				dist[x][y][z] = -1.0f; // initial distance -1
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz)
				 && obj[x+i][y+j][z+l]){
					if ( dist[x][y][z]== -1.0f) //check if it is inside
						dist[x][y][z] = 1.0f;
				}
			}
		}
		
		boolean stop = false;
		float count = 1.0f;
		while (!stop) {
			stop = true;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (dist[x][y][z]==count) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz) && ( dist[x+i][y+j][z+l]==-1.0f )) {
							dist[x+i][y+j][z+l] = (count+1.0f);//*(count+1.0f);
							stop = false;
						}
					}
				}
			}
			count = count+1.0f;
			if (count >= maxDist)
				stop = true;
				
		}
		return dist;
	}	
	
	
	final public float[][][] fastFloatDistance (boolean[][][] obj, byte maxDist,int nx, int ny, int nz){
		float[][][] dist = new float[nx][ny][nz];
		if (maxDist !=0){
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (obj[x][y][z])
				dist[x][y][z] = 0; //distance inside object set to zero
			else
				dist[x][y][z] = -1; // initial distance -1
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz)
				 && obj[x+i][y+j][z+l]){
					if ( dist[x][y][z]== -1) //check if it is inside
						dist[x][y][z] = 1;
				}
			}
		}
		
		boolean stop = false;
		byte count = 1;
		while (!stop) {
			stop = true;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (dist[x][y][z]==count) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz) && ( dist[x+i][y+j][z+l]==-1.0f )) {
							dist[x+i][y+j][z+l] = count+1.0f;
							stop = false;
						}
					}
				}
			}
			count = (byte)(count+1);
			if (count >= maxDist)
				stop = true;
				
		}
		}else
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) dist[x][y][z]=-1;
		
		return dist;
	}
	final public float[][][] fastNormalizedDistance (boolean[][][] obj, byte maxDist,float minDist, int nx, int ny, int nz){
		float[][][] dist = fastFloatDistance(obj, maxDist, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)  
			
			if (dist[x][y][z] == 0)
				dist[x][y][z] =minDist/(float)maxDist;
			else if (dist[x][y][z] !=-1)
				dist[x][y][z] = dist[x][y][z]/(float)maxDist;
		return dist;
		}	
	
	final public byte[][][] fastDistance (boolean[][][] obj, byte maxDist, int nx, int ny, int nz){
		byte[][][] dist = new byte[nx][ny][nz];

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (obj[x][y][z])
				dist[x][y][z] = 0; //distance inside object set to zero
			else
				dist[x][y][z] = -1; // initial distance -1
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz)
				 && obj[x+i][y+j][z+l]){
					if ( dist[x][y][z]== -1) //check if it is inside
						dist[x][y][z] = 1;
				}
			}
		}
		
		boolean stop = false;
		byte count = 1;
		while (!stop) {
			stop = true;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (dist[x][y][z]==count) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+l>=0) && (z+l<nz) && ( dist[x+i][y+j][z+l]==-1.0f )) {
							dist[x+i][y+j][z+l] = (byte)(count+1);//*(count+1.0f);
							stop = false;
						}
					}
				}
			}
			count = (byte)(count+1);
			if (count >= maxDist)
				stop = true;
				
		}
		
		return dist;
	}	
	
	final public void computeModifiedVentDistance(){
		
		boolean[][][] vent_mod_mask = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			
			if ( atlas.getNames()[max_mem_segmentation[x][y][z]].equals("Ventricles") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("LeftVentricle") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("RightVentricle") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("ThirdVentricle")) 
				vent_mod_mask[x][y][z] = true;
			else
				vent_mod_mask[x][y][z] = false;
		}
		
		vent_mod_mask = Morphology.erodeObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.erodeObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.dilateObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.dilateObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		ModifedVENTDist = fastNormalizedDistance(vent_mod_mask,MaxModifiedVentDist,1.0f,nx,ny,nz);
		if (debug) System.out.println("Modfied Ventricle Distance Computed. \n");
		
	}
		
	final public void computeVentirclsDistance(){
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			
			if ( atlas.getNames()[segmentation[x][y][z]].equals("Ventricles") 
				|| atlas.getNames()[segmentation[x][y][z]].equals("LeftVentricle") 
				|| atlas.getNames()[segmentation[x][y][z]].equals("RightVentricle") 
				|| atlas.getNames()[segmentation[x][y][z]].equals("ThirdVentricle")) 
				VENT[x][y][z] =true;
			else
				VENT[x][y][z] =false;
		}
		
		VentDist = fastDistance(VENT,MaxVentDist,nx,ny,nz);
		
		if (debug) System.out.println("Ventricle distance computed. \n");
		
	}
	
	
	
	final public void computeGMDistance(){
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){

			if ((atlas.getNames()[segmentation[x][y][z]].equals("Cerebrum-GM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("CerebralGM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("CorticalGM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Caudate"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("GM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Thalamus"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Putamen"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Hippocampus"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Amygdala"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Accumbens"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Sulcal-CSF"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("SulcalCSF"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("Subcortical-GM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("SubcorticalGM"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("BasalGanglia"))
				|| (atlas.getNames()[segmentation[x][y][z]].equals("GlobusPallidus")))
			
				GM[x][y][z] =true;
			else
				GM[x][y][z] =false;
		}
		
		GMDist = fastDistance(GM,MaxGMDist,nx,ny,nz);
		if (debug) System.out.println("GM distance computed. \n");
		
	}
	
	final public void computeBstemDistance(){
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			
			if (atlas.getNames()[segmentation[x][y][z]].equals("Brainstem")
			    || (atlas.getNames()[segmentation[x][y][z]].equals("InterVentricular-WM"))
			    || (atlas.getNames()[segmentation[x][y][z]].equals("FornixWM"))	) 
				BSTEM[x][y][z] =true;
			else
				BSTEM[x][y][z] =false;
		}
		
		BSTEMDist = fastDistance(BSTEM,MaxBstemDist,nx,ny,nz);
		if (debug) System.out.println("Outer Ventricular WM distance computed. \n");

		
	}
	final public void computeCSFDistance(){
		
		
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){

				
			if ( atlas.getNames()[segmentation[x][y][z]].equals("Ventricles") 
					|| atlas.getNames()[segmentation[x][y][z]].equals("LeftVentricle") 
					|| atlas.getNames()[segmentation[x][y][z]].equals("RightVentricle") 
					|| atlas.getNames()[segmentation[x][y][z]].equals("ThirdVentricle")) 
				CSF[x][y][z] =true;
			else
				CSF[x][y][z] =false;
		}
		
		CSFDist = fastDistance(CSF,MaxCSFDist,nx,ny,nz);
		
		
		if (debug) System.out.println("CSF Distance Computed. \n");
		
	}
	
final public void computeNormalizedModifiedVentDistance(){
		
		
		boolean[][][] vent_mod_mask = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			
			if ( atlas.getNames()[max_mem_segmentation[x][y][z]].equals("Ventricles") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("LeftVentricle") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("RightVentricle") 
					|| atlas.getNames()[max_mem_segmentation[x][y][z]].equals("ThirdVentricle")) 
				vent_mod_mask[x][y][z] = true;
			else
				vent_mod_mask[x][y][z] = false;
		}
		
		vent_mod_mask = Morphology.erodeObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.erodeObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.dilateObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		vent_mod_mask = Morphology.dilateObject6C(vent_mod_mask, nx,ny,nz, 1,1,1);
		ModifedVENTDist = fastNormalizedDistance(vent_mod_mask,MaxModifiedVentDist,spread,nx,ny,nz);
		
		
		if (debug) System.out.println("Normalized Modified Ventricle Distance Computed. \n");
	}
	
	/** 
	 *  compute the FCM membership functions given the centroids
	 *	with the different options (outliers, field, edges, MRF)
	 *	with q=2
	 */
    final public float computeMemberships() {
    	
    	float dist;
        float den,num;
        float numL,numO, numB;
        float ngb, ngbL, ngbS,ngbB;
		float[] shape = new float[classes];
		numL=0.0f;
		byte relationType;
		byte lesionRelationType;
		byte PRIOR = 3;
		byte DIRECT = 2;
		byte FARAWAY = 0;
		
		float[] energy = new float[classes];
        float lesion_energy =0.0f, blackHole_energy = 0.0f;
				
        // compute the centroid similarity coefficient and scaling factor
		float centroidDist;
		boolean csfFlag = true;
		if (t1==-1) csfFlag = false;
		for (int k=0;k<classes;k++) {
			for (int m=0;m<classes;m++) 
				if (m!=k) {
					centroidDist = 0;
					if ( csfFlag &&( (atlas.getNames()[k].equals("Sulcal-CSF")) 
							|| (atlas.getNames()[k].equals("SulcalCSF")) 
							|| (atlas.getNames()[k].equals("Ventricles"))
							|| (atlas.getNames()[m].equals("Sulcal-CSF")) 
							|| (atlas.getNames()[m].equals("SulcalCSF")) 
							|| (atlas.getNames()[m].equals("Ventricles"))
							)){
						centroidDist +=(centroid[t1][k]-centroid[t1][m])*(centroid[t1][k]-centroid[t1][m])
						 /(stddev[t1][k]+stddev[t1][m])/(stddev[t1][k]+stddev[t1][m]);
						
					}else{
					for (int c=0;c<nc;c++) centroidDist += channelWeight[c]* (centroid[c][k]-centroid[c][m])*(centroid[c][k]-centroid[c][m])
						 /(stddev[c][k]+stddev[c][m])/(stddev[c][k]+stddev[c][m]);
					}
					similarity[k][m] = 1.0f/(1.0f + (1-priorScale)/priorScale*centroidDist);
				}else {
					similarity[k][m] = 0.0f;
				}
		}

		//simlarity for lesion class
		for (int k =0; k<classes; k++){
			centroidDist = 0;
			for (int c=0;c<nc;c++) centroidDist += channelWeight[c]*(lesionCentroid[c]-centroid[c][k])*(lesionCentroid[c]-centroid[c][k])
			/(stddev[c][k]+lesionStDev[c])/(stddev[c][k]+lesionStDev[c]);	
		if ((objType[k] != LESION) && (objType[k] != BLACKHOLE) && (objType[k] != OPTIMIZED_LESION))
			lesionSimilarity[k] = 1.0f/(1.0f + (1-priorScale)/priorScale*centroidDist);
		else //set the similarity with WM to zero
			lesionSimilarity[k] = 0.0f;
		}
		
		
		// main loop
		dist = 0.0f;
		atlas.precomputeTransformMatrix(1.0f);
		
		long time = System.currentTimeMillis();
		if (hasLesion){
			if (MaxVentDist != 0.0f)
				computeVentirclsDistance();
			if (MaxBstemDist != 0.0f)
				computeBstemDistance();
			if (GMflag && (MaxGMDist != 0.0f))
				computeGMDistance();
		}
		
		if ( distanceFactor != 0){
			updateMaxmemClassification();
			//computeNormalizedVentirclsDistance();
			MaxModifiedVentDist = 2;
			computeModifiedVentDistance();
		}
		
		if (debug) System.out.println("Distance time: " + (System.currentTimeMillis() - time ));
		boolean useEverything =true;		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			den = 0;
			
			
			// precompute the priors
			if (priorCoefficient > 0) {
				shape = atlas.getTransformedShape(x,y,z);
			}
			
			for (int k=0;k<classes;k++) {
			
				relationType = PRIOR;
				lesionRelationType = FARAWAY;
				if (relationList[relationMap[x][y][z]][k]) 
					relationType = DIRECT;
				
				//update relationtype for lesions
				
				if ( (objType[k] == LESION) || (objType[k] == BLACKHOLE) || objType[k] == OPTIMIZED_LESION){
					lesionRelationType = relationType;
				}
				
				if (relationType != DIRECT && shape[k] == 0.0f) {
					mems[x][y][z][k] = 0.0f;
					if ((objType[k] == LESION) || objType[k] == BLACKHOLE || objType[k] == OPTIMIZED_LESION){
						lesionMembership[x][y][z] = 0.0f;
						blackHoleMem[x][y][z] = 0.0f;
						lesionRelation[x][y][z] = 0.0f;
					}
					
				} else {
									
					// data term
					num = 0.0f;
					numL = 0.0f;
					numB = 0.0f;
					numO = 0.0f;
										 
					if (objType[k]==OBJECT) {
						
							if (useField) {
								for (int c=0;c<nc;c++) 
									num += varweight[k][c] * channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
							}else {
								for (int c=0;c<nc;c++) 
									num += varweight[k][c] * channelWeight[c] * (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
							}
						
					} else if ( objType[k]==LESION ) {
						if (useField) {
							for (int c=0;c<nc;c++) {
								num += channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
								//num += varweight[c][k] * channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
							}
								if ((flair != -1) && useLesionWeight){
									numL += powerlesion[flair].lookup((field[flair][x][y][z]*images[flair][x][y][z]-lesionCentroid[flair])*(field[flair][x][y][z]*images[flair][x][y][z]-lesionCentroid[flair]),plesion[flair]);
								}else if ( (flair!=-1) && (t2!=-1) && useEverything){
									for (int c=0;c<nc;c++) if (c!=t1) 
										  numL += 0.5f * powerlesion[c].lookup((field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]),plesion[c]);
								}else	
									for (int c=0;c<nc;c++) 
									  numL += channelWeight[c] * powerlesion[c].lookup((field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]),plesion[c]);
							
						} else {
							for (int c=0;c<nc;c++) {
								num += channelWeight[c] * (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
								//num += varweight[c][k] *channelWeight[c] * (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
							}
								if ((flair != -1) && useLesionWeight){
									numL += powerlesion[flair].lookup((images[flair][x][y][z]-lesionCentroid[flair])*(images[flair][x][y][z]-lesionCentroid[flair]),plesion[flair]);
								}else
									for (int c=0;c<nc;c++)
										numL += channelWeight[c] * powerlesion[c].lookup((images[c][x][y][z]-lesionCentroid[c])*(images[c][x][y][z]-lesionCentroid[c]),plesion[c]);
						}
					}else if (objType[k]==BLACKHOLE) {
						if (useField) {
							for (int c=0;c<nc;c++) {
								num += channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
								numL +=channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]);
								numB +=channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-blackHoleCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-blackHoleCentroid[c]);
							}
						} else {
							for (int c=0;c<nc;c++) {
								num += channelWeight[c] * (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
								numL +=channelWeight[c] * (images[c][x][y][z]-lesionCentroid[c])*(images[c][x][y][z]-lesionCentroid[c]);
								numB +=channelWeight[c] * (images[c][x][y][z]-blackHoleCentroid[c])*(images[c][x][y][z]-blackHoleCentroid[c]);
							}
						}
					}/*else if (objType[k]==RAYLEIGH) {
						if (useField) {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*rayleighDistance(field[c][x][y][z]*images[c][x][y][z],centroid[c][k]);
						} else {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*rayleighDistance(images[c][x][y][z],centroid[c][k]);
						}
					}*/ else if (objType[k]==OPTIMIZED_CSF) {
						if (useField){
							if (hasMPRAGE)
								num += powercsf[t1].lookup( (field[t1][x][y][z]*images[t1][x][y][z]-centroid[t1][k])*(field[t1][x][y][z]*images[t1][x][y][z]-centroid[t1][k]), pcsf[t1]);
							else
								num += (field[t1][x][y][z]*images[t1][x][y][z]-centroid[t1][k])*(field[t1][x][y][z]*images[t1][x][y][z]-centroid[t1][k]);
						}else{
							if (hasMPRAGE)
								num += powercsf[t1].lookup( (images[t1][x][y][z]-centroid[t1][k])*(images[t1][x][y][z]-centroid[t1][k]), pcsf[t1]);
							else
								num += (images[t1][x][y][z]-centroid[t1][k])*(images[t1][x][y][z]-centroid[t1][k]);
						}
					}/* else if (objType[k]==OPTIMIZED_GM) {
						if (useField) {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*powergm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]), pgm[c]);
						} else {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*powergm[c].lookup( (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]), pgm[c]);
						}
					} else if (objType[k]==OPTIMIZED_WM) {
						if (useField) {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*powerwm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]), pwm[c]);
						} else {
							for (int c=0;c<nc;c++)
								num += channelWeight[c]*powerwm[c].lookup( (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]), pwm[c]);
						}
					}else if (objType[k]==OPTIMIZED_LESION) {
						if (useField) {
							for (int c=0;c<nc;c++){
								num += channelWeight[c]*powerwm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]), pwm[c]);
								numL += channelWeight[c]*powerwm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]), pwm[c]);
								sumOpt += Math.abs(channelWeight[c]*powerwm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]), pwm[c]) - channelWeight[c] * (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]));
								Opt += channelWeight[c]*powerwm[c].lookup( (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]), pwm[c]);
							}
						} else {
							for (int c=0;c<nc;c++){
								num += channelWeight[c]*powerwm[c].lookup( (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]), pwm[c]);
								numL += channelWeight[c]*powerwm[c].lookup( (images[c][x][y][z]-lesionCentroid[c])*(images[c][x][y][z]-lesionCentroid[c]), pwm[c]);
							}
						}
					}*/else if (objType[k]==MASK) {
						for (int c=0;c<nc;c++){
							num += channelWeight[c] * maskDistance(images[c][x][y][z],0.0f);
							//num += varweight[c][k] *channelWeight[c] * maskDistance(images[c][x][y][z],0.0f);
						}
					} else if (objType[k]==OUTMASK) {
						for (int c=0;c<nc;c++) {
							num  += channelWeight[c]*maskDistance(images[c][x][y][z],0.0f);
							numO += channelWeight[c]*outlier;
							//num  += varweight[c][k] *channelWeight[c]*maskDistance(images[c][x][y][z],0.0f);
							//numO += varweight[c][k] *channelWeight[c]*outlier;
						}
					}else {
						num = outlier;
					}
					if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) {
						energy[k] = num ;
						lesion_energy = numL;						
					}else if (objType[k]==BLACKHOLE ){
						energy[k] = num;
						lesion_energy = numL;
						blackHole_energy = numB;
					}else if (objType[k]==OUTMASK) {
						energy[k] = (1.0f-outlierMembership[x][y][z])*num + outlierMembership[x][y][z]*numO;
					}else {
						energy[k] = num;
					}
					
					// spatial smoothing
					if (smoothing > 0.0f) { 
						ngb = 0.0f;
						ngbS = 0.0f;
						ngbL = 0.0f;
						ngbB = 0.0f;
						// case by case	: X+
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x+1][y][z][m]*mems[x+1][y][z][m];
						}
						ngbL += lesionMembership[x+1][y][z]*lesionMembership[x+1][y][z];
						ngbB += blackHoleMem[x+1][y][z]*blackHoleMem[x+1][y][z];
						// case by case	: X-
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x-1][y][z][m]*mems[x-1][y][z][m];
						}
						ngbL += lesionMembership[x-1][y][z]*lesionMembership[x-1][y][z];
						ngbB += blackHoleMem[x-1][y][z]*blackHoleMem[x-1][y][z];
						// case by case	: Y+
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x][y+1][z][m]*mems[x][y+1][z][m];
						}
						ngbL += lesionMembership[x][y+1][z]*lesionMembership[x][y+1][z];
						ngbB += blackHoleMem[x][y+1][z]*blackHoleMem[x][y+1][z];
						// case by case	: Y-
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x][y-1][z][m]*mems[x][y-1][z][m];
						}
						ngbL += lesionMembership[x][y-1][z]*lesionMembership[x][y-1][z];
						ngbB += blackHoleMem[x][y-1][z]*blackHoleMem[x][y-1][z];
						// case by case	: Z+
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x][y][z+1][m]*mems[x][y][z+1][m];
						}
						ngbL += lesionMembership[x][y][z+1]*lesionMembership[x][y][z+1];
						ngbB += blackHoleMem[x][y][z+1]*blackHoleMem[x][y][z+1];
						// case by case	: Z-
						for (int m=0;m<classes;m++) if (m!=k) {
							ngb += smoothingFactor[m]*mems[x][y][z-1][m]*mems[x][y][z-1][m];
						}
						ngbL += lesionMembership[x][y][z-1]*lesionMembership[x][y][z-1];
						ngbB += blackHoleMem[x][y][z-1]*blackHoleMem[x][y][z-1];
						if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) {
							ngbS = mems[x-1][y][z][k]*mems[x-1][y][z][k] + mems[x+1][y][z][k]*mems[x+1][y][z][k] +
									mems[x][y-1][z][k]*mems[x][y-1][z][k] + mems[x][y+1][z][k]*mems[x][y+1][z][k] +
									mems[x][y][z-1][k]*mems[x][y][z-1][k] + mems[x][y][z+1][k]*mems[x][y][z+1][k] ;
							numL += smoothing*(ngb+ngbS+ngbB);
							lesion_energy += 0.5f*smoothing*(ngb+ngbS+ngbB);
						}else if (objType[k]==BLACKHOLE) {
							ngbS = mems[x-1][y][z][k]*mems[x-1][y][z][k] + mems[x+1][y][z][k]*mems[x+1][y][z][k] +
								mems[x][y-1][z][k]*mems[x][y-1][z][k] + mems[x][y+1][z][k]*mems[x][y+1][z][k] +
								mems[x][y][z-1][k]*mems[x][y][z-1][k] + mems[x][y][z+1][k]*mems[x][y][z+1][k] ;
							numL += smoothing*(ngb+ngbS+ngbB);
							lesion_energy += 0.5f*smoothing*(ngb+ngbS+ngbB);
							numB += smoothing*(ngb+ngbS+ngbL);
							blackHole_energy += 0.5f*smoothing*(ngb+ngbS+ngbL);
						}else if (objType[k]==OUTMASK) 
							numO += smoothing*(ngb+ngbL);
						
						num += smoothing*(ngb+ngbL+ngbB);
							
						energy[k] += 0.5f*smoothing*(ngb+ngbL+ngbB);
						

					}
					
					if (priorCoefficient > 0) {
						float invprior=0;
						float invpriorL=0;
						float invpriorB =0;
						for (int m=0;m<classes;m++) {
							if (m!=k) {
								invprior += similarity[k][m]*shape[m]*shape[m];
							}
							if (objType[k]==LESION || objType[k]==OPTIMIZED_LESION) invpriorL += lesionSimilarity[m]*shape[m]*shape[m];
							if (objType[k]==BLACKHOLE){ 
								invpriorL += lesionSimilarity[m]*shape[m]*shape[m];
								invpriorB += blackHoleSimilarity[m]*shape[m]*shape[m];
							}
						}
						if (hasLesion){ 
							invprior += lesionSimilarity[k]*shape[lesionClass]*shape[lesionClass];
						
							if (objType[lesionClass] == BLACKHOLE)
								invprior += blackHoleSimilarity[k]*shape[lesionClass]*shape[lesionClass];
						}
						num += priorCoefficient*invprior;
						if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) numL += priorCoefficient*invpriorL;
						else if (objType[k]==BLACKHOLE){
							numL += priorCoefficient*invpriorL;
							numB += priorCoefficient*invpriorB;
						}
						else if (objType[k]==OUTMASK) numO += priorCoefficient*invprior;

						energy[k] += 0.5f*priorCoefficient*invprior;
						if (objType[k] == LESION || objType[k] == OPTIMIZED_LESION) lesion_energy += 0.5f*priorCoefficient*invpriorL;
						if (objType[k] == BLACKHOLE) {
							blackHole_energy += 0.5f*priorCoefficient*invpriorB;
							lesion_energy += 0.5f*priorCoefficient*invpriorL;
						}
					}
			
					// invert the result
					if (num>ZERO) num = 1.0f/num;
					else num = INF;
					if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) {
						if (numL>ZERO) numL = 1.0f/numL;
						else numL = INF;
					}else if (objType[k] == BLACKHOLE ){
						if (numL>ZERO) numL = 1.0f/numL;
						else numL = INF;
						if (numB>ZERO) numB = 1.0f/numB;
						else numB = INF;
					}else if (objType[k]==OUTMASK) {
						if (numO>ZERO) numO = 1.0f/numO;
						else numO = INF;
					}
					
					if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) {
						if (relationType==DIRECT) {
							mems[x][y][z][k] = num;
						} else if (relationType == PRIOR) {
							mems[x][y][z][k] = shapeFactor(shape[k]) * num;
						}
						
						if (lesionRelationType == DIRECT){ 
							lesionMembership[x][y][z] = numL;
							lesionRelation[x][y][z] = 1.0f;
						}else if (lesionRelationType == PRIOR ) {
							lesionMembership[x][y][z] = shapeFactor(shape[k]) * numL;
							lesionRelation[x][y][z] = shapeFactor(shape[k]);
						
						}else{
							lesionMembership[x][y][z] = 0;
							lesionRelation[x][y][z] = 0.0f;
						}
						
						if ( (MaxVentDist != 0) && (VentDist[x][y][z] != -1) ) {
							if (VentDist[x][y][z] == 0){ 
								//VentDist[x][y][z] = spread;
								lesionMembership[x][y][z] = spread * lesionMembership[x][y][z] /MaxVentDist;
								lesionRelation[x][y][z] *= (spread / MaxVentDist);
							}else{
								lesionMembership[x][y][z] = VentDist[x][y][z] * lesionMembership[x][y][z] / MaxVentDist ;
								lesionRelation[x][y][z] *= (VentDist[x][y][z] / MaxVentDist);
							}
						}else if ( GMflag && (MaxGMDist != 0) && (GMDist[x][y][z] != -1)  ){	
							if (GMDist[x][y][z] == 0.0f){ 
								lesionMembership[x][y][z] = spread * lesionMembership[x][y][z] / MaxGMDist ;
								lesionRelation[x][y][z] *= (spread / MaxGMDist);
							}else{
								lesionMembership[x][y][z] = GMDist[x][y][z] * lesionMembership[x][y][z] / MaxGMDist ;
								lesionRelation[x][y][z] *= (GMDist[x][y][z] / MaxGMDist);
							}
						}
						
						if(  (MaxBstemDist != 0) && (BSTEMDist[x][y][z] != -1) ) {							
							lesionMembership[x][y][z] = BSTEMDist[x][y][z] * lesionMembership[x][y][z] / MaxBstemDist;
							lesionRelation[x][y][z] *= (BSTEMDist[x][y][z] / MaxBstemDist);
							
						}
						den += lesionMembership[x][y][z];
							
							
						
					} /*else if (objType[k] == BLACKHOLE){
						if (relationType==DIRECT) {
							mems[x][y][z][k] = num;
						} else if (relationType == PRIOR) {
							mems[x][y][z][k] = shapeFactor(shape[k]) * num;
						}
						
						if (lesionRelationType == DIRECT){ 
							lesionMembership[x][y][z] = numL;
							blackHoleMem[x][y][z] = numB;
							lesionRelation[x][y][z] = 1.0f;
						}else if (lesionRelationType == PRIOR ) {
							lesionMembership[x][y][z] = shapeFactor(shape[k]) * numL;
							blackHoleMem[x][y][z] = shapeFactor(shape[k]) * numB;
							lesionRelation[x][y][z] = shapeFactor(shape[k]);
						}else{
							lesionMembership[x][y][z] = 0.0f;
							blackHoleMem[x][y][z] = 0.0f;
							lesionRelation[x][y][z] = 0.0f;
						}
						
						if ( (MaxVentDist != 0) && (VentDist[x][y][z] != -1) ) {
							if (VentDist[x][y][z] == 0){ 
								//VentDist[x][y][z] = spread;
								lesionMembership[x][y][z] = spread * lesionMembership[x][y][z] /MaxVentDist;
								lesionRelation[x][y][z] *= (spread / MaxVentDist);
							}else{
								lesionMembership[x][y][z] = VentDist[x][y][z] * lesionMembership[x][y][z] / MaxVentDist ;
								lesionRelation[x][y][z] *= (VentDist[x][y][z] / MaxVentDist);
							}						
						}else if (GMflag && (MaxGMDist != 0) && (GMDist[x][y][z] != -1)  ){	
							if (GMDist[x][y][z] == 0){ 
								lesionMembership[x][y][z] = spread * lesionMembership[x][y][z] / MaxGMDist ;
								blackHoleMem[x][y][z] =  spread * blackHoleMem[x][y][z] / MaxGMDist ;
								lesionRelation[x][y][z] *= (spread / MaxGMDist);
							}else{
								lesionMembership[x][y][z] = GMDist[x][y][z] * lesionMembership[x][y][z] / MaxGMDist ;
								blackHoleMem[x][y][z] =  GMDist[x][y][z] * blackHoleMem[x][y][z] / MaxGMDist ;
								lesionRelation[x][y][z] *= (GMDist[x][y][z] / MaxGMDist);
							}
						}
						
						if( (MaxBstemDist != 0) && (BSTEMDist[x][y][z] != -1) ) {							
							lesionMembership[x][y][z] = BSTEMDist[x][y][z] * lesionMembership[x][y][z] / MaxBstemDist;
							blackHoleMem[x][y][z] = BSTEMDist[x][y][z] * blackHoleMem[x][y][z] / MaxBstemDist ;
							lesionRelation[x][y][z] *= (BSTEMDist[x][y][z] / MaxBstemDist);
							
						}
						
						den += lesionMembership[x][y][z] + blackHoleMem[x][y][z];
							
					}*/else if (objType[k]==OUTMASK) {
						if (relationType==DIRECT) {
							mems[x][y][z][k] = (1.0f-outlierMembership[x][y][z])*num + outlierMembership[x][y][z]*numO;
							outlierMembership[x][y][z] = outlierMembership[x][y][z]*numO;
						} else if (relationType == PRIOR){
							mems[x][y][z][k] = shapeFactor(shape[k])*((1.0f-outlierMembership[x][y][z])*num + outlierMembership[x][y][z]*numO);
							outlierMembership[x][y][z] = shapeFactor(shape[k])*outlierMembership[x][y][z]*numO;
						}
					
					}/*else if (  hasLesion  && (atlas.getNames()[k].equals("Cerebrum-GM") || atlas.getNames()[k].equals("GM")  )  
								   && (MaxVentDist != 0.0f)&& (VentDist[x][y][z] != -1.0f)  ) {
						if (relationType==DIRECT) {
							mems[x][y][z][k] = num*VentDist[x][y][z]/MaxVentDist;
						} else if (relationType == PRIOR){
							mems[x][y][z][k] = shapeFactor(shape[k])*num*VentDist[x][y][z]/MaxVentDist;
						}
						
					}*/else {
						if (relationType==DIRECT) {
							mems[x][y][z][k] = num;
						} else if (relationType == PRIOR) {
							mems[x][y][z][k] = shapeFactor(shape[k])*num;
						}
					}
					
					if (  atlas.getNames()[k].equals("Cerebrum-GM") || atlas.getNames()[k].equals("GM")  ){
						
						if (distanceFactor!=0){
							if ( ModifedVENTDist[x][y][z] !=-1)  {
								mems[x][y][z][k] *=ModifedVENTDist[x][y][z]/distanceFactor;
								//num += distanceFactor*(1/VentDist[x][y][z]);
								//energy[k] +=  0.5*distanceFactor*(1/ModifedVENTDist[x][y][z]);
								       
							}
						} else if (hasLesion  && (MaxVentDist != 0.0f)&& (VentDist[x][y][z] != -1.0f)  ) {
							
								mems[x][y][z][k] *= (VentDist[x][y][z]/MaxVentDist);

						}
					}
				}
				den += mems[x][y][z][k];
			} 
			// normalization
			if (den>0.0f) {
				for (int k=0;k<classes;k++) {
					mems[x][y][z][k] = mems[x][y][z][k]/den;
					mems_lesion_corrected[x][y][z][k] = mems[x][y][z][k];
					if (objType[k]==LESION || objType[k] == OPTIMIZED_LESION) {
						lesionMembership[x][y][z] = lesionMembership[x][y][z]/den;
						mems_lesion_corrected[x][y][z][k]+=lesionMembership[x][y][z];
					}else if (objType[k] == BLACKHOLE){
						lesionMembership[x][y][z] = lesionMembership[x][y][z]/den;
						blackHoleMem[x][y][z] = blackHoleMem[x][y][z]/den;
					}else if (objType[k]==OUTLIER) {
						outlierMembership[x][y][z] = outlierMembership[x][y][z]/den;
					}else if (objType[k]==OUTMASK) {
						outlierMembership[x][y][z] = outlierMembership[x][y][z]/den;
					}
					// compute energy
					dist += energy[k]*mems[x][y][z][k]*mems[x][y][z][k];
					
				}
			}
			dist += lesion_energy*lesionMembership[x][y][z]*lesionMembership[x][y][z]
			      + blackHole_energy*blackHoleMem[x][y][z]*blackHoleMem[x][y][z];                                                                
								
			
		}
		
		return dist;
    }// computeMembershipsWithPrior
    
    private final float shapeFactor(float sh) {
		return Numerics.min(1.0f, 4.0f*sh);
	}
    
    private final float rayleighDistance(float val, float sig) {
		if (val<=0) return 1.0f;
		else return val*val/sig*sig - 2.0f*(float)Math.log(val/sig) + LOGPI2MPI2;
	}
    	 
    
	
    public final void computeBoundaryDistances() {
		byte[][][] 		tmp = new byte[nx][ny][nz];

		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// find boundaries
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l==1) {
					if (segmentation[x][y][z]!=segmentation[x+i][y+j][z+l]) {
						boundaryDist[x+i][y+j][z+l][segmentation[x][y][z]] = 1;
						boundaryDist[x][y][z][segmentation[x+i][y+j][z+l]] = 1;
					}
				}
			}		
		}
		
		// diffuse the coefficients
		boolean change=true;
		for (int t=1;t<5*limit && change;t++) {
			for (int k=0;k<classes;k++) {
				change = false;
				for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
					tmp[x][y][z] = boundaryDist[x][y][z][k];
					if (segmentation[x][y][z]!=k && boundaryDist[x][y][z][k]==0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
							if (i*i+j*j+l*l==1) {
								if (boundaryDist[x+i][y+j][z+l][k]>0) {
									tmp[x][y][z] = (byte)Numerics.min(tmp[x][y][z],boundaryDist[x+i][y+j][z+l][k]+1);
								}
							}
						}
					}
				}
				for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
					boundaryDist[x][y][z][k] = tmp[x][y][z];
				}
			}
		}
		return;
	}
    /** pre-compute the smoothing and the relations */ 
	public final void computeRelations() {
		short[][][] 			tmp = new short[nx][ny][nz];
		float				num,den;
		int		best;
		String info;
		float factor;
	
		// start with the proper boundaries
		buildRelationList();
		
		// diffuse the coefficients
		for (int t=1;t<5*limit;t++) { 
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				if (relationMap[x][y][z]<classes) {
					     if (relationMap[x-1][y][z]>=classes) tmp[x][y][z] = relationMap[x-1][y][z];
					else if (relationMap[x+1][y][z]>=classes) tmp[x][y][z] = relationMap[x+1][y][z];
					else if (relationMap[x][y-1][z]>=classes) tmp[x][y][z] = relationMap[x][y-1][z];
					else if (relationMap[x][y+1][z]>=classes) tmp[x][y][z] = relationMap[x][y+1][z];
					else if (relationMap[x][y][z-1]>=classes) tmp[x][y][z] = relationMap[x][y][z-1];
					else if (relationMap[x][y][z+1]>=classes) tmp[x][y][z] = relationMap[x][y][z+1];
					else tmp[x][y][z] = relationMap[x][y][z];	
				} else tmp[x][y][z] = relationMap[x][y][z];
			}
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				relationMap[x][y][z] = tmp[x][y][z];
			}
		}
		return;
	}				
	
	/**
	 *  distances and other functions for computing memberships
	 */
	private final float maskDistance(float img, float mean) {
		if (img==mean) return 0.0f;
		else return 1.0f;
	}
	
	/**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean isHomeomorphicSegmentation(short x, short y, short z, short k) {
		// is the new object regular ? not the growing object, the original one!!
		
		// inside the original object ?
		if (segmentation[x][y][z]==k) return true;
		
		boolean [][][] obj = new boolean[3][3][3];
		
		// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[x+i][y+j][z+l]==k) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		// does it change the topology of the object it modifies ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[x+i][y+j][z+l]==segmentation[x][y][z]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = false;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;

		// does it change the topology of a relation between the modified object and its neighbors ?
		int  Nconfiguration = 0;
		short[] lb = new short[26];
		for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[x+i][y+j][z+l]!=k) 
				&& (segmentation[x+i][y+j][z+l]!=segmentation[x][y][z]) ) {
				boolean found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[x+i][y+j][z+l]==lb[n]) { found = true; break; }
				
				if (!found) {
					lb[Nconfiguration] = segmentation[x+i][y+j][z+l];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[x+i][y+j][z+l]==segmentation[x][y][z])
					|| (segmentation[x+i][y+j][z+l]==lb[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = false;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}
		for (int n=0;n<Nconfiguration;n++) {
			// in relation with new object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[x+i][y+j][z+l]==k)
					|| (segmentation[x+i][y+j][z+l]==lb[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = true;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}

		// triplets
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with previous object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[x+i][y+j][z+l]==segmentation[x][y][z])
						|| (segmentation[x+i][y+j][z+l]==lb[n])
						|| (segmentation[x+i][y+j][z+l]==lb[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = false;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with new object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[x+i][y+j][z+l]==k)
						|| (segmentation[x+i][y+j][z+l]==lb[n]) 
						|| (segmentation[x+i][y+j][z+l]==lb[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = true;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}

		// else, it works
		return true;
    }
	
	/**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean isSimpleSkeleton(short x, short y, short z) {
		boolean [][][] obj = new boolean[3][3][3];
		
		// does it change the topology of the skeleton ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[x+i][y+j][z+l]==skeleton[x][y][z]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;

		// else, it works
		return true;
    }
	
	
	/** global score for thinning object k */
	/*
	private final float highMembershipScore(int x, int y, int z, int k) {
		
		return membershipRatio*mems[x][y][z][k] 
				+ fullCurvatureRatio*(1.0f-convexityScore(skeleton,x,y,z,k))
				+ geometricRatio;
	}
	*/
	
	/** global score for growing object k */
	/*
	private final float lowMembershipScore(int x, int y, int z, int k) {
		
		return (membershipRatio*(1.0f-mems[x][y][z][k])
				+ fullCurvatureRatio*convexityScore(skeleton,x,y,z,k)
				+ geometricRatio)/(geometricRatio);
	}
	*/
	
	/** global score ratio for object k */
	/*
	private final float maxMembershipRatio(int x, int y, int z, int k) {
		float max = geometricRatio;
		float cur = mems[x][y][z][k];
		if (cur<geometricRatio) cur = geometricRatio;
		
		for (int m=0;m<classes;m++) if (m!=k) if (mems[x][y][z][m] > max)
			max = mems[x][y][z][m];
		
		return max/cur;
	}
	*/
	private final float maxMembershipRatio(int x, int y, int z, int k, int best) {
		float cur,max;
		if (k==best) return 1.0f;
		else {
			
			if ( (objType[k] == LESION) || (objType[k] == BLACKHOLE) ){ 
				
				if (useLesionMax)
					cur = Numerics.max(mems[x][y][z][k],lesionMembership[x][y][z],blackHoleMem[x][y][z]);
				else
					cur = mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMem[x][y][z];
				
			}else 
				cur = mems[x][y][z][k];
	
			if ((objType[best] == LESION) ||(objType[best] == BLACKHOLE)){
				if (useLesionMax)
					max = Numerics.max(mems[x][y][z][best],lesionMembership[x][y][z],blackHoleMem[x][y][z]);
				else
					max = mems[x][y][z][best]+lesionMembership[x][y][z]+blackHoleMem[x][y][z];
			}else 
				max = mems[x][y][z][best];
				
			// note the added ZERO value is important on both sides!
			return (max+ZERO)/(cur+ZERO);
		}
	}
	
	private final byte bestMembership(int x, int y, int z) {
		byte id = 0;
		float max = 0.0f;
		float mem = 0.0f;
		for (byte m=0;m<classes;m++){ 
			
			if ( (objType[m] == LESION) || ( objType[m] == BLACKHOLE) )
				if (useLesionMax)
					mem = Numerics.max(mems[x][y][z][m] , lesionMembership[x][y][z] , blackHoleMem[x][y][z]);
				else
					mem = mems[x][y][z][m] + lesionMembership[x][y][z] + blackHoleMem[x][y][z];
			else
				mem = mems[x][y][z][m];
				
			if (mem > max) {
				id = m;
				max = mems[x][y][z][m];
			}
		}
		//if (max < lesionMembership[x][y][z] )  id =lesionClass ;
		return id;
	}
		
	/**
	 *	check for convex points / concave points:
	 *  count the number of added boundaries
	 */
	/* 
	private final float convexityScore(short[][][] img, int x, int y, int z, int k) {
		//if (true) return 0.5f;
		
		float boundaries=0;
		float nb6 = 1.0f;
		float nb18 = 0.0f;
		float inb18 = 0.0f;
		float nb26 = 0.0f/3.0f;
		float inb26 = 0.0f/3.0f;
		float total = 6*nb6 + 12*nb18 + 8*nb26;
		
		// 6-C
		if (img[x-1][y][z]!=k) boundaries+=nb6;
		if (img[x+1][y][z]!=k) boundaries+=nb6;
		if (img[x][y-1][z]!=k) boundaries+=nb6;
		if (img[x][y+1][z]!=k) boundaries+=nb6;
		if (img[x][y][z-1]!=k) boundaries+=nb6;
		if (img[x][y][z+1]!=k) boundaries+=nb6;
		
		// 18-C
		if (img[x-1][y-1][z]!=k) boundaries+=nb18;
		else if ( (img[x-1][y][z]!=k) && (img[x][y-1][z]!=k) ) boundaries+=inb18;
		if (img[x-1][y+1][z]!=k) boundaries+=nb18;
		else if ( (img[x-1][y][z]!=k) && (img[x][y+1][z]!=k) ) boundaries+=inb18;
		if (img[x+1][y-1][z]!=k) boundaries+=nb18;
		else if ( (img[x+1][y][z]!=k) && (img[x][y-1][z]!=k) ) boundaries+=inb18;
		if (img[x+1][y+1][z]!=k) boundaries+=nb18;
		else if ( (img[x+1][y][z]!=k) && (img[x][y+1][z]!=k) ) boundaries+=inb18;
		
		if (img[x][y-1][z-1]!=k) boundaries+=nb18;
		else if ( (img[x][y-1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb18;
		if (img[x][y-1][z+1]!=k) boundaries+=nb18;
		else if ( (img[x][y-1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb18;
		if (img[x][y+1][z-1]!=k) boundaries+=nb18;
		else if ( (img[x][y+1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb18;
		if (img[x][y+1][z+1]!=k) boundaries+=nb18;
		else if ( (img[x][y+1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb18;
		
		if (img[x-1][y][z-1]!=k) boundaries+=nb18;
		else if ( (img[x][y][z-1]!=k) && (img[x-1][y][z]!=k) ) boundaries+=inb18;
		if (img[x+1][y][z-1]!=k) boundaries+=nb18;
		else if ( (img[x][y][z-1]!=k) && (img[x+1][y][z]!=k) ) boundaries+=inb18;
		if (img[x-1][y][z+1]!=k) boundaries+=nb18;
		else if ( (img[x][y][z+1]!=k) && (img[x-1][y][z]!=k) ) boundaries+=inb18;
		if (img[x+1][y][z+1]!=k) boundaries+=nb18;
		else if ( (img[x][y][z+1]!=k) && (img[x+1][y][z]!=k) ) boundaries+=inb18;
		
		// 26-C and implied 26-C
		if (img[x-1][y-1][z-1]!=k) boundaries+=nb26;
		else if ( (img[x-1][y-1][z]!=k) && (img[x][y-1][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y-1][z-1]!=k) && (img[x-1][y][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z-1]!=k) && (img[x-1][y-1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z]!=k) && (img[x][y-1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb26;

		if (img[x+1][y-1][z-1]!=k) boundaries+=nb26;
		else if ( (img[x+1][y-1][z]!=k) && (img[x][y-1][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y-1][z-1]!=k) && (img[x+1][y][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z-1]!=k) && (img[x+1][y-1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z]!=k) && (img[x][y-1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb26;

		if (img[x-1][y+1][z-1]!=k) boundaries+=nb26;
		else if ( (img[x-1][y+1][z]!=k) && (img[x][y+1][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y+1][z-1]!=k) && (img[x-1][y][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z-1]!=k) && (img[x-1][y+1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z]!=k) && (img[x][y+1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb26;

		if (img[x-1][y-1][z+1]!=k) boundaries+=nb26;
		else if ( (img[x-1][y-1][z]!=k) && (img[x][y-1][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y-1][z+1]!=k) && (img[x-1][y][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z+1]!=k) && (img[x-1][y-1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z]!=k) && (img[x][y-1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb26;

		if (img[x-1][y+1][z+1]!=k) boundaries+=nb26;
		else if ( (img[x-1][y+1][z]!=k) && (img[x][y+1][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y+1][z+1]!=k) && (img[x-1][y][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z+1]!=k) && (img[x-1][y+1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x-1][y][z]!=k) && (img[x][y+1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb26;

		if (img[x+1][y-1][z+1]!=k) boundaries+=nb26;
		else if ( (img[x+1][y-1][z]!=k) && (img[x][y-1][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y-1][z+1]!=k) && (img[x+1][y][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z+1]!=k) && (img[x+1][y-1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z]!=k) && (img[x][y-1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb26;

		if (img[x+1][y+1][z-1]!=k) boundaries+=nb26;
		else if ( (img[x+1][y+1][z]!=k) && (img[x][y+1][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y+1][z-1]!=k) && (img[x+1][y][z-1]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z-1]!=k) && (img[x+1][y+1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z]!=k) && (img[x][y+1][z]!=k) && (img[x][y][z-1]!=k) ) boundaries+=inb26;

		if (img[x+1][y+1][z+1]!=k) boundaries+=nb26;
		else if ( (img[x+1][y+1][z]!=k) && (img[x][y+1][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x][y+1][z+1]!=k) && (img[x+1][y][z+1]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z+1]!=k) && (img[x+1][y+1][z]!=k) ) boundaries+=inb26;
		else if ( (img[x+1][y][z]!=k) && (img[x][y+1][z]!=k) && (img[x][y][z+1]!=k) ) boundaries+=inb26;

		return (boundaries/total);
	}
	*/
	
	/**
	 *  propagate the skeleton out toward the boundaries
	 *  (geometric progression)
	 */
	public final void propagateGrowing() {	
		float   	val,order=0;
		int 		xi,yj,zl;
		short 		x,y,z,k;
		boolean		isRegular;
		
		// init: reset the boundary tree, the labels
		tree.reset();
		tree.setMinTree();
	
		initGrowingLabels();
		
		while ( isNotEmpty()) {
			
			// get the next value
			x = tree.getFirstX();
			y = tree.getFirstY();
			z = tree.getFirstZ();
			k = tree.getFirstK();
			val = tree.getFirst();
			tree.removeFirst();
			
			if (skeleton[x][y][z]==EMPTY) {
				// test for update
				
				isRegular = isHomeomorphicSegmentation(x,y,z,k);
				
				if (isRegular) {
					// regular point : set the distance
					skeleton[x][y][z] = (byte)k;
					for (int m=0;m<classes;m++) {
						available[x][y][z][m] = false;
					}
					
					// record ordering
					ordering[x][y][z] = val;
					
					// update the volume
					
					//#####LesionTOADS#####
					if (k==lesionClass){
						if (useLesionMax)
							volume[k] += Numerics.max(mems[x][y][z][k],lesionMembership[x][y][z],blackHoleMem[x][y][z]);
						else
							volume[k]+=mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMem[x][y][z];
					}else
						volume[k]+=mems[x][y][z][k];
					
					// update current classif
					segmentation[x][y][z] = (byte)k;
					
					// search for new neighbors *****
					for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						xi=x+i; yj=y+j; zl=z+l;
						if (available[xi][yj][zl][k]) {
							
							order = val + maxMembershipRatio(xi,yj,zl,k,bestMembership(xi,yj,zl));
							available[xi][yj][zl][k] = false;
							tree.addValue(order,xi,yj,zl,k);            
						}
					}
				} else {
					// critical point: back to background
					available[x][y][z][k] = true;
				}
			}
			

			
		}
		
		
		
		return;
	}//propagateGrowing

	
	/**
	 *  secondary function: propagate memberships with a fast marching technique
	 *  and maintain the topology; the distances are fixed but the memberships
	 *  can get lower if needed, introducing structural outliers
	 */
	/*public final void propagateThinning() {
		float   	val,order=0,factor;
		int 		xi,yj,zl;
		short 		x,y,z,k;
		boolean		isRegular;
		int			d2;
		float 		mem;
		int			rel = 0;
		boolean		changeBoundary;
		float[]		minVolume = new float[classes];
		
		// init: reset the boundary tree, the labels
		tree.reset();
		tree.setMaxTree();
		
		initThinningLabels();
		
		
		for (k=0;k<classes;k++) minVolume[k] = Numerics.max(100,volumeRatio*volume[k]);
		
		
		
		// thinning
		val = 0.0f;
		changeBoundary = true;
		while (isNotEmpty() && (val > -limit)) {


			// get the next value
			x = tree.getFirstX();
			y = tree.getFirstY();
			z = tree.getFirstZ();
			k = tree.getFirstK();
			val = tree.getFirst();
			tree.removeFirst();
				
			if (volume[k] > minVolume[k]) {
					
				//if ( (skeleton[x][y][z]==k) && (isSimpleSkeleton(x,y,z)) ) {
				if (skeleton[x][y][z]==k) {
						
					float best = 0.0f;
					byte id = bestMembership(x,y,z);
					boolean connected = false;
					boolean direct = false;
					float bestprior = 0.0f;
					byte next = EMPTY;
					boolean nextconnected = false;
					if (id!=k) {
						for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
							xi=x+i; yj=y+j; zl=z+l;
							if (segmentation[xi][yj][zl]==id) {
								direct = true;
								connected = (skeleton[xi][yj][zl]==id);
							} else if (segmentation[xi][yj][zl]!=k) {
								if (atlas.getTransformedShape(xi,yj,zl)[id]>bestprior) {
									next = segmentation[xi][yj][zl];
									nextconnected = (skeleton[xi][yj][zl]==next);
									bestprior = atlas.getTransformedShape(xi,yj,zl)[id];
								}
							}
						}
					}
					
					if (id==k) {
						// no change if it is already the best
						changeBoundary = false;
					} else if (direct) {
						// the best candidate is a local neighbor
						changeBoundary = true;
					} else if (next!=EMPTY) {
						// the best candidate is closest in this direction
						changeBoundary = true;
						id = next;
						connected = nextconnected;
					} else {
						// nothing is connected: no change possible
						changeBoundary = false;
					}
					
				
					// test for update
					if (changeBoundary) 
						changeBoundary = isHomeomorphicSegmentation(x,y,z,id);
				
					// much better (faster, and more flexible)
					isRegular = (volume[k]>minVolume[k]);
					
					if (isRegular) {
						// reset labels
						available[x][y][z][k] = false;	
					
						// record the ordering
						ordering[x][y][z] = val;
					
						// update the volume
						
						//#########LesionTOADS#########
						if (k==lesionClass){
							if (useLesionMax)
								volume[k] -= Numerics.max(mems[x][y][z][k], lesionMembership[x][y][z], blackHoleMem[x][y][z]);
							else
								volume[k]-= (mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMem[x][y][z]);
						}else
							volume[k]-=mems[x][y][z][k];
						
						// update skeleton
						skeleton[x][y][z] = EMPTY;
						
						// update the classification in good cases
						if (changeBoundary) segmentation[x][y][z] = id;
					
						// search for new or critical neighbors
						for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
							xi=x+i; yj=y+j; zl=z+l; 
							if (available[xi][yj][zl][k]) {
							
								order = val - 1.0f/maxMembershipRatio(xi,yj,zl,k,bestMembership(xi,yj,zl));
								if ( (volume[k] > minVolume[k]) && (order > -limit) ) {
									available[xi][yj][zl][k] = false;
									tree.addValue(order,xi,yj,zl,k);
								}
							}			
						}
					} else {
						// critical point: make available again
						available[x][y][z][k] = true;
					}
				}
			}
		}
		//System.out.print("final volume "+k+": "+volume[k]+"\n");	
		
	}*/
	
	public final void propagateThinning() {
		float   	val,order=0;
		int 		xi,yj,zl;
		short 		x,y,z,k;
		boolean		isRegular;
		boolean		changeBoundary;
		float[]		minVolume = new float[classes];
		// init: reset the boundary tree, the labels
		tree.reset();
		tree.setMaxTree();
		
		
		
		
		long time = System.currentTimeMillis(); 
		if (evolutionMode == NORMAL_DISTANCE){
			init1Ddistance();
			compute1Ddistances();
		}

	
		if (debug) System.out.println("Distances were updated in " + (System.currentTimeMillis() - time));
		initThinningLabels();
		
		
		for (k=0;k<classes;k++) minVolume[k] = Numerics.max(100,volumeRatio*volume[k]);
		
		// thinning
		val = 0.0f;
		changeBoundary = true;
		while (isNotEmpty() && (val > -limit)) {


			// get the next value
			x = tree.getFirstX();
			y = tree.getFirstY();
			z = tree.getFirstZ();
			k = tree.getFirstK();
			val = tree.getFirst();
			tree.removeFirst();
			if (volume[k] > minVolume[k]) {
					
				if (skeleton[x][y][z]==k) {
						
					float best = 0.0f;
					byte id = bestMembership(x,y,z);
					boolean connected = false;
					boolean direct = false;
					float bestprior = 0.0f;
					float closestNeighbor= -1;
					if (evolutionMode == NORMAL_DISTANCE) 
						closestNeighbor = structureDistances[id*nx*ny*nz+x*ny*nz+y*nz+z];
					byte next = EMPTY;
					boolean nextconnected = false;
					if (id!=k) {
						for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
							xi=x+i; yj=y+j; zl=z+l;
							if (segmentation[xi][yj][zl]==id) {
								direct = true;
								connected = (skeleton[xi][yj][zl]==id);
								break;
							} /*else if (segmentation[xi][yj][zl]!=k) {
								if (atlas.getTransformedShape(xi,yj,zl)[id]>bestprior) {
									next = segmentation[xi][yj][zl];
									nextconnected = (skeleton[xi][yj][zl]==next);
									bestprior = atlas.getTransformedShape(xi,yj,zl)[id];
								}*/
							else if (segmentation[xi][yj][zl]!=k){
								if ( evolutionMode ==NORMAL_DISTANCE && ( 
										(structureDistances[id*nx*ny*nz+xi*ny*nz+yj*nz+zl] !=-1 && structureDistances[id*nx*ny*nz+xi*ny*nz+yj*nz+zl] < closestNeighbor)
										|| (structureDistances[id*nx*ny*nz+xi*ny*nz+yj*nz+zl]!=-1 && closestNeighbor == -1))){
									next = segmentation[xi][yj][zl];
									nextconnected = (skeleton[xi][yj][zl]==next);
									closestNeighbor = structureDistances[id*nx*ny*nz+xi*ny*nz+yj*nz+zl];
								} else if (evolutionMode ==NONE && atlas.getTransformedShape(xi,yj,zl)[id]>bestprior) {
									next = segmentation[xi][yj][zl];
									nextconnected = (skeleton[xi][yj][zl]==next);
									bestprior = atlas.getTransformedShape(xi,yj,zl)[id];
								}
							}
						}
					}
					
					if (id==k || next == k) {
						// no change if it is already the best
						changeBoundary = false;
					} else if (direct) {
						// the best candidate is a local neighbor
						changeBoundary = true;
					} else if (next!=EMPTY) {
						// the best candidate is closest in this direction
						changeBoundary = true;
						id = next;
						connected = nextconnected;
					} else {
						// nothing is connected: no change possible
						changeBoundary = false;
					}
					
				
					// test for update
					if (changeBoundary) 
						changeBoundary = isHomeomorphicSegmentation(x,y,z,id);
				
					// much better (faster, and more flexible)
					isRegular = (volume[k]>minVolume[k]);
					
					if (isRegular) {
						// reset labels
						available[x][y][z][k] = false;	
					
						// record the ordering
						ordering[x][y][z] = val;
					
						// update the volume
						
						//#########LesionTOADS#########
						if (k==lesionClass){
							if (useLesionMax)
								volume[k] -= Numerics.max(mems[x][y][z][k], lesionMembership[x][y][z], blackHoleMem[x][y][z]);
							else
								volume[k]-= (mems[x][y][z][k]+lesionMembership[x][y][z]+blackHoleMem[x][y][z]);
						}else
							volume[k]-=mems[x][y][z][k];
						
						// update skeleton
						skeleton[x][y][z] = EMPTY;
						
						// update the classification in good cases
						if (changeBoundary) segmentation[x][y][z] = id;
					
						// search for new or critical neighbors
						for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
							xi=x+i; yj=y+j; zl=z+l; 
							if (available[xi][yj][zl][k]) {
							
								order = val - 1.0f/maxMembershipRatio(xi,yj,zl,k,bestMembership(xi,yj,zl));
								if ( (volume[k] > minVolume[k]) && (order > -limit) ) {
									available[xi][yj][zl][k] = false;
									tree.addValue(order,xi,yj,zl,k);
								}
							}			
						}
					} else {
						// critical point: make available again
						available[x][y][z][k] = true;
					}
					
					
				}
			}
		}
	}

	/**
	 *  propagate the skeleton out toward the boundaries
	 *  (geometric progression)
	 */
	public final void propagateCompetition() {	
		float   	val,order=0;
		int 		xi,yj,zl;
		short 		x,y,z,k;
		byte		m;
		boolean		isRegular;
		
		// init: reset the boundary tree, the labels
		tree.reset();
		tree.setMaxTree();
		
		initCompetitionLabels();
		
		// loop on the boundary
		while ( isNotEmpty()) {
			
			// get the next value
			x = tree.getFirstX();
			y = tree.getFirstY();
			z = tree.getFirstZ();
			// extract k,m: change is k -> m
			k = tree.getFirstK();
			m = (byte)Numerics.floor((byte)(k/classes));
			k = (short)(k-m*classes);
			val = tree.getFirst();
			tree.removeFirst();
			
			if (segmentation[x][y][z]==k) {
				// test for update
				isRegular = isHomeomorphicSegmentation(x,y,z,m);
				
				if (isRegular) {
					// regular point : set the distance
					segmentation[x][y][z] = m;
					skeleton[x][y][z] = m;
					
					// record ordering
					ordering[x][y][z] = val;
					
					// update the volume
					
					//####LesionToads#####
					
					if (m==lesionClass){
						
						volume[m]+=mems[x][y][z][m]+lesionMembership[x][y][z];
						if (objType[m]==BLACKHOLE) volume[m]+=blackHoleMem[x][y][z];
					}else
						volume[m]+=mems[x][y][z][m];
					if (k==lesionClass){
			
						volume[k]-=mems[x][y][z][k]-lesionMembership[x][y][z];
						if (objType[k]==BLACKHOLE) volume[k]-=blackHoleMem[x][y][z];
					}else
						volume[k]-=mems[x][y][z][k];
					
					// search for new neighbors
					for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						xi=x+i; yj=y+j; zl=z+l;
						if (segmentation[xi][yj][zl]==k) {
							//order = (1-curvatureRatio)*(membershipScore(xi,yj,zl,m)-membershipScore(xi,yj,zl,k))
							//		+ curvatureRatio*(1-2*convexityScore(xi,yj,zl,k));
							//order = highMembershipScore(xi,yj,zl,m)-highMembershipScore(xi,yj,zl,k);
							
							//####LesionTOADS
							if (k==lesionClass){
								//order = mems[xi][yj][zl][m]/(Math.max(mems[xi][yj][zl][k],lesionMembership[x][y][z])+ZERO);
								if (objType[k] == LESION)
									order = mems[xi][yj][zl][m]/(mems[xi][yj][zl][k]+lesionMembership[x][y][z]+ZERO);
								else if (objType[k] == BLACKHOLE) 
									order = mems[xi][yj][zl][m]/(mems[xi][yj][zl][k]+lesionMembership[x][y][z]+blackHoleMem[x][y][z]+ZERO);
							}else if (m==lesionClass){
								//order = Math.max(mems[xi][yj][zl][m],lesionMembership[x][y][z])/(mems[xi][yj][zl][k]+ZERO);
								if (objType[k] == LESION)
									order = (mems[xi][yj][zl][m]+lesionMembership[x][y][z])/(mems[xi][yj][zl][k]+ZERO);
								else if (objType[k] == BLACKHOLE)
									order = (mems[xi][yj][zl][m]+lesionMembership[x][y][z]+blackHoleMem[x][y][z])/(mems[xi][yj][zl][k]+ZERO);
							}else
								order = mems[xi][yj][zl][m]/(mems[xi][yj][zl][k]+ZERO);
									
							if (order>0) tree.addValue(order,xi,yj,zl,k);
						}
					}
				} else {
					// critical point: back to background
				}
			}
		}

		return;
	}//propagateCompetition

	/**
	 *	compute the cluster centroids given the distances 
	 */
    public final void computeCentroids() {
        float mem;
        float[] num,den;
		float numL,denL,numB,denB;
		float temp =0.0f;		 
		num = new float[classes];
		den = new float[classes];
        
		for (int c=0;c<nc;c++) {
			numL = denL = 0;
			numB = denB = 0;
			for (int k=0;k<classes;k++) {
				num[k] = 0;
				den[k] = 0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (segmentation[x][y][z] == k){
					
					if ((objType[k] == LESION) || objType[k] == BLACKHOLE) {
						mem = mems[x][y][z][k];
						temp = Numerics.max(mem,lesionMembership[x][y][z],blackHoleMem[x][y][z]);
						if (useField) {
							if (temp == mem){ 
								num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
								den[k] += mem*mem;
							}else if ( temp == lesionMembership[x][y][z]){
									float dist;
									float factor = 1.0f;
									dist = (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c])*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]);
									if (dist>ZERO) factor = (float)(powerlesion[c].lookup(dist, plesion[c])/dist);
									numL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
									denL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z];	
							}else{
								numB += blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
								denB += blackHoleMem[x][y][z]*blackHoleMem[x][y][z];
							}
						} else {
							if (temp == mem){ 
								num[k] += mem*mem*images[c][x][y][z];
								den[k] += mem*mem;
							}else if ( temp == lesionMembership[x][y][z]){
									float dist;
									float factor = 1.0f;
									dist = (images[c][x][y][z]-lesionCentroid[c])*(images[c][x][y][z]-lesionCentroid[c]);
									if (dist>ZERO) factor = (float)(powerlesion[c].lookup(dist, plesion[c])/dist);
									numL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z]*images[c][x][y][z];
									denL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z];	
							}else{
								numB += blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*images[c][x][y][z];
								denB += blackHoleMem[x][y][z]*blackHoleMem[x][y][z];
							}
						}
					}else if (objType[k]==OPTIMIZED_LESION || objType[k]==OPTIMIZED_BLACKHOLE) {
						mem = mems[x][y][z][k];
						temp = Numerics.max(mem,lesionMembership[x][y][z],blackHoleMem[x][y][z]);
						
						float dist=0.0f;
						float factor = 1.0f;
							
						if (useField){ 
							if (temp == mem){ 
								dist = mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
							}else if ( temp == lesionMembership[x][y][z]){
								dist = lesionMembership[x][y][z]*lesionMembership[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
							}else{
								dist = blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
							}
						}else{
							if (temp == mem){ 
								dist = mem*mem*images[c][x][y][z];
							}else if ( temp == lesionMembership[x][y][z]){
								dist = lesionMembership[x][y][z]*lesionMembership[x][y][z]*images[c][x][y][z];
							}else{
								dist = blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*images[c][x][y][z];
							}
						}
												
						if (dist>ZERO) factor = (float)(powerwm[c].lookup(dist, pwm[c])/dist);
							
						if (useField) {
							if (temp == mem){ 
								num[k] += factor*mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
								den[k] += factor*mem*mem;;
							}else if ( temp == lesionMembership[x][y][z]){
								
								numL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
								denL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z];
							}else{
								numB += factor*blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*(field[c][x][y][z]*images[c][x][y][z]);
								denB += factor*blackHoleMem[x][y][z]*blackHoleMem[x][y][z];
							}
						} else {
							if (temp == mem){ 
								num[k] += factor*mem*mem*images[c][x][y][z];
								den[k] += factor*mem*mem;
							}else if ( temp == lesionMembership[x][y][z]){
								numL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z]*images[c][x][y][z];
								denL += factor*lesionMembership[x][y][z]*lesionMembership[x][y][z];
							}else{
								numB += factor*blackHoleMem[x][y][z]*blackHoleMem[x][y][z]*images[c][x][y][z];
								denB += factor*blackHoleMem[x][y][z]*blackHoleMem[x][y][z];
							}
						}
					}else if (objType[k] == OBJECT){
						mem = mems[x][y][z][k];
						if (useField) 
							num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else 
							num[k] += mem*mem*(images[c][x][y][z]);
						den[k] += mem*mem;
					}else if (objType[k]==OUTLIER) {
						mem = mems[x][y][z][k]*(1.0f-outlierMembership[x][y][z]);
						
						if (useField) {
							num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
							den[k] += mem*mem;
						} else {
							num[k] += mem*mem*(images[c][x][y][z]);
							den[k] += mem*mem;
						}
					} else if (objType[k]==OPTIMIZED_CSF) {
							mem = mems[x][y][z][k];
						
							float dist;
							float factor = 1.0f;
							
							if (modality[c] == "T1_MPRAGE"){
								if (useField) dist = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
								else dist = (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
						
								if (dist>ZERO) factor = (float)(powercsf[c].lookup(dist, pcsf[c])/dist);
							
								if (useField) num[k] += factor*mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
								else num[k] += factor*mem*mem*(images[c][x][y][z]);
								den[k] += factor*mem*mem;
							}else {
								if (useField) 
									num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
								else 
									num[k] += mem*mem*(images[c][x][y][z]);
								den[k] += mem*mem;
							}
							
					} else if (objType[k]==OPTIMIZED_GM) {
						mem = mems[x][y][z][k];
						
						float dist;
						float factor = 1.0f;
							
						if (useField) dist = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
						else dist = (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
						
						if (dist>ZERO) factor = (float)(powergm[c].lookup(dist, pgm[c])/dist);
							
						if (useField) num[k] += factor*mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else num[k] += factor*mem*mem*(images[c][x][y][z]);
						den[k] += factor*mem*mem;	
					} else if (objType[k]==OPTIMIZED_WM) {
						mem = mems[x][y][z][k];
						
						float dist;
						float factor = 1.0f;
							
						if (useField) dist = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k])*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
						else dist = (images[c][x][y][z]-centroid[c][k])*(images[c][x][y][z]-centroid[c][k]);
						
						if (dist>ZERO) factor = (float)(powerwm[c].lookup(dist, pwm[c])/dist);
							
						if (useField) num[k] += factor*mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else num[k] += factor*mem*mem*(images[c][x][y][z]);
						den[k] += factor*mem*mem;	
					} else if (objType[k]==RAYLEIGH) {
						mem = mems[x][y][z][k];
						
						if (useField) num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z])*(field[c][x][y][z]*images[c][x][y][z]);
						else num[k] += mem*mem*(images[c][x][y][z])*(images[c][x][y][z]);
						den[k] += mem*mem;	
					} else {
						mem = mems[x][y][z][k];
						
						if (useField) num[k] += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else num[k] += mem*mem*(images[c][x][y][z]);
						den[k] += mem*mem;	
					}
				}
			}
			
			for (int k=0;k<classes;k++) {
				if (centroidMode.equals("none")) {
					if (den[k]>0.0) {
						if (objType[k]==RAYLEIGH) centroid[c][k] = (float)Math.sqrt(num[k]/(2.0f*den[k]));
						else centroid[c][k] = num[k]/den[k];
					} else {
						centroid[c][k] = 0.0f;
					}
				} else if (centroidMode.equals("linked")) {
					float nsum,dsum;
					nsum = num[k];
					dsum = den[k];
					for (int l=0;l<classes;l++) {
						if (l!=k && centroidLink[c][k]==centroidLink[c][l]) {
							nsum += centroidSmoothness*num[l];
							dsum += centroidSmoothness*den[l];
						}
					}
					if (dsum>0.0) {
						centroid[c][k] = nsum/dsum;	
					} else {
						centroid[c][k] = 0.0f;
					}
				} else if (centroidMode.equals("damped")) {
					if (den[k]>0.0) {
						centroid[c][k] = (1.0f-centroidSmoothness)*num[k]/den[k] + centroidSmoothness*centroid[c][k];
					} else {
						// no change: previous value
					}
				} else if (centroidMode.equals("prior")) {
					float nsum,dsum;
					nsum = num[k];
					dsum = den[k];
					
					nsum += centroidSmoothness*1e4f*priorCentroid[c][k];
					dsum += centroidSmoothness*1e4f;
					
					if (dsum>0.0) {
						centroid[c][k] = nsum/dsum;	
					} else {
						centroid[c][k] = 0.0f;
					}
				} else {
					if (den[k]>0.0) {
						centroid[c][k] = num[k]/den[k];
					} else {
						centroid[c][k] = 0.0f;
					}
				}
			}
			boolean linkedBlackHoleCentroid = false;
			if (linkedBlackHoleCentroid && (modality[c].equals("FLAIR")) ){
				
					if (denL>0.0)
						lesionCentroid[c] = (numL+centroidSmoothness*numB)/(denL+centroidSmoothness*denB);
					else
						lesionCentroid[c] = 0.0f;
					
					
					if (denB>0.0)
						blackHoleCentroid[c] = (numB + centroidSmoothness*numL)/(denB+centroidSmoothness*denL);
					else
						blackHoleCentroid[c] = 0.0f;
				
					
				
				
			}else{
				if (denL>0.0) {
			
				lesionCentroid[c] =  0.2f*(numL/denL) + 0.8f*lesionCentroid[c] ;
				//lesionCentroid[c] =  numL/denL;
			} else {
				lesionCentroid[c] = 0.0f;
			}
			if (denB>0.0) {
				blackHoleCentroid[c] =  (numB/denB) ;
			} else {
				blackHoleCentroid[c] = 0.0f;
			}
			}
		}
	        
		return;
    } // computeCentroids
    
	/**
	 *	compute the cluster centroids given the distances 
	 *  and averages with the previous ones for slower
	 *	centroid changes.
	 */
    public final void computeWeightedCentroids(float factor) {
        float num,den,mem;
        
        for (int c=0;c<nc;c++) {
			for (int k=0;k<classes;k++) {
				num = 0;
				den = 0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (segmentation[x][y][z] == k ) {
						mem = mems[x][y][z][k];
						if (useField) num += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else num += mem*mem*(images[c][x][y][z]);
						den += mem*mem;
					}
				}
				if (den>0.0) {
					centroid[c][k] = (1.0f-factor)*num/den + factor*centroid[c][k];
				} else {
					centroid[c][k] = 0.0f;
				}
			}
		}
        return;
    } // computeWeightedCentroids
    
	/**
	 *	compute the cluster centroids given the template
	 *	using template and images intensity, 
	 *	assuming increasing intensity
	 */
    public final void computeConstrainedCentroids() {
        float num,den,mem;
		float prev = -INF;
		float next;
        
		for (int c=0;c<nc;c++) {
			prev = -INF;
			next = INF;
			for (int k=0;k<classes;k++) {
				num = 0;
				den = 0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (segmentation[x][y][z] == k ) {
						mem = mems[x][y][z][k];
						if (useField) num += mem*mem*(field[c][x][y][z]*images[c][x][y][z]);
						else num += mem*mem*(images[c][x][y][z]);
						den += mem*mem;
					}
				}
				if (den>0.0) {
					next = num/den;
				} else {
					next = 0.0f;
				}
				// reset using the previous boundaries
				if (next<centroid[c][k]-0.45f*(centroid[c][k]-prev)) next = centroid[c][k]-0.45f*(centroid[c][k]-prev);
				else if ( (k+1<classes) && (next>centroid[c][k]+0.45f*(centroid[c][k+1]-centroid[c][k])) ) next = centroid[c][k]+0.45f*(centroid[c][k+1]-centroid[c][k]);
				// remember 
				prev = centroid[c][k];
				centroid[c][k] = next;
			}
		}
		return;
    } // computeConstrainedCentroids
    
    
    
	/**
	 *	compute the stddev of the cluster centroids given the distances 
	 */
    public final void computeVariances() {
        float num,den,mem,diff,lnum,lden,ldiff,bnum,bden,lmem;
        for (int c=0;c<nc;c++) {
        	lnum = lden = 0;
        	bnum = bden = 0;
			for (int k=0;k<classes;k++) {
				num = 0;
				den = 0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (segmentation[x][y][z]==k) {
						//lesiontodas
						if ( ( (objType[k] == LESION || objType[k] == OPTIMIZED_LESION) && (mems[x][y][z][k] < lesionMembership[x][y][z])) 
								|| ((objType[k]==BLACKHOLE)&&(mems[x][y][z][k] < blackHoleMem[x][y][z])) ){
							if (lesionMembership[x][y][z] >= blackHoleMem[x][y][z]){
								lmem=lesionMembership[x][y][z];
								if (useField) ldiff = field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c] ;
								else ldiff = images[c][x][y][z]-lesionCentroid[c];
								lnum += (lmem*lmem+1e-30f)*ldiff*ldiff;
								lden += (lmem*lmem+1e-30f);
							}else{
								lmem=blackHoleMem[x][y][z];
								if (useField) ldiff = field[c][x][y][z]*images[c][x][y][z]-blackHoleCentroid[c] ;
								else ldiff = images[c][x][y][z]-blackHoleCentroid[c];
								bnum += (lmem*lmem+1e-30f)*ldiff*ldiff;
								bden += (lmem*lmem+1e-30f);
							}
						}else{
							mem = mems[x][y][z][k];
							if (useField) diff = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
							else diff = (images[c][x][y][z]-centroid[c][k]);
							
							num += (mem*mem+1e-30f)*diff*diff;
							den += (mem*mem+1e-30f);
							//num += diff*diff;
							//den += 1.0f;
						}
					}
				}
				if (den>0.0) {
					stddev[c][k] = (float)Math.sqrt(num/den);
				} else {
					stddev[c][k] = 0.0f;
				}
				
			}
			if (lden>0)
				lesionStDev[c] = (float)Math.sqrt(lnum/lden);
			else 
				lesionStDev[c] = 0.0f;
			if (bden>0)
				blackHoleStdDev[c] = (float)Math.sqrt(bnum/bden);
			else
				blackHoleStdDev[c] = 0.0f;
		}
        
        
        return;
    } // computeVariances
    
    
    public final void computeChannelWeights() {
        float num,den,mem,diff,diffO,diffL,diffB,lmem,bmem;
		float sum;
        
		if (normType.equals("Diagonal Covariance")) { 
			if (debug) System.out.print("Diagonal Channel Weight Computation ....");
			sum = 0;
			for (int c=0;c<nc;c++) {
				num = 0;
				den = 0;
				for (int k=0;k<classes;k++) {
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						mem = mems[x][y][z][k];
						diff = 0; diffO = 0; diffL = 0; 
						if (objType[k]==OBJECT) {
							if (useField) diff = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);//*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
							else diff = (images[c][x][y][z]-centroid[c][k]);//*(images[c][x][y][z]-centroid[c][k]);
						
							num += (mem*mem+1e-30f)*diff*diff;
							den += (mem*mem+1e-30f);
						
						} else if (objType[k]==OUTLIER) {
							if (useField) {
								diff = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);//*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
								diffO = outlier;
							} else {
								diff = (images[c][x][y][z]-centroid[c][k]);//*(images[c][x][y][z]-centroid[c][k]);
								diffO = outlier;
							}
							float omem = outlierMembership[x][y][z];
							num += (mem*mem*(1.0f-omem)*(1.0f-omem)+1e-30f)*diff*diff
									+(mem*mem*omem*omem)*diffO*diffO;
							den += (mem*mem+1e-30f);
						
						} else if ( (objType[k]==LESION) || (objType[k] == OPTIMIZED_LESION)|| objType[k]==BLACKHOLE) {
							if (useField) {
								diff = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);//*(field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
								diffL = (field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]);//*(field[c][x][y][z]*images[c][x][y][z]-lesionCentroid[c]);
								diffB = (field[c][x][y][z]*images[c][x][y][z]-blackHoleCentroid[c]);//*(field[c][x][y][z]*images[c][x][y][z]-blackHoleCentroid[c]);
							} else {
								diff = (images[c][x][y][z]-centroid[c][k]);//*(images[c][x][y][z]-centroid[c][k]);
								diffL = (images[c][x][y][z]-lesionCentroid[c]);//*(images[c][x][y][z]-lesionCentroid[c]);
								diffB = (images[c][x][y][z]-blackHoleCentroid[c]);//*(images[c][x][y][z]-blackHoleCentroid[c]);
							}
							lmem = lesionMembership[x][y][z];
							bmem = blackHoleMem[x][y][z];
							num += (mem*mem+1e-30f)*diff//*diff 
								+(lmem*lmem+1e-30f)*diffL//*diffL
								+(bmem*bmem+1e-30f)*diffB;//*diffB;
							den += (mem*mem+1e-30f)+(lmem*lmem+1e-30f)+(bmem*bmem+1e-30f);
						
						}
					}
				}
				if (num>0.0) {
					channelWeight[c] = den/num;
				} else {
					channelWeight[c] = INF;
				}
				sum += channelWeight[c];
			}
			for (int c=0;c<nc;c++) {
				channelWeight[c] = channelWeight[c]/sum;
			}
		}else{
			for (int c =0; c<nc; c++)
        		channelWeight[c] = (float)1/nc;
		}
		return;
    } // computeImageVariances
    
	
    public final void computCovariance(){
    	
    	float den, denL;
    	for (int c1=0; c1<nc; c1++) for (int c2=0; c2<nc; c2++)
			lesionCov[c1][c2] = 0.0f;
    	for (int k=0; k<classes; k++){
    		for (int c1=0; c1<nc; c1++) for (int c2=0; c2<nc; c2++)
				covariance[k][c1][c2] = 0.0f;
    		den = denL = 0.0f;
    		for(int x=0; x<nx; x++) for(int y=0; y<ny; y++)for(int z=0; z<nz; z++) if (segmentation[x][y][z] == k){
    			if (objType[k] == LESION && mems[x][y][z][k] < lesionMembership[x][y][z]){
    				for (int c1=0; c1<nc; c1++) for (int c2=0; c2<nc; c2++){
        				if (useField)
        					lesionCov[c1][c2] += (lesionMembership[x][y][z]*lesionMembership[x][y][z]+1e-30f)
		                      *(field[c1][x][y][z]*images[c1][x][y][z] - lesionCentroid[c1]) * (field[c2][x][y][z]*images[c2][x][y][z] - lesionCentroid[c2]);
        				else
        					lesionCov[c1][c2] += (lesionMembership[x][y][z]*lesionMembership[x][y][z]+1e-30f)
		                      *(images[c1][x][y][z] - lesionCentroid[c1]) * (images[c2][x][y][z] - lesionCentroid[c2]);
        			}
    				denL += lesionMembership[x][y][z]*lesionMembership[x][y][z]+1e-30f;
    			}else{
    				for (int c1=0; c1<nc; c1++) for (int c2=0; c2<nc; c2++){
        				if (useField)
        					covariance[k][c1][c2] += (mems[x][y][z][k]*mems[x][y][z][k]+1e-30f)
		                      *(field[c1][x][y][z]*images[c1][x][y][z] - centroid[c1][k]) * (field[c2][x][y][z]*images[c2][x][y][z] - centroid[c2][k]);
        				else
        					covariance[k][c1][c2] += (mems[x][y][z][k]*mems[x][y][z][k]+1e-30f)
        				                      *(images[c1][x][y][z] - centroid[c1][k]) * (images[c2][x][y][z] - centroid[c2][k]);
        			}
    				den += mems[x][y][z][k]*mems[x][y][z][k]+1e-30f;
    			}
    		}
    		for (int c1=0; c1<nc; c1++) for (int c2=0; c2<nc; c2++){
				if (den != 0.0f)
					covariance[k][c1][c2] /= den;
				else
					covariance[k][c1][c2] = 0.0f;
				if (objType[k] == LESION)
					if (denL != 0.0f)
						lesionCov[c1][c2] /= denL;
					else 
						lesionCov[c1][c2] = 0.0f;
    		}
    	}
    	
    }
    
    /**
	 *	compute the stddev of the cluster centroids given the distances
	 *	estimate the noise variance (cluster independent)
	 */
    public final void computeNoiseVariance() {
        float num,den,mem,diff;
        
        for (int c=0;c<nc;c++) {
			num = 0;
			den = 0;
			for (int k=1;k<classes;k++) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (segmentation[x][y][z]==k) {
						mem = mems[x][y][z][k];
						
						if (useField) diff = (field[c][x][y][z]*images[c][x][y][z]-centroid[c][k]);
						else diff = (images[c][x][y][z]-centroid[c][k]);
						
						num += (mem*mem+1e-30f)*diff*diff;
						den += (mem*mem+1e-30f);
					
					}
				}
			}
			if (den>0.0) {
				for (int k=0;k<classes;k++)
					stddev[c][k] = (float)Math.sqrt(num/den);
			} else {
				for (int k=0;k<classes;k++)
					stddev[c][k] = 0.0f;
			}
		}
        return;
    } // computeVariances
    
	/**
	 *	output the centroids
	 */
    public final String displayCentroids() {
        String output = "centroids \n";
		for (int c=0;c<nc;c++) {
			output += modality[c] + " : " + centroid[c][0]*Iscale[c];
			for (int k=1;k<classes;k++) output +=" | "+centroid[c][k]*Iscale[c];
			output += "\n";
		}
		return output;
	}
	
	/**
	 *	output the variances
	 */
    public final String displayVariances() {
        String output = "stddevs \n";
		for (int c=0;c<nc;c++) {
			output += modality[c] + " : " + stddev[c][0];
			for (int k=1;k<classes;k++) output +=" | "+stddev[c][k]*stddev[c][k];//*Iscale[c];
			output += "\n";
		}
		return output;
	}
    
    /**
     * outout covariance Matrices
     */
    public final String displayCovariance(){
    	String output = "Covariance Matrices  \n  | " ;
    	float[][] inverse = new float[nc][nc];
    	float det;
    	for (int k =0; k<classes; k++){
       		inverse = covariance[k];
    		output += atlas.getNames()[k] + " : \n";
    		for (int c1=0; c1<nc; c1++){
    			for (int c2=0; c2<nc; c2++){
    				output += (float)inverse[c1][c2] + " ";
    			}
    		output += "\n";
    		}
    	}
    	return output;
    }
	
	/**
	 *	output the lesion centroids
	 */
    public final String displayLesionCentroid() {
        String output = "lesion centroid \n";
		for (int c=0;c<nc;c++) {
			output +=modality[c] +" : " +lesionCentroid[c]*Iscale[c] + " | ";
		}
		output += "\n";
		
		return output;
	}
    
    /**
	 *	output the black hole centroids
	 */
    public final String displayBlackHoleCentroid() {
        String output = "black hole centroid \n";
		for (int c=0;c<nc;c++) {
			output +=modality[c] +" : " +blackHoleCentroid[c]*Iscale[c] + " | ";
		}
		output += "\n";
		
		return output;
	}
    
    /**
     * output the lesion variances
     */
    
    public final String displayLesionVariance() {
        String output = "lesion std: \n";
		for (int c=0;c<nc;c++) {
			output +=modality[c] + " : "+lesionStDev[c]*Iscale[c] + " | ";
		}
		output += "\n";
		
		return output;
	}
    
    /**
     * output the black hole variances
     */
    
    public final String displayBlackHoleVariance() {
        String output = "black hole std: \n";
		for (int c=0;c<nc;c++) {
			output +=modality[c] + " : "+blackHoleStdDev[c]*Iscale[c] + " | ";
		}
		output += "\n";
		
		return output;
	}
    
    public final String displayChannelWeights(){
    	String output = "New Channel Weights: \n";
    	
    	
    	
    	for (int c=0; c<nc; c++)
    		output += modality[c] + " : " + channelWeight[c] + " | ";
    	output += "\n";
    	
    	return output;
    }
    /**
     * output lesion covariance matrix
     */
    
    /**
     * outout covariance Matrices
     */
    public final String displayLesionCovariance(){
    	String output = "Lesion Covariance Matrix  \n  | " ;
    	float[][] lesionInverse = new float[nc][nc];
    	float det;
    	lesionInverse = lesionCov;
    		for (int c1=0; c1<nc; c1++){
    			for (int c2=0; c2<nc; c2++)
    				output +=(float)lesionInverse[c1][c2]+ "   ";
    		output += "\n";
    		}
    		output += " |";
    	
    	return output;
    }
	
    public final float[][][][] generateCruiseInputs(){
    	
    	float[][][][]	Mems = new float[4][nx][ny][nz];
		
		// CSF
		for (int k=0;k<classes;k++) {
			if (atlas.getNames()[k].equals("SulcalCSF") 
				|| atlas.getNames()[k].equals("Sulcal-CSF")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) 
					Mems[0][x][y][z] = mems[x][y][z][k] ;
			}
		}
		
		// GM
		for (int k=0;k<classes;k++) {
			if (atlas.getNames()[k].equals("CorticalGM") 
				|| atlas.getNames()[k].equals("Cerebrum-GM") 
				|| atlas.getNames()[k].equals("CerebralGM")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) 
						Mems[1][x][y][z] = mems[x][y][z][k];
			}
		}

		// WM (everything else except background and cerebellum / brainstem
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			Mems[2][x][y][z] = (lesionMembership[x][y][z] + blackHoleMem[x][y][z]);
			Mems[3][x][y][z] = 0.0f;
		}
		for (int k=0;k<classes;k++) {
			if (!atlas.getNames()[k].equals("SulcalCSF") 
				&& !atlas.getNames()[k].equals("Sulcal-CSF")
				&& !atlas.getNames()[k].equals("CorticalGM") 
				&& !atlas.getNames()[k].equals("Cerebrum-GM") 
				&& !atlas.getNames()[k].equals("CerebralGM")
				&& !atlas.getNames()[k].startsWith("Cerebell")
				&& !atlas.getNames()[k].equals("Brainstem")
				&& !atlas.getNames()[k].equals("Background")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					Mems[2][x][y][z] += (mems[x][y][z][k]);
					if (segmentation[x][y][z]==k) Mems[3][x][y][z] = 1.0f;
				}
			}
		}
			
		return Mems;
	} // generateCruiseOutputs

    public final float[][][][] generateDuraRemovalOutputs(){
    	
    	float[][][][]	Mems = new float[4][nx][ny][nz];
    			
    			// CSF
    			for (int k=0;k<classes;k++) {
    				if (atlas.getNames()[k].equals("SulcalCSF") 
    					|| atlas.getNames()[k].equals("Sulcal-CSF")) {
    					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) 
    						Mems[0][x][y][z] = mems[x][y][z][k] ;
    				}
    			}
    			
    			
    			// WM (everything else except background and cerebellum-GM 
    			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    				Mems[1][x][y][z] = 0.0f;
    				Mems[2][x][y][z] = (lesionMembership[x][y][z] + blackHoleMem[x][y][z]);
    				Mems[3][x][y][z] = 0.0f;
    			}
    			for (int k=0;k<classes;k++) {
    				if (!atlas.getNames()[k].equals("SulcalCSF") 
    					&& !atlas.getNames()[k].equals("Sulcal-CSF")
    					&& !atlas.getNames()[k].equals("CorticalGM") 
    					&& !atlas.getNames()[k].equals("Cerebrum-GM") 
    					&& !atlas.getNames()[k].equals("CerebralGM")
    					&& !atlas.getNames()[k].equals("Cerebellum-GM")
    					&& !atlas.getNames()[k].equals("CerebellarGM")
    					&& !atlas.getNames()[k].equals("Background")) {
    					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    						Mems[2][x][y][z] += (mems[x][y][z][k]);
    						if (segmentation[x][y][z]==k) Mems[3][x][y][z] = 1.0f;
    					}
    				}
    			}
    			
    			// GM (Cerebellum and Cortical)
    			for (int k=0;k<classes;k++) {
    				if (atlas.getNames()[k].equals("CorticalGM") 
    					|| atlas.getNames()[k].equals("Cerebrum-GM") 
    					|| atlas.getNames()[k].equals("CerebralGM")
    					|| atlas.getNames()[k].equals("Cerebellum-GM")
    					|| atlas.getNames()[k].equals("CerebellarGM")) {
    					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    							Mems[1][x][y][z] += mems[x][y][z][k];
    							//include all the cerebellum in WM mask
    							if (atlas.getNames()[segmentation[x][y][z]].equals("Cerebellum-GM")
    								|| atlas.getNames()[segmentation[x][y][z]].equals("CerebellarGM)"))
    								Mems[3][x][y][z] = 1.0f;
    					}
    				}
    			}

    			return Mems;
    		} // generateDuraRemovalOutputs
    
    
	/**
	 *  analyze the template classes to make relation map 
	 */
	public final void buildRelationList() {
		// find the boundaries between labels
		int MAXLB=500;
		short[][] lbs = new short[MAXLB][];
		short Nlb=0;
		short[] newlb = new short[MAXLB];
		int Nb;
		boolean isFound,isFoundLine,newValue;
		// start with the single class relations
		for (short m=0;m<classes;m++) {
			lbs[Nlb] = new short[1];
			lbs[Nlb][0] = m;
			Nlb++;
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// find boundaries
			for (int n=0;n<MAXLB;n++) newlb[n] = EMPTY;
			Nb = 0;
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l==1) {
					if (segmentation[x][y][z]!=segmentation[x+i][y+j][z+l]) {
						newValue = true;
						for (int b=0;b<Nb;b++) 
							if (segmentation[x+i][y+j][z+l] == newlb[b])
								{ newValue = false; break; }
						if (newValue) {
							newlb[Nb] = segmentation[x+i][y+j][z+l];
							Nb++;
						}
					}
				}
			}
			if (newlb[0]!=EMPTY) {
				newlb[Nb] = segmentation[x][y][z]; 
				Nb++; 
				// check if boundary unique
				isFoundLine = false;
				for (short n=0;n<Nlb;n++) {
					if (lbs[n].length!=Nb) {
						isFoundLine = false;
					} else {
						isFoundLine = true;
						for (int b=0;b<Nb;b++) {
							isFound = false;
							for (int c=0;c<lbs[n].length;c++)
								if (newlb[b]==lbs[n][c]) { isFound = true; break; }
							if (!isFound) isFoundLine = false;
						}
					}
					if (isFoundLine) {
						// same labels
						relationMap[x][y][z] = n;
						break;
					}
				}
				if (!isFoundLine) {
					// create a new entry
					lbs[Nlb] = new short[Nb];
					for (int b=0;b<Nb;b++) lbs[Nlb][b] = newlb[b];
					relationMap[x][y][z] = Nlb;
					Nlb++;
				}
			} else {
				// not a boundary: find the object label
				relationMap[x][y][z] = segmentation[x][y][z];
			}
		}
		// convert into boolean relations
		Nrelations = Nlb;
		relationList = new boolean[Nlb][classes];
		for (int n=0;n<Nlb;n++) {
			for (int m=0;m<classes;m++) relationList[n][m] = false;
			
			for (int b=0;b<lbs[n].length;b++) {
				if (lbs[n][b]!=EMPTY) relationList[n][lbs[n][b]] = true;
			}
		}
		// for the first relations: list all related objects
		for (int k=0;k<classes;k++) {
			for (int n=classes;n<Nlb;n++) {
				if (relationList[n][k]) {
					for (int m=0;m<classes;m++) {
						if (relationList[n][m]) relationList[k][m] = true;
					}
				}
			}
		}
		lbs = null;
		newlb = null;

		return;
	}//updateRelationList

	/**
	 *  show the relation list
	 */
	public final void displayRelationList() {
		// output the list of boundaries
		System.out.println("relation list: \n");
		System.out.println(" : (");
		for (int m=0;m<classes;m++) System.out.println(""+templateLabel[m]+" ");
		System.out.println(")\n");
		for (int n=0;n<Nrelations;n++) {
			System.out.println(""+n+": (");
			//for (int m=0;m<classes+pairs;m++) 
			for (int m=0;m<classes;m++) 
				if (relationList[n][m]) System.out.println(" 1 ");
				else System.out.println(" 0 ");
			System.out.println(")\n");
		}

		System.out.println("intensity relations: \n");
		for (int k=0;k<classes;k++) {
			System.out.println(templateLabel[k]+" : ");
			for (int m=0;m<classes;m++) if (relationList[k][m]) System.out.println(templateLabel[m]+" ");
			System.out.println("\n");
		}
		System.out.println("relation labels: \n");
		for (int n=classes;n<Nrelations;n++) {
			System.out.println(n+" > ");
			for (int k=0;k<classes;k++) if (relationList[n][k]) System.out.println(templateLabel[k]+" ");
			System.out.println("\n");	
		}
		return;
	} //displayRelationList

	/**
	 *	compute the classification from the labels 
	 */
    public final int updateClassificationFromLabels() {
        int classDiff = 0;
		byte prev;
        
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			prev = classification[x][y][z];
			classification[x][y][z] = templateLabel[segmentation[x][y][z]];
			
			if (prev != classification[x][y][z]) classDiff++;		
       }
        return classDiff;
	}

	/**
	 *	compute the classification from the labels 
	 */
    public final int updateClassificationFromSkeleton() {
        int classDiff = 0;
		byte prev;
        
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			prev = classification[x][y][z];
			classification[x][y][z] = 0;
			if (skeleton[x][y][z]!=EMPTY)
				classification[x][y][z] = templateLabel[skeleton[x][y][z]];
			if (prev != classification[x][y][z]) classDiff++;		
       }
	   return classDiff;
	}


	/**
	 *	compute the classification from the memberships
	 */
    public final int updateClassificationFromMemberships() {
        int classDiff = 0;
		byte prev;
		int best;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			prev = classification[x][y][z];
			classification[x][y][z] = 0;
			best = 0;
			for (int k=1;k<classes;k++) if (mems[x][y][z][k]>mems[x][y][z][best]) {
				best = k;
			}
			if ( lesionMembership[x][y][z] > mems[x][y][z][best])
				classification[x][y][z] = templateLabel[lesionClass];
			else if (mems[x][y][z][best]>0) classification[x][y][z] = templateLabel[best];
			if (prev != classification[x][y][z]) classDiff++;		
       }
	   return classDiff;
	}
    
    public final int updateLesionClassifcation(){
    	
    	int lesionVolume =0;
    	   	
    	for (int x =0; x<nx; x++) for (int y =0; y<ny; y++) for (int z=0; z<nz; z++) if (hasLesion)
    		if (  (classification[x][y][z] == templateLabel[lesionClass])
    			 && ( (lesionMembership[x][y][z] > mems[x][y][z][lesionClass]) || (blackHoleMem[x][y][z] > mems[x][y][z][lesionClass]) ) ) {
    			
    			if (lesionMembership[x][y][z] < blackHoleMem[x][y][z])
    				lesionClassification[x][y][z] = (short)2;
    			else
    				lesionClassification[x][y][z] = (short)1;
    			lesionVolume++;
    			
    		}else
    			lesionClassification[x][y][z] = (short)0;
    	else
    		lesionClassification[x][y][z] = (short)0;
    	    	
    	return lesionVolume;
    }
    
 public final int updateLesionClassifcationformMemberships(){
    	
    	int lesionVolume =0;
    	   	
    	int gm = 0;
		for (int k=0; k<classes; k++)
			if (atlas.getNames()[k].equalsIgnoreCase("Cerebellum-GM")) gm =k;
    	for (int x =0; x<nx; x++) for (int y =0; y<ny; y++) for (int z=0; z<nz; z++) if (hasLesion)
    		if ( ( classification[x][y][z] == templateLabel[lesionClass]) && ( mems[x][y][z][lesionClass] < (lesionMembership[x][y][z] + mems[x][y][z][gm]))){
    			 
    			lesionClassification_mem[x][y][z] = (short)1;
    			lesionVolume++;
    		}else
    			lesionClassification_mem[x][y][z] = (short)0;
    	else
    		lesionClassification_mem[x][y][z] = (short)0;
    	return lesionVolume;
    }
 
 	public final void finalizeSegmentation(float voxelSize){
 		int minVoxelNum = Math.round(minLesionVolume/voxelSize);
 		boolean[][][] les = new boolean[nx][ny][nz];
 		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			les[x][y][z] = (lesionClassification[x][y][z]>0);
 		}
 		//Fill holes
 		les = ObjectProcessing.removeHoles(les, nx, ny, nz, 6);
 		//Prone Lesions
 		if (prone_lesions){
 		    int[][][] lb = ObjectProcessing.connected18Object3D(les, nx, ny, nz);
 		    int[] list = ObjectProcessing.listLabels(lb, nx, ny, nz);
 		    int[] vol = new int[list.length];
 		    for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (lb[x][y][z]>0) {
 		       for (int l=0;l<list.length;l++) if (list[l]==lb[x][y][z]) {
 		           vol[l]++;
 		       }
 		    }
 		    for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (lb[x][y][z]>0) {
 		       for (int l=0;l<list.length;l++) if (list[l]==lb[x][y][z]) {
 		           if (vol[l]<minVoxelNum) les[x][y][z] = false;
 		       }
 		    }
 		    if (debug) System.out.println("Removing lesions less than "+ minLesionVolume+" mm3 ...");
 		}
 		int wmClass=0;
 		for (int k=0; k<classes; k++) 
 			if (atlas.getNames()[k].equals("Cerebrum-WM")|| atlas.getNames()[k].equals("CerebralWM")|| atlas.getNames()[k].equals("OuterVentricular-WM")){
 				wmClass=k;
 				break;
 			}
 			
 		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			if (les[x][y][z])
				lesionClassification[x][y][z] = 1;
			else
				lesionClassification[x][y][z] = 0;
			if (lesionClassification[x][y][z]==1 & classification[x][y][z]!=templateLabel[wmClass]){
				lesionMembership[x][y][z]+=mems[x][y][z][segmentation[x][y][z]];
				mems[x][y][z][segmentation[x][y][z]]=0.0f;
				classification[x][y][z] = templateLabel[wmClass];
			}
		}
 	}
	/** 
	 *	export membership functions 
	 */
	public final float[][][][] exportMemberships() {
		float[][][][]	Mems = new float[classes][nx][ny][nz];
		
		for (int k=0;k<classes;k++) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				Mems[k][x][y][z] = mems[x][y][z][k];
			}
		}
		return Mems;
	} // exportMemberships
	
	public final byte[][][] exportClassification() {
		byte[][][]	Classif = new byte[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			//Classif[x][y][z] = templateLabel[segmentation[x][y][z]];
			Classif[x][y][z] = classification[x][y][z];
		}
		if (debug) System.out.println("Classification from Label");
		return Classif;
	} // exportClassification
	
	public final short[][][] exportLesionClassification() {
		
		return lesionClassification;
	} 
	
	public final short[][][] exportLesionClassificationformMembership() {
		
		return lesionClassification_mem;
	} 
	
	public final byte[][][] exportSkeleton() {
		byte[][][]	Classif = new byte[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (skeleton[x][y][z]!=EMPTY) Classif[x][y][z] = templateLabel[skeleton[x][y][z]];
			else Classif[x][y][z] = 0;
		}
		return Classif;
	} // exportSkeleton
	public final byte[][][] exportSegmentation() {
		byte[][][]	Classif = new byte[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (skeleton[x][y][z]!=EMPTY) Classif[x][y][z] = templateLabel[segmentation[x][y][z]];
			else Classif[x][y][z] = 0;
		}
		return Classif;
	}
	
	public final float[][][] exportOrdering() {
		float[][][]	Order = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			Order[x][y][z] = ordering[x][y][z];
		}
		return Order;
	} // exportOrdering
	
	public final float[][][] exportLesionMembership() {
		float[][][]	Lesion = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			Lesion[x][y][z] = lesionMembership[x][y][z];
		}
		return Lesion;
	} // exportLesionMembership
	
	public final float[][][] exportOutlierMembership() {
		float[][][]	Outlier = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			Outlier[x][y][z] = outlierMembership[x][y][z];
		}
		return Outlier;
	} // exportOutlierMembership
	
	public final byte[][][][] exportDistances(){
		byte[][][][] distances = new byte[3][nx][ny][nz];
		distances[0] = VentDist;distances[1]=GMDist;distances[2]=BSTEMDist;
		return distances;
	}
	public final byte[][][] exportRelationMap() {
		byte[][][]	Relations = new byte[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			Relations[x][y][z] = (byte)relationMap[x][y][z];
		}
		return Relations;
	} // exportRelationMap
	
	private final void updateMaxmemClassification(){
  			
  		
		byte prev;
		byte best;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			prev = max_mem_segmentation[x][y][z];
			max_mem_segmentation[x][y][z] = 0;
			best = 0;
			for (int k=1;k<classes;k++) if (mems[x][y][z][k]>mems[x][y][z][best]) {
				best = (byte)k;
			}
			if ( hasLesion && mems[x][y][z][lesionClass]+lesionMembership[x][y][z] > mems[x][y][z][best])
				max_mem_segmentation[x][y][z] = (byte)lesionClass;
			else if (mems[x][y][z][best]>0) 
				max_mem_segmentation[x][y][z] = best;
				
       }
	}
	
   
	/**  generate atlas image from information
	 */
    final public byte[][][] generateClassification() {
		float dist,max,count;
		int best=0;
		byte[][][] img = new byte[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			// compute each class probability : attribute the highest
			max = 0; best = -1;
			for (int k=0;k<classes;k++) {
				if (mems[x][y][z][k]>max) {
					best = k;
					max = mems[x][y][z][k];
				}
			}
			if ( lesionMembership[x][y][z] > max)
				img[x][y][z] = templateLabel[lesionClass];
			else if (best>-1) img[x][y][z] = templateLabel[best];
			else img[x][y][z] = 0;
		}
		return img;
	}
	

}
