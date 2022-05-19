package edu.jhmi.rad.medic.methods;
 
import java.io.*;
import java.util.*;

import edu.jhmi.rad.medic.libraries.*;
import edu.jhmi.rad.medic.utilities.*;
import edu.jhmi.rad.medic.structures.*;

/**
 *
 *  This class handles full structure atlas information:
 *	shape, topology, relations, etc.
 *
 *	@version    June 2006
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class DemonToadDeformableAtlas {
	private static final String cvsversion = "$Revision: 1.5 $"; 
	public static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", ""); 

	public static String get_version() {
	   return revnum;
	}


	// structures: basic information	
	private		int					classes;			// number of strcutures in the atlas
	private		String[]			name;				// their names
	private		byte[]				label;				// their labels
	private		String[]			topology;			// their topology type
	private		String				atlasFile;			// the atlas file
	
	// atlas quantities
	private 	int 				nix,niy,niz; 			// image dimensions
	private 	float 				rix,riy,riz; 			// image resolutions
	private		int					orient,orix,oriy,oriz;		// image orientations
	private		float				x0i,y0i,z0i;		// the center of the image
	
	// spatial transformations
	private		float[]				transform;		// the transform parameters to get into image space
	private		float[][]			rotation;		// the associated rotation matrix
	private		float[][][]			shapeTransform; // the global transform matrix (XI = A XP) for each shape
	private		float[][]			multitransform;		// the transform parameters to get into image space
	private		float[][][]			multirotation;		// the associated rotation matrix
	private		int					Nd;				// transform dimension
	private		ParametricTransform		transformModel;	// the type of transform (from possible ones below)
	private		int					transformMode;	// the type of transform (from possible ones below)
	private static final	int   	NONE = 0;
	private static final	int   	SINGLE = 1;
	private static final	int   	MULTIPLE = 2;
	private static final	int   	DEFORMABLE = 3;
	
	// shape maps
	private		boolean[]			hasShape;			// flag to notify which structures use a shape atlas
	private		float[][][][]		shape;				// the shape images
	private		String[]			shapeFile;			// the shape filenames
	private		float				shapeScale;			//	the slope of the sigmoid prior based on distance functions
	private		int					labelSamples;		// number of samples in the shape/distasnce/contact/direction model
	private		int[]				minx,miny,minz;		// the lowest image coordinate with non zero prior for each structure
	private		int[]				maxx,maxy,maxz;		// the highest image coordinate with non zero prior for each structure
	private		boolean[]			registeredShape;
	private		float				shapeSlope = 2.0f;	// slope of shape prior enhancement in deformable registration
	
	// topology template
	private		boolean				hasTopology;	// flag to notify if there is a topology template
	private		byte[][][]			template;		// the topology template for the segmentation
	private		String				templateFile;	// the topology template for the segmentation
	private		int					nax,nay,naz;	// the topology atlas dimensions
	private 	float 				rax,ray,raz; 	// image resolutions
	private		float				x0a,y0a,z0a;	// the center of the topology image
	
	// intensity models
	private		boolean[]			hasIntensity;	// which intensity models are available
	private		float[][]			intensity;		// the intensity models, normalised between 0 and 1
	private		int					intensitySamples;	// the number of samples used for the intensity model
	
	// intensity variance models
	private		boolean[]			hasIntensityVariance;	// which intensity models are available
	private		float[][]			intensityVariance;		// the intensity models, normalised between 0 and 1
	
	// intensity models in use..
	public		static final int	T1_SPGR = 0;
	public		static final int	T2 = 1;
	private		static final int	FLAIR = 2;
	public		static final int	T1_MPRAGE = 3;
	private		static final int	T1_RAW = 4;
	public		static final int	PD = 5;
	private		static final int	PDFSE = 6;
	private		static final int	DIR = 7;
	private		static final int	INTENSITY = 8;	// the number of possible intensity models
	
	// modality weighting
	private		float[][]			modweight;		// the weighting for modality and obj / lesion (temporary)
	private		static final int	OBJTYPES = 2;	// number of different object types (for now obj and lesion)


	private		float[][]			optimizedFactor;	// the factors for optimized distances of CSF/GM/WM
	private	static final 	int		OPTIMIZED = 4;	
	
	// lesions models
	private		boolean				hasLesions;		// whether there is a lesion model
	private		float[]				lesion;			// the lesion model for each intensity
	private		float[]             blackHole;		//the vlack hole model for each intensity
	
	// for registration
	private		float[][][][]		mems;			// memberships to align
	private		float[][]			centroid;		// corresponding intensity clusters
	private		float				priorCoefficient,priorScale;
	
	// non-linear registration
	private 	DemonsToadsAtlasWarping		demons;
	
	// Levenberg Marquardt parameters
	private		float		chisq,oldchisq;			// the chi square error value
	private		float		lambda	=	1.0f;		// the Levenberg-Marquardt parameter for Levenberg Marquardt estimation
	//private		float		cost, norm;
	//private		float[]		gradient, hessian;			// the coefficients to compute 
	//private		float[]		normgradient, normhessian;			// the coefficients to compute 
	private		float		chiPrecision = 1e-3f;	// the lower limit on the chi square value 
	private		float		lfactor = 1.5f;			// the multiplicative factor for the adaptive term
	private		int			itSupp = 10;				// maximum of steps if the cost function is not improving
	private		int			itMax,itPlus,Nturn;		// counters for various loops
	private static final	float   INIT_LAMBDA = 1;
	private		float		minEdiff = 1e-6f;		// the minimum variation of energy to require a better alignment
	private		float		minLambda = 0.001f;		// the minimum variation of energy to require a better alignment
	private		int			subsample = 3;			// scale for the registration: just subsample the volume 
	private		int 		levels = 1;				// number of image scales (pyramid)
	private		int			offset = 0;				// offset used in subsampling (cyclic)
	private 	boolean		precompute=true;
	
	private static final byte	MAXCORR = 1;
	private static final byte	SYMMAXCORR = 2;
	private				byte	costFunction = MAXCORR;
	
	private static final byte	X = 0;
	private static final byte	Y = 1;
	private static final byte	Z = 2;
	private static final byte	T = 3;
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=false;
	
	
	// numerics
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;

	/**
	 *	constructor: load the atlas information from a file.
	 *	<p>
	 *	The atlas files follow a certain template; 
	 *  the separator between numbers is a tab, not a space.
	 */
	public DemonToadDeformableAtlas(String filename) {
		
		transformMode = NONE;
		Nd = 0;
		transform = null;
		multitransform = null;
		
		labelSamples = 0;
		intensitySamples = 0;
		loadAtlas(filename);
		
	}
	
	/**
	 *	constructor: create an empty atlas.
	 */
	public DemonToadDeformableAtlas(int Nc_) {
		
		classes = Nc_;
		
		transformMode = NONE;
		Nd = 0;
		transform = null;
		multitransform = null;
		
		labelSamples = 0;
		intensitySamples = 0;
		
		// allocate everiything
		name = new String[classes];
		label = new byte[classes];
		topology = new String[classes];
			
		hasShape = new boolean[classes];
		for (int n=0;n<classes;n++) hasShape[n] = false;
		shape = new float[classes][][][];
		minx = new int[classes];
		miny = new int[classes];
		minz = new int[classes];
		maxx = new int[classes];
		maxy = new int[classes];
		maxz = new int[classes];
		hasIntensity = new boolean[INTENSITY];
		for (int i=0;i<INTENSITY;i++) hasIntensity[i] = false;
		intensity = new float[INTENSITY][classes];
		modweight = new float[INTENSITY][OBJTYPES];
		for (int i=0;i<INTENSITY;i++) for (int j=0;j<OBJTYPES;j++) modweight[i][j] = 1.0f;					
	}
	
	/** clean-up: destroy membership and centroid arrays */
	public final void finalize() {
		shape = null;
		template = null;
		intensity = null;
		System.gc();
	}
	
	/** link the variables */
	final public int 		getNumber() { return classes; }
	final public String[] 	getNames() { return name; }
	final public String 	getTemplateFile() { return templateFile; }
	final public String[] 	getTopology() { return topology; }
	final public byte[]		getLabels() { return label; }
	final public int getNameID (String txt){
		
		for (int k=0; k<classes; k++)
			if (name[k].equals(txt)) return k;
		return -1;
	}
final public int getNameLabel (String txt){
		
		for (int k=0; k<classes; k++)
			if (name[k].equals(txt)) return label[k];
		return -1;
	}
	final public void		setName(int id, String txt) { name[id] = txt; }
	final public void		setLabel(int id, byte val) { label[id] = val; }
	final public void		setTopology(int id, String txt) { topology[id] = txt; }
	
	final public void setShapeScale(float s_) { shapeScale = s_; }
	
	final public byte[][][] 	getTemplate() { return template; }
	final public float[][][] 	getShape(int n) { return shape[n]; }
	final public float[][][][] 	getShapes() { return shape; }
	
	/*final public void	setTemplate(byte[][][] t) { template = t; }
	final public void	setShapes(float[][][][] s) { shape = s; }
	final public void	setShape(int n, float[][][] s) { shape[n] = s; }*/
	final public void	setTemplate(byte[][][] t) { template = t; }
	final public void	setTemplate(byte[][][] t, int[] dim, float[] res) { 
		
		
		if (dim.length != 3)
			System.out.println("Warning: template dimension is unchanged!");
		else{
			nax=dim[0];
			nay=dim[1];
			naz=dim[2];
		}
		if (res.length != 3)
			System.out.println("Warning: template resolution is unchanged!");
		else{
			rax = res[0];
			ray = res[1];
			raz = res[2];
		}
		template = t; 
	}
	final public void	setShapes(float[][][][] s) { shape = s; }
	final public void	setShapes(float[][][][] s, int[] dim, float[] res) { 
		shape=null;
		if (dim.length != 3)
			System.out.println("Warning: shape atlas dimension is unchanged!");
		else{
			System.out.println("Atlas numbers are "+s.length+ " * "+s[0].length + " * " + s[0][0].length);
				/*nix= dim[0];
				niy = dim[1];
				niz = dim[2];
				x0i = nix/2.0f;
				y0i = niy/2.0f;
				z0i = niz/2.0f;*/
				nax= dim[0];
				nay = dim[1];
				naz = dim[2];
				x0a = nax/2.0f;
				y0a = nay/2.0f;
				z0a = naz/2.0f;
		}
		if (res.length != 3)
			System.out.println("Warning: shape atlas resolution is unchanged!");
		else{
			for (int i=0; i<s.length; i++){
				rax=res[0];
				ray=res[1];
				raz=res[2];
			}
		}
		shape = new float[classes][dim[0]][dim[1]][dim[2]];
		for (int k=0; k<classes; k++) for (int x=0; x<dim[0]; x++)for (int y=0; y<dim[1]; y++) for (int z=0; z<dim[2]; z++)
			shape[k][x][y][z] = s[k][x][y][z];
	}
	
	final public int 	getLabelSamples() { return labelSamples; }
	
	final public int 		getIntensitySamples() { return intensitySamples; }
	
	final public float[] getTransform() { 
		if (transformMode==MULTIPLE) buildAverageTransformFromMultiTransform();
		return transform; 
	}
	
	final public float[][][][] exportDeformationField() { 
		if (transformMode==DEFORMABLE) return demons.exportTransformField();
		else return exportRigidDeformationField(); 
	}
	final public void 		setTransform(float[] trans) { transform = trans; }
	
	final public float[] 	getIntensity(int id) { 
		if (hasIntensity[id]) return intensity[id]; 
		else return exportLabels();
	}
	
	final public float[][] 	getIntensityVariancePriors(String[] modality, int nc) { 
		float[][]	prior = new float[classes][nc];
		for (int n=0;n<nc;n++) {
			if (modalityId(modality[n])>-1) {
				for (int k=0;k<classes;k++) prior[k][n] = intensityVariance[modalityId(modality[n])][k];
			} else {
				for (int k=0;k<classes;k++) prior[k][n] = 1.0f;
			}	
		}
		return prior;
	}
	
	final public float[][] 	getIntensityPriors(String[] modality, int nc) {
		float[][]	prior = new float[nc][classes];
		for (int n=0;n<nc;n++) {
			if ( (modalityId(modality[n])>-1) && (hasIntensity[modalityId(modality[n])]) ) {
				for (int k=0;k<classes;k++) prior[n][k] = intensity[modalityId(modality[n])][k];
			} else {
				for (int k=0;k<classes;k++) prior[n][k] = label[k];
			}	
		}
		return prior;
	}
	
	final public float[] 	getLesionPriors(String[] modality, int nc) {
		float[]	prior = new float[nc];
		if (hasLesions) {
			for (int n=0;n<nc;n++) {
				if ( (modalityId(modality[n])>-1) && (hasIntensity[modalityId(modality[n])]) ) {
					prior[n] = lesion[modalityId(modality[n])];
				} else {
					prior[n] = 0.0f;
				}
			}	
		} else {
			for (int n=0;n<nc;n++) prior[n] = 0.0f;
		}
		return prior;
	}
	
	final public float[] 	getBlackHolePriors(String[] modality, int nc) {
		float[]	prior = new float[nc];
		if (hasLesions) {
			for (int n=0;n<nc;n++) {
				if (modality[n].equals("T1_SPGR")) prior[n] = blackHole[T1_SPGR];
				else if (modality[n].equals("T2")) prior[n] = blackHole[T2];
				else if (modality[n].equals("FLAIR")) prior[n] = blackHole[FLAIR];
				else if (modality[n].equals("T1_MPRAGE")) prior[n] = blackHole[T1_MPRAGE];
				//else if (modality[n].equals("T1_RAW")) prior[n] = blackHole[T1_RAW];
				else if (modality[n].equals("PD")) prior[n] = blackHole[PD];
			}	
		} else {
			for (int n=0;n<nc;n++) prior[n] = 0.0f;
		}
		return prior;
	}
	
	final public float[][] getModalityWeights(String[] modality, int nc) {
		float[][]	w = new float[nc][OBJTYPES];
		for (int n=0;n<nc;n++) {
			if ( (modalityId(modality[n])>-1) && (hasIntensity[modalityId(modality[n])]) ) {
				for (int k=0;k<OBJTYPES;k++) w[n][k] = modweight[modalityId(modality[n])][k];
			} else {
				for (int k=0;k<OBJTYPES;k++) w[n][k] = 1.0f;
			}	
		}
		if (debug) {
			for (int n=0;n<nc;n++) 
				System.out.print("modweight("+n+") = "+w[n][0]+", "+w[n][1]+"\n");
		}
		
		return w;
	}
	
	final public int[][] getIntensityGroups(String[] modality, int nc) {
		int[][]	group = new int[nc][classes];
		int lb = 1;
		for (int n=0;n<nc;n++) {
			for (int k=0;k<classes;k++) {
				group[n][k]=0;
			}
			if ( (modalityId(modality[n])>-1) && (hasIntensity[modalityId(modality[n])]) ) {
				int m = modalityId(modality[n]);
				for (int k=0;k<classes;k++) {
					for (int l=k;l<classes;l++) {
						if (intensity[m][l]==intensity[m][k] && group[n][l]==0) {
							group[n][l] = k+1;
						}
					}
				}
			} else {
				for (int k=0;k<classes;k++) group[n][k] = k+1;
			}	
		}
		return group;
	}
	
	final public int[][] getSortedIntensityGroups(String[] modality, int nc) {
		int[][]	group = new int[nc][classes];
		for (int n=0;n<nc;n++) {
			int ng = 0;
			for (int k=0;k<classes;k++) {
				group[n][k]=0;
			}
			int m = modalityId(modality[n]);
			if ( (m>-1) && (hasIntensity[m]) ) {
				for (int k=0;k<classes;k++) {
					for (int l=k;l<classes;l++) {
						if (intensity[m][l]==intensity[m][k] && group[n][l]==0) {
							if (l==k) {
								ng++;
								group[n][l] = ng;
							} else {
								group[n][l] = group[n][k];
							}
						}
					}
				}
			} else {
				for (int k=0;k<classes;k++) group[n][k] = (k+1);
			}	
		}
		
		return group;
	}
	
	final public float[] getOptimizedFactor(String modality) {
		float[]	pw = new float[OPTIMIZED];
		if ( (modalityId(modality)>-1) && (hasIntensity[modalityId(modality)]) ) {
			for (int k=0;k<OPTIMIZED;k++) pw[k] = optimizedFactor[modalityId(modality)][k];
		} else {
			for (int k=0;k<OPTIMIZED;k++) pw[k] = 1.0f;
		}	
		if (debug) {
			System.out.print("optimizedfactor = "+pw[0]+", "+pw[1]+", "+pw[2]+"\n");
		}
		
		return pw;
	}
	
	final public void setOptimizedFactor(String modality,float[] pw) {
		if ( (modalityId(modality)>-1) && (hasIntensity[modalityId(modality)]) ) {
			for (int k=0;k<OPTIMIZED;k++) optimizedFactor[modalityId(modality)][k]=pw[k];
		} else {
			for (int k=0;k<OPTIMIZED;k++) optimizedFactor[modalityId(modality)][k] = 1.0f;
		}	
	}
	
	final public float[]	exportLabels() {
		float[] lb = new float[classes];
		for (int k=0;k<classes;k++) lb[k] = (float)label[k];
		return lb;
	}
	
	final public void importShapes(float[][][][] img, int nax_, int nay_, int naz_) {
		for (int k=0;k<classes;k++) {
			nax = nax_;
			nay = nay_;
			naz = naz_;
			shape[k] = new float[nax][nay][naz];
		}
		for (int k=0;k<classes;k++) for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) 
			shape[k][x][y][z] = img[k][x][y][z];
	}
	
	final public int[] getAtlasDim() {
		int[] dim = new int[3];
		dim[0] = nax;
		dim[1] = nay;
		dim[2] = naz;
		return dim;
	}
	
    final public int[] getImageDim() {
		int[] dim = new int[3];
		dim[0] = nix;
		dim[1] = niy;
		dim[2] = niz;
		return dim;
	}
	
    final public float[] getImageRes() {
		float[] res = new float[3];
		res[0] = rix;
		res[1] = riy;
		res[2] = riz;
		return res;
	}
	
	final public int[] getImageOrient() {
		int[] ori = new int[4];
		ori[0] = orient;
		ori[1] = orix;
		ori[2] = oriy;
		ori[3] = oriz;
		return ori;
	}
	final public boolean hasTopology() { return hasTopology; }
	final public boolean hasShape(int id) { return hasShape[id]; }
	
	final public boolean isDeformable() { return (transformMode==DEFORMABLE); }
	
	/** keeping track of transforms inside the atlas */
	final public void buildAverageTransformFromMultiTransform() {
		
		if (transformMode==MULTIPLE) {
			transform = new float[Nd];
			for (int n=0;n<Nd;n++) {
				transform[n] = 0.0f;
			
				// simplistic average
				for (int k=1;k<classes;k++) 
					transform[n] += multitransform[k][n]/(float)(classes-1);
			}
			RotationMatrix R = new RotationMatrix();
			// set up rotation parameters
			R.setParameters(transform[0],transform[1],transform[2]);
			rotation = R.getMatrix();
			
		}
	}
		
    /** 
	 *  set image-related information for segmentation
	 */
	final public void setImageInfo(int nix_, int niy_, int niz_, float rix_, float riy_, float riz_, int orient_, int orix_, int oriy_, int oriz_) {
		nix = nix_; niy = niy_; niz = niz_;
		rix = rix_; riy = riy_; riz = riz_;
		orient = orient_;
		orix = orix_; oriy = oriy_; oriz = oriz_;
		
		x0i = nix/2.0f;
		y0i = niy/2.0f;
		z0i = niz/2.0f;
		
		if (debug) {
			System.out.print("dimensions: "+nix+", "+niy+", "+niz+"\n");
			System.out.print("resolutions: "+rix+", "+riy+", "+riz+"\n");
			System.out.print("orientation: "+orient+" | "+orix+", "+oriy+", "+oriz+"\n");
			System.out.print("Atlas\n");
			System.out.print("dimensions: "+nax+", "+nay+", "+naz+"\n");
			System.out.print("resolutions: "+rax+", "+ray+", "+raz+"\n");
		}
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public byte[][][] generateClassification() {
		float dist,max,count;
		int b=0;
		byte[][][] img = new byte[nix][niy][niz];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			// compute each class probability : attribute the highest
			max = 0; b = -1;
			for (int k=0;k<classes;k++) {
				if (shape[k][x][y][z]>max) {
					b = k;
					max = shape[k][x][y][z];
				}
			}
			if (b>-1) img[x][y][z] = label[b];
			else img[x][y][z] = 0;
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public byte[][][] generateTransformedClassification() {
		float dist,max,count,val;
		int b=0;
		byte[][][] img = new byte[nix][niy][niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			// compute each class probability : attribute the highest
			max = 0; b = -1;
			for (int k=0;k<classes;k++) {
				if (transformMode==DEFORMABLE)
					XP = demons.getCurrentMapping(x,y,z);
				else if (transformMode==MULTIPLE)
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				else
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation, 1.0f);
				
				val = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);	
				if (val>max) {
					b = k;
					max = val;
				}
			}
			if (b>-1) img[x][y][z] = label[b];
			else img[x][y][z] = 0;
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[][][][] generateTransformedShapes() {
		/*if (transformMode==DEFORMABLE) return demons.exportTransformedAtlas();
			
		float[][][][] img = new float[classes][nix][niy][niz];
		float[] XP=new float[3];
		
		for (int k=0;k<classes;k++) {
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				if (transformMode==MULTIPLE)
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				else
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation, 1.0f);
				
				img[k][x][y][z] = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);	
			}
		}
		return img;*/
    	float[][][][] img = new float[classes][nix][niy][niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			float[] point_atlas = getTransformedShape(x,y,z);
			for (int k=0; k<classes; k++)
				img[k][x][y][z] = point_atlas[k];
		}
		return img;
	}
    
    final public float[][][][] generateTransposeTransformedShapes() {
		float[][][][] img = new float[nix][niy][niz][classes];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			float[] point_atlas = getTransformedShape(x,y,z);
			for (int k=0; k<classes; k++)
				img[x][y][z][k] = point_atlas[k];
		}
		return img;
	}
    
    
    final public float[][][][] generateNonZeroTransformedShapes() {
		if (transformMode==DEFORMABLE) return demons.exportTransformedAtlas();
			
		float[][][][] img = new float[classes][nix][niy][niz];
		float[] XP=new float[3];
		
		for (int k=0;k<classes;k++) {
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				if (transformMode==MULTIPLE)
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				else
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation, 1.0f);
				
				img[k][x][y][z] = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)+0.001f;
			}
		}
		return img;
	}
	
    final public float[] generateTransformedShapes1D() {
		//if (transformMode==DEFORMABLE) return demons.exportTransformedAtlas();
			
		float[] img = new float[classes*nix*niy*niz];
		float[] XP=new float[3];
		
		for (int k=0;k<classes;k++) {
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				if (transformMode==MULTIPLE)
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				else
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation, 1.0f);
				
				img[k*nix*niy*niz+x*niy*niz+y*niz+z] = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);	
			}
		}
		return img;
	}
    /** 
	 *  generate atlas image from information
	 */
    final public float[][][] generateIntensityImage(int modality) {
		float dist,max,count;
		float[]	energy = new float[classes];
		int b=0;
		float[][][] img = new float[nix][niy][niz];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			// compute each class probability : attribute the highest
			max = 0; b = -1;
			for (int k=0;k<classes;k++) {
				if (shape[k][x][y][z]>max) {
					b = k;
					max = shape[k][x][y][z];
				}
			}
			if (b>-1) {
				if (hasIntensity[modality]) img[x][y][z] = intensity[modality][b];
				else img[x][y][z] = (float)label[b];
			} else img[x][y][z] = 0;
		}
		return img;
	}
	
	/** display the atlas data */
	final public String displayIntensity() {
		String output = "Intensity \n";
		
		output += "T1_SPGR ("+hasIntensity[0]+") : ";
		for (int k=0;k<classes;k++) output += intensity[0][k]+" ";
		output += "\n";
		output += "T2 ("+hasIntensity[1]+") : ";
		for (int k=0;k<classes;k++) output += intensity[1][k]+" ";
		output += "\n";
		output += "FLAIR ("+hasIntensity[2]+") : ";
		for (int k=0;k<classes;k++) output += intensity[2][k]+" ";
		output += "\n";
		output += "T1_MPRAGE ("+hasIntensity[3]+") : ";
		for (int k=0;k<classes;k++) output += intensity[3][k]+" ";
		output += "\n";
		output += "T1_RAW ("+hasIntensity[4]+") : ";
		for (int k=0;k<classes;k++) output += intensity[4][k]+" ";
		output += "\n";
		output += "PD ("+hasIntensity[5]+") : ";
		for (int k=0;k<classes;k++) output += intensity[5][k]+" ";
		output += "\n";
		output += "PDFSE ("+hasIntensity[6]+") : ";
		for (int k=0;k<classes;k++) output += intensity[6][k]+" ";
		output += "\n";
		
		return output;	
	}
	/** display the atlas data */
	final public String displayIntensityVariance() {
		String output = "Intensity Variance \n";
		
		output += "T1_SPGR ("+hasIntensityVariance[0]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[0][k]+" ";
		output += "\n";
		output += "T2 ("+hasIntensityVariance[1]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[1][k]+" ";
		output += "\n";
		output += "FLAIR ("+hasIntensityVariance[2]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[2][k]+" ";
		output += "\n";
		output += "T1_MPRAGE ("+hasIntensityVariance[3]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[3][k]+" ";
		output += "\n";
		output += "T1_RAW ("+hasIntensityVariance[4]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[4][k]+" ";
		output += "\n";
		output += "PD ("+hasIntensityVariance[5]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[5][k]+" ";
		output += "\n";
		output += "PDFSE ("+hasIntensityVariance[6]+") : ";
		for (int k=0;k<classes;k++) output += intensityVariance[6][k]+" ";
		output += "\n";
		
		return output;	
	}
	/** display the atlas data */
	final public String displayModalityWeights() {
		String output = "modality Weights \n";
		
		output += "T1 : ";
		for (int k=0;k<OBJTYPES;k++) output += modweight[0][k]+" ";
		output += "\n";
		output += "T2 : ";
		for (int k=0;k<OBJTYPES;k++) output += modweight[1][k]+" ";
		output += "\n";
		output += "FLAIR : ";
		for (int k=0;k<OBJTYPES;k++) output += modweight[2][k]+" ";
		output += "\n";
		output += "MPRAGE : ";
		for (int k=0;k<OBJTYPES;k++) output += modweight[3][k]+" ";
		output += "\n";
		output += "MPRAGE_RAW : ";
		for (int k=0;k<OBJTYPES;k++) output += modweight[4][k]+" ";
		output += "\n";
		
		return output;	
	}
	/** display the atlas data */
	final public String displayOptimizedFactors() {
		String output = "optimized factors \n";
		
		output += "T1 : ";
		for (int k=0;k<OPTIMIZED;k++) output += optimizedFactor[0][k]+" ";
		output += "\n";
		output += "T2 : ";
		for (int k=0;k<OPTIMIZED;k++) output += optimizedFactor[1][k]+" ";
		output += "\n";
		output += "FLAIR : ";
		for (int k=0;k<OPTIMIZED;k++) output += optimizedFactor[2][k]+" ";
		output += "\n";
		output += "MPRAGE : ";
		for (int k=0;k<OPTIMIZED;k++) output += optimizedFactor[3][k]+" ";
		output += "\n";
		output += "MPRAGE_RAW : ";
		for (int k=0;k<OPTIMIZED;k++) output += optimizedFactor[4][k]+" ";
		output += "\n";
		
		return output;	
	}
	/** display the atlas data */
	final public String displayLesionModel() {
		String output = "Lesion model ("+hasLesions+") : \n";
		
		output += "T1 "+lesion[0];
		output += ", T2 "+lesion[1];
		output += ", FLAIR "+lesion[2];
		output += ", MPRAGE "+lesion[3];
		output += ", MPRAGE_RAW "+lesion[4];
		output += ", PD "+lesion[5];
		output += ", PDFSE "+lesion[6];
		output += "\n";
		
		return output;	
	}
	/** display the atlas data */
	final public String displayNames() {
		String output = "Structures \n";
		
		for (int k=0;k<classes;k++) {
			output += name[k]+" ("+topology[k]+")	"+label[k]+"\n";
		}
		
		return output;	
	}
	/** display the atlas data */
	final public String displayRegisteredShapes() {
		String output = "Registered Shapes (0/1: true/false) \n";
		
		for (int k=0;k<classes;k++) output += k+"	";
		output += "\n";
		for (int k=0;k<classes;k++) {
			if (registeredShape[k]) output += "1	";
			else output += "0	";
		}
		output += "\n";
		
		return output;	
	}

	/**
	 *	read template image (the image must be in bytes)
	 */
	private byte[][][] loadTemplateImage(String filename, int Nx, int Ny, int Nz) {
		// read the raw data
		byte[] buffer = null;
		try {
          File f = new File( filename );
		  //System.out.println("exists ? "+f.exists());
          //System.out.println("can read ? "+f.canRead());
          FileInputStream fis = new FileInputStream( f );
            
		   buffer = new byte[Nx*Ny*Nz];
		   fis.read(buffer);
           fis.close();
		} catch (IOException io) {
           System.out.println("i/o pb: "+io.getMessage());
		}
        // convert to the image format
		byte [][][] img  = new byte[Nx][Ny][Nz];
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			img[x][y][z] = buffer[x + Nx*y + Nx*Ny*z];
		}
		buffer = null;
		
		return img;
	}
	/**
	 *	read shape image (the image must be in float, little endian)
	 */
	private final float[][][] loadShapeImage(String filename, int Nx, int Ny, int Nz) {
		// read the raw data
		byte[] buffer = null;
		try {
           File f = new File( filename );
           FileInputStream fis = new FileInputStream( f );
            
		   buffer = new byte[4*Nx*Ny*Nz];
		   fis.read(buffer);
           fis.close();
		} catch (IOException io) {
           System.out.println("i/o pb: "+io.getMessage());
		}
        // convert to the image format
		float [][][] img  = new float[Nx][Ny][Nz];
		
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			int b1 = buffer[4*(x+Nx*y+Nx*Ny*z)+0] & 0xff;
			int b2 = buffer[4*(x+Nx*y+Nx*Ny*z)+1] & 0xff;
			int b3 = buffer[4*(x+Nx*y+Nx*Ny*z)+2] & 0xff;
			int b4 = buffer[4*(x+Nx*y+Nx*Ny*z)+3] & 0xff;
			// big endian
			//int tmpInt = ((b1 << 24) | (b2 << 16) | (b3 << 8) | b4);
			// little endian
			int tmpInt = ((b4 << 24) | (b3 << 16) | (b2 << 8) | b1);

			img[x][y][z] = Float.intBitsToFloat(tmpInt);
		}
        buffer = null;
		
		return img;
	}
	
	/** 
	 *	load the atlas data from a file. 
	 *  All associated images are loaded at this time
	 */
	final public void loadAtlas(String filename) {
		if (verbose) System.out.println("loading atlas file: "+filename);
		try {
            File f = new File(filename);
			String dir = f.getParent();
            FileReader fr = new FileReader(f);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
			StringTokenizer st;
			String imageFile;
            // Exact corresponding template
            if (!line.equals("Structure Atlas File (edit at your own risks)")) {
                System.out.println("not a proper Structure Atlas file");
                br.close();
                fr.close();
                return;
            }
			line = br.readLine();
			while (line!=null) {
				if (line.startsWith("Structures")) {
					//System.out.println(line);
					// Structures:	classes	label	topology
					st = new StringTokenizer(line, "	");
					st.nextToken();
					classes = BasicInfo.getInt(st);
					name = new String[classes];
					label = new byte[classes];
					topology = new String[classes];
					for (int n=0;n<classes;n++) {
						// Name:label:topology
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						name[n] = st.nextToken();
						label[n] = (byte)BasicInfo.getInt(st);
						topology[n] = st.nextToken();
					}
					// allocate other quantities
					hasTopology = false;
					hasShape = new boolean[classes];
					for (int n=0;n<classes;n++) hasShape[n] = false;
					shape = new float[classes][][][];
					minx = new int[classes];
					miny = new int[classes];
					minz = new int[classes];
					maxx = new int[classes];
					maxy = new int[classes];
					maxz = new int[classes];
					hasIntensity = new boolean[INTENSITY];
					for (int i=0;i<INTENSITY;i++) hasIntensity[i] = false;
					intensity = new float[INTENSITY][classes];
					lesion = new float[INTENSITY];
					blackHole = new float[INTENSITY];
					hasIntensityVariance = new boolean[INTENSITY];
					for (int i=0;i<INTENSITY;i++) hasIntensityVariance[i] = false;
					intensityVariance = new float[INTENSITY][classes];
					for (int i=0;i<INTENSITY;i++) for (int j=0;j<classes;j++) intensityVariance[i][j] = 1.0f;
					modweight = new float[INTENSITY][OBJTYPES];
					for (int i=0;i<INTENSITY;i++) for (int j=0;j<OBJTYPES;j++) modweight[i][j] = 1.0f;
					optimizedFactor = new float[INTENSITY][OPTIMIZED];
					for (int i=0;i<INTENSITY;i++) for (int j=0;j<OPTIMIZED;j++) optimizedFactor[i][j] = 1.0f;
					lesion = new float[INTENSITY];
					shapeFile = new String[classes];
					for (int n=0;n<classes;n++) shapeFile[n] = null;
					registeredShape = new boolean[classes];
					for (int n=0;n<classes;n++) registeredShape[n] = true;
					registeredShape[0] = false;
					templateFile = null;
					if (debug) System.out.println(displayNames());
				} else
				if (line.startsWith("Topology Atlas")) {
					//System.out.println(line);
					// File: name
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					imageFile = dir+File.separator+st.nextToken();
					if (debug) System.out.print("file: "+imageFile+"\n");
					// Dimensions: nax nay naz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					nax = BasicInfo.getInt(st);
					nay = BasicInfo.getInt(st);
					naz = BasicInfo.getInt(st);
					x0a = nax/2.0f;
					y0a = nay/2.0f;
					z0a = naz/2.0f;
					if (debug) System.out.print("dims: "+nax+" "+nay+" "+naz+"\n");
					// Resolutions: rax ray raz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					rax = BasicInfo.getFloat(st);
					ray = BasicInfo.getFloat(st);
					raz = BasicInfo.getFloat(st);
					if (debug) System.out.print("res: "+rax+"x"+ray+"x"+raz+"\n");
					template = loadTemplateImage(imageFile, nax, nay, naz);
					hasTopology = true;
					templateFile = imageFile;
				} else
				if (line.startsWith("Shape Atlas")) {
					//if (debug) System.out.println(line);
					// Shape:	labelSamples
					st = new StringTokenizer(line, "	");
					st.nextToken();
					labelSamples = BasicInfo.getInt(st);
					// Dimensions: nix niy niz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					nix = BasicInfo.getInt(st);
					niy = BasicInfo.getInt(st);
					niz = BasicInfo.getInt(st);
					x0i = nix/2.0f;
					y0i = niy/2.0f;
					z0i = niz/2.0f;
					if (debug) System.out.print("image dim: "+nix+"x"+niy+"x"+niz+"\n");
					// Resolutions: rix riy riz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					rix = BasicInfo.getFloat(st);
					riy = BasicInfo.getFloat(st);
					riz = BasicInfo.getFloat(st);
					if (debug) System.out.print("image res: "+rix+"x"+riy+"x"+riz+"\n");
					// Orientations: orient orix oriy oriz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					orient = BasicInfo.getInt(st);
					orix = BasicInfo.getInt(st);
					oriy = BasicInfo.getInt(st);
					oriz = BasicInfo.getInt(st);
					if (debug) System.out.print("image orient: "+orient+"|"+orix+"x"+oriy+"x"+oriz+"\n");
					line = br.readLine();
					while (line.startsWith("Structure:")) {
						// find structure id
						st = new StringTokenizer(line, "	");
						st.nextToken();
						String title = st.nextToken();
						int id=-1;
						for (int n=0;n<classes;n++) {
							if (title.equals(name[n])) { id = n; break; }
						}
						if (debug) System.out.print("Shape: "+name[id]+"\n");
						if (id==-1) {
							line = br.readLine();
						} else {
							// File: name
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							st.nextToken();
							imageFile = dir+File.separator+st.nextToken();
							// min, max : initial values
							minx[id] = 0; miny[id] = 0; minz[id] = 0;
							maxx[id] = nax; maxy[id] = nay; maxz[id] = naz;
			
							shape[id] = loadShapeImage(imageFile, nax, nay, naz);
							hasShape[id] = true;
							shapeFile[id] = imageFile;
						}
						line = br.readLine();
					}
				} else
				if (line.startsWith("Intensity Atlas")) {
					//if (debug) System.out.println(line);
					// Intensity:	intensitySamples
					st = new StringTokenizer(line, "	");
					st.nextToken();
					intensitySamples = BasicInfo.getInt(st);
					// Type value value value...
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					String type = st.nextToken();
					//if (debug) System.out.println(type);
					int id = modalityId(type);
					
					while (id !=-1 ) {
						for (int n=0;n<classes;n++) {
							intensity[id][n] = BasicInfo.getFloat(st);
						}
						hasIntensity[id] = true;
						// search for next intensity profile
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						type = st.nextToken();
						if (debug) System.out.println(type);
						id = modalityId(type);
					}
					if (debug) System.out.println(displayIntensity());
				} else
				if (line.startsWith("Intensity Variance")) {
					// Type value value value...
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					String type = st.nextToken();
					//if (debug) System.out.println(type);
					int id = modalityId(type);
					
					while (id !=-1 ) {
						for (int n=0;n<classes;n++) {
							intensityVariance[id][n] = BasicInfo.getFloat(st);
						}
						hasIntensityVariance[id] = true;
						// search for next intensity profile
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						type = st.nextToken();
						if (debug) System.out.println(type);
						id = modalityId(type);
					}
					if (debug) System.out.println(displayIntensityVariance());
				} else
				if (line.startsWith("Lesions")) {
					// Lesions: 
					// T1 val T2 val FLAIR val MPRAGE val
					line = br.readLine();
					if (debug) System.out.println(line);
					st = new StringTokenizer(line, "	");
					while (st.hasMoreTokens()) {
						String type = st.nextToken();
						if (debug) System.out.print(type+" ");
						int id = modalityId(type);
						if (id>-1) {
							lesion[id] = BasicInfo.getFloat(st);
							if (debug) System.out.println(lesion[id]);
						}
					}
					if (debug) System.out.println("\n");
					hasLesions = true;
					if (debug) System.out.println(displayLesionModel());
				} else
				if (line.startsWith("Black Holes")){
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					while (st.hasMoreTokens()) {
						String type = st.nextToken();
						if (debug) System.out.print(type+" ");
						int id = -1;
						/*if (type.equals("T1_SPGR")) id = T1_SPGR;
						else if (type.equals("T2")) id = T2;
						else if (type.equals("FLAIR")) id = FLAIR;
						else if (type.equals("T1_MPRAGE")) id = T1_MPRAGE;
						else if (type.equals("T1_RAW")) id = T1_RAW;
						else if (type.equals("PD")) id = PD;*/
						id = modalityId(type);
						if (id>-1) {
							blackHole[id] = BasicInfo.getFloat(st);
							if (debug) System.out.println(blackHole[id]);
						}
					}
					if (debug) System.out.println("\n");
					hasLesions = true;
					if (debug) System.out.println(displayLesionModel());
				} else
				if (line.startsWith("Slope")) {
					// Slope: val 
					st = new StringTokenizer(line, "	");
					String text = st.nextToken();
					if (debug) System.out.println(text);
					shapeSlope = BasicInfo.getFloat(st);
					if (debug) System.out.println("shape slope = "+shapeSlope);
				} else
				if (line.startsWith("Modality weights")) {
					//if (debug) System.out.println(line);
					// Intensity:	intensitySamples
					st = new StringTokenizer(line, "	");
					st.nextToken();
					// Type value value
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					String type = st.nextToken();
					if (debug) System.out.println(type);
					int id = modalityId(type);
					
					while (id !=-1 ) {
						for (int n=0;n<OBJTYPES;n++) {
							modweight[id][n] = BasicInfo.getFloat(st);
						}
						// search for next intensity profile
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						type = st.nextToken();
						if (debug) System.out.println(type);
						id = modalityId(type);
					}
					if (debug) System.out.println(displayModalityWeights());
				} else
				if (line.startsWith("Registered Shapes")) {
					if (debug) System.out.println(line);
					// Type value value
					line = br.readLine();
					if (debug) System.out.println(line);
					st = new StringTokenizer(line, "	");
					for (int n=0;n<classes;n++) {
						registeredShape[n] = (BasicInfo.getInt(st)==1);
					}
					if (debug) System.out.println(displayRegisteredShapes());
				} else
				if (line.startsWith("Optimized factors")) {
					//if (debug) System.out.println(line);
					// Intensity:	intensitySamples
					st = new StringTokenizer(line, "	");
					st.nextToken();
					// Type value value
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					String type = st.nextToken();
					if (debug) System.out.println(type);
					int id = modalityId(type);
					
					while (id !=-1 ) {
						for (int n=0;n<OPTIMIZED;n++) {
							optimizedFactor[id][n] = BasicInfo.getFloat(st);
						}
						// search for next intensity profile
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						type = st.nextToken();
						if (debug) System.out.println(type);
						id = modalityId(type);
					}
					if (debug) System.out.println(displayOptimizedFactors());
				}
				line = br.readLine();
				if (debug) System.out.println(line);
			}		
			br.close();
            fr.close();
			atlasFile = filename;
        }
        catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
        catch (IOException e) {
            System.out.println(e.getMessage());
        } 
		catch (OutOfMemoryError e){
			System.out.println(e.getMessage());
		}
		catch (Exception e) {
			System.out.println(e.getMessage());
        }

		if (debug) System.out.println("initialisation\n");
	}

	public static final int modalityId(String type) {
		int id = -1;
		if (type.equals("T1") || type.equals("T1_SPGR")) id = T1_SPGR;
		else if (type.equals("T2")) id = T2;
		else if (type.equals("FLAIR")) id = FLAIR;
		else if (type.equals("MPRAGE") || type.equals("T1_MPRAGE")) id = T1_MPRAGE;
		else if (type.equals("MPRAGE_RAW") || type.equals("T1_RAW")) id = T1_RAW;
		else if (type.equals("PD")) id = PD;
		else if (type.equals("PDFSE")) id = PDFSE;
		else if (type.equals("DIR")) id = DIR;
		
		return id;
	}					
	
	public final float[] getTransformedShape(int x, int y, int z) {
		if (precompute) return getFastTransformedShape(x,y,z);
		else return getRegularTransformedShape(x,y,z);
	}
	
	/** transformations: how to get a transformed value 
	 *	simple hypotheses: start from the same system (same origin, same resolution)
	 */
	private final float[] getRegularTransformedShape(int x, int y, int z) {
		float[] val = new float[classes];
		if (transformMode==NONE) {
			for (int k=0;k<classes;k++) val[k] = shape[k][x][y][z];
			return val;
		}
		
		float[] XP=new float[3];
		boolean noCoordinates = true;
		for (int k=0;k<classes;k++) {
			if ( (x<minx[k]) || (x>=maxx[k]) || (y<miny[k]) || (y>=maxy[k]) || (z<minz[k]) || (z>=maxz[k]) ) {
				if (k==0) val[k] = 1.0f;
				else val[k] = 0.0f;
			} else {
				if (noCoordinates && transformMode==DEFORMABLE) {
					XP = demons.getCurrentMapping(x,y,z);
					noCoordinates = false;
				} else if (transformMode==MULTIPLE) {
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				} else if (noCoordinates) {
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation,1.0f);
					noCoordinates = false;
				}
				//val[k] = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);
				// post-processing ??
				val[k] = Numerics.bounded(shapeSlope*(ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
			}
		}		
		return val;
	}
	
	/** transformations: how to get a transformed value 
	*	simple hypotheses: start from the same system (same origin, same resolution)
	*/
	private final float[] getFastTransformedShape(int x, int y, int z) {
		float[] val = new float[classes];
		if (transformMode==NONE) {
			for (int k=0;k<classes;k++) val[k] = shape[k][x][y][z];
			return val;
		}
		
		float[] XP=new float[3];
		boolean noCoordinates = true;
		for (int k=0;k<classes;k++) {
			if ( (x<minx[k]) || (x>=maxx[k]) || (y<miny[k]) || (y>=maxy[k]) || (z<minz[k]) || (z>=maxz[k]) ) {
				if (k==0) val[k] = 1.0f;
				else val[k] = 0.0f;
			} else {
				if (noCoordinates) {
					XP = fastImageToShapeCoordinates(k,x,y,z);
					noCoordinates = false;
				}
				//val[k] = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);
				val[k] = Numerics.bounded(shapeSlope*(ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
			}
		}		
		return val;
	}
	
	public final float getTransformedShape(int x, int y, int z, int k) {
		if (precompute) return getFastTransformedShape(x,y,z,k);
		else return getRegularTransformedShape(x,y,z,k);
	}
	
	/** transformations: how to get a transformed value 
	 *	simple hypotheses: start from the same system (same origin, same resolution)
	 */
	private final float getRegularTransformedShape(int x, int y, int z, int k) {
		float val = 0.0f;
		if (transformMode==NONE) {
			return val;
		}
		
		float[] XP=new float[3];
		if ( (x<minx[k]) || (x>=maxx[k]) || (y<miny[k]) || (y>=maxy[k]) || (z<minz[k]) || (z>=maxz[k]) ) {
			if (k==0) val = 1.0f;
			else val = 0.0f;
		} else {
			if (transformMode==DEFORMABLE) {
				XP = demons.getCurrentMapping(x,y,z);
			} else if (transformMode==MULTIPLE) {
				XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
			} else {
				XP = transformModel.imageToTemplate(x,y,z,transform,rotation,1.0f);
			}
			//val = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);
			val = Numerics.bounded(shapeSlope*(ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
		}
		return val;
	}
	
	/** transformations: how to get a transformed value 
	*	simple hypotheses: start from the same system (same origin, same resolution)
	*/
	private final float getFastTransformedShape(int x, int y, int z, int k) {
		float val = 0.0f;
		if (transformMode==NONE) {
			return val;
		}
		
		float[] XP=new float[3];
		boolean noCoordinates = true;
		if ( (x<minx[k]) || (x>=maxx[k]) || (y<miny[k]) || (y>=maxy[k]) || (z<minz[k]) || (z>=maxz[k]) ) {
			if (k==0) val = 1.0f;
			else val = 0.0f;
		} else {
			if (noCoordinates) {
				XP = fastImageToShapeCoordinates(k,x,y,z);
				noCoordinates = false;
			}
			//val = ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz);
			val = Numerics.bounded(shapeSlope*(ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
		}		
		return val;
	}
	
	/** transformations: compute the shape bounding box in image with current transform
	 */
    public final void computeTransformedShapeBoundingBox() {
		float[] XP=new float[3];
		
		for (int k=0;k<classes;k++) {
			if (precompute) transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],transform,rotation,1.0f);
			minx[k] = nix; miny[k] = niy; minz[k] = niz; 
			maxx[k] = 0; maxy[k] = 0; maxz[k] = 0;
		
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				if (precompute) {
					XP = fastImageToShapeCoordinates(k,x,y,z);
				} else if (transformMode==DEFORMABLE) {
					XP = demons.getCurrentMapping(x,y,z);
				} else if (transformMode==MULTIPLE) {
					XP = transformModel.imageToTemplate(x,y,z, multitransform[k], multirotation[k], 1.0f);
				} else {
					XP = transformModel.imageToTemplate(x,y,z,transform,rotation,1.0f);
				}
				
				if (ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)>0) {
					if (x<minx[k]) minx[k] = x;
					if (y<miny[k]) miny[k] = y;
					if (z<minz[k]) minz[k] = z;
					if (x>maxx[k]) maxx[k] = x;
					if (y>maxy[k]) maxy[k] = y;
					if (z>maxz[k]) maxz[k] = z;
				}
			}
		}
		return;
	}
	
	/** transformations: compute the shape bounding box in image with current transform
	 */
	public final void computeSingleTransformedShapeBoundingBox(int k) {
		float[] XP=new float[3];
		
		if (precompute) 
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],multitransform[k],multirotation[k],1.0f);
		minx[k] = nix; miny[k] = niy; minz[k] = niz; 
		maxx[k] = 0; maxy[k] = 0; maxz[k] = 0;
	
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			if (precompute) XP = fastImageToShapeCoordinates(k,x,y,z);
			else XP = transformModel.imageToTemplate(x,y,z,multitransform[k],multirotation[k],1.0f);
			
			if (ImageFunctions.linearClosestInterpolation(shape[k],XP[0],XP[1],XP[2],nax,nay,naz)>0) {
				if (x<minx[k]) minx[k] = x;
				if (y<miny[k]) miny[k] = y;
				if (z<minz[k]) minz[k] = z;
				if (x>maxx[k]) maxx[k] = x;
				if (y>maxy[k]) maxy[k] = y;
				if (z>maxz[k]) maxz[k] = z;
			}
		}
		return;
	}
	
	/** transformations: re-compute the template using the transform
	 */
	public final void computeTransformedTemplate() {
		float[] XP=new float[3];
		byte[][][] tmp = new byte[nix][niy][niz];
		
		if (transformMode==MULTIPLE) 
			buildAverageTransformFromMultiTransform();
			
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			XP = transformModel.imageToTemplate(x,y,z,transform,rotation,1.0f);
			tmp[x][y][z] = ImageFunctions.nearestNeighborInterpolation(template,label[0],XP[0],XP[1],XP[2],nax,nay,naz);
			// debug
			boolean wrong=true;
			for (int k=0;k<classes;k++) {
				if (tmp[x][y][z]==label[k]) { wrong=false; }
			}
			if (wrong) {
				System.out.println(" bad interpolation!! : ("+x+", "+y+", "+z+") <- ("+XP[0]+", "+XP[1]+", "+XP[2]+") : "+tmp[x][y][z]);
			}
		}
		template = tmp;
		
		return;
	}
	
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeRegistrationCoefficients(float[] hessian, float[] gradient, float[] trans) {
        float vec,mat,dp, res,val,num=0,den,weight;
        float xP=0,yP=0,zP=0;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null,dRa = null,dRb = null,dRc = null;
        float[]         dprior = new float[Nd];
		float[] Xi;
        float[][] dXi;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		int				Npt;
		int				n,count;
		float			limit;
		float[]			mem = new float[classes];
		float 			cost, norm;
		float[]			ngrad = new float[Nd], nhess = new float[Nd];
		float			normalize = nix*niy*niz;
		
        // init the coefficients
		cost = 0.0f;
		norm = 0.0f;
		for (int i=0;i<Nd;i++) {
			hessian[i] = 0.0f;
			gradient[i] =  0.0f;
			ngrad[i] = 0.0f;
			nhess[i] = 0.0f;
		}
		
		// set up rotation parameters
        if (transformModel.useRotation()) {
		   R = transformModel.computeRotationMatrix(trans);
		   rot = R.getMatrix();
		   dRa = R.derivatives(1.0f, 0.0f, 0.0f);
		   dRb = R.derivatives(0.0f, 1.0f, 0.0f);
		   dRc = R.derivatives(0.0f, 0.0f, 1.0f);
		}
		
		if (precompute) for (int k=1;k<classes;k++) if (registeredShape[k]) 
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,1.0f);
			
		// main loop
		for (int x=offset;x<nix;x+=subsample) for (int y=offset;y<niy;y+=subsample) for (int z=offset;z<niz;z+=subsample) {
            // factor : classes
            vec = 0.0f; mat = 0.0f;
			for (int k=1;k<classes;k++) if (registeredShape[k]) {
				// compute the local position
				if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
				else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,1.0f);
				
				// compute interpolated values
				priorT = ImageFunctions.linearClosestInterpolation(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
				
				// check if the region is zero: no calculation needed then
				if (priorT>0) {
					// data term : function of the memberships alone
					//System.out.print(".");
					
					// derivatives
					dPx = ImageFunctions.linearInterpolationXderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
					dPy = ImageFunctions.linearInterpolationYderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
					dPz = ImageFunctions.linearInterpolationZderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
					
					// coordinate derivatives
					dXi = transformModel.imageToTemplateDerivatives(x,y,z,trans,rot,dRa,dRb,dRc,1.0f);
							
					// assemble everiything
					for (int i=0;i<Nd;i++) {
						dprior[i] = dPx*dXi[0][i] + dPy*dXi[1][i] + dPz*dXi[2][i];
					}
					cost += registrationCost(mems,priorT,x,y,z,k);
					norm += registrationNorm(mems,priorT,x,y,z,k);
					for (int i=0;i<Nd;i++) {
						gradient[i] += registrationCostGradient(mems,priorT,dprior,x,y,z,k,i);
						hessian[i] += registrationCostHessian(mems,priorT,dprior,x,y,z,k,i);
						ngrad[i] += registrationNormGradient(mems,priorT,dprior,x,y,z,k,i);
						nhess[i] += registrationNormHessian(mems,priorT,dprior,x,y,z,k,i);
					}
				}
			}
        }
		// assemble everything, derivatives first
		for (int i=0;i<Nd;i++) {
			hessian[i] = registrationMeasureHessian(hessian[i], nhess[i], norm, normalize);
			gradient[i] = registrationMeasureGradient(hessian[i], nhess[i], norm, normalize);
		}
		
		return registrationMeasure(cost, norm, normalize);
    } // computeRegistrationCoefficients
        
	/** transformations: how to update the transform using memberships */
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeRegistrationEnergy(float[] trans) {
        float weight;
        float[] Xi;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		float			cost, norm;
		float			normalize = nix*niy*niz;
		
        // set up rotation parameters
        if (transformModel.useRotation()) {
		   rot = transformModel.computeRotation(trans);
		}

		if (precompute) for (int k=1;k<classes;k++) if (registeredShape[k])
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,1.0f);
		
		// main loop
		cost = 0.0f;		
		norm = 0.0f;		
        for (int x=offset;x<nix;x+=subsample) for (int y=offset;y<niy;y+=subsample) for (int z=offset;z<niz;z+=subsample) {
            // factor : classes
            for (int k=1;k<classes;k++) if (registeredShape[k]) {
				// compute the local position
				if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
				else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,1.0f);
			
				// compute interpolated values
				priorT = ImageFunctions.linearClosestInterpolation(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
				
				// check if the region is zero: no calculation needed then
				if (priorT>0) {
					// data term : function of the memberships
					cost += registrationCost(mems,priorT,x,y,z,k);
					norm += registrationNorm(mems,priorT,x,y,z,k);
				}
			}
        }
 
        return registrationMeasure(cost, norm, normalize);
    } // computeRegistrationEnergy
    
	/** transformations: how to update the transform using memberships */
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeSingleRegistrationCoefficients(float[] hessian, float[] gradient, float[] trans, int k) {
        float vec,mat,dp, res,val,num=0,den,weight;
        float xP=0,yP=0,zP=0;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null,dRa = null,dRb = null,dRc = null;
        float[]         dprior = new float[Nd];
		float[] Xi;
        float[][] dXi;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		int				Npt;
		int				n,count;
		float			limit;
		float[]			mem = new float[classes];
		float			cost, norm;
		float[]			ngrad = new float[Nd], nhess = new float[Nd];
		float			normalize = nix*niy*niz;
		
        // init the coefficients
		cost = 0.0f;
		norm = 0.0f;
		for (int i=0;i<Nd;i++) {
			gradient[i] = 0.0f;
			hessian[i] = 0.0f;
			ngrad[i] = 0.0f;
			nhess[i] = 0.0f;
		}
		
        // set up rotation parameters
        if (transformModel.useRotation()) {
		   R = transformModel.computeRotationMatrix(trans);
		   rot = R.getMatrix();
		   dRa = R.derivatives(1.0f, 0.0f, 0.0f);
		   dRb = R.derivatives(0.0f, 1.0f, 0.0f);
		   dRc = R.derivatives(0.0f, 0.0f, 1.0f);
		}
		if (precompute)
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,1.0f);
		
		// main loop
		for (int x=offset;x<nix;x+=subsample) for (int y=offset;y<niy;y+=subsample) for (int z=offset;z<niz;z+=subsample) {
            // factor : classes
            vec = 0.0f; mat = 0.0f;
			
			// compute the local position
			if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
			else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,1.0f);
			
		
			// compute interpolated values
			priorT = ImageFunctions.linearClosestInterpolation(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
			
			// check if the region is zero: no calculation needed then
			if (priorT>0) {
				// data term : function of the memberships alone
				
				// derivatives
				dPx = ImageFunctions.linearInterpolationXderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
				dPy = ImageFunctions.linearInterpolationYderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
				dPz = ImageFunctions.linearInterpolationZderivative(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
				
				// coordinate derivatives
				dXi = transformModel.imageToTemplateDerivatives(x,y,z,trans,rot,dRa,dRb,dRc,1.0f);
						
				// assemble everything
				for (int i=0;i<Nd;i++) {
					dprior[i] = dPx*dXi[0][i] + dPy*dXi[1][i] + dPz*dXi[2][i];
				}
				
				cost += registrationCost(mems,priorT,x,y,z,k);
				norm += registrationNorm(mems,priorT,x,y,z,k);
				for (int i=0;i<Nd;i++) {
					gradient[i] += registrationCostGradient(mems,priorT,dprior,x,y,z,k,i);
					hessian[i] += registrationCostHessian(mems,priorT,dprior,x,y,z,k,i);
					ngrad[i] += registrationNormGradient(mems,priorT,dprior,x,y,z,k,i);
					nhess[i] += registrationNormHessian(mems,priorT,dprior,x,y,z,k,i);
				}
				
				// add coupling to selected priors
				/*
				for (int l=0;l<classes;l++) if (shapeCoupling[k][l]) {
					float priorC = getTransformedShape(x,y,z,l);
					cost += couplingFactor*couplingEnergy(mems,priorC,priorT,x,y,z,k);
					for (int i=0;i<Nd;i++) {
						hessian[i] += couplingFactor*couplingEnergyHessian(mems,priorC,priorT,dprior,x,y,z,k,i);
						gradient[i] += couplingFactor*couplingEnergyGradient(mems,priorC,priorT,dprior,x,y,z,k,i);
					}
				}
				*/
			}
        }
		// assemble everything, derivatives first
		for (int i=0;i<Nd;i++) {
			hessian[i] = registrationMeasureHessian(hessian[i],nhess[i],norm, normalize);
			gradient[i] = registrationMeasureGradient(gradient[i],ngrad[i],norm, normalize);
		}
		
        return registrationMeasure(cost,norm,normalize);
    } // computeRegistrationCoefficients
    
	/** transformations: how to update the transform using memberships */
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeSingleRegistrationEnergy(float[] trans, int k) {
        float weight;
        float[] Xi;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		float			cost, norm;
		float			normalize = nix*niy*niz;
		
        // set up rotation parameters
        if (transformModel.useRotation()) {
		   rot = transformModel.computeRotation(trans);
		}

		if (precompute) transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,1.0f);
		
		// main loop
		cost = 0.0f;
		norm = 0.0f;		
        for (int x=offset;x<nix;x+=subsample) for (int y=offset;y<niy;y+=subsample) for (int z=offset;z<niz;z+=subsample) {
            // factor : classes
            
			// compute the local position
			if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
			else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,1.0f);
			
			if ( Float.isNaN(Xi[0]) || Float.isNaN(Xi[1]) || Float.isNaN(Xi[2]) ) {
				System.out.print("Pb: X ("+x+","+y+","+z+") -> X' ("+Xi[0]+","+Xi[1]+","+Xi[2]+")\n");	
				System.out.print(displayTransform(trans));	
				System.out.print(displayMatrix(rot));	
				System.out.print(displayMatrix(shapeTransform[k]));	
			}
			
			// compute interpolated values
			priorT = ImageFunctions.linearClosestInterpolation(shape[k],Xi[0],Xi[1],Xi[2],nax,nay,naz);
			
			// check if the region is zero: no calculation needed then
			if (priorT>0) {
				// data term : function of the memberships
				cost += registrationCost(mems,priorT,x,y,z,k);
				norm += registrationNorm(mems,priorT,x,y,z,k);
				
				// add coupling to other memberships ?
				/*
				for (int l=0;l<classes;l++) if (shapeCoupling[k][l]) {
					cost += couplingFactor*couplingEnergy(mems,getTransformedShape(x,y,z,l),priorT,x,y,z,l);
				}
				*/
			}
        }
		
        return registrationMeasure(cost,norm,normalize);
    } // computeRegistrationEnergy
    
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeScaledRegistrationCoefficients(float[] hessian, float[] gradient, float[] trans, float[][][][] mms, float[][][][] shp, int nisx, int nisy, int nisz, int nasx, int nasy, int nasz, float scale) {
        float vec,mat,dp, res,val,num=0,den,weight;
        float xP=0,yP=0,zP=0;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null,dRa = null,dRb = null,dRc = null;
        float[]         dprior = new float[Nd];
		float[] Xi;
        float[][] dXi;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		int				Npt;
		int				n,count;
		float			limit;
		float[]			mem = new float[classes];
		float			cost, norm;
		float[]			ngrad = new float[Nd], nhess = new float[Nd];
		float			normalize = nisx*nisy*nisz;
		
		
        // init the coefficients
		cost = 0.0f;
		norm = 0.0f;
		for (int i=0;i<Nd;i++) {
			gradient[i] = 0.0f;
			hessian[i] = 0.0f;
			ngrad[i] = 0.0f;
			nhess[i] = 0.0f;
		}
		
        // set up rotation parameters
        if (transformModel.useRotation()) {
		   R = transformModel.computeRotationMatrix(trans);
		   rot = R.getMatrix();
		   dRa = R.derivatives(1.0f, 0.0f, 0.0f);
		   dRb = R.derivatives(0.0f, 1.0f, 0.0f);
		   dRc = R.derivatives(0.0f, 0.0f, 1.0f);
		}
		
		if (precompute) for (int k=1;k<classes;k++) if (registeredShape[k]) {
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,scale);
		
			//System.out.println("transfer matrix (shape "+k+")");
			//System.out.println(displayMatrix(shapeTransform[k]));
		}
			
		// main loop
		for (int x=offset;x<nisx;x+=subsample) for (int y=offset;y<nisy;y+=subsample) for (int z=offset;z<nisz;z+=subsample) {
			//System.out.print("-");
            // factor : classes
            vec = 0.0f; mat = 0.0f;
			for (int k=1;k<classes;k++) if (registeredShape[k]) {
				//System.out.print(":");
				// compute the local position
				if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
				else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,scale);
				
				// compute interpolated values
				priorT = ImageFunctions.linearClosestInterpolation(shp[k],Xi[0],Xi[1],Xi[2],nasx,nasy,nasz);
				//if (Xi[0]<0 || Xi[0]>nasx) System.out.print("+");
				
				// check if the region is zero: no calculation needed then
				if (priorT>0) {
					//System.out.print("!");
					// data term : function of the memberships
					
					// coefficient is 1 - coeff/K-1 *sum_l similarity_kl*u_jl^2
					
					// derivatives
					dPx = ImageFunctions.linearInterpolationXderivative(shp[k],Xi[0],Xi[1],Xi[2],nasx,nasy,nasz);
					dPy = ImageFunctions.linearInterpolationYderivative(shp[k],Xi[0],Xi[1],Xi[2],nasx,nasy,nasz);
					dPz = ImageFunctions.linearInterpolationZderivative(shp[k],Xi[0],Xi[1],Xi[2],nasx,nasy,nasz);
					
					// coordinate derivatives
					dXi = transformModel.imageToTemplateDerivatives(x,y,z,trans,rot,dRa,dRb,dRc,scale);
							
					// assemble everiything
					for (int i=0;i<Nd;i++) {
						dprior[i] = dPx*dXi[0][i] + dPy*dXi[1][i] + dPz*dXi[2][i];
					}
					cost += registrationCost(mms,priorT,x,y,z,k);
					norm += registrationNorm(mms,priorT,x,y,z,k);
					for (int i=0;i<Nd;i++) {
						gradient[i] += registrationCostGradient(mms,priorT,dprior,x,y,z,k,i);
						hessian[i] += registrationCostHessian(mms,priorT,dprior,x,y,z,k,i);
						ngrad[i] += registrationNormGradient(mms,priorT,dprior,x,y,z,k,i);
						nhess[i] += registrationNormHessian(mms,priorT,dprior,x,y,z,k,i);
					}
					//System.out.print("c: "+cost);
				}
			}
        }
		//System.out.print("\n");
		//System.out.print("hessian  "+displayVector(hessian));
		//System.out.print("gradient "+displayVector(gradient));
		
        // assemble everything, derivatives first
		for (int i=0;i<Nd;i++) {
			hessian[i] = registrationMeasureHessian(hessian[i],nhess[i],norm, normalize);
			gradient[i] = registrationMeasureGradient(gradient[i],ngrad[i],norm, normalize);
		}
		
        return registrationMeasure(cost,norm,normalize);
    } // computeRegistrationCoefficients
    
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeScaledRegistrationEnergy(float[] trans, float[][][][] mms, float[][][][] shp, int nisx, int nisy, int nisz, int nasx, int nasy, int nasz, float scale) {
        float weight;
        float[] Xi;
        float dPx,dPy,dPz,priorT;
        float[][]       rot = null;
		float[]          regMat,regVec,regParam;
		RotationMatrix  R;
		int				maskId;
		boolean			outside = false;
		float 			cost, norm;
		float			normalize = nisx*nisy*nisz;
		
        // set up rotation parameters
        if (transformModel.useRotation()) {
		   rot = transformModel.computeRotation(trans);
		}

		if (precompute) for (int k=1;k<classes;k++) if (registeredShape[k])
			transformModel.precomputeImageToTemplateMatrix(shapeTransform[k],trans,rot,scale);
		
		// main loop
		cost = 0.0f;
		norm = 0.0f;		
        for (int x=offset;x<nisx;x+=subsample) for (int y=offset;y<nisy;y+=subsample) for (int z=offset;z<nisz;z+=subsample) {
            // factor : classes
            for (int k=1;k<classes;k++) if (registeredShape[k]) {
				// compute the local position
				if (precompute) Xi = fastImageToShapeCoordinates(k,x,y,z);
				else Xi = transformModel.imageToTemplate(x,y,z,trans,rot,scale);
				
				// compute interpolated values
				priorT = ImageFunctions.linearClosestInterpolation(shp[k],Xi[0],Xi[1],Xi[2],nasx,nasy,nasz);
				
				// check if the region is zero: no calculation needed then
				if (priorT>0) {
					// data term : function of the memberships
					cost +=registrationCost(mms,priorT,x,y,z,k);
					norm +=registrationNorm(mms,priorT,x,y,z,k);
				}
			}
        }
		
        return registrationMeasure(cost,norm,normalize);
    } // computeRegistrationEnergy
    
	private final float registrationCost(float[][][][] mem, float pT, int x, int y, int z, int k) {
		return mem[x][y][z][k]*mem[x][y][z][k]*pT*pT;
	}
	
	private final float registrationCostGradient(float[][][][] mem, float pT, float[] dpT, int x, int y, int z, int k, int d) {
		return 2.0f*mem[x][y][z][k]*mem[x][y][z][k]*pT*dpT[d];
	}
	
	private final float registrationCostHessian(float[][][][] mem, float pT, float[] dpT, int x, int y, int z, int k, int d) {
		return 2.0f*mem[x][y][z][k]*mem[x][y][z][k]*dpT[d]*dpT[d];
	}
	
	private final float registrationNorm(float[][][][] mem, float pT, int x, int y, int z, int k) {
		return pT*pT;
	}
	
	private final float registrationNormGradient(float[][][][] mem, float pT, float[] dpT, int x, int y, int z, int k, int d) {
		return 2.0f*pT*dpT[d];
	}
	
	private final float registrationNormHessian(float[][][][] mem, float pT, float[] dpT, int x, int y, int z, int k, int d) {
		return 2.0f*dpT[d]*dpT[d];
	}
	
	private final float registrationMeasure(float cost, float norm, float scaling) {
		//return (cost-norm/(classes*classes))/scaling;
		return cost/norm;
		//return (cost-norm/(classes*classes))/norm;
	}
	
	private final float registrationMeasureGradient(float grad, float ngrad, float norm, float scaling) {
		//return (grad-ngrad/(classes*classes))/scaling;
		return grad/norm;
		//return (grad-ngrad/(classes*classes))/norm;
	}
	
	private final float registrationMeasureHessian(float hess, float nhess, float norm, float scaling) {
		//return (hess-nhess/(classes*classes))/scaling;
		return hess/norm;
		//return (hess-nhess/(classes*classes))/norm;
	}
	
	/**
	 * compute the new transform using gradient descent
	 * performs only one iteration, at scale l
	 */
    final private float registerGradientDescent() {
		float[]		trial = new float[Nd];
		float		E0,E,Eprev;
		boolean		stop = false;
		boolean		admissible = true;
		float[] 	hessian, gradient;
		
		gradient = new float[Nd];
		hessian  = new float[Nd];
		
		// compute the coefficients at current transform
		E0 = computeRegistrationCoefficients(hessian, gradient, transform);
		E = E0;
		
		if (debug) System.out.print( "H: "+displayVector(hessian)+"\n");
		if (debug) System.out.print( "G: "+displayVector(gradient)+"\n");

		// search along the line
		int iter = 0;
		while (!stop) {
			Eprev = E;
			
			// new values for the coeffs 
			for (int n=0; n<Nd; n++) 
				//trial[n] = transform[n] + lambda/(float)Numerics.max(ZERO,hessian[n])*(float)gradient[n];
				trial[n] = transform[n] + lambda/hessian[n]*gradient[n];
			
			// is the energy better ?
			E = computeRegistrationEnergy(trial);
			
			if (debug) System.out.print( "a: "+lambda+"("+E0+")->("+E+")\n");

			// test on energy value : maximisation
			if ( E > E0 ) {
				// better value: changing the transform and the scale
				for (int n=0; n<Nd; n++) 
					transform[n] = trial[n];
				
				lambda = lambda*lfactor;
				stop = true;
				if (debug) System.out.println(displayTransform(transform));
			} 
			else {
				lambda = lambda/lfactor;
				iter++;
				if (iter>itSupp) {
					stop = true;
					if (debug) System.out.print( "stop search\n");
				}
			}
		}
		
		return (E-E0)/E;
	} // registerGradientDescent

	/**
	 * compute the new transform using gradient descent
	 * performs only one iteration, at scale l
	 */
    final private float registerSingleGradientDescent(int k) {
		float[]		trial = new float[Nd];
		float		E0,E,Eprev;
		boolean		stop = false;
		boolean		admissible = true;
		float[] 	hessian, gradient;
		
		gradient = new float[Nd];
		hessian  = new float[Nd];
		
		// compute the coefficients at current transform
		E0 = computeSingleRegistrationCoefficients(hessian, gradient, multitransform[k], k);
		E = E0;
		
		// search along the line
		int iter = 0;
		while (!stop) {
			Eprev = E;
			
			// new values for the coeffs 
			for (int n=0; n<Nd; n++) 
				//trial[n] = multitransform[k][n] + lambda/(float)Numerics.max(ZERO,hessian[n])*(float)gradient[n];
				trial[n] = multitransform[k][n] + lambda/hessian[n]*gradient[n];
			
			// is the energy better ?
			E = computeSingleRegistrationEnergy(trial, k);
			
			if (debug) System.out.print( "a: "+lambda+"("+E0+")->("+E+")\n");

			// test on energy value : maximisation
			if ( E > E0 ) {
				// better value: changing the transform and the scale
				for (int n=0; n<Nd; n++) 
					multitransform[k][n] = trial[n];
				
				lambda = lambda*lfactor;
				stop = true;
				if (debug) {
					System.out.println("("+k+"): "+displayTransform(multitransform[k])+" (E:"+E+")\n");
					if (verbose) System.out.print("("+k+"): "+displayTransform(multitransform[k])+" (E:"+E+")\n");
				} /*else {
					System.out.println(".");
					if (verbose) System.out.print(".");
				}*/					
				
			} 
			else {
				lambda = lambda/lfactor;
				iter++;
				if (iter>itSupp) {
					stop = true;
					if (debug) System.out.print( "stop search\n");
				}
			}
		}
		
		return (E-E0)/E;
	} // registerGradientDescent

	/**
	 * compute the new transform using gradient descent
	 * performs only one iteration, at scale l
	 */
    final private float registerScaledGradientDescent(float[][][][] mms, float[][][][] shp, int nisx, int nisy, int nisz, int nasx, int nasy, int nasz, float scale) {
		float[]		trial = new float[Nd];
		float		E0,E,Eprev;
		boolean		stop = false;
		boolean		admissible = true;
		float[] 	hessian, gradient;
		
		gradient = new float[Nd];
		hessian  = new float[Nd];
			
		// compute the coefficients at current transform
		E0 = computeScaledRegistrationCoefficients(hessian, gradient, transform, mms, shp, nisx, nisy, nisz, nasx, nasy, nasz, scale);
		E = E0;
		
		if (debug) System.out.println( "gradient "+displayVector(gradient));
		if (debug) System.out.println( "hessian  "+displayVector(hessian));
		
		// search along the line
		int iter = 0;
		while (!stop) {
			Eprev = E;
			
			// new values for the coeffs 
			for (int n=0; n<Nd; n++) 
				//trial[n] = transform[n] + lambda/(float)Numerics.max(ZERO,hessian[n])*(float)gradient[n];
				trial[n] = transform[n] + lambda/hessian[n]*gradient[n];
			
			if (debug) System.out.println( "trial "+displayTransform(trial));
		
			// is the energy better ?
			E = computeScaledRegistrationEnergy(trial, mms, shp, nisx, nisy, nisz, nasx, nasy, nasz, scale);
			
			if (debug) System.out.print( "a: "+lambda+"("+E0+")->("+E+")\n");

			// test on energy value : maximisation
			if ( E > E0 ) {
				// better value: changing the transform and the scale
				for (int n=0; n<Nd; n++) 
					transform[n] = trial[n];
				
				lambda = lambda*lfactor;
				stop = true;
				if (debug) {
					System.out.println(displayTransform(transform)+" (E:"+E+")\n");
					if (verbose) System.out.print(displayTransform(transform)+" (E:"+E+")\n");
				} /*else {
					System.out.println(".");
					if (verbose) System.out.print(".");
				}					*/
			} 
			else {
				lambda = lambda/lfactor;
				iter++;
				if (iter>itSupp) {
					stop = true;
					if (debug) System.out.print( "stop search\n");
				}
			}
		}
		
		return (E-E0)/E;
	} // registerScaledGradientDescent

    /** 
	 *	runs a Levenberg-Marquardt step for registering shapes
	 */
	public final void registerShapes() {    
		if (transformMode==MULTIPLE) 
			registerMultiShapes();
		else if (transformMode==DEFORMABLE)
			for (int t=0;t<itMax;t++) demons.registerImageToTarget();
		else registerAllShapes();
	}
	
	private final void registerAllShapes() {	
		boolean stop;
		
        // one level
        if (debug) System.out.println("registration ");
        lambda = INIT_LAMBDA;
		//oldchisq = computeRegistrationCoefficients(hessian, gradient, transform);
		Nturn = 0; itPlus = 0; stop = false;
		//if (debug) System.out.print( "-first--->("+oldchisq+"\n");
		float diff = 1;
		for (int n=0;n<itMax && diff>minEdiff && lambda>minLambda;n++) {
			diff = registerGradientDescent();
		}
		// update the rotation coefficients
		if (transformModel.useRotation()) {
			rotation = transformModel.computeRotation(transform);
		}
		computeTransformedShapeBoundingBox();
		
		//if (changeTemplate) computeTransformedTemplate();
		
    }//registerShapes
 
    /** 
	 *	runs a Levenberg-Marquardt step for registering shapes
	 */
	public final void registerMultiShapes() {    
		boolean stop;
		
		for (int k=0;k<classes;k++) if (registeredShape[k]) {
			// one level
			if (debug) System.out.println("registration ");
			if (debug) System.out.println("\n init ("+k+"): "+displayTransform(multitransform[k]));
			lambda = INIT_LAMBDA;
			//oldchisq = computeSingleRegistrationCoefficients(hessian, gradient, multitransform[k], k);
			Nturn = 0; itPlus = 0; stop = false;
			//if (debug) System.out.print( "-first--->("+oldchisq+"\n");
			float diff = 1;
			for (int n=0;n<itMax && diff>minEdiff && lambda>minLambda;n++) {
				diff = registerSingleGradientDescent(k);
			}
			// update the rotation coefficients
			if (transformModel.useRotation()) {
				multirotation[k] = transformModel.computeRotation(multitransform[k]);
			}
			if (debug) System.out.println("\n final ("+k+"): "+displayTransform(multitransform[k]));
			
			computeSingleTransformedShapeBoundingBox(k);
			//if (changeTemplate) computeTransformedTemplate();
		}
    }//registerShapes
 
    /** 
	 *	runs a pyramid gradient step for registering shapes
	 */
	public final void registerShapesPyramid() {
		boolean stop;
		float[][][][] halfMems;
		float[][][][] halfShape = new float[classes][][][];
		float scale = 1;
		int nisx,nisy,nisz;
		int nasx,nasy,nasz;
		
		// pyramid: scale the image and shapes
		if (debug) System.out.println("registration \n");
		if (debug) System.out.print("structure atlas registration \n");
		
		scale = 1;
		for (int l=1;l<levels;l++) scale = 2*scale;
		for (int l=levels;l>1;l--) {
			if (debug) System.out.println("level "+l+"\n");
			if (debug) System.out.print("level "+l+"\n");
			
			// compute the half images
			nisx = Numerics.floor(nix/scale);
			nisy = Numerics.floor(niy/scale);
			nisz = Numerics.floor(niz/scale);
			halfMems = ImageFunctions.subsample(mems,nix,niy,niz,classes,(int)scale);
			nasx = Numerics.floor(nax/scale);
			nasy = Numerics.floor(nay/scale);
			nasz = Numerics.floor(naz/scale);
			for (int k=0;k<classes;k++) {
				//if (debug) System.out.println("DIMS "+nax+" "+nay+" "+naz+" "+scale);
				 
				halfShape[k] = ImageFunctions.subsample(shape[k],nax,nay,naz,(int)scale);
			}
			
			// perform the gradient descent
			//if (l==levels) subsample = 1; // no subsampling at highest level
			//else subsample = 2*subsample;
			lambda = INIT_LAMBDA;
			//oldchisq = computeExactScaledRegistrationCoefficients(hessian, gradient, transform, halfMems, halfShape, npx, npy, npz, nspx, nspy, nspz, scale);
			Nturn = 0; itPlus = 0; stop = false;
			//if (debug) System.out.print( "-first--->("+oldchisq+"\n");
			//while (!stop) stop = registerLevenbergMarquardt();
			float diff = 1;
			for (int n=0;n<itMax && diff>minEdiff && lambda>minLambda;n++) {
				diff = registerScaledGradientDescent(halfMems, halfShape, nisx, nisy, nisz, nasx, nasy, nasz, scale);
			}
			// update scaling
			scale = scale/2.0f;
		}
		halfMems = null;
		halfShape = null;
		
		// update the rotation coefficients
		if (transformModel.useRotation()) {
			rotation = transformModel.computeRotation(transform);
		}
		computeTransformedShapeBoundingBox();
		
		//if (changeTemplate) computeTransformedTemplate();
		
		System.out.println("\n");
		if (verbose) System.out.print("\n");
		
    }//registerShapesPyramid
 
	/**
	 *	initialize registration parameters
	 */
	public final void initShapeRegistration(float[][][][] mems_,  
											String transformType_, int iter_, int lvl_) {
		// memberships, iterations
		mems = mems_;
		
		itMax = iter_;
		levels = lvl_;
		subsample = 2;
		
		// image parameters
		x0i = nix/2.0f; 
		y0i = niy/2.0f; 
		z0i = niz/2.0f; 
			
		// transform
		if (transformType_.startsWith("multi_")) {
			transformType_ = transformType_.substring(6);
			transformMode = MULTIPLE;
			if (verbose) System.out.println("transform: "+transformType_+" (multi-object)\n");
		} else {
			transformMode = SINGLE;
			if (verbose) System.out.println("transform: "+transformType_+" (single transform)\n");
		}
		transformModel = new ParametricTransform(transformType_, x0i,y0i,z0i, rix,riy,riz, nix,niy,niz, x0a,y0a,z0a, rax,ray,raz, nax,nay,naz);
		
		Nd = transformModel.getDimension();
		
		// init parameters
		if (transformMode==MULTIPLE) {
			multitransform = new float[classes][Nd];
			multirotation = new float[classes][][];
			
			for (int k=0;k<classes;k++) 
				for (int n=0;n<Nd;n++) multitransform[k][n] = 0.0f;
			
			if (transformModel.useRotation())
				for (int k=0;k<classes;k++) 
					multirotation[k] = transformModel.computeRotation(multitransform[k]);
			
			transform = null;
			rotation = null;
		} else {
			transform = new float[Nd];
			
			for (int n=0;n<Nd;n++) transform[n] = 0.0f;
			if (transformModel.useRotation())
				rotation = transformModel.computeRotation(transform);
			
			multitransform = null;
			multirotation = null;
		}
		
		// quadratic scale: no pre-computing :(
		if (transformModel.isLinear()) {
			precompute = true;
			shapeTransform = new float[classes][3][4];
			precomputeTransformMatrix(1.0f);
		} else {
			precompute = false;
			shapeTransform = null;
		}		
	}
		
	/**
	 *	initialize registration parameters
	 */
	public final void updateShapeRegistration(float[][][][] mems_, byte[][][] seg_, String[] modality_,
											  String transformType_, int iter_, int lvl_, float dSmooth_, float dScale_) {
		// memberships, iterations
		mems = mems_;
		
		itMax = iter_;
		levels = lvl_;
		subsample = 2;
		
		// change the transformation type: compute the new parameters
		String oldType = transformModel.getTransformType();
		int newMode;
		if (transformType_.startsWith("deformable_")) {
			transformType_ = transformType_.substring(11);
			newMode = DEFORMABLE;
			if (debug) System.out.println("transform: "+transformType_+" (deformable)\n");
		} else if (transformType_.startsWith("multi_")) {
			transformType_ = transformType_.substring(6);
			newMode = MULTIPLE;
			if (debug) System.out.println("transform: "+transformType_+" (multi-object)\n");
			transformModel = new ParametricTransform(transformType_, x0i,y0i,z0i, rix,riy,riz, nix,niy,niz, x0a,y0a,z0a, rax,ray,raz, nax,nay,naz);
		} else {
			newMode = SINGLE;
			if (debug) System.out.println("transform: "+transformType_+" (single transform)\n");
			transformModel = new ParametricTransform(transformType_, x0i,y0i,z0i, rix,riy,riz, nix,niy,niz, x0a,y0a,z0a, rax,ray,raz, nax,nay,naz);
		}
		String newType = transformType_;
		
		
		// init parameters
		if (newMode==DEFORMABLE) {
			float[][] mat = new float[3][4];
			transformModel.precomputeImageToTemplateMatrix(mat, transform, transformModel.computeRotation(transform), 1.0f);
			if (debug) System.out.println("deformable transform init: "+transformModel.displayTransform(transform)+"\n");

			int type = DemonsToadsAtlasWarping.GAUSS_FLUID;
			if (transformType_.equals("gauss_diffusion")) type = DemonsToadsAtlasWarping.GAUSS_DIFFUSION;
			if (transformType_.equals("gauss_fluid")) type = DemonsToadsAtlasWarping.GAUSS_FLUID;
			if (transformType_.equals("gauss_mixed")) type = DemonsToadsAtlasWarping.GAUSS_MIXED;
			
			// start Demons
			demons = new DemonsToadsAtlasWarping(shape, classes, mems_, null, 0, classes, seg_,
												registeredShape,
												nax, nay,naz, rax, ray, raz, nix, niy, niz, rix, riy, riz,
												dSmooth_, 1.0f, false, 0.0f, dScale_, iter_, iter_,
												type,
												DemonsToadsAtlasWarping.MOVING,
												DemonsToadsAtlasWarping.COMPOSITIVE,
												shapeSlope,
												mat);
			
			demons.initializeTransform();
			
			transformMode = newMode;		
		} else {
			if (newMode==MULTIPLE && transformMode==MULTIPLE) {
				for (int k=0;k<classes;k++)
					multitransform[k] = transformModel.changeTransformType(multitransform[k], oldType, newType);
				transform = null;
			} else if (newMode==MULTIPLE && transformMode==SINGLE) {
				multitransform = new float[classes][Nd];
				for (int k=0;k<classes;k++) {
					for (int n=0;n<Nd;n++) multitransform[k][n] = transform[n];
					multitransform[k] = transformModel.changeTransformType(multitransform[k], oldType, newType);
				}
				transform = null;
			} else if (newMode==SINGLE && transformMode==SINGLE) {
				transform = transformModel.changeTransformType(transform, oldType, newType);
				multitransform = null;
			} else if (newMode==SINGLE && transformMode==MULTIPLE) {
				buildAverageTransformFromMultiTransform();
				transform = transformModel.changeTransformType(transform, oldType, newType);
				multitransform = null;
			} else {
				transform = new float[transformModel.getDimension()];
				for (int n=0;n<transformModel.getDimension();n++) transform[n] = 0.0f;
				multitransform = null;
			}
			// update other parameters
			Nd = transformModel.getDimension();
			transformMode = newMode;
		}
		
		// update rotation parameters
		if (transformMode==MULTIPLE) {
			multirotation = new float[classes][][];
			if (transformModel.useRotation())
				for (int k=0;k<classes;k++) 
					multirotation[k] = transformModel.computeRotation(multitransform[k]);
			
			rotation = null;
		} else if (transformMode==SINGLE) {
			if (transformModel.useRotation())
				rotation = transformModel.computeRotation(transform);
			
			multirotation = null;
		}
		
		// quadratic scale: no pre-computing :(
		if (transformMode!=DEFORMABLE && transformModel.isLinear()) {
			precompute = true;
			shapeTransform = new float[classes][3][4];
			precomputeTransformMatrix(1.0f);
		} else {
			precompute = false;
			shapeTransform = null;
		}
	}
		
	/** 
	 *	computes the transformed coordinates from image to template space
	 */
	private final float[] fastImageToShapeCoordinates(int s, int x,int y,int z) {
		float[] X = new float[3];
		X[0] = shapeTransform[s][0][0]*x + shapeTransform[s][0][1]*y + shapeTransform[s][0][2]*z + shapeTransform[s][0][3];
		X[1] = shapeTransform[s][1][0]*x + shapeTransform[s][1][1]*y + shapeTransform[s][1][2]*z + shapeTransform[s][1][3];
		X[2] = shapeTransform[s][2][0]*x + shapeTransform[s][2][1]*y + shapeTransform[s][2][2]*z + shapeTransform[s][2][3];
		
		return X;
	}
	public final void precomputeTransformMatrix(float scale) {
		if (!precompute) return;
		
		float[][] rot = null;
		if (transformMode==MULTIPLE) {
			for (int s=0;s<classes;s++) {
				if (transformModel.useRotation())
					rot = transformModel.computeRotation(multitransform[s]);
				transformModel.precomputeImageToTemplateMatrix(shapeTransform[s], multitransform[s], rot, scale);
			}
		} else {
			if (transformModel.useRotation())
				rot = transformModel.computeRotation(transform);
			for (int s=0;s<classes;s++) {
				transformModel.precomputeImageToTemplateMatrix(shapeTransform[s], transform, rot, scale);
			}
		}
	}
	
	final public float[][][][] exportRigidDeformationField() { 
		float[][][][] s = new float[3][nix][niy][niz];
		
		if (transformMode==MULTIPLE) buildAverageTransformFromMultiTransform();
		
		float[][] mat =new float[3][4];
		transformModel.precomputeImageToTemplateMatrix(mat, transform, transformModel.computeRotation(transform), 1.0f);
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			s[X][x][y][z] = mat[X][X]*x + mat[X][Y]*y + mat[X][Z]*z + mat[X][T] - x*rix/rax;
			s[Y][x][y][z] = mat[Y][X]*x + mat[Y][Y]*y + mat[Y][Z]*z + mat[Y][T] - y*riy/ray;
			s[Z][x][y][z] = mat[Z][X]*x + mat[Z][Y]*y + mat[Z][Z]*z + mat[Z][T] - z*riz/raz;
		}
		return s;
	}
	
	public final String displayTransform(float[] trans) {
		String info = "transform: (";
		for (int n=0;n<Nd-1;n++) info += trans[n]+", ";
		info += trans[Nd-1]+")\n";
		
		return info;
	}
	
	public final String displayMultiTransform(float[][] trans) {
		String info = "";
		for (int k=0;k<classes;k++) info += displayTransform(trans[k]);
		return info;
	}

	public final String displayVector(float[] vect) {
		String info = "vector: (";
		for (int n=0;n<vect.length-1;n++) info += vect[n]+", ";
		info += vect[vect.length-1]+")\n";
		
		return info;
	}
	
	public final String displayMatrix(float[][] mat) {
		String info = "matrix: (";
		for (int n=0;n<mat.length-1;n++) {
			for (int m=0;m<mat[n].length-1;m++) info += mat[n][m]+", ";
			info += mat[n][mat[n].length-1]+"),(";
		}
		for (int m=0;m<mat.length-1;m++) info += mat[mat.length-1][m]+", ";
		info += mat[mat.length-1][mat[mat.length-1].length-1]+")\n";
		
		return info;
	}
	
}
