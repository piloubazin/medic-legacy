package edu.jhmi.rad.medic.algorithms;

import edu.jhmi.rad.medic.methods.*;
import edu.jhmi.rad.medic.utilities.*;

import java.io.*;


/**
 *	LesionTOADS: TOpology-preserving Anatomy-Driven Segmentation adapted for lesions.
 *	<p>
 *	This is the main version of the segmentation. It includes inhomogeneity correction,
 *	template registration and adaptive coefficients.
 *
 *	@version    Dec 2011
 *	@author    Navid Shiee, Pierre-Louis Bazin 
 *  @see 		JDialogLesionToads
 *	@see 		LesionToads
 *	@see		RegisterTemplate
 *
*/
public class LesionToads {
	public LesionToads() {
	}

	private static final String cvsversion = "$Revision: 1.9 $";
	public static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");

	public String get_version() {
		return revnum;
	}

    // Fuzzy images require 1 image for each class
    // Hard images 1 image with assigned classes
    private float[][]  		destImage;
    private float[][]  		srcImage;
    private int         		destNum=0;
    private String      		resultOutput;
		
    // image size
	private int 		classes; 		// number of classes
	private int 		nx,ny,nz,nc;    // original image dimensions
	private int 		ntx,nty,ntz; 	// original template dimensions
	private	float		rx,ry,rz;		// image resolutions
	private	float		rtx,rty,rtz;	// template resolutions
	private int			orient;			// image orientation
	private int			orx,ory,orz;	// image axis orientation
	private	String[]	modality;		// image modality
	private	String		algorithmMode;		// variants of the algorithm
	private	String		connectivityType;		// choice of connectivity
	
    // segmentation parameters
    private	String				atlasName;
    private		DemonToadDeformableAtlas atlas =null;
   
	private float   	smoothing;
	private float		outlier;
	private int 		iterations;
	private float   	maxDistance;
	private float		bgthreshold;
	private	boolean		doGrowing;
		

	private float		firstLimit;
	private float		lastLimit;
	private float		spread;
	
	private short       maxGMDist;
	private short      maxBstemDist;
	private short       maxVentDist;
	private byte		lesionLabel =10;
	private byte		wmLabel = 25;
	

	private	float		atlasFirstCoefficient, atlasScale;
	private float		demonsSmoothing = 1.0f;
	private float		demonsScale =2.0f;
		

	
	private	String		centroidMode;
	private	float		centroidSmoothness = 0.5f;
	private boolean		lesionWeight = false;
	private String      normType;
	private	String		registrationMode = "rigid";
	
	private	boolean		register = false;
	private	int			levels;
	
	
	private	boolean		alignShapes = true;
	private	int			mainAlignIter = 1;
	private	int			initAlignIter = 10;
	
    private	boolean		correctInhomogeneity = false;
	private	int			polynomialDegree;
	private	float		splineKernel;
	private	int			correctionMethod;
	private float   	lowEratio = 0.1f;
	
	private boolean     includeLesionnClassification = false;
	private boolean     outputClassificationFromMembership = true;
	private boolean     outputField = false;
	private boolean		verbose = true;
	private boolean     debug   = false;
	
    /**
    *	Constructor for 3D images in which changes are placed in a predetermined destination image.
    *   @param destImg_      Image model where result image is to stored.
    *   @param srcImg_       Source image model.
    */
	public LesionToads(int nInput_,
			String[] imgModal_,
			String aName_, String segOutput_, 
			float smooth_, float out_, int nIterMax_, float distMax_, float bgth_,
			boolean outputMaxMem_, 
			float fLim_, float lLim_, float spread_,
			short maxGMDist_, short maxBstemDist_, short maxVentDist_,
			boolean inludeLesions_,
			float atlasCoeff_, float atlasScale_,
			String relMode_, 
			String centMode_, float smoothCentr_,
			boolean lesionWeight_,
			String regMode_,
			String nrmType_,
			boolean register_, 
			int lev_, int initAlignIter_, int mainAlignIter_,
			float dSmooth_, float dScale_,
			boolean correct_, boolean outputField_,
			String correctType_,
			int poly_, float kernel_,
			String algo_,
			String connect_) {
        
	    nc = nInput_;
		modality = imgModal_;
		srcImage = new float[nc][];
        		
		atlasName = aName_;
		atlas = new DemonToadDeformableAtlas(atlasName);
		classes = atlas.getNumber();
		
		resultOutput = segOutput_;
        
		includeLesionnClassification = inludeLesions_;
		smoothing = smooth_;
		outlier = out_;
		iterations = Numerics.abs(nIterMax_);
		doGrowing = (nIterMax_>0);
        maxDistance = distMax_;
		lowEratio = (float)Math.sqrt(maxDistance);
		bgthreshold = bgth_;
		
		outputClassificationFromMembership = outputMaxMem_;
		firstLimit = fLim_;
		lastLimit = lLim_;
		spread = spread_;
		
		maxGMDist = maxGMDist_;
		maxBstemDist = maxBstemDist_;
		maxVentDist = maxVentDist_;
		
		atlasFirstCoefficient = atlasCoeff_;
		atlasScale = atlasScale_;

		centroidMode = centMode_;
		centroidSmoothness = smoothCentr_;
		lesionWeight = lesionWeight_;
		normType = nrmType_;
		
		
		algorithmMode = algo_;
		registrationMode = regMode_;
		
		connectivityType = connect_;
		
		register = register_;
		levels = lev_;
		initAlignIter = initAlignIter_;
		mainAlignIter = mainAlignIter_;
		
		correctInhomogeneity = correct_;
		outputField = outputField_;
		if (correctInhomogeneity) {
			polynomialDegree = poly_;
			splineKernel = kernel_;
			if (correctType_.equals("Splines")) correctionMethod = InhomogeneityCorrection2.SPLINES;
			else correctionMethod = InhomogeneityCorrection2.CHEBYSHEV;
		} else {
			polynomialDegree = 0;
			splineKernel = 0.0f;
			correctionMethod = InhomogeneityCorrection2.CHEBYSHEV;
		}

	}

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; }
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setOrientations(int x, int y, int z, int main) { orx=x; ory=y; orz=z; orient=main; }
	
	public final void setSourceImageAt(int num, float[] val) { srcImage[num] = val; }
	

    /**
    *	Prepares this class for destruction.
    */
	public void finalize(){
	    destImage   = null;
	    srcImage    = null;
        System.gc();
	}

	
    /**
    *	produces the segmentation of the input image
    */
    private void calcSegmentation() {
		byte[][][]   	classification = null;
		byte[][][]   	classification_mem = null;
		float[][]		centroids = null;
		float[][][][] 	memberships = null;
        float[][][][] 	field = null;
        float[][][]   	lesions = null; 
		short[][][]   lesionClassification = null;
		float[]     buffer;
		byte[]		bytebuffer;
		short[]     shortbuffer;
		LesionToadSegmentation 	segmentation=null;
		LesionToadPreprocess		preprocess=null;
		InhomogeneityCorrection2[]	correction=null;
		float		energy, Eprev, Eratio;

		boolean stop;
		
		if (debug) {
			System.out.println("Lesion-TOADS\n");
			System.out.println(getParameterText());
		}
		
		
		//---------------------------------------/
		//-          MAIN SEGMENTATION          -/
		//---------------------------------------/
		
		// record time
		long start_time = System.currentTimeMillis();


        memberships = null;
        classification = null;
		
		// pre-processing
        if (verbose) System.out.println("Initialization\n");
		preprocess = new LesionToadPreprocess(srcImage, nc, nx, ny, nz, rx, ry, rz,
												orient, orx, ory, orz,
												atlas.getTemplate(), ntx, nty, ntz, rtx, rty, rtz,
												classes, (byte)1, 0.0f, 10);
		
		if (debug) System.out.println(preprocess.displayTransform());
		// put all intensities in [0,1]										
		preprocess.normalizeImages();
		// find unused regions
		preprocess.findCroppingBoundaries();
		// crop images and atlas
		preprocess.cropImages();
		if (debug) System.out.println("cropped: "+preprocess.displayTransform());
		int[] dim = preprocess.getCurrentImageDimensions();
		nx = dim[0]; ny = dim[1]; nz = dim[2];
		float[] res = preprocess.getSignedImageResolutions();
		rx = res[0]; ry = res[1]; rz = res[2];
		
		// first centroids, memberships : tissue types
		centroids = preprocess.initialCentroids(atlas.getIntensityPriors(modality,nc),classes);
		
		// align the atlas
		if (alignShapes) {
			if (verbose) System.out.println("Atlas alignment\n");
			memberships = new float[nx][ny][nz][classes];
			preprocess.initialMemberships(memberships,centroids, classes);
			
			atlas.setImageInfo(nx,ny,nz,rx,ry,rz,orient,orx,ory,orz);
			
			
			atlas.initShapeRegistration(memberships, "rigid", initAlignIter, levels);
			atlas.setTransform(preprocess.getTransform());
			
			if (debug) System.out.println("pre-process: "+preprocess.displayTransform());
			
			
			atlas.updateShapeRegistration(memberships, null, modality, "rigid", initAlignIter, levels, 1.0f, 1.0f);
			if (debug) System.out.println("atlas: "+atlas.displayTransform(atlas.getTransform()));
			
			atlas.registerShapesPyramid();
			if (debug) System.out.println("atlas: "+atlas.displayTransform(atlas.getTransform()));
			
			memberships = null;
			
			if (debug) System.out.println("aligned atlas: "+atlas.displayTransform(atlas.getTransform()));
			if (debug) System.out.println("aligned pre-process: "+preprocess.displayTransform());
			
			preprocess.uncropImages();
			preprocess.findCroppingBoundaries();
			preprocess.cropImages();
			
			if (debug) System.out.println("final pre-process: "+preprocess.displayTransform());
			if (debug) System.out.println("final aligned atlas: "+atlas.displayTransform(atlas.getTransform()));
			
			dim = preprocess.getCurrentImageDimensions();
			nx = dim[0]; ny = dim[1]; nz = dim[2];
			
			atlas.setImageInfo(nx,ny,nz,rx,ry,rz,orient,orx,ory,orz);
			
			atlas.computeTransformedTemplate();
			
			// re-compute the centroids to get better estimates
			preprocess.updateTransformedTemplate(atlas.getTemplate());
		
		}
		
		// create the segmentation algorithm
		segmentation = new LesionToadSegmentation(preprocess.getImages(), atlas, 
								nx, ny, nz, nc, 
								smoothing, outlier,
								firstLimit, lastLimit, 
								0, spread, 
								maxGMDist, maxBstemDist, maxVentDist,
								atlasFirstCoefficient, atlasScale,
								modality,
								algorithmMode, normType, connectivityType,
								centroidMode, centroidSmoothness);

		// get the centroids from initialization
		segmentation.setDistanceFactor(0);
		segmentation.initCentroids(centroids);
		segmentation.setLesionCentroid(preprocess.initialLesionCentroid(atlas.getLesionPriors(modality,nc)));
		segmentation.setBlackHoleCentroid(preprocess.initialBlackHoleCentroid(atlas.getBlackHolePriors(modality,nc)));
		if (lesionWeight) 
			segmentation.setLesionWeight();
		else
			segmentation.reSetLesionWeight();
		segmentation.setIntensityMax(preprocess.getIntensityMax());
		segmentation.setIntensityScale(preprocess.getIntensityScale());
		segmentation.reSetGMflag();
		if (verbose) System.out.println("initial "+segmentation.displayCentroids()); 
		if (verbose) System.out.println("initial "+segmentation.displayLesionCentroid());
		//if (verbose) System.out.println("initial "+segmentation.displayBlackHoleCentroid());
		
		// also use the initial memberships (re-computed)
		preprocess.initialMemberships(segmentation.getMemberships(), centroids, classes);

		// re-initialize the atlas alignment
		atlas.updateShapeRegistration(segmentation.getMemberships(), segmentation.getSegmentation(), 
						modality, registrationMode, mainAlignIter, 1, demonsSmoothing, demonsScale);
		// non-rigid shape alignment ?
		atlas.registerShapes();
		if (verbose) System.out.println("final aligned atlas: "+atlas.displayTransform(atlas.getTransform()));
			
		
		if (correctInhomogeneity) {
			//fireProgressStateChanged("initialization (inhomogeneity)");
				
			correction = new InhomogeneityCorrection2[nc];
			for (int c=0;c<nc;c++) {
				// create the inhomogeneity correction segmentation
				correction[c] = new InhomogeneityCorrection2(preprocess.getImages()[c],
						segmentation.getMemberships(),
						segmentation.getCentroids(c),
						classes,
						segmentation.getSegmentation(),atlas.getTopology(),
						polynomialDegree, splineKernel,
						InhomogeneityCorrection2.IMAGE,
						InhomogeneityCorrection2.GLOBAL,
						correctionMethod,
						InhomogeneityCorrection2.FIXED,
						InhomogeneityCorrection2.EQUAL_VARIANCE,
						2.0f, null, null,
						nx,ny,nz,rx,ry,rz);

				// integrate it into the segmentation
				segmentation.addInhomogeneityCorrection(correction[c].getField(),c);
			}
		} else correction = null;
		
		if (debug) System.out.print("pre-processing time: " + (System.currentTimeMillis()-start_time)+"\n"); 

		// relations
		segmentation.computeRelations();
		
		if (debug) System.out.print("relation initialization time: " + (System.currentTimeMillis()-start_time)+"\n"); 
		
		segmentation.computeVariances();
		if (debug) System.out.println(segmentation.displayVariances());
		segmentation.computeChannelWeights();
		energy = segmentation.computeMemberships();
		segmentation.setGMflag();
			
		if (debug) System.out.print("membership initialization time: " + (System.currentTimeMillis()-start_time)+"\n"); 
			
		Eprev = energy;
		Eratio = 1.0f;
			
		if (debug) System.out.print("final initialization time: " + (System.currentTimeMillis()-start_time)+"\n"); 
			
		if (debug) System.out.print("algorithm mode: " +algorithmMode+"\n"); 
		
		// 3. main loop
		//--------------------------------------------------------------------
		stop = false;
		int iter = 0;
		int deg = 0;
		//segmentation.setDistanceFactor(50.0f);
		if (iterations==0) stop=true;
		
		segmentation.setDistanceFactor(10.0f);
		while ((!stop) && segmentation.isWorking()) {
			long iterationTime = 0;
			if (debug) iterationTime = System.currentTimeMillis();
				
			iter++;
			if (verbose) System.out.println("iteration "+iter+"\n");
			
			// a. thinning
			//----------------------------------------------------
			
			if (algorithmMode.equals("narrow_band_competition")) {
				segmentation.propagateCompetition();
			} else {
				segmentation.propagateThinning();
			}
			
			if (debug) System.out.print("thinning time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
			if (debug) iterationTime = System.currentTimeMillis();
				
			// b. update parameter information
			//------------------------------------------------------------
			
			// inhomogeneity field
			if (correctInhomogeneity) {
				if ( (deg<polynomialDegree) && (Eratio < lowEratio) ) deg++;
				for (int c=0;c<nc;c++) correction[c].computeCorrectionField(deg);
							
				if (debug) System.out.print("inhomogeneity time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
				if (debug) iterationTime = System.currentTimeMillis();
			}
			
			
				
			// centroids

			segmentation.computeCentroids();
			if (verbose) System.out.println(segmentation.displayCentroids());
			if (verbose) System.out.println(segmentation.displayLesionCentroid());
			//if (verbose) System.out.println(segmentation.displayBlackHoleCentroid());
			
			// cluster-specific variations
			segmentation.computeVariances();
			if (debug) System.out.println(segmentation.displayVariances());
			if (debug) System.out.println(segmentation.displayLesionVariance());

			//Channel Weights 
			segmentation.computeChannelWeights();
			if (debug) System.out.println(segmentation.displayChannelWeights());
						
			if (debug) System.out.print("centroid time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
			if (debug) iterationTime = System.currentTimeMillis();
				
			// memberships & energy
			Eprev = energy;
			energy = segmentation.computeMemberships();
			Eratio = Numerics.abs( 2.0f*(Eprev-energy)/(Eprev+energy) );
			
			if (debug) System.out.print("membership time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
			if (debug) iterationTime = System.currentTimeMillis();
			
			
			
			// shape alignment
			if (alignShapes) atlas.registerShapes();
				
			if (debug) System.out.print("alignment time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
			if (debug) iterationTime = System.currentTimeMillis();
			
			// check for convergence
			if ( Eratio < maxDistance) stop = true;
			if (verbose) System.out.println("energy: "+energy+", " +"ratio: "+Eratio+"\n");
			
			if (iter >= iterations) stop = true;
				
			// c. growing
			//---------------------------------------------------------------
			if (algorithmMode.equals("narrow_band_competition")) {
				// do nothing
			} else {
				if ( (!stop) || doGrowing){
					segmentation.propagateGrowing();
				}
				
				if (debug) System.out.print("growing time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
				if (debug) iterationTime = System.currentTimeMillis();
			}
					
				
			// d. update classification only at the end
			//--------------------------------------------------------------
			if (stop) {
				segmentation.updateClassificationFromLabels();
				segmentation.updateLesionClassifcation();
				segmentation.finalizeSegmentation(rx*ry*rz);
			}
			
			// update the relation list
			if (!stop){
				segmentation.computeRelations();
			}

			if (debug) System.out.print("update time: " + (System.currentTimeMillis()-iterationTime)+"\n"); 
			if (debug) iterationTime = System.currentTimeMillis();
		}
		
		if (verbose){
			float secs = (float)((System.currentTimeMillis()-start_time)/1000.0f);
			int mins = Numerics.floor(secs/60.0f);
			secs -= 60*mins;  
			System.out.println("total processing time : " + mins + " minute(s) " + Numerics.round(secs) + " second(s) \n" ) ;
			System.out.println("Preparing output images\n");
		}
		segmentation.cleanUp();
		
		classification = segmentation.exportClassification();
		lesionClassification = segmentation.exportLesionClassification();
		if (outputClassificationFromMembership){
			segmentation.updateClassificationFromMemberships();
			classification_mem = segmentation.exportClassification();
		}
		if (resultOutput.equals("hard segmentation+memberships")) {
			memberships = segmentation.exportMemberships();
			lesions = segmentation.getLesionMemberships();
			
		}else if (resultOutput.equals("cruise inputs")){
			memberships = segmentation.generateCruiseInputs();
		} else if (resultOutput.equals("dura removal inputs")){
			memberships = segmentation.generateDuraRemovalOutputs();
		}
		if (correctInhomogeneity && outputField) {
			field = new float[nc][][][];
			for (int c=0;c<nc;c++) field[c] = correction[c].exportField();
		}
		
		
		// clean-up the algorithms
	
		if (segmentation!=null) segmentation.finalize();
		segmentation = null;
		if (correction!=null) {
			for (int c=0;c<nc;c++) if (correction[c]!=null) correction[c].finalize();
			correction = null;
		}
		for (int x =0; x<nx; x++) for (int y =0; y<ny; y++) for (int z =0; z<nz; z++){
			if (classification[x][y][z] > wmLabel) classification[x][y][z] =wmLabel;
			if (outputClassificationFromMembership)
				if ( classification_mem[x][y][z]>wmLabel) classification_mem[x][y][z]=wmLabel;
		}
		if (includeLesionnClassification)
			for (int x =0; x<nx; x++) for (int y =0; y<ny; y++) for (int z =0; z<nz; z++){
				if (lesionClassification[x][y][z]==1) classification[x][y][z] = lesionLabel;
				if (outputClassificationFromMembership)
					if (classification_mem[x][y][z]==wmLabel & lesionClassification[x][y][z]==1)
					classification_mem[x][y][z] = lesionLabel;
			}
		
		try {
			// Calculate the number of result images.
			if (resultOutput.equals("hard segmentation+memberships")) destNum = 3; // mems + hardClass + LesionHard  
			else if (resultOutput.equals("cruise inputs")||resultOutput.equals("dura removal inputs")) destNum = 6;
			else destNum = 2; // classification
			
			if (outputClassificationFromMembership) destNum++;//for classification from mem
			if (outputField) destNum++; //for inhomogeneity field
			
			destImage = new float[destNum][];
			
			// output dimensions
			int Nout=0;
            
			
			if (resultOutput.equals("all_images")
				|| resultOutput.equals("hard segmentation+memberships")
				) {
			       	destImage[Nout] = new float[(classes+1)*preprocess.getOriginalImageSize()];
			       	for (int k=0;k<classes;k++) {
			       	    buffer = preprocess.uncropAndBuffer(memberships[k],0.0f);
			       	    for (int n=0;n<buffer.length;n++) {
			       	        destImage[Nout][k*preprocess.getOriginalImageSize()+n] = buffer[n];
			       	    }
			       	}
            	buffer = preprocess.uncropAndBuffer(lesions,0.0f);
            	for (int n=0;n<buffer.length;n++) {
                    destImage[Nout][classes*preprocess.getOriginalImageSize()+n] = buffer[n];
                }
                memberships = null;
            	buffer = null;
            	Nout++;
			}
            
          	destImage[Nout] = new float[preprocess.getOriginalImageSize()];				

            bytebuffer = preprocess.uncropAndBuffer(classification,(byte)0);
            for (int n=0;n<bytebuffer.length;n++) {
                destImage[Nout][n] = bytebuffer[n];
            }	
            Nout++;
				
            if (outputClassificationFromMembership){
                destImage[Nout] = new float[preprocess.getOriginalImageSize()];				
                bytebuffer = preprocess.uncropAndBuffer(classification_mem,(byte)0);
				for (int n=0;n<bytebuffer.length;n++) {
                    destImage[Nout][n] = bytebuffer[n];
                }	
            	Nout++;
            }

//				lesion classification
            destImage[Nout] = new float[preprocess.getOriginalImageSize()];
            shortbuffer = preprocess.uncropAndBuffer(lesionClassification,(short)0);
            for (int n=0;n<shortbuffer.length;n++) {
                destImage[Nout][n] = shortbuffer[n];
            }
            Nout ++;
				
            buffer = null;
            bytebuffer = null;
            shortbuffer = null;
            classification = null;
            lesionClassification = null;
				
			

			if ( (resultOutput.equals("cruise inputs")) || (resultOutput.equals("dura removal inputs"))) {
				
				destImage[Nout] = new float[preprocess.getOriginalImageSize()];
				destImage[Nout+1] = new float[preprocess.getOriginalImageSize()];
				destImage[Nout+2] = new float[preprocess.getOriginalImageSize()];
				destImage[Nout+3] = new float[preprocess.getOriginalImageSize()];
				
				for (int i=0;i<4;i++){
					buffer = preprocess.uncropAndBuffer(memberships[i],0.0f);
					for (int n=0;n<buffer.length;n++) {
                        destImage[Nout+i][n] = buffer[n];
                    }
            	}
				Nout += 4;
				buffer = null;
			}
			
          
			if (outputField && correctInhomogeneity){
            	if (nc>1) {
					destImage[Nout] = new float[nc*preprocess.getOriginalImageSize()];
					for (int c=0;c<nc;c++) {
						buffer = preprocess.uncropAndBuffer(field[c],1.0f);
						for (int n=0;n<buffer.length;n++) {
                            destImage[Nout][c*preprocess.getOriginalImageSize()+n] = buffer[n];
                        }
					}
				} else {
					destImage[Nout] = new float[preprocess.getOriginalImageSize()];
					buffer = preprocess.uncropAndBuffer(field[0],1.0f);
					for (int n=0;n<buffer.length;n++) {
                        destImage[Nout][n] = buffer[n];
                    }
				}
        
		
			}
            
			buffer = null;
			
		} catch (OutOfMemoryError e) {
			buffer = null;
			System.err.println("Segmentation Lesion-TOADS: Out of memory creating hard classification");
			System.err.println(e.getMessage());
			finalize();
			
			return;
		} 
        if (debug) System.out.println("output...\n");

		// clean-up the algorithms
		if (preprocess!=null) preprocess.finalize();
		preprocess = null;
		
		memberships = null;
		classification = null;
		lesions = null;
		
		
    } // calcSegmentation
	
    public float[] getResultImageAt(int num){
    	return destImage[num];
    }
    
   
    /**
     * Construct a readable list of the parameters to this segmentation.
     * @return       the parameter string
     */
	public String getParameterText() {
		String delim = "\n";
        String str = new String();
		for (int c=0;c<nc;c++) {
			str += "image "+(c+1)+" modality: "+modality[c]+ delim;
		}
        str += "atlas name: " + atlasName + delim;
        str += "results: "+ resultOutput + " ("+destNum+") "+delim;
		
		str += "include lesion in classification: "+ includeLesionnClassification + delim;
		str += "smoothing: " + smoothing + delim;
        str += "outliers: "+ outlier + delim;
        str += "max difference: "+ maxDistance + delim;
        str += "iterations: "+ iterations + delim;
        str += "background threshold: "+ bgthreshold + delim;
        
        str += "max GM distance: " + maxGMDist + delim;
        str += "max Ventricle distance: "+ maxVentDist + delim;
        str += "max Brainstem distance: "+ maxBstemDist + delim;
		
		//str += "use topology: "+ useTopologyPrior + delim;
		str += "first limit: "+ firstLimit + delim;
        str += "last limit: "+ lastLimit + delim;
        str += "spread: "+ spread + delim;
		
		str += "output classification from membership: "+ outputClassificationFromMembership + delim;
		str += "atlas first coefficient: "+ atlasFirstCoefficient + delim;
		str += "atlas scale: "+ atlasScale + delim;
		
		//str += "use intensity prior: "+ useIntensityPrior + delim;
		
		str += "centroid mode: "+ centroidMode + delim;
		str += "centroid smoothing: "+ centroidSmoothness + delim;
		str += "use Lesion Weights: " + lesionWeight + delim;
		
		str += "Norm mode: " + normType + delim;
		
		str += "register: "+ register + delim;
		str += "registration iteration (init,main): " + initAlignIter +", "+ mainAlignIter + delim;
        str += "registration levels: "+ levels + delim;
        
		str += "inhomogeneity correction: "+ correctInhomogeneity + delim;
		str += "output inhomogeneity field: " + outputField + delim;
		str += "correction method: "+ correctionMethod + delim;
        str += "degree: "+ polynomialDegree + delim;
        str += "kernel: "+ splineKernel + delim;
        
        return str;
    }

	/** access the atlas information */
	public final DemonToadDeformableAtlas getAtlas() { return atlas; }
}
