/*****************************************************************************/
/*                                                                           */
/*       		RETINO_CORTICAL DYNAMICS MODEL                               */
/*                                                                           */
/*                           VERSION  1.1                                    */
/*                                                                           */
/*                        RETINAL  PARAMETERS                                */
/* HO 08/20/99                                                               */
/*****************************************************************************/


#define Pi 3.14159
#define numxon  580
//#define Zcell_x     556
//#define Zcell_y     (2*556)
#define numxonper       numxon
#define numxoff	 	0
#define numy 	 	numxon /*if set to 0, set also InAmp to 0*/
#define nums 	 	numxon

/******* Equation Parameters *********************/

/************** Transmitter equation constants For Y cells *****************/

#define alphay	8.4   //0.4   // When using this variable remember to comment out the variable alpha declared in the main program.
#define betay	40.0  //16.0  // if you increase this value you will get a response. As it is it does not produce any output.
#define gamay    1.53  //0.13
#define Jvaly	12.0  //12.0

/************** Transmitter equation constants For X cells *****************/

#define alpha	0.4   //0.4   // When using this variable remember to comment out the variable alpha declared in the main program.
#define beta	16.0  //16.0  // Beta controls the overshoot. The larger this value the more transient the overshoot and undershoot will be.
#define gama    0.13  //0.13  // Gamma controls the decay.. the lower this value the quicker the undershoot is going to decay.
#define Jval	12.0  //12.0


/*************** X ON cell equation constants *********************/

#define Axon	2.0    // 2.0
#define Bxon    250.0   // 250.0
#define Dxon    10.0    // 10.0
#define Jvalxon 6.0     // 6.0


       /**** These params are from Kaplan:  Vis Res, 1994,pp 7 **/
       /*** 23 arcsec is inter-receptor spacing Coletta-Williams'87
            Dacey'93 ***/

#define XGcVar	(6*0.03*60.0*60.0)/23.0   // 0.03  /*center extent  (7)*/
#define XGsVar	(6*0.1*60.0*60.0)/23.0    /* 0.18  surround extent (4)*/  //original --> (0.18*60.0*60.0)/23.0

      /* Amps are calculated as the ratio 4.4/325.0     */
#define XGcAmp	1.0 	/* amplitude of center Gaussian */
#define XGsAmp	0.0  /*  0.0135  amplitude of surround Gaussian */

#define halfx	28 	/* (28)		 */
#define fullx   (halfx*2)
//#define timescale1   100*0.0035  /*scales the timing of X with respect to cortex */
#define timescale   0.0035  /*scales the timing of X with respect to cortex */
#define xydelay   2.0  /*delay between the X and Y cell (X 8 gives ms) (was 8;
			was 6; now reduced to 2: p. 16) */
#define xcsdelay  2.0 /*delay between the center and surround of the X (2.0)*/

/**************X ON PER persistence parameter ********************/

//#define persistence 0.035 /* persistence parameter */

/*************** Y cell equation constants ******************/
#define	Ay 	2.0	/* (200.8) (1.0) */
#define By	600.0   /* (100.0)*/
#define delay	2.0    /* tau=1; duration of response (25.0) (6.0)*/ //I had to increase the delay to 2 units so as to generate the transient response
/* params from Kaplan */
#define YGcVar	(0.09*60*60)/23.0	/* Y center extent */ //0.1*60*60)/23.0
#define YGsVar	(0.72*60*60)/23.0	/* Y surround extent */ //(0.72*60*60)/23.0
#define YGcAmp	1.0	/* amplitude of the center Gaussian */
#define YGsAmp	0.00743	/* amplitude of the surround Gaussian */ // 0.00743
/*Because surround is very weak, only the center is used in the simulation*/
#define halfy	40	/* half receptive field */ // 40
#define fully 	(halfy*2)
#define yth     0.0     /*threshold for the y cell (1.0)*/
#define mgain	1.0     /* Gain in M's Michaelis-Menten eqn */
#define mhalfsat 50    /*Half-saturation in M's Michaelis-Menten eqn*/

/******* background luminance *********/

#define L0 10.0		/*background value (1.0)*/

#define Sigma 0.0

#define Kp 1.0

/******* two flash stimuli parameters *********/
#define dfM1 	60	 /*M*/
#define dfT	90	 /*T 135*/
#define dfM2	60	 /*M (180)*/
#define dfsizeT	9	 /*define spatial half extend of the T (21)*/
#define dfsizeM	9	 /*define spatial half extend of the M (21) */
#define dfdurationT	2.0          /*810.0*/
#define dfdurationM	2.0 	 /* duration of the flashes */
#define dfsoa	45.0	 /* stimulus onset asynchrony */
#define dfintT	1.0     /* intensity (1.0) */
#define dfintM	1.0


/******* three flash stimuli parameters *********/
#define tfM1 		210	 /*M1 (90;75;210)*/
#define tfT		250	 /*T 135 (150;200)*/
#define tfM2		65	 /*M2 180 (70;35)*/
#define tfsizeT		9	 /*define spatial half extend of the T 21 (38)*/
#define tfsizeM1	9	 /*define spatial half extend of the M1 21 (9)*/
#define tfsizeM2	9	 /*define spatial half extend of the M2 21 (9) */
#define tfdurationT	1.25 //1.25          /* (2) 810.0*/
#define tfdurationM1    1.25 //1.25 	 /* duration of the flashes (2) */
#define tfdurationM2	1.25 //1.25	 /* duration of the flashes */
#define tfNsoa          7    /*Number of SOAs in the simulation 13 17 */
#define NSF         7 // total number of Spatial frecuencies to be tested
/* stimulus onset asynchronies for M1*/


/* #define tfsoa1		0.0	/* was 10.0 */

// the values correspond to the temporal frequencies [0.13, 0.5, 1, 2, 5, 8, 16]
//double As[7] = {1.0058817606809214, 1.0109869426305931, 1.0179862099620915, 1.0474258731775667, 1.5, 1.9525741268224333, 2.0};

//double tfvsoa1[tfNsoa]	= {-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0};
/*double tfvsoa1[tfNsoa]	= {0.0,2.0,4.0,6.0,7.5,10.0,13.0,15.0,20.0,30.0};*/

//double SFs[11] = {0.0, 0.143, 0.429, 0.715, 1.57, 2.15, 3.29, 4.43, 6.72, 9.01, 18.16}; //{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; //{0.0, 0.143, 0.429, 0.715, 1.57, 2.15, 3.29, 4.43, 6.72, 9.01, 18.16};
//double TF[11] = {1.00, 2.0, 4.0 , 7.0, 10.0,12, 15.0, 20.0, 30.0, 40.0, 50.0};
//double B[13] = {0.0146, 0.0264, 0.0337, 0.0498, 0.0674, 0.0968, 0.127, 0.174, 0.253, 0.352, 0.5, 0.7, 1.0}; // contrast values for 1.22 hz for Lee's 1990 paper.
//double B[7] = {0.0146, 0.0337, 0.0674, 0.0968, 0.127, 0.174, 0.253};
//double count[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

//double alphas[2] = {0.4, 90}; //{0.4, 1, 3, 7, 10, 12, 15, 20, 25, 30, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150};


#define tfsoa2		0.0	 /* stimulus onset asynchrony for M2 */

/*double tfvsoa2[tfNsoa]	= {0.0,-4.0,-8.0,-11.0,-15.0,-20.0,-23.0,-30.0};*/



#define tfintT		1.0     /* intensity (1.0) */
#define tfintM1		1.0	/* (1.0) put to 0 for target only*/
#define tfintM2		0.0

#define Intensity 1.0


/******* flicker stimulus parameters *********/
/**use same params as 2 flash masking for location, duration, etc. **/
#define flickperiodT	33.32 /*period of the flicker for Target
                                           written as 1000/(8*freq(Hz)
                                           7.5 -> 16.67)
                                           use 1000/(4*freq(Hz))
                                           7.5 -> 33.32*/
#define flickperiodM	(1000.0)/(8.0*7.5) /*period of the flicker for Mask*/

#define myscale   0.35    /* facilitation equation parameters timesclae and spatial spread*/
#define halfs  40

/**********retcor map function****************/
//#define retthresh   246.36 //231.64  240.5  238.5
/*  #define thresh 98.05   */
/*  #define thresh 104.9   */


/******************************************************/

/** stimulus resulting from the convolution of an edge with a  **/
/** certain base blur with the line spread function of the eye **/
/** approximated by a Gaussian of std dev 2 arcmin. The values **/
/** of the input are stored in ConvolvedA[] to avoid recalculating **/
/** the input for every input() call.                          **/

#define dim 3*2*46*10   /** 3std dev*(2*max base blur )*
                          (# sampling pts per arcmin) **/
#define hdim (int)dim/2-1

/******************************************************/

