#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunmatrix/sunmatrix_dense.h>            /* access to dense SUNMatrix                    */
#include <sunlinsol/sunlinsol_dense.h>            /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h>             /* access to band SUNMatrix                     */
#include <sunlinsol/sunlinsol_band.h>             /* access to band SUNLinearSolver               */
#include <cvode/cvode_diag.h>                     /* access to CVDIAG linear solver               */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */
#include <sundials/sundials_types.h>
//#include </home/coglab/Documents/Retina/Retina/retprms.h>
#include </home/coglab/Documents/Retina/Demo_newCvode/retina.h>
            /* definition of realtype                       */

/* helpful macros */

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

/* Shared Problem Constants */

#define ATOL RCONST(1.0e-6)
#define RTOL RCONST(5.0e-8)

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define TWO    RCONST(2.0)
#define THIRTY RCONST(30.0)

/* Problem #1 Constants */

#define P1_NEQ        1160
#define P1_ETA        RCONST(3.0)
#define P1_NOUT       1000
#define P1_T0         RCONST(0.0)
#define P1_T1         RCONST(0.25)
#define P1_DTOUT      RCONST(0.25)
#define P1_TOL_FACTOR RCONST(1.0e4)

/* Problem #2 Constants */

#define P2_MESHX      5
#define P2_MESHY      5
#define P2_NEQ        P2_MESHX*P2_MESHY
#define P2_ALPH1      RCONST(1.0)
#define P2_ALPH2      RCONST(1.0)
#define P2_NOUT       5
#define P2_ML         5
#define P2_MU         0
#define P2_T0         RCONST(0.0)
#define P2_T1         RCONST(0.01)
#define P2_TOUT_MULT  RCONST(10.0)
#define P2_TOL_FACTOR RCONST(1.0e3)

#define rtimeinterval 0.25
#define numcell 580
#define rduration 250
#define retcorgain 1
#define Bp 1
#define Ap 1
#define Bt 1
#define At 10


/********** Global variable ***************/
//int sxzon,sxzons,sxon,sxonper,sxzoff,sxoff,sy,scell, syzon, syzons;
int xon, yon, Pch, Mch;
double gc[fullx+1],gs[fullx+1],gyc[fully+1], gys[fully+1], ConvolvedA[6000];
double RFyon[fully+1],RFyoff[fully+1];
double XONpthresh, XONnthresh,XOFFpthresh,XOFFnthresh,Cz,equxz,
equxon,equxoff,equy, equyz;/* equivalum point of x cell z */
double equcortex,tsum,tfsoa1, Sfreq, Sp, DeltaB, tf;
double gthresh();
/*Modification to allow -ve SOAs in 3 flash stimuli HO 8/8/01 */
/* start and end positions for 3 flash stimuli*/
int sfM11, efM11, sfM12, efM12, sfM21, efM21;
int sfM22, efM22, sfT, efT;
/* start and end positions for 2 flash and flicker stimuli*/
int sfdM1, efdM1, sfdM2, efdM2, sfdT, efdT;
double ttimeshift = 0.0; /*parameter to shift starting times to allow -ve SOAs*/

float StimM[(int) (rduration / rtimeinterval +1)][numcell];
float StimP[(int) (rduration / rtimeinterval +1)][numcell];
/* Linear Solver Options */
enum {FUNC, DENSE_USER, DENSE_DQ, DIAG, BAND_USER, BAND_DQ};

/* Private Helper Functions */

static int  Problem1(void);
static void PrintIntro1(void);
static void PrintHeader1(void);
static void PrintOutput1(realtype t, realtype y0, realtype y1, int qu, realtype hu);
static realtype MaxError(N_Vector y, realtype t);
static int PrepareNextRun(void *cvode_mem, int lmm, int miter, N_Vector y,
                          SUNMatrix* A, sunindextype mu, sunindextype ml,
                          SUNLinearSolver* LS, SUNNonlinearSolver* NLS);
static void PrintErrOutput(realtype tol_factor);
static void PrintFinalStats(void *cvode_mem, int miter, realtype ero);
static void PrintErrInfo(int nerr);

/* Functions Called by the Solver */

static int f1(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac1(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Implementation */

int main(void)
{
  int nerr;

  nerr = Problem1();
  PrintErrInfo(nerr);

  return(0);
}

static int Problem1(void)
{
  realtype reltol=RTOL, abstol=ATOL, t, tout, ero, er;
  int miter, retval, temp_retval, iout, nerr=0;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
  void *cvode_mem;
  booleantype firstrun;
  int qu;
  realtype hu;
  int i, j, k, stimuluspoint;
  double Zcellxeq, Zcellyeq, Xcelleq, Ycelleq;
  int Pch, Mch;
  double SumGe, SumGi, SumGey, SumGiy;
  double Attention(), Xthresh(), Ythresh(), inputS(), inputT();
  void newgauss();
  FILE *fopen(),*cortexPdat,*cortexMdat, *Pdatafile, *Mdatafile, *Stimuli;

  y = NULL;
  A = NULL;
  LS = NULL;
  NLS = NULL;
  cvode_mem = NULL;

  /*File names to store the data for the Spatial Resolution Exp.*/
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/positiveStep/gap6/Attention/";
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/positiveStep/gap6/Neutral/";

  /*File names to store the data for the Temporal Resolution Exp.*/
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp/isi3_/Attention/";
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp/isi3_/Neutral/";
  /*File names to store the data for the Temporal Resolution Exp. 3*/
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp3/isi3/Attention/";
  char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp3/isi3/Neutral/";
 /*File names to store the data for the Temporal Resolution Exp. 4*/
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp4/isi3/Attention/";
  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/TemporalExp4/isi3/Neutral/";

  char pathfile[100];
  //char pathfileint[100] = "./data";

  /***********INITIAL CONDITIONS ******/
  Pch = numcell;
  Mch = numcell + numcell;//numxon + numxon;
  xon = numxon + numxon + numxon;
  yon = numxon + numxon + numxon + numxon;
  printf("%i\n", xon);
// Gaussians
  newgauss(gc, XGcAmp, XGcVar, halfx); // C-S RFs for X cells
  newgauss(gs, XGsAmp, XGsVar, halfx);
  newgauss(gyc, YGcAmp, YGcVar, halfy); // C-S RFs for Y cells
  newgauss(gys, YGsAmp, YGsVar, halfy);
  SumGe=0; SumGi=0;
  for (i=0;i<=fullx;i++)
    {
		SumGe=SumGe+gc[i];
		SumGi=SumGi+gs[i];
		//printf("SumGey = %e\n",SumGe);
	}
  SumGey=0; SumGiy=0;
  for (i=0;i<=fully;i++)
    {
		SumGey=SumGey+gyc[i];
		SumGiy=SumGiy+gys[i];
		//printf("SumGey = %e\n",SumGe);
	}
 /******************open the output data files************************/

strcpy(pathfile,root);
strcat(pathfile, "cortexP.dat");
cortexPdat  = fopen(pathfile,"w");

strcpy(pathfile,root);
strcat(pathfile,"cortexM.dat");
cortexMdat = fopen(pathfile,"w");

strcpy(pathfile,root);
strcat(pathfile,"Stimuli.dat");
Stimuli = fopen(pathfile,"w");

/*********open the input data files and read inputs*******************/

strcpy(pathfile,root);
strcat(pathfile, "Mcell.dat");
Mdatafile = fopen(pathfile,"r");

strcpy(pathfile,root);
strcat(pathfile, "Pcell.dat");
Pdatafile = fopen(pathfile,"r");
    printf("reading the input files...");
	for (j=0;j<1001;j++){
		for(i=0;i<numcell;i++)
		{
		fscanf(Pdatafile,"%e",&StimP[j][i]); // row-major in C
		}
		}
	for (j=0;j<1001;j++) {
		for(i=0;i<numcell;i++)
		{
		fscanf(Mdatafile,"%e",&StimM[j][i]);
		}
		}
    printf(" reading completed.\n");
	fclose(Mdatafile);
	fclose(Pdatafile);


	/*** display and save stimuli ***/
	stimuluspoint=(int)(P1_NOUT);
	for (k=0;k<stimuluspoint;k++)
	{
	for (j=0; j<numxon;j++)
	{
	  fprintf(Stimuli, "%e\n", inputS(j,P1_DTOUT*k));
	  fflush(Stimuli);
	}
	}


  y = N_VNew_Serial(P1_NEQ);


  if(check_retval((void *)y, "N_VNew_Serial", 0)) return(1);


  PrintIntro1();

  cvode_mem = CVodeCreate(CV_BDF);

  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  miter=DIAG;
  ero = ZERO;
  //NV_Ith_S(y,0) = TWO;
  //NV_Ith_S(y,1) = ZERO;

  /*****************************
	initial values
  ******************************/
    printf("%i", Pch);

  for (i=0; i<Pch; i++)
    	NV_Ith_S(y,i) = 0.0;

  for (i=Pch; i<Mch; i++)
        NV_Ith_S(y,i) = 0.0;

    printf("came in!\n");

  /* initialize CVode */
  retval = CVodeInit(cvode_mem, f1, P1_T0, y);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  /* set scalar tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);


  retval = PrepareNextRun(cvode_mem, CV_BDF, miter, y, &A, 0, 0, &LS, &NLS);
  if(check_retval(&retval, "PrepareNextRun", 1)) return(1);

    PrintHeader1();

    for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      check_retval(&retval, "CVode", 1);
      temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
      if(check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
      temp_retval = CVodeGetLastStep(cvode_mem, &hu);
      if(check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
      PrintOutput1(t, NV_Ith_S(y,10), NV_Ith_S(y,250), qu, hu);
      if (retval != CV_SUCCESS) {
        nerr++;
        break;
      }


      /***************write to the files************************/

    for (j=0; j<Pch; j++)
      fprintf(cortexPdat,"%e\n", NV_Ith_S(y, j));
    for (j=Pch; j<Mch; j++)
      fprintf(cortexMdat,"%e\n", NV_Ith_S(y, j));

    fflush(cortexPdat);
    fflush(cortexMdat);

     if (iout%2 == 0) {
        er = fabs(NV_Ith_S(y,0)) / abstol;
        if (er > ero) ero = er;
        if (er > P1_TOL_FACTOR) {
          nerr++;
          //PrintErrOutput(P1_TOL_FACTOR);
        }
      }
    }

  /***************write to the files***********************

    for (j=0; j<numxon; j++){
    //realtype y1
    fprintf(Zon,"%10.5f  %e\n", t, NV_Ith_S(y, 250));
      //y1 = NV_Ith_S(y, j);
       // N_VPrintFile_Serial(y1, Zon);
}


    // N_VPrintFile_Serial(y, Zon);*/



  PrintFinalStats(cvode_mem, miter, ero);
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  N_VDestroy(y);

  return(nerr);
}

static void PrintIntro1(void)
{
  printf("Demonstration program for CVODE package - direct linear solvers\n");
  printf("\n\n");
  printf("Problem 1: Van der Pol oscillator\n");
  printf(" xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" neq = %d,  reltol = %.2Lg,  abstol = %.2Lg",
	 P1_NEQ, RTOL, ATOL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" neq = %d,  reltol = %.2g,  abstol = %.2g",
	 P1_NEQ, RTOL, ATOL);
#else
  printf(" neq = %d,  reltol = %.2g,  abstol = %.2g",
	 P1_NEQ, RTOL, ATOL);
#endif
}

static void PrintHeader1(void)
{
  printf("\n     t           x              xdot         qu     hu \n");

  return;
}

static void PrintOutput1(realtype t, realtype y0, realtype y1, int qu, realtype hu)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.5Lf    %12.5Le   %12.5Le   %2d    %6.4Le\n", t, y0, y1, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#else
  printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#endif

  return;
}

static int f1(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int i, k;
  double esum, isum, fbesum, fbisum;
  double Attention(), inputS(), inputT(), fb();
  //printf("came in!");


/******************** Cortex P cells ********************/

    for (i=0; i<=numcell-57; i++) { //81
    esum = 0; isum=0; fbesum=0, fbisum=0;
      for (k=0; k<57; k++){
        //esum=esum+(gc[57-k-1]*inputS(i+k, t) + 1*Attention(i+k,t)*0);  //gc[57-k-1];
        esum=esum+(1*inputS(i, t) + inputS(i, t)*Attention(i+k,t)*0); //gs[57-k-1]
        isum=isum+(gs[57-k-1]*inputS(i+k,t) + 1*(NV_Ith_S(y,i+k+numcell)));
        fbisum = fbisum + 0.75*fb(NV_Ith_S(y,i+0));  /*calculate the inhibitory feedback*/

         }
    /*calculate the excitatory feedback*/
     //fbisum = fbisum + 1*fb(NV_Ith_S(y,i));
     fbesum = fbesum + 1*fb(NV_Ith_S(y,i));
     NV_Ith_S(ydot,i)=1*(-(Ap+esum+isum+fbesum+fbisum)*
     NV_Ith_S(y,i)+Bp*(esum+fbesum)-0.0*(isum)); //Zcell_y
}

/******************** Cortex M cells ********************/

    for (i=0; i<=numcell-81; i++) {
    esum = 0; isum=0; fbesum=0;
      for (k=0; k<81; k++){
        esum=esum+(0.6*inputT(i, t)); // 0.5 for the max value
        isum=isum+(gys[81-k-1]*inputT(i+k, t) + 9*NV_Ith_S(y,i+0)); // 10 for the max value
        //fbesum = fbesum + 1*fb(NV_Ith_S(y,i+numcell+k));
         }
         //fbesum = fbesum + 5*fb(NV_Ith_S(y,i+numcell));
          NV_Ith_S(ydot,i+numcell)=1*(-(At+esum+isum+fbesum)*
     NV_Ith_S(y,i+numcell)+Bt*(esum+fbesum)-0.0*(isum));
}

  return(0);
}



static int PrepareNextRun(void *cvode_mem, int lmm, int miter, N_Vector y,
                          SUNMatrix* A, sunindextype mu, sunindextype ml,
                          SUNLinearSolver* LS, SUNNonlinearSolver* NLS)
{
  int retval = CV_SUCCESS;

  if (*NLS)
    SUNNonlinSolFree(*NLS);
  if (*LS)
    SUNLinSolFree(*LS);
  if (*A)
    SUNMatDestroy(*A);

  printf("\n\n-------------------------------------------------------------");

  printf("\n\nLinear Multistep Method : ");
  if (lmm == CV_ADAMS) {
    printf("ADAMS\n");
  } else {
    printf("BDF\n");
  }

  printf("Iteration               : ");
  if (miter == FUNC) {
    printf("FIXEDPOINT\n");

    /* create fixed point nonlinear solver object */
    *NLS = SUNNonlinSol_FixedPoint(y, 0);
    if(check_retval((void *)*NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

  } else {
    printf("NEWTON\n");

    /* create Newton nonlinear solver object */
    *NLS = SUNNonlinSol_Newton(y);
    if(check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

    printf("Linear Solver           : ");

    switch(miter) {

    case DENSE_USER :
      printf("Dense, User-Supplied Jacobian\n");

      /* Create dense SUNMatrix for use in linear solves */
      *A = SUNDenseMatrix(P1_NEQ, P1_NEQ);
      if(check_retval((void *)*A, "SUNDenseMatrix", 0)) return(1);

      /* Create dense SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Dense(y, *A);
      if(check_retval((void *)*LS, "SUNLinSol_Dense", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Set the user-supplied Jacobian routine Jac */
      //retval = CVodeSetJacFn(cvode_mem, Jac1);
      //if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;

    case DENSE_DQ :
      printf("Dense, Difference Quotient Jacobian\n");

      /* Create dense SUNMatrix for use in linear solves */
      *A = SUNDenseMatrix(P1_NEQ, P1_NEQ);
      if(check_retval((void *)*A, "SUNDenseMatrix", 0)) return(1);

      /* Create dense SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Dense(y, *A);
      if(check_retval((void *)*LS, "SUNLinSol_Dense", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Use a difference quotient Jacobian */
      retval = CVodeSetJacFn(cvode_mem, NULL);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;

    case DIAG :
      printf("Diagonal Jacobian\n");

      /* Call CVDiag to create/attach the CVODE-specific diagonal solver */
      retval = CVDiag(cvode_mem);
      if(check_retval(&retval, "CVDiag", 1)) return(1);
      break;

    case BAND_USER :
      printf("Band, User-Supplied Jacobian\n");

      /* Create band SUNMatrix for use in linear solves */
      *A = SUNBandMatrix(P2_NEQ, mu, ml);
      if(check_retval((void *)*A, "SUNBandMatrix", 0)) return(1);

      /* Create banded SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Band(y, *A);
      if(check_retval((void *)*LS, "SUNLinSol_Band", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);


    case BAND_DQ  :
      printf("Band, Difference Quotient Jacobian\n");

      /* Create band SUNMatrix for use in linear solves */
      *A = SUNBandMatrix(P2_NEQ, mu, ml);
      if(check_retval((void *)*A, "SUNBandMatrix", 0)) return(1);

      /* Create banded SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Band(y, *A);
      if(check_retval((void *)*LS, "SUNLinSol_Band", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Use a difference quotient Jacobian */
      retval = CVodeSetJacFn(cvode_mem, NULL);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;
    }
  }

  return(retval);
}

static void PrintErrOutput(realtype tol_factor)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\n\n Error exceeds %Lg * tolerance \n\n", tol_factor);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#else
  printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#endif

  return;
}

static void PrintFinalStats(void *cvode_mem, int miter, realtype ero)
{
  long int lenrw, leniw, lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf, nje, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  printf("\n Final statistics for this run:\n\n");
  printf(" CVode real workspace length              = %4ld \n",  lenrw);
  printf(" CVode integer workspace length           = %4ld \n",  leniw);
  printf(" Number of steps                          = %4ld \n",  nst);
  printf(" Number of f-s                            = %4ld \n",  nfe);
  printf(" Number of setups                         = %4ld \n",  nsetups);
  printf(" Number of nonlinear iterations           = %4ld \n",  nni);
  printf(" Number of nonlinear convergence failures = %4ld \n",  ncfn);
  printf(" Number of error test failures            = %4ld \n\n",netf);

  if (miter != FUNC) {
    if (miter != DIAG) {
      retval = CVodeGetNumJacEvals(cvode_mem, &nje);
      check_retval(&retval, "CVodeGetNumJacEvals", 1);
      retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);
      retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVodeGetLinWorkSpace", 1);
    } else {
      nje = nsetups;
      retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVDiagGetNumRhsEvals", 1);
      retval = CVDiagGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVDiagGetWorkSpace", 1);
    }
    printf(" Linear solver real workspace length      = %4ld \n", lenrwLS);
    printf(" Linear solver integer workspace length   = %4ld \n", leniwLS);
    printf(" Number of Jacobian evaluations           = %4ld \n", nje);
    printf(" Number of f evals. in linear solver      = %4ld \n\n", nfeLS);
  }

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" Error overrun = %.3Lf \n", ero);
#else
  printf(" Error overrun = %.3f \n", ero);
#endif
}

static void PrintErrInfo(int nerr)
{
  printf("\n\n-------------------------------------------------------------");
  printf("\n-------------------------------------------------------------");
  printf("\n\n Number of errors encountered = %d \n", nerr);

  return;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

/********************************** Functions ***********************************/


void newgauss(garray,Amp,var,halfsize)
double garray[], Amp, var;
int halfsize;
{
	int i, size;
	size=halfsize*2 + 1;
	for (i=0; i<=size; i++)
   	{
	 garray[i]=Amp*exp(-(i-halfsize)*(i-halfsize)/(var*var));
	}
}

double Attention(i, t)
int i;
realtype t;
{
    if (i > 225 && i < 325 && t > 1 && t < 250)
    {
        return(0.61); //0.61 for the temporal task // if you increase the pulse amplitude to a certain thershold the output will be negative due to the concalib function
        }            // 0.3 for the spatial task
    else
        return(0.0);
}


double inputS(i,t)
int i;
double t;
{
	double temp, fract;
	int timeindex;
	if(t < 0) return(0.0);
	else
	{
	temp = t / rtimeinterval;
	timeindex = temp;
	fract = temp - timeindex;
	return(retcorgain*(StimP[timeindex][i]
		+fract*(StimP[timeindex+1][i]-StimP[timeindex][i])));
	}
}


double inputT(i,t)
int i;
double t;
{
	double temp, fract;
	int timeindex;
	if(t < 0) return(0.0);
	else
	{
	temp = t / rtimeinterval;
	timeindex = temp;
	fract = temp - timeindex;
	return(StimM[timeindex][i]
		+fract*(StimM[timeindex+1][i]-StimM[timeindex][i]));
	}
}

double fb(a)
double a;
{
	if(a < 0.05)
		{
		return(10.0*a*(a+1)*(a+1)-10.0*a); //10.0*a*(a+1)*(a+1)-10.0*a
		}
	return(a*a+0.975*a); //a*a+0.975*a
}
/*
double concalib(a)
double a;
{
    //return (-9.26062 + 3.475*a - 0.164986*a*a);
    return(a);
}

double Xthresh(a)
double a;
{
	if(a < 231.64)
		{
		return(231.64);
		}
	return(a);
}

double Ythresh(a)
double a;
{
	if(a < 586.7856)
		{
		return(586.7856);
		}
	return(a);
}
*/
