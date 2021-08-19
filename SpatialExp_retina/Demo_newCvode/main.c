/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Demonstration program for CVODE - direct linear solvers.
 * Two separate problems are solved using both the CV_ADAMS and CV_BDF
 * linear multistep methods in combination with the
 * SUNNONLINSOL_FIXEDPOINT and SUNNONLINSOL_NEWTON nonlinear solver
 * modules:
 *
 * Problem 1: Van der Pol oscillator
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) dense, user-supplied, (2) dense, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * Problem 2: ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) band, user-supplied, (2) band, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * For each problem, in the series of eight runs, CVodeInit is
 * called only once, for the first run, whereas CVodeReInit is
 * called for each of the remaining seven runs.
 *
 * Notes: This program demonstrates the usage of the sequential
 * macros NV_Ith_S, SM_ELEMENT_D, SM_COLUMN_B, and
 * SM_COLUMN_ELEMENT_B. The NV_Ith_S macro is used to reference the
 * components of an N_Vector. It works for any size N=NEQ, but
 * due to efficiency concerns it should only by used when the
 * problem size is small. The Problem 1 right hand side and
 * Jacobian functions f1 and Jac1 both use NV_Ith_S. The
 * N_VGetArrayPointer function gives the user access to the
 * memory used for the component storage of an N_Vector. In the
 * sequential case, the user may assume that this is one contiguous
 * array of reals. The N_VGetArrayPointer function
 * gives a more efficient means (than the NV_Ith_S macro) to
 * access the components of an N_Vector and should be used when the
 * problem size is large. The Problem 2 right hand side function f2
 * uses the N_VGetArrayPointer function. The SM_ELEMENT_D macro
 * used in Jac1 gives access to an element of a dense SUNMatrix. It
 * should be used only when the problem size is small (the
 * size of a Dense SUNMatrix is NEQ x NEQ) due to efficiency concerns. For
 * larger problem sizes, the macro SM_COLUMN_D can be used in order
 * to work directly with a column of a Dense SUNMatrix. The SM_COLUMN_B and
 * SM_COLUMN_ELEMENT_B allow efficient columnwise access to the elements
 * of a Banded SUNMatix. These macros are used in the Jac2 function.
 * -----------------------------------------------------------------*/

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
#include "retina.h"
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

#define P1_NEQ        2900
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
#define retthresh   249.624 //247.96
#define persistence 0.105


/********** Global variable ***************/
//int sxzon,sxzons,sxon,sxonper,sxzoff,sxoff,sy,scell, syzon, syzons;
int xon, sxon, yon, Zcellx, Zcelly;
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

int L, R;
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
  double SumGe, SumGi, SumGey, SumGiy;
  double Impulse(), concalib(), Xthresh(), Ythresh(), retcor();
  void newgauss();
  FILE *fopen(),*stimuli,*xondat, *Zon, *Xcells, *SXcells, *Ycells;

  y = NULL;
  A = NULL;
  LS = NULL;
  NLS = NULL;
  cvode_mem = NULL;

  /*File names to store the data*/
  int gap = 6;
  double inc = 2.9;
  if (gap == 6) {
    L = 290 - (4 + 5.5*inc); //25
    R = 290 + (4 + 5.5*inc);
  }
  else if (gap == 5) {
    L = 290 - (4 + 4.5*inc); //20
    R = 290 + (4 + 4.5*inc);
  }
  else if (gap == 4) {
    L = 290 - (4 + 3.5*inc); //15
    R = 290 + (4 + 3.5*inc);
  }
  else if (gap == 3) {
    L = 290 - (4 + 2*inc); //11
    R = 290 + (4 + 2*inc);
  }
  else if (gap == 2) {
    L = 290 - (4 + inc); //7
    R = 290 + (4 + inc);
  }
  else {
    L = 290 - 4;
    R = 290 + 4;
  }

  //char root[100] = "/home/coglab/Documents/Results/Retina_p3/positiveStep/gap6/Neutral/";
  char root[100] = "/home/coglab/Documents/Results/Retina_p3/positiveStep/gap6/Attention/"; // the original values are in the folder positiveStep1

  char pathfile[100];
  //char pathfileint[100] = "./data";

  /***********INITIAL CONDITIONS ******/
  Zcellx = numxon;
  Zcelly = numxon + numxon;
  xon = numxon + numxon + numxon;
  sxon = numxon + numxon + numxon + numxon;
  yon = numxon + numxon + numxon + numxon + numxon;
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

  Zcellxeq =alpha*beta/(alpha+gama*(Jval+L0));
  Zcellyeq =alphay*betay/(alphay+gamay*(Jvaly+L0));

  Xcelleq =((Bxon*SumGe-Dxon*SumGi)*concalib((Jvalxon+L0))*Zcellxeq)/(Axon + (SumGe + SumGi)*concalib((Jvalxon+L0))*Zcellxeq);
  Ycelleq =((By*SumGey-Dxon*SumGiy)*concalib((Jvalxon+L0))*Zcellyeq)/(Axon + (SumGey + SumGiy)*concalib((Jvalxon+L0))*Zcellyeq);


  // Print intitial condition values
  printf("Zcells = %e\n", Zcellxeq);
  printf("Xcells = %e\n", Xcelleq);
  printf("Xcells = %e\n", alpha);


  strcpy(pathfile, root);
  strcat(pathfile, "Simulus.dat");
  stimuli  = fopen(pathfile,"w");


  strcpy(pathfile, root);
  strcat(pathfile, "Zcells.dat");
  Zon = fopen(pathfile,"w");

  strcpy(pathfile, root);
  strcat(pathfile, "Xcells.dat");
  Xcells = fopen(pathfile,"w");

  strcpy(pathfile, root);
  strcat(pathfile, "SXcells.dat");
  SXcells = fopen(pathfile,"w");

  strcpy(pathfile, root);
  strcat(pathfile, "Ycells.dat");
  Ycells = fopen(pathfile,"w");


  /*** display and save stimuli ***/
	stimuluspoint=(int)(P1_NOUT);
	for (k=0;k<stimuluspoint;k++)
	{
	for (j=0; j<numxon;j++)
	{
	  fprintf(stimuli, "%e\n", Impulse(j,P1_DTOUT*k));
	  fflush(stimuli);
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

  for (i=0; i<Zcellx; i++)
    	NV_Ith_S(y,i) = Zcellxeq;

  for (i=Zcellx; i<Zcelly; i++)
        NV_Ith_S(y,i) = Zcellyeq;

  for (i=Zcelly; i<xon; i++)
    	NV_Ith_S(y,i) = Xcelleq;

  for (i=xon; i<sxon; i++)
    	NV_Ith_S(y,i) = retcor(Xcelleq);

  for (i=sxon; i<yon; i++)
    	NV_Ith_S(y,i) = Ycelleq;


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

    for (j=0; j<numxon; j++)
      fprintf(Zon,"%e\n", NV_Ith_S(y, j));
    for (j=Zcelly; j<xon; j++)
      fprintf(Xcells,"%e\n", NV_Ith_S(y, j));
    for (j=xon; j<sxon; j++)
      fprintf(SXcells,"%e\n", NV_Ith_S(y, j));
    for (j=sxon; j<yon; j++)
      fprintf(Ycells,"%e\n", Ythresh(NV_Ith_S(y, j) ));



    fflush(Zon);
    fflush(Xcells);
    fflush(SXcells);
    fflush(Ycells);

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
  double esum, isum;
  double Impulse(), concalib(), retcor();
  realtype tau;

  tau = t - xydelay;

  /****transmitter eqns for the center X cell*****/
	for (i=0; i<Zcellx; i++)
	  NV_Ith_S(ydot,i)=timescale*(alpha*(beta-NV_Ith_S(y,i))
                             -gama*NV_Ith_S(y,i)*(Jval +Impulse(i,tau)));
  /****transmitter eqns for the center Y cell*****/
    for (i=0; i<Zcellx; i++)
	  NV_Ith_S(ydot,i+Zcellx)=timescale*(alphay*(betay-NV_Ith_S(y,i+Zcellx))
                             -gamay*NV_Ith_S(y,i+Zcellx)*(Jvaly +Impulse(i,t)));

/****Sustain eqns for the X cell*****/

    for (i=0; i<=Zcellx-57; i++) { //81
    esum = 0; isum=0;
      for (k=0; k<57; k++){
    esum=esum+gc[57-k-1]*concalib(Jvalxon + Impulse(i+k,tau))
         *NV_Ith_S(y,i+k);
        isum=isum+gs[57-k-1]*concalib(Jvalxon + Impulse(i+k,tau))
         *NV_Ith_S(y,i+k);
         }
          NV_Ith_S(ydot,i+Zcelly)=timescale*(-(Axon+esum+isum)*
     NV_Ith_S(y,i+Zcelly)+Bxon*esum-Dxon*(isum)); //Zcell_y
}


/*************************** Xon per equation ************************/

	for (i=0; i<Zcellx-57; i++)
	NV_Ith_S(ydot,i+xon)=persistence*(-NV_Ith_S(y,i+xon)+retcor(NV_Ith_S(y,i+Zcelly)));


/****Transient eqns for the Y cell*****/

    for (i=0; i<=numxon-81; i++) {
    esum = 0; isum=0;
      for (k=0; k<81; k++){
    esum=esum+gyc[81-k-1]*concalib(Jvalxon + Impulse(i+k,t))
         *NV_Ith_S(y,i+k+Zcellx);
        isum=isum+gys[81-k-1]*concalib(Jvalxon + Impulse(i+k,t))
         *NV_Ith_S(y,i+k+Zcellx);
         }
          NV_Ith_S(ydot,i+sxon)=timescale*(-(Axon+esum+isum)*
     NV_Ith_S(y,i+sxon)+By*esum-Dxon*(isum));
}

  return(0);
}

/* static int Jac1(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y0, y1;

  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  SM_ELEMENT_D(J,0,1) = ONE;
  SM_ELEMENT_D(J,1,0) = -TWO * P1_ETA * y0 * y1 - ONE;
  SM_ELEMENT_D(J,1,1) = P1_ETA * (ONE - SQR(y0));

  return(0);
}*/



/*static realtype MaxError(N_Vector y, realtype t)
{
  sunindextype i, j, k;
  realtype *ydata, er, ex=ZERO, yt, maxError=ZERO, ifact_inv, jfact_inv=ONE;

  if (t == ZERO) return(ZERO);

  ydata = N_VGetArrayPointer(y);
  if (t <= THIRTY) ex = exp(-TWO*t);

  for (j = 0; j < P2_MESHY; j++) {
    ifact_inv = ONE;
    for (i = 0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      yt = pow(t, i+j) * ex * ifact_inv * jfact_inv;
      er = fabs(ydata[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}*/

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

double Impulse(i, t)
int i;
realtype t;
{
    if (i > L && i < R && t > 1 && t < 50)
    {
        return((20)*Intensity + 1*L0);  // if you increase the pulse amplitude to a certain thershold the output will be negative due to the concalib function
        //return(0.0);
        }
    else
        return(L0);
}

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

double concalib(a)
double a;
{
    //return (-9.26062 + 3.475*a - 0.164986*a*a);
    return(a);
}

double Xthresh(a)
double a;
{
	if(a < 231.64) //231.64
		{
		return(231.64);
		}
	return(a);
}

double Ythresh(a)
double a;
{
	if(a < 587.4) // set value 587; except for plot use 587.4
		{
		return(0.0);
		}
	return(a-587.4); // a-587
}

double retcor(a)
double a;
{
	if(a < retthresh)
		{
		return(0.0); // 0.0
		}
/*      return(0.00125*50.0*pow(a-retthresh,2.0));   */
        return(1*pow(a-1*retthresh,1));
}

