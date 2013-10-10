/******************************************************************************/
/*                                                                            */
/*                                definitions                                 */
/*                                                                            */
/******************************************************************************/

#define ERR_VAR       265
#define CROSS_FCN_W   266
int RES_NORM_IN_FCN=L2_NORM;      /* L1_NORM, L2_NORM, LP_NORM, RES_FCN,   */
double P_IN_LP_NORM=2.;           /*                          RES_X_CROSS  */
double RES_SCALING=0.1;
double MAX_CROSS_DER=1.;
double CROSS_DER_WEIGHT=0.;
double RES_NORM_WEIGHT=1.;
#define TAU_SPACE     P0           /*  P0, P1_NC, P1C, ...  */
int TAU_VARIABLE=1;                /*  in which index is tau saved, changeable
                                                     during the program run  */
#define A_IN_FCN 1.
#define SCALED_IND    NO           /*  YES or NO  */
#define PAR_TYPE      Q_VE
#define USE_GRAD_IND  NO
#define ELEMENT_TYPE  SIMPLEX      /*  SIMPLEX or CUBE                        */
#define SC_EXAMPLE    9            /*  8, 81, 9, 20, 55, 66, 67, 800  */
#define COARSE_GRID   cube         /*  cube or step2_xy  */
//#define COARSE_GRID   step2_xy      /*  cube or step2_xy  */
#define REFINEMENTS   0
#define U_SPACE       P2C          /*  P1_NC, P1_MOD, P1C, P2C  */
#define SDFEM         YES          /*  YES or NO  */
#define DIR_DIAM      YES          /*  YES or NO .. diameter influenced by b  */
#define USE_QUADRATURE NO          /*  YES or NO  */
#define PW_CONST_PAR   YES
#define LOC_PROJECT   NO           /*  YES or NO  */
#define CONV_LOC_PROJECT   NO      /*  YES or NO  */
#define USE_BM        YES          /*  YES or NO  */
#define USE_BEL       NO           /*  YES or NO  */
#define METHOD        UMFPACK      /*  BASIC_GMRES, DIAG_GMRES, ILU_GMRES,    */
                                   /*  MULTI_GRID or UMFPACK                  */
#define TAU0 25.
#define DAMPING_PARAM 1.
#define SOLVE_BY_IMPL_EULER  NO
#define IE_TAU 0.001

int COMPUTE_SC=NO;
#define SC_BETA     0.4
#define SC_BETA2    0.
#define SC_TYPE    SCM_JSW
#define REG_FABS   NO              /*  YES or NO  */   /* for K4 */
#define REG_EPS    1.
#define DELTA_TYPE    OUTFLOW_D    /*  BASIC_D1, BASIC_D2, CLASSIC_D or
                                                                    OUTFLOW_D */
int RESTRICT_FCN=NO;
int RESTRICT_DER=NO;
#define FCN_AT_INFLOW NO
#define PAR_OPT       YES

/*----------------------------------------------------------------------------*/

#define DIM           2            /*  space dimension                        */
#define SOLVED_PDE    CONV_DIFF
#define USE_COARSE_F  NO           /*  YES or NO  -  coarse grid from a file  */
#define COARSE_FILE   "tr2360"     /*  file containing the coarse grid        */
#define NS            142      
#define REF_TYPE      UNIFORM      /*  UNIFORM, BISECTION, RED_GREEN          */
#define U_STRUCTURE   SCALAR       /*  SCALAR or VECTOR  */
#define P_SPACE       NONE
#define RHS_INTEGR    CUBIC        /*  LINEAR  or  CUBIC  */
#define CONV_DISCR    CONV         /*  CONV  or  SKEW  */
#define DATA_STR      0
#define PROBLEM       1
#define INITIAL       ZERO         /*  ZERO  or  SOLUTION  */
#define DELTA_FACTOR  0.25         /*  delta = hT*DELTA_FACTOR  */
#define SAVE_TRIANG   NO           /*  YES or  NO  */
#define ERR_ON_SQUARE NO           /*  YES or  NO  */
#define X_SQUARE      0.801
#define Y_SQUARE      0.801
#define MOVING_BOUNDARY NO
#include "sc_examples.h"
double  TNU=VALUE_OF_EPS;
double  TAU_K;
double AUX1;
INT TOPLEVEL;
//#include "bc-cd/bc1.c"
//#include "bc-cd/bc888.h"
#define T_FOR_U         121
#if METHOD == UMFPACK
#define INCLUDE_UMFPACK  YES
#endif

/*----------------------------------------------------------------------------*/

#define NMASK           2
#define PNMASK          4
#define OUT_MASK        16
#define FMASK           1024

#define IS_DIR(n,t)     ( t != 0 && (NTYPE(n) & NMASK))  /*  for t=0, we have  \
                         IS_DIR = 0, otherwise IS_DIR = NTYPE(n) & NMASK.      \
                         t=0 is used for the MINI element                     */
/* #define IS_DIR(n,t)     ( t==1 && ( NTYPE(n) & GW || NTYPE(n) & GLS ) ||    \
                          t==2 && ( NTYPE(n) & GW || NTYPE(n) & GLG ) ) */

#define BSC             1.0E-10
#define EPS             1.0E-18
#define EPSA            1.0E-12
#define EPSC            1.0E-10
#define EPSAL           1.0E-8
#define EPST            1.0E-13
#define EPSE            1.0E-4
#define EPSUP           0.0
#define EPS_CGAR        5.0E-2
#define EPS_GMRES       1.0E-13

int     NV=             9;      /* number of vertices in one direction        */
int     NVX=            9;
int     NVY=            9;
int     TNF;                    /* total number of inner faces                */

#define V_MATR          24 
#define F_MATR          24
#define B_MATR          24
#define S_MATR          24
#define D_MATR          24
#define NS_VECT         272
#define NV_VECT         24
#define NM_VECT         24
#define FS_VECT         272
#define FV_VECT         24
#define FM_VECT         24
#define FDV_VECT        24
#define SNS_VECT        24
#define SFS_VECT        24
#define LG_VECT         24
#define EL_VECT         272
#define EL_VECTD        272
#define EL_VECTN        272
#define EM_VECT         24

double NU  = 1.;
double FMULT3  = 1.;
double ALPHA1 = 0.;
double ALPHA3 = 0.;
double KKK = 1.0;
double MAX_SIGMA = 4.;

#define FMULT           1.        /* multiplier for faces                     */
#define FMULT2          1.        /* FMULT*FMULT                              */
#define GMN             20        /* for GMRES                                */
#define N_GGEM          100
#define N_SGGEM         100
