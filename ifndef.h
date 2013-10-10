/******************************************************************************/
/*                                                                            */
/*                                 macros 2                                   */
/*                                                                            */
/******************************************************************************/

#ifndef T_FOR_P
#define T_FOR_P         0
#endif

#ifndef A_LOW
#define A_LOW           1
#endif

#ifndef B_LOW
#define B_LOW           1
#endif

#ifndef LOW_AND_HIGH
#define LOW_AND_HIGH    NO
#endif

#ifndef U_SPACE_LOW
#define U_SPACE_LOW     -1
#endif

#ifndef P_SPACE_LOW
#define P_SPACE_LOW     -1
#endif

#ifndef U_TYPE_LOW
#define U_TYPE_LOW      -1
#endif

#ifndef P_TYPE_LOW
#define P_TYPE_LOW      -1
#endif

#ifndef C_STRUCT    
#define C_STRUCT         0
#endif

#ifndef A_STRUCT_LOW
#define A_STRUCT_LOW    -1
#endif

#ifndef B_STRUCT_LOW
#define B_STRUCT_LOW    -1
#endif

#ifndef C_STRUCT_LOW
#define C_STRUCT_LOW     0
#endif

#ifndef ILU_TYPE_LOW
#define ILU_TYPE_LOW    -1
#endif

#ifndef ILU_STRUCT_LOW
#define ILU_STRUCT_LOW  -1
#endif

#ifndef USE_DIV_STAB
#define USE_DIV_STAB     NO
#endif

#ifndef USE_DIV_STAB_LOW
#define USE_DIV_STAB_LOW NO
#endif

#ifndef KORN_LAPLACE  
#define KORN_LAPLACE    NO
#endif

#ifndef REF_TYPE
#define REF_TYPE        UNIFORM
#endif

#ifndef F_NAME
#define F_NAME          "Aug"
#endif

#ifndef DIR_FOR_RESULTS
#define DIR_FOR_RESULTS   "VYSLEDKY/NU_1/"
#endif

#ifndef PROFILE_LEVEL
#define PROFILE_LEVEL 0.5
#endif

#ifndef DIRCRITERION
#define DIRCRITERION(x,y,type)  (type && fabs(y)<EPS)
#endif

#ifndef ZNCRITERION
#define ZNCRITERION(x,y,type)   (type && (fabs(x)<EPS || fabs(1.-x)<EPS))
#endif

#ifndef FBCRITERION
#define FBCRITERION(x,y,type)   (type && y > PROFILE_LEVEL - EPS)
#endif

#ifndef CC0
#define CC0   0.
#endif

#ifndef DELTA_TYPE  
#define DELTA_TYPE      BASIC_D1
#endif

#ifndef SDFEM
#define SDFEM           NO
#endif

#ifndef CONV_LOC_PROJECT
#define CONV_LOC_PROJECT NO
#endif

#ifndef LOC_PROJECT
#define LOC_PROJECT     NO
#endif

#ifndef STABILIZATION
#define STABILIZATION   NONE
#endif

#ifndef DELTA_FACTOR
#define DELTA_FACTOR    1.
#endif

#ifndef PECLET_FACTOR
#define PECLET_FACTOR  1.
#endif

#ifndef DAMPING_PARAM
#define DAMPING_PARAM  1.
#endif

#ifndef MAKE_AVERAGES
#define MAKE_AVERAGES  NO
#endif

#ifndef INCLUDE_UMFPACK
#define INCLUDE_UMFPACK  NO
#endif

#ifndef Q1ROT_T
#define Q1ROT_T(x)   ((x)*(x))
#endif

#ifndef Q1ROT_TD
#define Q1ROT_TD(x)  (2.*(x))
#endif

#ifndef Q1ROT_K
#define Q1ROT_K      4.
#endif

#ifndef Q1ROT_E
#define Q1ROT_E      (1./3.)
#endif

#ifndef DOF
#define DOF          MIDPOINT
#endif

#ifndef RHS_INTEGR
#define RHS_INTEGR      LINEAR
#endif

#ifndef METHOD
#define METHOD      BASIC_CG
#endif

#ifndef INITIAL
#define INITIAL     ZERO
#endif

#ifndef NMASK_FOR_SF
#define NMASK_FOR_SF  4096
#endif

#ifndef NMASK_FOR_FB
#define NMASK_FOR_FB     4
#endif

#ifndef NMASK_ZN
#define NMASK_ZN         8
#endif

#ifndef FMASK_ZN
#define FMASK_ZN       512
#endif

#ifndef SNMASK
#define SNMASK      4096
#endif

#ifndef SFMASK
#define SFMASK      4096
#endif

#ifndef SURFACE
#define SURFACE        1.
#endif

#ifndef COORDS_FOR_FACES
#define COORDS_FOR_FACES   NO
#endif

#ifndef LOWER_LEFT_X
#define LOWER_LEFT_X  -0.5
#endif

#ifndef LOWER_LEFT_Y
#define LOWER_LEFT_Y  -0.5
#endif

#ifndef LENGTH_X
#define LENGTH_X  1.
#endif

#ifndef LENGTH_Y
#define LENGTH_Y  1.
#endif

#ifndef N_OF_NODE_FUNC
#define N_OF_NODE_FUNC  1
#endif

#ifndef N_OF_FACE_FUNC
#define N_OF_FACE_FUNC  1
#endif

#ifndef N_OF_ELEM_FUNC
#define N_OF_ELEM_FUNC  1
#endif

#ifndef NM_VECT
#define NM_VECT         1
#endif

#ifndef FM_VECT
#define FM_VECT         1
#endif

#ifndef EM_VECT
#define EM_VECT         1
#endif

#ifndef MOVING_BOUNDARY
#define MOVING_BOUNDARY NO
#endif

#ifndef SLIDING_POINTS
#define SLIDING_POINTS NO
#endif

#ifndef USE_QUADRATURE
#define USE_QUADRATURE NO
#endif

#ifndef PW_CONST_PAR
#define PW_CONST_PAR   YES
#endif

#ifndef USE_BM
#define USE_BM         NO
#endif

#ifndef USE_BEL
#define USE_BEL        YES
#endif

#ifndef SOLVE_BY_IMPL_EULER
#define SOLVE_BY_IMPL_EULER  NO
#endif

#ifndef TAU0
#define TAU0    1.
#endif

#ifndef IE_TAU
#define IE_TAU  1.
#endif

#ifndef MH_VER
#define MH_VER  1
#endif

#ifndef SC_EXAMPLE
#define SC_EXAMPLE 0
#endif

#ifndef COMPUTE_SC
#define COMPUTE_SC NO
#endif

#ifndef CUBE_TYPE
#define CUBE_TYPE 1
#endif

#ifndef DIR_DIAM
#define DIR_DIAM  NO
#endif

#ifndef SC_BETA
#define SC_BETA     0.
#endif

#ifndef SC_BETA2
#define SC_BETA2    0.
#endif

#ifndef SC_TYPE
#define SC_TYPE     0
#endif

#ifndef REG_FABS
#define REG_FABS   NO
#endif

#ifndef REG_EPS
#define REG_EPS    1.
#endif

#ifndef PAR_OPT
#define PAR_OPT    NO
#endif

#ifndef PNMASK
#define PNMASK          32
#endif

#ifndef PFMASK
#define PFMASK          2048
#endif

#ifndef PERIODIC_BC
#define PERIODIC_BC     NO
#endif

#ifndef DIRECTION_OF_EDGES
#define DIRECTION_OF_EDGES   YES
#endif

#ifndef ELEM
#define ELEM            q1c_element
#endif

#ifndef REF_MAP
#if ELEMENT_TYPE == SIMPLEX
#define REF_MAP      P1_REF_MAP
#elif ELEMENT_TYPE == CUBE
#define REF_MAP      Q1_REF_MAP
#else
#define REF_MAP      P1_REF_MAP
#endif
#endif

#ifndef L2_Q_RULE
#define L2_Q_RULE    NULL
#endif

#ifndef H10_Q_RULE
#define H10_Q_RULE   NULL
#endif

#ifndef LAPL_Q_RULE
#define LAPL_Q_RULE  NULL
#endif

#ifndef CONV_Q_RULE
#define CONV_Q_RULE  NULL
#endif

#ifndef REACT_Q_RULE
#define REACT_Q_RULE NULL
#endif

#ifndef SD_Q_RULE
#define SD_Q_RULE    NULL
#endif

#ifndef SC_Q_RULE
#define SC_Q_RULE    NULL
#endif

#ifndef RHS_Q_RULE
#define RHS_Q_RULE   NULL
#endif

#ifndef SDRHS_Q_RULE
#define SDRHS_Q_RULE NULL
#endif

#ifndef N_FOR_EVAL
#define N_FOR_EVAL   10001
#endif


