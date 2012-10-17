/******************************************************************************/
/*                                                                            */
/*                                 macros 1                                   */
/*                                                                            */
/******************************************************************************/

int NUMBER_OF_ERRORS=0;

#define Q29             2.0/9.0 

/* #define NULL 0 */

#define SHORT  short
#define INT    int
#define UNS    unsigned short
#define FLOAT  double
#define DOUBLE double

/*  for DATA_STR:
    -------------                                                             */
#define LG_DATA         2
#define KORN_MATRIX     16 /* (DATA_STR & KORN_MATRIX) == 0 ... velocity matrix 
               corresponds to the bilinear form from Friedrichs' inequality   */
#define BOUNDARY_APPROX 64
#define REDUCED         128
#define NEW_BUBBLES     256
#define BD_DATA         512

/*  for PROBLEM:
    ------------                                                              */
#define WITH_F       1  /* (PROBLEM & WITH_F) == 1 ... rhs for faces computed */
#define BUBBLE_BC    2  /* (PROBLEM & BUBBLE_BC) == 1 ... b.c. with bubbles   */

/*  for DATA_S:
    -----------                                                               */
#define N_LINK_TO_FACES   1
#define F_LINK_TO_FACES   2
#define F_LINK_TO_NODES   4
#define N_LINK_TO_NODES   8
#define PREVIOUS_NODE    16
#define PREVIOUS_FACE    32
#define SPECIAL_NODES_AND_FACES    64
#define PERIODIC_NODES     128
#define PERIODIC_FACES     256
#define N_LINK_TO_ELEMENTS 512

/*  for N_DATA:
    -----------                                                               */
#define ONE_NODE_MATR        1
#define IxD_NODE_MATR        2
#define DxD_NODE_MATR        4
#define SCALAR_NODE_DATA     8
#define VECTOR_NODE_DATA    16
#define Dx1_NODE_FACE_MATR  32
#define ONE_NODE_FACE_MATR  64
#define DxD_NODE_FACE_MATR 128
#define NODE_IFLAG         256
#define NUMBER_OF_N_NEL    512  /*  # neighbouring elements for each node  */
#define Ix2D_NODE_FACE_MATR  1024
#define NODE_BD_DATA         2048
#define NxN_NODE_MATR        4096
#define NxM_NODE_FACE_MATR   8192
#define MVECTOR_NODE_DATA   16384
#define NODE_ITYPE          32768

/*  for F_DATA:
    -----------                                                               */
#define ONE_FACE_MATR          1
#define IxD_FACE_MATR          2
#define Ix2D_FACE_MATR         4
#define DxD_FACE_MATR          8
#define SCALAR_FACE_DATA      16
#define VECTOR_FACE_DATA      32
#define DVECTOR_FACE_DATA     64
#define IxD_FACE_NODE_MATR   128
#define ONE_FACE_NODE_MATR   256
#define DxD_FACE_NODE_MATR   512
#define FACE_IFLAG          1024
#define FACE_MIDDLE         2048
#define CURVED_FACE_MIDDLE  4096
#define FACE_BD_DATA        8192
#define NxN_FACE_MATR      16384
#define MxN_FACE_NODE_MATR 32768
#define MVECTOR_FACE_DATA  65536

/*  for E_DATA:
    -----------                                                               */
#define SCALAR_ELEMENT_DATA            1
#define VECTOR_ELEMENT_DATA            2
#define SCALAR_DATA_IN_ELEMENT_NODES   4
#define ExDN_MATR                      8
#define ExF_MATR                      16
#define ExDF_MATR                     32
#define ExFx2DF_MATR                  64
#define E_E_NEIGHBOURS               128
#define FxE_MATR                     256
#define ExN_MATR                     512
#define NxE_MATR                    1024
#define ExE_MATR                    2048
#define ExFxDN_MATR                 4096
#define ExFxDF_MATR                 8192
#define MxM_E_E_MATR               16384
#define MxN_E_N_MATR               32768
#define NxM_N_E_MATR               65536
#define MxN_E_F_MATR              131072
#define NxM_F_E_MATR              262144
#define MVECTOR_ELEMENT_DATA      524288
#define E_E_FNEIGHBOURS          1048576
#define E_TAU                    2097152
#define ELEM_ITYPE               4194304
#define E_TAU_SOLD               8388608

/*  for T: 
    ------                                                                    */
#define ONLY_INNER                     1  /*  only inner nodes/faces          */
#define STOP_IS_FIRST_INNER            2  /*  nodes/faces before 1st inner    */
#define USE_IS_DIR                     4  /*  IS_DIR determines inner/bound   */
#define NSTART_FROM_INNER              8  /*  only inner neighbouring nodes   */
#define NFSTART_FROM_INNER            16  /*  only inner neighbouring faces   */
#define FSTART_FROM_INNER             32  /*  only inner neighbouring faces   */
#define FNSTART_FROM_INNER            64  /*  only inner neighbouring nodes   */
#define WITHOUT_FIRST                512
#define GO_THROUGH_ALL              1024
#define ZERO_MEAN                   2048
#define ADD_ZN_CONDITION            4096  /*  zero normal comp. of velocity   */

#define NONE            0
#define ALL        100000
#define SIMPLEX         1
#define CUBE            2
#define P1_MOD          1
#define P0              2
#define P1C             3
#define P1_DISC         4
#define P1C_FBUB        5
#define P1C_NEW_FBUB    6
#define P1C_FBUB_LG     7
#define P1C_ELBUB       8
#define P1_NC           9
#define P2C            10
#define P2C_ELBUB      11
#define Q1ROT          12
#define IP1C           13
#define IP2C           14
#define MINI_L_DIV_FR  15 /* div-free basis of MINI el. obtained from lin. fcn*/
#define IP2C_ELBUB     16
#define IP1_DISC       17
#define P1C_P2L        18 /* P2 Lagrange multipliers on some (boundary) edges */
#define P2L            19 /* P2 Lagrange multipliers on some (boundary) edges */
#define Q1C            20
#define Q2C            21
#define Q3C            22
#define Q4C            23
#define P3C            33
#define P4C            34
#define GP1C          503
#define GP1X3C        504
#define GP1C_ELBUB    508
#define GP2C          510
#define GP2X3C        511
#define GP2C_3ELBUB   518
#define GP2C_6ELBUB   519
#define GQ1C          520
#define GQ1X4C        521
#define GQ1C_ELBUB    525
#define GQ2C          540
#define GQ2X4C        541
#define GQ2B3C        542 /* Q3 on elements along the boundary */
#define GQ2C_2ELBUB   545
#define GQ2C_3ELBUB   546
#define GQ3C          550
#define GQ4C          560
#define GP3C          570
#define GP4C          580
#define P1_REF_MAP      1
#define P2_REF_MAP      2
#define Q1_REF_MAP      3
#define MIDPOINT        1
#define MEAN_VALUE      2
#define NO              0
#define YES             1
#define ZERO            0
#define SOLUTION        1
#define UNIFORM         1
#define REFINE_TO_MID   2
#define BISECTION       3
#define RED_GREEN       4
#define BASIC_GMRES     0
#define DIAG_GMRES      1
#define ILU_GMRES       2
#define BASIC_CG        3
#define BS_SMOOTHER     6
#define VANKA_SMOOTHER  7
#define MULTI_GRID      8
#define AUGM_LAGR       9
#define BS_SMOOTHER_EX 10
#define ILU_IT         11
#define SOR_F          12
#define UMFPACK        30
#define LINEAR          1
#define QUADRATIC       2
#define CUBIC           3
#define Q_FORMULA     100
#define CONV            0
#define SKEW            1
#define SCALAR          1
#define P_SCALAR        2
#define VECTOR          3
#define V_CYCLE         1
#define W_CYCLE         2
#define NONZERO_INIT    16384
#define IDENT_PR        1
#define DIAG_PR         2
#define MDIAG_PR        4
#define ILU_PR          8
#define SSOR_PR        16
#define UMF_PR         32
#define POISSON         1
#define CONV_DIFF       2
#define STOKES          3
#define OSEEN           4
#define NAVIER_STOKES   5
#define BASIC_D1        1
#define BASIC_D2        2
#define CLASSIC_D       3
#define SE1_D           4   /* Shih, Elman (1999) */
#define SE2_D           5   /* Shih, Elman (1999) */
#define OUTFLOW_D       6
#define SDFEM_STAB      1
#define UPWIND_STAB     2
#define CONV_NORM       1
#define H10_NORM        2
#define L2_NORM         3
#define H10_CONV_NORM   4
#define L1_NORM         5
#define LP_NORM         6
#define RES_FCN         7
#define RES_X_CROSS     8

#define SCM_MH         51   /* Mizukami, Hughes (1985) */
#define SCM_HMM        52   /* Hughes, Mallet, Mizukami (1986) */
#define SCM_TP1        53   /* Tezduyar, Park (1986) */
#define SCM_TP2        54   /* Tezduyar, Park (1986) */
#define SCM_JSW        55   /* Johnson, Schatz, Wahlbin (1987) */
#define SCM_LIN        56
#define SCM_GC1        58   /* Galeao, do Carmo (1988) */
#define SCM_GC2        59   /* Galeao, do Carmo (1988) */
#define SCM_CA         60   /* do Carmo, Alvarez (2003) */
#define SCM_AS         61   /* Almeida, Silva (1997) */
#define SCM_C          62   /* Codina (1993) */
#define SCM_KLR1       63   /* Knopp, Lube, Rapin (2002) */
#define SCM_KLR2       64   /* Knopp, Lube, Rapin (2002) */
#define SCM_BE1        65   /* Burman, Ern (2002) */
#define SCM_BH         66   /* Burman, Hansbo (2004) */
#define SCM_BE2        67   /* Burman, Ern (2005) */
#define SCM_BE3        68   /* Burman, Ern (2005) */
#define SCM_J          69   /* Johnson (1990) */
#define SCM_LP         70   /* Layton, Polman (1996) */
#define SCM_CS         71   /* Codina, Soto (1999) - anisotropic part */
#define SCM_CS1        72   /* Codina, Soto (1999) - isotropic part */
#define SCM_K3         83   /* modified KLR1 */
#define SCM_K4         84   /* modified BE1 */
#define SCM_K6         86   /* modified C */
#define SCM_K6red      87   /* SCM_K6 with eps=0 */
#define SCM_K6red_iso  88   /* isotropic artficial diffusion with SCM_K6red */
#define SCM_K8         90
#define SCM_K9         91

/* array entries in node data */
#define U                       0
#define UU                      1
#define F                       2
#define D                       3
#define Q                       4
#define R                       5
#define B                       6
#define W                       7
#define X                       8
#define Y                       9
#define P                       10
#define G                       11
#define J                       12
#define K                       13
#define L                       14
#define N                       15

/* array entries in element data */
#define EU			0
#define EB			1
#define EF			2
#define EQ			3
#define ER			4
#define EH                      5
#define EW			6
#define EX			7
#define EY			8
#define EP                      9
#define EG                      10
#define EI                      11
#define EJ                      12

/* array entries in link/node diag data */
#define A       	0
#define A1              1
#define A2              2
#define A3              3
#define A4              4
#define H1              5
#define H2              6
#define H3              7
#define H4              8
#define BB      	0

/* types of vectors */
#define Q_SN        1  /*  scalars in nodes  */
#define Q_VN        2  /*  vectors in nodes  */
#define Q_DVN       4  /*  dvectors in nodes  */
#define Q_SF        8  /*  scalars in faces  */
#define Q_VF       16  /*  vectors in faces  */
#define Q_DVF      32  /*  dvectors in faces  */
#define Q_SE       64  /*  scalars in elements  */
#define Q_SNE     128  /*  scalars in nodes of elements  */
#define Q_VE      256  /*  vectors in elements  */
#define Q_DVE     512  /*  dvectors in elements  */
#define Q_VNSFLG 1024  /*  vectors in nodes, scalars in faces, LG data  */
#define Q_SSN    2048  /*  scalars in special nodes  */
#define Q_SSF    4096  /*  scalars in special faces  */
#define Q_MVN    8192  /*  general vectors in nodes  */
#define Q_MVF   16384  /*  general vectors in faces  */
#define Q_MVE   32768  /*  general vectors in elements  */
#define Q_SNSF      9  /*  scalars in nodes, scalars in faces  */
#define Q_SNSE     65  /*  scalars in nodes, scalars in elements  */
#define Q_SFSE     72  /*  scalars in faces, scalars in elements  */
#define Q_SNSFSE   73  /*  scalars in nodes, faces and elements  */
#define Q_VNSF     10  /*  vectors in nodes, scalars in faces  */
#define Q_VNVF     18  /*  vectors in nodes, vectors in faces  */
#define Q_VNVE    258  /*  vectors in nodes, vectors in elements  */
#define Q_VNVFVE  274  /*  vectors in nodes, faces and elements  */
#define Q_SSNSSF 6144  /*  scalars in special nodes, scalars in special faces */
#define Q_SNSSNSSF 6145  /* scalars in nodes, special nodes and special faces */
#define Q_MVNMVF  24576  /*  general vectors in nodes and faces  */  
#define Q_GENERAL 57344  /*  general vectors in nodes, faces and elements  */  

/* structures of matrices */
#define Q_TRANS   1  /*  transposed matrix  */
#define Q_FULL    2  /*  all couplings  */
#define Q_NEBDIAG 4  /*  node matrix is block diagonal with identical blocks  */
#define Q_FEBDIAG 8  /*  face matrix is block diagonal with identical blocks  */
#define Q_DVFEBDIAG       16
#define Q_NEBDIAG_FDIAG   32
#define Q_NEBDIAG_NRED    64
#define Q_NEBDIAG_NFRED  128
#define Q_NEBDIAG_EDIAG  256
#define Q_FIRST_DV       512 /* only 1st component of a dvector multiplied    */
#define Q_EBDIAG        1024 /* matrix is block diagonal with identical blocks*/
#define Q_NE_ZERO       2048
#define Q_ZERO_DIAG     4096
#define Q_BDAUGMENT     8192 /* augmented B matrix */

/* NTYPE */
#define BNDRY_TYPE_BITS 126
#define GW              2
#define GLG             4
#define GLS             8
#define GWW             16

#define NTYPE(n)        (n)->myvertex->type


/* FTYPE */
#define TOP_LEVEL_BITS  31
#define ORIENT          32 /* (FTYPE(f) & ORIENT) == 0 ... orientation of the
                           normal vector to f according to vertices numbering */
#define OR_SET          64 /* (FTYPE(f) & OR_SET == 0) ... orientation not set*/
#define BNDRY_FACE_BIT  128
#define IBNDRY_FACE_BIT 256
#define AUX_FACE_BIT    512

#define FTYPE(f)             (f)->type

#define BND_EDGE                 1
#define EDGE_TO_MIDDLE           2
#define EDGE_TO_OLD              4
#define EDGE_TO_DIVIDE           8
#define DIVIDED_EDGE             16
#define TREATED_EDGE             32

#define PI              3.1415926535897932384626433832795
#define PI2             6.283185307179586476925286766559

#define FIFO_NAME "/tmp/myfifo"
INT PIPE_FOR_GNUPLOT_OPENED=NO, PIPE_FD_FOR_GNUPLOT;

