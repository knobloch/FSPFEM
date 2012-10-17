/******************************************************************************/
/*                                                                            */
/*                              data structures                               */
/*                                                                            */
/******************************************************************************/

   struct vertex {

        UNS type;
        FLOAT x[DIM];                 /* vertex coordinates                   */
        struct node *topnode;
#if SLIDING_POINTS == YES
        FLOAT ref_x[DIM];
        FLOAT phi;
#endif
   } ;

   struct bpoint {                    /* point on the boundary                */

        UNS type;
        FLOAT x[DIM];                 /* coordinates                          */
   } ;

/*----------------------------------------------------------------------------*/

#if MOVING_BOUNDARY == YES
#define NODE_NEW_VERTEX   struct vertex *newvertex;
#define FACE_NEW_VERTEX   struct vertex *newvertex;  /* if the face lies on the                                        boundary,otherwise the pointer is null */
#define FACE_VERTEX       struct vertex *myvertex;   /* coordinates of the 
                                         point on the profile of the boundary */
#else
#define NODE_NEW_VERTEX
#define FACE_NEW_VERTEX
#define FACE_VERTEX
#endif

#if DATA_S & PREVIOUS_NODE
#define PREV_NODE    struct node *prev;
#else
#define PREV_NODE  
#endif
        
#if DATA_S & PREVIOUS_FACE
#define PREV_FACE    struct face *prev;
#else
#define PREV_FACE  
#endif
        
#if DATA_S & N_LINK_TO_NODES                  /* list of neighbour nodes  */
#define NN_LINK      struct link *start;                                       \
                     struct link *tstart;
#else
#define NN_LINK    
#endif
        
#if DATA_S & N_LINK_TO_FACES                  /* list of neighbour faces  */
#define NF_LINK      struct nflink *nfstart;                                   \
                     struct nflink *tnfstart;
#else
#define NF_LINK    
#endif
        
#if DATA_S & N_LINK_TO_ELEMENTS               /* list of neighbour elements */
#define NE_LINK      struct nelink *nestart;
#else
#define NE_LINK    
#endif
        
#if DATA_S & F_LINK_TO_FACES                  /* list of neighbour faces  */
#define FF_LINK      struct flink *fstart;                                     \
                     struct flink *tfstart;
#else
#define FF_LINK    
#endif
        
#if DATA_S & F_LINK_TO_NODES                  /* list of neighbour nodes  */
#define FN_LINK      struct fnlink *fnstart;                                   \
                     struct fnlink *tfnstart;
#else
#define FN_LINK    
#endif

#if DATA_S & SPECIAL_NODES_AND_FACES
#define SPEC_NODE    struct snode *s_node;
#define SPEC_FACE    struct sface *s_face;
#else
#define SPEC_NODE    
#define SPEC_FACE    
#endif

#if DATA_STR & LG_DATA
#define NLG_LINK     struct nlglink *lgstart;
#define NLG_DATA     struct lgdata *lgd;
#define FLG_LINK     struct flglink *lgstart;
#else
#define NLG_LINK     
#define NLG_DATA      
#define FLG_LINK     
#endif

#if N_DATA & NODE_BD_DATA
#define NBD_DATA                FLOAT *sbddata;
#define NBD_MATR_1xD    FLOAT (*ddiag)[DIM];
#else
#define NBD_DATA
#define NBD_MATR_1xD
#endif

#if N_DATA & NODE_IFLAG
#define N_IFLAG      INT iflag;
#else
#define N_IFLAG
#endif

#if N_DATA & NODE_ITYPE
#define N_ITYPE           UNS itype;
#define N_ITYPE_INDICES   INT nn_i; UNS nn_f;
#else
#define N_ITYPE
#define N_ITYPE_INDICES
#endif

#if N_DATA & NUMBER_OF_N_NEL
#define N_NEL        UNS nel;
#else
#define N_NEL  
#endif

#if N_DATA & ONE_NODE_MATR 
#define N_MATR_1x1   FLOAT diag[S_MATR];
#else
#define N_MATR_1x1      
#endif

#if N_DATA & IxD_NODE_MATR 
#define N_MATR_1xD   FLOAT bdiag[B_MATR][DIM];
#else
#define N_MATR_1xD      
#endif

#if N_DATA & DxD_NODE_MATR 
#define N_MATR_DxD   FLOAT a[V_MATR][DIM][DIM];
#else
#define N_MATR_DxD      
#endif

#if N_DATA & NxN_NODE_MATR 
#define N_MATR_NxN   FLOAT ann[V_MATR][N_OF_NODE_FUNC][N_OF_NODE_FUNC];
#else
#define N_MATR_NxN      
#endif

#if N_DATA & ONE_NODE_FACE_MATR 
#define NF_MATR_1x1  FLOAT diag[S_MATR];
#else
#define NF_MATR_1x1      
#endif

#if N_DATA & Dx1_NODE_FACE_MATR 
#define NF_MATR_Dx1  FLOAT vdiag[V_MATR][DIM];
#else
#define NF_MATR_Dx1      
#endif

#if N_DATA & DxD_NODE_FACE_MATR 
#define NF_MATR_DxD  FLOAT a[V_MATR][DIM][DIM];
#else
#define NF_MATR_DxD      
#endif

#if N_DATA & NxM_NODE_FACE_MATR 
#define NF_MATR_NxM  FLOAT ann[V_MATR][N_OF_NODE_FUNC][N_OF_FACE_FUNC];
#else
#define NF_MATR_NxM      
#endif

#if N_DATA & Ix2D_NODE_FACE_MATR 
#define NF_MATR_Ix2D  FLOAT bddiag[V_MATR][2][DIM];
#else
#define NF_MATR_Ix2D
#endif

#if N_DATA & SCALAR_NODE_DATA 
#define N_SCALAR_DATA   FLOAT sdata[NS_VECT];
#else
#define N_SCALAR_DATA      
#endif

#if N_DATA & VECTOR_NODE_DATA 
#define N_VECTOR_DATA   FLOAT vdata[NV_VECT][DIM];
#else
#define N_VECTOR_DATA      
#endif

#if N_DATA & MVECTOR_NODE_DATA 
#define N_MVECTOR_DATA  FLOAT mvdata[NM_VECT][N_OF_NODE_FUNC];
#else
#define N_MVECTOR_DATA      
#endif

#if F_DATA & FACE_BD_DATA
#define FBD_DATA                FLOAT *sbddata;
#define FBD_MATR_1xD    FLOAT (*ddiag)[DIM];
#else
#define FBD_DATA
#define FBD_MATR_1xD
#endif

#if F_DATA & FACE_IFLAG
#define F_IFLAG      INT iflag;
#else
#define F_IFLAG
#endif

#if F_DATA & FACE_MIDDLE
#define F_MIDDLE     struct node *midpoint;
#else
#define F_MIDDLE
#endif

#if F_DATA & CURVED_FACE_MIDDLE
#define F_CURVED_MIDDLE     struct bpoint *c_midpoint;
#else
#define F_CURVED_MIDDLE
#endif

#if F_DATA & ONE_FACE_MATR 
#define F_MATR_1x1   FLOAT diag[V_MATR];
#else
#define F_MATR_1x1      
#endif

#if F_DATA & IxD_FACE_MATR 
#define F_MATR_1xD   FLOAT bdiag[F_MATR][DIM];
#else
#define F_MATR_1xD      
#endif

#if F_DATA & Ix2D_FACE_MATR 
#define F_MATR_1x2D   FLOAT bddiag[F_MATR][2][DIM];
#else
#define F_MATR_1x2D      
#endif

#if F_DATA & DxD_FACE_MATR 
#define F_MATR_DxD   FLOAT a[F_MATR][DIM][DIM];
#else
#define F_MATR_DxD      
#endif

#if F_DATA & NxN_FACE_MATR 
#define F_MATR_NxN   FLOAT ann[F_MATR][N_OF_FACE_FUNC][N_OF_FACE_FUNC];
#else
#define F_MATR_NxN      
#endif

#if F_DATA & ONE_FACE_NODE_MATR 
#define FN_MATR_1x1  FLOAT diag[V_MATR];
#else
#define FN_MATR_1x1      
#endif

#if F_DATA & IxD_FACE_NODE_MATR 
#define FN_MATR_1xD  FLOAT vdiag[V_MATR][DIM];
#else
#define FN_MATR_1xD      
#endif

#if F_DATA & DxD_FACE_NODE_MATR 
#define FN_MATR_DxD  FLOAT a[F_MATR][DIM][DIM];
#else
#define FN_MATR_DxD      
#endif

#if F_DATA & MxN_FACE_NODE_MATR 
#define FN_MATR_MxN  FLOAT ann[F_MATR][N_OF_FACE_FUNC][N_OF_NODE_FUNC];
#else
#define FN_MATR_MxN      
#endif

#if F_DATA & SCALAR_FACE_DATA
#define F_SCALAR_DATA   FLOAT sdata[FS_VECT];
#else
#define F_SCALAR_DATA
#endif

#if F_DATA & VECTOR_FACE_DATA
#define F_VECTOR_DATA   FLOAT vdata[FV_VECT][DIM];
#else
#define F_VECTOR_DATA
#endif

#if F_DATA & MVECTOR_FACE_DATA 
#define F_MVECTOR_DATA  FLOAT mvdata[FM_VECT][N_OF_FACE_FUNC];
#else
#define F_MVECTOR_DATA      
#endif

#if F_DATA & DVECTOR_FACE_DATA
#define F_DVECTOR_DATA   FLOAT dvdata[FDV_VECT][2][DIM];
#else
#define F_DVECTOR_DATA
#endif

   struct node {                      /* level dependent part of a vertex     */

        struct node *succ;            /* linked list of nodes per level       */
        INT index;
        INT index2;
        N_IFLAG
        N_ITYPE
        N_ITYPE_INDICES
        N_NEL
        struct node *father;          /* father node                          */
        struct node *son;             /* node on finer level (NULL if none)   */
        struct vertex *myvertex;      /* corresponding vertex structure       */
        NODE_NEW_VERTEX
        PREV_NODE
        SPEC_NODE    
        NN_LINK
        NF_LINK
        NE_LINK
        NLG_LINK
        NLG_DATA
        N_SCALAR_DATA
        N_VECTOR_DATA
        N_MVECTOR_DATA      
        N_MATR_1x1
        N_MATR_1xD
        N_MATR_DxD
        N_MATR_NxN

        NBD_DATA
        NBD_MATR_1xD
   } ;

   struct link {

        struct node *nbnode;          /* neighbour node                       */
        struct link *next;
        INT flag;
        N_MATR_1x1
        N_MATR_1xD
        N_MATR_DxD
        N_MATR_NxN
        FLOAT flux;                   /* needed for upwinding                 */

        NBD_MATR_1xD
   } ;

   struct nflink{

        struct face *nbface;
        struct nflink *next;
        NF_MATR_1x1      
        NF_MATR_Dx1      
        NF_MATR_DxD      
        NF_MATR_NxM      
        NF_MATR_Ix2D

        NBD_MATR_1xD
   };

   struct nelink{

        struct element *nbel;
        struct nelink *next;
   };

   struct face{
        
        struct face *succ;
        UNS type;
        INT index;
        INT index2;
        FACE_VERTEX
        FACE_NEW_VERTEX
        F_IFLAG
        struct face *father;          /* father face on coarser grid          */
  	struct face *sons[FSONS];     /* face tree                            */
        PREV_FACE
        SPEC_FACE    
        F_MIDDLE
        F_CURVED_MIDDLE
        FF_LINK
        FN_LINK
        FLG_LINK
        F_SCALAR_DATA
        F_VECTOR_DATA
        F_MVECTOR_DATA      
        F_DVECTOR_DATA
        F_MATR_1x1
        F_MATR_1xD
        F_MATR_1x2D      
        F_MATR_DxD
        F_MATR_NxN      

        FBD_DATA
        FBD_MATR_1xD
   };

   struct flink{

        struct face *nbface;
        struct flink *next;
        F_MATR_1x1
        F_MATR_1xD
        F_MATR_1x2D      
        F_MATR_DxD
        F_MATR_NxN      
   };


   struct fnlink{

        struct node *nbnode;
        struct fnlink *next;
        FN_MATR_1x1
        FN_MATR_1xD
        FN_MATR_DxD
        FN_MATR_MxN

        FBD_MATR_1xD
   };

   struct snode{

        struct snode *succ;           /*  linked list of snodes per level     */
        struct node *n0, *n1, *n2;    /*  n0 is a common node of f1, f2;      */
        struct face *f1, *f2;         /*  n0, n1 are nodes of f1;             */
        struct element *pel1, *pel2;  /*  n0, n2 are nodes of f2;             */
        double ssdata[SNS_VECT];      /*  pel1 is adjacent to f1, pel2 to f2  */
        double ann0[D_MATR][DIM], ann1[D_MATR][DIM], ann2[D_MATR][DIM];
        double anf1[D_MATR][DIM], anf2[D_MATR][DIM];
   };

   struct sface{

        struct sface *succ;     /*  linked list of sfaces per level  */
        struct node *n1, *n2;
        struct face *f;         /*  n1, n2 are nodes of f  */
        struct element *pel;    /*  pel is adjacent to f   */
        double ssdata[SFS_VECT];
        double aff[D_MATR][DIM], afn1[D_MATR][DIM], afn2[D_MATR][DIM];
   };

   struct p_node{

        struct p_node *succ;          /*  linked list of p_nodes per level    */
        struct node *n1, *n2;         /*  corresponding nodes                 */
   };

   struct p_face{

        struct p_face *succ;          /*  linked list of p_faces per level    */
        struct face *f1, *f2;         /*  corresponding faces                 */
   };


/*----------------------------------------------------------------------------*/

#if E_DATA & SCALAR_ELEMENT_DATA
#define E_SCALAR_DATA   FLOAT sdata[EL_VECT];
#else
#define E_SCALAR_DATA        
#endif

#if E_DATA & VECTOR_ELEMENT_DATA
#define E_VECTOR_DATA   FLOAT vdata[EL_VECTD][DIM];
#else
#define E_VECTOR_DATA        
#endif

#if E_DATA & MVECTOR_ELEMENT_DATA
#define E_MVECTOR_DATA   FLOAT mvdata[EM_VECT][N_OF_ELEM_FUNC];
#else
#define E_MVECTOR_DATA        
#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES
#define E_SCALAR_NDATA   FLOAT sndata[EL_VECTN][NVERT];
#else
#define E_SCALAR_NDATA        
#endif

#if E_DATA & ExE_MATR
#define E_MATR_1x1     FLOAT ee_diag[S_MATR];      
#else
#define E_MATR_1x1
#endif

#if E_DATA & ExN_MATR
#define EN_MATR_1x1     FLOAT en_diag[S_MATR][NVERT];      
#else
#define EN_MATR_1x1
#endif

#if E_DATA & NxE_MATR
#define NE_MATR_1x1     FLOAT ne_diag[S_MATR][NVERT];      
#else
#define NE_MATR_1x1
#endif

#if E_DATA & ExDN_MATR
#define EN_MATR_1xD     FLOAT bn[S_MATR][NVERT][DIM]; 
#else
#define EN_MATR_1xD
#endif

#if E_DATA & ExF_MATR
#define EF_MATR_1x1     FLOAT ef_diag[S_MATR][SIDES];      
#else
#define EF_MATR_1x1
#endif

#if E_DATA & FxE_MATR
#define FE_MATR_1x1     FLOAT fe_diag[S_MATR][SIDES];      
#else
#define FE_MATR_1x1
#endif

#if E_DATA & ExDF_MATR
#define EF_MATR_1xD     FLOAT bdf[S_MATR][SIDES][DIM];      
#else
#define EF_MATR_1xD
#endif

#if E_DATA & MxM_E_E_MATR
#define E_MATR_MxM     FLOAT ee_mm[S_MATR][N_OF_ELEM_FUNC][N_OF_ELEM_FUNC];
#else
#define E_MATR_MxM
#endif

#if E_DATA & MxN_E_N_MATR
#define EN_MATR_MxN  FLOAT en_mn[S_MATR][NVERT][N_OF_ELEM_FUNC][N_OF_NODE_FUNC];
#else
#define EN_MATR_MxN
#endif

#if E_DATA & NxM_N_E_MATR
#define NE_MATR_NxM  FLOAT ne_nm[S_MATR][NVERT][N_OF_NODE_FUNC][N_OF_ELEM_FUNC];
#else
#define NE_MATR_NxM
#endif

#if E_DATA & MxN_E_F_MATR
#define EF_MATR_MxN  FLOAT ef_mn[S_MATR][NVERT][N_OF_ELEM_FUNC][N_OF_FACE_FUNC];
#else
#define EF_MATR_MxN
#endif

#if E_DATA & NxM_F_E_MATR
#define FE_MATR_NxM  FLOAT fe_nm[S_MATR][NVERT][N_OF_FACE_FUNC][N_OF_ELEM_FUNC];
#else
#define FE_MATR_NxM
#endif

#if E_DATA & ExFxDN_MATR
#define EFN_MATR_1x1xD FLOAT bfdn[S_MATR][SIDES][NVERT][DIM];      
#else
#define EFN_MATR_1x1xD
#endif

#if E_DATA & ExFxDF_MATR
#define EFF_MATR_1x1xD FLOAT bfdf[S_MATR][SIDES][SIDES][DIM];      
#else
#define EFF_MATR_1x1xD
#endif

#if E_DATA & ExFx2DF_MATR
#define EFF_MATR_1x1x2D FLOAT bdfs[S_MATR][SIDES][SIDES][2][DIM];      
#else
#define EFF_MATR_1x1x2D
#endif

#if E_DATA & E_E_FNEIGHBOURS
#define NB_ELEMS        struct element *nbel[SIDES];
#else
#define NB_ELEMS
#endif

#if E_DATA & E_E_NEIGHBOURS
#define E_LINK          struct elink *estart;
#else
#define E_LINK          
#endif

#if E_DATA & ELEM_ITYPE
#define E_ITYPE      UNS itype;
#else
#define E_ITYPE
#endif

#if E_DATA & E_TAU
#define ELEM_TAU        DOUBLE tau;
#else
#define ELEM_TAU
#endif

#if E_DATA & E_TAU_SOLD
#define ELEM_TAU_SOLD   DOUBLE tau_sold;
#else
#define ELEM_TAU_SOLD
#endif

   struct element {                   

        struct element *succ; /* list of elements                             */
        INT status;           /* 0 normal, 1 further refinement not possible
	                                 (green element), -1 fictive element  */
        INT eflag;    /* 0 normal, 1 new green elements, 2 primary refinement,
	                 3 uniform refinement of green elements, 
	                 4 double refinement of green elements                */

        INT index2;
        struct node *n[NVERT];          /* nodes of that element "n[i+1]>n[i] */
        struct face *f[SIDES]; /* faces of element; n[i] not contained in f[i]*/
        struct element *father;         /* father element on coarser grid     */
  	struct element *sons[SONS];     /* element tree                       */
        E_ITYPE
        ELEM_TAU
        ELEM_TAU_SOLD
        NB_ELEMS
        E_LINK
        E_SCALAR_DATA
        E_VECTOR_DATA
        E_MVECTOR_DATA        
        E_SCALAR_NDATA        
        E_MATR_1x1
        EN_MATR_1x1
        NE_MATR_1x1
        EF_MATR_1x1
        FE_MATR_1x1
        EN_MATR_1xD
        EF_MATR_1xD
        E_MATR_MxM
        EN_MATR_MxN
        NE_MATR_NxM
        EF_MATR_MxN
        FE_MATR_NxM
        EFN_MATR_1x1xD
        EFF_MATR_1x1xD
        EFF_MATR_1x1x2D
   } ;

   struct elink{

        struct element *nbel;
        struct elink *next;
   } ;

/*----------------------------------------------------------------------------*/

   struct lgdata {                          
	
        struct lgnlink *nstart;
        struct lgflink *fstart;
        struct lglglink *lgstart;
        FLOAT  n[DIM];  /* discrete normal vector  */
        FLOAT t1[DIM];  /*  t1 \perp n  */
        FLOAT t2[DIM];  /*  t2 \perp n  and  t2 \perp t1  */
        FLOAT tdata[LG_VECT][DIM1];
        FLOAT a[V_MATR][DIM1][DIM1];
   } ;

   struct nlglink {                          
	
        struct node *nbnode;            
        struct nlglink *next;  
        FLOAT a[V_MATR][DIM][DIM1];
   } ;

   struct lgnlink {                          
	
        struct node *nbnode;            
        struct lgnlink *next;  
        FLOAT a[V_MATR][DIM1][DIM];
   } ;

   struct lgflink {                          
	
        struct face *nbface;            
        struct lgflink *next;  
        FLOAT a[V_MATR][DIM1];
   } ;

   struct lglglink {                          
	
        struct node *nbnode;            
        struct lglglink *next;  
        FLOAT a[V_MATR][DIM1][DIM1];
   } ;

   struct flglink {                          
	
        struct node *nbnode;            
        struct flglink *next;  
        FLOAT a[V_MATR][DIM1];
   } ;

/*----------------------------------------------------------------------------*/

typedef double (*DOUBLE_FUNC)();
 
struct finite_element {

   UNS domain;   /* SIMPLEX, CUBE */
   UNS k, n, m;  /* k=0,1 ... #basis f. at a vertex
                               n     ... #basis f. inside an edge
                               m     ... #basis f. inside the element */
   DOUBLE_FUNC *nbasis;  /* basis functions at vertices */
   DOUBLE_FUNC *fbasis;  /* basis functions assigned to nodes on edges */
   DOUBLE_FUNC *ebasis;  /* basis functions assigned to nodes inside elements */
   DOUBLE_FUNC *nbasis_0;
   DOUBLE_FUNC *fbasis_0;
   DOUBLE_FUNC *ebasis_0;
   DOUBLE_FUNC *nbasis_1;
   DOUBLE_FUNC *fbasis_1;
   DOUBLE_FUNC *ebasis_1;
   DOUBLE_FUNC *nbasis_00;
   DOUBLE_FUNC *fbasis_00;
   DOUBLE_FUNC *ebasis_00;
   DOUBLE_FUNC *nbasis_01;
   DOUBLE_FUNC *fbasis_01;
   DOUBLE_FUNC *ebasis_01;
   DOUBLE_FUNC *nbasis_11;
   DOUBLE_FUNC *fbasis_11;
   DOUBLE_FUNC *ebasis_11;
   DOUBLE_FUNC *ndofs;   /* degrees of freedom assigned to vertices */
   DOUBLE_FUNC *fdofs;   /* degrees of freedom assigned to nodes on edges */
   DOUBLE_FUNC *edofs;   /* dofs assigned to nodes inside elements */

} ;

struct quadr_rule_1d {
   int    n;                        /*  number of points  */
   double *points, *weights;
} ;

struct quadrature_rule {
   int    n;                        /*  number of points  */
   double (*points)[DIM], *weights;
} ;

struct ref_mapping {
   UNS type;     /*  type of ref_mapping  */ 
   UNS domain;   /* SIMPLEX, CUBE         */
   FLOAT area;

   FLOAT p1_b[DIM][DIM], p1_a[DIM][DIM], p1_c[DIM], p1_jac;  /* see 
    P1_reference_mapping0 and P1_reference_mapping0_with_inverse in triang.h  */

   FLOAT q1_b[2][2][DIM2], q1_a[DIM][DIM], q1_c[DIM], q1_alpha[DIM], 
         q1_jac[DIM2];  /* see Q1_reference_mapping_with_inverse in triang.h  */
} ;

typedef struct finite_element   FINITE_ELEMENT;
typedef struct quadr_rule_1d    QUADR_RULE_1D;
typedef struct quadrature_rule  QUADRATURE_RULE;
typedef struct ref_mapping      REF_MAPPING;

/*----------------------------------------------------------------------------*/

#define N_NN_ITYPE 10000
#define MAX_ITYPES 10  /* DO NOT CHANGE WITHOUT CHANGING THE PROGRAMM */
struct itypes_positions {
   int n;            /*  number of itypes  */
   int n_e;          /*  number of element itypes  */
   int nn;           /*  number of used entries of nn_link_itypes */
   int it[MAX_ITYPES];/* all itypes of nodes; it[0] ... itype of boundary nodes,
                        it[1] > it[2] > ... > it[n-1] */
   int it_e[MAX_ITYPES];     /*  element itypes  */
   struct node *firstN, *ldbn, *firstNode;
   struct node *pn1[MAX_ITYPES], *pn2[MAX_ITYPES];  /* pn1[i]/pn2[i] ... 
                                            first/last node with itype it[i] */
   struct element *pe1[MAX_ITYPES], *pe2[MAX_ITYPES];
} ;

typedef struct itypes_positions ITYPES_POSITIONS;
ITYPES_POSITIONS itypes_pos;
struct link *nn_link_itypes[N_NN_ITYPE];

/*----------------------------------------------------------------------------*/
INT Ap[MAX_ROW], Ai[MAX_ENT], Nj;
DOUBLE Ax[MAX_ENT], Rhs[MAX_ROW], Sol[MAX_ROW];
/*----------------------------------------------------------------------------*/

   struct grid {

        UNS level;    
        INT minNodeIndex, maxNodeIndex; 
        INT minFaceIndex, maxFaceIndex; 
        struct node *firstN;                 /*  first node in the list       */
        struct node *fdbn, *ldbn;
        struct node *firstNode, *lastNode;
        struct face *firstF;                 /*  first face in the list       */
        struct face *fdbf, *ldbf;
        struct face *firstFace, *lastFace;
        struct snode *firstSN;
        struct sface *firstSF;
        struct p_node *firstPN;              /*  for periodic b.c.            */
        struct p_face *firstPF;              /*  for periodic b.c.            */
        struct element *firstElement;
        struct element *fdbe; 
        struct grid *coarser, *finer, *first_grid, *top_grid;
        FLOAT alpha;                    /*  estimate of spectral radius of A  */
   } ;

   struct multigrid {
 
        UNS toplevel;
        INT nodeNumber;                 /* number of nodes                    */
        INT faceNumber;                 /* number of all faces                */
        INT elemNumber;                 /* number of all elements             */
        INT fictNumber;                 /* number of fictive faces            */
        INT fictElNumber;               /* number of fictive elements         */
        struct grid *grids[MAXLEVEL];
   };

typedef struct vertex      VERTEX;
typedef struct bpoint      BPOINT;
typedef struct node        NODE;
typedef struct snode       SNODE;
typedef struct p_node      P_NODE;
typedef struct element     ELEMENT;
typedef struct link        LINK;
typedef struct elink       ELINK;
typedef struct face        FACE;
typedef struct sface       SFACE;
typedef struct p_face      P_FACE;
typedef struct lgdata      LGDATA;
typedef struct nlglink     NLGLINK;
typedef struct lgnlink     LGNLINK;
typedef struct lgflink     LGFLINK;
typedef struct lglglink    LGLGLINK;
typedef struct flink       FLINK;
typedef struct nflink      NFLINK;
typedef struct nelink      NELINK;
typedef struct fnlink      FNLINK; 
typedef struct flglink     FLGLINK; 
typedef struct grid        GRID; 
typedef struct multigrid   MULTIGRID;

