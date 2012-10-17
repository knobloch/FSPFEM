/******************************************************************************/
/*                                                                            */
/*                                 macros 2                                   */
/*                                                                            */
/******************************************************************************/

#if P_SPACE == NONE

#define P_TYPE          0
#define B_STRUCT        0

#if (U_SPACE == P1C) || (U_SPACE == Q1C) || (U_SPACE == MINI_L_DIV_FR)

#if U_STRUCTURE == SCALAR

#define U_TYPE          Q_SN
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SN

#else

#define U_TYPE          Q_VN
#define ILU_TYPE        Q_VN
#if DATA_STR & KORN_MATRIX
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#else
#define A_STRUCT        Q_NEBDIAG
#define ILU_STRUCT      Q_NEBDIAG
#endif

#endif

#if U_SPACE == P1C
#define F_NAME          "Aug_p1c"
#elif U_SPACE == Q1C
#define F_NAME          "Aug_q1c"
#elif U_SPACE == MINI_L_DIV_FR
#define F_NAME          "Aug_mini_l_div_fr"
#endif
#endif /*  U_SPACE == P1C || U_SPACE == Q1C || U_SPACE == MINI_L_DIV_FR  */

#if (U_SPACE == Q2C) && (U_STRUCTURE == SCALAR)

#define U_TYPE          Q_SNSFSE
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SNSFSE

#endif

#if (U_SPACE == GP1C) || (U_SPACE == GP1X3C) || (U_SPACE == GQ1C) || (U_SPACE == GP1C_ELBUB) || (U_SPACE == GQ1C_ELBUB) || (U_SPACE == GQ1X4C) || (U_SPACE == GP2C) || (U_SPACE == GP2X3C) || (U_SPACE == GP2C_3ELBUB) || (U_SPACE == GP2C_6ELBUB) || (U_SPACE == GQ2C) || (U_SPACE == GQ2X4C) || (U_SPACE == GQ2B3C) || (U_SPACE == GQ2C_2ELBUB) || (U_SPACE == GQ2C_3ELBUB)

#if U_STRUCTURE == SCALAR

#define U_TYPE          Q_GENERAL
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_GENERAL

#else

#endif

#if U_SPACE == GP1C
#define F_NAME          "Aug_gp1c"
#elif U_SPACE == GQ1C
#define F_NAME          "Aug_gq1c"
#elif U_SPACE == GP2C
#define F_NAME          "Aug_gp2c"
#endif
#endif /*  U_SPACE == GP1C || U_SPACE == GQ1C || U_SPACE == GP2C */

#if U_SPACE == P1_NC
#define U_TYPE          Q_SF
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SF
#define F_NAME          "Aug_nc"
#endif

#if U_SPACE == P1_MOD
#define U_TYPE          Q_VF
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#define ILU_TYPE        Q_VF
#define F_NAME          "Aug_mod"
#endif

#if (U_SPACE == P1C_ELBUB) && (U_STRUCTURE == SCALAR)
#define U_TYPE          Q_SNSE
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SNSE
#define F_NAME          "Aug_p1cb"
#endif

#if U_SPACE == P2C
#define U_TYPE          Q_SNSF
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SNSF
#define F_NAME          "Aug_p2c"
#endif

#if U_SPACE == Q1ROT
#if U_STRUCTURE == SCALAR
#define U_TYPE          Q_SF
#define A_STRUCT        0
#define ILU_STRUCT      0
#define ILU_TYPE        Q_SF
#elif !(DATA_STR & KORN_MATRIX) 
#define U_TYPE          Q_VF
#define A_STRUCT        Q_FEBDIAG
#define ILU_STRUCT      Q_FEBDIAG
#define ILU_TYPE        Q_SF
#else
#define U_TYPE          Q_VF
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#define ILU_TYPE        Q_VF
#endif
#define F_NAME          "Aug_q1rot"
#endif

#endif

#if U_SPACE == P1C && P_SPACE == P0

#define U_TYPE          Q_VN
#define P_TYPE          Q_SE

#if !(DATA_STR & KORN_MATRIX) 
#define A_STRUCT        Q_NEBDIAG
#define ILU_STRUCT      Q_NEBDIAG
#else  /*  (DATA_STR & KORN_MATRIX)  */
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#endif

#define B_STRUCT        0
#define C_STRUCT        1
#define ILU_TYPE        Q_VN
#define F_NAME          "Aug"

#endif

#if (U_SPACE == P1C_FBUB || U_SPACE == P1C_NEW_FBUB) && (P_SPACE == P0)

#define U_TYPE          Q_VNSF
#define P_TYPE          Q_SE

/*------------------ not Korn matrix  -----------------------*/

#if !(DATA_STR & KORN_MATRIX) 

#if (DATA_STR & REDUCED) && (U_SPACE == P1C_NEW_FBUB)
#define A_STRUCT        Q_NEBDIAG_NFRED
#define ILU_STRUCT      Q_NEBDIAG_NFRED
#elif (DATA_STR & REDUCED) && !(U_SPACE == P1C_NEW_FBUB)
#define A_STRUCT        Q_NEBDIAG_NRED
#define ILU_STRUCT      Q_NEBDIAG_NRED
#else
#define A_STRUCT        Q_NEBDIAG
#define ILU_STRUCT      Q_NEBDIAG
#endif

/*-------------------- Korn matrix --------------------------*/

#else  /*  (DATA_STR & KORN_MATRIX)  */

#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL

#endif  /*  (DATA_STR & KORN_MATRIX)  */

#define B_STRUCT        0
#define ILU_TYPE        Q_VNSF
#define F_NAME          "Aug"

#endif

#if (U_SPACE == P1C_FBUB_LG) && (P_SPACE == P0)

#define U_TYPE          Q_VNSFLG
#define P_TYPE          Q_SE

#if DATA_STR & KORN_MATRIX 

#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL

#else

#define A_STRUCT        Q_NEBDIAG
#define ILU_STRUCT      Q_NEBDIAG

#endif

#define B_STRUCT        0
#define ILU_TYPE        0
#define F_NAME          "Aug"

#endif

#if (U_SPACE == P1C_ELBUB) && (P_SPACE == P1C)

#define U_TYPE          Q_VNVE
#define P_TYPE          Q_SN
#define A_STRUCT        Q_NEBDIAG_EDIAG
#define ILU_STRUCT      Q_NEBDIAG_EDIAG
#define B_STRUCT        0
#define ILU_TYPE        Q_VNVE
#define F_NAME          "Aug_mini"

#endif

#if (U_SPACE == P2C || U_SPACE == IP2C) && (P_SPACE == P0 || P_SPACE == P1C || P_SPACE == P1C_P2L || P_SPACE == IP1C || P_SPACE == P1_DISC || P_SPACE == IP1_DISC)

#if P_SPACE == P0
#define P_TYPE          Q_SE
#elif P_SPACE == P1C || P_SPACE == IP1C
#define P_TYPE          Q_SN
#elif P_SPACE == P1_DISC || P_SPACE == IP1_DISC
#define P_TYPE          Q_SNE
#elif P_SPACE == P1C_P2L
#define P_TYPE          Q_SNSSNSSF
#endif

#define U_TYPE          Q_VNVF
#define F_NAME          "Aug_TaylorHood"
#define B_STRUCT        0
#define ILU_TYPE        Q_VNVF

#if DATA_STR & KORN_MATRIX 
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#else
#define A_STRUCT        Q_EBDIAG
#define ILU_STRUCT      Q_EBDIAG
#endif

#endif

#if (U_SPACE == P2C_ELBUB || U_SPACE == IP2C_ELBUB) && (P_SPACE == P0 || P_SPACE == P1C || P_SPACE == IP1C || P_SPACE == P1_DISC || P_SPACE == IP1_DISC)

#if P_SPACE == P0
#define P_TYPE          Q_SE
#define B_STRUCT        Q_NE_ZERO
#elif P_SPACE == P1C || P_SPACE == IP1C
#define P_TYPE          Q_SN
#define B_STRUCT        0
#elif P_SPACE == P1_DISC || P_SPACE == IP1_DISC
#define P_TYPE          Q_SNE
#define B_STRUCT        0
#elif P_SPACE == P1C_P2L
#define P_TYPE          Q_SNSSNSSF
#define B_STRUCT        0
#endif

#define U_TYPE          Q_VNVFVE
#define A_STRUCT        Q_EBDIAG
#define ILU_STRUCT      Q_EBDIAG
#define ILU_TYPE        Q_VNVFVE
#define F_NAME          "Aug_P2Cbubble"

#endif

#if (U_SPACE == P1_NC || U_SPACE == Q1ROT) && (P_SPACE == P0)

#define U_TYPE          Q_VF
#define P_TYPE          Q_SE
#define B_STRUCT        0
#if U_SPACE == P1_NC
#define F_NAME          "Aug_nc"
#else
#define F_NAME          "Aug_q1rot"
#endif

#if DATA_STR & KORN_MATRIX 
#define A_STRUCT        Q_FULL
#define ILU_STRUCT      Q_FULL
#define ILU_TYPE        Q_VF
#else
#define A_STRUCT        Q_FEBDIAG
#define ILU_STRUCT      Q_FEBDIAG
#define ILU_TYPE        Q_SF
#endif

#endif

#if (U_SPACE_LOW == P1_NC) && (P_SPACE_LOW == P0)

#define T_FOR_U_LOW     121
#define T_FOR_P_LOW     0
#define U_TYPE_LOW      Q_VF
#define P_TYPE_LOW      Q_SE
#define B_STRUCT_LOW    0

#if DATA_STR & KORN_MATRIX
#define A_STRUCT_LOW    Q_FULL
#define ILU_STRUCT_LOW  Q_FULL
#define ILU_TYPE_LOW    Q_VF
#else
#define A_STRUCT_LOW    Q_FEBDIAG
#define ILU_STRUCT_LOW  Q_FEBDIAG
#define ILU_TYPE_LOW    Q_SF
#endif

#endif

#if U_SPACE_LOW == P1C && P_SPACE_LOW == P0

#define T_FOR_U_LOW     121
#define T_FOR_P_LOW     0
#define U_TYPE_LOW      Q_VN
#define P_TYPE_LOW      Q_SE

#if !(DATA_STR & KORN_MATRIX) 
#define A_STRUCT_LOW    Q_NEBDIAG
#define ILU_STRUCT_LOW  Q_NEBDIAG
#else  /*  (DATA_STR & KORN_MATRIX)  */
#define A_STRUCT_LOW    Q_FULL
#define ILU_STRUCT_LOW  Q_FULL
#endif

#define B_STRUCT_LOW    0
#define C_STRUCT_LOW    1
#define ILU_TYPE_LOW    Q_VN

#endif

#if (U_SPACE_LOW == P2C) && (P_SPACE_LOW == P0)

#define T_FOR_U_LOW     121
#define T_FOR_P_LOW     0
#define U_TYPE_LOW      Q_VNVF
#define P_TYPE_LOW      Q_SE
#define A_STRUCT_LOW    Q_EBDIAG
#define ILU_STRUCT_LOW  Q_EBDIAG
#define B_STRUCT_LOW    0
#define ILU_TYPE_LOW    Q_VNVF

#endif

#if (U_SPACE_LOW == P2C) && (P_SPACE_LOW == P1C)

#define T_FOR_U_LOW     T_FOR_U
#define T_FOR_P_LOW     T_FOR_P
#define U_TYPE_LOW      Q_VNVF
#define P_TYPE_LOW      Q_SN
#define B_STRUCT_LOW    0
#define ILU_TYPE_LOW    Q_VNVF

#if DATA_STR & KORN_MATRIX 
#define A_STRUCT_LOW    Q_FULL
#define ILU_STRUCT_LOW  Q_FULL
#else
#define A_STRUCT_LOW    Q_EBDIAG
#define ILU_STRUCT_LOW  Q_EBDIAG
#endif

#endif

#if (U_SPACE == P1_MOD) && (P_SPACE == P0)

#define U_TYPE          Q_DVF
#define P_TYPE          Q_SE
#define A_STRUCT        Q_DVFEBDIAG
#define B_STRUCT        Q_FIRST_DV
#define ILU_TYPE        Q_VF
#define ILU_STRUCT      Q_FULL
#define F_NAME          "Aug_p1mod"

#endif

#if (U_SPACE == P1_MOD) && (P_SPACE == P1_DISC)

#define U_TYPE          Q_DVF
#define P_TYPE          Q_SNE
#define A_STRUCT        Q_DVFEBDIAG
#define B_STRUCT        0
#define ILU_TYPE        Q_VF
#define ILU_STRUCT      Q_FULL
#define F_NAME          "Aug_p1mod"

#endif

#if (U_SPACE == Q1C) || (U_SPACE == GQ1C) || (U_SPACE == GQ1C_ELBUB)
#if (U_SPACE == Q1C) || (U_SPACE == GQ1C)
#define ELEM         q1c_element
#elif U_SPACE == GQ1C_ELBUB
#define ELEM         q1cb_element
#endif
#define REF_MAP      Q1_REF_MAP
#define L2_Q_RULE    qs9
#define H10_Q_RULE   qs9
#define LAPL_Q_RULE  qs9
#define CONV_Q_RULE  qs9
#define REACT_Q_RULE qs9
#define SD_Q_RULE    qs9
#define SC_Q_RULE    qs9
#define RHS_Q_RULE   qs9
#define SDRHS_Q_RULE qs9
#endif

#if (U_SPACE == GQ1X4C) || (U_SPACE == GQ2X4C)
#if U_SPACE == GQ1X4C
#define ELEM         q1x4c_element
#elif U_SPACE == GQ2X4C
#define ELEM         q2x4c_element
#endif
#define REF_MAP      Q1_REF_MAP
#define L2_Q_RULE    qs4x9
#define H10_Q_RULE   qs4x9
#define LAPL_Q_RULE  qs4x9
#define CONV_Q_RULE  qs4x9
#define REACT_Q_RULE qs4x9
#define SD_Q_RULE    qs4x9
#define SC_Q_RULE    qs4x9
#define RHS_Q_RULE   qs4x9
#define SDRHS_Q_RULE qs4x9
#endif

#if (U_SPACE == P1C) || (U_SPACE == GP1C) || (U_SPACE == P1C_ELBUB) || (U_SPACE == GP1C_ELBUB)
#if (U_SPACE == P1C) || (U_SPACE == GP1C)
#define ELEM         p1c_element
#elif (U_SPACE == P1C_ELBUB) || (U_SPACE == GP1C_ELBUB)
#define ELEM         p1cb_element
#endif
#define REF_MAP      P1_REF_MAP
#define L2_Q_RULE    qr8
#define H10_Q_RULE   qr8
#define LAPL_Q_RULE  qr8
#define CONV_Q_RULE  qr8
#define REACT_Q_RULE qr8
#define SD_Q_RULE    qr8
#define SC_Q_RULE    qr8
#define RHS_Q_RULE   qr8
#define SDRHS_Q_RULE qr8
#endif

#if (U_SPACE == P2C) || (U_SPACE == GP2C) 
#define ELEM         p2c_element
#define REF_MAP      P1_REF_MAP
#define L2_Q_RULE    qr8
#define H10_Q_RULE   qr8
#define LAPL_Q_RULE  qr8
#define CONV_Q_RULE  qr8
#define REACT_Q_RULE qr8
#define SD_Q_RULE    qr8
#define SC_Q_RULE    qr8
#define RHS_Q_RULE   qr8
#define SDRHS_Q_RULE qr8
#endif

#if (U_SPACE == GP2C_3ELBUB) 
#define ELEM         p2c3b_element
#define REF_MAP      P1_REF_MAP
#define L2_Q_RULE    qr8
#define H10_Q_RULE   qr8
#define LAPL_Q_RULE  qr8
#define CONV_Q_RULE  qr8
#define REACT_Q_RULE qr8
#define SD_Q_RULE    qr8
#define SC_Q_RULE    qr8
#define RHS_Q_RULE   qr8
#define SDRHS_Q_RULE qr8
#endif

#if (U_SPACE == GP2C_6ELBUB) 
#define ELEM         p2c6b_element
#define REF_MAP      P1_REF_MAP
#define L2_Q_RULE    qr8
#define H10_Q_RULE   qr8
#define LAPL_Q_RULE  qr8
#define CONV_Q_RULE  qr8
#define REACT_Q_RULE qr8
#define SD_Q_RULE    qr8
#define SC_Q_RULE    qr8
#define RHS_Q_RULE   qr8
#define SDRHS_Q_RULE qr8
#endif

#if (U_SPACE == GP1X3C) || (U_SPACE == GP2X3C)
#if U_SPACE == GP1X3C
#define ELEM         p1x3c_element
#elif U_SPACE == GP2X3C
#define ELEM         p2x3c_element
#endif
#define REF_MAP      P1_REF_MAP
#define L2_Q_RULE    qr3x8
#define H10_Q_RULE   qr3x8
#define LAPL_Q_RULE  qr3x8
#define CONV_Q_RULE  qr3x8
#define REACT_Q_RULE qr3x8
#define SD_Q_RULE    qr3x8
#define SC_Q_RULE    qr3x8
#define RHS_Q_RULE   qr3x8
#define SDRHS_Q_RULE qr3x8
#endif

#if (U_SPACE == Q2C) || (U_SPACE == GQ2C) || (U_SPACE == GQ2B3C) || (U_SPACE == GQ2C_2ELBUB) || (U_SPACE == GQ2C_3ELBUB)
#if (U_SPACE == Q2C) || (U_SPACE == GQ2C)
#define ELEM         q2c_element
#elif U_SPACE == GQ2B3C
#define ELEM         q2b3c_element
#elif U_SPACE == GQ2C_2ELBUB
#define ELEM         q2c2b_element
#elif U_SPACE == GQ2C_3ELBUB
#define ELEM         q2c3b_element
#endif
#define REF_MAP      Q1_REF_MAP
#define L2_Q_RULE    qs9
#define H10_Q_RULE   qs9
#define LAPL_Q_RULE  qs9
#define CONV_Q_RULE  qs9
#define REACT_Q_RULE qs9
#define SD_Q_RULE    qs9
#define SC_Q_RULE    qs9
#define RHS_Q_RULE   qs9
#define SDRHS_Q_RULE qs9
#endif

#define IS_BN(n)        (NTYPE(n) & 1)        /* boundary node                */
#define NOT_BN(n)       (!(NTYPE(n) & 1))
#define IS_FN(n)        (!(NTYPE(n) & NMASK)) /* 3 dofs for vel. in n         */
#define NOT_FN(n)       (NTYPE(n) & NMASK)
#define IS_DN(n)        (NTYPE(n) & NMASK)    /* Dirichlet b.c. for vel. in n */
#define IS_FBN(n)       (NTYPE(n) & NMASK_FOR_FB)       /* free boundary node */
#define NOT_FBN(n)      (!(NTYPE(n) & NMASK_FOR_FB))
#define IS_ZNN(n)       ((NTYPE(n) & NMASK_ZN) && IS_FN(n)) /* zero normal 
                                                               condition in n */
#define NOT_ZNN(n)      (!IS_ZNN(n))
#define IS_ZNF(f)       (FTYPE(f) & FMASK_ZN)   /* zero normal condition in f */
#define NOT_ZNF(f)      (!(FTYPE(f) & FMASK_ZN)) 
#define INCR(n,m)       ((NOT_FN(n) || IS_FN(m)) &&                            \
                        ((NOT_FN(n) && IS_FN(m)) || INDEX(n) < INDEX(m)))
#define IS_IN_GLG(n)    (  NTYPE(n) & GLG &&                                   \
                        !(NTYPE(n) & GW || NTYPE(n) & GLS) )

#define FTOPLEVEL(f)         ((INT)(TOP_LEVEL_BITS & (f)->type))
#define SET_FTOPLEVEL(f,l)   FTYPE(f) = (FTYPE(f) & ~TOP_LEVEL_BITS) | (l)
#define SET_ORIENT(f,o)      FTYPE(f) = (FTYPE(f) & ~ORIENT) | (o) | OR_SET
#define IS_BF(f)         (FTYPE(f) & BNDRY_FACE_BIT)  /* boundary face        */
#define NOT_BF(f)        (!(FTYPE(f) & BNDRY_FACE_BIT))
#define IS_FF(f)         (!(FTYPE(f) & FMASK))        /* dofs for vel. in f   */
#define NOT_FF(f)        (FTYPE(f) & FMASK)        /* NOT_FF(f) -> IS_BF(f)   */
#define IFLAG(p)         ((p)->iflag)
#define ITYPE(p)         ((p)->itype)
#define IS_LN(n)        ((n)->s_node)
#define IS_LF(f)        ((f)->s_face)

#define MARK_BND_EDGE(l)         (l)->flag |= BND_EDGE
#define MARK_EDGE_TO_MIDDLE(l)   (l)->flag |= EDGE_TO_MIDDLE /* l is a link 
                                                     to the middle of an edge */
#define MARK_EDGE_TO_OLD(l)      (l)->flag |= EDGE_TO_OLD /* l is a link to 
                                              a neighbour from an older level */
#define MARK_EDGE_TO_DIVIDE(l)   (l)->flag |= EDGE_TO_DIVIDE
#define MARK_DIVIDED_EDGE(l)     (l)->flag  = ((l)->flag & ~EDGE_TO_DIVIDE) |  \
                                              DIVIDED_EDGE
#define MARK_TREATED_EDGE(l)     (l)->flag |= TREATED_EDGE
#define MARK_TYPE_OF_EDGE(l,r)   (l)->flag |= r

#define IS_BND_EDGE(l)           ((l)->flag & BND_EDGE)
#define IS_EDGE_TO_MIDDLE(l)     ((l)->flag & EDGE_TO_MIDDLE)
#define IS_EDGE_TO_OLD(l)        ((l)->flag & EDGE_TO_OLD)
#define IS_EDGE_TO_DIVIDE(l)     ((l)->flag & EDGE_TO_DIVIDE)
#define IS_DIVIDED_EDGE(l)       ((l)->flag & DIVIDED_EDGE)
#define IS_TREATED_EDGE(l)       ((l)->flag & TREATED_EDGE)

#define PREV(p)         (p)->prev
#define SUCC(p) 	(p)->succ
#define INDEX(p) 	(p)->index
#define START(p) 	(p)->start
#define TSTART(p) 	(p)->tstart
#define LGSTART(p)      (p)->lgstart
#define NFSTART(p) 	(p)->nfstart
#define TNFSTART(p) 	(p)->tnfstart
#define FSTART(p) 	(p)->fstart
#define TFSTART(p) 	(p)->tfstart
#define FNSTART(p) 	(p)->fnstart
#define TFNSTART(p) 	(p)->tfnstart
#define NESTART(p) 	(p)->nestart
#define MYVERTEX(p)     (p)->myvertex

#define ND(p,i,j) ((p)->vdata)[(i)][(j)]
#define NDD(p,i)  ((p)->vdata)[(i)]
#define NDMV(p,i,j)  ((p)->mvdata)[(i)][(j)]
#define NDMVP(p,i)   ((p)->mvdata)[(i)]
#define NDLG(p,i,j)  ((p)->lgd->tdata)[(i)][(j)]
#define NDLGP(p,i,j) ((p)->lgd->tdata)[(i)]
#define NDS(p,i)     ((p)->sdata)[(i)]
#define SNDS(p,i)    ((p)->ssdata)[(i)]
#define FD(p,i)      ((p)->sdata)[(i)]
#define FDV(p,i,j)   ((p)->vdata)[(i)][(j)]
#define FDVP(p,i)    ((p)->vdata)[(i)]
#define FDDV(p,i,j,k)  ((p)->dvdata)[(i)][(j)][(k)]
#define FDDVP(p,i,j)   ((p)->dvdata)[(i)][(j)]
#define FDDVPP(p,i)    ((p)->dvdata)[(i)]
#define FDMV(p,i,j)    ((p)->mvdata)[(i)][(j)]
#define FDMVP(p,i)     ((p)->mvdata)[(i)]
#define SFDS(p,i)      ((p)->ssdata)[(i)]
#define ED(p,i)        ((p)->sdata)[(i)]
#define EDV(p,i,j)     ((p)->vdata)[(i)][(j)]
#define EDVP(p,i)      ((p)->vdata)[(i)]
#define EDMV(p,i,j)    ((p)->mvdata)[(i)][(j)]
#define EDMVP(p,i)     ((p)->mvdata)[(i)]
#define EDSN(p,i,j)    ((p)->sndata)[(i)][(j)]
#define EDSNP(p,i)     ((p)->sndata)[(i)]
#define COEFFN(p,i)        ((p)->diag)[(i)]
#define COEFF_FF(p,i)      ((p)->diag)[(i)]
#define COEFFL(p,i)        ((p)->diag)[(i)]
#define COEFF_FL(p,i)      ((p)->diag)[(i)]
#define COEFFS(p,i)        ((p)->diag)[(i)]
#define COEFFLS(p,i)       ((p)->diag)[(i)]
#define COEFF_NF(p,i,j)    ((p)->vdiag)[(i)][(j)]
#define COEFF_FN(p,i,j)    ((p)->vdiag)[(i)][(j)]
#define COEFF_NFP(p,i)     ((p)->vdiag)[(i)]
#define COEFF_FNP(p,i)     ((p)->vdiag)[(i)]
#define COEFF_EE(p,i)      ((p)->ee_diag)[(i)]
#define COEFFB(p,i,j)      ((p)->bdiag)[(i)][(j)]
#define COEFFBL(p,i,j)     ((p)->bdiag)[(i)][(j)]
#define COEFFBP(p,i)       ((p)->bdiag)[(i)]
#define COEFFBLP(p,i)      ((p)->bdiag)[(i)]
#define COEFFBD(p,i,j,k)   ((p)->bddiag)[(i)][(j)][(k)]
#define COEFFBDP(p,i,j)    ((p)->bddiag)[(i)][(j)]
#define COEFFBDPP(p,i)     ((p)->bddiag)[(i)]
#define COEFFLL(p,k,i,j)   ((p)->a)[(k)][(i)][(j)]
#define COEFFNN(p,k,i,j)   ((p)->a)[(k)][(i)][(j)]
#define COEFFLLP(p,k)      ((p)->a)[(k)]
#define COEFFNNP(p,k)      ((p)->a)[(k)]
#define COEFFNNPP(p,k,i)   ((p)->a)[(k)][(i)]
#define COEFF_NN(p,k,i,j)  ((p)->ann)[(k)][(i)][(j)]
#define COEFF_NNP(p,k)     ((p)->ann)[(k)]
#define COEFF_EE_MM(p,k,i,j)    ((p)->ee_mm)[(k)][(i)][(j)]
#define COEFF_EN_MN(p,k,l,i,j)  ((p)->en_mn)[(k)][(l)][(i)][(j)]
#define COEFF_NE_NM(p,k,l,i,j)  ((p)->ne_nm)[(k)][(l)][(i)][(j)]
#define COEFF_EF_MN(p,k,l,i,j)  ((p)->ef_mn)[(k)][(l)][(i)][(j)]
#define COEFF_FE_NM(p,k,l,i,j)  ((p)->fe_nm)[(k)][(l)][(i)][(j)]
#define COEFF_EE_MMP(p,k)       ((p)->ee_mm)[(k)]
#define COEFF_EN_MNP(p,k,l)     ((p)->en_mn)[(k)][(l)]
#define COEFF_NE_NMP(p,k,l)     ((p)->ne_nm)[(k)][(l)]
#define COEFF_EF_MNP(p,k,l)     ((p)->ef_mn)[(k)][(l)]
#define COEFF_FE_NMP(p,k,l)     ((p)->fe_nm)[(k)][(l)]
#define COEFF_NLG(p,k,i,j) ((p)->a)[(k)][(i)][(j)]
#define COEFF_FLG(p,k,j)   ((p)->a)[(k)][(j)]
#define COEFF_LGN(p,k,i,j) ((p)->a)[(k)][(i)][(j)]
#define COEFF_LGF(p,k,j)   ((p)->a)[(k)][(j)]
#define COEFF_LGL(p,k,i,j) ((p)->a)[(k)][(i)][(j)]
#define COEFF_LG(p,k,i,j)  ((p)->lgd->a)[(k)][(i)][(j)]
#define COEFF_NLGP(p,k)    ((p)->a)[(k)]
#define COEFF_FLGP(p,k)    ((p)->a)[(k)]
#define COEFF_LGNP(p,k)    ((p)->a)[(k)]
#define COEFF_LGFP(p,k)    ((p)->a)[(k)]
#define COEFF_LGLP(p,k)    ((p)->a)[(k)]
#define COEFF_LGP(p,k)     ((p)->lgd->a)[(k)]
#define COEFF_BN(p,i,n,j)    ((p)->bn)[(i)][(n)][(j)]
#define COEFF_BNP(p,i,n)     ((p)->bn)[(i)][(n)]
#define COEFF_BF(p,i,j)      ((p)->ef_diag)[(i)][(j)]
#define COEFF_EN(p,i,j)      ((p)->en_diag)[(i)][(j)]
#define COEFF_NE(p,i,j)      ((p)->ne_diag)[(i)][(j)]
#define COEFF_EF(p,i,j)      ((p)->ef_diag)[(i)][(j)]
#define COEFF_FE(p,i,j)      ((p)->fe_diag)[(i)][(j)]
#define COEFF_ENP(p,i)       ((p)->en_diag)[(i)]
#define COEFF_NEP(p,i)       ((p)->ne_diag)[(i)]
#define COEFF_EFP(p,i)       ((p)->ef_diag)[(i)]
#define COEFF_FEP(p,i)       ((p)->fe_diag)[(i)]
#define COEFF_BDF(p,i,j,k)   ((p)->bdf)[(i)][(j)][(k)]
#define COEFF_BDFP(p,i,j)    ((p)->bdf)[(i)][(j)]
#define COEFF_BFDN(p,i,j,k,l)    ((p)->bfdn)[(i)][(j)][(k)][(l)]
#define COEFF_BFDNP(p,i,j,k)     ((p)->bfdn)[(i)][(j)][(k)]
#define COEFF_BFDF(p,i,j,k,l)    ((p)->bfdf)[(i)][(j)][(k)][(l)]
#define COEFF_BFDFP(p,i,j,k)     ((p)->bfdf)[(i)][(j)][(k)]
#define COEFF_BDFS(p,i,j,k,l,m)  ((p)->bdfs)[(i)][(j)][(k)][(l)][(m)]
#define COEFF_BDFSP(p,i,j,k,l)   ((p)->bdfs)[(i)][(j)][(k)][(l)]

#define NBNODE(p)	(p)->nbnode
#define NBFACE(p)	(p)->nbface
#define NBELEM(p)	(p)->nbel
#define NB_EL(p,i)	((p)->nbel)[(i)]
#define NEXT(p)		(p)->next

#define FIRSTN(p)       (p)->firstN
#define FIRSTF(p)       (p)->firstF
#define FIRSTNODE(p)    (p)->firstNode
#define FIRSTFACE(p)    (p)->firstFace
#define FIRSTSN(p)      (p)->firstSN
#define FIRSTSF(p)      (p)->firstSF
#define FIRSTPN(p)      (p)->firstPN
#define FIRSTPF(p)      (p)->firstPF
#define FIRSTELEMENT(p) (p)->firstElement
#define FIRSTGRID(p)    (p)->grids[1]
#define FDBN(p)         (p)->fdbn
#define FDBF(p)         (p)->fdbf
#define FDBE(p)         (p)->fdbe
#define LASTNODE(p)     (p)->lastNode
#define LASTFACE(p)     (p)->lastFace
#define IS_CN(p,t)      ((t)&STOP_IS_FIRST_INNER ? NOT_FN(p) : IS_FN(p))
#define IS_CF(p,t)      ((t)&STOP_IS_FIRST_INNER ? NOT_FF(p) : IS_FF(p))
#define FIRST_FACE(p,t) ((t)&ONLY_INNER  ? FIRSTFACE(p) : FIRSTF(p))
#define FIRST_NODE(p,t) ((t)&ONLY_INNER  ? FIRSTNODE(p) : FIRSTN(p))
#define STOP_FACE(p,t)  ((t)&STOP_IS_FIRST_INNER ? FIRSTFACE(p) : NULL)
#define STOP_NODE(p,t)  ((t)&STOP_IS_FIRST_INNER ? FIRSTNODE(p) : NULL)
#define STOP_NN(p,t)    ((t)&STOP_IS_FIRST_INNER ? START(p) : NULL)
#define STOP_NF(p,t)    ((t)&STOP_IS_FIRST_INNER ? NFSTART(p) : NULL)
#define STOP_FF(p,t)    ((t)&STOP_IS_FIRST_INNER ? FSTART(p) : NULL)
#define STOP_FN(p,t)    ((t)&STOP_IS_FIRST_INNER ? FNSTART(p) : NULL)
#define T_START(p,t)    ((t)&  NSTART_FROM_INNER ? START(p)   : TSTART(p))
#define NF_START(p,t) 	((t)& NFSTART_FROM_INNER ? NFSTART(p) : TNFSTART(p))
#define F_START(p,t) 	((t)&  FSTART_FROM_INNER ? FSTART(p)  : TFSTART(p))
#define FN_START(p,t) 	((t)& FNSTART_FROM_INNER ? FNSTART(p) : TFNSTART(p))

#define TOP_NODE(n)     (n)->myvertex->topnode
#define TOP_GRID(mg)    (mg->grids)[mg->toplevel]
#define IS_TOP_NODE(n)    ((n)->son == NULL) 
#define IS_TOP_FACE(f)    ((f)->sons[0] == NULL) 
#define IS_TOP_ELEMENT(e) ((e)->sons[0] == NULL)
#define IS_LTOP_NODE(n,g) ((n)->son == NULL ||                                 \
                                           INDEX((n)->son) > (g)->maxNodeIndex) 
#define IS_LTOP_FACE(f,g)    (FTOPLEVEL(f) >= (INT)(g)->level)
#define IS_LTOP_ELEMENT(e,g) ((e)->sons[0] == NULL ||                          \
                               INDEX(((e)->sons[0])->n[0]) > (g)->maxNodeIndex)
#define FIND_LTOP_FACE(f,pfl,g)                                                \
               f = NBFACE(pfl);                                                \
               if (f->sons[0] != f)                                            \
                  while (f->sons[0] && INDEX(f->sons[0]) <= (g)->maxFaceIndex  \
                                    && INDEX(f->sons[0]) > 0) f = f->sons[0];
#define ON_FINEST_GRID(n,mg) (n->index >=                                      \
                                       (mg->grids)[mg->toplevel]->minNodeIndex)
#define NO_NODE_ON_NEW_LEVEL(n,mg) (TOP_NODE(n)->index <= mg->nodeNumber)
#define NODE_ON_NEW_LEVEL(n,mg) ((n)->index > mg->nodeNumber)
#define FACE_ON_NEW_LEVEL(f,mg) ((f)->index > mg->faceNumber)
#define F_ON_FINEST_GRID(f,mg) (f->index >=                                    \
                                       (mg->grids)[mg->toplevel]->minFaceIndex)
#define NOT_ON_NEW_LEVEL(f,mg) ((f)->index > 0 && (f)->index <= mg->faceNumber)
#define MAX(a,b)          ((a) > (b) ? (a) : (b))
#define MIN(a,b)          ((a) < (b) ? (a) : (b))
#define SGN(a)            ((a) > 0 ? 1 : ( (a) < 0 ? -1 : 0 ))
#define TMIN(a,b)         ((a) && (b) ? MIN(a,b) : MAX(a,b))
#define INDI(a,b)         ((a) > (b) ? (1.) : (-1.))
#define NINDI(a,b)        ((a)->index > (b)->index ? (1.) : (-1.))
#define IINDI(a,b)        ((a)->index > (b)->index ? 1 : 0)
#define EXCHANGE(a,b,z)   {z = a; a = b; b = z;}
#define SAME_COUPLES(i1,i2,j1,j2) (((i1==j1)&&(i2==j2))||((i1==j2)&&(i2==j1)))
#if !(PROBLEM & WITH_F)
#define RHS_FOR_FACES    ;
#define RHS_FOR_FACESN   ;
#endif
#define SET2x2(a,a00,a01,a10,a11)  { a[0][0]=a00; a[0][1]=a01;                 \
                                     a[1][0]=a10; a[1][1]=a11; }
#define SET3x3(a,a00,a01,a02,a10,a11,a12,a20,a21,a22)                          \
                                 { a[0][0]=a00; a[0][1]=a01; a[0][2]=a02;      \
                                   a[1][0]=a10; a[1][1]=a11; a[1][2]=a12;      \
                                   a[2][0]=a20; a[2][1]=a21; a[2][2]=a22; } 
#define MFILL_2x2(a,b,i,j,in,jn) { a[i ][j] = b[0][0]; a[i ][jn] = b[0][1];    \
                                   a[in][j] = b[1][0]; a[in][jn] = b[1][1]; }
#define MFILL_SUBTR_2x2(a0,a1,b,i,j,z,p)                                   \
                                { a0[i][j] = b[0][0];  a1[i][j] = b[1][1];     \
                                  z[0] -= b[0][1]*p[1]; z[1] -= b[1][0]*p[0]; }
#define polynom4(a0,a1,a2,a3,a4,x) (a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4)))))
#define POINT_VALUES_3(x1,x2,x3,v1,v2,v3,v)  {  v1=v(x1); v2=v(x2); v3=v(x3); }
#define POINT_VALUES_6(x1,x2,x3,x4,x5,x6,v1,v2,v3,v4,v5,v6,v)                  \
              {  v1=v(x1); v2=v(x2); v3=v(x3); v4=v(x4); v5=v(x5); v6=v(x6);  }
#define POINT_VALUES_7(x1,x2,x3,x4,x5,x6,x7,v1,v2,v3,v4,v5,v6,v7,v)            \
     {  v1=v(x1); v2=v(x2); v3=v(x3); v4=v(x4); v5=v(x5); v6=v(x6); v7=v(x7);  }
#define POINT_VALUES_9(x1,x2,x3,x4,x5,x6,x7,x8,x9,                             \
                       v1,v2,v3,v4,v5,v6,v7,v8,v9,v)                           \
              {  v1=v(x1); v2=v(x2); v3=v(x3); v4=v(x4); v5=v(x5); v6=v(x6);   \
                 v7=v(x7); v8=v(x8); v9=v(x9); }
#define POINT_VALUES_10(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,                        \
                        v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v)                      \
              {  v1=v(x1); v2=v(x2); v3=v(x3); v4=v(x4); v5=v(x5); v6=v(x6);   \
                 v7=v(x7); v8=v(x8); v9=v(x9); v10=v(x10); }
#define QUADR_VALUES_10(u0,u1,u2,u01,u02,u12,                                  \
                        f0,f1,f2,f001,f002,f110,f112,f220,f221,f012)           \
              {  f0=u0; f1=u1; f2=u2; f001=((u0+u0+u1)*3.+u01+u01)/9.;         \
     f002=((u0+u0+u2)*3.+u02+u02)/9.; f110=((u1+u1+u0)*3.+u01+u01)/9.;         \
     f112=((u1+u1+u2)*3.+u12+u12)/9.; f220=((u2+u2+u0)*3.+u02+u02)/9.;         \
     f221=((u2+u2+u1)*3.+u12+u12)/9.; f012=((u0+u1+u2)*3.+u01+u02+u12)/9.;  }
#define SUBTR_CONST_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,c)                       \
                     { v1 -= c; v2 -= c; v3 -= c; v4 -= c; v5 -= c;            \
                       v6 -= c; v7 -= c; v8 -= c; v9 -= c; v10 -= c; }
#define SUBTR_10_FROM_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,                       \
                         x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)                       \
                     {  v1-=x1; v2-=x2; v3-=x3; v4-=x4; v5-=x5;                \
                        v6-=x6; v7-=x7; v8-=x8; v9-=x9; v10-=x10; }
#define ADD_10_TO_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,                           \
                     x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)                           \
                     {  v1+=x1; v2+=x2; v3+=x3; v4+=x4; v5+=x5;                \
                        v6+=x6; v7+=x7; v8+=x8; v9+=x9; v10+=x10; }
#define MULT_AND_ADD_10_TO_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,                  \
                              x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,r)                \
                     {  v1+=x1*r; v2+=x2*r; v3+=x3*r; v4+=x4*r; v5+=x5*r;      \
                        v6+=x6*r; v7+=x7*r; v8+=x8*r; v9+=x9*r; v10+=x10*r; }
#define SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)                                   \
                   {n11 = n1->son; n22 = n2->son; n33 = n3->son;}
#define SONS_OF_NODES4(n1,n2,n3,n4,n11,n22,n33,n44)                            \
                   {n11 = n1->son; n22 = n2->son; n33 = n3->son; n44 = n4->son;}
#define AREA_OF_TRIANGLE(x0,x1,x2)                                             \
             (0.5*fabs((x0[0]-x2[0])*(x1[1]-x2[1])-(x1[0]-x2[0])*(x0[1]-x2[1])))

#if DIM == 3
#define DIM1            2        /* DIM - 1                           */
#define DIM2            4        /* DIM + 1                           */
#define DIM3            6        /* 2*DIM                             */
#define DIM4            7        /* 2*DIM + 1                         */
#if ELEMENT_TYPE == SIMPLEX
#define NVERT		4	 /* max number of vertices of an element */
#define SIDES		4	 /* max number of sides of an element */
#define EDGES		6	 /* max number of edges of an element */
#define IS_FTYPE(e,mask) (FTYPE(e->f[0]) & mask || FTYPE(e->f[1]) & mask ||    \
                          FTYPE(e->f[2]) & mask || FTYPE(e->f[3]) & mask)
#define VOLUME(e)       volume(e->n[0]->myvertex->x,e->n[1]->myvertex->x,      \
                               e->n[2]->myvertex->x,e->n[3]->myvertex->x)
#define IS_B_EL(e)      (IS_DN(e->n[0]) || IS_DN(e->n[1]) ||                   \
                         IS_DN(e->n[2]) || IS_DN(e->n[3]))
#elif ELEMENT_TYPE == CUBE
#define NVERT		8	 /* max number of vertices of an element */
#define SIDES		6	 /* max number of sides of an element */
#define EDGES	       12        /* max number of edges of an element */
#define IS_FTYPE(e,mask) (FTYPE(e->f[0]) & mask || FTYPE(e->f[1]) & mask ||    \
                          FTYPE(e->f[2]) & mask || FTYPE(e->f[3]) & mask)||    \
                          FTYPE(e->f[4]) & mask || FTYPE(e->f[5]) & mask)
#define VOLUME(e)        (-1.)
#define IS_B_EL(e)      (IS_DN(e->n[0]) || IS_DN(e->n[1]) || IS_DN(e->n[2]) || \
                         IS_DN(e->n[3]) || IS_DN(e->n[4]) || IS_DN(e->n[5]) || \
                         IS_DN(e->n[6]) || IS_DN(e->n[7]))
#endif
#define SONS		8        /* max number of sons of an element  */
#define FSONS		4        /* max number of sons of a face      */
#define AVER_N(e,p,i)   (ND(e->n[0],p,i) + ND(e->n[1],p,i) +                   \
                         ND(e->n[2],p,i) + ND(e->n[3],p,i))/4.
#define AVER_NS(e,p)    (NDS(e->n[0],p) + NDS(e->n[1],p) +                     \
                         NDS(e->n[2],p) + NDS(e->n[3],p))/4.
#define AVER_FS(e,p)    (FD(e->f[0],p) + FD(e->f[1],p) +                       \
                         FD(e->f[2],p) + FD(e->f[3],p))/4.
#define AVER_ESN(e,p)   (EDSN(e,p,0) +EDSN(e,p,1) +EDSN(e,p,2) +EDSN(e,p,3))/4.
#define VERTICES_OF_ELEMENT(x1,x2,x3,x4,e) x1 = e->n[0]->myvertex->x;          \
                                           x2 = e->n[1]->myvertex->x;          \
                                           x3 = e->n[2]->myvertex->x;          \
                                           x4 = e->n[3]->myvertex->x
#define NODES_OF_ELEMENT(n1,n2,n3,n4,e)    n1=e->n[0]; n2=e->n[1];             \
                                           n3=e->n[2]; n4=e->n[3]
#define INODES_OF_ELEMENT(n1,n2,n3,n4,i1,i2,i3,i4,e)                           \
                                           n1=e->n[i1]; n2=e->n[i2];           \
                                           n3=e->n[i3]; n4=e->n[i4]
#define FACES_OF_ELEMENT(f1,f2,f3,f4,e)    f1=e->f[0]; f2=e->f[1];             \
                                           f3=e->f[2]; f4=e->f[3]
#define TOPNODES_OF_ELEMENT(n1,n2,n3,n4,e) n1=TOP_NODE(e->n[0]);               \
              n2=TOP_NODE(e->n[1]); n3=TOP_NODE(e->n[2]); n4=TOP_NODE(e->n[3])
#define TOPFACES_OF_ELEMENT(f1,f2,f3,f4,e) f1=top_face(e->f[0]);               \
              f2=top_face(e->f[1]); f3=top_face(e->f[2]); f4=top_face(e->f[3])
#define LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,e,g) n1=ltop_node(e->n[0],g);         \
     n2=ltop_node(e->n[1],g); n3=ltop_node(e->n[2],g); n4=ltop_node(e->n[3],g)
#define LTOPFACES_OF_ELEMENT(f1,f2,f3,f4,e,g) f1=ltop_face(e->f[0],g);         \
     f2=ltop_face(e->f[1],g); f3=ltop_face(e->f[2],g); f4=ltop_face(e->f[3],g)
#define S_NODE_VALUES(n0,n1,n2,n3,f,g0,g1,g2,g3) { g0 = f(n0->myvertex->x);    \
                                                   g1 = f(n1->myvertex->x);    \
                                                   g2 = f(n2->myvertex->x);    \
                                                   g3 = f(n3->myvertex->x); }
#define V_NODE_VALUES(n0,n1,n2,n3,bb0,bb1,bb2,bb_0,bb_1,bb_2,bb_3)             \
         { bb_0[0] = bb0(n0->myvertex->x);   bb_0[1] = bb1(n0->myvertex->x);   \
                                             bb_0[2] = bb2(n0->myvertex->x);   \
           bb_1[0] = bb0(n1->myvertex->x);   bb_1[1] = bb1(n1->myvertex->x);   \
                                             bb_1[2] = bb2(n1->myvertex->x);   \
           bb_2[0] = bb0(n2->myvertex->x);   bb_2[1] = bb1(n2->myvertex->x);   \
                                             bb_2[2] = bb2(n2->myvertex->x);   \
           bb_3[0] = bb0(n3->myvertex->x);   bb_3[1] = bb1(n3->myvertex->x);   \
                                             bb_3[2] = bb2(n3->myvertex->x); }
#define NORMAL_VECTORS(nn1,nn2,nn3,nn4,e)                                      \
                    normal_vector(e->n[1]->myvertex->x,e->n[2]->myvertex->x,   \
                                  e->n[3]->myvertex->x,e->f[0],nn1);           \
                    normal_vector(e->n[0]->myvertex->x,e->n[2]->myvertex->x,   \
                                  e->n[3]->myvertex->x,e->f[1],nn2);           \
                    normal_vector(e->n[0]->myvertex->x,e->n[1]->myvertex->x,   \
                                  e->n[3]->myvertex->x,e->f[2],nn3);           \
                    normal_vector(e->n[0]->myvertex->x,e->n[1]->myvertex->x,   \
                                  e->n[2]->myvertex->x,e->f[3],nn4)
#define ARE_NEIGHBOURS(p,q)  (p->n[0]==q->n[0] || p->n[0]==q->n[1] ||          \
                              p->n[0]==q->n[2] || p->n[0]==q->n[3] ||          \
                              p->n[1]==q->n[0] || p->n[1]==q->n[1] ||          \
                              p->n[1]==q->n[2] || p->n[1]==q->n[3] ||          \
                              p->n[2]==q->n[0] || p->n[2]==q->n[1] ||          \
                              p->n[2]==q->n[2] || p->n[2]==q->n[3] ||          \
                              p->n[3]==q->n[0] || p->n[3]==q->n[1] ||          \
                              p->n[3]==q->n[2] || p->n[3]==q->n[3])
#define ARE_FNEIGHBOURS(p,q) (p->f[0]==q->f[0] || p->f[0]==q->f[1] ||          \
                              p->f[0]==q->f[2] || p->f[0]==q->f[3] ||          \
                              p->f[1]==q->f[0] || p->f[1]==q->f[1] ||          \
                              p->f[1]==q->f[2] || p->f[1]==q->f[3] ||          \
                              p->f[2]==q->f[0] || p->f[2]==q->f[1] ||          \
                              p->f[2]==q->f[2] || p->f[2]==q->f[3] ||          \
                              p->f[3]==q->f[0] || p->f[3]==q->f[1] ||          \
                              p->f[3]==q->f[2] || p->f[3]==q->f[3])
#define FACE_CONTAINED(fa,p) (fa==p->f[0] || fa==p->f[1] ||                    \
                              fa==p->f[2] || fa==p->f[3])
#define GIVE_INDEX(p,n,i)  { if (p[0] == n) i = 0;                             \
                             else if (p[1] == n) i = 1;                        \
                             else if (p[2] == n) i = 2;                        \
                             else if (p[3] == n) i = 3;                        \
                             else i = -1; }
#define IS_IN_SQUARE(e,xm,ym)  ( e->n[0]->myvertex->x[0] <= xm &&              \
                                 e->n[1]->myvertex->x[0] <= xm &&              \
                                 e->n[2]->myvertex->x[0] <= xm &&              \
                                 e->n[3]->myvertex->x[0] <= xm &&              \
                                 e->n[0]->myvertex->x[1] <= ym &&              \
                                 e->n[1]->myvertex->x[1] <= ym &&              \
                                 e->n[2]->myvertex->x[1] <= ym &&              \
                                 e->n[3]->myvertex->x[1] <= ym )
#define IS_SON_FACE(f,g)       ( f == g->sons[0] || f == g->sons[1] ||         \
                                 f == g->sons[2] || f == g->sons[3] )
#if PROBLEM & WITH_F
#define RHS_FOR_FACES    FD(fa0,f) +=  integr3(n0,n1,n2,n3,nn0,rdetB)*FMULT;
#define RHS_FOR_FACESN   FD(fa0,f) += integr3c(n0,n1,n2,n3,nn0,rdetB)*FMULT;
#endif
#define IS_IN_SIMPLEX(b,x,eps)    (DOT(b[0],x)+b[0][3] >= eps &&               \
                                   DOT(b[1],x)+b[1][3] >= eps &&               \
                                   DOT(b[2],x)+b[2][3] >= eps &&               \
                                   DOT(b[3],x)+b[3][3] >= eps)
#define IS_IN(n,n1,n2,n3,n4) ( n == n1 || n == n2 || n == n3 || n == n4 )
#define NOT_ALL_FN(p)    (NOT_FN(p->n[0]) || NOT_FN(p->n[1]) ||                \
                          NOT_FN(p->n[2]) || NOT_FN(p->n[3]))
#define SUBTR(a,b,c)      c[0]=a[0]-b[0]; c[1]=a[1]-b[1]; c[2]=a[2]-b[2]
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define DOTS(a,b)        (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])
#define DOT2(a,b,c)      (a[0]*(b[0]+4.*c[0]) + a[1]*(b[1]+4.*c[1]) +          \
                                                a[2]*(b[2]+4.*c[2]))
#define DOT_DIFF(a,b,c,d)  ( (a[0]-b[0])*(c[0]-d[0]) + (a[1]-b[1])*(c[1]-d[1]) \
                           + (a[2]-b[2])*(c[2]-d[2]) )
#define DISTANCE(a,b)    (sqrt((a[0]-b[0])*(a[0]-b[0])+                        \
                               (a[1]-b[1])*(a[1]-b[1])+                        \
                               (a[2]-b[2])*(a[2]-b[2])))
#define SET1(a,b)        {a[0]=b[0]; a[1]=b[1]; a[2]=b[2];}
#define SSET1(a,b)       {a[0]=b[0]; a[1]=b[1]; a[2]=b[2]; a[3]=b[3];}
#define SET2(a,b,p)      {a[0]=b[0]*(p); a[1]=b[1]*(p); a[2]=b[2]*(p);}
#define SET3(a,b,p,q,r)  {r=(p)*(q); a[0]=b[0]*r; a[1]=b[1]*r; a[2]=b[2]*r;}
#define SET4(a,b,p)      {a[0]+=b[0]*(p); a[1]+=b[1]*(p); a[2]+=b[2]*(p);}
#define SET5(a,b)        {a[0]+=b[0];   a[1]+=b[1];   a[2]+=b[2];  }
#define SET6(a,b)        {a[0]-=b[0];   a[1]-=b[1];   a[2]-=b[2];  }
#define SET7(a,p)        {a[0]=p; a[1]=p; a[2]=p;}
#define SSET7(a,p)       {a[0]=p; a[1]=p; a[2]=p; a[3]=p;}
#define SET8(a,b)        {a[0]=-b[0]; a[1]=-b[1]; a[2]=-b[2];}
#define SET9(a,b,p)      {a[0]-=b[0]*(p); a[1]-=b[1]*(p); a[2]-=b[2]*(p);}
#define SET10(a,b,c)     {a[0]=b[0]+c[0]; a[1]=b[1]+c[1]; a[2]=b[2]+c[2];}
#define SET11(a,b,c)     {a[0]=b[0]-c[0]; a[1]=b[1]-c[1]; a[2]=b[2]-c[2];}
#define SET12(a,b,c,p)   {a[0]=(b[0]-c[0])/(p); a[1]=(b[1]-c[1])/(p);          \
                          a[2]=(b[2]-c[2])/(p);}
#define SET13(a,p)       {a[0]/=(p); a[1]/=(p); a[2]/=(p);}
#define SET14(z,a,b,c)   {z[0]=a[0]+b[0]+c[0]; z[1]=a[1]+b[1]+c[1];            \
                                               z[2]=a[2]+b[2]+c[2];}
#define SET15(z,a,b,c,d) {z[0]=a[0]+b[0]+c[0]+d[0]; z[1]=a[1]+b[1]+c[1]+d[1];  \
                                                    z[2]=a[2]+b[2]+c[2]+d[2];}
#define SET16(z,a,b,c,d,e) {z[0]=a[0]+b[0]+c[0]+d[0]+e[0];                     \
                            z[1]=a[1]+b[1]+c[1]+d[1]+e[1];                     \
                            z[2]=a[2]+b[2]+c[2]+d[2]+e[2];}
#define SET17(a,b,c)     {a[0]=b[0]+b[0]+c[0]; a[1]=b[1]+b[1]+c[1];            \
                                               a[2]=b[2]+b[2]+c[2];}
#define SET18(a,b,c)     {a[0]=b[0]+b[0]-c[0]; a[1]=b[1]+b[1]-c[1];            \
                                               a[2]=b[2]+b[2]-c[2];}
#define SET19(a,b,c)     {a[0]=b[0]+b[0]+b[0]-c[0]-c[0];                       \
                          a[1]=b[1]+b[1]+b[1]-c[1]-c[1];                       \
                          a[2]=b[2]+b[2]+b[2]-c[2]-c[2];}
#define SET20(z,a,p,b,q) {z[0]=a[0]*(p) + b[0]*(q);  z[1]=a[1]*(p) + b[1]*(q); \
                          z[2]=a[2]*(p) + b[2]*(q);}
#define SET21(z,a,p,b,q) {z[0] = a[0]/(p) + b[0]/(q);                          \
                          z[1] = a[1]/(p) + b[1]/(q);                          \
                          z[2] = a[2]/(p) + b[2]/(q);}
#define SET22(a,b,p)     {a[0]=b[0]/(p); a[1]=b[1]/(p); a[2]=b[2]/(p);}
#define SET23(a,b,c,p)   {a[0]=b[0]+c[0]*(p); a[1]=b[1]+c[1]*(p);              \
                          a[2]=b[2]+c[2]*(p);}
#define SET24(a,b,c,d)   {a[0]=2.*(b[0]+b[0]-c[0]-d[0]);                       \
                          a[1]=2.*(b[1]+b[1]-c[1]-d[1]);                       \
                          a[2]=2.*(b[2]+b[2]-c[2]-d[2]);}
#define SET25(a,b,c)     {a[0]=fabs(b[0]-c[0]); a[1]=fabs(b[1]-c[1]);          \
                          a[2]=fabs(b[2]-c[2]);}
#define SET26(a,p)       {a[0]+=p; a[1]+=p; a[2]+=p; a[3]+=p;}
#define SET27(z,a,p,b,q,c,r) {z[0]=a[0]*(p) + b[0]*(q) + c[0]*(r);             \
                              z[1]=a[1]*(p) + b[1]*(q) + c[1]*(r);             \
                              z[2]=a[2]*(p) + b[2]*(q) + c[2]*(r);}
#define SET28(a,b,c,p)   {a[0]=(b[0]-c[0])*(p); a[1]=(b[1]-c[1])*(p);          \
                          a[2]=(b[2]-c[2])*(p);}
#define SET29(a,b,c,d,p) {a[0]=b[0]+c[0]+d[0]*(p); a[1]=b[1]+c[1]+d[1]*(p);    \
                                                   a[2]=b[2]+c[2]+d[2]*(p);}
#define DOTLG(a,b)       (a[0]*b[0] + a[1]*b[1])
#define SET1LG(a,b)      {a[0]=b[0]; a[1]=b[1];}
#define SET4LG(a,b,p)    {a[0]+=b[0]*(p); a[1]+=b[1]*(p);}
#define SET7LG(a,p)      {a[0]=p; a[1]=p;}
#define MSET1(a,b)       {a[0][0]=b[0][0]; a[0][1]=b[0][1]; a[0][2]=b[0][2];   \
                          a[1][0]=b[1][0]; a[1][1]=b[1][1]; a[1][2]=b[1][2];   \
                          a[2][0]=b[2][0]; a[2][1]=b[2][1]; a[2][2]=b[2][2];}
#define MSET2(a,b,p)     {a[0] = b[0][0]*p[0] + b[0][1]*p[1] + b[0][2]*p[2];   \
                          a[1] = b[1][0]*p[0] + b[1][1]*p[1] + b[1][2]*p[2];   \
                          a[2] = b[2][0]*p[0] + b[2][1]*p[1] + b[2][2]*p[2];}
#define MSET3(a,b,p)     {a[0] = b[0][0]*p[0]; a[1] = b[1][1]*p[1];            \
                          a[2] = b[2][2]*p[2];}
#define MSET33(a,b,p)    {a[0] = p[0]/b[0][0]; a[1] = p[1]/b[1][1];            \
                          a[2] = p[2]/b[2][2];}
#define MSET4(a,b,p)     {a[0] += b[0][0]*p[0] + b[0][1]*p[1] + b[0][2]*p[2];  \
                          a[1] += b[1][0]*p[0] + b[1][1]*p[1] + b[1][2]*p[2];  \
                          a[2] += b[2][0]*p[0] + b[2][1]*p[1] + b[2][2]*p[2];}
#define MSET44(a,b,p)    {a[0][0] += b[0][0]*(p); a[0][1] += b[0][1]*(p);      \
                                                  a[0][2] += b[0][2]*(p);      \
                          a[1][0] += b[1][0]*(p); a[1][1] += b[1][1]*(p);      \
                                                  a[1][2] += b[1][2]*(p);      \
                          a[2][0] += b[2][0]*(p); a[2][1] += b[2][1]*(p);      \
                                                  a[2][2] += b[2][2]*(p);}
#define MSET5(a,b)       {a[0][0]+=b[0][0]; a[0][1]+=b[0][1]; a[0][2]+=b[0][2];\
                          a[1][0]+=b[1][0]; a[1][1]+=b[1][1]; a[1][2]+=b[1][2];\
                          a[2][0]+=b[2][0]; a[2][1]+=b[2][1]; a[2][2]+=b[2][2];}
#define MSET7(a,p)       {a[0][0]=p; a[0][1]=p; a[0][2]=p;                     \
                          a[1][0]=p; a[1][1]=p; a[1][2]=p;                     \
                          a[2][0]=p; a[2][1]=p; a[2][2]=p;}
#define MSET9(a,b,p)     {a[0] -= b[0][0]*p[0] + b[0][1]*p[1] + b[0][2]*p[2];  \
                          a[1] -= b[1][0]*p[0] + b[1][1]*p[1] + b[1][2]*p[2];  \
                          a[2] -= b[2][0]*p[0] + b[2][1]*p[1] + b[2][2]*p[2];}
#define MSET2LG(a,b,p)   {a[0] = b[0][0]*p[0] + b[0][1]*p[1];                  \
                          a[1] = b[1][0]*p[0] + b[1][1]*p[1];}
#define MSET4LG(a,b,p)   {a[0] += b[0][0]*p[0] + b[0][1]*p[1];                 \
                          a[1] += b[1][0]*p[0] + b[1][1]*p[1];} 
#define MSET4NLG(a,b,p)  {a[0] += b[0][0]*p[0] + b[0][1]*p[1];                 \
                          a[1] += b[1][0]*p[0] + b[1][1]*p[1];                 \
                          a[2] += b[2][0]*p[0] + b[2][1]*p[1];}
#define MSET4LGN(a,b,p)  {a[0] += b[0][0]*p[0] + b[0][1]*p[1] + b[0][2]*p[2];  \
                          a[1] += b[1][0]*p[0] + b[1][1]*p[1] + b[1][2]*p[2];}
#define MSET7LG(a,p)     {a[0][0]=p; a[0][1]=p;                                \
                          a[1][0]=p; a[1][1]=p;}
#define MSET7NLG(a,p)    {a[0][0]=p; a[0][1]=p;                                \
                          a[1][0]=p; a[1][1]=p;                                \
                          a[2][0]=p; a[2][1]=p;}            
#define MSET7LGN(a,p)    {a[0][0]=p; a[0][1]=p; a[0][2]=p;                     \
                          a[1][0]=p; a[1][1]=p; a[1][2]=p;}
#define AVERAGE(x,y,z)   {z[0] = (x[0]+y[0])/2.; z[1] = (x[1]+y[1])/2.;        \
                          z[2] = (x[2]+y[2])/2.;}
#define ADDMULT(a,x,b,y,z) z[0] = (a)*x[0]+(b)*y[0]; z[1] = (a)*x[1]+(b)*y[1]; \
                           z[2] = (a)*x[2]+(b)*y[2]
/*  z=(a*x+y)/c  */
#define POINT2(x,y,z,a,c) z[0] = ((a)*x[0]+y[0])/(c);                          \
                          z[1] = ((a)*x[1]+y[1])/(c);                          \
                          z[2] = ((a)*x[2]+y[2])/(c)
/*  z=(w+x+y)/3.0  */
#define POINT3(w,x,y,z)   z[0]=(w[0]+x[0]+y[0])/3.0; z[1]=(w[1]+x[1]+y[1])/3.0;\
                          z[2]=(w[2]+x[2]+y[2])/3.0
#define POINT4(v,w,x,y,z) z[0]=(v[0]+w[0]+x[0]+y[0])/4.0;                      \
                          z[1]=(v[1]+w[1]+x[1]+y[1])/4.0;                      \
                          z[2]=(v[2]+w[2]+x[2]+y[2])/4.0
#define POINT5(w,x,y,z,a,c)     z[0] = ((a)*w[0]+x[0]+y[0])/(c);               \
                                z[1] = ((a)*w[1]+x[1]+y[1])/(c);               \
                                z[2] = ((a)*w[2]+x[2]+y[2])/(c)  
#define MIDPOINTS(x0,x1,x2,x3,x01,x02,x03,x12,x13,x23)                         \
           {  AVERAGE(x0,x1,x01); AVERAGE(x0,x2,x02); AVERAGE(x0,x3,x03);      \
              AVERAGE(x1,x2,x12); AVERAGE(x1,x3,x13); AVERAGE(x2,x3,x23);  }
#define MAX_DIFF(a,b)  MAX(MAX(fabs(a[0]-b[0]),fabs(a[1]-b[1])),fabs(a[2]-b[2]))
#define SET_VALUE(p,r,z)  ND(p,z,0) = ND(p,z,1) = ND(p,z,2) = r; 
#define COPY(p,x,z)    {  ND(p,z,0) = ND(p,x,0);                               \
                          ND(p,z,1) = ND(p,x,1);                               \
                          ND(p,z,2) = ND(p,x,2);  } 
#define INV(p,x,z)     {  ND(p,z,0) = -ND(p,x,0);                              \
                          ND(p,z,1) = -ND(p,x,1);                              \
                          ND(p,z,2) = -ND(p,x,2);  } 
#define ADD(p,x,y,z)   {  ND(p,z,0) = ND(p,x,0) + ND(p,y,0);                   \
                          ND(p,z,1) = ND(p,x,1) + ND(p,y,1);                   \
                          ND(p,z,2) = ND(p,x,2) + ND(p,y,2); } 
#define NSUBTR(p,x,y,z) {  ND(p,z,0) = ND(p,x,0) - ND(p,y,0);                  \
                           ND(p,z,1) = ND(p,x,1) - ND(p,y,1);                  \
                           ND(p,z,2) = ND(p,x,2) - ND(p,y,2); }
#define NASUBTR(p,w,x,y,z) { ND(p,z,0) = ND(p,w,0) + ND(p,x,0) - ND(p,y,0);    \
                             ND(p,z,1) = ND(p,w,1) + ND(p,x,1) - ND(p,y,1);    \
                             ND(p,z,2) = ND(p,w,2) + ND(p,x,2) - ND(p,y,2); }
#define MADD(p,r,x,y,z) {  ND(p,z,0) = (r)*ND(p,x,0) + ND(p,y,0);              \
                           ND(p,z,1) = (r)*ND(p,x,1) + ND(p,y,1);              \
                           ND(p,z,2) = (r)*ND(p,x,2) + ND(p,y,2);} 
#define MSUBTR(p,r,x,y,z) {  ND(p,z,0) = (r)*ND(p,x,0) - ND(p,y,0);            \
                             ND(p,z,1) = (r)*ND(p,x,1) - ND(p,y,1);            \
                             ND(p,z,2) = (r)*ND(p,x,2) - ND(p,y,2);} 
#define ADD_MULT(p,r,x,z) {ND(p,z,0) += (r)*ND(p,x,0);                         \
                           ND(p,z,1) += (r)*ND(p,x,1);                         \
                           ND(p,z,2) += (r)*ND(p,x,2);} 
#define SUBTR_MULT(p,r,x,z) {ND(p,z,0) -= (r)*ND(p,x,0);                       \
                             ND(p,z,1) -= (r)*ND(p,x,1);                       \
                             ND(p,z,2) -= (r)*ND(p,x,2);} 
#define DAMP(p,r,x,y,z) { ND(p,z,0) = ND(p,x,0) + (r)*(ND(p,y,0) - ND(p,x,0)); \
                          ND(p,z,1) = ND(p,x,1) + (r)*(ND(p,y,1) - ND(p,x,1)); \
                          ND(p,z,2) = ND(p,x,2) + (r)*(ND(p,y,2) - ND(p,x,2));} 
#define MULT(p,r,x,z)   {  ND(p,z,0) = ND(p,x,0)*(r);                          \
                           ND(p,z,1) = ND(p,x,1)*(r);                          \
                           ND(p,z,2) = ND(p,x,2)*(r);  }
#define MULT_S(p,r,x,z) {  ND(p,z,0) = ND(p,x,0)*r[0];                         \
                           ND(p,z,1) = ND(p,x,1)*r[1];                         \
                           ND(p,z,2) = ND(p,x,2)*r[2];  }
#define DIVIDE(p,x,y,z) {  ND(p,z,0) = ND(p,x,0)/ND(p,y,0);                    \
                           ND(p,z,1) = ND(p,x,1)/ND(p,y,1);                    \
                           ND(p,z,2) = ND(p,x,2)/ND(p,y,2);  }
#define MULTIPLY(p,x,y,z) {ND(p,z,0) = ND(p,x,0)*ND(p,y,0);                    \
                           ND(p,z,1) = ND(p,x,1)*ND(p,y,1);                    \
                           ND(p,z,2) = ND(p,x,2)*ND(p,y,2);  }
#define MSQRT(p,x,z)    {  ND(p,z,0) = sqrt(ND(p,x,0));                        \
                           ND(p,z,1) = sqrt(ND(p,x,1));                        \
                           ND(p,z,2) = sqrt(ND(p,x,2));  }
#define ADD_VMULT(n,y,a,b)   { ND(n[0],y,0) += a[0]*b[0];                      \
                               ND(n[1],y,0) += a[1]*b[0];                      \
                               ND(n[2],y,0) += a[2]*b[0];                      \
                               ND(n[3],y,0) += a[3]*b[0];                      \
                               ND(n[0],y,1) += a[0]*b[1];                      \
                               ND(n[1],y,1) += a[1]*b[1];                      \
                               ND(n[2],y,1) += a[2]*b[1];                      \
                               ND(n[3],y,1) += a[3]*b[1];                      \
                               ND(n[0],y,2) += a[0]*b[2];                      \
                               ND(n[1],y,2) += a[1]*b[2];                      \
                               ND(n[2],y,2) += a[2]*b[2];                      \
                               ND(n[3],y,2) += a[3]*b[2]; }
#define SUBTR_VMULT(n,y,a,b) { ND(n[0],y,0) -= a[0]*b[0];                      \
                               ND(n[1],y,0) -= a[1]*b[0];                      \
                               ND(n[2],y,0) -= a[2]*b[0];                      \
                               ND(n[3],y,0) -= a[3]*b[0];                      \
                               ND(n[0],y,1) -= a[0]*b[1];                      \
                               ND(n[1],y,1) -= a[1]*b[1];                      \
                               ND(n[2],y,1) -= a[2]*b[1];                      \
                               ND(n[3],y,1) -= a[3]*b[1];                      \
                               ND(n[0],y,2) -= a[0]*b[2];                      \
                               ND(n[1],y,2) -= a[1]*b[2];                      \
                               ND(n[2],y,2) -= a[2]*b[2];                      \
                               ND(n[3],y,2) -= a[3]*b[2]; }
#define ADD_SUM(n,x,a,b)     { a[0] += b[0]*ND(n[0],x,0) +                     \
                                       b[1]*ND(n[1],x,0) +                     \
                                       b[2]*ND(n[2],x,0) +                     \
                                       b[3]*ND(n[3],x,0);                      \
                               a[1] += b[0]*ND(n[0],x,1) +                     \
                                       b[1]*ND(n[1],x,1) +                     \
                                       b[2]*ND(n[2],x,1) +                     \
                                       b[3]*ND(n[3],x,1);                      \
                               a[2] += b[0]*ND(n[0],x,2) +                     \
                                       b[1]*ND(n[1],x,2) +                     \
                                       b[2]*ND(n[2],x,2) +                     \
                                       b[3]*ND(n[3],x,2); }
#define SUBTR_SUM(n,x,a,b)   { a[0] -= b[0]*ND(n[0],x,0) +                     \
                                       b[1]*ND(n[1],x,0) +                     \
                                       b[2]*ND(n[2],x,0) +                     \
                                       b[3]*ND(n[3],x,0);                      \
                               a[1] -= b[0]*ND(n[0],x,1) +                     \
                                       b[1]*ND(n[1],x,1) +                     \
                                       b[2]*ND(n[2],x,1) +                     \
                                       b[3]*ND(n[3],x,1);                      \
                               a[2] -= b[0]*ND(n[0],x,2) +                     \
                                       b[1]*ND(n[1],x,2) +                     \
                                       b[2]*ND(n[2],x,2) +                     \
                                       b[3]*ND(n[3],x,2); }
#define E_SUM(p,z)       (EDSN(p,z,0) + EDSN(p,z,1) + EDSN(p,z,2) + EDSN(p,z,3))
#define E_SET_VALUE(p,r,z) {EDSN(p,z,0) = EDSN(p,z,1) = r;                     \
                            EDSN(p,z,2) = EDSN(p,z,3) = r; }
#define E_ADD_VALUE(p,r,z) {EDSN(p,z,0) += r;   EDSN(p,z,1) += r;              \
                            EDSN(p,z,2) += r;   EDSN(p,z,3) += r; }
#define E_SUBTR_C(p,f,z) {  EDSN(p,z,0) -= f;                                  \
                            EDSN(p,z,1) -= f;                                  \
                            EDSN(p,z,2) -= f;                                  \
                            EDSN(p,z,3) -= f;  } 
#define E_COPY(p,x,z)    {  EDSN(p,z,0) = EDSN(p,x,0);                         \
                            EDSN(p,z,1) = EDSN(p,x,1);                         \
                            EDSN(p,z,2) = EDSN(p,x,2);                         \
                            EDSN(p,z,3) = EDSN(p,x,3);  } 
#define E_INV(p,x,z)     {  EDSN(p,z,0) = -EDSN(p,x,0);                        \
                            EDSN(p,z,1) = -EDSN(p,x,1);                        \
                            EDSN(p,z,2) = -EDSN(p,x,2);                        \
                            EDSN(p,z,3) = -EDSN(p,x,3);  } 
#define E_ADD(p,x,y,z)   {  EDSN(p,z,0) = EDSN(p,x,0) + EDSN(p,y,0);           \
                            EDSN(p,z,1) = EDSN(p,x,1) + EDSN(p,y,1);           \
                            EDSN(p,z,2) = EDSN(p,x,2) + EDSN(p,y,2);           \
                            EDSN(p,z,3) = EDSN(p,x,3) + EDSN(p,y,3); } 
#define E_NSUBTR(p,x,y,z) { EDSN(p,z,0) = EDSN(p,x,0) - EDSN(p,y,0);           \
                            EDSN(p,z,1) = EDSN(p,x,1) - EDSN(p,y,1);           \
                            EDSN(p,z,2) = EDSN(p,x,2) - EDSN(p,y,2);           \
                            EDSN(p,z,3) = EDSN(p,x,3) - EDSN(p,y,3); }
#define E_NASUBTR(p,w,x,y,z) { EDSN(p,z,0) =                                   \
                               EDSN(p,w,0) + EDSN(p,x,0) - EDSN(p,y,0);        \
                               EDSN(p,z,1) =                                   \
                               EDSN(p,w,1) + EDSN(p,x,1) - EDSN(p,y,1);        \
                               EDSN(p,z,2) =                                   \
                               EDSN(p,w,2) + EDSN(p,x,2) - EDSN(p,y,2);        \
                               EDSN(p,z,3) =                                   \
                               EDSN(p,w,3) + EDSN(p,x,3) - EDSN(p,y,3); }
#define E_MADD(p,r,x,y,z) {  EDSN(p,z,0) = (r)*EDSN(p,x,0) + EDSN(p,y,0);      \
                             EDSN(p,z,1) = (r)*EDSN(p,x,1) + EDSN(p,y,1);      \
                             EDSN(p,z,2) = (r)*EDSN(p,x,2) + EDSN(p,y,2);      \
                             EDSN(p,z,3) = (r)*EDSN(p,x,3) + EDSN(p,y,3);} 
#define E_ADD_MULT(p,r,x,z) {EDSN(p,z,0) += (r)*EDSN(p,x,0);                   \
                             EDSN(p,z,1) += (r)*EDSN(p,x,1);                   \
                             EDSN(p,z,2) += (r)*EDSN(p,x,2);                   \
                             EDSN(p,z,3) += (r)*EDSN(p,x,3);} 
#define E_SUBTR_MULT(p,r,x,z) {EDSN(p,z,0) -= (r)*EDSN(p,x,0);                 \
                               EDSN(p,z,1) -= (r)*EDSN(p,x,1);                 \
                               EDSN(p,z,2) -= (r)*EDSN(p,x,2);                 \
                               EDSN(p,z,3) -= (r)*EDSN(p,x,3);} 
#define E_DAMP(p,r,x,y,z) {  EDSN(p,z,0) =                                     \
                             EDSN(p,x,0) + (r)*(EDSN(p,y,0) - EDSN(p,x,0));    \
                             EDSN(p,z,1) =                                     \
                             EDSN(p,x,1) + (r)*(EDSN(p,y,1) - EDSN(p,x,1));    \
                             EDSN(p,z,2) =                                     \
                             EDSN(p,x,2) + (r)*(EDSN(p,y,2) - EDSN(p,x,2));    \
                             EDSN(p,z,3) =                                     \
                             EDSN(p,x,3) + (r)*(EDSN(p,y,3) - EDSN(p,x,3));} 
#define E_MULT(p,r,x,z)   {  EDSN(p,z,0) = EDSN(p,x,0)*(r);                    \
                             EDSN(p,z,1) = EDSN(p,x,1)*(r);                    \
                             EDSN(p,z,2) = EDSN(p,x,2)*(r);                    \
                             EDSN(p,z,3) = EDSN(p,x,3)*(r);  }
#define E_DIVIDE(p,x,y,z) {  EDSN(p,z,0) = EDSN(p,x,0)/EDSN(p,y,0);            \
                             EDSN(p,z,1) = EDSN(p,x,1)/EDSN(p,y,1);            \
                             EDSN(p,z,2) = EDSN(p,x,2)/EDSN(p,y,2);            \
                             EDSN(p,z,3) = EDSN(p,x,3)/EDSN(p,y,3);  }
#define E_MULTIPLY(p,x,y,z) {EDSN(p,z,0) = EDSN(p,x,0)*EDSN(p,y,0);            \
                             EDSN(p,z,1) = EDSN(p,x,1)*EDSN(p,y,1);            \
                             EDSN(p,z,2) = EDSN(p,x,2)*EDSN(p,y,2);            \
                             EDSN(p,z,3) = EDSN(p,x,3)*EDSN(p,y,3);  }
#define E_MSQRT(p,x,z)    {  EDSN(p,z,0) = sqrt(EDSN(p,x,0));                  \
                             EDSN(p,z,1) = sqrt(EDSN(p,x,1));                  \
                             EDSN(p,z,2) = sqrt(EDSN(p,x,2));                  \
                             EDSN(p,z,3) = sqrt(EDSN(p,x,3));  }
#define D_SET_VALUE(p,r,z) {FDDV(p,z,0,0) = FDDV(p,z,0,1) = FDDV(p,z,0,2) = r; \
                            FDDV(p,z,1,0) = FDDV(p,z,1,1) = FDDV(p,z,1,2) = r;} 
#define D_COPY(p,x,z)    {  FDDV(p,z,0,0) = FDDV(p,x,0,0);                     \
                            FDDV(p,z,0,1) = FDDV(p,x,0,1);                     \
                            FDDV(p,z,0,2) = FDDV(p,x,0,2);                     \
                            FDDV(p,z,1,0) = FDDV(p,x,1,0);                     \
                            FDDV(p,z,1,1) = FDDV(p,x,1,1);                     \
                            FDDV(p,z,1,2) = FDDV(p,x,1,2);  } 
#define D_INV(p,x,z)     {  FDDV(p,z,0,0) = -FDDV(p,x,0,0);                    \
                            FDDV(p,z,0,1) = -FDDV(p,x,0,1);                    \
                            FDDV(p,z,0,2) = -FDDV(p,x,0,2);                    \
                            FDDV(p,z,1,0) = -FDDV(p,x,1,0);                    \
                            FDDV(p,z,1,1) = -FDDV(p,x,1,1);                    \
                            FDDV(p,z,1,2) = -FDDV(p,x,1,2);  } 
#define D_ADD(p,x,y,z)   {  FDDV(p,z,0,0) = FDDV(p,x,0,0) + FDDV(p,y,0,0);     \
                            FDDV(p,z,0,1) = FDDV(p,x,0,1) + FDDV(p,y,0,1);     \
                            FDDV(p,z,0,2) = FDDV(p,x,0,2) + FDDV(p,y,0,2);     \
                            FDDV(p,z,1,0) = FDDV(p,x,1,0) + FDDV(p,y,1,0);     \
                            FDDV(p,z,1,1) = FDDV(p,x,1,1) + FDDV(p,y,1,1);     \
                            FDDV(p,z,1,2) = FDDV(p,x,1,2) + FDDV(p,y,1,2); } 
#define D_NSUBTR(p,x,y,z) {  FDDV(p,z,0,0) = FDDV(p,x,0,0) - FDDV(p,y,0,0);    \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1) - FDDV(p,y,0,1);    \
                             FDDV(p,z,0,2) = FDDV(p,x,0,2) - FDDV(p,y,0,2);    \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0) - FDDV(p,y,1,0);    \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1) - FDDV(p,y,1,1);    \
                             FDDV(p,z,1,2) = FDDV(p,x,1,2) - FDDV(p,y,1,2); }
#define D_NASUBTR(p,w,x,y,z) { FDDV(p,z,0,0) =                                 \
                               FDDV(p,w,0,0) + FDDV(p,x,0,0) - FDDV(p,y,0,0);  \
                               FDDV(p,z,0,1) =                                 \
                               FDDV(p,w,0,1) + FDDV(p,x,0,1) - FDDV(p,y,0,1);  \
                               FDDV(p,z,0,2) =                                 \
                               FDDV(p,w,0,2) + FDDV(p,x,0,2) - FDDV(p,y,0,2);  \
                               FDDV(p,z,1,0) =                                 \
                               FDDV(p,w,1,0) + FDDV(p,x,1,0) - FDDV(p,y,1,0);  \
                               FDDV(p,z,1,1) =                                 \
                               FDDV(p,w,1,1) + FDDV(p,x,1,1) - FDDV(p,y,1,1);  \
                               FDDV(p,z,1,2) =                                 \
                               FDDV(p,w,1,2) + FDDV(p,x,1,2) - FDDV(p,y,1,2); }
#define D_MADD(p,r,x,y,z) {  FDDV(p,z,0,0) = (r)*FDDV(p,x,0,0) + FDDV(p,y,0,0);\
                             FDDV(p,z,0,1) = (r)*FDDV(p,x,0,1) + FDDV(p,y,0,1);\
                             FDDV(p,z,0,2) = (r)*FDDV(p,x,0,2) + FDDV(p,y,0,2);\
                             FDDV(p,z,1,0) = (r)*FDDV(p,x,1,0) + FDDV(p,y,1,0);\
                             FDDV(p,z,1,1) = (r)*FDDV(p,x,1,1) + FDDV(p,y,1,1);\
                             FDDV(p,z,1,2) = (r)*FDDV(p,x,1,2) + FDDV(p,y,1,2);}
#define D_ADD_MULT(p,r,x,z) {FDDV(p,z,0,0) += (r)*FDDV(p,x,0,0);               \
                             FDDV(p,z,0,1) += (r)*FDDV(p,x,0,1);               \
                             FDDV(p,z,0,2) += (r)*FDDV(p,x,0,2);               \
                             FDDV(p,z,1,0) += (r)*FDDV(p,x,1,0);               \
                             FDDV(p,z,1,1) += (r)*FDDV(p,x,1,1);               \
                             FDDV(p,z,1,2) += (r)*FDDV(p,x,1,2);} 
#define D_SUBTR_MULT(p,r,x,z) {FDDV(p,z,0,0) -= (r)*FDDV(p,x,0,0);             \
                               FDDV(p,z,0,1) -= (r)*FDDV(p,x,0,1);             \
                               FDDV(p,z,0,2) -= (r)*FDDV(p,x,0,2);             \
                               FDDV(p,z,1,0) -= (r)*FDDV(p,x,1,0);             \
                               FDDV(p,z,1,1) -= (r)*FDDV(p,x,1,1);             \
                               FDDV(p,z,1,2) -= (r)*FDDV(p,x,1,2);} 
#define D_DAMP(p,r,x,y,z) {FDDV(p,z,0,0) =                                     \
                           FDDV(p,x,0,0) + (r)*(FDDV(p,y,0,0) - FDDV(p,x,0,0));\
                           FDDV(p,z,0,1) =                                     \
                           FDDV(p,x,0,1) + (r)*(FDDV(p,y,0,1) - FDDV(p,x,0,1));\
                           FDDV(p,z,0,2) =                                     \
                           FDDV(p,x,0,2) + (r)*(FDDV(p,y,0,2) - FDDV(p,x,0,2));\
                           FDDV(p,z,1,0) =                                     \
                           FDDV(p,x,1,0) + (r)*(FDDV(p,y,1,0) - FDDV(p,x,1,0));\
                           FDDV(p,z,1,1) =                                     \
                           FDDV(p,x,1,1) + (r)*(FDDV(p,y,1,1) - FDDV(p,x,1,1));\
                           FDDV(p,z,1,2) =                                     \
                           FDDV(p,x,1,2) + (r)*(FDDV(p,y,1,2) - FDDV(p,x,1,2));}
#define D_MULT(p,r,x,z)   {  FDDV(p,z,0,0) = FDDV(p,x,0,0)*(r);                \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)*(r);                \
                             FDDV(p,z,0,2) = FDDV(p,x,0,2)*(r);                \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)*(r);                \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)*(r);                \
                             FDDV(p,z,1,2) = FDDV(p,x,1,2)*(r);  }
#define D_DIVIDE(p,x,y,z) {  FDDV(p,z,0,0) = FDDV(p,x,0,0)/FDDV(p,y,0,0);      \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)/FDDV(p,y,0,1);      \
                             FDDV(p,z,0,2) = FDDV(p,x,0,2)/FDDV(p,y,0,2);      \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)/FDDV(p,y,1,0);      \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)/FDDV(p,y,1,1);      \
                             FDDV(p,z,1,2) = FDDV(p,x,1,2)/FDDV(p,y,1,2);  }
#define D_MULTIPLY(p,x,y,z) {FDDV(p,z,0,0) = FDDV(p,x,0,0)*FDDV(p,y,0,0);      \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)*FDDV(p,y,0,1);      \
                             FDDV(p,z,0,2) = FDDV(p,x,0,2)*FDDV(p,y,0,2);      \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)*FDDV(p,y,1,0);      \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)*FDDV(p,y,1,1);      \
                             FDDV(p,z,1,2) = FDDV(p,x,1,2)*FDDV(p,y,1,2);  }
#define D_MSQRT(p,x,z)    {  FDDV(p,z,0,0) = sqrt(FDDV(p,x,0,0));              \
                             FDDV(p,z,0,1) = sqrt(FDDV(p,x,0,1));              \
                             FDDV(p,z,0,2) = sqrt(FDDV(p,x,0,2));              \
                             FDDV(p,z,1,0) = sqrt(FDDV(p,x,1,0));              \
                             FDDV(p,z,1,1) = sqrt(FDDV(p,x,1,1));              \
                             FDDV(p,z,1,2) = sqrt(FDDV(p,x,1,2));  }
#define SDIVIDE(p,x,y,z) { ND(p,z,0) = ND(p,x,0)/NDS(p,y);                     \
                           ND(p,z,1) = ND(p,x,1)/NDS(p,y);                     \
                           ND(p,z,2) = ND(p,x,2)/NDS(p,y); }
#define SMULTIPLY(p,x,y,z) { ND(p,z,0) = ND(p,x,0)*NDS(p,y);                   \
                             ND(p,z,1) = ND(p,x,1)*NDS(p,y);                   \
                             ND(p,z,2) = ND(p,x,2)*NDS(p,y); }
#define SET_COEFFNN(p,k,a,b,c,d,e,f,g,h,i)  COEFFNN(p,k,0,0) += a;             \
                                            COEFFNN(p,k,0,1) += b;             \
                                            COEFFNN(p,k,0,2) += c;             \
                                            COEFFNN(p,k,1,0) += d;             \
                                            COEFFNN(p,k,1,1) += e;             \
                                            COEFFNN(p,k,1,2) += f;             \
                                            COEFFNN(p,k,2,0) += g;             \
                                            COEFFNN(p,k,2,1) += h;             \
                                            COEFFNN(p,k,2,2) += i
#define LG_SET_VALUE(p,r,z) NDLG(p,z,0) = NDLG(p,z,1) = r;
#define LG_COPY(p,x,z)  {  NDLG(p,z,0) = NDLG(p,x,0);                          \
                           NDLG(p,z,1) = NDLG(p,x,1);  }
#define LG_INV(p,x,z)   {  NDLG(p,z,0) = -NDLG(p,x,0);                         \
                           NDLG(p,z,1) = -NDLG(p,x,1);  }
#define LG_ADD(p,x,y,z) {  NDLG(p,z,0) = NDLG(p,x,0) + NDLG(p,y,0);            \
                           NDLG(p,z,1) = NDLG(p,x,1) + NDLG(p,y,1);  }
#define LG_SUBTR(p,x,y,z)  {  NDLG(p,z,0) = NDLG(p,x,0) - NDLG(p,y,0);         \
                              NDLG(p,z,1) = NDLG(p,x,1) - NDLG(p,y,1);  }
#define LG_MADD(p,r,x,y,z) {  NDLG(p,z,0) = (r)*NDLG(p,x,0) + NDLG(p,y,0);     \
                              NDLG(p,z,1) = (r)*NDLG(p,x,1) + NDLG(p,y,1);  }
#define LG_DAMP(p,r,x,y,z) {  NDLG(p,z,0) = NDLG(p,x,0) +                      \
                                             (r)*(NDLG(p,y,0) - NDLG(p,x,0));  \
                              NDLG(p,z,1) = NDLG(p,x,1) +                      \
                                             (r)*(NDLG(p,y,1) - NDLG(p,x,1)); }
#define LG_MULT(p,r,x,z)   {  NDLG(p,z,0) = NDLG(p,x,0)*(r);                   \
                              NDLG(p,z,1) = NDLG(p,x,1)*(r);  }
#define LG_DOT(p,x,y)         NDLG(p,x,0)*NDLG(p,y,0) + NDLG(p,x,1)*NDLG(p,y,1);
#define LG_DIVIDE(p,x,y,z) {  NDLG(p,z,0) = NDLG(p,x,0)/NDLG(p,y,0);           \
                              NDLG(p,z,1) = NDLG(p,x,1)/NDLG(p,y,1);  }
#define LG_MULTIPLY(p,x,y,z) {  NDLG(p,z,0) = NDLG(p,x,0)*NDLG(p,y,0);         \
                                NDLG(p,z,1) = NDLG(p,x,1)*NDLG(p,y,1);  }
#define LG_MSQRT(p,x,z)      {  NDLG(p,z,0) = sqrt(NDLG(p,x,0));               \
                                NDLG(p,z,1) = sqrt(NDLG(p,x,1));  }

#else
#define DIM1            1        /* DIM - 1                           */
#define DIM2            3        /* DIM + 1                           */
#define DIM3            4        /* 2*DIM                             */
#define DIM4            5        /* 2*DIM + 1                         */
#if ELEMENT_TYPE == SIMPLEX
#define NVERT		3	 /* max number of vertices of an element */
#define SIDES		3	 /* max number of sides of an element */
#define EDGES	        3        /* max number of edges of an element */
#define IS_FTYPE(e,mask) (FTYPE(e->f[0]) & mask || FTYPE(e->f[1]) & mask ||    \
                          FTYPE(e->f[2]) & mask)
#define VOLUME(e)       volume(e->n[0]->myvertex->x,e->n[1]->myvertex->x,      \
                                                    e->n[2]->myvertex->x)
#define BAR_COORD(p,b)   barycentric_coordinates((p)->n[0]->myvertex->x,       \
                          (p)->n[1]->myvertex->x,(p)->n[2]->myvertex->x,b)
#define IS_B_EL(e)      (IS_DN(e->n[0]) || IS_DN(e->n[1]) || IS_DN(e->n[2]))
#define FACES_AND_NODES_OF_ELEMENT(e,ef,n0,n1)                                 \
                         { ef[0] = e->f[0]; n0[0] = e->n[1]; n1[0] = e->n[2];  \
                           ef[1] = e->f[1]; n0[1] = e->n[2]; n1[1] = e->n[0];  \
                           ef[2] = e->f[2]; n0[2] = e->n[0]; n1[2] = e->n[1]; }
#elif ELEMENT_TYPE == CUBE
#define NVERT		4	 /* max number of vertices of an element */
#define SIDES		4	 /* max number of sides of an element */
#define EDGES	        4        /* max number of edges of an element */
#define IS_FTYPE(e,mask) (FTYPE(e->f[0]) & mask || FTYPE(e->f[1]) & mask ||    \
                          FTYPE(e->f[2]) & mask || FTYPE(e->f[3]) & mask)
#define VOLUME(e)       volume_q(e->n[0]->myvertex->x,e->n[1]->myvertex->x,    \
                                 e->n[2]->myvertex->x,e->n[3]->myvertex->x)
#define IS_B_EL(e)      (IS_DN(e->n[0]) || IS_DN(e->n[1]) ||                   \
                         IS_DN(e->n[2]) || IS_DN(e->n[3]))
#define FACES_AND_NODES_OF_ELEMENT(e,ef,n0,n1)                                 \
                         { ef[0] = e->f[0]; n0[0] = e->n[0]; n1[0] = e->n[1];  \
                           ef[1] = e->f[1]; n0[1] = e->n[1]; n1[1] = e->n[2];  \
                           ef[2] = e->f[2]; n0[2] = e->n[2]; n1[2] = e->n[3];  \
                           ef[3] = e->f[3]; n0[3] = e->n[3]; n1[3] = e->n[0]; }
#endif
#define SONS		4        /* max number of sons of an element  */
#define FSONS		2        /* max number of sons of a face      */
#define AVER_N(e,p,i)   (ND(e->n[0],p,i) + ND(e->n[1],p,i) + ND(e->n[2],p,i))/3.
#define AVER_NS(e,p)    (NDS(e->n[0],p) + NDS(e->n[1],p) + NDS(e->n[2],p))/3.
#define AVER_FS(e,p)    (FD(e->f[0],p) + FD(e->f[1],p) + FD(e->f[2],p))/3.
#define AVER_ESN(e,p)   (EDSN(e,p,0) + EDSN(e,p,1) + EDSN(e,p,2))/3.
#define VERTICES_OF_ELEMENT(x1,x2,x3,e) x1 = e->n[0]->myvertex->x;             \
                                        x2 = e->n[1]->myvertex->x;             \
                                        x3 = e->n[2]->myvertex->x
#define NODES_OF_ELEMENT(n1,n2,n3,e)    n1=e->n[0]; n2=e->n[1]; n3=e->n[2]
#define INODES_OF_ELEMENT(n1,n2,n3,i1,i2,i3,e)                                 \
                                        n1=e->n[i1]; n2=e->n[i2]; n3=e->n[i3]; 
#define FACES_OF_ELEMENT(f1,f2,f3,e)    f1=e->f[0]; f2=e->f[1]; f3=e->f[2]
#define TOPNODES_OF_ELEMENT(n1,n2,n3,e) n1=TOP_NODE(e->n[0]);                  \
                  n2=TOP_NODE(e->n[1]); n3=TOP_NODE(e->n[2])
#define TOPFACES_OF_ELEMENT(f1,f2,f3,e) f1=top_face(e->f[0]);                  \
                   f2=top_face(e->f[1]); f3=top_face(e->f[2])
#define LTOPNODES_OF_ELEMENT(n1,n2,n3,e,g) n1=ltop_node(e->n[0],g);            \
          n2=ltop_node(e->n[1],g); n3=ltop_node(e->n[2],g)
#define LTOPFACES_OF_ELEMENT(f1,f2,f3,e,g) f1=ltop_face(e->f[0],g);            \
          f2=ltop_face(e->f[1],g); f3=ltop_face(e->f[2],g)
#define S_NODE_VALUES(n0,n1,n2,f,g0,g1,g2)       { g0 = f(n0->myvertex->x);    \
                                                   g1 = f(n1->myvertex->x);    \
                                                   g2 = f(n2->myvertex->x); }
#define V_NODE_VALUES(n0,n1,n2,bb0,bb1,bb_0,bb_1,bb_2)                         \
         { bb_0[0] = bb0(n0->myvertex->x);   bb_0[1] = bb1(n0->myvertex->x);   \
           bb_1[0] = bb0(n1->myvertex->x);   bb_1[1] = bb1(n1->myvertex->x);   \
           bb_2[0] = bb0(n2->myvertex->x);   bb_2[1] = bb1(n2->myvertex->x); }
#define V_MID_VALUES(x01,x02,x12,bb0,bb1,bb_0,bb_1,bb_2)                       \
                                 { bb_0[0] = bb0(x12);   bb_0[1] = bb1(x12);   \
                                   bb_1[0] = bb0(x02);   bb_1[1] = bb1(x02);   \
                                   bb_2[0] = bb0(x01);   bb_2[1] = bb1(x01); }
#define NORMAL_VECTORS(nn1,nn2,nn3,e)                                          \
                    normal_vector(e->n[1]->myvertex->x,e->n[2]->myvertex->x,   \
                           e->f[0],nn1);                                       \
                    normal_vector(e->n[0]->myvertex->x,e->n[2]->myvertex->x,   \
                           e->f[1],nn2);                                       \
                    normal_vector(e->n[0]->myvertex->x,e->n[1]->myvertex->x,   \
                           e->f[2],nn3)
#define ARE_NEIGHBOURS(p,q)  (p->n[0]==q->n[0] || p->n[0]==q->n[1] ||          \
                              p->n[0]==q->n[2] ||                              \
                              p->n[1]==q->n[0] || p->n[1]==q->n[1] ||          \
                              p->n[1]==q->n[2] ||                              \
                              p->n[2]==q->n[0] || p->n[2]==q->n[1] ||          \
                              p->n[2]==q->n[2])
#define ARE_FNEIGHBOURS(p,q) (p->f[0]==q->f[0] || p->f[0]==q->f[1] ||          \
                              p->f[0]==q->f[2] ||                              \
                              p->f[1]==q->f[0] || p->f[1]==q->f[1] ||          \
                              p->f[1]==q->f[2] ||                              \
                              p->f[2]==q->f[0] || p->f[2]==q->f[1] ||          \
                              p->f[2]==q->f[2])
#define FACE_CONTAINED(fa,p) (fa==p->f[0] || fa==p->f[1] || fa==p->f[2])
#define GIVE_INDEX(p,n,i)  { if (p[0] == n) i = 0;                             \
                             else if (p[1] == n) i = 1;                        \
                             else if (p[2] == n) i = 2;                        \
                             else i = -1; }
#if ELEMENT_TYPE == SIMPLEX
#define IS_IN_SQUARE(e,xm,ym)  ( e->n[0]->myvertex->x[0] <= xm &&              \
                                 e->n[1]->myvertex->x[0] <= xm &&              \
                                 e->n[2]->myvertex->x[0] <= xm &&              \
                                 e->n[0]->myvertex->x[1] <= ym &&              \
                                 e->n[1]->myvertex->x[1] <= ym &&              \
                                 e->n[2]->myvertex->x[1] <= ym )
#define VERTEX_ON_LINE(e,k,d)  ( fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 ||   \
                                 fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 ||   \
                                 fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 )
#define EDGE_ON_LINE(e,k,d) ( (fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 )  ||  \
                              (fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 )  ||  \
                              (fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 ) )
#define IS_IN_G(pel,out_mask)               (NTYPE(pel->n[0]) & out_mask ||    \
              NTYPE(pel->n[1]) & out_mask || NTYPE(pel->n[2]) & out_mask)
#define IS_FROM_G2(pel,out_mask)                                               \
        (((NTYPE(pel->n[0]) & out_mask) && (NTYPE(pel->n[1]) & out_mask)) ||   \
         ((NTYPE(pel->n[1]) & out_mask) && (NTYPE(pel->n[2]) & out_mask)) ||   \
         ((NTYPE(pel->n[2]) & out_mask) && (NTYPE(pel->n[0]) & out_mask)))
#define HAS_DN(pel) (IS_DN(pel->n[0]) && IS_DN(pel->n[1]) && IS_DN(pel->n[2]))
#elif ELEMENT_TYPE == CUBE
#define IS_IN_SQUARE(e,xm,ym)  ( e->n[0]->myvertex->x[0] <= xm &&              \
                                 e->n[1]->myvertex->x[0] <= xm &&              \
                                 e->n[2]->myvertex->x[0] <= xm &&              \
                                 e->n[3]->myvertex->x[0] <= xm &&              \
                                 e->n[0]->myvertex->x[1] <= ym &&              \
                                 e->n[1]->myvertex->x[1] <= ym &&              \
                                 e->n[2]->myvertex->x[1] <= ym &&              \
                                 e->n[3]->myvertex->x[1] <= ym )
#define VERTEX_ON_LINE(e,k,d)  ( fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 ||   \
                                 fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 ||   \
                                 fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 ||   \
                                 fabs(e->n[3]->myvertex->x[k]-d) < 1.e-15 )
#define EDGE_ON_LINE(e,k,d) ( (fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 )  ||  \
                              (fabs(e->n[1]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 )  ||  \
                              (fabs(e->n[2]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[3]->myvertex->x[k]-d) < 1.e-15 )  ||  \
                              (fabs(e->n[3]->myvertex->x[k]-d) < 1.e-15 &&     \
                               fabs(e->n[0]->myvertex->x[k]-d) < 1.e-15 ) )
#else
#define IS_IN_SQUARE(e,xm,ym)  0
#endif

#if F_DATA & CURVED_FACE_MIDDLE

#define C_DATA_OF_ELEMENT(n0,n1,n2,f0,f1,f2,x0,x1,x2,e)                        \
{ if (e->f[0]->c_midpoint){ NODES_OF_ELEMENT(n0,n1,n2,e);                      \
                            FACES_OF_ELEMENT(f0,f1,f2,e);                      \
                            VERTICES_OF_ELEMENT(x0,x1,x2,e); }                 \
  else if (e->f[1]->c_midpoint){ NODES_OF_ELEMENT(n2,n0,n1,e);                 \
                                 FACES_OF_ELEMENT(f2,f0,f1,e);                 \
                                 VERTICES_OF_ELEMENT(x2,x0,x1,e); }            \
  else{ NODES_OF_ELEMENT(n1,n2,n0,e);                                          \
        FACES_OF_ELEMENT(f1,f2,f0,e);                                          \
        VERTICES_OF_ELEMENT(x1,x2,x0,e); }} 

#define LTOP_C_DATA_OF_ELEMENT(n0,n1,n2,f0,f1,f2,x0,x1,x2,e,g)                 \
{ if (e->f[0]->c_midpoint){ LTOPNODES_OF_ELEMENT(n0,n1,n2,e,g);                \
                            LTOPFACES_OF_ELEMENT(f0,f1,f2,e,g);                \
                            VERTICES_OF_ELEMENT(x0,x1,x2,e); }                 \
  else if (e->f[1]->c_midpoint){ LTOPNODES_OF_ELEMENT(n2,n0,n1,e,g);           \
                                 LTOPFACES_OF_ELEMENT(f2,f0,f1,e,g);           \
                                 VERTICES_OF_ELEMENT(x2,x0,x1,e); }            \
  else{ LTOPNODES_OF_ELEMENT(n1,n2,n0,e,g);                                    \
        LTOPFACES_OF_ELEMENT(f1,f2,f0,e,g);                                    \
        VERTICES_OF_ELEMENT(x1,x2,x0,e); }} 

#endif

#define VERTICES_OF_4ELEMENT(x1,x2,x3,x4,e) x1 = e->n[0]->myvertex->x;         \
                                            x2 = e->n[1]->myvertex->x;         \
                                            x3 = e->n[2]->myvertex->x;         \
                                            x4 = e->n[3]->myvertex->x
#define NODES_OF_4ELEMENT(n1,n2,n3,n4,e)    n1=e->n[0]; n2=e->n[1];            \
                                            n3=e->n[2]; n4=e->n[3]
#define FACES_OF_4ELEMENT(f1,f2,f3,f4,e)    f1=e->f[0]; f2=e->f[1];            \
                                            f3=e->f[2]; f4=e->f[3]
#define TOPNODES_OF_4ELEMENT(n1,n2,n3,n4,e) n1=TOP_NODE(e->n[0]);              \
              n2=TOP_NODE(e->n[1]); n3=TOP_NODE(e->n[2]); n4=TOP_NODE(e->n[3])
#define TOPFACES_OF_4ELEMENT(f1,f2,f3,f4,e) f1=top_face(e->f[0]);              \
              f2=top_face(e->f[1]); f3=top_face(e->f[2]); f4=top_face(e->f[3])
#define LTOPNODES_OF_4ELEMENT(n1,n2,n3,n4,e,g) n1=ltop_node(e->n[0],g);        \
     n2=ltop_node(e->n[1],g); n3=ltop_node(e->n[2],g); n4=ltop_node(e->n[3],g)
#define LTOPFACES_OF_4ELEMENT(f1,f2,f3,f4,e,g) f1=ltop_face(e->f[0],g);        \
     f2=ltop_face(e->f[1],g); f3=ltop_face(e->f[2],g); f4=ltop_face(e->f[3],g)
#define S_NODE_4VALUES(n0,n1,n2,n3,f,g0,g1,g2,g3)   { g0 = f(n0->myvertex->x); \
                                                      g1 = f(n1->myvertex->x); \
                                                      g2 = f(n2->myvertex->x); \
                                                      g3 = f(n3->myvertex->x); }

#define IS_SON_FACE(f,g)       ( f == g->sons[0] || f == g->sons[1] )
#if PROBLEM & WITH_F
#define RHS_FOR_FACES   FD(fa0,f) +=  integr3(n0,n1,n2,nn0,rdetB,f01,f02)*FMULT;
#define RHS_FOR_FACESN  FD(fa0,f) += integr3c(n0,n1,n2,nn0,rdetB,f01,f02)*FMULT;
#endif
#define IS_IN_SIMPLEX(b,x,eps)    (DOT(b[0],x)+b[0][2] >= eps &&               \
                                   DOT(b[1],x)+b[1][2] >= eps &&               \
                                   DOT(b[2],x)+b[2][2] >= eps)
#define IS_IN(n,n1,n2,n3) ( n == n1 || n == n2 || n == n3 )
#define NOT_ALL_FN(p)    (NOT_FN(p->n[0]) || NOT_FN(p->n[1]) || NOT_FN(p->n[2]))
#define SUBTR(a,b,c)      c[0]=a[0]-b[0]; c[1]=a[1]-b[1]
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1])
#define DOTS(a,b)        (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define DOT2(a,b,c)      (a[0]*(b[0]+4.*c[0]) + a[1]*(b[1]+4.*c[1]))
#define DOT_DIFF(a,b,c,d)  ( (a[0]-b[0])*(c[0]-d[0]) + (a[1]-b[1])*(c[1]-d[1]) )
#define DISTANCE(a,b)    (sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])))
#define ORT_VECT(a,b)    {a[0]=-b[1]; a[1]=b[0];}
#define THREE_ORT_VECT(a1,a2,a3,b1,b2,b3)  {ORT_VECT(a1,b1);                   \
                                            ORT_VECT(a2,b2);  ORT_VECT(a3,b3);}
#define SET1(a,b)        {a[0]=b[0]; a[1]=b[1];}
#define SSET1(a,b)       {a[0]=b[0]; a[1]=b[1]; a[2]=b[2];}
#define SET2(a,b,p)      {a[0]=b[0]*(p); a[1]=b[1]*(p);}
#define SET3(a,b,p,q,r)  {r=(p)*(q); a[0]=b[0]*r; a[1]=b[1]*r;}
#define SET4(a,b,p)      {a[0]+=b[0]*(p); a[1]+=b[1]*(p);}
#define SET5(a,b)        {a[0]+=b[0];   a[1]+=b[1];}
#define SET6(a,b)        {a[0]-=b[0];   a[1]-=b[1];}
#define SET7(a,p)        {a[0]=p; a[1]=p;}
#define SSET7(a,p)       {a[0]=p; a[1]=p; a[2]=p;}
#define SET8(a,b)        {a[0]=-b[0]; a[1]=-b[1];}
#define SET9(a,b,p)      {a[0]-=b[0]*(p); a[1]-=b[1]*(p);}
#define SET10(a,b,c)     {a[0]=b[0]+c[0]; a[1]=b[1]+c[1];}
#define SET11(a,b,c)     {a[0]=b[0]-c[0]; a[1]=b[1]-c[1];}
#define SET12(a,b,c,p)   {a[0]=(b[0]-c[0])/(p); a[1]=(b[1]-c[1])/(p);}
#define SET13(a,p)       {a[0]/=(p); a[1]/=(p);}
#define SET14(z,a,b,c)   {z[0]=a[0]+b[0]+c[0]; z[1]=a[1]+b[1]+c[1];}
#define SET15(z,a,b,c,d) {z[0]=a[0]+b[0]+c[0]+d[0]; z[1]=a[1]+b[1]+c[1]+d[1];}
#define SET16(z,a,b,c,d,e) {z[0]=a[0]+b[0]+c[0]+d[0]+e[0];                     \
                            z[1]=a[1]+b[1]+c[1]+d[1]+e[1];}
#define SET17(a,b,c)     {a[0]=b[0]+b[0]+c[0]; a[1]=b[1]+b[1]+c[1];}
#define SET18(a,b,c)     {a[0]=b[0]+b[0]-c[0]; a[1]=b[1]+b[1]-c[1];}
#define SET19(a,b,c)     {a[0]=b[0]+b[0]+b[0]-c[0]-c[0];                       \
                          a[1]=b[1]+b[1]+b[1]-c[1]-c[1];}
#define SET20(z,a,p,b,q) {z[0]=a[0]*(p) + b[0]*(q);  z[1]=a[1]*(p) + b[1]*(q);}
#define SET21(z,a,p,b,q) {z[0]=a[0]/(p) + b[0]/(q);  z[1]=a[1]/(p) + b[1]/(q);}
#define SET22(a,b,p)     {a[0]=b[0]/(p); a[1]=b[1]/(p);}
#define SET23(a,b,c,p)   {a[0]=b[0]+c[0]*(p); a[1]=b[1]+c[1]*(p);}
#define SET24(a,b,c,d)   {a[0]=2.*(b[0]+b[0]-c[0]-d[0]);                       \
                          a[1]=2.*(b[1]+b[1]-c[1]-d[1]);}
#define SET25(a,b,c)     {a[0]=fabs(b[0]-c[0]); a[1]=fabs(b[1]-c[1]);}
#define SET26(a,p)       {a[0]+=p; a[1]+=p; a[2]+=p;}
#define SET27(z,a,p,b,q,c,r) {z[0]=a[0]*(p) + b[0]*(q) + c[0]*(r);             \
                              z[1]=a[1]*(p) + b[1]*(q) + c[1]*(r);}
#define SET28(a,b,c,p)   {a[0]=(b[0]-c[0])*(p); a[1]=(b[1]-c[1])*(p);}
#define SET29(a,b,c,d,p) {a[0]=b[0]+c[0]+d[0]*(p); a[1]=b[1]+c[1]+d[1]*(p);}
#define MSET1(a,b)       {a[0][0]=b[0][0]; a[0][1]=b[0][1];                    \
                          a[1][0]=b[1][0]; a[1][1]=b[1][1];}
#define MSET2(a,b,p)     {a[0] = b[0][0]*p[0] + b[0][1]*p[1];                  \
                          a[1] = b[1][0]*p[0] + b[1][1]*p[1];}
#define MSET3(a,b,p)     {a[0] = b[0][0]*p[0]; a[1] = b[1][1]*p[1];}
#define MSET33(a,b,p)    {a[0] = p[0]/b[0][0]; a[1] = p[1]/b[1][1];}
#define MSET4(a,b,p)     {a[0] += b[0][0]*p[0] + b[0][1]*p[1];                 \
                          a[1] += b[1][0]*p[0] + b[1][1]*p[1];}
#define MSET44(a,b,p)    {a[0][0] += b[0][0]*(p); a[0][1] += b[0][1]*(p);      \
                          a[1][0] += b[1][0]*(p); a[1][1] += b[1][1]*(p);}
#define MSET5(a,b)       {a[0][0] += b[0][0]; a[0][1] += b[0][1];              \
                          a[1][0] += b[1][0]; a[1][1] += b[1][1];}
#define MSET7(a,p)       {a[0][0]=p; a[0][1]=p;                                \
                          a[1][0]=p; a[1][1]=p;}
#define MSET9(a,b,p)     {a[0] -= b[0][0]*p[0] + b[0][1]*p[1];                 \
                          a[1] -= b[1][0]*p[0] + b[1][1]*p[1];}
#define V_LIN_VALUE(y,a,c,x0,x1)  {y[0] = a[0][0]*(x0) + a[0][1]*(x1) + c[0];  \
                                   y[1] = a[1][0]*(x0) + a[1][1]*(x1) + c[1];}
#define V_BILIN_VALUE(y,a,c,alpha,x0,x1)                                       \
             {y[0] = a[0][0]*(x0) + a[0][1]*(x1) + c[0] + alpha[0]*(x0)*(x1);  \
              y[1] = a[1][0]*(x0) + a[1][1]*(x1) + c[1] + alpha[1]*(x0)*(x1);}
#define REF_LIN(x,u0,u1,u2)      (u0*(1-x[0]-x[1])+u1*x[0]+u2*x[1])
#define REF_BILIN(x,u0,u1,u2,u3) (u0*r_q1_0(x)+u1*r_q1_1(x)+                   \
                                  u2*r_q1_2(x)+u3*r_q1_3(x))
#define DX_REF_BILIN(x,u0,u1,u2,u3) (u0*r_q1_0_0(x)+u1*r_q1_1_0(x)+            \
                                     u2*r_q1_2_0(x)+u3*r_q1_3_0(x))
#define DY_REF_BILIN(x,u0,u1,u2,u3) (u0*r_q1_0_1(x)+u1*r_q1_1_1(x)+            \
                                     u2*r_q1_2_1(x)+u3*r_q1_3_1(x))
#define REF_Q1ROT(x,u0,u1,u2,u3) (u0*r_phi0(x)+u1*r_phi1(x)+                   \
                                  u2*r_phi2(x)+u3*r_phi3(x))
#define DX_REF_Q1ROT(x,u0,u1,u2,u3) (u0*r_phi0_0(x)+u1*r_phi1_0(x)+            \
                                     u2*r_phi2_0(x)+u3*r_phi3_0(x))
#define DY_REF_Q1ROT(x,u0,u1,u2,u3) (u0*r_phi0_1(x)+u1*r_phi1_1(x)+            \
                                     u2*r_phi2_1(x)+u3*r_phi3_1(x))
#define REF_QUADR(x,u0,u1,u2,u01,u02,u12)                                      \
             (u0*(1-x[0]-x[1])+u1*x[0]+u2*x[1]+                                \
              u01*x[0]*(1-x[0]-x[1])+u02*x[1]*(1-x[0]-x[1])+u12*x[0]*x[1])
#define DX_REF_QUADR(x,u0,u1,u2,u01,u02,u12)                                   \
                             (-u0 + u1 + u01*(1-2.*x[0]-x[1]) + (u12-u02)*x[1])
#define DY_REF_QUADR(x,u0,u1,u2,u01,u02,u12)                                   \
                             (-u0 + u2 + u02*(1-x[0]-2.*x[1]) + (u12-u01)*x[0])
#define AVERAGE(x,y,z)    {z[0] = (x[0]+y[0])/2.; z[1] = (x[1]+y[1])/2.;}
#define ADDMULT(a,x,b,y,z) z[0] = (a)*x[0]+(b)*y[0]; z[1] = (a)*x[1]+(b)*y[1]
#define POINT2(x,y,z,a,c) z[0] = ((a)*x[0]+y[0])/(c); z[1] = ((a)*x[1]+y[1])/(c)
/*  z=(w+x+y)/3.0  */
#define POINT3(w,x,y,z)   z[0]=(w[0]+x[0]+y[0])/3.0; z[1]=(w[1]+x[1]+y[1])/3.0
#define POINT4(v,w,x,y,z) z[0]=(v[0]+w[0]+x[0]+y[0])/4.0;                      \
                          z[1]=(v[1]+w[1]+x[1]+y[1])/4.0
#define POINT5(w,x,y,z,a,c)     z[0] = ((a)*w[0]+x[0]+y[0])/(c);               \
                                z[1] = ((a)*w[1]+x[1]+y[1])/(c)  
#define MIDPOINTS(x0,x1,x2,x01,x02,x12) {  AVERAGE(x0,x1,x01);                 \
                                           AVERAGE(x0,x2,x02);                 \
                                           AVERAGE(x1,x2,x12);  }
#define V_POINT_VALUES_6(x1,x2,x3,x4,x5,x6,v1,v2,v3,v4,v5,v6,f0,f1)            \
                                {  v1[0]=f0(x1); v2[0]=f0(x2); v3[0]=f0(x3);   \
                                   v4[0]=f0(x4); v5[0]=f0(x5); v6[0]=f0(x6);   \
                                   v1[1]=f1(x1); v2[1]=f1(x2); v3[1]=f1(x3);   \
                                   v4[1]=f1(x4); v5[1]=f1(x5); v6[1]=f1(x6);  }
#define MAX_DIFF(a,b)     MAX(fabs(a[0]-b[0]),fabs(a[1]-b[1]))
#define INTEGR1_2(f0,f1,f2,f01,f02,f12)  ( (2.*(f0) - (f1) - (f2))/60. +       \
                                           (2.*((f01)+(f02)) + (f12))/15. )
#define INTEGR3_2(f0,f01,f02,f12)   ((4.*((f01)+(f02)) + 8.*(f12) - (f0))/180.)
#define INTEGR4_2(f0,f001,f002,z)   ((z + (f0) + 9.*((f001)+(f002)))/120.)
#define INTEGR5_2(f0,f112,f221,f001,f002,z) ((z + (f0) + 6.*((f112)+(f221))    \
                                                 -9.*((f001)+(f002)))/840.)
#define SET_VALUE(p,r,z)  ND(p,z,0) = ND(p,z,1) = r; 
#define COPY(p,x,z)    {  ND(p,z,0) = ND(p,x,0);                               \
                          ND(p,z,1) = ND(p,x,1);  }
#define INV(p,x,z)     {  ND(p,z,0) = -ND(p,x,0);                              \
                          ND(p,z,1) = -ND(p,x,1); }
#define ADD(p,x,y,z)   {  ND(p,z,0) = ND(p,x,0) + ND(p,y,0);                   \
                          ND(p,z,1) = ND(p,x,1) + ND(p,y,1); }
#define NSUBTR(p,x,y,z) {  ND(p,z,0) = ND(p,x,0) - ND(p,y,0);                  \
                           ND(p,z,1) = ND(p,x,1) - ND(p,y,1); }
#define NASUBTR(p,w,x,y,z) { ND(p,z,0) = ND(p,w,0) + ND(p,x,0) - ND(p,y,0);    \
                             ND(p,z,1) = ND(p,w,1) + ND(p,x,1) - ND(p,y,1); }
#define MADD(p,r,x,y,z) {  ND(p,z,0) = (r)*ND(p,x,0) + ND(p,y,0);              \
                           ND(p,z,1) = (r)*ND(p,x,1) + ND(p,y,1); }
#define MSUBTR(p,r,x,y,z) {  ND(p,z,0) = (r)*ND(p,x,0) - ND(p,y,0);            \
                             ND(p,z,1) = (r)*ND(p,x,1) - ND(p,y,1); }
#define ADD_MULT(p,r,x,z) {ND(p,z,0) += (r)*ND(p,x,0);                         \
                           ND(p,z,1) += (r)*ND(p,x,1); }
#define SUBTR_MULT(p,r,x,z) {ND(p,z,0) -= (r)*ND(p,x,0);                       \
                             ND(p,z,1) -= (r)*ND(p,x,1);}
#define DAMP(p,r,x,y,z) {  ND(p,z,0) = ND(p,x,0) + (r)*(ND(p,y,0) - ND(p,x,0));\
                           ND(p,z,1) = ND(p,x,1) + (r)*(ND(p,y,1) - ND(p,x,1)); }
#define MULT(p,r,x,z)   {  ND(p,z,0) = ND(p,x,0)*(r);                          \
                           ND(p,z,1) = ND(p,x,1)*(r); }
#define MULT_S(p,r,x,z) {  ND(p,z,0) = ND(p,x,0)*r[0];                         \
                           ND(p,z,1) = ND(p,x,1)*r[1]; }
#define DIVIDE(p,x,y,z) {  ND(p,z,0) = ND(p,x,0)/ND(p,y,0);                    \
                           ND(p,z,1) = ND(p,x,1)/ND(p,y,1); } 
#define MULTIPLY(p,x,y,z) {ND(p,z,0) = ND(p,x,0)*ND(p,y,0);                    \
                           ND(p,z,1) = ND(p,x,1)*ND(p,y,1); } 
#define MSQRT(p,x,z)    {  ND(p,z,0) = sqrt(ND(p,x,0));                        \
                           ND(p,z,1) = sqrt(ND(p,x,1)); } 
#define ADD_VMULT(n,y,a,b)   { ND(n[0],y,0) += a[0]*b[0];                      \
                               ND(n[1],y,0) += a[1]*b[0];                      \
                               ND(n[2],y,0) += a[2]*b[0];                      \
                               ND(n[0],y,1) += a[0]*b[1];                      \
                               ND(n[1],y,1) += a[1]*b[1];                      \
                               ND(n[2],y,1) += a[2]*b[1]; }
#define SUBTR_VMULT(n,y,a,b) { ND(n[0],y,0) -= a[0]*b[0];                      \
                               ND(n[1],y,0) -= a[1]*b[0];                      \
                               ND(n[2],y,0) -= a[2]*b[0];                      \
                               ND(n[0],y,1) -= a[0]*b[1];                      \
                               ND(n[1],y,1) -= a[1]*b[1];                      \
                               ND(n[2],y,1) -= a[2]*b[1]; }
#define ADD_SUM(n,x,a,b)     { a[0] += b[0]*ND(n[0],x,0) +                     \
                                       b[1]*ND(n[1],x,0) +                     \
                                       b[2]*ND(n[2],x,0);                      \
                               a[1] += b[0]*ND(n[0],x,1) +                     \
                                       b[1]*ND(n[1],x,1) +                     \
                                       b[2]*ND(n[2],x,1); }
#define SUBTR_SUM(n,x,a,b)   { a[0] -= b[0]*ND(n[0],x,0) +                     \
                                       b[1]*ND(n[1],x,0) +                     \
                                       b[2]*ND(n[2],x,0);                      \
                               a[1] -= b[0]*ND(n[0],x,1) +                     \
                                       b[1]*ND(n[1],x,1) +                     \
                                       b[2]*ND(n[2],x,1); }
#define E_SUM(p,z)         (EDSN(p,z,0) + EDSN(p,z,1) + EDSN(p,z,2))
#define E_SET_VALUE(p,r,z)  EDSN(p,z,0) =    EDSN(p,z,1) =    EDSN(p,z,2) = r; 
#define E_ADD_VALUE(p,r,z) {EDSN(p,z,0)+= r; EDSN(p,z,1)+= r; EDSN(p,z,2)+= r;}
#define E_SUBTR_C(p,f,z) {  EDSN(p,z,0) -= f;                                  \
                            EDSN(p,z,1) -= f;                                  \
                            EDSN(p,z,2) -= f;  } 
#define E_COPY(p,x,z)    {  EDSN(p,z,0) = EDSN(p,x,0);                         \
                            EDSN(p,z,1) = EDSN(p,x,1);                         \
                            EDSN(p,z,2) = EDSN(p,x,2);  } 
#define E_INV(p,x,z)     {  EDSN(p,z,0) = -EDSN(p,x,0);                        \
                            EDSN(p,z,1) = -EDSN(p,x,1);                        \
                            EDSN(p,z,2) = -EDSN(p,x,2);  } 
#define E_ADD(p,x,y,z)   {  EDSN(p,z,0) = EDSN(p,x,0) + EDSN(p,y,0);           \
                            EDSN(p,z,1) = EDSN(p,x,1) + EDSN(p,y,1);           \
                            EDSN(p,z,2) = EDSN(p,x,2) + EDSN(p,y,2); } 
#define E_NSUBTR(p,x,y,z) { EDSN(p,z,0) = EDSN(p,x,0) - EDSN(p,y,0);           \
                            EDSN(p,z,1) = EDSN(p,x,1) - EDSN(p,y,1);           \
                            EDSN(p,z,2) = EDSN(p,x,2) - EDSN(p,y,2); }
#define E_NASUBTR(p,w,x,y,z) { EDSN(p,z,0) =                                   \
                               EDSN(p,w,0) + EDSN(p,x,0) - EDSN(p,y,0);        \
                               EDSN(p,z,1) =                                   \
                               EDSN(p,w,1) + EDSN(p,x,1) - EDSN(p,y,1);        \
                               EDSN(p,z,2) =                                   \
                               EDSN(p,w,2) + EDSN(p,x,2) - EDSN(p,y,2); }
#define E_MADD(p,r,x,y,z) {  EDSN(p,z,0) = (r)*EDSN(p,x,0) + EDSN(p,y,0);      \
                             EDSN(p,z,1) = (r)*EDSN(p,x,1) + EDSN(p,y,1);      \
                             EDSN(p,z,2) = (r)*EDSN(p,x,2) + EDSN(p,y,2);} 
#define E_ADD_MULT(p,r,x,z) {EDSN(p,z,0) += (r)*EDSN(p,x,0);                   \
                             EDSN(p,z,1) += (r)*EDSN(p,x,1);                   \
                             EDSN(p,z,2) += (r)*EDSN(p,x,2);} 
#define E_SUBTR_MULT(p,r,x,z) {EDSN(p,z,0) -= (r)*EDSN(p,x,0);                 \
                               EDSN(p,z,1) -= (r)*EDSN(p,x,1);                 \
                               EDSN(p,z,2) -= (r)*EDSN(p,x,2);} 
#define E_DAMP(p,r,x,y,z) {  EDSN(p,z,0) =                                     \
                             EDSN(p,x,0) + (r)*(EDSN(p,y,0) - EDSN(p,x,0));    \
                             EDSN(p,z,1) =                                     \
                             EDSN(p,x,1) + (r)*(EDSN(p,y,1) - EDSN(p,x,1));    \
                             EDSN(p,z,2) =                                     \
                             EDSN(p,x,2) + (r)*(EDSN(p,y,2) - EDSN(p,x,2));} 
#define E_MULT(p,r,x,z)   {  EDSN(p,z,0) = EDSN(p,x,0)*(r);                    \
                             EDSN(p,z,1) = EDSN(p,x,1)*(r);                    \
                             EDSN(p,z,2) = EDSN(p,x,2)*(r);  }
#define E_DIVIDE(p,x,y,z) {  EDSN(p,z,0) = EDSN(p,x,0)/EDSN(p,y,0);            \
                             EDSN(p,z,1) = EDSN(p,x,1)/EDSN(p,y,1);            \
                             EDSN(p,z,2) = EDSN(p,x,2)/EDSN(p,y,2);  }
#define E_MULTIPLY(p,x,y,z) {EDSN(p,z,0) = EDSN(p,x,0)*EDSN(p,y,0);            \
                             EDSN(p,z,1) = EDSN(p,x,1)*EDSN(p,y,1);            \
                             EDSN(p,z,2) = EDSN(p,x,2)*EDSN(p,y,2);  }
#define E_MSQRT(p,x,z)    {  EDSN(p,z,0) = sqrt(EDSN(p,x,0));                  \
                             EDSN(p,z,1) = sqrt(EDSN(p,x,1));                  \
                             EDSN(p,z,2) = sqrt(EDSN(p,x,2));  }
#define D_SET_VALUE(p,r,z) {FDDV(p,z,0,0) = FDDV(p,z,0,1) = r;                 \
                            FDDV(p,z,1,0) = FDDV(p,z,1,1) = r; }
#define D_COPY(p,x,z)    {  FDDV(p,z,0,0) = FDDV(p,x,0,0);                     \
                            FDDV(p,z,0,1) = FDDV(p,x,0,1);                     \
                            FDDV(p,z,1,0) = FDDV(p,x,1,0);                     \
                            FDDV(p,z,1,1) = FDDV(p,x,1,1);  }
#define D_INV(p,x,z)     {  FDDV(p,z,0,0) = -FDDV(p,x,0,0);                    \
                            FDDV(p,z,0,1) = -FDDV(p,x,0,1);                    \
                            FDDV(p,z,1,0) = -FDDV(p,x,1,0);                    \
                            FDDV(p,z,1,1) = -FDDV(p,x,1,1); }
#define D_ADD(p,x,y,z)   {  FDDV(p,z,0,0) = FDDV(p,x,0,0) + FDDV(p,y,0,0);     \
                            FDDV(p,z,0,1) = FDDV(p,x,0,1) + FDDV(p,y,0,1);     \
                            FDDV(p,z,1,0) = FDDV(p,x,1,0) + FDDV(p,y,1,0);     \
                            FDDV(p,z,1,1) = FDDV(p,x,1,1) + FDDV(p,y,1,1); }
#define D_NSUBTR(p,x,y,z) {  FDDV(p,z,0,0) = FDDV(p,x,0,0) - FDDV(p,y,0,0);    \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1) - FDDV(p,y,0,1);    \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0) - FDDV(p,y,1,0);    \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1) - FDDV(p,y,1,1); }
#define D_NASUBTR(p,w,x,y,z) { FDDV(p,z,0,0) =                                 \
                               FDDV(p,w,0,0) + FDDV(p,x,0,0) - FDDV(p,y,0,0);  \
                               FDDV(p,z,0,1) =                                 \
                               FDDV(p,w,0,1) + FDDV(p,x,0,1) - FDDV(p,y,0,1);  \
                               FDDV(p,z,1,0) =                                 \
                               FDDV(p,w,1,0) + FDDV(p,x,1,0) - FDDV(p,y,1,0);  \
                               FDDV(p,z,1,1) =                                 \
                               FDDV(p,w,1,1) + FDDV(p,x,1,1) - FDDV(p,y,1,1); }
#define D_MADD(p,r,x,y,z) {  FDDV(p,z,0,0) = (r)*FDDV(p,x,0,0) + FDDV(p,y,0,0);\
                             FDDV(p,z,0,1) = (r)*FDDV(p,x,0,1) + FDDV(p,y,0,1);\
                             FDDV(p,z,1,0) = (r)*FDDV(p,x,1,0) + FDDV(p,y,1,0);\
                             FDDV(p,z,1,1) = (r)*FDDV(p,x,1,1) + FDDV(p,y,1,1);}
#define D_ADD_MULT(p,r,x,z) {FDDV(p,z,0,0) += (r)*FDDV(p,x,0,0);               \
                             FDDV(p,z,0,1) += (r)*FDDV(p,x,0,1);               \
                             FDDV(p,z,1,0) += (r)*FDDV(p,x,1,0);               \
                             FDDV(p,z,1,1) += (r)*FDDV(p,x,1,1); }
#define D_SUBTR_MULT(p,r,x,z) {FDDV(p,z,0,0) -= (r)*FDDV(p,x,0,0);             \
                               FDDV(p,z,0,1) -= (r)*FDDV(p,x,0,1);             \
                               FDDV(p,z,1,0) -= (r)*FDDV(p,x,1,0);             \
                               FDDV(p,z,1,1) -= (r)*FDDV(p,x,1,1);}
#define D_DAMP(p,r,x,y,z) {FDDV(p,z,0,0) =                                     \
                           FDDV(p,x,0,0) + (r)*(FDDV(p,y,0,0) - FDDV(p,x,0,0));\
                           FDDV(p,z,0,1) =                                     \
                           FDDV(p,x,0,1) + (r)*(FDDV(p,y,0,1) - FDDV(p,x,0,1));\
                           FDDV(p,z,1,0) =                                     \
                           FDDV(p,x,1,0) + (r)*(FDDV(p,y,1,0) - FDDV(p,x,1,0));\
                           FDDV(p,z,1,1) =                                     \
                           FDDV(p,x,1,1) + (r)*(FDDV(p,y,1,1) - FDDV(p,x,1,1));}
#define D_MULT(p,r,x,z)   {  FDDV(p,z,0,0) = FDDV(p,x,0,0)*(r);                \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)*(r);                \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)*(r);                \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)*(r); }
#define D_DIVIDE(p,x,y,z) {  FDDV(p,z,0,0) = FDDV(p,x,0,0)/FDDV(p,y,0,0);      \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)/FDDV(p,y,0,1);      \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)/FDDV(p,y,1,0);      \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)/FDDV(p,y,1,1); } 
#define D_MULTIPLY(p,x,y,z) {FDDV(p,z,0,0) = FDDV(p,x,0,0)*FDDV(p,y,0,0);      \
                             FDDV(p,z,0,1) = FDDV(p,x,0,1)*FDDV(p,y,0,1);      \
                             FDDV(p,z,1,0) = FDDV(p,x,1,0)*FDDV(p,y,1,0);      \
                             FDDV(p,z,1,1) = FDDV(p,x,1,1)*FDDV(p,y,1,1); } 
#define D_MSQRT(p,x,z)    {  FDDV(p,z,0,0) = sqrt(FDDV(p,x,0,0));              \
                             FDDV(p,z,0,1) = sqrt(FDDV(p,x,0,1));              \
                             FDDV(p,z,1,0) = sqrt(FDDV(p,x,1,0));              \
                             FDDV(p,z,1,1) = sqrt(FDDV(p,x,1,1));  }
#define SDIVIDE(p,x,y,z) { ND(p,z,0) = ND(p,x,0)/NDS(p,y);                     \
                           ND(p,z,1) = ND(p,x,1)/NDS(p,y); } 
#define SMULTIPLY(p,x,y,z) { ND(p,z,0) = ND(p,x,0)*NDS(p,y);                   \
                             ND(p,z,1) = ND(p,x,1)*NDS(p,y); } 

#define SET_COEFFNN(p,k,a,b,c,d)   COEFFNN(p,k,0,0) += a;                      \
                                   COEFFNN(p,k,0,1) += b;                      \
                                   COEFFNN(p,k,1,0) += c;                      \
                                   COEFFNN(p,k,1,1) += d
#define SET_COEFFNN2(p,k,i0,i1,a,b,c,d)   COEFFNN(p,k,i0,i0) += a;             \
                                          COEFFNN(p,k,i0,i1) += b;             \
                                          COEFFNN(p,k,i1,i0) += c;             \
                                          COEFFNN(p,k,i1,i1) += d
#define LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n) {if (IS_FN(n0)) in0 = n++;  \
                                                    if (IS_FN(n1)) in1 = n++;  \
                                                    if (IS_FN(n2)) in2 = n++;}
#define LOCAL_FACE_INDICES(f0,f1,f2,if0,if1,if2,n) {if (IS_FF(f0)) if0 = n++;  \
                                                    if (IS_FF(f1)) if1 = n++;  \
                                                    if (IS_FF(f2)) if2 = n++;}
#define FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)    a[i][i] = COEFFN(n0,ZA);   \
                    for (pli = START(n0); pli; pli = NEXT(pli)){               \
                       if (NBNODE(pli) == n1)      a[i][in1] = COEFFL(pli,ZA); \
                       else if (NBNODE(pli) == n2) a[i][in2] = COEFFL(pli,ZA); \
                       else SET9(u_rhs[i],NDD(NBNODE(pli),u),COEFFL(pli,ZA)) }
#define FILL_AND_SUBTR_NF_MATR(n0)                                             \
                for (pnfl = NFSTART(n0); pnfl; pnfl = NEXT(pnfl)){             \
                   if (NBFACE(pnfl) == fa0)      a[i][ifa0] = COEFFL(pnfl,ZA); \
                   else if (NBFACE(pnfl) == fa1) a[i][ifa1] = COEFFL(pnfl,ZA); \
                   else if (NBFACE(pnfl) == fa2) a[i][ifa2] = COEFFL(pnfl,ZA); \
                   else SET9(u_rhs[i],FDVP(NBFACE(pnfl),u),COEFFL(pnfl,ZA)) }
#define FILL_AND_SUBTR_FN_MATR(fa0)                                            \
                  for (pfnl = FNSTART(fa0); pfnl; pfnl = NEXT(pfnl)){          \
                     if (NBNODE(pfnl) == n0)      a[i][in0] = COEFFL(pfnl,ZA); \
                     else if (NBNODE(pfnl) == n1) a[i][in1] = COEFFL(pfnl,ZA); \
                     else if (NBNODE(pfnl) == n2) a[i][in2] = COEFFL(pfnl,ZA); \
                     else SET9(u_rhs[i],NDD(NBNODE(pfnl),u),COEFFL(pfnl,ZA)) }
#define FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)                          \
                  a[i][i] = COEFF_FF(fa0,ZA);                                  \
                  for (pfl = FSTART(fa0); pfl; pfl = NEXT(pfl)){               \
                     if (NBFACE(pfl) == fa1)      a[i][ifa1] = COEFFL(pfl,ZA); \
                     else if (NBFACE(pfl) == fa2) a[i][ifa2] = COEFFL(pfl,ZA); \
                     else SET9(u_rhs[i],FDVP(NBFACE(pfl),u),COEFFL(pfl,ZA)) }
#define FILL_P1DISC_B_MATR(p,i,j,macro) { b0[0][i] = macro(p,ZB,0,j,0);        \
                                          b0[1][i] = macro(p,ZB,1,j,0);        \
                                          b0[2][i] = macro(p,ZB,2,j,0);        \
                                          b1[0][i] = macro(p,ZB,0,j,1);        \
                                          b1[1][i] = macro(p,ZB,1,j,1);        \
                                          b1[2][i] = macro(p,ZB,2,j,1); }
#define SUBTR_P1DISC_B_MATR(a,q,i,macro) { SET9(a,macro(q,ZB,0,i),EDSN(q,p,0)) \
                                           SET9(a,macro(q,ZB,1,i),EDSN(q,p,1)) \
                                           SET9(a,macro(q,ZB,2,i),EDSN(q,p,2)) }
#define COMPUTE_MAX_ND_DIFF                                                    \
             if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][in0]),               \
                                     fabs(ND(n0,u,1)-y[1][in0]))) > s)  s = q; \
             if (IS_FN(n1) && (q=MAX(fabs(ND(n1,u,0)-y[0][in1]),               \
                                     fabs(ND(n1,u,1)-y[1][in1]))) > s)  s = q; \
             if (IS_FN(n2) && (q=MAX(fabs(ND(n2,u,0)-y[0][in2]),               \
                                     fabs(ND(n2,u,1)-y[1][in2]))) > s)  s = q;
#define COMPUTE_MAX_FDV_DIFF                                                   \
         if (IS_FF(fa0) && (q=MAX(fabs(FDV(fa0,u,0)-y[0][ifa0]),               \
                                  fabs(FDV(fa0,u,1)-y[1][ifa0]))) > s)  s = q; \
         if (IS_FF(fa1) && (q=MAX(fabs(FDV(fa1,u,0)-y[0][ifa1]),               \
                                  fabs(FDV(fa1,u,1)-y[1][ifa1]))) > s)  s = q; \
         if (IS_FF(fa2) && (q=MAX(fabs(FDV(fa2,u,0)-y[0][ifa2]),               \
                                  fabs(FDV(fa2,u,1)-y[1][ifa2]))) > s)  s = q;
#define UPDATE_FDV_VALUES                                                      \
         if (IS_FF(fa0)){                                                      \
            FDV(fa0,u_new,0) = FDV(fa0,u,0) + om*(y[0][ifa0] - FDV(fa0,u,0));  \
            FDV(fa0,u_new,1) = FDV(fa0,u,1) + om*(y[1][ifa0] - FDV(fa0,u,1));} \
         if (IS_FF(fa1)){                                                      \
            FDV(fa1,u_new,0) = FDV(fa1,u,0) + om*(y[0][ifa1] - FDV(fa1,u,0));  \
            FDV(fa1,u_new,1) = FDV(fa1,u,1) + om*(y[1][ifa1] - FDV(fa1,u,1));} \
         if (IS_FF(fa2)){                                                      \
            FDV(fa2,u_new,0) = FDV(fa2,u,0) + om*(y[0][ifa2] - FDV(fa2,u,0));  \
            FDV(fa2,u_new,1) = FDV(fa2,u,1) + om*(y[1][ifa2] - FDV(fa2,u,1));}
#define UPDATE_ND_VALUES                                                       \
         if (IS_FN(n0)){                                                       \
            ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][in0] - ND(n0,u,0));         \
            ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][in0] - ND(n0,u,1)); }       \
         if (IS_FN(n1)){                                                       \
            ND(n1,u_new,0) = ND(n1,u,0) + om*(y[0][in1] - ND(n1,u,0));         \
            ND(n1,u_new,1) = ND(n1,u,1) + om*(y[1][in1] - ND(n1,u,1)); }       \
         if (IS_FN(n2)){                                                       \
            ND(n2,u_new,0) = ND(n2,u,0) + om*(y[0][in2] - ND(n2,u,0));         \
            ND(n2,u_new,1) = ND(n2,u,1) + om*(y[1][in2] - ND(n2,u,1)); }
#endif

#define GDOTS(a,b,s,i,m)     { for(i=0; i < (m); i++) s += a[i]*b[i]; }
#define GSET1(a,b,i,m)       { for(i=0; i < (m); i++) a[i] = b[i]; }
#define GSET2(a,b,p,i,m)     { for(i=0; i < (m); i++) a[i] = b[i]*(p); }
#define GSET4(a,b,p,i,m)     { for(i=0; i < (m); i++) a[i] += b[i]*(p); }
#define GSET7(a,p,i,m)       { for(i=0; i < (m); i++) a[i] = p; }
#define GSET8(a,b,i,m)       { for(i=0; i < (m); i++) a[i] = -b[i]; }
#define GSET9(a,b,p,i,m)     { for(i=0; i < (m); i++) a[i] -= b[i]*(p); }
#define GSET10(a,b,c,i,m)    { for(i=0; i < (m); i++) a[i] = b[i] + c[i]; }
#define GSET11(a,b,c,i,m)    { for(i=0; i < (m); i++) a[i] = b[i] - c[i]; }
#define GSET23(a,b,c,p,i,m)  { for(i=0; i < (m); i++) a[i] =  b[i] + c[i]*(p); }
#define GSET24(a,b,c,p,i,m)  { for(i=0; i < (m); i++) a[i] = -b[i] + c[i]*(p); }
#define GSET30(a,b,c,d,i,m)  { for(i=0; i < (m); i++) a[i]=b[i]+c[i]-d[i]; }
#define GSET31(a,b,c,r,i,m)  { for(i=0; i < (m); i++)                          \
                                                   a[i]=b[i]+(r)*(c[i]-b[i]); }
#define GSET32(a,b,c,i,m)    { for(i=0; i < (m); i++) a[i] = b[i]/c[i]; }
#define GSET33(a,b,c,i,m)    { for(i=0; i < (m); i++) a[i] = b[i]*c[i]; }
#define GSET34(a,b,i,m)      { for(i=0; i < (m); i++) a[i] = sqrt(b[i]); }
#define MMSET1(a,b,i,j,m,n)  { for(i=0; i < (m); i++)                       \
                               for(j=0; j < (n); j++) a[i][j] = b[i][j]; }
#define MMSET2(a,b,p,i,j,m,n)   { for(i=0; i < (m); i++){  a[i] = 0;           \
                            for(j=0; j < (n); j++) a[i] += b[i][j]*p[j]; } }
#define MMSET4(a,b,p,i,j,m,n)   { for(i=0; i < (m); i++)                       \
                                  for(j=0; j < (n); j++) a[i] += b[i][j]*p[j]; }
#define MMSET7(a,p,i,j,m,n)     { for(i=0; i < (m); i++)                       \
                                  for(j=0; j < (n); j++) a[i][j] = p; }
#define MMSET9(a,b,p,i,j,m,n)   { for(i=0; i < (m); i++)                       \
                                  for(j=0; j < (n); j++) a[i] -= b[i][j]*p[j]; }

#define NODE_INDEX(p,pn,i)   { for (i = 0; i < NVERT && (p)->n[i] != pn; i++); \
                               if ((p)->n[i] != pn) i = -1;  }
#define FACE_INDEX(p,pf,i)   { for (i = 0; i < SIDES && (p)->f[i] != pf; i++); \
                               if ((p)->f[i] != pf) i = -1;  }

#define LINV(a,x) (DOT(a,x)+a[DIM])
#define UNIT_VECTOR(a,b,p)  {p = sqrt(DOT(b,b)); SET22(a,b,p)}
#define ORT_VECT_DIR(a,b,c) { ORT_VECT(a,b) if (DOT(a,c) < 0.) SET8(a,a) }
#define CHECK_TAU(r,r1,a,amin,hK) { if (fabs(r) < 1.e-20) r = -1.;             \
                  else{  r = r1; if (r > 0.){ a *= amin/hK; a = MAX(a,r); } } }
#define SET_TAU(pel,a,r)   { if (a < 1.e-20 || r < 0.)  pel->tau = 0.;         \
                             else                       pel->tau = 1./a; }

#define ADD_TO_FILE1(fp,name,text)                                             \
                 {fp = fopen(name,"a"); fprintf(fp,text);         fclose(fp);}
#define ADD_TO_FILE2(fp,name,text,m)                                           \
                 {fp = fopen(name,"a"); fprintf(fp,text,m);       fclose(fp);}
#define ADD_TO_FILE3(fp,name,text,m,n)                                         \
                 {fp = fopen(name,"a"); fprintf(fp,text,m,n);     fclose(fp);}
#define ADD_TO_FILE4(fp,name,text,m,n,o)                                       \
                 {fp = fopen(name,"a"); fprintf(fp,text,m,n,o);   fclose(fp);}
#define ADD_TO_FILE5(fp,name,text,m,n,o,p)                                     \
                 {fp = fopen(name,"a"); fprintf(fp,text,m,n,o,p); fclose(fp);}

#define ADD_NEW_LINK(start,link,pli)  {if (start == NULL)                      \
                                          start = link;                        \
                                       else{                                   \
                                          for (pli = start; pli->next != NULL; \
                                                            pli=pli->next);    \
                                          pli->next = link;                    \
                                       }}

#define ADD_FACE(nf)   { if(IS_FF(nf)) ADD_NEW_LINK(n->nfstart, *pnflink,pnf)  \
                         else          ADD_NEW_LINK(n->tnfstart,*pnflink,pnf)  \
                         (*pnflink)->nbface = nf;                              \
                         ((*pnflink)++)->next = NULL; }
                                        
#define ADD_FIRST_FACE1(nf) NFLINK *pnf;                                       \
                            ADD_NEW_LINK(n->nfstart,*pnflink,pnf)              \
                            (*pnflink)->nbface = nf;
                            
#define ADD_FIRST_FACE2(nf) NFLINK *pnf;                                       \
                            for (pnf = n->nfstart; pnf->next != NULL;          \
                                                               pnf=pnf->next); \
                            pnf->next = *pnflink; (*pnflink)->nbface = nf; 
#define ADD_NEXT_FACE(nf) { (*pnflink)->next = *pnflink + 1;                   \
                            (++(*pnflink))->nbface = nf; }

#if DATA_S & N_LINK_TO_FACES
#define SET_NF for (pnflink=pnode->tnfstart;pnflink!=NULL;pnflink=pnflink->next)\
                  SET7(COEFF_NFP(pnflink,Z),r)
#else
#define SET_NF
#endif

#if DATA_S & F_LINK_TO_FACES
#define SET_FF   for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)  \
                    COEFF_FL(pflink,Z) = r;
#else
#define SET_FF
#endif

#if DATA_S & F_LINK_TO_NODES
#define SET_FN for (pfnlink=pface->tfnstart;pfnlink!=NULL;pfnlink=pfnlink->next)\
                  SET7(COEFF_FNP(pfnlink,Z),r)
#else
#define SET_FN
#endif

#define GET_ITYPE_INDEX(i,itype,it) for (i = 0; it[i] != itype; i++);

FLOAT lin_coeff0_0, lin_coeff0_1, lin_coeff0_2, 
      lin_coeff1_0, lin_coeff1_1, lin_coeff1_2;

double lin_func0(x)
double x[2];
{
   return(lin_coeff0_0*x[0] + lin_coeff0_1*x[1] + lin_coeff0_2);
}

double lin_func1(x)
double x[2];
{
   return(lin_coeff1_0*x[0] + lin_coeff1_1*x[1] + lin_coeff1_2);
}

double zero_func(x)
double x[2];
{
   return(0.);
}

/*   principal lattices  */

#define F1_6 0.166666666666666666666666666666666666666666
#define F1_3 0.333333333333333333333333333333333333333333
#define F2_3 0.666666666666666666666666666666666666666666
#define F5_6 0.833333333333333333333333333333333333333333

double LP1[ 3][DIM]={{0.,0.}, {1.,0.}, {0.,1.}};
double LP2[ 6][DIM]={{0.,0.}, {1.,0.}, {0.,1.}, {.5,.5}, {.5,0.}, {0.,.5}};
double LP3[10][DIM]={{0.,0.  }, {F1_3,0.  }, {F2_3,0.  }, {1.,0.}, 
                     {0.,F1_3}, {F1_3,F1_3}, {F2_3,F1_3}, 
                     {0.,F2_3}, {F1_3,F2_3}, 
                     {0.,1.  }};
double LP4[15][DIM]={{0.,0. }, {.25,0. }, {.5,0. }, {.75,0. },  {1.,0.},
                     {0.,.25}, {.25,.25}, {.5,.25}, {.75,.25},
                     {0.,.5 }, {.25,.5 }, {.5,.5 }, 
                     {0.,.75}, {.25,.75}, 
                     {0.,1.  }};
double LP5[21][DIM]={{0.,0.}, {.2,0.}, {.4,0.}, {.6,0.}, {.8,0.}, {1.,0.},
                     {0.,.2}, {.2,.2}, {.4,.2}, {.6,.2}, {.8,.2},
                     {0.,.4}, {.2,.4}, {.4,.4}, {.6,.4},
                     {0.,.6}, {.2,.6}, {.4,.6}, 
                     {0.,.8}, {.2,.8},
                     {0.,1.}};
double LP6[28][DIM]={{0.,0.  }, {F1_6,0.  }, {F1_3,0.  }, {.5,0.  }, {F2_3,0.  }, {F5_6,0.  }, {1.,0.},
                     {0.,F1_6}, {F1_6,F1_6}, {F1_3,F1_6}, {.5,F1_6}, {F2_3,F1_6}, {F5_6,F1_6},
                     {0.,F1_3}, {F1_6,F1_3}, {F1_3,F1_3}, {.5,F1_3}, {F2_3,F1_3},
                     {0.,0.5 }, {F1_6,0.5 }, {F1_3,0.5 }, {.5,0.5 },
                     {0.,F2_3}, {F1_6,F2_3}, {F1_3,F2_3},
                     {0.,F5_6}, {F1_6,F5_6},
                     {0.,1.  }};

double LQ1[ 4][DIM]={{-1.,-1.}, { 1.,-1.}, {-1., 1.}, { 1., 1.}};
double LQ2[ 9][DIM]={{-1.,-1.}, {0.,-1.}, {1.,-1.},
                     {-1., 0.}, {0., 0.}, {1., 0.},
                     {-1., 1.}, {0., 1.}, {1., 1.}};
double LQ3[16][DIM]={{-1.,-1.  }, {-F1_3,-1.  }, {F1_3,-1.  }, {1.,-1.}, 
                     {-1.,-F1_3}, {-F1_3,-F1_3}, {F1_3,-F1_3}, {1.,-F1_3},
                     {-1., F1_3}, {-F1_3, F1_3}, {F1_3, F1_3}, {1., F1_3},
                     {-1., 1.  }, {-F1_3, 1.  }, {F1_3, 1.  }, {1., 1.  }};
double LQ4[25][DIM]={{-1.,-1.}, {-.5,-1.}, {0.,-1.}, {.5,-1.},  {1.,-1.},
                     {-1.,-.5}, {-.5,-.5}, {0.,-.5}, {.5,-.5},  {1.,-.5},
                     {-1., 0.}, {-.5, 0.}, {0., 0.}, {.5, 0.},  {1., 0.},
                     {-1., .5}, {-.5, .5}, {0., .5}, {.5, .5},  {1., .5},
                     {-1., 1.}, {-.5, 1.}, {0., 1.}, {.5, 1.},  {1., 1.}};
double LQ5[36][DIM]={{-1.,-1.}, {-.6,-1.}, {-.2,-1.}, { .2,-1.}, { .6,-1.}, {1.,-1.},
                     {-1.,-.6}, {-.6,-.6}, {-.2,-.6}, { .2,-.6}, { .6,-.6}, {1.,-.6},
                     {-1.,-.2}, {-.6,-.2}, {-.2,-.2}, { .2,-.2}, { .6,-.2}, {1.,-.2},
                     {-1., .2}, {-.6, .2}, {-.2, .2}, { .2, .2}, { .6, .2}, {1., .2},
                     {-1., .6}, {-.6, .6}, {-.2, .6}, { .2, .6}, { .6, .6}, {1., .6},
                     {-1., 1.}, {-.6, 1.}, {-.2, 1.}, { .2, 1.}, { .6, 1.}, {1., 1.}};
double LQ6[49][DIM]={{-1.,-1.  }, {-F2_3,-1.  }, {-F1_3,-1.  }, {0.,-1.  }, {F1_3,-1.  }, {F2_3,-1.  }, {1.,-1.  },
                     {-1.,-F2_3}, {-F2_3,-F2_3}, {-F1_3,-F2_3}, {0.,-F2_3}, {F1_3,-F2_3}, {F2_3,-F2_3}, {1.,-F2_3},
                     {-1.,-F1_3}, {-F2_3,-F1_3}, {-F1_3,-F1_3}, {0.,-F1_3}, {F1_3,-F1_3}, {F2_3,-F1_3}, {1.,-F1_3},
                     {-1., 0.  }, {-F2_3, 0.  }, {-F1_3, 0.  }, {0., 0.  }, {F1_3, 0.  }, {F2_3, 0.  }, {1., 0.  },
                     {-1., F1_3}, {-F2_3, F1_3}, {-F1_3, F1_3}, {0., F1_3}, {F1_3, F1_3}, {F2_3, F1_3}, {1., F1_3},
                     {-1., F2_3}, {-F2_3, F2_3}, {-F1_3, F2_3}, {0., F2_3}, {F1_3, F2_3}, {F2_3, F2_3}, {1., F2_3},
                     {-1., 1.  }, {-F2_3, 1.  }, {-F1_3, 1.  }, {0., 1.  }, {F1_3, 1.  }, {F2_3, 1.  }, {1., 1.  }};

double LP1S[ 4][DIM]={{0.,0.}, {1.,0.}, {0.,1.}, {F1_3,F1_3}};
double LP2S[10][DIM]={{0.,0.}, {1.,0.}, {0.,1.}, {.5,.5}, {.5,0.}, {0.,.5},
                      {F1_6,F1_6}, {F2_3,F1_6}, {F1_6,F2_3}, {F1_3,F1_3}};

#define ALLOCATED_MEMORY    (MAXNODE*sizeof(NODE) + MAXFACE*sizeof(FACE)       \
    + MAXVERT*sizeof(VERTEX) + MAXBPOINT*sizeof(BPOINT)                        \
    + MAXELEM*sizeof(ELEMENT) + MAXLINK*sizeof(LINK) + MAXELINK*sizeof(ELINK)  \
    + MAXFNLINK*sizeof(FNLINK) + MAXNFLINK*sizeof(NFLINK)                      \
    + MAXFLINK*sizeof(FLINK) + MAXSNODE*sizeof(SNODE)                          \
    + MAXSFACE*sizeof(SFACE) + MAXNLGLINK*sizeof(NLGLINK)                      \
    + MAXLGNLINK*sizeof(LGNLINK) + MAXLGFLINK*sizeof(LGFLINK)                  \
    + MAXLGLGLINK*sizeof(LGLGLINK) + MAXFLGLINK*sizeof(FLGLINK)                \
    + MAXLGDATA*sizeof(LGDATA) + MAXLEVEL*sizeof(GRID))

