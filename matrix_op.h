/*

    POZOR NA VOLBU q!!!!!!!!
    NESMI SE SHODOVAT S JINYM ARGUMENTEM!!!!!!!!!!!

*/


/******************************************************************************/
/*                                                                            */
/*                             matrix operations                              */
/*                                                                            */
/******************************************************************************/

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) 

void sset_mat_value(tGrid,Z,r)  /*  Z := r */
GRID *tGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   LINK *plink;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = SUCC(pnode)){
         COEFFS(pnode,Z) = r;
         for (plink = TSTART(pnode); plink != NULL; plink = NEXT(plink))
               COEFFLS(plink,Z) = r;
   }
}

#else

void sset_mat_value(tGrid,Z,r)  /*  Z := r */
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA)

void smultA(tGrid,Z,x,y,t)  /* y := Z x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid;
INT Z, x, y, t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   LINK *pl, *stop;
   DOUBLE sum;
	
   if (t & USE_IS_DIR){
      for (theNode=FIRSTN(tGrid); theNode; theNode=SUCC(theNode))
         if (!IS_DIR(theNode,t)){
            sum = COEFFS(theNode,Z)*NDS(theNode,x);
            for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
               pnode = NBNODE(pl);
               if (!IS_DIR(pnode,t))
	          sum += COEFFLS(pl,Z)*NDS(pnode,x);
            }
            NDS(theNode,y) = sum;
         }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         for (theNode=FIRSTN(theGrid); theNode; theNode=SUCC(theNode)) 
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t)){
               sum = COEFFS(theNode,Z)*NDS(theNode,x);
               for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
                  pnode = ltop_node(NBNODE(pl),tGrid);
                  if (!IS_DIR(pnode,t))
	             sum += COEFFLS(pl,Z)*NDS(pnode,x);
               }
               NDS(theNode,y) = sum;
            }
      }
   }
   else{
      for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
         sum = COEFFS(theNode,Z)*NDS(theNode,x);
         stop = STOP_NN(theNode,t);
         for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
            sum += COEFFLS(pl,Z)*NDS(NBNODE(pl),x);
         NDS(theNode,y) = sum;
      }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         for (theNode=FIRSTNODE(theGrid); theNode; theNode=SUCC(theNode)) 
            if (IS_LTOP_NODE(theNode,tGrid)){
               sum = COEFFS(theNode,Z)*NDS(theNode,x);
               stop = STOP_NN(theNode,t);
               for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
	          sum += COEFFLS(pl,Z)*NDS(ltop_node(NBNODE(pl),tGrid),x);
               NDS(theNode,y) = sum;
            }
      }
   }
}

void smultAT(tGrid,Z,x,y,t)  /* y = Z^T x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid;
INT Z, x, y, t;
{
   NODE *theNode, *pnode;
   LINK *pl, *pli;
   DOUBLE sum;
	
   for (theNode=FIRSTN(tGrid); theNode!= NULL; theNode=SUCC(theNode))
      if (!IS_DIR(theNode,t)){
         sum = COEFFS(theNode,Z)*NDS(theNode,x);
         for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
           pnode = NBNODE(pl);
           if (!IS_DIR(pnode,t)){
             for (pli = TSTART(pnode); pli->nbnode != theNode; pli = pli->next);
	     sum += COEFFLS(pli,Z)*NDS(pnode,x);
	   }
         }
         NDS(theNode,y) = sum;
      }
}

void sSOR_step_forward_sn(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   NODE *theNode;
   LINK *pl;
   DOUBLE sum;
	
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      if (fabs(COEFFS(theNode,Z)) > EPSA){
         sum = NDS(theNode,b);
         for (pl=START(theNode); pl; pl=NEXT(pl))
            sum -= COEFFLS(pl,Z)*NDS(NBNODE(pl),x);
         NDS(theNode,x) += om*(sum/COEFFS(theNode,Z)-NDS(theNode,x));
      }
}

void mso_smultA(tGrid,Z,x,y,t)  /* y := Z x; */
GRID *tGrid;
INT Z, x, y, t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   LINK *pl, *stop;
   DOUBLE sum;
	
   if (t & USE_IS_DIR){
      for (theNode=FIRSTN(tGrid); theNode; theNode=SUCC(theNode))
         if (!IS_DIR(theNode,t)){
            sum = COEFFS(theNode,Z)*NDS(theNode,x);
            for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
               pnode = NBNODE(pl);
			   if (t & STOP_IS_FIRST_INNER) {
				  if (IS_DIR(pnode,t))
					 sum += COEFFLS(pl,Z)*NDS(pnode,x);
			   } else {
				  if (!IS_DIR(pnode,t))
					 sum += COEFFLS(pl,Z)*NDS(pnode,x);
			   }
            }
            NDS(theNode,y) = sum;
         }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         for (theNode=FIRSTN(theGrid); theNode; theNode=SUCC(theNode)) 
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t)){
               sum = COEFFS(theNode,Z)*NDS(theNode,x);
               for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
                  pnode = ltop_node(NBNODE(pl),tGrid);
                  if (t & STOP_IS_FIRST_INNER) {
  				     if (IS_DIR(pnode,t))
					    sum += COEFFLS(pl,Z)*NDS(pnode,x);
				  } else {
				     if (!IS_DIR(pnode,t))
					    sum += COEFFLS(pl,Z)*NDS(pnode,x);
				  }
               }
               NDS(theNode,y) = sum;
            }
      }
   }
   else{
      for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
         sum = COEFFS(theNode,Z)*NDS(theNode,x);
         stop = STOP_NN(theNode,t);
         for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
            sum += COEFFLS(pl,Z)*NDS(NBNODE(pl),x);
         NDS(theNode,y) = sum;
      }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         for (theNode=FIRSTNODE(theGrid); theNode; theNode=SUCC(theNode)) 
            if (IS_LTOP_NODE(theNode,tGrid)){
               sum = COEFFS(theNode,Z)*NDS(theNode,x);
               stop = STOP_NN(theNode,t);
               for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
	          sum += COEFFLS(pl,Z)*NDS(ltop_node(NBNODE(pl),tGrid),x);
               NDS(theNode,y) = sum;
            }
      }
   }
}


#else  /*  if !((N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && 
                (N_DATA & SCALAR_NODE_DATA))  */

void smultA(tGrid,Z,x,y,t)  /* y := Z x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: smultA not available.\n");  }

void smultAT(tGrid,Z,x,y,t)  /* y = Z^T x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: smultAT not available.\n");  }

void sSOR_step_forward_sn(tGrid,Z,b,x,om)
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: sSOR_step_forward_sn not available.\n");  }

void mso_smultA(tGrid,Z,x,y,t)  /* y := Z x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: mso_smultA not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES)

FLOAT snorm_of_A(tGrid,Z,t)  /* := |Z|;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid;
INT Z, t;
{
   GRID *theGrid;
   NODE *pnode;
   LINK *plink;
   FLOAT sum, max=0.0;
   
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser)
      for (pnode = FIRSTN(theGrid); pnode != NULL; pnode = pnode->succ)
         if (IS_LTOP_NODE(pnode,tGrid) && !IS_DIR(pnode,t)){
            sum = fabs(COEFFS(pnode,Z));
            for (plink = TSTART(pnode); plink != NULL; plink = plink->next)
               if (!IS_DIR(NBNODE(plink),t))
                  sum += fabs(COEFFS(plink,Z));
            if (sum > max) max = sum;
         }
   return(max);
}

void smakeAT(tGrid,Z,a,t)  /* a := Z^T ;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid;
INT Z, a, t;
{
   NODE *theNode, *pnode;
   LINK *pl, *pli;

   for (theNode=FIRSTN(tGrid); theNode!= NULL; theNode=SUCC(theNode))
      if (!IS_DIR(theNode,t)){
         COEFFS(theNode,a) = COEFFS(theNode,Z);
         for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
           pnode = NBNODE(pl);
           if (!IS_DIR(pnode,t)){
             for (pli = TSTART(pnode); pli->nbnode != theNode; pli = pli->next);
             COEFFLS(pl,a) = COEFFLS(pli,Z);
           }
         }
      }
}

void sadd_matr(mg,Z,a,t)  /*  Z := Z + A;  t=1 ... temp., t=2 ... conc. */
MULTIGRID *mg;
INT Z, a, t;
{
   NODE *pnode;
   LINK *plink;

   for (pnode=FIRSTN(TOP_GRID(mg)); pnode!= NULL; pnode=SUCC(pnode))
      if (!IS_DIR(pnode,t)){
         COEFFS(pnode,Z) += COEFFS(pnode,a);
         for (plink=TSTART(pnode); plink!=NULL; plink=NEXT(plink))
            if (!IS_DIR(NBNODE(plink),t))
               COEFFLS(plink,Z) += COEFFLS(plink,a);
      }
} 
 
void ssubtr_matr(mg,Z,a,t)  /*  Z := Z - A;  t=1 ... temp., t=2 ... conc. */
MULTIGRID *mg;
INT Z, a, t;
{
   NODE *pnode;
   LINK *plink;

   for (pnode=FIRSTN(TOP_GRID(mg)); pnode!= NULL; pnode=SUCC(pnode))
      if (!IS_DIR(pnode,t)){
         COEFFS(pnode,Z) -= COEFFS(pnode,a);
         for (plink=TSTART(pnode); plink!=NULL; plink=NEXT(plink))
            if (!IS_DIR(NBNODE(plink),t))
               COEFFLS(plink,Z) -= COEFFLS(plink,a);
      }
} 

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES))       */

FLOAT snorm_of_A(tGrid,Z,t)  /* := |Z|;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT Z, t;
{  eprintf("Error: snorm_of_A not available.\n"); return(0.);  }

void smakeAT(tGrid,Z,a,t)  /* a := Z^T ;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT Z, a, t;
{  eprintf("Error: smakeAT not available.\n");  }

void sadd_matr(mg,Z,a,t)  /*  Z := Z + A;  t=1 ... temp., t=2 ... conc. */
MULTIGRID *mg; INT Z, a, t;
{  eprintf("Error: sadd_matr not available.\n");  }
 
void ssubtr_matr(mg,Z,a,t)  /*  Z := Z - A;  t=1 ... temp., t=2 ... conc. */
MULTIGRID *mg; INT Z, a, t;
{  eprintf("Error: ssubtr_matr not available.\n");  }

#endif
 
#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA)

void smultA_vn(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   LINK *pl, *stop;
   DOUBLE sum[DIM];
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,Z))
      SET1(NDD(theNode,y),sum)
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               SET4(sum,NDD(pnode,x),COEFFL(pl,Z))
            }
            SET1(NDD(theNode,y),sum)
         }
   }
}

#else  /*  if !(N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) &&
               (N_DATA & VECTOR_NODE_DATA))  */

void smultA_vn(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vn not available.\n");  }

#endif  /*  !(N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) &&
             (N_DATA & VECTOR_NODE_DATA))  */

#if (N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA)

void multA_sn_X_vn(tGrid,Z,x,y,t)  /*  y := A x  */
GRID *tGrid;           /*  rows ... scalars in nodes                          */
INT Z, x, y, t;        /*  columns ... vectors in nodes                       */
{
   NODE *pnode, *pn;
   LINK *pli, *stop;

   if (!(t & STOP_IS_FIRST_INNER))
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ){
         NDS(pnode,y) = 0.;
         for (pli = pnode->start; pli != NULL; pli = pli->next){
            pn=NBNODE(pli);
            NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
         }
      }
   else
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ){
         NDS(pnode,y) = DOT(COEFFBP(pnode,Z),NDD(pnode,x));
         for (pli = pnode->tstart; pli != NULL; pli = pli->next)
            if (NOT_FN(pn=NBNODE(pli)))
               NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
      }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ){
      NDS(pnode,y) = DOT(COEFFBP(pnode,Z),NDD(pnode,x));
      stop = STOP_NN(pnode,t);
      for (pli = T_START(pnode,t); pli != stop; pli = pli->next){
         pn = NBNODE(pli);
         NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
      }
   }
}

void multA_vn_X_sn(tGrid,Z,x,y,t)  /*  y = A x  */
GRID *tGrid;           /*  rows ... vectors in nodes                          */
INT Z, x, y, t;        /*  columns ... scalars in nodes                       */
{
   NODE *pnode, *pn;
   LINK *pli;

   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      SET2(NDD(pnode,y),COEFFBP(pnode,Z),NDS(pnode,x))
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn=NBNODE(pli);
         SET4(NDD(pn,y),COEFFBLP(pli,Z),NDS(pnode,x))
      }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ){
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn = NBNODE(pli);
         SET4(NDD(pn,y),COEFFBLP(pli,Z),NDS(pnode,x))
      }
   }
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && 
                (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA)    */

void multA_sn_X_vn(tGrid,Z,x,y,t)  /*  y := A x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sn_X_vn not available.\n");  }

void multA_vn_X_sn(tGrid,Z,x,y,t)  /*  y = A x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_vn_X_sn not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA)

void vset_mat_value_vn(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   LINK *plink;
    
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      MSET7(COEFFNNP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         MSET7(COEFFLLP(plink,Z),r)
   }
}

void vmultA_vn(tGrid,Z,x,y,t) /*  y := Z x  */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   LINK *pl, *stop;
   DOUBLE sum[DIM];

   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         MSET4(sum,COEFFLLP(pl,Z),NDD(NBNODE(pl),x))
      SET1(NDD(theNode,y),sum)
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               MSET4(sum,COEFFLLP(pl,Z),NDD(pnode,x))
            }
            SET1(NDD(theNode,y),sum)
         }
   }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && 
                (N_DATA & VECTOR_NODE_DATA))  */

void vset_mat_value_vn(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: vset_mat_value_vn not available.\n");  }

void vmultA_vn(tGrid,Z,x,y,t) /*  y := Z x  */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: vmultA_vn not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)

void sset_mat_value_sn_sf(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
  
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ) {
      COEFFN(pnode,Z) = r;
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         COEFFL(plink,Z) = r;
      for (pnflink=pnode->tnfstart;pnflink!=NULL;pnflink=pnflink->next)
         COEFFL(pnflink,Z) = r;
   }
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ) {
      COEFF_FF(pface,Z) = r;
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next) 
         COEFF_FL(pflink,Z) = r;
      for (pfnlink=pface->tfnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         COEFF_FL(pfnlink,Z) = r;
   }
}

void smakeAT_sn_sf(tGrid,Z,a)
GRID *tGrid;
INT Z, a;
{
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *pli;
   NFLINK *pnfl;
   FLINK *pfl, *pfli;
   FNLINK *pfnl;
	
   for (theNode = FIRSTN(tGrid); theNode; theNode = SUCC(theNode)){
      COEFFN(theNode,a) = COEFFN(theNode,Z);
      for (pl = TSTART(theNode); pl; pl = NEXT(pl)){
         pnode = NBNODE(pl);
         for (pli = TSTART(pnode); pli->nbnode != theNode; pli = NEXT(pli));
         COEFFL(pl,a) = COEFFL(pli,Z);
      }
      for (pnfl = TNFSTART(theNode); pnfl; pnfl = NEXT(pnfl)){
         pface = NBFACE(pnfl);
         for (pfnl = TFNSTART(pface); pfnl->nbnode != theNode; pfnl=NEXT(pfnl));
         COEFFL(pnfl,a) = COEFFL(pfnl,Z);
         COEFFL(pfnl,a) = COEFFL(pnfl,Z);
      } 
   }
   for (theFace = FIRSTF(tGrid); theFace; theFace = SUCC(theFace)){
      COEFF_FF(theFace,a) = COEFF_FF(theFace,Z);
      for (pfl = TFSTART(theFace); pfl; pfl = NEXT(pfl)){
         pface = NBFACE(pfl);
         for (pfli = TFSTART(pface); pfli->nbface != theFace; pfli =NEXT(pfli));
         COEFF_FL(pfl,a) = COEFF_FL(pfli,Z);
      }
   }
}

#else

void sset_mat_value_sn_sf(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_sn_sf not available.\n");  }

void smakeAT_sn_sf(tGrid,Z,a)
GRID *tGrid;
INT Z, a;
{  eprintf("Error: smakeAT_sn_sf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void smultA_sn_sf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   LINK *pl, *stop;
   NFLINK *pnfl, *nfstop;
   FLINK *pfl, *fstop;
   FNLINK *pfnl, *fnstop;
   DOUBLE sum;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      sum = NDS(theNode,x)*COEFFN(theNode,Z);
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         sum += NDS(NBNODE(pl),x)*COEFFL(pl,Z);
      nfstop = STOP_NF(theNode,t);
      for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
         sum += COEFFL(pnfl,Z)*FD(NBFACE(pnfl),x);
      NDS(theNode,y) = sum;
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sum = COEFF_FF(theFace,Z)*FD(theFace,x);
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         sum += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      fnstop = STOP_FN(theFace,t);
      for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
         sum += COEFF_FL(pfnl,Z)*NDS(NBNODE(pfnl),x);
      FD(theFace,y) = sum;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            sum = NDS(theNode,x)*COEFFN(theNode,Z);
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
               sum += NDS(ltop_node(NBNODE(pl),tGrid),x)*COEFFL(pl,Z);
            nfstop = STOP_NF(theNode,t);
            for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
               sum += COEFFL(pnfl,Z)*FD(ltop_face(NBFACE(pnfl),tGrid),x);
            NDS(theNode,y) = sum;
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sum = COEFF_FF(theFace,Z)*FD(theFace,x);
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
               sum += COEFF_FL(pfl,Z)*FD(ltop_face(NBFACE(pfl),tGrid),x);
            fnstop = STOP_FN(theFace,t);
            for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
               sum += COEFF_FL(pfnl,Z)*NDS(ltop_node(NBNODE(pfnl),tGrid),x);
            FD(theFace,y) = sum;
         }
   }
}

void sSOR_step_forward_sn_sf(tGrid,Z,b,x,om)    /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   NODE *theNode;
   FACE *theFace;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   DOUBLE sum;

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      if (fabs(COEFFN(theNode,Z)) > EPSA){
         sum = NDS(theNode,b);
         for (pl=START(theNode); pl; pl=NEXT(pl))
            sum -= COEFFL(pl,Z)*NDS(NBNODE(pl),x);
         for (pnfl=NFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
            sum -= COEFFL(pnfl,Z)*FD(NBFACE(pnfl),x);
         NDS(theNode,x) += om*(sum/COEFFN(theNode,Z)-NDS(theNode,x));
      }
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      if (fabs(COEFF_FF(theFace,Z)) > EPSA){
         sum = FD(theFace,b);
         for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
            sum -= COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
         for (pfnl=FNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
            sum -= COEFF_FL(pfnl,Z)*NDS(NBNODE(pfnl),x);
         FD(theFace,x) += om*(sum/COEFF_FF(theFace,Z)-FD(theFace,x));
      }
}

#else

void smultA_sn_sf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_sn_sf not available.\n");  }

void sSOR_step_forward_sn_sf(tGrid,Z,b,x,om)
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: sSOR_step_forward_sn_sf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void sset_mat_value_vn_sf(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
  
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ) {
      COEFFN(pnode,Z) = r;
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         COEFFL(plink,Z) = r;
      SET_NF
   }
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ) {
      COEFF_FF(pface,Z) = r;
      SET_FF
      SET_FN
   }
}

void smultA_vn_sf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *stop;
   NFLINK *pnfl, *nfstop;
   FLINK *pfl, *fstop;
   FNLINK *pfnl, *fnstop;
   DOUBLE sum[DIM], sumf;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,Z))
      nfstop = STOP_NF(theNode,t);
      for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
         SET4(sum,COEFF_NFP(pnfl,Z),FD(NBFACE(pnfl),x))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      fnstop = STOP_FN(theFace,t);
      for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
         sumf += DOT(COEFF_FNP(pfnl,Z),NDD(NBNODE(pfnl),x));
      FD(theFace,y) = sumf;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               SET4(sum,NDD(pnode,x),COEFFL(pl,Z))
            }
            nfstop = STOP_NF(theNode,t);
            for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl)){
               pface = ltop_face(NBFACE(pnfl),tGrid);
               SET4(sum,COEFF_NFP(pnfl,Z),FD(pface,x));
            } 
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
               sumf += COEFF_FL(pfl,Z)*FD(ltop_face(NBFACE(pfl),tGrid),x);
            fnstop = STOP_FN(theFace,t);
            for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl)){
               pnode = ltop_node(NBNODE(pfnl),tGrid);
               sumf += DOT(COEFF_FNP(pfnl,Z),NDD(pnode,x));
            } 
            FD(theFace,y) = sumf;
         }
   }
}

void smultAT_vn_sf(tGrid,Z,x,y) /*  y := Z^T x  */
GRID *tGrid;
INT Z,x,y;
{
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *pli;
   NFLINK *pnfl;
   FLINK *pfl, *pfli;
   FNLINK *pfnl;
   DOUBLE sum[DIM], sumf;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      for (pl=START(theNode); pl!=NULL; pl=NEXT(pl)){
         pnode = NBNODE(pl);
         for (pli = START(pnode); pli->nbnode != theNode; pli = pli->next);
         SET4(sum,NDD(pnode,x),COEFFL(pli,Z))
      }
      for (pnfl=NFSTART(theNode); pnfl != NULL; pnfl=NEXT(pnfl)){
         pface = NBFACE(pnfl);
         for (pfnl=FNSTART(pface); pfnl->nbnode != theNode; pfnl=NEXT(pfnl));
         SET4(sum,COEFF_FNP(pfnl,Z),FD(pface,x))
      } 
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
   
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
         pface = NBFACE(pfl);
         for (pfli=FSTART(pface); pfli->nbface != theFace; pfli=NEXT(pfli));
         sumf += COEFF_FL(pfli,Z)*FD(pface,x);
      }
      for (pfnl=FNSTART(theFace); pfnl != NULL; pfnl=NEXT(pfnl)){
         pnode = NBNODE(pfnl);
         for (pnfl=NFSTART(pnode); pnfl->nbface != theFace; pnfl=NEXT(pnfl));
         sumf += DOT(COEFF_NFP(pnfl,Z),NDD(pnode,x));
      }
      FD(theFace,y) = sumf;
   }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
               (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && 
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void sset_mat_value_vn_sf(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_vn_sf not available.\n");  }

void smultA_vn_sf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vn_sf not available.\n");  }

void smultAT_vn_sf(tGrid,Z,x,y) /*  y := Z^T x  */
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: smultAT_vn_sf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void smultA_vn_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *stop;
   NFLINK *pnfl, *nfstop;
   FLINK *pfl, *fstop;
   FNLINK *pfnl, *fnstop;
   DOUBLE sum[DIM];
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,Z))
      nfstop = STOP_NF(theNode,t);
      for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
         SET4(sum,FDVP(NBFACE(pnfl),x),COEFFL(pnfl,Z))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      SET2(sum,FDVP(theFace,x),COEFF_FF(theFace,Z))
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         SET4(sum,FDVP(NBFACE(pfl),x),COEFF_FL(pfl,Z))
      fnstop = STOP_FN(theFace,t);
      for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
         SET4(sum,NDD(NBNODE(pfnl),x),COEFF_FL(pfnl,Z))
      SET1(FDVP(theFace,y),sum)
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               SET4(sum,NDD(pnode,x),COEFFL(pl,Z))
            }
            nfstop = STOP_NF(theNode,t);
            for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl)){
               pface = ltop_face(NBFACE(pnfl),tGrid);
               SET4(sum,FDVP(pface,x),COEFFL(pnfl,Z))
            } 
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            SET2(sum,FDVP(theFace,x),COEFF_FF(theFace,Z))
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl)){
               pface = ltop_face(NBFACE(pfl),tGrid);
               SET4(sum,FDVP(pface,x),COEFF_FL(pfl,Z))
            }
            fnstop = STOP_FN(theFace,t);
            for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl)){
               pnode = ltop_node(NBNODE(pfnl),tGrid);
               SET4(sum,NDD(pnode,x),COEFF_FL(pfnl,Z))
            } 
            SET1(FDVP(theFace,y),sum)
         }
   }
}

#else

void smultA_vn_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vn_vf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void smult_inv_diag_vn_vf(tGrid,Z,x,y)  /*  y := diag(Z)^{-1} x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode))
      SET22(NDD(theNode,y),NDD(theNode,x),COEFFN(theNode,Z))
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      SET22(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid))
            SET22(NDD(theNode,y),NDD(theNode,x),COEFFN(theNode,Z))
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            SET22(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))
   }
}

#else

void smult_inv_diag_vn_vf(tGrid,Z,x,y)  /*  y := diag(Z)^{-1} x */
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: smult_inv_diag_vn_vf not available.\n");  }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExE_MATR)

void sset_mat_value_sn_se(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pel;
  
   sset_mat_value(tGrid,Z,r);
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      SSET7(COEFF_NEP(pel,Z),r)
      SSET7(COEFF_ENP(pel,Z),r)
      COEFF_EE(pel,Z) = r;
   }
}
   
#else

void sset_mat_value_sn_se(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_sn_se not available.\n");  }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExE_MATR) && (N_DATA & SCALAR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (DIM==2)

void smultA_sn_se(tGrid,Z,x,y,q,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,q,t;
{
   ELEMENT *pel;
   NODE *theNode;
  
   if (!(t & STOP_IS_FIRST_INNER)){
      b_copy(tGrid,x,q,1,Q_SN);
      for (theNode = FDBN(tGrid); theNode != FIRSTNODE(tGrid); theNode = SUCC(theNode))
         NDS(theNode,x) = 0.;
   }
   smultA(tGrid,Z,x,y,t);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NDS(pel->n[0],y) += COEFF_NE(pel,Z,0)*ED(pel,x);
      NDS(pel->n[1],y) += COEFF_NE(pel,Z,1)*ED(pel,x);
      NDS(pel->n[2],y) += COEFF_NE(pel,Z,2)*ED(pel,x);
      ED(pel,y) = COEFF_EN(pel,Z,0)*NDS(pel->n[0],x) + COEFF_EN(pel,Z,1)*NDS(pel->n[1],x) +
                  COEFF_EN(pel,Z,2)*NDS(pel->n[2],x) + COEFF_EE(pel,Z)*ED(pel,x);
   }
   if (!(t & STOP_IS_FIRST_INNER))
      b_copy(tGrid,q,x,1,Q_SN);
}
   
#else

void smultA_sn_se(tGrid,Z,x,y,q,t)
GRID *tGrid; INT Z,x,y,q,t;
{  eprintf("Error: smultA_sn_se not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (N_DATA & SCALAR_NODE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

#if DIM == 3

void multA_sn_X_ve(tGrid,Z,xe,y,t)  /*  y := y + A xe  */
GRID *tGrid;           /*  rows ... scalars in nodes                          */
INT Z, xe, y, t;       /*  columns ... vectors in elements                    */
{
   ELEMENT *pel;

   if (!(t & STOP_IS_FIRST_INNER))
      for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
         NDS(pel->n[0],y) += DOT(COEFF_BNP(pel,Z,0),EDVP(pel,xe));
         NDS(pel->n[1],y) += DOT(COEFF_BNP(pel,Z,1),EDVP(pel,xe));
         NDS(pel->n[2],y) += DOT(COEFF_BNP(pel,Z,2),EDVP(pel,xe));
         NDS(pel->n[3],y) += DOT(COEFF_BNP(pel,Z,3),EDVP(pel,xe));
      }
}

void multA_ve_X_sn(tGrid,Z,x,ye,t)  /*  ye = A x  */
GRID *tGrid;           /*  rows ... vectors in elements                       */
INT Z, x, ye, t;       /*  columns ... scalars in nodes                       */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      EDV(pel,ye,0) = COEFF_BN(pel,Z,0,0)*NDS(pel->n[0],x)+
                      COEFF_BN(pel,Z,1,0)*NDS(pel->n[1],x)+
                      COEFF_BN(pel,Z,2,0)*NDS(pel->n[2],x)+
                      COEFF_BN(pel,Z,3,0)*NDS(pel->n[3],x);
      EDV(pel,ye,1) = COEFF_BN(pel,Z,0,1)*NDS(pel->n[0],x)+
                      COEFF_BN(pel,Z,1,1)*NDS(pel->n[1],x)+
                      COEFF_BN(pel,Z,2,1)*NDS(pel->n[2],x)+
                      COEFF_BN(pel,Z,3,1)*NDS(pel->n[3],x);
      EDV(pel,ye,2) = COEFF_BN(pel,Z,0,2)*NDS(pel->n[0],x)+
                      COEFF_BN(pel,Z,1,2)*NDS(pel->n[1],x)+
                      COEFF_BN(pel,Z,2,2)*NDS(pel->n[2],x)+
                      COEFF_BN(pel,Z,3,2)*NDS(pel->n[3],x);
   }
}

#else  /*  if DIM != 3  */

void multA_sn_X_ve(tGrid,Z,xe,y,t)  /*  y := y + A xe  */
GRID *tGrid;           /*  rows ... scalars in nodes                          */
INT Z, xe, y, t;       /*  columns ... vectors in elements                    */
{
   ELEMENT *pel;

   if (!(t & STOP_IS_FIRST_INNER))
      for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
         NDS(pel->n[0],y) += DOT(COEFF_BNP(pel,Z,0),EDVP(pel,xe));
         NDS(pel->n[1],y) += DOT(COEFF_BNP(pel,Z,1),EDVP(pel,xe));
         NDS(pel->n[2],y) += DOT(COEFF_BNP(pel,Z,2),EDVP(pel,xe));
      }
}

void multA_ve_X_sn(tGrid,Z,x,ye,t)  /*  ye = A x  */
GRID *tGrid;           /*  rows ... vectors in elements                       */
INT Z, x, ye, t;       /*  columns ... scalars in nodes                       */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      EDV(pel,ye,0) = COEFF_BN(pel,Z,0,0)*NDS(pel->n[0],x)+
                      COEFF_BN(pel,Z,1,0)*NDS(pel->n[1],x)+
                      COEFF_BN(pel,Z,2,0)*NDS(pel->n[2],x);
      EDV(pel,ye,1) = COEFF_BN(pel,Z,0,1)*NDS(pel->n[0],x)+
                      COEFF_BN(pel,Z,1,1)*NDS(pel->n[1],x)+
                      COEFF_BN(pel,Z,2,1)*NDS(pel->n[2],x);
   }
}

#endif

#else  /*  if !((E_DATA & ExDN_MATR) && (N_DATA & SCALAR_NODE_DATA) && 
                (E_DATA & VECTOR_ELEMENT_DATA))  */

void multA_sn_X_ve(tGrid,Z,xe,y,t)  /*  y := y + A xe  */
GRID *tGrid; INT Z, xe, y, t;
{  eprintf("Error: multA_sn_X_ve not available.\n");  }

void multA_ve_X_sn(tGrid,Z,x,ye,t)  /*  ye = A x  */
GRID *tGrid; INT Z, x, ye, t;
{  eprintf("Error: multA_ve_X_sn not available.\n");  }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR)

void sset_mat_value_sn_sf_se(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pel;

   sset_mat_value_sn_sf(theGrid,Z,r);
   for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
      SSET7(COEFF_NEP(pel,Z),r)
      SSET7(COEFF_ENP(pel,Z),r)
      SSET7(COEFF_FEP(pel,Z),r)
      SSET7(COEFF_EFP(pel,Z),r)
      COEFF_EE(pel,Z) = r;
   }
}

#else

void sset_mat_value_sn_sf_se(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_sn_sf_se not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (E_DATA & ExDN_MATR)

void set_mat_value_sn_X_vn_vf_ve(tGrid,Z,r)  /*  A_ij := r  */
GRID *tGrid;             /*  rows ... scalars in nodes                        */
INT Z;                   /*  columns ... vectors in nodes, faces and elements */
FLOAT r;
{
   ELEMENT *pelem;
   NODE *pnode;
   LINK *pli;
   NFLINK *pnf; 
   INT i;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ){
      SET7(COEFFBP(pnode,Z),r)
      for (pli = TSTART(pnode); pli != NULL; pli = pli->next)
         SET7(COEFFBLP(pli,Z),r)
      for (pnf = TNFSTART(pnode); pnf != NULL; pnf = pnf->next)
         SET7(COEFF_NFP(pnf,Z),r)
   }
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ)
      for (i=0; i < NVERT; i++)
         SET7(COEFF_BNP(pelem,Z,i),r)
}

#else

void set_mat_value_sn_X_vn_vf_ve(tGrid,Z,r)  /*  A_ij := r  */
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_sn_X_vn_vf_ve not available.\n");  }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void smultA_vn_vf_ve(tGrid,Z,x,y,q,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,q,t;
{
   ELEMENT *pel;
   NODE *theNode;
   FACE *theFace;
  
   if (!(t & STOP_IS_FIRST_INNER)){
      b_copy(tGrid,x,q,1,Q_VNVF);
      for (theNode=FDBN(tGrid); theNode!=FIRSTNODE(tGrid);theNode=SUCC(theNode))
         SET_VALUE(theNode,0.,x)
      for (theFace=FDBF(tGrid); theFace!=FIRSTFACE(tGrid);theFace=SUCC(theFace))
         SET_VALUE(theFace,0.,x)
   }
   smultA_vn_vf(tGrid,Z,x,y,t);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      ADD_VMULT(pel->n,y,COEFF_NEP(pel,Z),EDVP(pel,x))
      ADD_VMULT(pel->f,y,COEFF_FEP(pel,Z),EDVP(pel,x))
      SET2(EDVP(pel,y),EDVP(pel,x),COEFF_EE(pel,Z))
      ADD_SUM(pel->n,x,EDVP(pel,y),COEFF_ENP(pel,Z))
      ADD_SUM(pel->f,x,EDVP(pel,y),COEFF_EFP(pel,Z))
   }
   if (!(t & STOP_IS_FIRST_INNER))
      b_copy(tGrid,q,x,1,Q_VNVF);
}
   
#else

void smultA_vn_vf_ve(tGrid,Z,x,y,q,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,q,t;
{  eprintf("Error: smultA_vn_vf_ve not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vset_mat_value_vn_vf(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
    
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      MSET7(COEFFNNP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         MSET7(COEFFLLP(plink,Z),r)
      for (pnflink=pnode->tnfstart;pnflink!=NULL;pnflink=pnflink->next)
         MSET7(COEFFLLP(pnflink,Z),r)
   }
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ){
      MSET7(COEFFNNP(pface,Z),r)
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)
         MSET7(COEFFLLP(pflink,Z),r)
      for (pfnlink=pface->tfnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         MSET7(COEFFLLP(pfnlink,Z),r)
   }
}

void vmultA_vn_vf(tGrid,Z,x,y,t) /*  y := Z x  */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *stop;
   NFLINK *pnfl, *nfstop;
   FLINK *pfl, *fstop;
   FNLINK *pfnl, *fnstop;
   DOUBLE sum[DIM];

   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         MSET4(sum,COEFFLLP(pl,Z),NDD(NBNODE(pl),x))
      nfstop = STOP_NF(theNode,t);
      for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
         MSET4(sum,COEFFLLP(pnfl,Z),FDVP(NBFACE(pnfl),x))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      MSET2(sum,COEFFNNP(theFace,Z),FDVP(theFace,x))
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         MSET4(sum,COEFFLLP(pfl,Z),FDVP(NBFACE(pfl),x))
      fnstop = STOP_FN(theFace,t);
      for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
         MSET4(sum,COEFFLLP(pfnl,Z),NDD(NBNODE(pfnl),x))
      SET1(FDVP(theFace,y),sum)
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               MSET4(sum,COEFFLLP(pl,Z),NDD(pnode,x))
            }
            nfstop = STOP_NF(theNode,t);
            for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl)){
               pface = ltop_face(NBFACE(pnfl),tGrid);
               MSET4(sum,COEFFLLP(pnfl,Z),FDVP(pface,x))
            }
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            MSET2(sum,COEFFNNP(theFace,Z),FDVP(theFace,x))
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl)){
               pface = ltop_face(NBFACE(pfl),tGrid);
               MSET4(sum,COEFFLLP(pfl,Z),FDVP(pface,x))
            }
            fnstop = STOP_FN(theFace,t);
            for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl)){
               pnode = ltop_node(NBNODE(pfnl),tGrid);
               MSET4(sum,COEFFLLP(pfnl,Z),NDD(pnode,x))
            }
            SET1(FDVP(theFace,y),sum)
         }
   }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && 
               (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) && 
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

void vset_mat_value_vn_vf(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: vset_mat_value_vn_vf not available.\n");  }

void vmultA_vn_vf(tGrid,Z,x,y,t) /*  y := Z x  */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: vmultA_vn_vf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & F_LINK_TO_FACES)

void sset_mat_value_sn_sf_nred(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
  
   for (pnode = FIRSTN(theGrid); pnode; pnode = pnode->succ){
      COEFFN(pnode,Z) = r;
      for (plink = pnode->tstart; plink; plink = plink->next)
         COEFFL(plink,Z) = r;
   }
  
   for (pface = FIRSTF(theGrid); pface; pface = pface->succ) {
      COEFF_FF(pface,Z) = r;
      for (pflink = pface->tfstart; pflink; pflink = pflink->next) 
         COEFF_FL(pflink,Z) = r;
   }
}

#else

void sset_mat_value_sn_sf_nred(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_sn_sf_nred not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES)

void sset_mat_value_sn_sf_nfred(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
  
   for (pnode = FIRSTN(theGrid); pnode; pnode = pnode->succ){
      COEFFN(pnode,Z) = r;
      for (plink = pnode->tstart; plink; plink = plink->next)
         COEFFL(plink,Z) = r;
   }
  
   for (pface = FIRSTF(theGrid); pface; pface = pface->succ) 
      COEFF_FF(pface,Z) = r;
}

#else

void sset_mat_value_sn_sf_nfred(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_sn_sf_nfred not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DATA_STR & LG_DATA)

void sset_mat_value_vn_sf_lg(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
   NLGLINK *pnlg;
   FLGLINK *pflg;
   LGNLINK *plgn;
   LGFLINK *plgf;
   LGLGLINK *plglg;
  
   for (pnode=FIRSTN(theGrid); pnode!=NULL; pnode=pnode->succ) {
      COEFFN(pnode,Z) = r;
      for (plink=pnode->tstart; plink!=NULL; plink=plink->next)
         COEFFL(plink,Z) = r;
      for (pnflink=pnode->tnfstart; pnflink!=NULL; pnflink=pnflink->next)
         SET7(COEFF_NFP(pnflink,Z),r)
      for (pnlg = pnode->lgstart; pnlg != NULL; pnlg = pnlg->next)
         MSET7NLG(COEFF_NLGP(pnlg,Z),r)
   }
   for (pface=FIRSTF(theGrid); pface!=NULL; pface=pface->succ) {
      COEFF_FF(pface,Z) = r;
      for (pflink=pface->tfstart; pflink!=NULL; pflink=pflink->next)
         COEFF_FL(pflink,Z) = r;
      for (pfnlink=pface->tfnstart; pfnlink!=NULL; pfnlink=pfnlink->next)
         SET7(COEFF_FNP(pfnlink,Z),r)
      for (pflg = pface->lgstart; pflg != NULL; pflg = pflg->next)
         SET7LG(COEFF_FLGP(pflg,Z),r)
   }
   for (pnode=FDBN(theGrid); pnode!=FIRSTNODE(theGrid); pnode=pnode->succ)
      if (pnode->lgd){
         MSET7LG(COEFF_LGP(pnode,Z),r)
         for (plglg = pnode->lgd->lgstart; plglg != NULL; plglg = plglg->next)
            MSET7LG(COEFF_LGLP(plglg,Z),r)
         for (plgn = pnode->lgd->nstart; plgn != NULL; plgn = plgn->next)
            MSET7LGN(COEFF_LGNP(plgn,Z),r)
         for (plgf = pnode->lgd->fstart; plgf != NULL; plgf = plgf->next)
            SET7LG(COEFF_LGFP(plgf,Z),r)
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
               (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) &&
               (DATA_STR & LG_DATA))  */

void sset_mat_value_vn_sf_lg(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: sset_mat_value_vn_sf_lg not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (DIM == 3)

FLOAT norm_of_A_vn_sf(tGrid,Z)
GRID *tGrid;
INT Z;
{
   GRID *theGrid;
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
   FLOAT sum, sum0, sum1, sum2, max=0.0;
   
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser){
      for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ)
         if (IS_LTOP_NODE(pnode,tGrid)){
            sum = fabs(COEFFN(pnode,Z));
            for (plink = pnode->start; plink != NULL; plink = plink->next)
               sum += fabs(COEFFL(plink,Z));
            sum0 = sum1 = sum2 = 0.0;
            for (pnflink = pnode->nfstart; pnflink != NULL; 
                                                      pnflink = pnflink->next){
               sum0 += fabs(COEFF_NF(pnflink,Z,0));
               sum1 += fabs(COEFF_NF(pnflink,Z,1));
               sum2 += fabs(COEFF_NF(pnflink,Z,2));
            }
            if (sum0 > sum1) sum1 = sum0;
            if (sum1 > sum2) sum2 = sum1;
            sum += sum2;
            if (sum > max) max = sum;
         }
         
      for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
         if (IS_LTOP_FACE(pface,tGrid)){
            sum = fabs(COEFF_FF(pface,Z));
            for (pflink = pface->fstart; pflink != NULL; pflink = pflink->next)
               sum += fabs(COEFF_FL(pflink,Z));
            for (pfnlink = pface->fnstart; pfnlink != NULL; 
                                                       pfnlink = pfnlink->next)
               sum += fabs(COEFF_FN(pfnlink,Z,0)) + fabs(COEFF_FN(pfnlink,Z,1)) 
                                                  + fabs(COEFF_FN(pfnlink,Z,2));
            if (sum > max) max = sum;
         }
   }
   return(max);
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
     (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
     (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
     (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (DIM == 3))  */

FLOAT norm_of_A_vn_sf(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: norm_of_A_vn_sf not available.\n"); return(0.);  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)

void smakeAT_vn_sf(tGrid,Z,a)
GRID *tGrid;
INT Z, a;
{
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl, *pli;
   NFLINK *pnfl;
   FLINK *pfl, *pfli;
   FNLINK *pfnl;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      COEFFN(theNode,a) = COEFFN(theNode,Z);
      	   
      for (pl=START(theNode); pl!=NULL; pl=NEXT(pl)){
         pnode = NBNODE(pl);
         for (pli = START(pnode); pli->nbnode != theNode; pli = pli->next);
         COEFFL(pl,a) = COEFFL(pli,Z);
      }
      for (pnfl=NFSTART(theNode); pnfl != NULL; pnfl=NEXT(pnfl)){
         pface = NBFACE(pnfl);
         for (pfnl=FNSTART(pface); pfnl->nbnode != theNode; pfnl=NEXT(pfnl));
         SET1(COEFF_NFP(pnfl,a),COEFF_FNP(pfnl,Z))
         SET1(COEFF_FNP(pfnl,a),COEFF_NFP(pnfl,Z))
      } 
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      COEFF_FF(theFace,a) = COEFF_FF(theFace,Z);
   
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
         pface = NBFACE(pfl);
         for (pfli=FSTART(pface); pfli->nbface != theFace; pfli=NEXT(pfli));
         COEFF_FL(pfl,a) = COEFF_FL(pfli,Z);
      }
   }
}

void sadd_matr_vn_sf(mg,Z,a)  /*  Z := Z + A  */
MULTIGRID *mg;
INT Z, a;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
   
   for (pnode=FIRSTNODE(TOP_GRID(mg));pnode!=NULL;pnode=pnode->succ) {
      COEFFN(pnode,Z) += COEFFN(pnode,a);
      for (plink=pnode->start;plink!=NULL;plink=plink->next)
         COEFFL(plink,Z) += COEFFL(plink,a);
      for (pnflink=pnode->nfstart;pnflink!=NULL;pnflink=pnflink->next)
         SET5(COEFF_NFP(pnflink,Z),COEFF_NFP(pnflink,a))
   }

   for (pface=FIRSTFACE(TOP_GRID(mg));pface!=NULL;pface=pface->succ) {
      COEFF_FF(pface,Z) += COEFF_FF(pface,a);
      for (pflink=pface->fstart;pflink!=NULL;pflink=pflink->next)
         COEFF_FL(pflink,Z) += COEFF_FL(pflink,a);
      for (pfnlink=pface->fnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         SET5(COEFF_FNP(pfnlink,Z),COEFF_FNP(pfnlink,a))
   }
} 
 
void ssubtr_matr_vn_sf(mg,Z,a)  /*  Z := Z - A  */
MULTIGRID *mg;
INT Z, a;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;

   for (pnode=FIRSTNODE(TOP_GRID(mg));pnode!=NULL;pnode=pnode->succ) {
      COEFFN(pnode,Z) -= COEFFN(pnode,a);
      for (plink=pnode->start;plink!=NULL;plink=plink->next)
         COEFFL(plink,Z) -= COEFFL(plink,a);
      for (pnflink=pnode->nfstart;pnflink!=NULL;pnflink=pnflink->next)
         SET6(COEFF_NFP(pnflink,Z),COEFF_NFP(pnflink,a))
   }
  
   for (pface=FIRSTFACE(TOP_GRID(mg));pface!=NULL;pface=pface->succ) {
      COEFF_FF(pface,Z) -= COEFF_FF(pface,a);
      for (pflink=pface->fstart;pflink!=NULL;pflink=pflink->next)
         COEFF_FL(pflink,Z) -= COEFF_FL(pflink,a);
      for (pfnlink=pface->fnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         SET6(COEFF_FNP(pfnlink,Z),COEFF_FNP(pfnlink,a))
   }
} 
 
#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
                (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)  */ 

void smakeAT_vn_sf(tGrid,Z,a)
GRID *tGrid; INT Z, a;
{  eprintf("Error: smakeAT_vn_sf not available.\n");  }

void sadd_matr_vn_sf(mg,Z,a)  /*  Z := Z + A  */
MULTIGRID *mg; INT Z, a;
{  eprintf("Error: sadd_matr_vn_sf not available.\n");  }

void ssubtr_matr_vn_sf(mg,Z,a)  /*  Z := Z - A  */
MULTIGRID *mg; INT Z, a;
{  eprintf("Error: ssubtr_matr_vn_sf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & F_LINK_TO_FACES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void smultA_vn_sf_nred(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace;
   LINK *pl, *stop;
   FLINK *pfl, *fstop;
   DOUBLE sum[DIM], sumf;
	
   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,Z))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      FD(theFace,y) = sumf;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               SET4(sum,NDD(pnode,x),COEFFL(pl,Z))
            }
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
               sumf += COEFF_FL(pfl,Z)*FD(ltop_face(NBFACE(pfl),tGrid),x);
            FD(theFace,y) = sumf;
         }
   }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & F_LINK_TO_FACES) && 
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void smultA_vn_sf_nred(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vn_sf_nred not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void smultA_vn_sf_nfred(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace;
   LINK *pl, *stop;
   DOUBLE sum[DIM];

   for (theNode=FIRSTNODE(tGrid); theNode!= NULL; theNode=SUCC(theNode)){
      SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,Z))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,y) = COEFF_FF(theFace,Z)*FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            SET2(sum,NDD(theNode,x),COEFFN(theNode,Z))
            stop = STOP_NN(theNode,t);
            for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               SET4(sum,NDD(pnode,x),COEFFL(pl,Z))
            }
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,y) = COEFF_FF(theFace,Z)*FD(theFace,x);
   }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && 
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void smultA_vn_sf_nfred(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vn_sf_nfred not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void vset_mat_value_vn_sf(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
    
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      MSET7(COEFFNNP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         MSET7(COEFFLLP(plink,Z),r);
      for (pnflink=pnode->tnfstart;pnflink!=NULL;pnflink=pnflink->next)
         SET7(COEFF_NFP(pnflink,Z),r)
   }
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ){
      COEFF_FF(pface,Z)=r;
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)
         COEFF_FL(pflink,Z)=r;
      for (pfnlink=pface->tfnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         SET7(COEFF_FNP(pfnlink,Z),r)
   }
}

void vmultA_vn_sf(tGrid,Z,x,y) /*  y := Z x  */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   FACE *theFace, *pface;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   DOUBLE sum[DIM], sumf;
   	
   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
      for (pl=START(theNode); pl!=NULL; pl=NEXT(pl))
         MSET4(sum,COEFFLLP(pl,Z),NDD(NBNODE(pl),x))
      for (pnfl=NFSTART(theNode); pnfl != NULL; pnfl=NEXT(pnfl))
         SET4(sum,COEFF_NFP(pnfl,Z),FD(NBFACE(pnfl),x))
      SET1(NDD(theNode,y),sum)
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl))
         sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      for (pfnl=FNSTART(theFace); pfnl != NULL; pfnl=NEXT(pfnl))
         sumf += DOT(COEFF_FNP(pfnl,Z),NDD(NBNODE(pfnl),x));
      FD(theFace,y) = sumf;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid)){
            MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
            for (pl=START(theNode); pl!=NULL; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               MSET4(sum,COEFFLLP(pl,Z),NDD(pnode,x))
            }
            for (pnfl=NFSTART(theNode); pnfl != NULL; pnfl=NEXT(pnfl)){
               pface = ltop_face(NBFACE(pnfl),tGrid);
               SET4(sum,COEFF_NFP(pnfl,Z),FD(pface,x))
            }
            SET1(NDD(theNode,y),sum)
         }
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
            for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl))
               sumf += COEFF_FL(pfl,Z)*FD(ltop_face(NBFACE(pfl),tGrid),x);
            for (pfnl=FNSTART(theFace); pfnl != NULL; pfnl=NEXT(pfnl)){
               pnode = ltop_node(NBNODE(pfnl),tGrid);
               sumf += DOT(COEFF_FNP(pfnl,Z),NDD(pnode,x));
            } 
            FD(theFace,y) = sumf;
         }
   }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
               (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void vset_mat_value_vn_sf(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: vset_mat_value_vn_sf not available.\n");  }

void vmultA_vn_sf(tGrid,Z,x,y) /*  y := Z x  */
GRID *tGrid;
INT Z,x,y;
{  eprintf("Error: vmultA_vn_sf not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DIM == 3)

FLOAT vnorm_of_A_vn_sf(tGrid,Z)
GRID *tGrid;
INT Z;
{
   GRID *theGrid;
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
   FLOAT sum, sum0, sum1, sum2, max=0.0;
   
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser){
      for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ)
         if (IS_LTOP_NODE(pnode,tGrid)){
            sum0 = fabs(COEFFNN(pnode,Z,0,0)) + fabs(COEFFNN(pnode,Z,0,1))
                                              + fabs(COEFFNN(pnode,Z,0,2));
            sum1 = fabs(COEFFNN(pnode,Z,1,0)) + fabs(COEFFNN(pnode,Z,1,1))
                                              + fabs(COEFFNN(pnode,Z,1,2));
            sum2 = fabs(COEFFNN(pnode,Z,2,0)) + fabs(COEFFNN(pnode,Z,2,1))
                                              + fabs(COEFFNN(pnode,Z,2,2));
            for (plink = pnode->start; plink != NULL; plink = plink->next){
               sum0 += fabs(COEFFLL(plink,Z,0,0)) + fabs(COEFFLL(plink,Z,0,1))
                                                  + fabs(COEFFLL(plink,Z,0,2));
               sum1 += fabs(COEFFLL(plink,Z,1,0)) + fabs(COEFFLL(plink,Z,1,1))
                                                  + fabs(COEFFLL(plink,Z,1,2));
               sum2 += fabs(COEFFLL(plink,Z,2,0)) + fabs(COEFFLL(plink,Z,2,1))
                                                  + fabs(COEFFLL(plink,Z,2,2));
            }
            for (pnflink = pnode->nfstart; pnflink != NULL; 
                                                      pnflink = pnflink->next){
               sum0 += fabs(COEFF_NF(pnflink,Z,0));
               sum1 += fabs(COEFF_NF(pnflink,Z,1));
               sum2 += fabs(COEFF_NF(pnflink,Z,2));
            }
            if (sum0 > max) max = sum0;
            if (sum1 > max) max = sum1;
            if (sum2 > max) max = sum2;
         }
         
      for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
         if (IS_LTOP_FACE(pface,tGrid)){
            sum = fabs(COEFF_FF(pface,Z));
            for (pflink = pface->fstart; pflink != NULL; pflink = pflink->next)
               sum += fabs(COEFF_FL(pflink,Z));
            for (pfnlink = pface->fnstart; pfnlink != NULL; 
                                                       pfnlink = pfnlink->next)
               sum += fabs(COEFF_FN(pfnlink,Z,0)) + fabs(COEFF_FN(pfnlink,Z,1)) 
                                                  + fabs(COEFF_FN(pfnlink,Z,2));
            if (sum > max) max = sum;
         }
   }
   return(max);
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
               (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
               (DIM == 3))   */

FLOAT vnorm_of_A_vn_sf(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: vnorm_of_A_vn_sf not available.\n"); return(0.);  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DATA_STR & LG_DATA)

void vset_mat_value_vn_sf_lg(theGrid,Z,r)
GRID *theGrid; 
INT Z;
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   NFLINK *pnflink;
   FNLINK *pfnlink;
   NLGLINK *pnlg;
   FLGLINK *pflg;
   LGNLINK *plgn;
   LGFLINK *plgf;
   LGLGLINK *plglg;
   
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      MSET7(COEFFNNP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         MSET7(COEFFLLP(plink,Z),r)
      for (pnflink=pnode->tnfstart;pnflink!=NULL;pnflink=pnflink->next)
         SET7(COEFF_NFP(pnflink,Z),r)
      for (pnlg = pnode->lgstart; pnlg != NULL; pnlg = pnlg->next)
         MSET7NLG(COEFF_NLGP(pnlg,Z),r)
   }
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ){
      COEFF_FF(pface,Z)=r;
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)
         COEFF_FL(pflink,Z)=r;
      for (pfnlink=pface->tfnstart;pfnlink!=NULL;pfnlink=pfnlink->next)
         SET7(COEFF_FNP(pfnlink,Z),r)
      for (pflg = pface->lgstart; pflg != NULL; pflg = pflg->next)
         SET7LG(COEFF_FLGP(pflg,Z),r)
   }
   for (pnode=FDBN(theGrid); pnode != FIRSTNODE(theGrid); pnode=pnode->succ)
      if (pnode->lgd){
         MSET7LG(COEFF_LGP(pnode,Z),r)
         for (plglg = pnode->lgd->lgstart; plglg != NULL; plglg = plglg->next)
            MSET7LG(COEFF_LGLP(plglg,Z),r)
         for (plgn = pnode->lgd->nstart; plgn != NULL; plgn = plgn->next)
            MSET7LGN(COEFF_LGNP(plgn,Z),r)
         for (plgf = pnode->lgd->fstart; plgf != NULL; plgf = plgf->next)
            SET7LG(COEFF_LGFP(plgf,Z),r)
      }
}

INT vmultA_vn_sf_lg(tGrid,Z,x,y)
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   NODE *theNode, *pn;
   FACE *theFace, *pface;
   LINK *pl;
   NLGLINK *pnlg;
   NFLINK *pnfl;
   LGNLINK *plgn;
   LGFLINK *plgf;
   LGLGLINK *plg;
   FLINK *pfl;
   FNLINK *pfnl;
   FLGLINK *pflg;
   DOUBLE sum[DIM], sumf;
   	
   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      MSET2(sum,COEFFNNP(theNode,Z),NDD(theNode,x))
      for (pl=START(theNode); pl!=NULL; pl=NEXT(pl))
         MSET4(sum,COEFFLLP(pl,Z),NDD(NBNODE(pl),x))
      for (pnlg = LGSTART(theNode); pnlg != NULL; pnlg = pnlg->next)
         MSET4NLG(sum,COEFF_NLGP(pnlg,Z),NDLGP(NBNODE(pnlg),x))
      for (pnfl=NFSTART(theNode); pnfl != NULL; pnfl=NEXT(pnfl))
         SET4(sum,COEFF_NFP(pnfl,Z),FD(NBFACE(pnfl),x))
      SET1(NDD(theNode,y),sum)
   }
   for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd){
         MSET2LG(sum,COEFF_LGP(theNode,Z),NDLGP(theNode,x))
         for (plg = theNode->lgd->lgstart; plg != NULL; plg = plg->next)
            MSET4LG(sum,COEFF_LGLP(plg,Z),NDLGP(NBNODE(plg),x))
         for (plgn = theNode->lgd->nstart; plgn != NULL; plgn = plgn->next){
            pn = NBNODE(plgn);
            MSET4LGN(sum,COEFF_LGNP(plgn,Z),NDD(pn,x))
         }
         for (plgf = theNode->lgd->fstart; plgf != NULL; plgf = plgf->next)
            SET4LG(sum,COEFF_LGFP(plgf,Z),FD(NBFACE(plgf),x))
         SET1LG(NDLGP(theNode,y),sum)
      }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
   
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl))
         sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      for (pflg = LGSTART(theFace); pflg != NULL; pflg = NEXT(pflg))
         sumf += DOT(COEFF_FLGP(pflg,Z),NDLGP(NBNODE(pflg),x));
      for (pfnl=FNSTART(theFace); pfnl != NULL; pfnl=NEXT(pfnl))
         sumf += DOT(COEFF_FNP(pfnl,Z),NDD(NBNODE(pfnl),x));
      FD(theFace,y) = sumf;
   }
   return(0);
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
               (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
               (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
               (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
               (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) &&
               (DATA_STR & LG_DATA))  */

void vset_mat_value_vn_sf_lg(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: vset_mat_value_vn_sf_lg not available.\n");  }

INT vmultA_vn_sf_lg(tGrid,Z,x,y)
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: vmultA_vn_sf_lg not available.\n"); return(0);  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES)

void set_mat_value_snXvn(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   NODE *pnode;
   LINK *plink;
  
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      SET7(COEFFBP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         SET7(COEFFBLP(plink,Z),r)
   }
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES))  */

void set_mat_value_snXvn(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_snXvn not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (E_DATA & ExDN_MATR)

void set_mat_value_sn_X_vn_ve(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pelem;
   NODE *pnode;
   LINK *plink;
   INT i;
  
   for (pnode=FIRSTN(theGrid);pnode!=NULL;pnode=pnode->succ){
      SET7(COEFFBP(pnode,Z),r)
      for (plink=pnode->tstart;plink!=NULL;plink=plink->next)
         SET7(COEFFBLP(plink,Z),r)
   }
   for (pelem = FIRSTELEMENT(theGrid);pelem != NULL;pelem = pelem->succ)
      for (i=0; i < NVERT; i++)
         SET7(COEFF_BNP(pelem,Z,i),r)
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && 
                (E_DATA & ExDN_MATR))  */

void set_mat_value_sn_X_vn_ve(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_sn_X_vn_ve not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (F_DATA & SCALAR_FACE_DATA)

void mult_by_inverse_face_diagonal_sf(tGrid,Z,x,y)/* y := face_diag(Z)^{-1} x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace;
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,y) = FD(theFace,x)/COEFF_FF(theFace,Z);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
	    FD(theFace,y) = FD(theFace,x)/COEFF_FF(theFace,Z);
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (F_DATA & SCALAR_FACE_DATA))  */

void mult_by_inverse_face_diagonal_sf(tGrid,Z,x,y)/* y := face_diag(Z)^{-1} x */
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: mult_by_inverse_face_diagonal_sf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA)

void multA_f(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   FACE *theFace;
   FLINK *pfl, *stop;
   DOUBLE sumf;
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
      stop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != stop; pfl=NEXT(pfl))
         sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      FD(theFace,y) = sumf;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sumf = COEFF_FF(theFace,Z)*FD(theFace,x);
            stop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != stop; pfl=NEXT(pfl))
               sumf += COEFF_FL(pfl,Z)*FD(ltop_face(NBFACE(pfl),tGrid),x);
            FD(theFace,y) = sumf;
         }
   }
}

void sSOR_step_forward_sf(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   FACE *theFace;
   FLINK *pfl;
   DOUBLE sum;
	
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      if (fabs(COEFF_FF(theFace,Z)) > EPSA){
         sum = FD(theFace,b);
         for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
            sum -= COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
         FD(theFace,x) += om*(sum/COEFF_FF(theFace,Z)-FD(theFace,x));
      }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & SCALAR_FACE_DATA))  */

void multA_f(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{ eprintf("Error: multA_f not available.\n"); }

void sSOR_step_forward_sf(tGrid,Z,b,x,om)
GRID *tGrid; INT Z, b, x; FLOAT om;
{ eprintf("Error: sSOR_step_forward_sf not available.\n"); }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)

void set_mat_value_f(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   FACE *pface;
   FLINK *pflink;
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ) {
      COEFF_FF(pface,Z) = r;
      SET_FF
   }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES))  */

void set_mat_value_f(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_f not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA)

void smultA_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; 
INT Z,x,y,t;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl, *fstop;
   DOUBLE sum[DIM];
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      SET2(sum,FDVP(theFace,x),COEFF_FF(theFace,Z))
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         SET4(sum,FDVP(NBFACE(pfl),x),COEFF_FL(pfl,Z))
      SET1(FDVP(theFace,y),sum)
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            SET2(sum,FDVP(theFace,x),COEFF_FF(theFace,Z))
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl)){
               pface = ltop_face(NBFACE(pfl),tGrid);
               SET4(sum,FDVP(pface,x),COEFF_FL(pfl,Z))
            }
            SET1(FDVP(theFace,y),sum)
         }
   }
}

void smultAL_vf(tGrid,Z,x,y,q) /*  y := Z_lower_triangle x  (without diag)  */
GRID *tGrid;
INT Z, x, y, q;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i;

   for (theGrid = tGrid; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < DIM; i++){
               sumf = 0.;
               for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                  FIND_LTOP_FACE(pface,pfl,tGrid)
                  if (FDV(pface,q,i) < 0.)
                     sumf += FDV(pface,x,i)*COEFF_FL(pfl,Z);
               }  
               FDV(theFace,y,i) = sumf;
               FDV(theFace,q,i) = -1.;
            }
}

void smultAU_vf(tGrid,Z,x,y,q) /*  y := Z_upper_triangle x  (without diag)  */
GRID *tGrid;
INT Z, x, y, q;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i;

   for (theGrid = tGrid; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < DIM; i++){
               sumf = 0.;
               for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                  FIND_LTOP_FACE(pface,pfl,tGrid)
                  if (FDV(pface,q,i) > 0.)
                     sumf += FDV(pface,x,i)*COEFF_FL(pfl,Z);
               }  
               FDV(theFace,y,i) = sumf;
               FDV(theFace,q,i) = -1.;
            }
}

void sSOR_step_forward_vf(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i;

   for (theGrid = tGrid; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < DIM; i++){
               sumf = FDV(theFace,b,i);
               for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                  FIND_LTOP_FACE(pface,pfl,tGrid)
                  sumf -= FDV(pface,x,i)*COEFF_FL(pfl,Z);
               }  
               FDV(theFace,x,i) += 
                                 om*(sumf/COEFF_FF(theFace,Z)-FDV(theFace,x,i));
            }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & VECTOR_FACE_DATA))  */

void smultA_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: smultA_vf not available.\n");  }

void smultAL_vf(tGrid,Z,x,y,q) /*  y := Z_lower_triangle x  (without diag)  */
GRID *tGrid; INT Z, x, y, q;
{  eprintf("Error: smultAL_vf not available.\n");  }

void smultAU_vf(tGrid,Z,x,y,q) /*  y := Z_upper_triangle x  (without diag)  */
GRID *tGrid; INT Z, x, y, q;
{  eprintf("Error: smultAU_vf not available.\n");  }

void sSOR_step_forward_vf(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: sSOR_step_forward_vf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)

void print_smultA_vf(tGrid,Z)
GRID *tGrid;
INT Z;
{
   INT i=0;
   FACE *theFace;
   FLINK *pfl;
   FILE *fp;

   fp = fopen("maticeA_P1NC","w");

   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      theFace->index2 = ++i;
   TNF = i;

   fprintf(fp,"%i\n",2*i);
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      fprintf(fp,"%i %i %e\n",theFace->index2,theFace->index2,COEFF_FF(theFace,Z));
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",theFace->index2,NBFACE(pfl)->index2,COEFF_FL(pfl,Z));
   }
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      fprintf(fp,"%i %i %e\n",theFace->index2+i,theFace->index2+i,COEFF_FF(theFace,Z));
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",theFace->index2+i,NBFACE(pfl)->index2+i,COEFF_FL(pfl,Z));
   }
   fclose(fp);
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES))  */

void print_smultA_vf(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: print_smultA_vf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA) && (DATA_S & PREVIOUS_FACE)

void sSOR_step_backward_vf(tGrid,Z,b,x,om)      /*  one backward step of SOR  */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   GRID *theGrid;
   FACE *theFace, *pface, *stop;
   FLINK *pfl;
   DOUBLE sumf;
   INT i;

   for (theGrid = tGrid->first_grid; theGrid; theGrid = theGrid->finer){
      stop = PREV(FIRSTFACE(theGrid));
      for (theFace=LASTFACE(theGrid); theFace != stop; theFace=PREV(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=DIM1; i > -1; i--){
               sumf = FDV(theFace,b,i);
               for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                  FIND_LTOP_FACE(pface,pfl,tGrid)
                  sumf -= FDV(pface,x,i)*COEFF_FL(pfl,Z);
               }  
               FDV(theFace,x,i) += 
                                 om*(sumf/COEFF_FF(theFace,Z)-FDV(theFace,x,i));
            }
   }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & VECTOR_FACE_DATA) && (DATA_S & PREVIOUS_FACE))  */

void sSOR_step_backward_vf(tGrid,Z,b,x,om)     /*  one backward step of SOR   */
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: sSOR_step_backward_vf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA)

void smult_diagA_vf(tGrid,Z,x,y)  /*  y := diag(Z) x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      SET2(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            SET2(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))
}

void smult_inv_diagA_vf(tGrid,Z,x,y)  /*  y := diag(Z)^{-1} x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      SET22(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            SET22(FDVP(theFace,y),FDVP(theFace,x),COEFF_FF(theFace,Z))
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA))  */

void smult_diagA_vf(tGrid,Z,x,y)
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: smult_diagA_vf not available.\n");  }

void smult_inv_diagA_vf(tGrid,Z,x,y)
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: smult_inv_diagA_vf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA)

void vmultA_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl, *stop;
   DOUBLE sumf[DIM];
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      MSET2(sumf,COEFFNNP(theFace,Z),FDVP(theFace,x))
      stop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != stop; pfl=NEXT(pfl))
         MSET4(sumf,COEFFNNP(pfl,Z),FDVP(NBFACE(pfl),x))
      SET1(FDVP(theFace,y),sumf)
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            MSET2(sumf,COEFFNNP(theFace,Z),FDVP(theFace,x))
            stop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != stop; pfl=NEXT(pfl))
               MSET4(sumf,COEFFNNP(pfl,Z),FDVP(ltop_face(NBFACE(pfl),tGrid),x))
            SET1(FDVP(theFace,y),sumf)
         }
   }
}

void vSOR_step_forward_vf(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sum;
   INT i, j;
	
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      for (i=0; i < 2; i++)
         if (fabs(COEFFNN(theFace,Z,i,i)) > EPSA){
            j = 1 - i;
            sum = FDV(theFace,b,i) - COEFFNN(theFace,Z,i,j)*FDV(theFace,x,j);
            for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
               sum -= DOT(COEFFNNPP(pfl,Z,i),FDVP(NBFACE(pfl),x));
            FDV(theFace,x,i) += om*(sum/COEFFNN(theFace,Z,i,i)-FDV(theFace,x,i));
         }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & VECTOR_FACE_DATA))  */

void vmultA_vf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: vmultA_vf not available.\n");  }

void vSOR_step_forward_vf(tGrid,Z,b,x,om)       /*  one forward step of SOR   */
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: vSOR_step_forward_vf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)

void vset_mat_value_vf(theGrid,Z,r)
GRID *theGrid;
INT Z;
FLOAT r;
{
   FACE *pface;
   FLINK *pflink;
  
   for (pface=FIRSTF(theGrid);pface!=NULL;pface=pface->succ){
      MSET7(COEFFNNP(pface,Z),r)
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)
         MSET7(COEFFNNP(pflink,Z),r)
   }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES))  */

void vset_mat_value_vf(theGrid,Z,r)
GRID *theGrid; INT Z; FLOAT r;
{  eprintf("Error: vset_mat_value_vf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA)

void vmult_diagA_vf(tGrid,Z,x,y)  /*  y := diag(Z) x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      MSET3(FDVP(theFace,y),COEFFNNP(theFace,Z),FDVP(theFace,x))

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MSET3(FDVP(theFace,y),COEFFNNP(theFace,Z),FDVP(theFace,x))
}

void vmult_inv_diagA_vf(tGrid,Z,x,y)  /*  y := diag(Z)^{-1} x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      MSET33(FDVP(theFace,y),COEFFNNP(theFace,Z),FDVP(theFace,x))

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MSET33(FDVP(theFace,y),COEFFNNP(theFace,Z),FDVP(theFace,x))
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA))  */

void vmult_diagA_vf(tGrid,Z,x,y)
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: vmult_diagA_vf not available.\n");  }

void vmult_inv_diagA_vf(tGrid,Z,x,y)
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: vmult_inv_diagA_vf not available.\n");  }

#endif

#if (DIM == 2) && (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)

void set_mat_value_dvf(tGrid,Z,r)
GRID *tGrid;
FLOAT r;
INT Z;
{
   FACE *pface;
   FLINK *pflink;

   for (pface=FIRSTF(tGrid);pface!=NULL;pface=pface->succ){
      COEFFNN(pface,Z,0,0) = COEFFNN(pface,Z,0,1) = 0.;
      COEFFNN(pface,Z,1,0) = COEFFNN(pface,Z,1,1) = 0.;
      for (pflink=pface->tfstart;pflink!=NULL;pflink=pflink->next)
         COEFFNN(pflink,Z,0,0) = COEFFNN(pflink,Z,0,1) =
         COEFFNN(pflink,Z,1,0) = COEFFNN(pflink,Z,1,1) = 0.;
   }
}

#else

void set_mat_value_dvf(tGrid,Z,r)
GRID *tGrid; FLOAT r; INT Z;
{  eprintf("Error: set_mat_value_dvf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & DVECTOR_FACE_DATA)

void vmultA_dvf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid;
INT Z,x,y,t;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl, *fstop;
   DOUBLE sumf[2][DIM];
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sumf[0][0] = COEFFNN(theFace,Z,0,0)*FDDV(theFace,x,0,0) +
                   COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,0);
      sumf[1][0] = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,0) +
                   COEFFNN(theFace,Z,1,1)*FDDV(theFace,x,1,0);
      sumf[0][1] = COEFFNN(theFace,Z,0,0)*FDDV(theFace,x,0,1) +
                   COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,1);
      sumf[1][1] = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,1) +
                   COEFFNN(theFace,Z,1,1)*FDDV(theFace,x,1,1);
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl)){
         pface = NBFACE(pfl);
         sumf[0][0] += COEFFNN(pfl,Z,0,0)*FDDV(pface,x,0,0) +
                       COEFFNN(pfl,Z,0,1)*FDDV(pface,x,1,0);
         sumf[1][0] += COEFFNN(pfl,Z,1,0)*FDDV(pface,x,0,0) +
                       COEFFNN(pfl,Z,1,1)*FDDV(pface,x,1,0);
         sumf[0][1] += COEFFNN(pfl,Z,0,0)*FDDV(pface,x,0,1) +
                       COEFFNN(pfl,Z,0,1)*FDDV(pface,x,1,1);
         sumf[1][1] += COEFFNN(pfl,Z,1,0)*FDDV(pface,x,0,1) +
                       COEFFNN(pfl,Z,1,1)*FDDV(pface,x,1,1);
      } 
      FDDV(theFace,y,0,0) = sumf[0][0];
      FDDV(theFace,y,1,0) = sumf[1][0];
      FDDV(theFace,y,0,1) = sumf[0][1];
      FDDV(theFace,y,1,1) = sumf[1][1];
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sumf[0][0] = COEFFNN(theFace,Z,0,0)*FDDV(theFace,x,0,0) +
                         COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,0);
            sumf[1][0] = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,0) +
                         COEFFNN(theFace,Z,1,1)*FDDV(theFace,x,1,0);
            sumf[0][1] = COEFFNN(theFace,Z,0,0)*FDDV(theFace,x,0,1) +
                         COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,1);
            sumf[1][1] = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,1) +
                         COEFFNN(theFace,Z,1,1)*FDDV(theFace,x,1,1);
            fstop = STOP_FF(theFace,t);
            for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl)){
               pface = ltop_face(NBFACE(pfl),tGrid);
               sumf[0][0] += COEFFNN(pfl,Z,0,0)*FDDV(pface,x,0,0) +
                             COEFFNN(pfl,Z,0,1)*FDDV(pface,x,1,0);
               sumf[1][0] += COEFFNN(pfl,Z,1,0)*FDDV(pface,x,0,0) +
                             COEFFNN(pfl,Z,1,1)*FDDV(pface,x,1,0);
               sumf[0][1] += COEFFNN(pfl,Z,0,0)*FDDV(pface,x,0,1) +
                             COEFFNN(pfl,Z,0,1)*FDDV(pface,x,1,1);
               sumf[1][1] += COEFFNN(pfl,Z,1,0)*FDDV(pface,x,0,1) +
                             COEFFNN(pfl,Z,1,1)*FDDV(pface,x,1,1);
            }
            FDDV(theFace,y,0,0) = sumf[0][0];
            FDDV(theFace,y,1,0) = sumf[1][0];
            FDDV(theFace,y,0,1) = sumf[0][1];
            FDDV(theFace,y,1,1) = sumf[1][1];
         }
   }
}

void vmultAL_dvf(tGrid,Z,x,y,q) /*  y := Z_lower_triangle x  (without diag)  */
GRID *tGrid;
INT Z, x, y, q;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i, j;
	
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < 2; i++)
               for (j=0; j < DIM; j++){
                  sumf = 0.;
                  if (i==0 && FDDV(theFace,q,1,j) < 0.)
                     sumf = COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,j);
                  if (i==1 && FDDV(theFace,q,0,j) < 0.)
                     sumf = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,j);
                  for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                     FIND_LTOP_FACE(pface,pfl,tGrid)
                     if (FDDV(pface,q,0,j) < 0.)
                        sumf += COEFFNN(pfl,Z,i,0)*FDDV(pface,x,0,j);
                     if (FDDV(pface,q,1,j) < 0.)
                        sumf += COEFFNN(pfl,Z,i,1)*FDDV(pface,x,1,j);
                  }
                  FDDV(theFace,y,i,j) = sumf;
                  FDDV(theFace,q,i,j) = -1.;
               }
}

void vmultAU_dvf(tGrid,Z,x,y,q) /*  y := Z_upper_triangle x  (without diag)  */
GRID *tGrid;
INT Z, x, y, q;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i, j;
	
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < 2; i++)
               for (j=0; j < DIM; j++){
                  sumf = 0.;
                  if (i==0 && FDDV(theFace,q,1,j) > 0.)
                     sumf = COEFFNN(theFace,Z,0,1)*FDDV(theFace,x,1,j);
                  if (i==1 && FDDV(theFace,q,0,j) > 0.)
                     sumf = COEFFNN(theFace,Z,1,0)*FDDV(theFace,x,0,j);
                  for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                     FIND_LTOP_FACE(pface,pfl,tGrid)
                     if (FDDV(pface,q,0,j) > 0.)
                        sumf += COEFFNN(pfl,Z,i,0)*FDDV(pface,x,0,j);
                     if (FDDV(pface,q,1,j) > 0.)
                        sumf += COEFFNN(pfl,Z,i,1)*FDDV(pface,x,1,j);
                  }
                  FDDV(theFace,y,i,j) = sumf;
                  FDDV(theFace,q,i,j) = -1.;
               }
}

void vSOR_step_forward_dvf(tGrid,Z,b,x,om)      /*  one forward step of SOR   */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   FLINK *pfl;
   DOUBLE sumf;
   INT i, j;
	
   for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=0; i < 2; i++)
               for (j=0; j < DIM; j++){
                  sumf = FDDV(theFace,b,i,j) -
                                  COEFFNN(theFace,Z,i,0)*FDDV(theFace,x,0,j) -
                                  COEFFNN(theFace,Z,i,1)*FDDV(theFace,x,1,j) +
                                  COEFFNN(theFace,Z,i,i)*FDDV(theFace,x,i,j);
                  for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                     FIND_LTOP_FACE(pface,pfl,tGrid)
                     sumf -= COEFFNN(pfl,Z,i,0)*FDDV(pface,x,0,j) +
                             COEFFNN(pfl,Z,i,1)*FDDV(pface,x,1,j);
                  }
                  FDDV(theFace,x,i,j) += 
                           om*(sumf/COEFFNN(theFace,Z,i,i)-FDDV(theFace,x,i,j));
               }
}

DOUBLE vnorm_of_A_dvf(tGrid,Z) /* L infinity norm of Z */
GRID *tGrid;
INT Z;
{
   GRID *theGrid;
   FACE *theFace;
   FLINK *pfl;
   DOUBLE sum0, sum1, max=0.;
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      sum0 = fabs(COEFFNN(theFace,Z,0,0)) + fabs(COEFFNN(theFace,Z,0,1));
      sum1 = fabs(COEFFNN(theFace,Z,1,0)) + fabs(COEFFNN(theFace,Z,1,1));
      for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
         sum0 += fabs(COEFFNN(pfl,Z,0,0)) + fabs(COEFFNN(pfl,Z,0,1));
         sum1 += fabs(COEFFNN(pfl,Z,1,0)) + fabs(COEFFNN(pfl,Z,1,1));
      } 
      if (sum0 > max) max = sum0;
      if (sum1 > max) max = sum1;
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            sum0 = fabs(COEFFNN(theFace,Z,0,0)) + fabs(COEFFNN(theFace,Z,0,1));
            sum1 = fabs(COEFFNN(theFace,Z,1,0)) + fabs(COEFFNN(theFace,Z,1,1));
            for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
               sum0 += fabs(COEFFNN(pfl,Z,0,0)) + fabs(COEFFNN(pfl,Z,0,1));
               sum1 += fabs(COEFFNN(pfl,Z,1,0)) + fabs(COEFFNN(pfl,Z,1,1));
            }
            if (sum0 > max) max = sum0;
            if (sum1 > max) max = sum1;
         }
   }
   return(max);
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & VECTOR_FACE_DATA))  */

void vmultA_dvf(tGrid,Z,x,y,t) /*  y := Z x */
GRID *tGrid; INT Z,x,y,t;
{  eprintf("Error: vmultA_dvf not available.\n");  }

void vmultAL_dvf(tGrid,Z,x,y,q) /*  y := Z_lower_triangle x  (without diag)  */
GRID *tGrid; INT Z, x, y, q;
{  eprintf("Error: vmultAL_dvf not available.\n");  }

void vmultAU_dvf(tGrid,Z,x,y,q) /*  y := Z_upper_triangle x  (without diag)  */
GRID *tGrid; INT Z, x, y, q;
{  eprintf("Error: vmultAU_dvf not available.\n");  }

void vSOR_step_forward_dvf(tGrid,Z,b,x,om)      /*  one forward step of SOR   */
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: vSOR_step_forward_dvf not available.\n");  }

DOUBLE vnorm_of_A_dvf(tGrid,Z) /* L infinity norm of Z */
GRID *tGrid; INT Z;
{  eprintf("Error: vnorm_of_A_dvf not available.\n"); return(0.);  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)

void printf_vmultA_dvf(tGrid,Z)
GRID *tGrid;
INT Z;
{
   INT i=0, j;
   GRID *theGrid;
   FACE *tF, *pf;
   FLINK *pfl;
   FILE *fp;

   fp = fopen("maticeA_P1MOD","w");

   for (tF=FIRSTFACE(tGrid); tF != NULL; tF=SUCC(tF))
      tF->index2 = ++i;
   TNF = j = 2*i;

   fprintf(fp,"%i\n",4*i);

   for (tF=FIRSTFACE(tGrid); tF != NULL; tF=SUCC(tF)){
      fprintf(fp,"%i %i %e\n",tF->index2  ,tF->index2  ,COEFFNN(tF,Z,0,0));
      fprintf(fp,"%i %i %e\n",tF->index2  ,tF->index2+i,COEFFNN(tF,Z,0,1));
      fprintf(fp,"%i %i %e\n",tF->index2+i,tF->index2  ,COEFFNN(tF,Z,1,0));
      fprintf(fp,"%i %i %e\n",tF->index2+i,tF->index2+i,COEFFNN(tF,Z,1,1));
      for (pfl=FSTART(tF); pfl != NULL; pfl=NEXT(pfl)){
         pf = NBFACE(pfl);
         fprintf(fp,"%i %i %e\n",tF->index2  ,pf->index2  ,COEFFNN(pfl,Z,0,0));
         fprintf(fp,"%i %i %e\n",tF->index2  ,pf->index2+i,COEFFNN(pfl,Z,0,1));
         fprintf(fp,"%i %i %e\n",tF->index2+i,pf->index2  ,COEFFNN(pfl,Z,1,0));
         fprintf(fp,"%i %i %e\n",tF->index2+i,pf->index2+i,COEFFNN(pfl,Z,1,1));
      }
   }
   for (tF=FIRSTFACE(tGrid); tF != NULL; tF=SUCC(tF)){
      fprintf(fp,"%i %i %e\n",j+tF->index2  ,j+tF->index2  ,COEFFNN(tF,Z,0,0));
      fprintf(fp,"%i %i %e\n",j+tF->index2  ,j+tF->index2+i,COEFFNN(tF,Z,0,1));
      fprintf(fp,"%i %i %e\n",j+tF->index2+i,j+tF->index2  ,COEFFNN(tF,Z,1,0));
      fprintf(fp,"%i %i %e\n",j+tF->index2+i,j+tF->index2+i,COEFFNN(tF,Z,1,1));
      for (pfl=FSTART(tF); pfl != NULL; pfl=NEXT(pfl)){
        pf = NBFACE(pfl);
        fprintf(fp,"%i %i %e\n",j+tF->index2  ,j+pf->index2  ,COEFFNN(pfl,Z,0,0));
        fprintf(fp,"%i %i %e\n",j+tF->index2  ,j+pf->index2+i,COEFFNN(pfl,Z,0,1));
        fprintf(fp,"%i %i %e\n",j+tF->index2+i,j+pf->index2  ,COEFFNN(pfl,Z,1,0));
        fprintf(fp,"%i %i %e\n",j+tF->index2+i,j+pf->index2+i,COEFFNN(pfl,Z,1,1));
      }
   }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES))  */

void printf_vmultA_dvf(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: printf_vmultA_dvf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & DVECTOR_FACE_DATA) && (DATA_S & PREVIOUS_FACE)

void vSOR_step_backward_dvf(tGrid,Z,b,x,om)     /*  one backward step of SOR  */
GRID *tGrid;                                    /*  for the solution of Zx=b  */
INT Z, b, x;                                    /*  x ... initial guess       */
FLOAT om;
{
   GRID *theGrid;
   FACE *theFace, *pface, *stop;
   FLINK *pfl;
   DOUBLE sumf;
   INT i, j;
	
   for (theGrid = tGrid->first_grid; theGrid; theGrid = theGrid->finer){
      stop = PREV(FIRSTFACE(theGrid));
      for (theFace=LASTFACE(theGrid); theFace != stop; theFace=PREV(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            for (i=1; i > -1; i--)
               for (j=DIM1; j > -1; j--){
                  sumf = FDDV(theFace,b,i,j) -
                         COEFFNN(theFace,Z,i,0)*FDDV(theFace,x,0,j) -
                         COEFFNN(theFace,Z,i,1)*FDDV(theFace,x,1,j) +
                         COEFFNN(theFace,Z,i,i)*FDDV(theFace,x,i,j);
                  for (pfl=FSTART(theFace); pfl != NULL; pfl=NEXT(pfl)){
                     FIND_LTOP_FACE(pface,pfl,tGrid)
                     sumf -= COEFFNN(pfl,Z,i,0)*FDDV(pface,x,0,j) +
                             COEFFNN(pfl,Z,i,1)*FDDV(pface,x,1,j);
                  }
                  FDDV(theFace,x,i,j) += 
                     om*(sumf/COEFFNN(theFace,Z,i,i)-FDDV(theFace,x,i,j));
               }
   }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & VECTOR_FACE_DATA) && (DATA_S & PREVIOUS_FACE))  */

void vSOR_step_backward_dvf(tGrid,Z,b,x,om)    /*  one backward step of SOR   */
GRID *tGrid; INT Z, b, x; FLOAT om;
{  eprintf("Error: vSOR_step_backward_dvf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (F_DATA & DVECTOR_FACE_DATA)

void vmult_diagA_dvf(tGrid,Z,x,y) /*  y := diag(Z) x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   DOUBLE sumf[2][DIM];
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      FDDV(theFace,y,0,0) = FDDV(theFace,x,0,0)*COEFFNN(theFace,Z,0,0);
      FDDV(theFace,y,1,0) = FDDV(theFace,x,1,0)*COEFFNN(theFace,Z,1,1);
      FDDV(theFace,y,0,1) = FDDV(theFace,x,0,1)*COEFFNN(theFace,Z,0,0);
      FDDV(theFace,y,1,1) = FDDV(theFace,x,1,1)*COEFFNN(theFace,Z,1,1);
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            FDDV(theFace,y,0,0) = FDDV(theFace,x,0,0)*COEFFNN(theFace,Z,0,0);
            FDDV(theFace,y,1,0) = FDDV(theFace,x,1,0)*COEFFNN(theFace,Z,1,1);
            FDDV(theFace,y,0,1) = FDDV(theFace,x,0,1)*COEFFNN(theFace,Z,0,0);
            FDDV(theFace,y,1,1) = FDDV(theFace,x,1,1)*COEFFNN(theFace,Z,1,1);
         }
   }
}

void vmult_inv_diagA_dvf(tGrid,Z,x,y) /*  y := diag(Z)^{-1} x */
GRID *tGrid;
INT Z,x,y;
{
   GRID *theGrid;
   FACE *theFace, *pface;
   DOUBLE sumf[2][DIM];
	
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      FDDV(theFace,y,0,0) = FDDV(theFace,x,0,0)/COEFFNN(theFace,Z,0,0);
      FDDV(theFace,y,1,0) = FDDV(theFace,x,1,0)/COEFFNN(theFace,Z,1,1);
      FDDV(theFace,y,0,1) = FDDV(theFace,x,0,1)/COEFFNN(theFace,Z,0,0);
      FDDV(theFace,y,1,1) = FDDV(theFace,x,1,1)/COEFFNN(theFace,Z,1,1);
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theFace=FIRSTFACE(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            FDDV(theFace,y,0,0) = FDDV(theFace,x,0,0)/COEFFNN(theFace,Z,0,0);
            FDDV(theFace,y,1,0) = FDDV(theFace,x,1,0)/COEFFNN(theFace,Z,1,1);
            FDDV(theFace,y,0,1) = FDDV(theFace,x,0,1)/COEFFNN(theFace,Z,0,0);
            FDDV(theFace,y,1,1) = FDDV(theFace,x,1,1)/COEFFNN(theFace,Z,1,1);
         }
   }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (F_DATA & DVECTOR_FACE_DATA))  */

void vmult_diagA_dvf(tGrid,Z,x,y) /*  y := diag(Z) x */
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: vmult_diagA_dvf not available.\n");  }

void vmult_inv_diagA_dvf(tGrid,Z,x,y) /*  y := diag(Z)^{-1} x */
GRID *tGrid; INT Z,x,y;
{  eprintf("Error: vmult_inv_diagA_dvf not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR)

#if DIM == 3

void set_mat_value_e_X_vn_sf(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      SET7(COEFF_BNP(pelem,Z,0),r)
      SET7(COEFF_BNP(pelem,Z,1),r)
      SET7(COEFF_BNP(pelem,Z,2),r)
      SET7(COEFF_BNP(pelem,Z,3),r)
      COEFF_BF(pelem,Z,0) = r;
      COEFF_BF(pelem,Z,1) = r;
      COEFF_BF(pelem,Z,2) = r;
      COEFF_BF(pelem,Z,3) = r;
   } 
}

#else

void set_mat_value_e_X_vn_sf(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      SET7(COEFF_BNP(pelem,Z,0),r)
      SET7(COEFF_BNP(pelem,Z,1),r)
      SET7(COEFF_BNP(pelem,Z,2),r)
      COEFF_BF(pelem,Z,0) = r;
      COEFF_BF(pelem,Z,1) = r;
      COEFF_BF(pelem,Z,2) = r;
   } 
}

#endif

#else

void set_mat_value_e_X_vn_sf(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_e_X_vn_sf not available.\n");  }

#endif

#if (E_DATA & E_E_FNEIGHBOURS) && (E_DATA & ExE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA)

void multA_se_X_se(tGrid,Z,x,y)
GRID *tGrid;
INT Z, x, y;
{
   ELEMENT *pelem, *pel;
   INT i;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      ED(pelem,y) = COEFF_EE(pelem,Z)*ED(pelem,x);
      for (i = 0; i < SIDES; i++)
         if (pel=NB_EL(pelem,i))
            ED(pelem,y) += COEFF_BF(pelem,Z,i)*ED(pel,x);
   }
}

#else

void multA_se_X_se(tGrid,Z,x,y)
GRID *tGrid; INT Z, x, y;
{  eprintf("Error: multA_se_X_se not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

#if DIM == 3

void multA_e_X_vn_sf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in nodes, scalars in faces */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_CN(pel->n[0],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
     if (IS_CN(pel->n[1],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
     if (IS_CN(pel->n[2],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
     if (IS_CN(pel->n[3],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,3),NDD(pel->n[3],x));
     if (IS_CF(pel->f[0],t)) ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
     if (IS_CF(pel->f[1],t)) ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
     if (IS_CF(pel->f[2],t)) ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
     if (IS_CF(pel->f[3],t)) ED(pel,ye) += COEFF_BF(pel,Z,3)*FD(pel->f[3],x);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_CN(pel->n[0],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
            if (IS_CN(pel->n[1],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
            if (IS_CN(pel->n[2],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
            if (IS_CN(pel->n[3],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,3),NDD(pel->n[3],x));
            if (IS_CF(pel->f[0],t))
               ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
            if (IS_CF(pel->f[1],t))
               ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
            if (IS_CF(pel->f[2],t))
               ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
            if (IS_CF(pel->f[3],t))
               ED(pel,ye) += COEFF_BF(pel,Z,3)*FD(pel->f[3],x);
         }
}
 
void multA_vn_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;            /*  columns ... scalars in elements                   */
                        /*  rows ... vectors in nodes, scalars in faces       */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   vs_set_value_nf(tGrid,0.,y,1);
   vs_b_copy_nf(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      SET4(NDD(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
      FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
      FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
      FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
      FD(pel->f[3],y) += COEFF_BF(pel,Z,3)*ED(pel,xe);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
            SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
            SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
            SET4(NDD(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
            FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
            FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
            FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
            FD(pel->f[3],y) += COEFF_BF(pel,Z,3)*ED(pel,xe);
         }
   vs_b_copy_nf(tGrid,q,y);
}

#else  /*  if DIM != 3  */

void multA_e_X_vn_sf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in nodes, scalars in faces */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_CN(pel->n[0],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
     if (IS_CN(pel->n[1],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
     if (IS_CN(pel->n[2],t)) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
     if (IS_CF(pel->f[0],t)) ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
     if (IS_CF(pel->f[1],t)) ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
     if (IS_CF(pel->f[2],t)) ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_CN(pel->n[0],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
            if (IS_CN(pel->n[1],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
            if (IS_CN(pel->n[2],t))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
            if (IS_CF(pel->f[0],t))
               ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
            if (IS_CF(pel->f[1],t))
               ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
            if (IS_CF(pel->f[2],t))
               ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
         }
}
 
void multA_vn_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;            /*  columns ... scalars in elements                   */
                        /*  rows ... vectors in nodes, scalars in faces       */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   vs_set_value_nf(tGrid,0.,y,1);
   vs_b_copy_nf(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
      FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
      FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
            SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
            SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
            FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
            FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
            FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
         }
   vs_b_copy_nf(tGrid,q,y);
}

#endif

#else  /*  if !((E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && 
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                (E_DATA & SCALAR_ELEMENT_DATA))  */

void multA_e_X_vn_sf(tGrid,Z,x,ye,t) /* ye := A x */
GRID *tGrid; INT Z, x, ye, t; 
{  eprintf("Error: multA_e_X_vn_sf not available.\n");  }

void multA_vn_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; INT Z, xe, y, q, t;
{  eprintf("Error: multA_vn_sf_X_e not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (E_DATA & ExDF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (DIM == 2)

void set_mat_value_e_X_vn_vf(tGrid,Z,r)
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z;                     /*  columns ... vectors in nodes, vectors in faces */
FLOAT r;
{
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET7(COEFF_BNP(pel,Z,0),0.)
      SET7(COEFF_BNP(pel,Z,1),0.)
      SET7(COEFF_BNP(pel,Z,2),0.)
      SET7(COEFF_BDFP(pel,Z,0),0.)
      SET7(COEFF_BDFP(pel,Z,1),0.)
      SET7(COEFF_BDFP(pel,Z,2),0.)
   }
}

void multA_e_X_vn_vf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in nodes, vectors in faces */
{
   ELEMENT *pel;
  
   if (!(t & STOP_IS_FIRST_INNER))
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (IS_FN(pel->n[0])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
         if (IS_FN(pel->n[1])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
         if (IS_FN(pel->n[2])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
         if (IS_FF(pel->f[0])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (IS_FF(pel->f[1])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (IS_FF(pel->f[2])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
      }
   else
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (NOT_FN(pel->n[0])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
         if (NOT_FN(pel->n[1])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
         if (NOT_FN(pel->n[2])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
         if (NOT_FF(pel->f[0])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (NOT_FF(pel->f[1])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (NOT_FF(pel->f[2])) ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
      }
}
 
void multA_vn_vf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;            /*  columns ... scalars in elements                   */
                        /*  rows ... vectors in nodes, vectors in faces       */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   ELEMENT *pel;
   
   set_value(tGrid,0.,y,1,Q_VNVF);
   b_copy(tGrid,y,q,1,Q_VNVF);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      SET4(FDVP(pel->f[0],y),COEFF_BDFP(pel,Z,0),ED(pel,xe))
      SET4(FDVP(pel->f[1],y),COEFF_BDFP(pel,Z,1),ED(pel,xe))
      SET4(FDVP(pel->f[2],y),COEFF_BDFP(pel,Z,2),ED(pel,xe))
   }
   b_copy(tGrid,q,y,1,Q_VNVF);
}

#else

void set_mat_value_e_X_vn_vf(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_e_X_vn_vf not available.\n");  }

void multA_e_X_vn_vf(tGrid,Z,x,ye,t)
GRID *tGrid; INT Z, x, ye, t;
{  eprintf("Error: multA_e_X_vn_vf not available.\n");  }

void multA_vn_vf_X_e(tGrid,Z,xe,y,q,t) 
GRID *tGrid; INT Z, xe, y, q, t;   
{  eprintf("Error: multA_vn_vf_X_e not available.\n");  }

#endif

#if E_DATA & ExDN_MATR

#if DIM == 3

void set_mat_value_e_X_vn(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      SET7(COEFF_BNP(pelem,Z,0),r)
      SET7(COEFF_BNP(pelem,Z,1),r)
      SET7(COEFF_BNP(pelem,Z,2),r)
      SET7(COEFF_BNP(pelem,Z,3),r)
   } 
}

#else

void set_mat_value_e_X_vn(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      SET7(COEFF_BNP(pelem,Z,0),r)
      SET7(COEFF_BNP(pelem,Z,1),r)
      SET7(COEFF_BNP(pelem,Z,2),r)
   } 
}

#endif

#else

void set_mat_value_e_X_vn(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_e_X_vn not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (N_DATA & VECTOR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

#if DIM == 3

void multA_e_X_vn(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in nodes                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_FN(pel->n[0])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
     if (IS_FN(pel->n[1])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
     if (IS_FN(pel->n[2])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
     if (IS_FN(pel->n[3])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,3),NDD(pel->n[3],x));
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_FN(pel->n[0])) 
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
            if (IS_FN(pel->n[1]))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
            if (IS_FN(pel->n[2]))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
            if (IS_FN(pel->n[3]))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,3),NDD(pel->n[3],x));
         }
}

void multA_vn_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;  /*  columns ... scalars in elements, rows ... vectors in nodes  */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   v_set_value_n(tGrid,0.,y,1);
   v_b_copy_n(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      SET4(NDD(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
            SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
            SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
            SET4(NDD(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
         }
   v_b_copy_n(tGrid,q,y);
}

#else  /*  if DIM != 3  */

void multA_e_X_vn(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in nodes                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_FN(pel->n[0])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
     if (IS_FN(pel->n[1])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
     if (IS_FN(pel->n[2])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_FN(pel->n[0])) 
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
            if (IS_FN(pel->n[1]))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
            if (IS_FN(pel->n[2]))
               ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
         }
}

void multA_vn_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;  /*  columns ... scalars in elements, rows ... vectors in nodes  */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   v_set_value_n(tGrid,0.,y,1);
   v_b_copy_n(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
            SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
            SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
         }
   v_b_copy_n(tGrid,q,y);
}

#endif

#else  /*  if !((E_DATA & ExDN_MATR) && (N_DATA & VECTOR_NODE_DATA) &&
                (E_DATA & SCALAR_ELEMENT_DATA))  */

void multA_e_X_vn(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid; INT Z, x, ye, t;
{  eprintf("Error: multA_e_X_vn not available.\n");  }

void multA_vn_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; INT Z, xe, y, q, t;
{  eprintf("Error: multA_vn_X_e not available.\n");  }

#endif

#if (E_DATA & ExF_MATR) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

#if DIM == 3

void multA_e_X_sf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... scalars in faces                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_FF(pel->f[0]))
        ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
     if (IS_FF(pel->f[1]))
        ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
     if (IS_FF(pel->f[2]))
        ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
     if (IS_FF(pel->f[3]))
        ED(pel,ye) += COEFF_BF(pel,Z,3)*FD(pel->f[3],x);
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_FF(pel->f[0]))
               ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
            if (IS_FF(pel->f[1]))
               ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
            if (IS_FF(pel->f[2]))
               ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
            if (IS_FF(pel->f[3]))
               ED(pel,ye) += COEFF_BF(pel,Z,3)*FD(pel->f[3],x);
         }
}

void multA_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;  /*  columns ... scalars in elements, rows ... scalars in faces  */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   set_value_f(tGrid,0.,y,1);
   b_copy_f(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
      FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
      FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
      FD(pel->f[3],y) += COEFF_BF(pel,Z,3)*ED(pel,xe);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
            FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
            FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
            FD(pel->f[3],y) += COEFF_BF(pel,Z,3)*ED(pel,xe);
         }
   b_copy_f(tGrid,q,y);
}

#else  /*  if DIM != 3  */

void multA_e_X_sf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... scalars in faces                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_FF(pel->f[0]))
        ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
     if (IS_FF(pel->f[1]))
        ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
     if (IS_FF(pel->f[2]))
        ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_FF(pel->f[0]))
               ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
            if (IS_FF(pel->f[1]))
               ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
            if (IS_FF(pel->f[2]))
               ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
         }
}

void multA_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid;  /*  columns ... scalars in elements, rows ... scalars in faces  */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   set_value_f(tGrid,0.,y,1);
   b_copy_f(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
      FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
      FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
            FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
            FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
         }
   b_copy_f(tGrid,q,y);
}

#endif

#else  /*  if !((E_DATA & ExF_MATR) && (F_DATA & SCALAR_FACE_DATA) &&
                (E_DATA & SCALAR_ELEMENT_DATA))  */

void multA_e_X_sf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid; INT Z, x, ye, t;
{  eprintf("Error: multA_e_X_sf not available.\n");  }

void multA_sf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; INT Z, xe, y, q, t;
{  eprintf("Error: multA_sf_X_e not available.\n");  }

#endif

#if (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

void set_mat_value_e_X_vf(tGrid,Z,r)
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z;                     /*  columns ... vectors in faces                   */
FLOAT r;
{
   ELEMENT *pel;
 
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET7(COEFF_BDFP(pel,Z,0),0.)
      SET7(COEFF_BDFP(pel,Z,1),0.)
      SET7(COEFF_BDFP(pel,Z,2),0.)
   }
}

void multA_e_X_vf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in faces                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   if (!(t & STOP_IS_FIRST_INNER)){
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (IS_FF(pel->f[0]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (IS_FF(pel->f[1]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (IS_FF(pel->f[2]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
      }
      for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
         for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
            if (IS_LTOP_ELEMENT(pel,tGrid)){
               ED(pel,ye) = 0.0;
               if (IS_FF(pel->f[0]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
               if (IS_FF(pel->f[1]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
               if (IS_FF(pel->f[2]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
            }
   }
   else{
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (NOT_FF(pel->f[0]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (NOT_FF(pel->f[1]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (NOT_FF(pel->f[2]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
      }
      for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
         for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
            if (IS_LTOP_ELEMENT(pel,tGrid)){
               ED(pel,ye) = 0.0;
               if (NOT_FF(pel->f[0]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
               if (NOT_FF(pel->f[1]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
               if (NOT_FF(pel->f[2]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
            }
   }
}

void multA_vf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; /*  columns ... scalars in elements, rows ... vectors in faces   */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   vset_value_f(tGrid,0.,y,1);
   vb_copy_f(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(FDVP(pel->f[0],y),COEFF_BDFP(pel,Z,0),ED(pel,xe));
      SET4(FDVP(pel->f[1],y),COEFF_BDFP(pel,Z,1),ED(pel,xe));
      SET4(FDVP(pel->f[2],y),COEFF_BDFP(pel,Z,2),ED(pel,xe));
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(FDVP(pel->f[0],y),COEFF_BDFP(pel,Z,0),ED(pel,xe));
            SET4(FDVP(pel->f[1],y),COEFF_BDFP(pel,Z,1),ED(pel,xe));
            SET4(FDVP(pel->f[2],y),COEFF_BDFP(pel,Z,2),ED(pel,xe));
         }
   vb_copy_f(tGrid,q,y);
}

#elif (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == CUBE)

void set_mat_value_e_X_vf(tGrid,Z,r)
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z;                     /*  columns ... vectors in faces                   */
FLOAT r;
{
   ELEMENT *pel;
 
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET7(COEFF_BDFP(pel,Z,0),0.)
      SET7(COEFF_BDFP(pel,Z,1),0.)
      SET7(COEFF_BDFP(pel,Z,2),0.)
      SET7(COEFF_BDFP(pel,Z,3),0.)
   }
}

void multA_e_X_vf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... vectors in faces                   */
{
   GRID *theGrid;
   ELEMENT *pel;
  
   if (!(t & STOP_IS_FIRST_INNER)){
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (IS_FF(pel->f[0]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (IS_FF(pel->f[1]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (IS_FF(pel->f[2]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
         if (IS_FF(pel->f[3]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,3),FDVP(pel->f[3],x));
      }
      for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
         for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
            if (IS_LTOP_ELEMENT(pel,tGrid)){
               ED(pel,ye) = 0.0;
               if (IS_FF(pel->f[0]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
               if (IS_FF(pel->f[1]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
               if (IS_FF(pel->f[2]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
               if (IS_FF(pel->f[3]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,3),FDVP(pel->f[3],x));
            }
   }
   else{
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         ED(pel,ye) = 0.0;
         if (NOT_FF(pel->f[0]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
         if (NOT_FF(pel->f[1]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
         if (NOT_FF(pel->f[2]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
         if (NOT_FF(pel->f[3]))
            ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,3),FDVP(pel->f[3],x));
      }
      for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
         for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
            if (IS_LTOP_ELEMENT(pel,tGrid)){
               ED(pel,ye) = 0.0;
               if (NOT_FF(pel->f[0]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDVP(pel->f[0],x));
               if (NOT_FF(pel->f[1]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDVP(pel->f[1],x));
               if (NOT_FF(pel->f[2]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDVP(pel->f[2],x));
               if (NOT_FF(pel->f[3]))
                  ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,3),FDVP(pel->f[3],x));
            }
   }
}

void multA_vf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; /*  columns ... scalars in elements, rows ... vectors in faces   */
INT Z, xe, y, q, t;           /* only boundary components of q used */
{
   GRID *theGrid;
   ELEMENT *pel;
   
   vset_value_f(tGrid,0.,y,1);
   vb_copy_f(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(FDVP(pel->f[0],y),COEFF_BDFP(pel,Z,0),ED(pel,xe));
      SET4(FDVP(pel->f[1],y),COEFF_BDFP(pel,Z,1),ED(pel,xe));
      SET4(FDVP(pel->f[2],y),COEFF_BDFP(pel,Z,2),ED(pel,xe));
      SET4(FDVP(pel->f[3],y),COEFF_BDFP(pel,Z,3),ED(pel,xe));
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(FDVP(pel->f[0],y),COEFF_BDFP(pel,Z,0),ED(pel,xe));
            SET4(FDVP(pel->f[1],y),COEFF_BDFP(pel,Z,1),ED(pel,xe));
            SET4(FDVP(pel->f[2],y),COEFF_BDFP(pel,Z,2),ED(pel,xe));
            SET4(FDVP(pel->f[3],y),COEFF_BDFP(pel,Z,3),ED(pel,xe));
         }
   vb_copy_f(tGrid,q,y);
}

#else  /*  if !((E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && 
                (F_DATA & VECTOR_FACE_DATA))  */

void set_mat_value_e_X_vf(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{ eprintf("Error: set_mat_value_e_X_vf not available.\n"); }
 
void multA_e_X_vf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid; INT Z, x, ye, t;
{ eprintf("Error: multA_e_X_vf not available.\n"); }

void multA_vf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; INT Z, xe, y, q, t;
{ eprintf("Error: multA_vf_X_e not available.\n"); }

#endif

#if (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA)

void print_multA_e_X_vf(tGrid,Z)
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z;                     /*  columns ... vectors in faces                   */
{
 INT i=0;
 GRID *theGrid;
 ELEMENT *pel;
 FILE *fp;

 fp = fopen("maticeB_P1NC_P0","w");
 
 for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
    pel->index2 = ++i;

 fprintf(fp,"%i\n",i);
 for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
  if (IS_FF(pel->f[0])){
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[0]->index2,    COEFF_BDF(pel,Z,0,0));
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[0]->index2+TNF,COEFF_BDF(pel,Z,0,1));
  }
  if (IS_FF(pel->f[1])){
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[1]->index2,    COEFF_BDF(pel,Z,1,0));
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[1]->index2+TNF,COEFF_BDF(pel,Z,1,1));
  }
  if (IS_FF(pel->f[2])){
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[2]->index2,    COEFF_BDF(pel,Z,2,0));
   fprintf(fp,"%i %i %e\n",pel->index2,pel->f[2]->index2+TNF,COEFF_BDF(pel,Z,2,1));
  }
 }
 fclose(fp);
}

#else  /*  if !((E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA))  */

void print_multA_e_X_vf(tGrid,Z)
GRID *tGrid; INT Z;
{ eprintf("Error: print_multA_e_X_vf not available.\n"); }

#endif

#if (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void multA_e_X_dvf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid;               /*  rows ... scalars in elements                   */
INT Z, x, ye, t;           /*  columns ... dvectors in faces                  */
{                          /*  only first component of dvector multiplied     */
   GRID *theGrid;
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_CF(pel->f[0],t))
        ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDDVP(pel->f[0],x,0));
     if (IS_CF(pel->f[1],t))
        ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDDVP(pel->f[1],x,0));
     if (IS_CF(pel->f[2],t))
        ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDDVP(pel->f[2],x,0));
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            ED(pel,ye) = 0.0;
            if (IS_CF(pel->f[0],t))
               ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,0),FDDVP(pel->f[0],x,0));
            if (IS_CF(pel->f[1],t))
               ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,1),FDDVP(pel->f[1],x,0));
            if (IS_CF(pel->f[2],t))
               ED(pel,ye) += DOT(COEFF_BDFP(pel,Z,2),FDDVP(pel->f[2],x,0));
         }
}

void multA_dvf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe                            */
GRID *tGrid; /*  columns ... scalars in elements, rows ... dvectors in faces  */
INT Z, xe, y, q, t;           /* only first comp. of dvectors used  */
{                             /* only boundary components of q used */
   GRID *theGrid;
   ELEMENT *pel;
   
   dvset_value_f(tGrid,0.,y,1);
   dvb_copy_f(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(FDDVP(pel->f[0],y,0),COEFF_BDFP(pel,Z,0),ED(pel,xe));
      SET4(FDDVP(pel->f[1],y,0),COEFF_BDFP(pel,Z,1),ED(pel,xe));
      SET4(FDDVP(pel->f[2],y,0),COEFF_BDFP(pel,Z,2),ED(pel,xe));
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            SET4(FDDVP(pel->f[0],y,0),COEFF_BDFP(pel,Z,0),ED(pel,xe));
            SET4(FDDVP(pel->f[1],y,0),COEFF_BDFP(pel,Z,1),ED(pel,xe));
            SET4(FDDVP(pel->f[2],y,0),COEFF_BDFP(pel,Z,2),ED(pel,xe));
         }
   dvb_copy_f(tGrid,q,y);
}

#else  /*  if !((E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && 
                (F_DATA & DVECTOR_FACE_DATA))  */

void multA_e_X_dvf(tGrid,Z,x,ye,t)  /* ye := A x */
GRID *tGrid; INT Z, x, ye, t;
{ eprintf("Error: multA_e_X_dvf not available.\n"); }

void multA_dvf_X_e(tGrid,Z,xe,y,q,t)  /* y := A xe */
GRID *tGrid; INT Z, xe, y, q, t;
{ eprintf("Error: multA_dvf_X_e not available.\n"); }

#endif

#if (F_DATA & ONE_FACE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA)

void diag_BAmBT(tGrid,ZA,Z,ze)  /* ze := diag(B diag(A^{-1}) B^T) */
GRID *tGrid;
INT ZA,Z,ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      ED(pel,ze) = 0.;
      for (i = 0; i < DIM2; i++){
         if (IS_FN(pel->n[i]))
            ED(pel,ze) += DOT(COEFF_BNP(pel,Z,i),COEFF_BNP(pel,Z,i))
                                                          /COEFFN(pel->n[i],ZA);
         if (IS_FF(pel->f[i]))
            ED(pel,ze) += COEFF_BF(pel,Z,i)*COEFF_BF(pel,Z,i)/COEFF_FF(pel->f[i],ZA);
      }
   }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (E_DATA & ExDN_MATR) && 
                (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA))  */

void diag_BAmBT(tGrid,ZA,Z,ze)  /* ze := diag(B diag(A^{-1}) B^T) */
GRID *tGrid; INT ZA,Z,ze;
{  eprintf("Error: diag_BAmBT not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (N_DATA & SCALAR_NODE_DATA)

void diag_BAmBT_P2C_P1C(tGrid,ZA,Z,d)  
GRID *tGrid;                             /* d := diag(B diag(A^{-1}) B^T) */
INT ZA, Z, d;
{
   NODE *pnode, *pn;
   FACE *pf;
   LINK *pli;
   NFLINK *pnf;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ){
      NDS(pnode,d) = 0.;
      if (IS_FN(pnode))
         NDS(pnode,d)+=DOT(COEFFBP(pnode,Z),COEFFBP(pnode,Z))/COEFFN(pnode,ZA);
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn=NBNODE(pli);
         NDS(pnode,d)+=DOT(COEFFBLP(pli,Z),COEFFBLP(pli,Z))/COEFFN(pn,ZA);
      }
      for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next){
         pf=NBFACE(pnf);
         NDS(pnode,d)+=DOT(COEFF_NFP(pnf,Z),COEFF_NFP(pnf,Z))/COEFF_FF(pf,ZA);
      }
   }
//   NDS(FIRSTN(tGrid),d) = 1.;
}

#else 

void diag_BAmBT_P2C_P1C(tGrid,ZA,Z,d)  
GRID *tGrid; INT ZA, Z, d;
{  eprintf("Error: diag_BAmBT_P2C_P1C not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

#if DIM == 3

void diag_rBTB(tGrid,r,ZA,Z,x)   /*  x := diag(r B^T B)  */
GRID *tGrid;
FLOAT r;
INT ZA, Z, x;
{
   ELEMENT *pel;
   NODE *theNode;
   FACE *theFace;
   	
/* for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode))
      ND(theNode,x,0) = ND(theNode,x,1) = ND(theNode,x,2) = COEFFN(theNode,ZA);
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,x) = COEFF_FF(theFace,ZA);*/
   vs_set_value_nf(tGrid,0.,x,1);
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ){
     if (IS_FN(pel->n[0])){ 
        ND(pel->n[0],x,0) += r*COEFF_BN(pel,Z,0,0)*COEFF_BN(pel,Z,0,0); 
        ND(pel->n[0],x,1) += r*COEFF_BN(pel,Z,0,1)*COEFF_BN(pel,Z,0,1); 
        ND(pel->n[0],x,2) += r*COEFF_BN(pel,Z,0,2)*COEFF_BN(pel,Z,0,2);
     }
     if (IS_FN(pel->n[1])){ 
        ND(pel->n[1],x,0) += r*COEFF_BN(pel,Z,1,0)*COEFF_BN(pel,Z,1,0); 
        ND(pel->n[1],x,1) += r*COEFF_BN(pel,Z,1,1)*COEFF_BN(pel,Z,1,1); 
        ND(pel->n[1],x,2) += r*COEFF_BN(pel,Z,1,2)*COEFF_BN(pel,Z,1,2);
     }
     if (IS_FN(pel->n[2])){ 
        ND(pel->n[2],x,0) += r*COEFF_BN(pel,Z,2,0)*COEFF_BN(pel,Z,2,0); 
        ND(pel->n[2],x,1) += r*COEFF_BN(pel,Z,2,1)*COEFF_BN(pel,Z,2,1); 
        ND(pel->n[2],x,2) += r*COEFF_BN(pel,Z,2,2)*COEFF_BN(pel,Z,2,2);
     }
     if (IS_FN(pel->n[3])){ 
        ND(pel->n[3],x,0) += r*COEFF_BN(pel,Z,3,0)*COEFF_BN(pel,Z,3,0); 
        ND(pel->n[3],x,1) += r*COEFF_BN(pel,Z,3,1)*COEFF_BN(pel,Z,3,1); 
        ND(pel->n[3],x,2) += r*COEFF_BN(pel,Z,3,2)*COEFF_BN(pel,Z,3,2);
     }
     if (IS_FF(pel->f[0]))
        FD(pel->f[0],x) += r*COEFF_BF(pel,Z,0)*COEFF_BF(pel,Z,0);
     if (IS_FF(pel->f[1]))
        FD(pel->f[1],x) += r*COEFF_BF(pel,Z,1)*COEFF_BF(pel,Z,1);
     if (IS_FF(pel->f[2]))
        FD(pel->f[2],x) += r*COEFF_BF(pel,Z,2)*COEFF_BF(pel,Z,2);
     if (IS_FF(pel->f[3]))
        FD(pel->f[3],x) += r*COEFF_BF(pel,Z,3)*COEFF_BF(pel,Z,3);
   }
}

#else  /*  if DIM != 3  */

void diag_rBTB(tGrid,r,ZA,Z,x)   /*  x := diag(r B^T B)  */
GRID *tGrid;
FLOAT r;
INT ZA, Z, x;
{
   ELEMENT *pel;
   NODE *theNode;
   FACE *theFace;

/* for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode))
      ND(theNode,x,0) = ND(theNode,x,1) = COEFFN(theNode,ZA);
   for (theFace=FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,x) = COEFF_FF(theFace,ZA);*/
   vs_set_value_nf(tGrid,0.,x,1);
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ){
     if (IS_FN(pel->n[0])){
        ND(pel->n[0],x,0) += r*COEFF_BN(pel,Z,0,0)*COEFF_BN(pel,Z,0,0);
        ND(pel->n[0],x,1) += r*COEFF_BN(pel,Z,0,1)*COEFF_BN(pel,Z,0,1);
     }
     if (IS_FN(pel->n[1])){
        ND(pel->n[1],x,0) += r*COEFF_BN(pel,Z,1,0)*COEFF_BN(pel,Z,1,0);
        ND(pel->n[1],x,1) += r*COEFF_BN(pel,Z,1,1)*COEFF_BN(pel,Z,1,1);
     }
     if (IS_FN(pel->n[2])){
        ND(pel->n[2],x,0) += r*COEFF_BN(pel,Z,2,0)*COEFF_BN(pel,Z,2,0);
        ND(pel->n[2],x,1) += r*COEFF_BN(pel,Z,2,1)*COEFF_BN(pel,Z,2,1);
     }
     if (IS_FF(pel->f[0]))
        FD(pel->f[0],x) += r*COEFF_BF(pel,Z,0)*COEFF_BF(pel,Z,0);
     if (IS_FF(pel->f[1]))
        FD(pel->f[1],x) += r*COEFF_BF(pel,Z,1)*COEFF_BF(pel,Z,1);
     if (IS_FF(pel->f[2]))
        FD(pel->f[2],x) += r*COEFF_BF(pel,Z,2)*COEFF_BF(pel,Z,2);
   }
}

#endif

#else  /*  if !((E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && 
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)  */

void diag_rBTB(tGrid,r,ZA,Z,x)   /*  x := diag(r B^T B)  */
GRID *tGrid; FLOAT r; INT ZA, Z, x;
{  eprintf("Error: diag_rBTB not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (DIM == 3)

void add_BTB(tGrid,r,ZA,Z)
GRID *tGrid;
FLOAT r;
INT ZA,Z;
{
   ELEMENT *pel;
   LINK *pli;
   NFLINK *pnf;
   FNLINK *pfn;
   FLINK *pfl;
   INT i, j;
   
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel != NULL; pel = pel->succ)
      for (i = 0; i < 4; i++){
         if (IS_FN(pel->n[i])){
            COEFFNN(pel->n[i],ZA,0,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,i,0);
            COEFFNN(pel->n[i],ZA,0,1) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,i,1);
            COEFFNN(pel->n[i],ZA,0,2) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,i,2);
            COEFFNN(pel->n[i],ZA,1,0) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,i,0);
            COEFFNN(pel->n[i],ZA,1,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,i,1);
            COEFFNN(pel->n[i],ZA,1,2) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,i,2);
            COEFFNN(pel->n[i],ZA,2,0) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,i,0);
            COEFFNN(pel->n[i],ZA,2,1) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,i,1);
            COEFFNN(pel->n[i],ZA,2,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,i,2);
            for (j = 0; j < 4; j++){
               if (j != i && IS_FN(pel->n[j])){
                  for (pli=pel->n[i]->start; pli->nbnode!=pel->n[j]; 
                                                               pli=pli->next);
                  COEFFLL(pli,ZA,0,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,j,0);
                  COEFFLL(pli,ZA,0,1) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,j,1);
                  COEFFLL(pli,ZA,0,2) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,j,2);
                  COEFFLL(pli,ZA,1,0) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,j,0);
                  COEFFLL(pli,ZA,1,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,j,1);
                  COEFFLL(pli,ZA,1,2) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,j,2);
                  COEFFLL(pli,ZA,2,0) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,j,0);
                  COEFFLL(pli,ZA,2,1) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,j,1);
                  COEFFLL(pli,ZA,2,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,j,2);
               }
               if (IS_FF(pel->f[j])){
                  for (pnf=pel->n[i]->nfstart; pnf->nbface!=pel->f[j]; 
                                                                pnf=pnf->next);
                  COEFF_NF(pnf,ZA,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BF(pel,Z,j);
                  COEFF_NF(pnf,ZA,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BF(pel,Z,j);
                  COEFF_NF(pnf,ZA,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BF(pel,Z,j);
               }
            }
         }
         if (IS_FF(pel->f[i])){
            COEFF_FF(pel->f[i],ZA) += r*COEFF_BF(pel,Z,i)*COEFF_BF(pel,Z,i);
            for (j = 0; j < 4; j++){
               if (j != i && IS_FF(pel->f[j])){
                  for (pfl=pel->f[i]->fstart; pfl->nbface!=pel->f[j]; 
                                                                pfl=pfl->next);
                  COEFF_FL(pfl,ZA) += r*COEFF_BF(pel,Z,i)*COEFF_BF(pel,Z,j);
               }
               if (IS_FN(pel->n[j])){
                  for (pfn=pel->f[i]->fnstart; pfn->nbnode!=pel->n[j]; 
                                                                pfn=pfn->next);
                  COEFF_FN(pfn,ZA,0) += r*COEFF_BF(pel,Z,i)*COEFF_BN(pel,Z,j,0);
                  COEFF_FN(pfn,ZA,1) += r*COEFF_BF(pel,Z,i)*COEFF_BN(pel,Z,j,1);
                  COEFF_FN(pfn,ZA,2) += r*COEFF_BF(pel,Z,i)*COEFF_BN(pel,Z,j,2);
               }
            }
         }
      }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
                (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
                (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
                (DIM == 3))  */

void add_BTB(tGrid,r,ZA,Z)
GRID *tGrid; FLOAT r; INT ZA,Z;
{  eprintf("Error: add_BTB not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (DATA_S & N_LINK_TO_NODES) && (DIM == 3)

void add_BTB_prep(tGrid,r,ZA,Z)
GRID *tGrid;
FLOAT r;
INT ZA,Z;
{
   ELEMENT *pel;
   LINK *pli;
   INT i, j;
   
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel != NULL; pel = pel->succ)
      for (i = 0; i < 4; i++){
         if (IS_FN(pel->n[i])){
            COEFFNN(pel->n[i],ZA,0,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,i,0);
            COEFFNN(pel->n[i],ZA,1,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,i,1);
            COEFFNN(pel->n[i],ZA,2,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,i,2);
            for (j = 0; j < 4; j++){
               if (j != i && IS_FN(pel->n[j])){
                  for (pli=pel->n[i]->start; pli->nbnode!=pel->n[j]; 
                                                               pli=pli->next);
                  COEFFLL(pli,ZA,0,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,j,0);
                  COEFFLL(pli,ZA,1,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,j,1);
                  COEFFLL(pli,ZA,2,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,j,2);
               }
            }
         }
         if (IS_FF(pel->f[i])){
            COEFF_FF(pel->f[i],ZA) += r*COEFF_BF(pel,Z,i)*COEFF_BF(pel,Z,i);
         }
      }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && 
                (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DIM == 3))  */

void add_BTB_prep(tGrid,r,ZA,Z)
GRID *tGrid; FLOAT r; INT ZA,Z;
{  eprintf("Error: add_BTB_prep not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (DIM == 3)

void add_BTB_prep2(tGrid,r,ZA,Z)
GRID *tGrid;
FLOAT r;
INT ZA,Z;
{
   ELEMENT *pel;
   INT i, j;

   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel != NULL; pel = pel->succ)
      for (i = 0; i < 4; i++){
         if (IS_FN(pel->n[i])){
            COEFFNN(pel->n[i],ZA,0,0) += r*COEFF_BN(pel,Z,i,0)*COEFF_BN(pel,Z,i,0);
            COEFFNN(pel->n[i],ZA,1,1) += r*COEFF_BN(pel,Z,i,1)*COEFF_BN(pel,Z,i,1);
            COEFFNN(pel->n[i],ZA,2,2) += r*COEFF_BN(pel,Z,i,2)*COEFF_BN(pel,Z,i,2);
         }
         if (IS_FF(pel->f[i])){
            COEFF_FF(pel->f[i],ZA) += r*COEFF_BF(pel,Z,i)*COEFF_BF(pel,Z,i);
         }
      }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && 
                (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (DIM == 3))  */

void add_BTB_prep2(tGrid,r,ZA,Z)
GRID *tGrid; FLOAT r; INT ZA,Z;
{  eprintf("Error: add_BTB_prep2 not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES)

void setA_sn_X_vn_vf(tGrid,Z,r)  /*  A_ij := r  */
GRID *tGrid;              /*  rows ... scalars in nodes                       */
INT Z;                    /*  columns ... vectors in nodes, vectors in faces  */
FLOAT r;
{
   NODE *pnode;
   LINK *pli;
   NFLINK *pnf;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ){
      SET7(COEFFBP(pnode,Z),r)
      for (pli = TSTART(pnode); pli != NULL; pli = pli->next)
         SET7(COEFFBLP(pli,Z),r)
      for (pnf = TNFSTART(pnode); pnf != NULL; pnf = pnf->next)
         SET7(COEFF_NFP(pnf,Z),r)
   }
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES))  */

void setA_sn_X_vn_vf(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: setA_sn_X_vn_vf not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void multA_sn_X_vn_vf(tGrid,Z,x,y,t)  /*  y := A x  */
GRID *tGrid;              /*  rows ... scalars in nodes                       */
INT Z, x, y, t;           /*  columns ... vectors in nodes, vectors in faces  */
{
   NODE *pnode, *pn;
   FACE *pf;
   LINK *pli, *stop;
   NFLINK *pnf, *nfstop;

   if (!(t & STOP_IS_FIRST_INNER))
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ){
         NDS(pnode,y) = 0.;
         for (pli = pnode->start; pli != NULL; pli = pli->next){
            pn=NBNODE(pli);
            NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
         }
         for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next){
            pf=NBFACE(pnf);
            NDS(pnode,y) += DOT(COEFF_NFP(pnf,Z),FDVP(pf,x));
         }
      }
   else
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ){
         NDS(pnode,y) = DOT(COEFFBP(pnode,Z),NDD(pnode,x));
         for (pli = pnode->tstart; pli != NULL; pli = pli->next)
            if (NOT_FN(pn=NBNODE(pli)))
               NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
         for (pnf = pnode->tnfstart; pnf != NULL; pnf = pnf->next)
            if (NOT_FF(pf=NBFACE(pnf)))
               NDS(pnode,y) += DOT(COEFF_NFP(pnf,Z),FDVP(pf,x));
      }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ){
      NDS(pnode,y) = DOT(COEFFBP(pnode,Z),NDD(pnode,x));
      stop = STOP_NN(pnode,t);
      for (pli = T_START(pnode,t); pli != stop; pli = pli->next){
         pn = NBNODE(pli);
         NDS(pnode,y) += DOT(COEFFBLP(pli,Z),NDD(pn,x));
      }
      nfstop = STOP_NF(pnode,t);
      for (pnf = NF_START(pnode,t); pnf != nfstop; pnf = pnf->next){
         pf = NBFACE(pnf);
         NDS(pnode,y) += DOT(COEFF_NFP(pnf,Z),FDVP(pf,x));
      }
   }
}

void multA_vn_vf_X_sn(tGrid,Z,x,y)  /*  y = A x  */
GRID *tGrid;                /*  rows ... vectors in nodes, vectors in faces   */
INT Z, x, y;                /*  columns ... scalars in nodes                  */
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pli;
   NFLINK *pnf;

   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      SET_VALUE(pnode,0.,y)
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      SET_VALUE(pface,0.,y)
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ){
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn=NBNODE(pli);
         SET4(NDD(pn,y),COEFFBLP(pli,Z),NDS(pnode,x))
      }
      for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next){
         pf=NBFACE(pnf);
         SET4(FDVP(pf,y),COEFF_NFP(pnf,Z),NDS(pnode,x))
      }
   }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ){
      SET4(NDD(pnode,y),COEFFBP(pnode,Z),NDS(pnode,x))
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn = NBNODE(pli);
         SET4(NDD(pn,y),COEFFBLP(pli,Z),NDS(pnode,x))
      }
      for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next){
         pf = NBFACE(pnf);
         SET4(FDVP(pf,y),COEFF_NFP(pnf,Z),NDS(pnode,x))
      }
   }
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && 
                (F_DATA & VECTOR_FACE_DATA))  */

void multA_sn_X_vn_vf(tGrid,Z,x,y,t)
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sn_X_vn_vf not available.\n");  }

void multA_vn_vf_X_sn(tGrid,Z,x,y)
GRID *tGrid; INT Z, x, y;
{  eprintf("Error: multA_vn_vf_X_sn not available.\n");  }

#endif

void set_mat_value_ssn_ssf_X_vn_vf(tGrid,Z,r)  /*  Z := r  */
GRID *tGrid;      /*  rows ... vectors in nodes, vectors in faces             */
INT Z;            /*  columns ... scalars in special nodes and special faces  */
FLOAT r;
{
   SNODE *psn;
   SFACE *psf;

   for (psn = FIRSTSN(tGrid); psn; psn = psn->succ)
      psn->ann0[Z][0] = psn->ann1[Z][0] = psn->ann2[Z][0] =
      psn->ann0[Z][1] = psn->ann1[Z][1] = psn->ann2[Z][1] =
      psn->anf1[Z][0] = psn->anf2[Z][0] =
      psn->anf1[Z][1] = psn->anf2[Z][1] = r;
   for (psf = FIRSTSF(tGrid); psf; psf = psf->succ)
      psf->aff[Z][0] = psf->afn1[Z][0] = psf->afn2[Z][0] =
      psf->aff[Z][1] = psf->afn1[Z][1] = psf->afn2[Z][1] = r;
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t)  /*  y := A x  */
GRID *tGrid;         /*  rows ... scalars in special nodes and special faces  */
INT Z, x, y, t;      /*  columns ... vectors in nodes, vectors in faces       */
{
   SNODE *psn;
   SFACE *psf;

   for (psn = FIRSTSN(tGrid); psn; psn = psn->succ)
      SNDS(psn,y) = 0.;
   for (psf = FIRSTSF(tGrid); psf; psf = psf->succ)
      SFDS(psf,y) = 0.;
   if (!(t & ONLY_INNER)){
      for (psn = FIRSTSN(tGrid); psn; psn = psn->succ){
         if (NOT_FN(psn->n0)) SNDS(psn,y) += DOT(psn->ann0[Z],NDD(psn->n0,x));
         if (NOT_FN(psn->n1)) SNDS(psn,y) += DOT(psn->ann1[Z],NDD(psn->n1,x));
         if (NOT_FN(psn->n2)) SNDS(psn,y) += DOT(psn->ann2[Z],NDD(psn->n2,x));
         if (NOT_FF(psn->f1)) SNDS(psn,y) += DOT(psn->anf1[Z],FDVP(psn->f1,x));
         if (NOT_FF(psn->f2)) SNDS(psn,y) += DOT(psn->anf2[Z],FDVP(psn->f2,x));
      }
      for (psf = FIRSTSF(tGrid); psf; psf = psf->succ){
         if (NOT_FF(psf->f )) SFDS(psf,y) += DOT(psf->aff[Z],FDVP(psf->f,x));
         if (NOT_FN(psf->n1)) SFDS(psf,y) += DOT(psf->afn1[Z],NDD(psf->n1,x));
         if (NOT_FN(psf->n2)) SFDS(psf,y) += DOT(psf->afn2[Z],NDD(psf->n2,x));
      }
   }
   if (!(t & STOP_IS_FIRST_INNER)){
      for (psn = FIRSTSN(tGrid); psn; psn = psn->succ){
         if (IS_FN(psn->n0)) SNDS(psn,y) += DOT(psn->ann0[Z],NDD(psn->n0,x));
         if (IS_FN(psn->n1)) SNDS(psn,y) += DOT(psn->ann1[Z],NDD(psn->n1,x));
         if (IS_FN(psn->n2)) SNDS(psn,y) += DOT(psn->ann2[Z],NDD(psn->n2,x));
         if (IS_FF(psn->f1)) SNDS(psn,y) += DOT(psn->anf1[Z],FDVP(psn->f1,x));
         if (IS_FF(psn->f2)) SNDS(psn,y) += DOT(psn->anf2[Z],FDVP(psn->f2,x));
      }
      for (psf = FIRSTSF(tGrid); psf; psf = psf->succ){
         if (IS_FF(psf->f )) SFDS(psf,y) += DOT(psf->aff[Z],FDVP(psf->f,x));
         if (IS_FN(psf->n1)) SFDS(psf,y) += DOT(psf->afn1[Z],NDD(psf->n1,x));
         if (IS_FN(psf->n2)) SFDS(psf,y) += DOT(psf->afn2[Z],NDD(psf->n2,x));
      }
   }
}

void add_multA_vn_vf_X_ssn_ssf(tGrid,Z,x,y)  /*  y := y + A x  */
GRID *tGrid;      /*  rows ... vectors in nodes, vectors in faces             */
INT Z, x, y;      /*  columns ... scalars in special nodes and special faces  */
{
   SNODE *psn;
   SFACE *psf;

   for (psn = FIRSTSN(tGrid); psn; psn = psn->succ){
      SET4(NDD(psn->n0,y),psn->ann0[Z],SNDS(psn,x))
      SET4(NDD(psn->n1,y),psn->ann1[Z],SNDS(psn,x))
      SET4(NDD(psn->n2,y),psn->ann2[Z],SNDS(psn,x))
      SET4(FDVP(psn->f1,y),psn->anf1[Z],SNDS(psn,x))
      SET4(FDVP(psn->f2,y),psn->anf2[Z],SNDS(psn,x))
   }
   for (psf = FIRSTSF(tGrid); psf; psf = psf->succ){
      SET4(FDVP(psf->f,y),psf->aff[Z],SFDS(psf,x))
      SET4(NDD(psf->n1,y),psf->afn1[Z],SFDS(psf,x))
      SET4(NDD(psf->n2,y),psf->afn2[Z],SFDS(psf,x))
   }
}

#else

void multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t)  /*  y := A x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_ssn_ssf_X_vn_vf not available.\n");  }

void add_multA_vn_vf_X_ssn_ssf(tGrid,Z,x,y)  /*  y := y + A x  */
GRID *tGrid; INT Z, x, y;
{  eprintf("Error: add_multA_vn_vf_X_ssn_ssf not available.\n");  }

#endif

#if (F_DATA & FACE_BD_DATA) && (N_DATA & NODE_BD_DATA) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES)

void setA_sn_X_vn_vf_BD(tGrid,Z,r)  /*  A_ij := r  */
GRID *tGrid;              /*  rows ... scalars in nodes                       */
INT Z;                    /*  columns ... vectors in nodes, vectors in faces  */
FLOAT r;
{
   NODE *pnode;
   FACE *pface;
   LINK *pli;
   NFLINK *pnf;
   FNLINK *pfn;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ){
      SET7(COEFFBP(pnode,Z),r)
      for (pli = TSTART(pnode); pli != NULL; pli = pli->next)
         SET7(COEFFBLP(pli,Z),r)
      for (pnf = TNFSTART(pnode); pnf != NULL; pnf = pnf->next)
         SET7(COEFF_NFP(pnf,Z),r)
//matrix D:
      if (IS_ZNN(pnode)) {
          SET7(COEFFDP(pnode,Z),r)
          for (pli = TSTART(pnode); pli != NULL; pli = pli->next)
              if (IS_ZNN(pli->nbnode))
				  SET7(COEFFDP(pli,Z),r)
          for (pnf = TNFSTART(pnode); pnf != NULL; pnf = pnf->next)
              if (IS_ZNF(pnf->nbface))
				  SET7(COEFFDP(pnf,Z),r)
      }
   }
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
       if (IS_ZNF(pface)) {
           SET7(COEFFDP(pface,Z),r)
           for (pfn = TFNSTART(pface); pfn != NULL; pfn = pfn->next)
               if (IS_ZNN(pfn->nbnode))
				   SET7(COEFFDP(pfn,Z),r)
       }
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES))  */

void setA_sn_X_vn_vf_BD(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: setA_sn_X_vn_vf_BD not available.\n");  }

#endif

#if (N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void multD_sn_sf_X_vn_vf(tGrid,Z,x,y)	// y = Z x
GRID *tGrid;
INT Z, x, y;
{
	NODE *pnode, *pn;
	FACE *pface, *pf;
	LINK *pli;
	NFLINK *pnf;
	FNLINK *pfn;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode)) {
			NDSBD(pnode,y) = DOT(COEFFDP(pnode,Z),NDD(pnode,x));
			for (pli = START(pnode); pli; pli = pli->next)
				if (IS_ZNN(pn=pli->nbnode))
					NDSBD(pnode,y) += DOT(COEFFDP(pli,Z),NDD(pn,x));
			for (pnf = NFSTART(pnode); pnf; pnf = pnf->next)
				if (IS_ZNF(pf=pnf->nbface))
					NDSBD(pnode,y) += DOT(COEFFDP(pnf,Z),FDVP(pf,x));
		}
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface)) {
			FDBD(pface,y) = DOT(COEFFDP(pface,Z),FDVP(pface,x));
			for (pfn = FNSTART(pface); pfn; pfn = pfn->next)
				if (IS_ZNN(pn=pfn->nbnode))
					FDBD(pface,y) += DOT(COEFFDP(pfn,Z),NDD(pn,x));
		}
}

void multD_vn_vf_X_sn_sf(tGrid,Z,x,y)	// y += Z x
GRID *tGrid;
INT Z, x, y;
{
	NODE *pnode, *pn;
	FACE *pface, *pf;
	LINK *pli;
	NFLINK *pnf;
	FNLINK *pfn;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode)) {
			SET4(NDD(pnode,y),COEFFDP(pnode,Z),NDSBD(pnode,x))
			for (pli = START(pnode); pli; pli = pli->next)
				if (IS_ZNN(pn=pli->nbnode))
					SET4(NDD(pn,y),COEFFDP(pli,Z),NDSBD(pnode,x))
			for (pnf = NFSTART(pnode); pnf; pnf = pnf->next)
				if (IS_ZNF(pf=pnf->nbface))
					SET4(FDVP(pf,y),COEFFDP(pnf,Z),NDSBD(pnode,x))
		}
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface)) {
			SET4(FDVP(pface,y),COEFFDP(pface,Z),FDBD(pface,x))
			for (pfn = FNSTART(pface); pfn; pfn = pfn->next)
				if (IS_ZNN(pn=pfn->nbnode))
					SET4(NDD(pn,y),COEFFDP(pfn,Z),FDBD(pface,x))
		}
}

#else

void multD_sn_sf_X_vn_vf(tGrid,Z,x,y)
GRID *tGrid; INT Z, x, y;
{ eprintf("Error: multD_sn_sf_X_vn_vf not available.\n"); }

void multD_vn_vf_X_sn_sf(tGrid,Z,x,y)
GRID *tGrid; INT Z, x, y;
{ eprintf("Error: multD_vn_vf_X_sn_sf not available.\n"); }

#endif

#if (N_DATA & IxD_NODE_MATR) && (E_DATA & ExDN_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void diag_BAmBT_for_mini(tGrid,ZA,Z,d)  
GRID *tGrid;                             /* d := diag(B diag(A^{-1}) B^T) */
INT ZA, Z, d;
{
   ELEMENT *pel;
   NODE *pnode, *pn;
   LINK *pli;
   INT i;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ){
      NDS(pnode,d) = 0.;
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         pn=NBNODE(pli);
         NDS(pnode,d)+=DOT(COEFFBLP(pli,Z),COEFFBLP(pli,Z))/COEFFN(pn,ZA);
      }
   }
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i=0; i < NVERT; i++)
         NDS(pel->n[i],d) += DOT(COEFF_BNP(pel,Z,i),COEFF_BNP(pel,Z,i))/ED(pel,ZA);
   NDS(FIRSTN(tGrid),d) = 1.;
}

#else  /*  if !((N_DATA & IxD_NODE_MATR) && (E_DATA & ExDN_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && 
                (E_DATA & SCALAR_ELEMENT_DATA))  */

void diag_BAmBT_for_mini(tGrid,Z,Zbn,d)  
GRID *tGrid; INT Z, Zbn, d;
{  eprintf("Error: diag_BAmBT_for_mini not available.\n");  }

#endif

#if (E_DATA & ExFxDN_MATR) && (E_DATA & ExFxDF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)

void set_mat_value_sne_X_vn_vf(tGrid,Z,r)  /*  Z := r */
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      for (i = 0; i < NVERT; i++){
         for (j = 0; j < NVERT; j++)
            SET7(COEFF_BFDNP(pel,Z,i,j),r)
         for (j = 0; j < SIDES; j++)
            SET7(COEFF_BFDFP(pel,Z,i,j),r)
      }
}

void multA_sef_X_vn_vf(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... scalars in element faces                      */
INT Z, x, y, t;    /*  columns ... vectors in nodes and faces                 */
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      for (i = 0; i < NVERT; i++){
         EDSN(pel,y,i) = 0.;
         for (j = 0; j < NVERT; j++)
            if (IS_CN(pel->n[j],t))
               EDSN(pel,y,i) += DOT(COEFF_BFDNP(pel,Z,i,j),NDD(pel->n[j],x));
         for (j = 0; j < SIDES; j++)
            if (IS_CF(pel->f[j],t))
               EDSN(pel,y,i) += DOT(COEFF_BFDFP(pel,Z,i,j),FDVP(pel->f[j],x));
      }
}

void multA_vn_vf_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... vectors in nodes and faces                    */
INT Z, x, y, t;    /*  columns ... scalars in element faces                   */
{
   ELEMENT *pel;
   INT i, j;

   set_value(tGrid,0.,y,1,Q_VNVF);
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      for (i = 0; i < SIDES; i++){
         for (j = 0; j < NVERT; j++)
            SET4(NDD(pel->n[j],y),COEFF_BFDNP(pel,Z,i,j),EDSN(pel,x,i))
         for (j = 0; j < SIDES; j++)
            SET4(FDVP(pel->f[j],y),COEFF_BFDFP(pel,Z,i,j),EDSN(pel,x,i))
      }
}

#else

void set_mat_value_sne_X_vn_vf(tGrid,Z,r)  /*  Z := r */
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_sne_X_vn_vf not available.\n");  }

void multA_sef_X_vn_vf(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sef_X_vn_vf not available.\n");  }

void multA_vn_vf_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_vn_vf_X_sef not available.\n");  }

#endif

#if (E_DATA & ExFxDN_MATR) && (E_DATA & ExFxDF_MATR) && (E_DATA & ExDF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)

void set_mat_value_sne_X_vn_vf_ve(tGrid,Z,r)  /*  Z := r */
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < NVERT; i++){
         SET7(COEFF_BDFP(pel,Z,i),r)
         for (j = 0; j < NVERT; j++)
            SET7(COEFF_BFDNP(pel,Z,i,j),r)
         for (j = 0; j < SIDES; j++)
            SET7(COEFF_BFDFP(pel,Z,i,j),r)
      }
}

void multA_sef_X_vn_vf_ve(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... scalars in element faces                      */
INT Z, x, y, t;    /*  columns ... vectors in nodes, faces and elements       */
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < NVERT; i++){
         if (!(t & STOP_IS_FIRST_INNER))
            EDSN(pel,y,i) = DOT(COEFF_BDFP(pel,Z,i),EDVP(pel,x));
         else
            EDSN(pel,y,i) = 0.;
         for (j = 0; j < NVERT; j++)
            if (IS_CN(pel->n[j],t))
               EDSN(pel,y,i) += DOT(COEFF_BFDNP(pel,Z,i,j),NDD(pel->n[j],x));
         for (j = 0; j < SIDES; j++)
            if (IS_CF(pel->f[j],t))
               EDSN(pel,y,i) += DOT(COEFF_BFDFP(pel,Z,i,j),FDVP(pel->f[j],x));
      }
}

void multA_vn_vf_ve_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... vectors in nodes, faces and elements          */
INT Z, x, y, t;    /*  columns ... scalars in element faces                   */
{
   ELEMENT *pel;
   INT i, j;

   set_value(tGrid,0.,y,1,Q_VNVF);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET7(EDVP(pel,y),0.)
      for (i = 0; i < SIDES; i++){
         for (j = 0; j < NVERT; j++)
            SET4(NDD(pel->n[j],y),COEFF_BFDNP(pel,Z,i,j),EDSN(pel,x,i))
         for (j = 0; j < SIDES; j++)
            SET4(FDVP(pel->f[j],y),COEFF_BFDFP(pel,Z,i,j),EDSN(pel,x,i))
         SET4(EDVP(pel,y),COEFF_BDFP(pel,Z,i),EDSN(pel,x,i))
      }
   }
}

#else

void set_mat_value_sne_X_vn_vf_ve(tGrid,Z,r)  /*  Z := r */
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_sne_X_vn_vf_ve not available.\n");  }

void multA_sef_X_vn_vf_ve(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sef_X_vn_vf_ve not available.\n");  }

void multA_vn_vf_ve_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_vn_vf_ve_X_sef not available.\n");  }

#endif

#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (DATA_STR & LG_DATA) && (DIM == 3)

void multB_lg(tGrid,Z,x,ye)  /* ye := B x */
GRID *tGrid;
INT Z, x, ye;
{
   ELEMENT *pel;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
     ED(pel,ye) = 0.0;
     if (IS_FN(pel->n[0])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,0),NDD(pel->n[0],x));
     else if (pel->n[0]->lgd)
                       ED(pel,ye) += DOTLG(COEFF_BNP(pel,Z,0),NDLGP(pel->n[0],x));
     if (IS_FN(pel->n[1])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,1),NDD(pel->n[1],x));
     else if (pel->n[1]->lgd)
                       ED(pel,ye) += DOTLG(COEFF_BNP(pel,Z,1),NDLGP(pel->n[1],x));
     if (IS_FN(pel->n[2])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,2),NDD(pel->n[2],x));
     else if (pel->n[2]->lgd)
                       ED(pel,ye) += DOTLG(COEFF_BNP(pel,Z,2),NDLGP(pel->n[2],x));
     if (IS_FN(pel->n[3])) ED(pel,ye) += DOT(COEFF_BNP(pel,Z,3),NDD(pel->n[3],x));
     else if (pel->n[3]->lgd)
                       ED(pel,ye) += DOTLG(COEFF_BNP(pel,Z,3),NDLGP(pel->n[3],x));
     if (IS_FF(pel->f[0]))
        ED(pel,ye) += COEFF_BF(pel,Z,0)*FD(pel->f[0],x);
     if (IS_FF(pel->f[1]))
        ED(pel,ye) += COEFF_BF(pel,Z,1)*FD(pel->f[1],x);   
     if (IS_FF(pel->f[2]))
        ED(pel,ye) += COEFF_BF(pel,Z,2)*FD(pel->f[2],x);
     if (IS_FF(pel->f[3]))
        ED(pel,ye) += COEFF_BF(pel,Z,3)*FD(pel->f[3],x);
   }
}
 
void multBT_lg(tGrid,Z,xe,y,q)  /* y := B^T xe */
GRID *tGrid;
INT Z, xe, y, q;              /* only boundary components of q are used */
{
   ELEMENT *pel;
   
   vs_set_value_nflg(tGrid,0.,y,1);
   vs_b_copy_nf(tGrid,y,q);
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      SET4(NDD(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      SET4(NDD(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      SET4(NDD(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      SET4(NDD(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
      if (pel->n[0]->lgd)
         SET4LG(NDLGP(pel->n[0],y),COEFF_BNP(pel,Z,0),ED(pel,xe))
      if (pel->n[1]->lgd)
         SET4LG(NDLGP(pel->n[1],y),COEFF_BNP(pel,Z,1),ED(pel,xe))
      if (pel->n[2]->lgd)
         SET4LG(NDLGP(pel->n[2],y),COEFF_BNP(pel,Z,2),ED(pel,xe))
      if (pel->n[3]->lgd)
         SET4LG(NDLGP(pel->n[3],y),COEFF_BNP(pel,Z,3),ED(pel,xe))
      FD(pel->f[0],y) += COEFF_BF(pel,Z,0)*ED(pel,xe);
      FD(pel->f[1],y) += COEFF_BF(pel,Z,1)*ED(pel,xe);
      FD(pel->f[2],y) += COEFF_BF(pel,Z,2)*ED(pel,xe);
      FD(pel->f[3],y) += COEFF_BF(pel,Z,3)*ED(pel,xe);
   }
   vs_b_copy_nf(tGrid,q,y);
}

void diag_rBTB_lg(tGrid,r,Z,x)   /*  x := diag(r B^T B)  */
GRID *tGrid;
FLOAT r;
INT Z, x;
{
   ELEMENT *pel;
   	
   vs_set_value_nflg(tGrid,0.,x,1);
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ){
      ND(pel->n[0],x,0) += r*COEFF_BN(pel,Z,0,0)*COEFF_BN(pel,Z,0,0); 
      ND(pel->n[0],x,1) += r*COEFF_BN(pel,Z,0,1)*COEFF_BN(pel,Z,0,1); 
      ND(pel->n[0],x,2) += r*COEFF_BN(pel,Z,0,2)*COEFF_BN(pel,Z,0,2);
      ND(pel->n[1],x,0) += r*COEFF_BN(pel,Z,1,0)*COEFF_BN(pel,Z,1,0); 
      ND(pel->n[1],x,1) += r*COEFF_BN(pel,Z,1,1)*COEFF_BN(pel,Z,1,1); 
      ND(pel->n[1],x,2) += r*COEFF_BN(pel,Z,1,2)*COEFF_BN(pel,Z,1,2);
      ND(pel->n[2],x,0) += r*COEFF_BN(pel,Z,2,0)*COEFF_BN(pel,Z,2,0); 
      ND(pel->n[2],x,1) += r*COEFF_BN(pel,Z,2,1)*COEFF_BN(pel,Z,2,1); 
      ND(pel->n[2],x,2) += r*COEFF_BN(pel,Z,2,2)*COEFF_BN(pel,Z,2,2);
      ND(pel->n[3],x,0) += r*COEFF_BN(pel,Z,3,0)*COEFF_BN(pel,Z,3,0); 
      ND(pel->n[3],x,1) += r*COEFF_BN(pel,Z,3,1)*COEFF_BN(pel,Z,3,1); 
      ND(pel->n[3],x,2) += r*COEFF_BN(pel,Z,3,2)*COEFF_BN(pel,Z,3,2);
      if (pel->n[0]->lgd){
         NDLG(pel->n[0],x,0) += r*COEFF_BN(pel,Z,0,0)*COEFF_BN(pel,Z,0,0);
         NDLG(pel->n[0],x,1) += r*COEFF_BN(pel,Z,0,1)*COEFF_BN(pel,Z,0,1);
      }
      if (pel->n[1]->lgd){
         NDLG(pel->n[1],x,0) += r*COEFF_BN(pel,Z,1,0)*COEFF_BN(pel,Z,1,0);
         NDLG(pel->n[1],x,1) += r*COEFF_BN(pel,Z,1,1)*COEFF_BN(pel,Z,1,1);
      }
      if (pel->n[2]->lgd){
         NDLG(pel->n[2],x,0) += r*COEFF_BN(pel,Z,2,0)*COEFF_BN(pel,Z,2,0);
         NDLG(pel->n[2],x,1) += r*COEFF_BN(pel,Z,2,1)*COEFF_BN(pel,Z,2,1);
      }
      if (pel->n[3]->lgd){
         NDLG(pel->n[3],x,0) += r*COEFF_BN(pel,Z,3,0)*COEFF_BN(pel,Z,3,0);
         NDLG(pel->n[3],x,1) += r*COEFF_BN(pel,Z,3,1)*COEFF_BN(pel,Z,3,1);
      }
      FD(pel->f[0],x) += r*COEFF_BF(pel,Z,0)*COEFF_BF(pel,Z,0);
      FD(pel->f[1],x) += r*COEFF_BF(pel,Z,1)*COEFF_BF(pel,Z,1);
      FD(pel->f[2],x) += r*COEFF_BF(pel,Z,2)*COEFF_BF(pel,Z,2);
      FD(pel->f[3],x) += r*COEFF_BF(pel,Z,3)*COEFF_BF(pel,Z,3);
   }
}

#else  /*  if !((E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && 
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                (E_DATA & SCALAR_ELEMENT_DATA)&& (DATA_STR & LG_DATA) && 
                (DIM == 3))  */

void multB_lg(tGrid,Z,x,ye)  /* ye := B x */
GRID *tGrid; INT Z, x, ye;
{  eprintf("Error: multB_lg not available.\n");  }
 
void multBT_lg(tGrid,Z,xe,y,q)  /* y := B^T xe */
GRID *tGrid; INT xe, y, q;          /* only boundary components of q are used */
{  eprintf("Error: multBT_lg not available.\n");  }

void diag_rBTB_lg(tGrid,r,Z,x)   /*  x := diag(r B^T B)  */
GRID *tGrid; FLOAT r; INT Z, x;
{  eprintf("Error: diag_rBTB_lg not available.\n");  }

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void mass_matrix(tGrid,ze)  /* ze_i := volume(T_i) */
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      ED(pel,ze) = VOLUME(pel);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            ED(pel,ze) = VOLUME(pel);
}

void inv_mass_matrix(tGrid,ze)  /* ze_i := 1./volume(T_i) */
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      ED(pel,ze) = 1./VOLUME(pel);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
        if (IS_LTOP_ELEMENT(pel,tGrid))
          ED(pel,ze) = 1./VOLUME(pel);
}

#else  /*  if !(E_DATA & SCALAR_ELEMENT_DATA)  */

void mass_matrix(tGrid,ze)  /* ze_i := volume(T_i) */
GRID *tGrid; INT ze;
{  eprintf("Error: mass_matrix not available.\n");  }

void inv_mass_matrix(tGrid,ze)  /* ze_i := 1./volume(T_i) */
GRID *tGrid; INT ze;
{  eprintf("Error: inv_mass_matrix not available.\n");  }

#endif

#if (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void multA_sn_X_dvf0(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... scalars in nodes                              */
INT Z, x, y, t;    /*  columns ... vectors in faces                           */
{                  /*  matrix multiplied by the first component of dvector    */
   NODE *pnode, *stop;
   NFLINK *pnf, *nfstop;
   FNLINK *pfn;

   stop = STOP_NODE(tGrid,t);
   for (pnode = FIRST_NODE(tGrid,t); pnode != stop; pnode = SUCC(pnode)){
      NDS(pnode,y) = 0.;
      nfstop = STOP_NF(pnode,t);
      for (pnf = NF_START(pnode,t); pnf != nfstop; pnf=NEXT(pnf)){
         for (pfn = FN_START(NBFACE(pnf),t); NBNODE(pfn) != pnode; pfn=NEXT(pfn));
         NDS(pnode,y) += DOT(COEFF_FNP(pfn,Z),FDDVP(NBFACE(pnf),x,0));
      }
   }
}

#else  /*  if !((F_DATA & IxD_FACE_NODE_MATR) && 
                (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
                (N_DATA & SCALAR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void multA_sn_X_dvf0(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sn_X_dvf0 not available.\n");  }

#endif

#if (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void multA_dvf0_X_sn(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... vectors in faces                              */
INT Z, x, y, t;    /*  columns ... scalars in nodes                           */
{                  /*  result is the first component of dvector               */
   FACE *pface, *stop;
   FNLINK *pfn, *fnstop;

   stop = STOP_FACE(tGrid,t);
   for (pface = FIRST_FACE(tGrid,t); pface != stop; pface = SUCC(pface)){
      SET7(FDDVP(pface,y,0),0.)
      fnstop = STOP_FN(pface,t);
      for (pfn = FN_START(pface,t); pfn != fnstop; pfn=NEXT(pfn))
         SET4(FDDVP(pface,y,0),COEFF_FNP(pfn,Z),NDS(NBNODE(pfn),x))
   }
}

#else  /*  if !((F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & F_LINK_TO_NODES) && 
                (N_DATA & SCALAR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void multA_dvf0_X_sn(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_dvf0_X_sn not available.\n");  }

#endif

#if (N_DATA & Ix2D_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_FACES) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void multA_sn_X_dvf(tGrid,Z,x,y,t_row,t_column)  /*  y := Z x  */
GRID *tGrid;                     /*  rows ... scalars in nodes                */
INT Z, x, y, t_row, t_column;    /*  columns ... dvectors in faces            */
{
   NODE *pnode, *stop;
   NFLINK *pnf, *nfstop;

   stop = STOP_NODE(tGrid,t_row);
   for (pnode = FIRST_NODE(tGrid,t_row); pnode != stop; pnode = SUCC(pnode)){
      NDS(pnode,y) = 0.;
      nfstop = STOP_NF(pnode,t_column);
      for (pnf = NF_START(pnode,t_column); pnf != nfstop; pnf=NEXT(pnf))
         NDS(pnode,y) += DOT(COEFFBDP(pnf,Z,0),FDDVP(NBFACE(pnf),x,0)) +
                         DOT(COEFFBDP(pnf,Z,1),FDDVP(NBFACE(pnf),x,1));
   }
}

void multA_dvf_X_sn(tGrid,Z,x,y) /*  y := Z x  */
GRID *tGrid;                     /*  rows ... dvectors in faces               */
INT Z, x, y;                     /*  columns ... scalars in nodes             */
{
   NODE *pnode;
   FACE *pf;
   NFLINK *pnf;

   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      D_SET_VALUE(pface,0.,y)
   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ)
      for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next){
         pf = NBFACE(pnf);
         SET4(FDDVP(pf,y,0),COEFFBDP(pnf,Z,0),NDS(pnode,x))
         SET4(FDDVP(pf,y,1),COEFFBDP(pnf,Z,1),NDS(pnode,x))
      }
}

#else

void multA_sn_X_dvf(tGrid,Z,x,y,t_row,t_column)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t_row, t_column;
{  eprintf("Error: multA_sn_X_dvf not available.\n");  }

void multA_dvf_X_sn(tGrid,Z,x,y) /*  y := Z x  */
GRID *tGrid; INT Z, x, y;
{  eprintf("Error: multA_dvf_X_sn not available.\n");  }

#endif

#if (F_DATA & Ix2D_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void multA_sf_X_dvf_zero_diag(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... scalars in faces                              */
INT Z, x, y, t;    /*  columns ... dvectors in faces                          */
{
   FACE *pface;
   FLINK *pfl, *fstop;

   for (pface = FIRSTF(tGrid); pface != NULL; pface = SUCC(pface)){
/*      FD(pface,y) = DOT(COEFFBDP(pface,Z,0),FDDVP(pface,x,0));  */
      FD(pface,y) = 0.;
      fstop = STOP_FF(pface,t);
      for (pfl = F_START(pface,t); pfl != fstop; pfl = NEXT(pfl))
         FD(pface,y) += DOT(COEFFBDP(pfl,Z,0),FDDVP(NBFACE(pfl),x,0)) +
                        DOT(COEFFBDP(pfl,Z,1),FDDVP(NBFACE(pfl),x,1));
   }
}

void multA_dvf_X_sf_zero_diag(tGrid,Z,x,y)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... dvectors in faces                             */
INT Z, x, y;       /*  columns ... scalars in faces                           */
{
   FACE *pface;
   FLINK *pfl;

   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      D_SET_VALUE(pface,0.,y)
   for (pface = FIRSTF(tGrid); pface != NULL; pface = SUCC(pface))
      for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl)){
         SET4(FDDVP(NBFACE(pfl),y,0),COEFFBDP(pfl,Z,0),FD(pface,x))
         SET4(FDDVP(NBFACE(pfl),y,1),COEFFBDP(pfl,Z,1),FD(pface,x))
      }
/*
      SET2(FDDVP(pface,y,0),COEFFBDP(pface,Z,0),FD(pface,x))
      SET2(FDDVP(pface,y,1),COEFFBDP(pface,Z,1),FD(pface,x))
*/
}

#else  /*  if !((F_DATA & Ix2D_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (F_DATA & SCALAR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void multA_sf_X_dvf_zero_diag(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sf_X_dvf_zero_diag not available.\n");  }

void multA_dvf_X_sf_zero_diag(tGrid,Z,x,y)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y;
{  eprintf("Error: multA_dvf_X_sf_zero_diag not available.\n");  }

#endif

#if (E_DATA & ExFx2DF_MATR) && (DIM == 2)

void set_mat_value_sne_X_dvf(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < NVERT; i++)
         for (j = 0; j < SIDES; j++){
             SET7(COEFF_BDFSP(pel,Z,i,j,0),r)
             SET7(COEFF_BDFSP(pel,Z,i,j,1),r)
         }
}

#else

void set_mat_value_sne_X_dvf(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_sne_X_dvf not available.\n");  }

#endif

#if (E_DATA & ExFx2DF_MATR) && (F_DATA & DVECTOR_FACE_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)

void multA_sef_X_dvf(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... scalars in element faces                      */
INT Z, x, y, t;    /*  columns ... dvectors in faces                          */
{
   ELEMENT *pel;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < NVERT; i++){
         EDSN(pel,y,i) = 0.;
         for (j = 0; j < SIDES; j++)
            if (IS_CF(pel->f[j],t))
             EDSN(pel,y,i) += DOT(COEFF_BDFSP(pel,Z,i,j,0),FDDVP(pel->f[j],x,0)) +
                              DOT(COEFF_BDFSP(pel,Z,i,j,1),FDDVP(pel->f[j],x,1));
       }
}

void multA_dvf_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid;       /*  rows ... dvectors in faces                             */
INT Z, x, y, t;    /*  columns ... scalars in element faces                   */
{
   ELEMENT *pel;
   INT i, j;

   dvset_value_f(tGrid,0.,y,1);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < NVERT; i++)
         for (j = 0; j < SIDES; j++){
            SET4(FDDVP(pel->f[j],y,0),COEFF_BDFSP(pel,Z,i,j,0),EDSN(pel,x,i))
            SET4(FDDVP(pel->f[j],y,1),COEFF_BDFSP(pel,Z,i,j,1),EDSN(pel,x,i))
         }
}

#else  /*  if !((E_DATA & ExFx2DF_MATR) && (F_DATA & DVECTOR_FACE_DATA) && 
                (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES))                      */

void multA_sef_X_dvf(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_sef_X_dvf not available.\n");  }

void multA_dvf_X_sef(tGrid,Z,x,y,t)  /*  y := Z x  */
GRID *tGrid; INT Z, x, y, t;
{  eprintf("Error: multA_dvf_X_sef not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (F_DATA & SCALAR_FACE_DATA)

void copy_sdiag_to_vector_sf(theGrid,a,u)
GRID *theGrid;
INT a, u;
{
   FACE *pface;
 
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      FD(pface,u) = COEFF_FF(pface,a);
}

#else

void copy_sdiag_to_vector_sf(theGrid,a,u)
GRID *theGrid; INT a, u;
{  eprintf("Error: copy_sdiag_to_vector_sf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA)

void copy_sdiag_to_vector_vf(theGrid,a,u)
GRID *theGrid;
INT a, u;
{
   FACE *pface;
 
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      SET7(FDVP(pface,u),COEFF_FF(pface,a))
}

#else

void copy_sdiag_to_vector_vf(theGrid,a,u)
GRID *theGrid; INT a, u;
{  eprintf("Error: copy_sdiag_to_vector_vf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (F_DATA & VECTOR_FACE_DATA) && (DIM == 2)

void copy_vdiag_to_vector_vf(theGrid,a,u)
GRID *theGrid;
INT a, u;
{  
   FACE *pface;
      
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      FDV(pface,u,0) = COEFFNN(pface,a,0,0);
      FDV(pface,u,1) = COEFFNN(pface,a,1,1);
   }
}

#else

void copy_vdiag_to_vector_vf(theGrid,a,u)
GRID *theGrid; INT a, u;
{  eprintf("Error: copy_vdiag_to_vector_vf not available.\n");  }

#endif

#if (N_DATA & NxN_NODE_MATR) && (N_DATA & NxM_NODE_FACE_MATR) && (F_DATA & NxN_FACE_MATR) && (F_DATA & MxN_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & MxM_E_E_MATR) && (E_DATA & MxN_E_N_MATR) && (E_DATA & NxM_N_E_MATR) && (E_DATA & MxN_E_F_MATR) && (E_DATA & NxM_F_E_MATR)

void set_mat_value_general(tGrid,Z,r)
GRID *tGrid;
INT Z;
FLOAT r;
{
   NODE *theNode;
   FACE *theFace;
   ELEMENT *pel;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   INT i, j, k, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC, mm=N_OF_ELEM_FUNC;

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      MMSET7(COEFF_NNP(theNode,Z),r,i,j,kk,kk)
      for (pl=TSTART(theNode); pl; pl=NEXT(pl))
         MMSET7(COEFF_NNP(pl,Z),r,i,j,kk,kk)
      for (pnfl=TNFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
         MMSET7(COEFF_NNP(pnfl,Z),r,i,j,kk,nn)
   }
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      MMSET7(COEFF_NNP(theFace,Z),r,i,j,nn,nn)
      for (pfl=TFSTART(theFace); pfl; pfl=NEXT(pfl))
         MMSET7(COEFF_NNP(pfl,Z),r,i,j,nn,nn)
      for (pfnl=TFNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
         MMSET7(COEFF_NNP(pfnl,Z),r,i,j,nn,kk)
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      MMSET7(COEFF_EE_MMP(pel,Z),r,i,j,mm,mm)
      for (k=0; k < NVERT; k++){
         MMSET7(COEFF_NE_NMP(pel,Z,k),r,i,j,kk,mm)
         MMSET7(COEFF_FE_NMP(pel,Z,k),r,i,j,nn,mm)
         MMSET7(COEFF_EF_MNP(pel,Z,k),r,i,j,mm,nn)
         MMSET7(COEFF_EN_MNP(pel,Z,k),r,i,j,mm,kk)
      }
   }
}

void set_mat_to_zero_for_q2b3(tGrid,Z)
GRID *tGrid;
INT Z;
{
   NODE *theNode;
   FACE *theFace;
   ELEMENT *pel;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   DOUBLE q;
   INT i, j, k, mm=N_OF_ELEM_FUNC;

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      for (pnfl=TNFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
         COEFF_NNP(pnfl,Z)[0][1] = 0.;
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      COEFF_NNP(theFace,Z)[0][1] = COEFF_NNP(theFace,Z)[1][0] = 0.;
      COEFF_NNP(theFace,Z)[1][1] = 1.;
      for (pfl=TFSTART(theFace); pfl; pfl=NEXT(pfl))
         COEFF_NNP(pfl,Z)[0][1] = COEFF_NNP(pfl,Z)[1][0] 
                                = COEFF_NNP(pfl,Z)[1][1] = 0.;
      for (pfnl=TFNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
         COEFF_NNP(pfnl,Z)[1][0] = 0.;
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      q = COEFF_EE_MMP(pel,Z)[0][0];
      for(i=0; i < mm; i++)
      for(j=0; j < mm; j++) COEFF_EE_MMP(pel,Z)[i][j] = 0.;
      COEFF_EE_MMP(pel,Z)[0][0] = q;
      for(i=1; i < mm; i++) COEFF_EE_MMP(pel,Z)[i][i] = 1.;
      for (k=0; k < NVERT; k++){
         for(i=1; i < mm; i++){
            COEFF_EN_MNP(pel,Z,k)[i][0] = 0.;
            COEFF_NE_NMP(pel,Z,k)[0][i] = 0.;
            COEFF_EF_MNP(pel,Z,k)[i][0] = 0.;
            COEFF_FE_NMP(pel,Z,k)[0][i] = 0.;
         }
         for(i=0; i < mm; i++)
            COEFF_EF_MNP(pel,Z,k)[i][1] = COEFF_FE_NMP(pel,Z,k)[1][i] = 0.;
      }
   }
}

void partial_set_mat_to_zero_for_q2b3(tGrid,Z)
GRID *tGrid;
INT Z;
{
   NODE *n1, *n2;
   FACE *theFace;
   ELEMENT *pel;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   DOUBLE q;
   INT i, j, k, mm=N_OF_ELEM_FUNC;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
   if (IS_B_EL(pel)){
      for (k=0; k < NVERT; k++){
         theFace = pel->f[k];
         n1 = pel->n[k];
         n2 = pel->n[(k+1)%NVERT];
         if (IS_BF(theFace) || (IS_FN(n1) && IS_FN(n2))){
            for(i=0; i < mm; i++)
               COEFF_EF_MNP(pel,Z,k)[i][1] = COEFF_FE_NMP(pel,Z,k)[1][i] = 0.;
            for (j=0; j < NVERT; j++){
               for (pnfl=TNFSTART(pel->n[j]); pnfl->nbface != theFace; 
                                                               pnfl=NEXT(pnfl));
               COEFF_NNP(pnfl,Z)[0][1] = 0.;
               if (j != k){
                  for (pfl=TFSTART(pel->f[j]); pfl->nbface != theFace; 
                                                                pfl=NEXT(pfl));
                  COEFF_NNP(pfl,Z)[0][1] = COEFF_NNP(pfl,Z)[1][1] = 0.;
               }
            }
            for (pfnl=TFNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
               COEFF_NNP(pfnl,Z)[1][0] = 0.;
            COEFF_NNP(theFace,Z)[0][1] = COEFF_NNP(theFace,Z)[1][0] = 0.;
            COEFF_NNP(theFace,Z)[1][1] = 1.;
            for (pfl=TFSTART(theFace); pfl; pfl=NEXT(pfl))
               COEFF_NNP(pfl,Z)[1][0] = COEFF_NNP(pfl,Z)[1][1] = 0.;
         }
      }
   }
   else
   {
      q = COEFF_EE_MMP(pel,Z)[0][0];
      for(i=0; i < mm; i++)
      for(j=0; j < mm; j++) COEFF_EE_MMP(pel,Z)[i][j] = 0.;
      COEFF_EE_MMP(pel,Z)[0][0] = q;
      for(i=1; i < mm; i++) COEFF_EE_MMP(pel,Z)[i][i] = 1.;
      for (k=0; k < NVERT; k++){
         for(i=1; i < mm; i++){
            COEFF_EN_MNP(pel,Z,k)[i][0] = 0.;
            COEFF_NE_NMP(pel,Z,k)[0][i] = 0.;
            COEFF_EF_MNP(pel,Z,k)[i][0] = 0.;
            COEFF_FE_NMP(pel,Z,k)[0][i] = 0.;
         }
         for(i=0; i < mm; i++)
            COEFF_EF_MNP(pel,Z,k)[i][1] = COEFF_FE_NMP(pel,Z,k)[1][i] = 0.;
         for (pnfl=TNFSTART(pel->n[k]); pnfl; pnfl=NEXT(pnfl))
            if (pnfl->nbface == pel->f[0] || pnfl->nbface == pel->f[1] ||
                pnfl->nbface == pel->f[2] || pnfl->nbface == pel->f[3])
               COEFF_NNP(pnfl,Z)[0][1] = 0.;
         theFace = pel->f[k];
         COEFF_NNP(theFace,Z)[0][1] = COEFF_NNP(theFace,Z)[1][0] = 0.;
         COEFF_NNP(theFace,Z)[1][1] = 1.;
         for (pfl=TFSTART(theFace); pfl; pfl=NEXT(pfl))
            if (pfl->nbface == pel->f[0] || pfl->nbface == pel->f[1] ||
                pfl->nbface == pel->f[2] || pfl->nbface == pel->f[3])
               COEFF_NNP(pfl,Z)[0][1] = COEFF_NNP(pfl,Z)[1][0] 
                                      = COEFF_NNP(pfl,Z)[1][1] = 0.;
         for (pfnl=TFNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
            COEFF_NNP(pfnl,Z)[1][0] = 0.;
      }
   }
}

#else

void set_mat_value_general(tGrid,Z,r)
GRID *tGrid; INT Z; FLOAT r;
{  eprintf("Error: set_mat_value_general not available.\n");  }

void set_mat_to_zero_for_q2b3(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: set_mat_to_zero_for_q2b3 not available.\n");  }

void partial_set_mat_to_zero_for_q2b3(tGrid,Z)
GRID *tGrid; INT Z;
{  eprintf("Error: partial_set_mat_to_zero_for_q2b3 not available.\n");  }

#endif

#if (N_DATA & NxN_NODE_MATR) && (N_DATA & NxM_NODE_FACE_MATR) && (F_DATA & NxN_FACE_MATR) && (F_DATA & MxN_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & MxM_E_E_MATR) && (E_DATA & MxN_E_N_MATR) && (E_DATA & NxM_N_E_MATR) && (E_DATA & MxN_E_F_MATR) && (E_DATA & NxM_F_E_MATR) && (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)

void multA_general(tGrid,Z,x,y,q,t) /*  y := Z x  */
GRID *tGrid;
INT Z, x, y, q, t;
{
   NODE *theNode;
   FACE *theFace;
   ELEMENT *pel;
   LINK *pl, *stop;
   NFLINK *pnfl, *nfstop;
   FLINK *pfl, *fstop;
   FNLINK *pfnl, *fnstop;
   INT i, j, k, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC, mm=N_OF_ELEM_FUNC;

   if (!(t & STOP_IS_FIRST_INNER)){
      for (theNode=FDBN(tGrid);theNode!=FIRSTNODE(tGrid);theNode=SUCC(theNode)){
         GSET1(NDMVP(theNode,q),NDMVP(theNode,x),i,kk)
         GSET7(NDMVP(theNode,x),0.,i,kk)
      }
      for (theFace=FDBF(tGrid);theFace!=FIRSTFACE(tGrid);theFace=SUCC(theFace)){
         GSET1(FDMVP(theFace,q),FDMVP(theFace,x),i,nn)
         GSET7(FDMVP(theFace,x),0.,i,nn)
      }
   }
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      MMSET2(NDMVP(theNode,y),COEFF_NNP(theNode,Z),NDMVP(theNode,x),i,j,kk,kk)
      stop = STOP_NN(theNode,t);
      for (pl=T_START(theNode,t); pl!=stop; pl=NEXT(pl))
         MMSET4(NDMVP(theNode,y),COEFF_NNP(pl,Z),NDMVP(NBNODE(pl),x),i,j,kk,kk)
      nfstop = STOP_NF(theNode,t);
      for (pnfl=NF_START(theNode,t); pnfl != nfstop; pnfl=NEXT(pnfl))
         MMSET4(NDMVP(theNode,y),COEFF_NNP(pnfl,Z),FDMVP(NBFACE(pnfl),x),i,j,kk,nn)
   }
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      MMSET2(FDMVP(theFace,y),COEFF_NNP(theFace,Z),FDMVP(theFace,x),i,j,nn,nn)
      fstop = STOP_FF(theFace,t);
      for (pfl=F_START(theFace,t); pfl != fstop; pfl=NEXT(pfl))
         MMSET4(FDMVP(theFace,y),COEFF_NNP(pfl,Z),FDMVP(NBFACE(pfl),x),i,j,nn,nn)
      fnstop = STOP_FN(theFace,t);
      for (pfnl=FN_START(theFace,t); pfnl != fnstop; pfnl=NEXT(pfnl))
         MMSET4(FDMVP(theFace,y),COEFF_NNP(pfnl,Z),NDMVP(NBNODE(pfnl),x),i,j,nn,kk)
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      MMSET2(EDMVP(pel,y),COEFF_EE_MMP(pel,Z),EDMVP(pel,x),i,j,mm,mm)
      for (k=0; k < NVERT; k++){
         MMSET4(NDMVP(pel->n[k],y),COEFF_NE_NMP(pel,Z,k),EDMVP(pel,x),i,j,kk,mm)
         MMSET4(FDMVP(pel->f[k],y),COEFF_FE_NMP(pel,Z,k),EDMVP(pel,x),i,j,nn,mm)
         MMSET4(EDMVP(pel,y),COEFF_EF_MNP(pel,Z,k),FDMVP(pel->f[k],x),i,j,mm,nn)
         MMSET4(EDMVP(pel,y),COEFF_EN_MNP(pel,Z,k),NDMVP(pel->n[k],x),i,j,mm,kk)
      }
   }
   if (!(t & STOP_IS_FIRST_INNER)){
      for (theNode=FDBN(tGrid);theNode!=FIRSTNODE(tGrid);theNode=SUCC(theNode))
         GSET1(NDMVP(theNode,x),NDMVP(theNode,q),i,kk)
      for (theFace=FDBF(tGrid);theFace!=FIRSTFACE(tGrid);theFace=SUCC(theFace))
         GSET1(FDMVP(theFace,x),FDMVP(theFace,q),i,nn)
   }
}

#else

void multA_general(tGrid,Z,x,y,q,t)
GRID *tGrid; INT Z, x, y, q, t;
{  eprintf("Error: multA_general not available.\n");  }

#endif

void mult_A(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;                                                      /* y := Z x */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
INT p0,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT r1,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   if (PERIODIC_BC == YES)
      copy_1st_to_2nd_periodic_boundary(tGrid,x,row_type);
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_SN: if (structure & Q_TRANS) 
                         smultAT(tGrid,Z,x,y,t_row);
                      else
                         smultA(tGrid,Z,x,y,t_row);
                break;
           case Q_VN: multA_sn_X_vn(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVF: multA_sn_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVE: multA_sn_X_vn(tGrid,Z,x,y,t_row);
                        multA_sn_X_ve(tGrid,Z,x,y,t_row);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                          multA_sn_X_dvf0(tGrid,Z,x,y,t_row);
                       else
                          multA_sn_X_dvf(tGrid,Z,x,y,t_row,t_column);
                break;
           case Q_VE: set_value(tGrid,0.,y,t_row,row_type);
                      multA_sn_X_ve(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VN: 
        switch(column_type){
           case Q_SN: multA_vn_X_sn(tGrid,Z,x,y,t_row);
                break;
           case Q_VN: if (structure & Q_FULL)
                         vmultA_vn(tGrid,Z,x,y,t_row);
                      else if (structure & Q_NEBDIAG)
                         smultA_vn(tGrid,Z,x,y,t_row);
                      else
                         eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: multA_vn_X_e(tGrid,Z,x,y,q,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_SF: 
        switch(column_type){
           case Q_SF: multA_f(tGrid,Z,x,y,t_row);
                break;
           case Q_DVF: if (structure & Q_ZERO_DIAG)
                          multA_sf_X_dvf_zero_diag(tGrid,Z,x,y,t_row);
                       else
                          eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: multA_sf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         vmultA_vf(tGrid,Z,x,y,t_row);
                      else if (structure & Q_FEBDIAG)
                         smultA_vf(tGrid,Z,x,y,t_row);
                      else
                         eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: multA_vf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vmultA_dvf(tGrid,Z,x,y,t_row);
                       else
                          eprintf("Error: mult_A not available.\n");
                break;
           case Q_SN: if (structure & Q_FIRST_DV)
                         multA_dvf0_X_sn(tGrid,Z,x,y,t_row);
                      else
                         multA_dvf_X_sn(tGrid,Z,x,y);
                break;
           case Q_SF: if (structure & Q_ZERO_DIAG)
                         multA_dvf_X_sf_zero_diag(tGrid,Z,x,y);
                      else
                         eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: if (structure & Q_FIRST_DV)
                         multA_dvf_X_e(tGrid,Z,x,y,q,t_row);
                      else
                         eprintf("Error: mult_A not available.\n");
                break;
           case Q_SNE: multA_dvf_X_sef(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_SE: 
           switch(column_type){
           case Q_VN: multA_e_X_vn(tGrid,Z,x,y,t_row);
                break;
           case Q_SF: multA_e_X_sf(tGrid,Z,x,y,t_row);
                break;
           case Q_VF: multA_e_X_vf(tGrid,Z,x,y,t_row);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                         multA_e_X_dvf(tGrid,Z,x,y,t_row);
                      else
                         eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: multA_se_X_se(tGrid,Z,x,y);
                break;
           case Q_VNSF: multA_e_X_vn_sf(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVF: multA_e_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNSFLG: multB_lg(tGrid,Z,x,y);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
           }
        break;
   case Q_SNE: 
           switch(column_type){
           case Q_DVF: multA_sef_X_dvf(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVF: multA_sef_X_vn_vf(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVFVE: multA_sef_X_vn_vf_ve(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
           }
        break;
   case Q_VE: 
           switch(column_type){
           case Q_SN: multA_ve_X_sn(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
           }
        break;
   case Q_SNSF: 
        switch(column_type){
           case Q_SNSF: smultA_sn_sf(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_SNSE: 
        switch(column_type){
           case Q_SNSE: smultA_sn_se(tGrid,Z,x,y,q,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VNSF: 
        switch(column_type){
           case Q_VNSF: if (structure & Q_FULL)
                           vmultA_vn_sf(tGrid,Z,x,y);
                        else if ((structure & Q_NEBDIAG) && 
                                 (structure & Q_TRANS))
                           smultAT_vn_sf(tGrid,Z,x,y);
                        else if (structure & Q_NEBDIAG)
                           smultA_vn_sf(tGrid,Z,x,y,t_row);
                        else if (structure & Q_NEBDIAG_NRED)
                           smultA_vn_sf_nred(tGrid,Z,x,y,t_row);
                        else if (structure & Q_NEBDIAG_NFRED)
                           smultA_vn_sf_nfred(tGrid,Z,x,y,t_row);
                        else
                           eprintf("Error: mult_A not available.\n");
                break;
           case Q_SE: multA_vn_sf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VNVF: 
        switch(column_type){
           case Q_VNVF: if (structure & Q_FULL)
                           vmultA_vn_vf(tGrid,Z,x,y,t_row);
                        else if (structure & Q_EBDIAG)
                           smultA_vn_vf(tGrid,Z,x,y,t_row);
                        else
                           eprintf("Error: mult_A not available.\n");
                break;
           case Q_SN: multA_vn_vf_X_sn(tGrid,Z,x,y);
                break;
           case Q_SE: multA_vn_vf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_SNE: multA_vn_vf_X_sef(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VNVE:
        switch(column_type){
           case Q_VNVE: if (structure & Q_NEBDIAG_EDIAG){
                           smultA_vn(tGrid,Z,x,y,t_row);
                           if (!(t_row & STOP_IS_FIRST_INNER))
                              vmultiply_1e(tGrid,x,Z,y);
                        }
                        else
                           eprintf("Error: mult_A not available.\n");
                break;
           case Q_SN: multA_vn_X_sn(tGrid,Z,x,y,t_row);
                      multA_ve_X_sn(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VNVFVE: 
        switch(column_type){
           case Q_VNVFVE: if (structure & Q_EBDIAG)
                             smultA_vn_vf_ve(tGrid,Z,x,y,q,t_row);
                          else
                             eprintf("Error: mult_A not available.\n");
                break;
           case Q_SNE: multA_vn_vf_ve_X_sef(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_SSNSSF: 
        switch(column_type){
           case Q_VNVF: multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVFVE: multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_VNSFLG: 
        switch(column_type){
           case Q_VNSFLG: if (structure & Q_FULL)
                              vmultA_vn_sf_lg(tGrid,Z,x,y);
                           else
                              eprintf("Error: mult_A not available.\n");
                break;
           default:
                eprintf("Error: mult_A not available.\n");
                break;
        }
        break;
   case Q_GENERAL:
        if (column_type == Q_GENERAL)
           multA_general(tGrid,Z,x,y,q,t_row);
        else
           eprintf("Error: mult_A not available.\n");
        break;
   default:
        eprintf("Error: mult_A not available.\n");
        break;
   }
   if (PERIODIC_BC == YES){
      add_2nd_to_1st_periodic_boundary(tGrid,y,column_type);
      set_zero_on_2nd_periodic_boundary(tGrid,y,column_type);
      set_zero_on_2nd_periodic_boundary(tGrid,x,row_type);
   }
}

void subtract_local_proj();

void mult_A_lp(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
               p0,p1,p2,p3,p4,p5,p6,p7,p8,space,r1,r2,r3)
GRID *tGrid;                                                     /* y := Z x */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
INT p0, p1, p2, p3, p4, p5, p6, p7, p8, space;
FLOAT r1, r2, r3;
{
   mult_A(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
          p0,p1,p2,p3,p4,p5,p6,p7,p8,0,r1,r2,r3);
   copy(tGrid,x,q,STOP_IS_FIRST_INNER,column_type);
   set_value(tGrid,0.,x,STOP_IS_FIRST_INNER,column_type);
   subtract_local_proj(tGrid,TNU,bb0,bb1,x,y,space);
   copy(tGrid,q,x,STOP_IS_FIRST_INNER,column_type);
}

void mult_AL(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;                       /* y := Z_lower_triangle x  (without diag) */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
INT p0,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT r1,r2,r3;
{
   set_value(tGrid,1.,q,t_row,row_type);
   switch(row_type){
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         eprintf("Error: mult_AL not available.\n");
                      else if (structure & Q_FEBDIAG)
                         smultAL_vf(tGrid,Z,x,y,q);
                      else
                         eprintf("Error: mult_AL not available.\n");
                break;
           default:
                eprintf("Error: mult_AL not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vmultAL_dvf(tGrid,Z,x,y,q);
                       else
                          eprintf("Error: mult_AL not available.\n");
                break;
           default:
                eprintf("Error: mult_AL not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_AL not available.\n");
        break;
   }
}

void mult_AU(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;                       /* y := Z_upper_triangle x  (without diag) */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
INT p0,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT r1,r2,r3;
{
   set_value(tGrid,1.,q,t_row,row_type);
   switch(row_type){
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         eprintf("Error: mult_AU not available.\n");
                      else if (structure & Q_FEBDIAG)
                         smultAU_vf(tGrid,Z,x,y,q);
                      else
                         eprintf("Error: mult_AU not available.\n");
                break;
           default:
                eprintf("Error: mult_AU not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vmultAU_dvf(tGrid,Z,x,y,q);
                       else
                          eprintf("Error: mult_AU not available.\n");
                break;
           default:
                eprintf("Error: mult_AU not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_AU not available.\n");
        break;
   }
}

void copy_diag_to_vector(tGrid,Z,x,t,type,structure)  /*  x := diag Z  */
GRID *tGrid;
INT Z, x, t, type, structure;
{
   switch(type){
   case Q_SF: copy_sdiag_to_vector_sf(tGrid,Z,x);
        break;
   case Q_VF: if (structure & Q_FULL)
                 copy_vdiag_to_vector_vf(tGrid,Z,x);
              else if (structure & Q_FEBDIAG)
                 copy_sdiag_to_vector_vf(tGrid,Z,x);
              else
                 eprintf("Error: copy_diag_to_vector not available.\n");
        break;
   default:
        eprintf("Error: copy_diag_to_vector not available.\n");
        break;
   }
}

void mult_diag(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure)
GRID *tGrid;                                               /* y := (diag Z) x */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
{
   switch(row_type){
   case Q_VF:
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         vmult_diagA_vf(tGrid,Z,x,y);
                      else if (structure & Q_FEBDIAG)
                         smult_diagA_vf(tGrid,Z,x,y);
                      else
                         eprintf("Error: mult_diag not available.\n");
                break;
           default:
                eprintf("Error: mult_diag not available.\n");
                break;
        }
        break;
   case Q_DVF:
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vmult_diagA_dvf(tGrid,Z,x,y);
                       else
                          eprintf("Error: mult_diag not available.\n");
                break;
           default:
                eprintf("Error: mult_diag not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_diag not available.\n");
        break;
   }
}

void mult_inv_diag(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure)
GRID *tGrid;                                          /* y := (diag Z)^(-1) x */INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
{
   switch(row_type){
   case Q_VF:
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         vmult_inv_diagA_vf(tGrid,Z,x,y);
                      else if (structure & Q_FEBDIAG)
                         smult_inv_diagA_vf(tGrid,Z,x,y);
                      else
                         eprintf("Error: mult_inv_diag not available.\n");
                break;
           default:
                eprintf("Error: mult_inv_diag not available.\n");
                break;
        }
        break;
   case Q_DVF:
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vmult_inv_diagA_dvf(tGrid,Z,x,y);
                       else
                          eprintf("Error: mult_inv_diag not available.\n");
                break;
           default:
                eprintf("Error: mult_inv_diag not available.\n");
                break;
        }
        break;
   case Q_VNVF: 
        switch(column_type){
           case Q_VNVF: if (structure & Q_FULL)
                           eprintf("Error: mult_inv_diag not available.\n");
                        else if (structure & Q_EBDIAG)
                           smult_inv_diag_vn_vf(tGrid,Z,x,y);
                        else
                           eprintf("Error: mult_inv_diag not available.\n");
                break;
           default:
                eprintf("Error: mult_inv_diag not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_inv_diag not available.\n");
        break;
   }
}

FLOAT norm_of_A(tGrid,Z,t_row,t_column,row_type,column_type,structure)
GRID *tGrid;                                        /*  L infinity norm of Z  */
INT Z, t_row, t_column, row_type, column_type, structure;
{
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_SN: return(snorm_of_A(tGrid,Z,t_row));
                break;
           default:
                eprintf("Error: norm_of_A not available.\n");
                return(0.);
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG) 
                          return(vnorm_of_A_dvf(tGrid,Z));
                       else{
                          eprintf("Error: norm_of_A not available.\n");
                          return(0.);
                       }
                break;
           default:
                eprintf("Error: norm_of_A not available.\n");
                return(0.);
                break;
        }
        break;
   case Q_VNSF: 
        switch(column_type){
           case Q_VNSF: if (structure & Q_FULL)
                           return(vnorm_of_A_vn_sf(tGrid,Z));
                        else if (structure & Q_NEBDIAG)
                           return(norm_of_A_vn_sf(tGrid,Z));
                        else{
                           eprintf("Error: norm_of_A not available.\n");
                           return(0.);
                        }
                break;
           default:
                eprintf("Error: norm_of_A not available.\n");
                return(0.);
                break;
        }
        break;
   default:
        eprintf("Error: norm_of_A not available.\n");
        return(0.);
        break;
   }
}

void set_mat_value(tGrid,Z,r,t_row,t_column,row_type,column_type,structure)
GRID *tGrid;
INT Z, t_row, t_column, row_type, column_type, structure;
FLOAT r;
{
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_SN: sset_mat_value(tGrid,Z,r);
                break;
           case Q_VN: set_mat_value_snXvn(tGrid,Z,r);
                break;
           case Q_VNVF: if (structure & Q_BDAUGMENT)
                            setA_sn_X_vn_vf_BD(tGrid,Z,r);
                        else
                            setA_sn_X_vn_vf(tGrid,Z,r);
                break;
           case Q_VNVE: set_mat_value_sn_X_vn_ve(tGrid,Z,r);
                break;
           case Q_VNVFVE: set_mat_value_sn_X_vn_vf_ve(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VN:
        switch(column_type){
           case Q_VN: if (structure & Q_FULL)
                         vset_mat_value_vn(tGrid,Z,r);
                      else if (structure & Q_NEBDIAG)
                         sset_mat_value(tGrid,Z,r);
                      else
                         eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SF: 
        switch(column_type){
           case Q_SF: set_mat_value_f(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         vset_mat_value_vf(tGrid,Z,r);
                      else if (structure & Q_FEBDIAG)
                         set_mat_value_f(tGrid,Z,r);
                      else
                         eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          set_mat_value_dvf(tGrid,Z,r);
                       else
                          eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SE:
        switch(column_type){
           case Q_VN: set_mat_value_e_X_vn(tGrid,Z,r);
                break;
           case Q_VF: set_mat_value_e_X_vf(tGrid,Z,r);
                break;
           case Q_DVF: set_mat_value_e_X_vf(tGrid,Z,r);
                break;
           case Q_VNSF: set_mat_value_e_X_vn_sf(tGrid,Z,r);
                break;
           case Q_VNVF: set_mat_value_e_X_vn_vf(tGrid,Z,r);
                break;
           case Q_VNVFVE: if (structure & Q_NE_ZERO)
                        set_mat_value_e_X_vn_vf(tGrid,Z,r);
                       else
                        eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SNSF:
        switch(column_type){
           case Q_SNSF: sset_mat_value_sn_sf(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SNSE:
        switch(column_type){
           case Q_SNSE: sset_mat_value_sn_se(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SNE: 
           switch(column_type){
           case Q_DVF: set_mat_value_sne_X_dvf(tGrid,Z,r);
                break;
           case Q_VNVF: set_mat_value_sne_X_vn_vf(tGrid,Z,r);
                break;
           case Q_VNVFVE: set_mat_value_sne_X_vn_vf_ve(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
           }
        break;
   case Q_VNSF: 
        switch(column_type){
           case Q_VNSF: if (structure & Q_FULL)
                           vset_mat_value_vn_sf(tGrid,Z,r);
                        else if (structure & Q_NEBDIAG)
                           sset_mat_value_vn_sf(tGrid,Z,r);
                        else if (structure & Q_NEBDIAG_NRED)
                           sset_mat_value_sn_sf_nred(tGrid,Z,r);
                        else if (structure & Q_NEBDIAG_NFRED)
                           sset_mat_value_sn_sf_nfred(tGrid,Z,r);
                        else
                           eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VNVF:
        switch(column_type){
           case Q_VNVF: if (structure & Q_FULL)
                           vset_mat_value_vn_vf(tGrid,Z,r);
                        else if (structure & Q_EBDIAG)
                           sset_mat_value_sn_sf(tGrid,Z,r);
                        else
                           eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VNVE: 
        switch(column_type){
           case Q_VNVE: if (structure & Q_NEBDIAG_EDIAG){
                           sset_mat_value(tGrid,Z,r);
                           set_value_e(tGrid,r,Z);
                        }
                        else
                           eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VNVFVE:
        switch(column_type){
           case Q_VNVFVE: if (structure & Q_EBDIAG)
                             sset_mat_value_sn_sf_se(tGrid,Z,r);
                          else
                             eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_VNSFLG: 
        switch(column_type){
           case Q_VNSFLG: if (structure & Q_FULL)
                             vset_mat_value_vn_sf_lg(tGrid,Z,r);
                          else if (structure & Q_NEBDIAG)
                             sset_mat_value_vn_sf_lg(tGrid,Z,r);
                          else
                             eprintf("Error: set_mat_value not available.\n");
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_SNSSNSSF: 
        switch(column_type){
           case Q_VNVF: setA_sn_X_vn_vf(tGrid,Z,r);
                        set_mat_value_ssn_ssf_X_vn_vf(tGrid,Z,r);
                break;
           case Q_VNVFVE: set_mat_value_sn_X_vn_vf_ve(tGrid,Z,r);
                          set_mat_value_ssn_ssf_X_vn_vf(tGrid,Z,r);
                break;
           default:
                eprintf("Error: set_mat_value not available.\n");
                break;
        }
        break;
   case Q_GENERAL:
        if (column_type == Q_GENERAL)
           set_mat_value_general(tGrid,Z,r);
        else
           eprintf("Error: set_mat_value not available.\n");
        break;
   default:
        eprintf("Error: set_mat_value not available.\n");
        break;
   }
}

void make_AT(tGrid,Z,a,row_type,column_type,structure)
GRID *tGrid;                                                      /* a := Z^T */
INT Z, a, row_type, column_type, structure;
{
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_SN: smakeAT(tGrid,Z,a,0);
                break;
           default:
                eprintf("Error: make_AT not available.\n");
                break;
        }
        break;
   case Q_SNSF: 
        switch(column_type){
           case Q_SNSF: smakeAT_sn_sf(tGrid,Z,a);
                break;
           default:
                eprintf("Error: make_AT not available.\n");
                break;
        }
        break;
   case Q_VNSF: 
        switch(column_type){
           case Q_VNSF: if (structure & Q_NEBDIAG)
                           smakeAT_vn_sf(tGrid,Z,a);
                        else
                           eprintf("Error: make_AT not available.\n");
                break;
           default:
                eprintf("Error: make_AT not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: make_AT not available.\n");
        break;
   }
}

void SOR_step_forward(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,
                      structure,b,p1,p2,p3,p4,p5,p6,p7,p8,p9,om,r2,r3)
GRID *tGrid;       /*  one forward step of SOR for the solution of Zx=b       */
INT Z, x, y, q;    /*  x ... initial guess, y ... next iterate, y=x possible  */
INT t_row, t_column, row_type, column_type, structure;      /*  row = column  */
INT b,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT om,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   if (x != y)
      copy(tGrid,x,y,t_row,row_type);
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_SN: sSOR_step_forward_sn(tGrid,Z,b,x,om);
                break;
           default:
                eprintf("Error: SOR_step_forward not available.\n");
                break;
        }
        break;
   case Q_SF: 
        switch(column_type){
           case Q_SF: sSOR_step_forward_sf(tGrid,Z,b,x,om);
                break;
           default:
                eprintf("Error: SOR_step_forward not available.\n");
                break;
        }
        break;
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                         vSOR_step_forward_vf(tGrid,Z,b,x,om);
                      else if (structure & Q_FEBDIAG)
                         sSOR_step_forward_vf(tGrid,Z,b,x,om);
                      else
                         eprintf("Error: SOR_step_forward not available.\n");
                break;
           default:
                eprintf("Error: SOR_step_forward not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vSOR_step_forward_dvf(tGrid,Z,b,x,om);
                       else
                          eprintf("Error: SOR_step_forward not available.\n");
                break;
           default:
                eprintf("Error: SOR_step_forward not available.\n");
                break;
        }
        break;
   case Q_SNSF: 
        switch(column_type){
           case Q_SNSF: sSOR_step_forward_sn_sf(tGrid,Z,b,x,om);
                break;
           default:
                eprintf("Error: SOR_step_forward not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: SOR_step_forward not available.\n");
        break;
   }
}

void SOR_step_backward(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,
                       structure,b,p1,p2,p3,p4,p5,p6,p7,p8,p9,om,r2,r3)
GRID *tGrid;       /*  one backward step of SOR for the solution of Zx=b      */
INT Z, x, y, q;    /*  x ... initial guess, y ... next iterate, y=x possible  */
INT t_row, t_column, row_type, column_type, structure;      /*  row = column  */
INT b,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT om,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   if (x != y)
      copy(tGrid,x,y,t_row,row_type);
   switch(row_type){
   case Q_VF: 
        switch(column_type){
           case Q_VF: if (structure & Q_FULL)
                      /* vSOR_step_backward_vf(tGrid,Z,b,x,om); */
                         eprintf("Error: SOR_step_backward not available.\n");
                      else if (structure & Q_FEBDIAG)
                         sSOR_step_backward_vf(tGrid,Z,b,x,om);
                      else
                         eprintf("Error: SOR_step_backward not available.\n");
                break;
           default:
                eprintf("Error: SOR_step_backward not available.\n");
                break;
        }
        break;
   case Q_DVF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_DVFEBDIAG)
                          vSOR_step_backward_dvf(tGrid,Z,b,x,om);
                       else
                          eprintf("Error: SOR_step_backward not available.\n");
                break;
           default:
                eprintf("Error: SOR_step_backward not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: SOR_step_backward not available.\n");
        break;
   }
}

void SSOR_step(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure,
               b,p1,p2,p3,p4,p5,p6,p7,p8,p9,om,r2,r3)
GRID *tGrid;       /*  one step of SSOR for the solution of Zx=b              */
INT Z, x, y, q;    /*  x ... initial guess, y ... next iterate, y=x possible  */
INT t_row, t_column, row_type, column_type, structure;      /*  row = column  */
INT b,p1,p2,p3,p4,p5,p6,p7,p8,p9;
FLOAT om,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   SOR_step_forward(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,
                    structure,b,p1,p2,p3,p4,p5,p6,p7,p8,p9,om,r2,r3);
   SOR_step_backward(tGrid,Z,y,y,q,t_row,t_column,row_type,column_type,
                     structure,b,p1,p2,p3,p4,p5,p6,p7,p8,p9,om,r2,r3);
}

void multAW(tGrid,Z,x,y,q,q2,t,type,A_struct) /*  y := Z_symmGS x  */
GRID *tGrid;
INT Z, x, y, q, q2, t, type, A_struct;
{
   mult_AU(tGrid,Z,x,y,q,t,t,type,type,A_struct,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_inv_diag(tGrid,Z,y,q,q,t,t,type,type,A_struct);
   add(tGrid,x,q,y,t,type);
   mult_AL(tGrid,Z,y,q,q2,t,t,type,type,A_struct,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_diag(tGrid,Z,y,q2,q2,t,t,type,type,A_struct);
   add(tGrid,q,q2,y,t,type);
}

void mult_B(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure)  
GRID *tGrid;                                                    /*  y := B x  */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
/* q ... auxiliary vector of the type y (i.e., row_type); function modifies 
   values of q which do not correspond to degrees of freedom of y;
   not used here */
{
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_VNVF: multA_sn_X_vn_vf(tGrid,Z,x,y,t_column);
                        if (structure & Q_BDAUGMENT)
                           multD_sn_sf_X_vn_vf(tGrid,Z,x,y);
                break;
           case Q_VNVE: multA_sn_X_vn(tGrid,Z,x,y,t_column);
                        multA_sn_X_ve(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVFVE: multA_sn_X_vn_vf(tGrid,Z,x,y,t_column);
                          multA_sn_X_ve(tGrid,Z,x,y,t_column);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                          multA_sn_X_dvf0(tGrid,Z,x,y,t_column);
                       else
                          multA_sn_X_dvf(tGrid,Z,x,y,t_row,t_column);
                break;
           default:
                eprintf("Error: mult_B not available.\n");
                break;
        }
        break;
   case Q_SF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_ZERO_DIAG)
                          multA_sf_X_dvf_zero_diag(tGrid,Z,x,y,t_row);
                       else
                          eprintf("Error: mult_B not available.\n");
                break;
           default:
                eprintf("Error: mult_B not available.\n");
                break;
        }
        break;
   case Q_SE: 
        switch(column_type){
           case Q_VN: multA_e_X_vn(tGrid,Z,x,y,t_row);
                break;
           case Q_SF: multA_e_X_sf(tGrid,Z,x,y,t_row);
                break;
           case Q_VF: multA_e_X_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                         multA_e_X_dvf(tGrid,Z,x,y,t_column);
                      else
                         eprintf("Error: mult_B not available.\n");
                break;
           case Q_VNSF: multA_e_X_vn_sf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVF: multA_e_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVFVE: if (structure & Q_NE_ZERO)
                        multA_e_X_vn_vf(tGrid,Z,x,y,t_column);
                       else
                          eprintf("Error: mult_B not available.\n");
                break;
           case Q_VNSFLG: multB_lg(tGrid,Z,x,y);
                break;
           default:
                eprintf("Error: mult_B not available.\n");
                break;
        }
        break;
   case Q_SNE: 
        switch(column_type){
           case Q_DVF: multA_sef_X_dvf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVF: multA_sef_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVFVE: multA_sef_X_vn_vf_ve(tGrid,Z,x,y,t_column);
                break;
           default:
                eprintf("Error: mult_B not available.\n");
                break;
        }
        break;
   case Q_SNSSNSSF: 
        switch(column_type){
           case Q_VNVF: multA_sn_X_vn_vf(tGrid,Z,x,y,t_column);
                        multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           case Q_VNVFVE: multA_sn_X_vn_vf(tGrid,Z,x,y,t_column);
                          multA_sn_X_ve(tGrid,Z,x,y,t_column);
                          multA_ssn_ssf_X_vn_vf(tGrid,Z,x,y,t_column);
                break;
           default:
                eprintf("Error: mult_B not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_B not available.\n");
        break;
   }
}

void mult_BT(tGrid,Z,x,y,q,t_row,t_column,row_type,column_type,structure)
GRID *tGrid;                                                  /*  y := B^T x  */
INT Z, x, y, q, t_row, t_column, row_type, column_type, structure;
/* q ... auxiliary vector of the type y (i.e., column_type); function modifies 
   values of q which do not correspond to degrees of freedom of y */
{
   if (t_row & WITHOUT_FIRST)
      set_first(tGrid,0.,x,t_row,row_type);
   switch(row_type){
   case Q_SN: 
        switch(column_type){
           case Q_VNVF: multA_vn_vf_X_sn(tGrid,Z,x,y);
                        if (structure & Q_BDAUGMENT)
                           multD_vn_vf_X_sn_sf(tGrid,Z,x,y);
                break;
           case Q_VNVE: multA_vn_X_sn(tGrid,Z,x,y,t_row);
                        multA_ve_X_sn(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVFVE: multA_vn_vf_X_sn(tGrid,Z,x,y);
                          multA_ve_X_sn(tGrid,Z,x,y,t_row);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                          multA_dvf0_X_sn(tGrid,Z,x,y,t_row);
                       else
                          multA_dvf_X_sn(tGrid,Z,x,y);
                break;
           default:
                eprintf("Error: mult_BT not available.\n");
                break;
        }
        break;
   case Q_SF: 
        switch(column_type){
           case Q_DVF: if (structure & Q_ZERO_DIAG)
                          multA_dvf_X_sf_zero_diag(tGrid,Z,x,y);
                       else
                          eprintf("Error: mult_BT not available.\n");
                break;
           default:
                eprintf("Error: mult_BT not available.\n");
                break;
        }
        break;
   case Q_SE: 
        switch(column_type){
           case Q_VN: multA_vn_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_SF: multA_sf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_VF: multA_vf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_DVF: if (structure & Q_FIRST_DV)
                         multA_dvf_X_e(tGrid,Z,x,y,q,t_row);
                      else
                         eprintf("Error: mult_BT not available.\n");
                break;
           case Q_VNSF: multA_vn_sf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_VNVF: multA_vn_vf_X_e(tGrid,Z,x,y,q,t_row);
                break;
           case Q_VNVFVE: if (structure & Q_NE_ZERO){
                        multA_vn_vf_X_e(tGrid,Z,x,y,q,t_row);
                        set_value(tGrid,0.,y,t_column,Q_VE);
                      }
                      else
                         eprintf("Error: mult_BT not available.\n");
                break;
           case Q_VNSFLG: multBT_lg(tGrid,Z,x,y,q);
                break;
           default:
                eprintf("Error: mult_BT not available.\n");
                break;
        }
        break;
   case Q_SNE: 
        switch(column_type){
           case Q_DVF: multA_dvf_X_sef(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVF: multA_vn_vf_X_sef(tGrid,Z,x,y,t_row);
                break;
           case Q_VNVFVE: multA_vn_vf_ve_X_sef(tGrid,Z,x,y,t_row);
                break;
           default:
                eprintf("Error: mult_BT not available.\n");
                break;
        }
        break;
   case Q_SNSSNSSF: 
        switch(column_type){
           case Q_VNVF: multA_vn_vf_X_sn(tGrid,Z,x,y);
                        add_multA_vn_vf_X_ssn_ssf(tGrid,Z,x,y);
                break;
           case Q_VNVFVE: multA_vn_vf_X_sn(tGrid,Z,x,y);
                          multA_ve_X_sn(tGrid,Z,x,y,t_row);
                          add_multA_vn_vf_X_ssn_ssf(tGrid,Z,x,y);
                break;
           default:
                eprintf("Error: mult_BT not available.\n");
                break;
        }
        break;
   default:
        eprintf("Error: mult_BT not available.\n");
        break;
   }
}

/* y := (Z + r B^T B) x */
void mult_Ar(tGrid,matrix_a,x,y,q,t_row,t_column,row_type,column_type,a_struct,
             matrix_b,b_struct,w,v,t_for_p,p_type,p6,p7,p8,p9,r,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_u; row_type = column_type = u_type */
INT x, y, q, w;  /*  u_type; w is auxilaiary vector; for q see below  */
INT v;           /*  p_type; auxiliary vector  */
INT matrix_a,t_row,t_column,row_type,column_type,a_struct,
    matrix_b,b_struct,t_for_p,p_type,p6,p7,p8,p9;
FLOAT r,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   mult_B(tGrid,matrix_b,x,v,v,t_for_p,t_row,p_type,row_type,b_struct);
   mult_BT(tGrid,matrix_b,v,w,q,t_for_p,t_row,p_type,row_type,b_struct);
   mult_A(tGrid,matrix_a,x,y,q,t_row,t_row,row_type,row_type,a_struct,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_and_add(tGrid,r,w,y,y,t_row,row_type);
}

/* y := d^{-1} (Z + r B^T B) d^{-1} x */
void Pmult_Ar(tGrid,matrix_a,x,y,q,t_row,t_column,row_type,column_type,a_struct,
             matrix_b,b_struct,w,v,t_for_p,p_type,d,p7,p8,p9,r,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_u; row_type = column_type = u_type */
INT x, y, q, w;  /*  u_type; w is auxilaiary vector; for q see below  */
INT v;           /*  p_type; auxiliary vector  */
INT d;           /*  u_type; diagonal preconditioner  */
INT matrix_a,t_row,t_column,row_type,column_type,a_struct,
    matrix_b,b_struct,t_for_p,p_type,p7,p8,p9;
FLOAT r,r2,r3;
/* q ... auxiliary vector of the type y; function modifies values of q which
   do not correspond to degrees of freedom of y */
{
   divide(tGrid,x,d,w,t_row,row_type);
   mult_A(tGrid,matrix_a,w,y,q,t_row,t_row,row_type,row_type,a_struct,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_B(tGrid,matrix_b,w,v,v,t_for_p,t_row,p_type,row_type,b_struct);
   mult_BT(tGrid,matrix_b,v,w,q,t_for_p,t_row,p_type,row_type,b_struct);
   mult_and_add(tGrid,r,w,y,y,t_row,row_type);
   divide(tGrid,y,d,y,t_row,row_type);
}

void subtract_Acontribution_from_bc(tGrid,Z,u,f,q,r,s,t,type,a_struct)
GRID *tGrid;
INT Z, u, f, q, r, s, t, type, a_struct;
{
 if (MOVING_BOUNDARY == YES && (t & USE_IS_DIR)){
   mso_scopy(tGrid,u,q,USE_IS_DIR | STOP_IS_FIRST_INNER);
   set_value(tGrid,0.,q,t,type);
   mult_A(tGrid,Z,q,r,s,USE_IS_DIR | T_FOR_BC,USE_IS_DIR | T_FOR_BC,type,type,a_struct,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   subtr(tGrid,f,r,f,t,type);
 }
 else{
   copy(tGrid,u,q,STOP_IS_FIRST_INNER,type);
   set_value(tGrid,0.,q,t,type);
   mult_A(tGrid,Z,q,r,s,T_FOR_BC,T_FOR_BC,type,type,a_struct,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   subtr(tGrid,f,r,f,t,type);
 }
}

void make_Bcontribution_from_bc(tGrid,Z,u,f,q,r,t_u,t_p,u_type,p_type,b_struct)
GRID *tGrid;
INT Z, u, f, q, r, t_u, t_p, u_type, p_type, b_struct;
{
   FLOAT flux;
   INT n;

   copy(tGrid,u,q,STOP_IS_FIRST_INNER,u_type);
   set_value(tGrid,0.,q,t_u,u_type);
   mult_B(tGrid,Z,q,f,r,t_p,T_FOR_BC,p_type,u_type,b_struct);
   inv(tGrid,f,f,t_p,p_type);
   if (t_p & ZERO_MEAN){
      flux = sum_n(tGrid,f,&n,t_p,p_type);
      printf("Flux of the b.c.: %e\n",flux);
      if (fabs(flux) > EPSC){
         eprintf("Error: flux of the b.c. > EPSC. RHS changed.\n");
         add_value(tGrid,-flux/n,f,t_p,p_type); 
         printf("     changed rhs: %e\n",sum(tGrid,f,t_p,p_type));
      }
   }
}

void defect(tGrid,Z,f,u,d,q,t,type,structure)  /*  d  = A u - f  */
GRID *tGrid;
INT Z, f, u, d, q, t, type, structure;
{
   mult_A(tGrid,Z,u,d,q,t,t,type,type,structure,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   subtr(tGrid,d,f,d,t,type);
}

void Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,
                   t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct)
GRID *tGrid;
INT ZA, ZB, f, fe, u, ue, d, de, q, 
    t_u, t_p, u_type, p_type, A_struct, B_struct, C_struct;
{                                                  /*  d  = A u + B^T ue - f  */
                                                   /*  de = B u - C ue - fe   */
   mult_A(tGrid,ZA,u,d,q,t_u,t_u,u_type,u_type,A_struct,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_BT(tGrid,ZB,ue,q,q,t_p,t_u,p_type,u_type,B_struct);
   mult_B(tGrid,ZB,u,de,q,t_p,t_u,p_type,u_type,B_struct);
   add_and_subtr(tGrid,d,q,f,d,t_u,u_type);
   subtr(tGrid,de,fe,de,t_p,p_type);
   if (C_struct){
      mult_A(tGrid,ZA,ue,q,q,t_p,t_p,p_type,p_type,C_struct,
             0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
      subtr(tGrid,de,q,de,t_p,p_type);
   }
}

void mult_by_Stokes_matrix(tGrid,ZA,ZB,u,ue,d,de,q,t_u,t_p,u_type,p_type,A_struct,B_struct)
GRID *tGrid;
INT ZA, ZB, u, ue, d, de, q, t_u, t_p, u_type, p_type, A_struct, B_struct;
{                                                  /*  d  = A u + B^T ue  */
                                                   /*  de = B u           */
   mult_A(tGrid,ZA,u,d,q,t_u,t_u,u_type,u_type,A_struct,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_BT(tGrid,ZB,ue,q,q,t_p,t_u,p_type,u_type,B_struct);
   mult_B(tGrid,ZB,u,de,q,t_p,t_u,p_type,u_type,B_struct);
   add(tGrid,d,q,d,t_u,u_type);
}

/******************************************************************************/
/******************************************************************************/

#if F_DATA & DVECTOR_FACE_DATA

FLOAT max_eigenvalue(tGrid,Z,x,y,t,type,structure) /*  Z is a square matrix  */
GRID *tGrid;                                           
INT Z, x, y, t, type, structure;
{
   FACE *theFace;
   FLOAT c, lambda, lambda_old, lambda_max=0.;
   INT i, j;

   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      for (i=0; i < 2; i++)
         for (j=0; j < 2; j++){
            set_value(tGrid,0.,x,t,type);
            FDDV(theFace,x,i,j) = 1.;
            lambda_old = 0.;
            lambda = 1.;
            while( fabs((lambda-lambda_old)/lambda) > 1.e-6 ){
               lambda_old = lambda;
               mult_A(tGrid,Z,x,y,y,t,t,type,type,structure,
                      0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
               c = dot(tGrid,x,x,t,type);
               lambda = dot(tGrid,x,y,t,type)/c;
               if (c > 1.e20)
                  mult(tGrid,1.e-10,y,x,t,type);
               else
                  copy(tGrid,y,x,t,type);
            }
            printf("%e\n",lambda);
            if (lambda > lambda_max)
               lambda_max = lambda;
         }
   printf("lambda_max = %e\n",lambda_max);
   return(lambda_max);
}

#else

FLOAT max_eigenvalue(tGrid,Z,x,y,t,type,structure)
GRID *tGrid; INT Z, x, y, t, type, structure;
{ eprintf("Error: max_eigenvalue not available.\n"); return(0.); }

#endif
