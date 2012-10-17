/******************************************************************************/
/*                                                                            */
/*                             ILU decomposition                              */
/*                                                                            */
/******************************************************************************/

void fill_Ax();
void make_vector_from_grid_data();
void make_grid_data_from_vector();
void solve_system_using_umfpack();

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA)

void scopy_mat_to_mat(theGrid,a,h)
GRID *theGrid;
INT a, h;
{
   NODE *pnode;
   LINK *plink;

   for (pnode=FIRSTN(theGrid); pnode!= NULL; pnode=SUCC(pnode)){
         COEFFS(pnode,h) = COEFFS(pnode,a);
         for (plink=TSTART(pnode); plink!=NULL; plink=NEXT(plink))
            COEFFLS(plink,h) = COEFFLS(plink,a);
      }
} 

void sILU(theGrid,a,h,t)
GRID *theGrid;
INT a, h, t;
{
   NODE *pnr, *pni, *pnj;
   LINK *plri, *plir, *plrj, *plij;
   FLOAT pivot;
    
   scopy_mat_to_mat(theGrid,a,h);

   if (t & USE_IS_DIR){
      for (pnr = FIRSTN(theGrid); pnr != NULL; pnr = pnr->succ)
         if (!IS_DIR(pnr,t))
            for (plri = TSTART(pnr); plri != NULL; plri = plri->next){
               pni = NBNODE(plri);
               if ( !IS_DIR(pni,t) && INCR(pnr,pni) ){
                  for (plir = TSTART(pni); NBNODE(plir) != pnr; plir = plir->next);
                  pivot = COEFFLS(plir,h) /= COEFFS(pnr,h);
                  COEFFS(pni,h) -= pivot*COEFFLS(plri,h);
                  for (plrj = TSTART(pnr); plrj != NULL; plrj=plrj->next){
                     pnj = NBNODE(plrj);
                     if (!IS_DIR(pnj,t) && INCR(pnr,pnj) && pnj != pni){
                        for (plij = TSTART(pni); plij != NULL && 
                                        NBNODE(plij) != pnj; plij = plij->next);
                        if (plij != NULL)
                           COEFFLS(plij,h) -= pivot*COEFFLS(plrj,h); 
                     }
                  }
               }
            }
   }
   else
      for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ)
         for (plri = START(pnr); plri != NULL; plri = plri->next){
            pni = NBNODE(plri);
            if (INCR(pnr,pni)){
               for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
               pivot = COEFFLS(plir,h) /= COEFFS(pnr,h);
               COEFFS(pni,h) -= pivot*COEFFLS(plri,h);
               for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){
                  pnj = NBNODE(plrj);
                  if (INCR(pnr,pnj) && pnj != pni){
                     for (plij = START(pni); plij != NULL && 
                                        NBNODE(plij) != pnj; plij = plij->next);
                     if (plij != NULL)
                        COEFFLS(plij,h) -= pivot*COEFFLS(plrj,h); 
                  }
               }
            }
         }
}
 
/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void ssolveILU(theGrid,h,x,y,t) 
GRID *theGrid;
INT h, x, y, t;
{    
   NODE *pni, *pnj;
   LINK *plij;
   FLOAT sum;

   if (t & USE_IS_DIR){
      /* find x: L x = y */
      for (pni = FIRSTN(theGrid); pni != NULL; pni = pni->succ)
         if (!IS_DIR(pni,t)){
            sum = 0.0;
            for (plij = TSTART(pni); plij != NULL; plij = plij->next){
               pnj = NBNODE(plij);
               if ( !IS_DIR(pnj,t) && INCR(pnj,pni) )
                  sum += COEFFLS(plij,h)*NDS(pnj,x);
            }
            NDS(pni,x) = NDS(pni,y) - sum;
         }

      /* find solution of U . = x and store it in x */
      for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
         if (IS_FN(pni) && !IS_DIR(pni,t)){
            sum = 0.0;
            for (plij = START(pni); plij != NULL; plij = plij->next){
               pnj = NBNODE(plij);
               if ( !IS_DIR(pnj,t) && INDEX(pnj) > INDEX(pni) )
                  sum += COEFFLS(plij,h)*NDS(pnj,x);
            }
            NDS(pni,x) = (NDS(pni,x)-sum)/COEFFS(pni,h);
         }

      for (pni = theGrid->ldbn; pni >= FIRSTN(theGrid); pni--)
         if (NOT_FN(pni) && !IS_DIR(pni,t)){
            sum = 0.0;
            for (plij = TSTART(pni); plij != NULL; plij = plij->next){
               pnj = NBNODE(plij);
               if ( !IS_DIR(pnj,t) && INCR(pni,pnj) )
                  sum += COEFFLS(plij,h)*NDS(pnj,x);
            }
            NDS(pni,x) = (NDS(pni,x)-sum)/COEFFS(pni,h);
         }
      }
   else{
      /* find x: L x = y */
      for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
         sum = 0.;
         for (plij = START(pni); plij != NULL; plij = plij->next){
            pnj = NBNODE(plij);
            if (INCR(pnj,pni) )
               sum += COEFFLS(plij,h)*NDS(pnj,x);
         }
         NDS(pni,x) = NDS(pni,y) - sum;
      }

      /* find solution of U . = x and store it in x */
      for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
         if (IS_FN(pni)){
            sum = 0.;
            for (plij = START(pni); plij != NULL; plij = plij->next){
               pnj = NBNODE(plij);
               if (INDEX(pnj) > INDEX(pni) )
                  sum += COEFFLS(plij,h)*NDS(pnj,x);
            }
            NDS(pni,x) = (NDS(pni,x)-sum)/COEFFS(pni,h);
         }
   }
}

void ILU_n_beta_sparse(theGrid,a,h,q,beta)
GRID *theGrid;
INT a, h, q;
FLOAT beta;
{
   NODE *pni, *pnj, *pnk;
   LINK *plij, *plik, *plkj, *plki;
   FLOAT s;
   INT ind;
   
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      s = COEFFN(pni,h)-COEFFN(pni,a);
      ind = INDEX(pni)-1; 
      for (plik = START(pni); plik != NULL; plik = plik->next)
         if ( INDEX(pnk = NBNODE(plik)) <= ind){
            for (plki = START(pnk); NBNODE(plki) != pni; plki= plki->next);
            s += COEFFL(plik,h)*COEFFL(plki,h);
         }
      NDS(pni,q) = fabs(s);
      for (plij = START(pni); plij != NULL; plij = plij->next){
         pnj = NBNODE(plij); 
         s = -COEFFL(plij,a);
         if (INDEX(pni) > INDEX(pnj))
            ind = INDEX(pnj);
         else{
            ind = INDEX(pni)-1; 
            s += COEFFL(plij,h); 
         }
         for (plik = START(pni); plik != NULL; plik = plik->next)
            if ( INDEX(pnk = NBNODE(plik)) <= ind){
               for (plkj = START(pnk); plkj != NULL && NBNODE(plkj) != pnj;
                                                           plkj = plkj->next);
               if (plkj != NULL)
                  s += COEFFL(plik,h)*COEFFL(plkj,h);
            }
         NDS(pni,q) += fabs(s);
      }
   }
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      s = COEFFN(pni,h) + NDS(pni,q)*beta;
      NDS(pni,q) = COEFFN(pni,h)/s;
      COEFFN(pni,h) = s;
   }
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ)
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            COEFFL(plij,h) *= NDS(pnj,q);
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) &&
                (N_DATA & SCALAR_NODE_DATA))  */

void sILU(theGrid,a,h,t)
GRID *theGrid; INT a, h, t;
{  eprintf("Error: sILU not available.\n");  }

void ssolveILU(theGrid,h,x,y,t) 
GRID *theGrid; INT h, x, y, t;
{  eprintf("Error: ssolveILU not available.\n");  }

void ILU_n_beta_sparse(theGrid,a,h,q,beta)
GRID *theGrid; INT a, h, q; FLOAT beta;
{  eprintf("Error: ILU_n_beta_sparse not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES)

void nn_copy_mat_to_mat(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   LINK *plink;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      COEFFN(pnode,h) = COEFFN(pnode,a);
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         COEFFL(plink,h) = COEFFL(plink,a);
   }
}
 
void nn_ILU(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   LINK *plri, *plir, *plrj, *plij;
   FLOAT pivot;
   
   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni = NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot = COEFFL(plir,h) /= COEFFN(pnr,h);
            COEFFN(pni,h) -= pivot*COEFFL(plri,h);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
               if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni ){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj; 
                                                            plij = plij->next);
                  if (plij != NULL)
                     COEFFL(plij,h) -= pivot*COEFFL(plrj,h); 
              }
         }
   }         
}   

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES))  */

void nn_copy_mat_to_mat(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: nn_copy_mat_to_mat not available.\n");  }

void nn_ILU(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: nn_ILU not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA)

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void nn_solveILU(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{    
   NODE *pni, *pnj;
   LINK *plij;
   FLOAT sum[DIM];
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      SET7(sum,0.)
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            SET4(sum,NDD(pnj,x),COEFFL(plij,h))
      SET11(NDD(pni,x),NDD(pni,y),sum)
   }
   
   /* find solution of U . = x and store it in x */
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         SET7(sum,0.)
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               SET4(sum,NDD(pnj,x),COEFFL(plij,h))
         SET12(NDD(pni,x),NDD(pni,x),sum,COEFFN(pni,h))
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && 
                (N_DATA & VECTOR_NODE_DATA)  */

void nn_solveILU(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{  eprintf("Error: nn_solveILU not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)

void scopy_mat_to_mat_sn_sf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pflink;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      COEFFN(pnode,h) = COEFFN(pnode,a);
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         COEFFL(plink,h) = COEFFL(plink,a);
      for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl=pnfl->next)
         COEFFL(pnfl,h) = COEFFL(pnfl,a);
   }
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         COEFFL(pflink,h) = COEFFL(pflink,a);
      for (pfnl = FNSTART(pface); pfnl != NULL; pfnl=pfnl->next)
         COEFFL(pfnl,h) = COEFFL(pfnl,a);
   }
}
 
void sILU_sn_sf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   NFLINK *pnfri,*pnfrj, *pnfij;
   FNLINK *pfnir, *pfnij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot;
   
   scopy_mat_to_mat_sn_sf(theGrid,a,h);

   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni = NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot = COEFFL(plir,h) /= COEFFN(pnr,h);
            COEFFN(pni,h) -= pivot*COEFFL(plri,h);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
               if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni ){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj; 
                                                           plij = plij->next);
                  if (plij != NULL)
                     COEFFL(plij,h) -= pivot*COEFFL(plrj,h); 
              }
            for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
               pfj = NBFACE(pnfrj);
               for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) !=pfj;
                                                          pnfij = pnfij->next);
               if (pnfij != NULL)
                  COEFFL(pnfij,h) -= pivot*COEFFL(pnfrj,h);
            }
         }
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         pfi = NBFACE(pnfri);
         for (pfnir = FNSTART(pfi); NBNODE(pfnir) !=  pnr; pfnir = pfnir->next);
         pivot = COEFFL(pfnir,h) /= COEFFN(pnr,h);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
            if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) ){
               for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj;
                                                           pfnij = pfnij->next);
               if (pfnij != NULL)
                  COEFFL(pfnij,h) -= pivot*COEFFL(plrj,h);
            }
         COEFF_FF(pfi,h) -= pivot*COEFFL(pnfri,h);
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                          pflij = pflij->next);
               if (pflij != NULL)
                  COEFFL(pflij,h) -= pivot*COEFFL(pnfrj,h);
            } 
      }        
   }         
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFFL(pflir,h) /= COEFF_FF(pfr,h);
           COEFF_FF(pfi,h) -= pivot*COEFF_FL(pflri,h);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL)
                   COEFF_FL(pflij,h) -= pivot*COEFFL(pflrj,h);
             }
         }         
}   

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && 
                (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES))  */

void sILU_sn_sf(theGrid,a,h)
GRID *theGrid; INT a,h;
{   eprintf("Error: sILU_sn_sf not available.\n");   }

#endif 

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void ssolveILU_sn_sf(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{    
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum;
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = 0.;
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFL(plij,h)*NDS(pnj,x);
      NDS(pni,x) = NDS(pni,y) - sum;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sum = 0.;
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         sum += COEFFL(pfnij,h)*NDS(NBNODE(pfnij),x);
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sum += COEFFL(pflij,h)*FD(pfj,x);
      FD(pfi,x) = FD(pfi,y) - sum;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sum = 0.;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sum += COEFFL(pflij,h)*FD(pfj,x);
         FD(pfi,x) = (FD(pfi,x) - sum)/COEFF_FF(pfi,h);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = 0.;
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFFL(pnfij,h)*FD(NBFACE(pnfij),x);
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFL(plij,h)*NDS(pnj,x);
         NDS(pni,x) = (NDS(pni,x) - sum)/COEFFN(pni,h);
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && 
                (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void ssolveILU_sn_sf(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{   eprintf("Error: ssolveILU_sn_sf not available.\n");   }

#endif 
 
#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

/* Matrix h contains an ILU-decomposition; L x = y is solved.  */           
void ssolveIL_vn_vf(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{    
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum[DIM];
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      SET7(sum,0.)
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            SET4(sum,NDD(pnj,x),COEFFL(plij,h))
      SET11(NDD(pni,x),NDD(pni,y),sum)
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      SET7(sum,0.)
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         SET4(sum,NDD(NBNODE(pfnij),x),COEFFL(pfnij,h))
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            SET4(sum,FDVP(pfj,x),COEFFL(pflij,h))
      SET11(FDVP(pfi,x),FDVP(pfi,y),sum)
   }
}

/* Matrix h contains an ILU-decomposition; U . = x is solved.  */           
void ssolveIU_vn_vf(theGrid,h,x) 
GRID *theGrid;
INT h,x;
{    
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FLINK *pflij;
   FLOAT sum[DIM];
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         SET7(sum,0.)
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               SET4(sum,FDVP(pfj,x),COEFFL(pflij,h))
         SET12(FDVP(pfi,x),FDVP(pfi,x),sum,COEFF_FF(pfi,h))
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         SET7(sum,0.)
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            SET4(sum,FDVP(NBFACE(pnfij),x),COEFFL(pnfij,h))
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               SET4(sum,NDD(pnj,x),COEFFL(plij,h))
         SET12(NDD(pni,x),NDD(pni,x),sum,COEFFN(pni,h))
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && 
                (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

void ssolveIL_vn_vf(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{   eprintf("Error: ssolveIL_vn_vf not available.\n");   }
    
void ssolveIU_vn_vf(theGrid,h,x) 
GRID *theGrid; INT h,x;
{   eprintf("Error: ssolveIU_vn_vf not available.\n");   }

#endif 
 
#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR)

void sILU_sn_sf_se(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   ELEMENT *pel;
   NODE *pnr;
   FACE *pfr;
   LINK *plri, *plir;
   NFLINK *pnfri;
   FNLINK *pfnir;
   FLINK *pflri, *pflir;
   FLOAT pivot;
   INT i, r;

   for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
      SSET1(COEFF_NEP(pel,h),COEFF_NEP(pel,a))
      SSET1(COEFF_ENP(pel,h),COEFF_ENP(pel,a))
      SSET1(COEFF_FEP(pel,h),COEFF_FEP(pel,a))
      SSET1(COEFF_EFP(pel,h),COEFF_EFP(pel,a))
      COEFF_EE(pel,h) = COEFF_EE(pel,a);
   }
   sILU_sn_sf(theGrid,a,h);
   
   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
         GIVE_INDEX(pel->n,pnr,r);
         if (r > -1){
            pivot = COEFF_EN(pel,h,r) /= COEFFN(pnr,h);
            COEFF_EE(pel,h) -= pivot*COEFF_NE(pel,h,r);
            for (i = 0; i < NVERT; i++){
               if (INDEX(pel->n[i]) > INDEX(pnr) && IS_FN(pel->n[i])){
                  for (plir = START(pel->n[i]); NBNODE(plir) != pnr; 
                                                             plir = plir->next);
                  for (plri = START(pnr); NBNODE(plri) != pel->n[i];
                                                             plri = plri->next);
                  COEFF_NE(pel,h,i) -= COEFFL(plir,h)*COEFF_NE(pel,h,r);
                  COEFF_EN(pel,h,i) -= pivot*COEFFL(plri,h); 
               }
               if (IS_FF(pel->f[i])){
                  for (pfnir = FNSTART(pel->f[i]); NBNODE(pfnir) != pnr; 
                                                           pfnir = pfnir->next);
                  for (pnfri = NFSTART(pnr); NBFACE(pnfri) != pel->f[i];
                                                           pnfri = pnfri->next);
                  COEFF_FE(pel,h,i) -= COEFFL(pfnir,h)*COEFF_NE(pel,h,r);
                  COEFF_EF(pel,h,i) -= pivot*COEFFL(pnfri,h); 
               } 
            } 
         } 
      }
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
         GIVE_INDEX(pel->f,pfr,r);
         if (r > -1){
            pivot = COEFF_EF(pel,h,r) /= COEFF_FF(pfr,h);
            COEFF_EE(pel,h) -= pivot*COEFF_FE(pel,h,r);
            for (i = 0; i < SIDES; i++)
               if (INDEX(pel->f[i]) > INDEX(pfr) && IS_FF(pel->f[i])){
                  for (pflir = FSTART(pel->f[i]); NBFACE(pflir) != pfr; 
                                                           pflir = pflir->next);
                  for (pflri = FSTART(pfr); NBFACE(pflri) != pel->f[i];
                                                           pflri = pflri->next);
                  COEFF_FE(pel,h,i) -= COEFFL(pflir,h)*COEFF_FE(pel,h,r);
                  COEFF_EF(pel,h,i) -= pivot*COEFFL(pflri,h); 
               } 
         } 
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && 
                (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && 
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && 
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && 
                (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && 
                (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && 
                (E_DATA & ExE_MATR))  */

void sILU_sn_sf_se(theGrid,a,h)
GRID *theGrid; INT a,h;
{   eprintf("Error: sILU_sn_sf_se not available.\n");   }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void ssolveILU_vn_vf_ve(theGrid,h,x,y)
GRID *theGrid;
INT h,x,y;
{
   ELEMENT *pel;

   set_value(theGrid,0.,x,0,Q_VNVF);
   ssolveIL_vn_vf(theGrid,h,x,y);
   for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
      SET1(EDVP(pel,x),EDVP(pel,y))
      SUBTR_SUM(pel->n,x,EDVP(pel,x),COEFF_ENP(pel,h))
      SUBTR_SUM(pel->f,x,EDVP(pel,x),COEFF_EFP(pel,h))
      SET13(EDVP(pel,x),COEFF_EE(pel,h))
   }
   for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ){
      SUBTR_VMULT(pel->n,x,COEFF_NEP(pel,h),EDVP(pel,x))
      SUBTR_VMULT(pel->f,x,COEFF_FEP(pel,h),EDVP(pel,x))
   }
   ssolveIU_vn_vf(theGrid,h,x);
}

#else  /*  if !((E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && 
                (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && 
                (E_DATA & ExE_MATR) && (N_DATA & VECTOR_NODE_DATA) &&
                (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA))*/

void ssolveILU_vn_vf_ve(theGrid,h,x,y)
GRID *theGrid; INT h,x,y;
{   eprintf("Error: ssolveILU_vn_vf_ve not available.\n");   }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void scopy_mat_to_mat_vn_sf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pflink;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      COEFFN(pnode,h) = COEFFN(pnode,a);
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         COEFFL(plink,h) = COEFFL(plink,a);
      for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl=pnfl->next)
         SET1(COEFF_NFP(pnfl,h),COEFF_NFP(pnfl,a));
   }
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         COEFF_FL(pflink,h) = COEFF_FL(pflink,a);
      for (pfnl = FNSTART(pface); pfnl != NULL; pfnl=pfnl->next)
         SET1(COEFF_FNP(pfnl,h),COEFF_FNP(pfnl,a));
   }
}
 
void sILU_vn_sf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   NFLINK *pnfri,*pnfrj, *pnfij;
   FNLINK *pfnir, *pfnij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot, pivots[DIM];
   
   scopy_mat_to_mat_vn_sf(theGrid,a,h);

   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni = NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot = COEFFL(plir,h) /= COEFFN(pnr,h);
            COEFFN(pni,h) -= pivot*COEFFL(plri,h);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
               if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni ){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj; 
                                                           plij = plij->next);
                  if (plij != NULL)
                     COEFFL(plij,h) -= pivot*COEFFL(plrj,h); 
              }
            for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
               pfj = NBFACE(pnfrj);
               for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) !=pfj;
                                                          pnfij = pnfij->next);
               if (pnfij != NULL)
                  SET9(COEFF_NFP(pnfij,h),COEFF_NFP(pnfrj,h),pivot)
            }
         }
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         pfi = NBFACE(pnfri);
         for (pfnir = FNSTART(pfi); NBNODE(pfnir) !=  pnr; pfnir = pfnir->next);
         SET13(COEFF_FNP(pfnir,h),COEFFN(pnr,h))
         SET1(pivots,COEFF_FNP(pfnir,h));
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
            if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) ){
               for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj;
                                                           pfnij = pfnij->next);
               if (pfnij != NULL)
                  SET9(COEFF_FNP(pfnij,h),pivots,COEFFL(plrj,h)) 
            }
         COEFF_FF(pfi,h) -= DOT(pivots,COEFF_NFP(pnfri,h));
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                          pflij = pflij->next);
               if (pflij != NULL)
                  COEFF_FL(pflij,h) -= DOT(pivots,COEFF_NFP(pnfrj,h));
            } 
      }        
   }         
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFF_FL(pflir,h) /= COEFF_FF(pfr,h);
           COEFF_FF(pfi,h) -= pivot*COEFF_FL(pflri,h);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL)
                   COEFF_FL(pflij,h) -= pivot*COEFF_FL(pflrj,h);
             }
         }         
}   

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void ssolveILU_vn_sf(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{    
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum[DIM], sumf;
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      SET7(sum,0.)
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            SET4(sum,NDD(pnj,x),COEFFL(plij,h))
      SET11(NDD(pni,x),NDD(pni,y),sum)
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sumf = 0.0;
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         sumf += DOT(COEFF_FNP(pfnij,h),NDD(NBNODE(pfnij),x));
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sumf += COEFF_FL(pflij,h)*FD(pfj,x);
      FD(pfi,x) = FD(pfi,y) - sumf;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sumf = 0.0;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sumf += COEFF_FL(pflij,h)*FD(pfj,x);
         FD(pfi,x) = (FD(pfi,x) - sumf)/COEFF_FF(pfi,h);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         SET7(sum,0.)
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            SET4(sum,COEFF_NFP(pnfij,h),FD(NBFACE(pnfij),x))
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               SET4(sum,NDD(pnj,x),COEFFL(plij,h))
         SET12(NDD(pni,x),NDD(pni,x),sum,COEFFN(pni,h))
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
                (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) &&
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void sILU_vn_sf(theGrid,a,h)
GRID *theGrid; INT a,h;
{   eprintf("Error: sILU_vn_sf not available.\n");   }

void ssolveILU_vn_sf(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{   eprintf("Error: ssolveILU_vn_sf not available.\n");   }

#endif 
 
#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & F_LINK_TO_FACES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void scopy_mat_to_mat_vn_sf_nred(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLINK *pflink;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      COEFFN(pnode,h) = COEFFN(pnode,a);
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         COEFFL(plink,h) = COEFFL(plink,a);
   }
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         COEFF_FL(pflink,h) = COEFF_FL(pflink,a);
   }
}
 
void sILU_vn_sf_nred(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot;
   
   scopy_mat_to_mat_vn_sf_nred(theGrid,a,h);

   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni = NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot = COEFFL(plir,h) /= COEFFN(pnr,h);
            COEFFN(pni,h) -= pivot*COEFFL(plri,h);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
               if ( INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni ){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj; 
                                                           plij = plij->next);
                  if (plij != NULL)
                     COEFFL(plij,h) -= pivot*COEFFL(plrj,h); 
              }
         }
   }         
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFF_FL(pflir,h) /= COEFF_FF(pfr,h);
           COEFF_FF(pfi,h) -= pivot*COEFF_FL(pflri,h);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL)
                   COEFF_FL(pflij,h) -= pivot*COEFF_FL(pflrj,h);
             }
         }         
}   

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void ssolveILU_vn_sf_nred(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   FLINK *pflij;
   FLOAT sum[DIM], sumf;
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      SET7(sum,0.)
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            SET4(sum,NDD(pnj,x),COEFFL(plij,h))
      SET11(NDD(pni,x),NDD(pni,y),sum)
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sumf = 0.0;
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sumf += COEFF_FL(pflij,h)*FD(pfj,x);
      FD(pfi,x) = FD(pfi,y) - sumf;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sumf = 0.0;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sumf += COEFF_FL(pflij,h)*FD(pfj,x);
         FD(pfi,x) = (FD(pfi,x) - sumf)/COEFF_FF(pfi,h);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         SET7(sum,0.)
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               SET4(sum,NDD(pnj,x),COEFFL(plij,h))
         SET12(NDD(pni,x),NDD(pni,x),sum,COEFFN(pni,h))
      }
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & F_LINK_TO_FACES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void sILU_vn_sf_nred(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: sILU_vn_sf_nred not available.\n");  }

void ssolveILU_vn_sf_nred(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{  eprintf("Error: ssolveILU_vn_sf_nred not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void scopy_mat_to_mat_vn_sf_nfred(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   FACE *pface;
   
   nn_copy_mat_to_mat(theGrid,a,h);
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
}

void sILU_vn_sf_nfred(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   scopy_mat_to_mat_vn_sf_nfred(theGrid,a,h);
   nn_ILU(theGrid,a,h);
}   
 
/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void ssolveILU_vn_sf_nfred(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{    
   FACE *pfi;
   
   nn_solveILU(theGrid,h,x,y);
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ)
      FD(pfi,x) = FD(pfi,y)/COEFF_FF(pfi,h);
}

#else  /*  if !((N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void sILU_vn_sf_nfred(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: sILU_vn_sf_nfred not available.\n");  }

void ssolveILU_vn_sf_nfred(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{  eprintf("Error: ssolveILU_vn_sf_nfred not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void vcopy_mat_to_mat_vn_sf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pflink;
    
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      MSET1(COEFFNNP(pnode,h),COEFFNNP(pnode,a))
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         MSET1(COEFFLLP(plink,h),COEFFLLP(plink,a))
      for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl=pnfl->next)
         SET1(COEFF_NFP(pnfl,h),COEFF_NFP(pnfl,a))
   }
    
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         COEFF_FL(pflink,h) = COEFF_FL(pflink,a);
      for (pfnl = FNSTART(pface); pfnl != NULL; pfnl=pfnl->next)
         SET1(COEFF_FNP(pfnl,h),COEFF_FNP(pfnl,a))
   }
}   

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
                (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) &&
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void vcopy_mat_to_mat_vn_sf(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: vcopy_mat_to_mat_vn_sf not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DIM == 3)

void vILU_vn_sf_3d(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   NFLINK *pnfri,*pnfrj, *pnfij;
   FNLINK *pfnir, *pfnij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot0, pivot1, pivot2, pivot;
   
   vcopy_mat_to_mat_vn_sf(theGrid,a,h);

   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      pivot1 = COEFFNN(pnr,h,1,0) /= COEFFNN(pnr,h,0,0);
      pivot2 = COEFFNN(pnr,h,2,0) /= COEFFNN(pnr,h,0,0);
      COEFFNN(pnr,h,1,1) -= pivot1*COEFFNN(pnr,h,0,1);
      COEFFNN(pnr,h,1,2) -= pivot1*COEFFNN(pnr,h,0,2);
      COEFFNN(pnr,h,2,1) -= pivot2*COEFFNN(pnr,h,0,1);
      COEFFNN(pnr,h,2,2) -= pivot2*COEFFNN(pnr,h,0,2);
      for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){   
         if ( INDEX(NBNODE(plrj)) > INDEX(pnr) ){
            COEFFLL(plrj,h,1,0) -= pivot1*COEFFLL(plrj,h,0,0);
            COEFFLL(plrj,h,2,0) -= pivot2*COEFFLL(plrj,h,0,0);
         }
         COEFFLL(plrj,h,1,1) -= pivot1*COEFFLL(plrj,h,0,1);
         COEFFLL(plrj,h,1,2) -= pivot1*COEFFLL(plrj,h,0,2);
         COEFFLL(plrj,h,2,1) -= pivot2*COEFFLL(plrj,h,0,1);
         COEFFLL(plrj,h,2,2) -= pivot2*COEFFLL(plrj,h,0,2);
      }
      for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
         COEFF_NF(pnfrj,h,1) -= pivot1*COEFF_NF(pnfrj,h,0);
         COEFF_NF(pnfrj,h,2) -= pivot2*COEFF_NF(pnfrj,h,0);   
      }
      for (plri = START(pnr); plri != NULL; plri = plri->next){
         for (plir = START(pni=NBNODE(plri)); NBNODE(plir) != pnr; 
                                                                 plir = plir->next);
         pivot1 = COEFFLL(plir,h,1,0) /= COEFFNN(pnr,h,0,0);
         pivot2 = COEFFLL(plir,h,2,0) /= COEFFNN(pnr,h,0,0);
         if ( INDEX(pni) > INDEX(pnr) ){
            pivot0 = COEFFLL(plir,h,0,0) /= COEFFNN(pnr,h,0,0);
            COEFFNN(pni,h,0,0) -= pivot0*COEFFLL(plri,h,0,0);
            COEFFNN(pni,h,0,1) -= pivot0*COEFFLL(plri,h,0,1);
            COEFFNN(pni,h,0,2) -= pivot0*COEFFLL(plri,h,0,2);
            COEFFNN(pni,h,1,0) -= pivot1*COEFFLL(plri,h,0,0);
            COEFFNN(pni,h,2,0) -= pivot2*COEFFLL(plri,h,0,0);
            COEFFLL(plir,h,0,1) -= pivot0*COEFFNN(pnr,h,0,1);
            COEFFLL(plir,h,0,2) -= pivot0*COEFFNN(pnr,h,0,2);
         }
         COEFFNN(pni,h,1,1) -= pivot1*COEFFLL(plri,h,0,1);
         COEFFNN(pni,h,1,2) -= pivot1*COEFFLL(plri,h,0,2);
         COEFFNN(pni,h,2,1) -= pivot2*COEFFLL(plri,h,0,1);
         COEFFNN(pni,h,2,2) -= pivot2*COEFFLL(plri,h,0,2);
         COEFFLL(plir,h,1,1) -= pivot1*COEFFNN(pnr,h,0,1);
         COEFFLL(plir,h,1,2) -= pivot1*COEFFNN(pnr,h,0,2);
         COEFFLL(plir,h,2,1) -= pivot2*COEFFNN(pnr,h,0,1);
         COEFFLL(plir,h,2,2) -= pivot2*COEFFNN(pnr,h,0,2);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)   
            if ((pnj=NBNODE(plrj)) != pni){
               for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj;
                                                                 plij = plij->next);
               if (plij != NULL){
                  if ( INDEX(pni) > INDEX(pnr) ){
                     if ( INDEX(pnj) > INDEX(pnr) )
                        COEFFLL(plij,h,0,0) -= pivot0*COEFFLL(plrj,h,0,0);
                     COEFFLL(plij,h,0,1) -= pivot0*COEFFLL(plrj,h,0,1);
                     COEFFLL(plij,h,0,2) -= pivot0*COEFFLL(plrj,h,0,2);
                  }
                  if ( INDEX(pnj) > INDEX(pnr) ){
                     COEFFLL(plij,h,1,0) -= pivot1*COEFFLL(plrj,h,0,0);
                     COEFFLL(plij,h,2,0) -= pivot2*COEFFLL(plrj,h,0,0);
                  }     
                  COEFFLL(plij,h,1,1) -= pivot1*COEFFLL(plrj,h,0,1);   
                  COEFFLL(plij,h,1,2) -= pivot1*COEFFLL(plrj,h,0,2);   
                  COEFFLL(plij,h,2,1) -= pivot2*COEFFLL(plrj,h,0,1);   
                  COEFFLL(plij,h,2,2) -= pivot2*COEFFLL(plrj,h,0,2);   
               }
            }
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
            pfj = NBFACE(pnfrj);
            for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) != pfj;
                                                               pnfij = pnfij->next);
            if (pnfij != NULL){
               if ( INDEX(pni) > INDEX(pnr) )
                  COEFF_NF(pnfij,h,0) -= pivot0*COEFF_NF(pnfrj,h,0);
               COEFF_NF(pnfij,h,1) -= pivot1*COEFF_NF(pnfrj,h,0);
               COEFF_NF(pnfij,h,2) -= pivot2*COEFF_NF(pnfrj,h,0);
            }
         }
      }   
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         for (pfnir = FNSTART(pfi=NBFACE(pnfri)); NBNODE(pfnir) !=  pnr; 
                                                               pfnir = pfnir->next);
         pivot = COEFF_FN(pfnir,h,0) /= COEFFNN(pnr,h,0,0);
         COEFF_FN(pfnir,h,1) -= pivot*COEFFNN(pnr,h,0,1); 
         COEFF_FN(pfnir,h,2) -= pivot*COEFFNN(pnr,h,0,2); 
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){
            pnj=NBNODE(plrj);
            for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj; 
                                                               pfnij = pfnij->next);
            if (pfnij != NULL){
               if ( INDEX(pnj) > INDEX(pnr) )
                  COEFF_FN(pfnij,h,0) -= pivot*COEFFLL(plrj,h,0,0);
               COEFF_FN(pfnij,h,1) -= pivot*COEFFLL(plrj,h,0,1); 
               COEFF_FN(pfnij,h,2) -= pivot*COEFFLL(plrj,h,0,2); 
            }
         }
         COEFF_FF(pfi,h) -= pivot*COEFF_NF(pnfri,h,0); 
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                               pflij = pflij->next);
               if (pflij != NULL)
                  COEFF_FL(pflij,h) -= pivot*COEFF_NF(pnfrj,h,0);         
            }
      }
   }
      
   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      pivot2 = COEFFNN(pnr,h,2,1) /= COEFFNN(pnr,h,1,1);
      COEFFNN(pnr,h,2,2) -= pivot2*COEFFNN(pnr,h,1,2);
      for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){   
         if ( INDEX(NBNODE(plrj)) > INDEX(pnr) )
            COEFFLL(plrj,h,2,1) -= pivot2*COEFFLL(plrj,h,1,1);
         COEFFLL(plrj,h,2,2) -= pivot2*COEFFLL(plrj,h,1,2);
      }
      for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
         COEFF_NF(pnfrj,h,2) -= pivot2*COEFF_NF(pnfrj,h,1);   
      for (plri = START(pnr); plri != NULL; plri = plri->next){
         for (plir = START(pni=NBNODE(plri)); NBNODE(plir) != pnr; 
                                                                 plir = plir->next);
         pivot2 = COEFFLL(plir,h,2,1) /= COEFFNN(pnr,h,1,1);
         if ( INDEX(pni) > INDEX(pnr) ){
            pivot1 = COEFFLL(plir,h,1,1) /= COEFFNN(pnr,h,1,1);
            COEFFNN(pni,h,1,1) -= pivot1*COEFFLL(plri,h,1,1);
            COEFFNN(pni,h,1,2) -= pivot1*COEFFLL(plri,h,1,2);
            COEFFNN(pni,h,2,1) -= pivot2*COEFFLL(plri,h,1,1);
            COEFFLL(plir,h,1,2) -= pivot1*COEFFNN(pnr,h,1,2);
         }
         COEFFNN(pni,h,2,2) -= pivot2*COEFFLL(plri,h,1,2);
         COEFFLL(plir,h,2,2) -= pivot2*COEFFNN(pnr,h,1,2);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)   
            if ((pnj=NBNODE(plrj)) != pni){
               for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj;
                                                                 plij = plij->next);
               if (plij != NULL){
                  if ( INDEX(pni) > INDEX(pnr) ){
                     if ( INDEX(pnj) > INDEX(pnr) )
                        COEFFLL(plij,h,1,1) -= pivot1*COEFFLL(plrj,h,1,1);
                     COEFFLL(plij,h,1,2) -= pivot1*COEFFLL(plrj,h,1,2);
                  }
                  if ( INDEX(pnj) > INDEX(pnr) )
                     COEFFLL(plij,h,2,1) -= pivot2*COEFFLL(plrj,h,1,1);   
                  COEFFLL(plij,h,2,2) -= pivot2*COEFFLL(plrj,h,1,2);   
               }
            }
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
            pfj = NBFACE(pnfrj);
            for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) != pfj;
                                                               pnfij = pnfij->next);
            if (pnfij != NULL){
               if ( INDEX(pni) > INDEX(pnr) )
                  COEFF_NF(pnfij,h,1) -= pivot1*COEFF_NF(pnfrj,h,1);
               COEFF_NF(pnfij,h,2) -= pivot2*COEFF_NF(pnfrj,h,1);
            }
         }
      }   
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         for (pfnir = FNSTART(pfi=NBFACE(pnfri)); NBNODE(pfnir) != pnr; 
                                                               pfnir = pfnir->next);
         pivot = COEFF_FN(pfnir,h,1) /= COEFFNN(pnr,h,1,1);
         COEFF_FN(pfnir,h,2) -= pivot*COEFFNN(pnr,h,1,2); 
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){
            pnj=NBNODE(plrj);
            for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj; 
                                                               pfnij = pfnij->next);
            if (pfnij != NULL){
               if ( INDEX(pnj) > INDEX(pnr) )
                  COEFF_FN(pfnij,h,1) -= pivot*COEFFLL(plrj,h,1,1); 
               COEFF_FN(pfnij,h,2) -= pivot*COEFFLL(plrj,h,1,2); 
            }
         }
         COEFF_FF(pfi,h) -= pivot*COEFF_NF(pnfri,h,1);
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                               pflij = pflij->next);
               if (pflij != NULL)
                  COEFF_FL(pflij,h) -= pivot*COEFF_NF(pnfrj,h,1);         
            }
      }
   }
      
   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni=NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot2 = COEFFLL(plir,h,2,2) /= COEFFNN(pnr,h,2,2);
            COEFFNN(pni,h,2,2) -= pivot2*COEFFLL(plri,h,2,2);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)   
               if (INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj;
                                                                 plij = plij->next);
                  if (plij != NULL)
                     COEFFLL(plij,h,2,2) -= pivot2*COEFFLL(plrj,h,2,2);
               }
            for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
               pfj = NBFACE(pnfrj);
               for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) != pfj;
                                                               pnfij = pnfij->next);
               if (pnfij != NULL)
                  COEFF_NF(pnfij,h,2) -= pivot2*COEFF_NF(pnfrj,h,2);
            }
         }   
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         pfi = NBFACE(pnfri);
         for (pfnir = FNSTART(pfi); NBNODE(pfnir) !=  pnr; pfnir = pfnir->next);
         pivot = COEFF_FN(pfnir,h,2) /= COEFFNN(pnr,h,2,2);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
            if (INDEX(pnj=NBNODE(plrj)) > INDEX(pnr)){
               for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj;
                                                               pfnij = pfnij->next);
               if (pfnij != NULL)
                  COEFF_FN(pfnij,h,2) -= pivot*COEFFLL(plrj,h,2,2);
            }
         COEFF_FF(pfi,h) -= pivot*COEFF_NF(pnfri,h,2); 
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                               pflij = pflij->next);
               if (pflij != NULL)
                  COEFF_FL(pflij,h) -= pivot*COEFF_NF(pnfrj,h,2);         
            }
      }
   }
             
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if (INDEX(pfi=NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFF_FL(pflir,h) /= COEFF_FF(pfr,h);
           COEFF_FF(pfi,h) -= pivot*COEFF_FL(pflri,h);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
              if (INDEX(pfj=NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                 for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                               pflij = pflij->next);
                 if (pflij != NULL)
                    COEFF_FL(pflij,h) -= pivot*COEFF_FL(pflrj,h);
              }
         }         
}
 
/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void vsolveILU_vn_sf_3d(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum;
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = 0.0;
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFLL(plij,h,0,0)*ND(pnj,x,0);
      ND(pni,x,0) = ND(pni,y,0) - sum;
   }
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = COEFFNN(pni,h,1,0)*ND(pni,x,0);
      for (plij = START(pni); plij != NULL; plij = plij->next){
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFLL(plij,h,1,1)*ND(pnj,x,1);
         sum += COEFFLL(plij,h,1,0)*ND(pnj,x,0);
      }
      ND(pni,x,1) = ND(pni,y,1) - sum;
   }
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = COEFFNN(pni,h,2,0)*ND(pni,x,0) + COEFFNN(pni,h,2,1)*ND(pni,x,1);
      for (plij = START(pni); plij != NULL; plij = plij->next){
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFLL(plij,h,2,2)*ND(pnj,x,2);
         sum += COEFFLL(plij,h,2,0)*ND(pnj,x,0) + 
                COEFFLL(plij,h,2,1)*ND(pnj,x,1);
      }
      ND(pni,x,2) = ND(pni,y,2) - sum;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sum = 0.0;
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         sum += COEFF_FN(pfnij,h,0)*ND(NBNODE(pfnij),x,0) + 
                COEFF_FN(pfnij,h,1)*ND(NBNODE(pfnij),x,1) + 
                COEFF_FN(pfnij,h,2)*ND(NBNODE(pfnij),x,2); 
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sum += COEFF_FL(pflij,h)*FD(pfj,x);
      FD(pfi,x) = FD(pfi,y) - sum;
   }    
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sum = 0.0;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sum += COEFF_FL(pflij,h)*FD(pfj,x);
         FD(pfi,x) = (FD(pfi,x) - sum)/COEFF_FF(pfi,h);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = 0.0;
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFF_NF(pnfij,h,2)*FD(NBFACE(pnfij),x);
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFLL(plij,h,2,2)*ND(pnj,x,2);
         ND(pni,x,2) = (ND(pni,x,2)-sum)/COEFFNN(pni,h,2,2);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = COEFFNN(pni,h,1,2)*ND(pni,x,2);
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFF_NF(pnfij,h,1)*FD(NBFACE(pnfij),x);
         for (plij = START(pni); plij != NULL; plij = plij->next){
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFLL(plij,h,1,1)*ND(pnj,x,1);
            sum += COEFFLL(plij,h,1,2)*ND(pnj,x,2);
         }
         ND(pni,x,1) = (ND(pni,x,1)-sum)/COEFFNN(pni,h,1,1);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = COEFFNN(pni,h,0,1)*ND(pni,x,1) + COEFFNN(pni,h,0,2)*ND(pni,x,2);
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFF_NF(pnfij,h,0)*FD(NBFACE(pnfij),x);
         for (plij = START(pni); plij != NULL; plij = plij->next){
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFLL(plij,h,0,0)*ND(pnj,x,0);
            sum += COEFFLL(plij,h,0,1)*ND(pnj,x,1) + 
                   COEFFLL(plij,h,0,2)*ND(pnj,x,2);
         }
         ND(pni,x,0) = (ND(pni,x,0)-sum)/COEFFNN(pni,h,0,0);
      }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) &&
                (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) &&
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) &&
                (DIM == 3))  */

void vILU_vn_sf_3d(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: vILU_vn_sf_3d not available.\n");  }

void vsolveILU_vn_sf_3d(theGrid,h,x,y) 
GRID *theGrid; INT h,x,y;
{  eprintf("Error: vsolveILU_vn_sf_3d not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vcopy_mat_to_mat_vn_vf(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pflink;
    
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      MSET1(COEFFNNP(pnode,h),COEFFNNP(pnode,a))
      for (plink = START(pnode); plink != NULL; plink=plink->next)
         MSET1(COEFFLLP(plink,h),COEFFLLP(plink,a))
      for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl=pnfl->next)
         MSET1(COEFFLLP(pnfl,h),COEFFLLP(pnfl,a))
   }
    
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      MSET1(COEFFNNP(pface,h),COEFFNNP(pface,a))
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         MSET1(COEFFLLP(pflink,h),COEFFLLP(pflink,a))
      for (pfnl = FNSTART(pface); pfnl != NULL; pfnl=pfnl->next)
         MSET1(COEFFLLP(pfnl,h),COEFFLLP(pfnl,a))
   }
}   

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) &&
                (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) &&
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

void vcopy_mat_to_mat_vn_vf(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: vcopy_mat_to_mat_vn_vf not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (DIM == 2)

void vILU_vn_vf_2d(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   NFLINK *pnfri,*pnfrj, *pnfij;
   FNLINK *pfnir, *pfnij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot0, pivot1;
   
   vcopy_mat_to_mat_vn_vf(theGrid,a,h);

   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      pivot1 = COEFFNN(pnr,h,1,0) /= COEFFNN(pnr,h,0,0);
      COEFFNN(pnr,h,1,1) -= pivot1*COEFFNN(pnr,h,0,1);
      for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){   
         if ( INDEX(NBNODE(plrj)) > INDEX(pnr) )
            COEFFLL(plrj,h,1,0) -= pivot1*COEFFLL(plrj,h,0,0);
         COEFFLL(plrj,h,1,1) -= pivot1*COEFFLL(plrj,h,0,1);
      }
      for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
         COEFFLL(pnfrj,h,1,0) -= pivot1*COEFFLL(pnfrj,h,0,0);
         COEFFLL(pnfrj,h,1,1) -= pivot1*COEFFLL(pnfrj,h,0,1);
      }
      for (plri = START(pnr); plri != NULL; plri = plri->next){
         for (plir = START(pni=NBNODE(plri)); NBNODE(plir) != pnr; 
                                                             plir = plir->next);
         pivot1 = COEFFLL(plir,h,1,0) /= COEFFNN(pnr,h,0,0);
         if ( INDEX(pni) > INDEX(pnr) ){
            pivot0 = COEFFLL(plir,h,0,0) /= COEFFNN(pnr,h,0,0);
            COEFFNN(pni,h,0,0) -= pivot0*COEFFLL(plri,h,0,0);
            COEFFNN(pni,h,0,1) -= pivot0*COEFFLL(plri,h,0,1);
            COEFFNN(pni,h,1,0) -= pivot1*COEFFLL(plri,h,0,0);
            COEFFLL(plir,h,0,1) -= pivot0*COEFFNN(pnr,h,0,1);
         }
         COEFFNN(pni,h,1,1) -= pivot1*COEFFLL(plri,h,0,1);
         COEFFLL(plir,h,1,1) -= pivot1*COEFFNN(pnr,h,0,1);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)   
            if ((pnj=NBNODE(plrj)) != pni){
               for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj;
                                                             plij = plij->next);
               if (plij != NULL){
                  if ( INDEX(pni) > INDEX(pnr) ){
                     if ( INDEX(pnj) > INDEX(pnr) )
                        COEFFLL(plij,h,0,0) -= pivot0*COEFFLL(plrj,h,0,0);
                     COEFFLL(plij,h,0,1) -= pivot0*COEFFLL(plrj,h,0,1);
                  }
                  if ( INDEX(pnj) > INDEX(pnr) )
                     COEFFLL(plij,h,1,0) -= pivot1*COEFFLL(plrj,h,0,0);   
                  COEFFLL(plij,h,1,1) -= pivot1*COEFFLL(plrj,h,0,1);   
               }
            }
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
            pfj = NBFACE(pnfrj);
            for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) != pfj;
                                                           pnfij = pnfij->next);
            if (pnfij != NULL){
               if ( INDEX(pni) > INDEX(pnr) ){
                  COEFFLL(pnfij,h,0,0) -= pivot0*COEFFLL(pnfrj,h,0,0);
                  COEFFLL(pnfij,h,0,1) -= pivot0*COEFFLL(pnfrj,h,0,1);
               }
               COEFFLL(pnfij,h,1,0) -= pivot1*COEFFLL(pnfrj,h,0,0);
               COEFFLL(pnfij,h,1,1) -= pivot1*COEFFLL(pnfrj,h,0,1);
            }
         }
      }   
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         for (pfnir = FNSTART(pfi=NBFACE(pnfri)); NBNODE(pfnir) != pnr; 
                                                           pfnir = pfnir->next);
         pivot0 = COEFFLL(pfnir,h,0,0) /= COEFFNN(pnr,h,0,0);
         pivot1 = COEFFLL(pfnir,h,1,0) /= COEFFNN(pnr,h,0,0);
         COEFFLL(pfnir,h,0,1) -= pivot0*COEFFNN(pnr,h,0,1); 
         COEFFLL(pfnir,h,1,1) -= pivot1*COEFFNN(pnr,h,0,1); 
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next){
            pnj=NBNODE(plrj);
            for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj; 
                                                           pfnij = pfnij->next);
            if (pfnij != NULL){
               if ( INDEX(pnj) > INDEX(pnr) ){
                  COEFFLL(pfnij,h,0,0) -= pivot0*COEFFLL(plrj,h,0,0); 
                  COEFFLL(pfnij,h,1,0) -= pivot1*COEFFLL(plrj,h,0,0); 
               }
               COEFFLL(pfnij,h,0,1) -= pivot0*COEFFLL(plrj,h,0,1); 
               COEFFLL(pfnij,h,1,1) -= pivot1*COEFFLL(plrj,h,0,1); 
            }
         }
         COEFFNN(pfi,h,0,0) -= pivot0*COEFFLL(pnfri,h,0,0);
         COEFFNN(pfi,h,0,1) -= pivot0*COEFFLL(pnfri,h,0,1);
         COEFFNN(pfi,h,1,0) -= pivot1*COEFFLL(pnfri,h,0,0);
         COEFFNN(pfi,h,1,1) -= pivot1*COEFFLL(pnfri,h,0,1);
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
               if (pflij != NULL){
                  COEFFLL(pflij,h,0,0) -= pivot0*COEFFLL(pnfrj,h,0,0);         
                  COEFFLL(pflij,h,0,1) -= pivot0*COEFFLL(pnfrj,h,0,1);         
                  COEFFLL(pflij,h,1,0) -= pivot1*COEFFLL(pnfrj,h,0,0);         
                  COEFFLL(pflij,h,1,1) -= pivot1*COEFFLL(pnfrj,h,0,1);         
               }
            }
      }
   }
 
   for (pnr = FIRSTNODE(theGrid); pnr != NULL; pnr = pnr->succ){ 
      for (plri = START(pnr); plri != NULL; plri = plri->next)
         if ( INDEX(pni=NBNODE(plri)) > INDEX(pnr) ){
            for (plir = START(pni); NBNODE(plir) != pnr; plir = plir->next);
            pivot1 = COEFFLL(plir,h,1,1) /= COEFFNN(pnr,h,1,1);
            COEFFNN(pni,h,1,1) -= pivot1*COEFFLL(plri,h,1,1);
            for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)   
               if (INDEX(pnj=NBNODE(plrj)) > INDEX(pnr) && pnj != pni){
                  for (plij = START(pni); plij != NULL && NBNODE(plij) != pnj;
                                                             plij = plij->next);
                  if (plij != NULL)
                     COEFFLL(plij,h,1,1) -= pivot1*COEFFLL(plrj,h,1,1);
               }
            for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next){
               pfj = NBFACE(pnfrj);
               for (pnfij = NFSTART(pni); pnfij != NULL && NBFACE(pnfij) != pfj;
                                                           pnfij = pnfij->next);
               if (pnfij != NULL){
                  COEFFLL(pnfij,h,1,0) -= pivot1*COEFFLL(pnfrj,h,1,0);
                  COEFFLL(pnfij,h,1,1) -= pivot1*COEFFLL(pnfrj,h,1,1);
               }
            }
         }   
      for (pnfri = NFSTART(pnr); pnfri != NULL; pnfri=pnfri->next){
         pfi = NBFACE(pnfri);
         for (pfnir = FNSTART(pfi); NBNODE(pfnir) !=  pnr; pfnir = pfnir->next);
         pivot0 = COEFFLL(pfnir,h,0,1) /= COEFFNN(pnr,h,1,1);
         pivot1 = COEFFLL(pfnir,h,1,1) /= COEFFNN(pnr,h,1,1);
         for (plrj = START(pnr); plrj != NULL; plrj=plrj->next)
            if (INDEX(pnj=NBNODE(plrj)) > INDEX(pnr)){
               for (pfnij = FNSTART(pfi); pfnij != NULL && NBNODE(pfnij) != pnj;
                                                           pfnij = pfnij->next);
               if (pfnij != NULL){
                  COEFFLL(pfnij,h,0,1) -= pivot0*COEFFLL(plrj,h,1,1);
                  COEFFLL(pfnij,h,1,1) -= pivot1*COEFFLL(plrj,h,1,1);
               }
            }
         COEFFNN(pfi,h,0,0) -= pivot0*COEFFLL(pnfri,h,1,0); 
         COEFFNN(pfi,h,0,1) -= pivot0*COEFFLL(pnfri,h,1,1); 
         COEFFNN(pfi,h,1,0) -= pivot1*COEFFLL(pnfri,h,1,0); 
         COEFFNN(pfi,h,1,1) -= pivot1*COEFFLL(pnfri,h,1,1); 
         for (pnfrj = NFSTART(pnr); pnfrj != NULL; pnfrj = pnfrj->next)
            if ( (pfj = NBFACE(pnfrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
               if (pflij != NULL){
                  COEFFLL(pflij,h,0,0) -= pivot0*COEFFLL(pnfrj,h,1,0);         
                  COEFFLL(pflij,h,0,1) -= pivot0*COEFFLL(pnfrj,h,1,1);         
                  COEFFLL(pflij,h,1,0) -= pivot1*COEFFLL(pnfrj,h,1,0);         
                  COEFFLL(pflij,h,1,1) -= pivot1*COEFFLL(pnfrj,h,1,1);         
               }
            }
      }
   }

   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ){
      pivot1 = COEFFNN(pfr,h,1,0) /= COEFFNN(pfr,h,0,0);
      COEFFNN(pfr,h,1,1) -= pivot1*COEFFNN(pfr,h,0,1);
      for (pflrj = FSTART(pfr); pflrj != NULL; pflrj=pflrj->next){   
         if ( INDEX(NBFACE(pflrj)) > INDEX(pfr) )
            COEFFLL(pflrj,h,1,0) -= pivot1*COEFFLL(pflrj,h,0,0);
         COEFFLL(pflrj,h,1,1) -= pivot1*COEFFLL(pflrj,h,0,1);
      }
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next){
         for (pflir = FSTART(pfi=NBFACE(pflri)); NBFACE(pflir) != pfr; 
                                                           pflir = pflir->next);
         pivot1 = COEFFLL(pflir,h,1,0) /= COEFFNN(pfr,h,0,0);
         if ( INDEX(pfi) > INDEX(pfr) ){
            pivot0 = COEFFLL(pflir,h,0,0) /= COEFFNN(pfr,h,0,0);
            COEFFNN(pfi,h,0,0) -= pivot0*COEFFLL(pflri,h,0,0);
            COEFFNN(pfi,h,0,1) -= pivot0*COEFFLL(pflri,h,0,1);
            COEFFNN(pfi,h,1,0) -= pivot1*COEFFLL(pflri,h,0,0);
            COEFFLL(pflir,h,0,1) -= pivot0*COEFFNN(pfr,h,0,1);
         }
         COEFFNN(pfi,h,1,1) -= pivot1*COEFFLL(pflri,h,0,1);
         COEFFLL(pflir,h,1,1) -= pivot1*COEFFNN(pfr,h,0,1);
         for (pflrj = FSTART(pfr); pflrj != NULL; pflrj=pflrj->next)   
            if ((pfj=NBFACE(pflrj)) != pfi){
               for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
               if (pflij != NULL){
                  if ( INDEX(pfi) > INDEX(pfr) ){
                     if ( INDEX(pfj) > INDEX(pfr) )
                        COEFFLL(pflij,h,0,0) -= pivot0*COEFFLL(pflrj,h,0,0);
                     COEFFLL(pflij,h,0,1) -= pivot0*COEFFLL(pflrj,h,0,1);
                  }
                  if ( INDEX(pfj) > INDEX(pfr) )
                     COEFFLL(pflij,h,1,0) -= pivot1*COEFFLL(pflrj,h,0,0);   
                  COEFFLL(pflij,h,1,1) -= pivot1*COEFFLL(pflrj,h,0,1);   
               }
            }
      }   
   }
      
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ){ 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi=NBFACE(pflri)) > INDEX(pfr) ){
            for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir =pflir->next);
            pivot1 = COEFFLL(pflir,h,1,1) /= COEFFNN(pfr,h,1,1);
            COEFFNN(pfi,h,1,1) -= pivot1*COEFFLL(pflri,h,1,1);
            for (pflrj = FSTART(pfr); pflrj != NULL; pflrj=pflrj->next)   
               if (INDEX(pfj=NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                  for (pflij = FSTART(pfi); pflij != NULL && 
                                     NBFACE(pflij) != pfj; pflij = pflij->next);
                  if (pflij != NULL)
                     COEFFLL(pflij,h,1,1) -= pivot1*COEFFLL(pflrj,h,1,1);
               }
         }   
   }
}
 
/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void vsolveILU_vn_vf_2d(theGrid,h,x,y)
GRID *theGrid;
INT h,x,y;
{
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum;
   
   /* find x: L x = y */
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = 0.;
      for (plij = START(pni); plij != NULL; plij = plij->next)
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFLL(plij,h,0,0)*ND(pnj,x,0);
      ND(pni,x,0) = ND(pni,y,0) - sum;
   }
   for (pni = FIRSTNODE(theGrid); pni != NULL; pni = pni->succ){
      sum = COEFFNN(pni,h,1,0)*ND(pni,x,0);
      for (plij = START(pni); plij != NULL; plij = plij->next){
         if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
            sum += COEFFLL(plij,h,1,1)*ND(pnj,x,1);
         sum += COEFFLL(plij,h,1,0)*ND(pnj,x,0);
      }
      ND(pni,x,1) = ND(pni,y,1) - sum;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sum = 0.;
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         sum += COEFFLL(pfnij,h,0,0)*ND(NBNODE(pfnij),x,0) + 
                COEFFLL(pfnij,h,0,1)*ND(NBNODE(pfnij),x,1);
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi) )
            sum += COEFFLL(pflij,h,0,0)*FDV(pfj,x,0);
      FDV(pfi,x,0) = FDV(pfi,y,0) - sum;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sum = COEFFNN(pfi,h,1,0)*FDV(pfi,x,0);
      for (pfnij = FNSTART(pfi); pfnij != NULL; pfnij = pfnij->next)
         sum += COEFFLL(pfnij,h,1,0)*ND(NBNODE(pfnij),x,0) + 
                COEFFLL(pfnij,h,1,1)*ND(NBNODE(pfnij),x,1);
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi) )
            sum += COEFFLL(pflij,h,1,1)*FDV(pfj,x,1);
         sum += COEFFLL(pflij,h,1,0)*FDV(pfj,x,0);
      }
      FDV(pfi,x,1) = FDV(pfi,y,1) - sum;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sum = 0.;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sum += COEFFLL(pflij,h,1,1)*FDV(pfj,x,1);
         FDV(pfi,x,1) = (FDV(pfi,x,1)-sum)/COEFFNN(pfi,h,1,1);
      }
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){
         sum = COEFFNN(pfi,h,0,1)*FDV(pfi,x,1);
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sum += COEFFLL(pflij,h,0,0)*FDV(pfj,x,0);
            sum += COEFFLL(pflij,h,0,1)*FDV(pfj,x,1);
         }
         FDV(pfi,x,0) = (FDV(pfi,x,0)-sum)/COEFFNN(pfi,h,0,0);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = 0.;
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFFLL(pnfij,h,1,0)*FDV(NBFACE(pnfij),x,0) +
                   COEFFLL(pnfij,h,1,1)*FDV(NBFACE(pnfij),x,1);
         for (plij = START(pni); plij != NULL; plij = plij->next)
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFLL(plij,h,1,1)*ND(pnj,x,1);
         ND(pni,x,1) = (ND(pni,x,1)-sum)/COEFFNN(pni,h,1,1);
      }
   for (pni = theGrid->lastNode; pni >= FIRSTNODE(theGrid); pni--)
      if (IS_FN(pni)){
         sum = COEFFNN(pni,h,0,1)*ND(pni,x,1);
         for (pnfij = NFSTART(pni); pnfij != NULL; pnfij = pnfij->next)
            sum += COEFFLL(pnfij,h,0,0)*FDV(NBFACE(pnfij),x,0) +
                   COEFFLL(pnfij,h,0,1)*FDV(NBFACE(pnfij),x,1);
         for (plij = START(pni); plij != NULL; plij = plij->next){
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFFLL(plij,h,0,0)*ND(pnj,x,0);
            sum += COEFFLL(plij,h,0,1)*ND(pnj,x,1);
         }
         ND(pni,x,0) = (ND(pni,x,0)-sum)/COEFFNN(pni,h,0,0);
      }
}

#else  /*  if !((N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) &&
                (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) &&
                (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) &&
                (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) &&
                (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) &&
                (DIM == 2))  */

void vILU_vn_vf_2d(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: vILU_vn_vf_2d not available.\n");  }

void vsolveILU_vn_vf_2d(theGrid,h,x,y)
GRID *theGrid; INT h,x,y;
{  eprintf("Error: vsolveILU_vn_vf_2d not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) 

void copy_mat_to_mat_f(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   FACE *pface;
   FLINK *pflink;
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      COEFF_FF(pface,h) = COEFF_FF(pface,a);
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         COEFF_FL(pflink,h) = COEFF_FL(pflink,a);
   }
}
 
void ILU_f(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   FACE *pfr, *pfi, *pfj;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot;
   
   copy_mat_to_mat_f(theGrid,a,h);

   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFF_FL(pflir,h) /= COEFF_FF(pfr,h);
           COEFF_FF(pfi,h) -= pivot*COEFF_FL(pflri,h);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL)
                   COEFF_FL(pflij,h) -= pivot*COEFF_FL(pflrj,h);
             }
         }         
}   

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES)  */

void ILU_f(theGrid,a,h)
GRID *theGrid; INT a, h;
{  eprintf("Error: ILU_f not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA)

void ILU_f_beta_full(theGrid,a,h,q,beta)
GRID *theGrid;
INT a, h, q;
FLOAT beta;
{
   FACE *pfi, *pfj, *pfk;
   FLINK *pflij, *pflik, *pflkj;
   FLOAT s, hij;
   INT ind;
   
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      FD(pfi,q) = 0.;
      for (pfj = FIRSTFACE(theGrid); pfj != NULL; pfj = pfj->succ){
         if (pfi == pfj){
            s = -COEFF_FF(pfi,a);
            hij = COEFF_FF(pfi,h);
         }
         else{ 
            for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
            if (pflij != NULL){
               s = -COEFF_FL(pflij,a);
               hij = COEFF_FL(pflij,h);
            }
            else
               s = hij = 0.;
         }
         if (INDEX(pfi) > INDEX(pfj))
            ind = INDEX(pfj);
         else{
            ind = INDEX(pfi)-1; 
            s += hij; 
         }
         for (pflik = FSTART(pfi); pflik != NULL; pflik = pflik->next)
            if ( INDEX(pfk = NBFACE(pflik)) <= ind){
               for (pflkj = FSTART(pfk); pflkj != NULL && NBFACE(pflkj) != pfj;
                                                           pflkj = pflkj->next);
               if (pflkj != NULL)
                  s += COEFF_FL(pflik,h)*COEFF_FL(pflkj,h);
            }
         FD(pfi,q) += fabs(s);
      }
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      s = COEFF_FF(pfi,h) + FD(pfi,q)*beta;
      FD(pfi,q) = COEFF_FF(pfi,h)/s;
      COEFF_FF(pfi,h) = s;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ)
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi) )
            COEFF_FL(pflij,h) *= FD(pfj,q);
}

void ILU_f_beta_sparse(theGrid,a,h,q,beta)
GRID *theGrid;
INT a, h, q;
FLOAT beta;
{
   FACE *pfi, *pfj, *pfk;
   FLINK *pflij, *pflik, *pflkj, *pflki;
   FLOAT s;
   INT ind;
   
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      s = COEFF_FF(pfi,h)-COEFF_FF(pfi,a);
      ind = INDEX(pfi)-1; 
      for (pflik = FSTART(pfi); pflik != NULL; pflik = pflik->next)
         if ( INDEX(pfk = NBFACE(pflik)) <= ind){
            for (pflki = FSTART(pfk); NBFACE(pflki) != pfi; pflki= pflki->next);
            s += COEFF_FL(pflik,h)*COEFF_FL(pflki,h);
         }
      FD(pfi,q) = fabs(s);
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
         pfj = NBFACE(pflij); 
         s = -COEFF_FL(pflij,a);
         if (INDEX(pfi) > INDEX(pfj))
            ind = INDEX(pfj);
         else{
            ind = INDEX(pfi)-1; 
            s += COEFF_FL(pflij,h); 
         }
         for (pflik = FSTART(pfi); pflik != NULL; pflik = pflik->next)
            if ( INDEX(pfk = NBFACE(pflik)) <= ind){
               for (pflkj = FSTART(pfk); pflkj != NULL && NBFACE(pflkj) != pfj;
                                                           pflkj = pflkj->next);
               if (pflkj != NULL)
                  s += COEFF_FL(pflik,h)*COEFF_FL(pflkj,h);
            }
         FD(pfi,q) += fabs(s);
      }
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      s = COEFF_FF(pfi,h) + FD(pfi,q)*beta;
      FD(pfi,q) = COEFF_FF(pfi,h)/s;
      COEFF_FF(pfi,h) = s;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ)
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi) )
            COEFF_FL(pflij,h) *= FD(pfj,q);
}

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void solveILU_f(theGrid,h,x,y) 
GRID *theGrid;
INT h,x,y;
{
   FACE *pfi, *pfj;
   FLINK *pflij;
   FLOAT sumf;
   
   /* find x: L x = y */
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sumf = 0.0;
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sumf += COEFF_FL(pflij,h)*FD(pfj,x);
      FD(pfi,x) = FD(pfi,y) - sumf;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         sumf = 0.0;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sumf += COEFF_FL(pflij,h)*FD(pfj,x);
         FD(pfi,x) = (FD(pfi,x) - sumf)/COEFF_FF(pfi,h);
      }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) &&
                (F_DATA & SCALAR_FACE_DATA))  */

void ILU_f_beta_full(theGrid,a,h,q,beta)
GRID *theGrid; INT a, h, q; FLOAT beta;
{  eprintf("Error: ILU_f_beta_full not available.\n");  }

void ILU_f_beta_sparse(theGrid,a,h,q,beta)
GRID *theGrid; INT a, h, q; FLOAT beta;
{  eprintf("Error: ILU_f_beta_sparse not available.\n");  }

void solveILU_f(theGrid,h,x,y)
GRID *theGrid; INT h, x, y;
{  eprintf("Error: solveILU_f not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) 

void vcopy_mat_to_mat_f(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   FACE *pface;
   FLINK *pflink;
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ){
      MSET1(COEFFNNP(pface,h),COEFFNNP(pface,a))
      for (pflink = FSTART(pface); pflink != NULL; pflink=pflink->next)
         MSET1(COEFFNNP(pflink,h),COEFFNNP(pflink,a))
   }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES))  */

void vcopy_mat_to_mat_f(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: vcopy_mat_to_mat_f not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (DIM == 2)

void vILU_vf(theGrid,a,h)
GRID *theGrid;
INT a, h;
{
   FACE *pfr, *pfi, *pfj;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   FLOAT pivot;
   
   vcopy_mat_to_mat_f(theGrid,a,h);

   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ){
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFFNN(pflir,h,0,0) /= COEFFNN(pfr,h,0,0);
           COEFFNN(pflir,h,0,1) -= pivot*COEFFNN(pfr,h,0,1);
           COEFFNN(pfi,h,0,0) -= pivot*COEFFNN(pflri,h,0,0);
           COEFFNN(pfi,h,0,1) -= pivot*COEFFNN(pflri,h,0,1);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( (pfj = NBFACE(pflrj)) != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL){
                   if ( INDEX(pfj) > INDEX(pfr))
                      COEFFNN(pflij,h,0,0) -= pivot*COEFFNN(pflrj,h,0,0);
                   COEFFNN(pflij,h,0,1) -= pivot*COEFFNN(pflrj,h,0,1);
                }
             }
         }
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next){
           pfi = NBFACE(pflri);
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFFNN(pflir,h,1,0) /= COEFFNN(pfr,h,0,0);
           COEFFNN(pflir,h,1,1) -= pivot*COEFFNN(pfr,h,0,1);
           if ( INDEX(pfi) > INDEX(pfr))
              COEFFNN(pfi,h,1,0) -= pivot*COEFFNN(pflri,h,0,0);
           COEFFNN(pfi,h,1,1) -= pivot*COEFFNN(pflri,h,0,1);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( (pfj = NBFACE(pflrj)) != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL){
                   if ( INDEX(pfj) > INDEX(pfr) )
                      COEFFNN(pflij,h,1,0) -= pivot*COEFFNN(pflrj,h,0,0);
                   COEFFNN(pflij,h,1,1) -= pivot*COEFFNN(pflrj,h,0,1);
                }
             }
         }
      pivot = COEFFNN(pfr,h,1,0) /= COEFFNN(pfr,h,0,0);
      COEFFNN(pfr,h,1,1) -= pivot*COEFFNN(pfr,h,0,1);
      for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next){
         if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) )
            COEFFNN(pflrj,h,1,0)  -= pivot*COEFFNN(pflrj,h,0,0);
         COEFFNN(pflrj,h,1,1)  -= pivot*COEFFNN(pflrj,h,0,1);
      }
   }
   for (pfr = FIRSTFACE(theGrid); pfr != NULL; pfr = pfr->succ) 
      for (pflri = FSTART(pfr); pflri != NULL; pflri = pflri->next)
         if ( INDEX(pfi = NBFACE(pflri)) > INDEX(pfr)){
           for (pflir = FSTART(pfi); NBFACE(pflir) != pfr; pflir = pflir->next);
           pivot = COEFFNN(pflir,h,1,1) /= COEFFNN(pfr,h,1,1);
           COEFFNN(pfi,h,1,1) -= pivot*COEFFNN(pflri,h,1,1);
           for (pflrj = FSTART(pfr); pflrj != NULL; pflrj = pflrj->next)
             if ( INDEX(pfj = NBFACE(pflrj)) > INDEX(pfr) && pfj != pfi){
                for (pflij = FSTART(pfi); pflij != NULL && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                if (pflij != NULL)
                   COEFFNN(pflij,h,1,1) -= pivot*COEFFNN(pflrj,h,1,1);
             }
         }         
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && 
                (DIM==2))  */

void vILU_vf(theGrid,a,h)
GRID *theGrid; INT a, h;
{  eprintf("Error: vILU_vf not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA)

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */
void ssolveILU_vf(theGrid,h,x,y)
GRID *theGrid;
INT h,x,y;
{
   FACE *pfi, *pfj;
   FLINK *pflij;
   FLOAT sumf[DIM];
   
   /* find x: L x = y */
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      SET7(sumf,0.)
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            SET4(sumf,FDVP(pfj,x),COEFF_FL(pflij,h))
      SET11(FDVP(pfi,x),FDVP(pfi,y),sumf)
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         SET7(sumf,0.)
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               SET4(sumf,FDVP(pfj,x),COEFF_FL(pflij,h))
         SET12(FDVP(pfi,x),FDVP(pfi,x),sumf,COEFF_FF(pfi,h))
      }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) &&
                (F_DATA & VECTOR_FACE_DATA))  */

void ssolveILU_vf(theGrid,h,x,y)
GRID *theGrid; INT h,x,y;
{  eprintf("Error: ssolveILU_vf not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA) && (DIM==2)

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void vsolveILU_vf_2d(theGrid,h,x,y)
GRID *theGrid;
INT h, x, y;
{
   FACE *pfi, *pfj;
   FLINK *pflij;
   FLOAT sumf;
   
   /* find x: L x = y */
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sumf = 0.0;
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sumf += COEFFNN(pflij,h,0,0)*FDV(pfj,x,0);
      FDV(pfi,x,0) = FDV(pfi,y,0) - sumf;
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      sumf = COEFFNN(pfi,h,1,0)*FDV(pfi,x,0);
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            sumf += COEFFNN(pflij,h,1,1)*FDV(pfj,x,1);
         sumf += COEFFNN(pflij,h,1,0)*FDV(pfj,x,0);
      }
      FDV(pfi,x,1) = FDV(pfi,y,1) - sumf;
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         sumf = 0.0;
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sumf += COEFFNN(pflij,h,1,1)*FDV(pfj,x,1);
         FDV(pfi,x,1) = (FDV(pfi,x,1) - sumf)/COEFFNN(pfi,h,1,1);
      }
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         sumf = COEFFNN(pfi,h,0,1)*FDV(pfi,x,1);
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sumf += COEFFNN(pflij,h,0,0)*FDV(pfj,x,0);
            sumf += COEFFNN(pflij,h,0,1)*FDV(pfj,x,1);
         }
         FDV(pfi,x,0) = (FDV(pfi,x,0) - sumf)/COEFFNN(pfi,h,0,0);
      }
}

#else  /*  if !((F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) &&
                (F_DATA & VECTOR_FACE_DATA) && (DIM==2))  */

void vsolveILU_vf_2d(theGrid,h,x,y)
GRID *theGrid; INT h, x, y;
{  eprintf("Error: vsolveILU_vf_2d not available.\n");  }

#endif 

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & DVECTOR_FACE_DATA)

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void dvsolveILU_f(theGrid,h,x,y)
GRID *theGrid;
INT h, x, y;
{
   FACE *pfi, *pfj;
   FLINK *pflij;
   FLOAT sumf[DIM];
   
   /* find x: L x = y */
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      SET7(sumf,0.)
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            SET4(sumf,FDDVP(pfj,x,0),COEFFNN(pflij,h,0,0))
      SET11(FDDVP(pfi,x,0),FDDVP(pfi,y,0),sumf)
   }
   for (pfi = FIRSTFACE(theGrid); pfi != NULL; pfi = pfi->succ){
      SET2(sumf,FDDVP(pfi,x,0),COEFFNN(pfi,h,1,0))
      for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
         if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi))
            SET4(sumf,FDDVP(pfj,x,1),COEFFNN(pflij,h,1,1))
         SET4(sumf,FDDVP(pfj,x,0),COEFFNN(pflij,h,1,0))
      }
      SET11(FDDVP(pfi,x,1),FDDVP(pfi,y,1),sumf)
   }
   
   /* find solution of U . = x and store it in x */
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         SET7(sumf,0.)
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next)
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               SET4(sumf,FDDVP(pfj,x,1),COEFFNN(pflij,h,1,1))
         SET12(FDDVP(pfi,x,1),FDDVP(pfi,x,1),sumf,COEFFNN(pfi,h,1,1))
      }
   for (pfi = theGrid->lastFace; pfi >= FIRSTFACE(theGrid); pfi--)
      if (IS_FF(pfi) && INDEX(pfi) > 0){  
         SET2(sumf,FDDVP(pfi,x,1),COEFFNN(pfi,h,0,1))
         for (pflij = FSTART(pfi); pflij != NULL; pflij = pflij->next){
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               SET4(sumf,FDDVP(pfj,x,0),COEFFNN(pflij,h,0,0))
            SET4(sumf,FDDVP(pfj,x,1),COEFFNN(pflij,h,0,1))
         }
         SET12(FDDVP(pfi,x,0),FDDVP(pfi,x,0),sumf,COEFFNN(pfi,h,0,0))
      }
}

#else  /*  if !((F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) &&
                (F_DATA & VECTOR_FACE_DATA))  */

void dvsolveILU_f(theGrid,h,x,y)
GRID *theGrid; INT h, x, y;
{  eprintf("Error: dvsolveILU_f not available.\n");  }

#endif

#if (N_DATA & NxN_NODE_MATR) && (N_DATA & NxM_NODE_FACE_MATR) && (F_DATA & NxN_FACE_MATR) && (F_DATA & MxN_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (DATA_S & N_LINK_TO_ELEMENTS) && (DATA_S & PREVIOUS_NODE) && (DATA_S & PREVIOUS_FACE) && (E_DATA & MxM_E_E_MATR) && (E_DATA & MxN_E_N_MATR) && (E_DATA & NxM_N_E_MATR) && (E_DATA & MxN_E_F_MATR) && (E_DATA & NxM_F_E_MATR) && (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) 

void copy_mat_to_mat_general(theGrid,a,h)
GRID *theGrid;
INT a, h;
{
   NODE *theNode;
   FACE *theFace;
   ELEMENT *pel;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   INT i, j, k, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC, mm=N_OF_ELEM_FUNC;

   for (theNode=FIRSTNODE(theGrid); theNode; theNode=SUCC(theNode)){
      MMSET1(COEFF_NNP(theNode,h),COEFF_NNP(theNode,a),i,j,kk,kk)
      for (pl=TSTART(theNode); pl; pl=NEXT(pl))
         MMSET1(COEFF_NNP(pl,h),COEFF_NNP(pl,a),i,j,kk,kk)
      for (pnfl=TNFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
         MMSET1(COEFF_NNP(pnfl,h),COEFF_NNP(pnfl,a),i,j,kk,nn)
   }
   for (theFace=FIRSTFACE(theGrid); theFace; theFace=SUCC(theFace)){
      MMSET1(COEFF_NNP(theFace,h),COEFF_NNP(theFace,a),i,j,nn,nn)
      for (pfl=TFSTART(theFace); pfl; pfl=NEXT(pfl))
         MMSET1(COEFF_NNP(pfl,h),COEFF_NNP(pfl,a),i,j,nn,nn)
      for (pfnl=TFNSTART(theFace); pfnl; pfnl=NEXT(pfnl))
         MMSET1(COEFF_NNP(pfnl,h),COEFF_NNP(pfnl,a),i,j,nn,kk)
   }
   for (pel = FIRSTELEMENT(theGrid); pel; pel=pel->succ){
      MMSET1(COEFF_EE_MMP(pel,h),COEFF_EE_MMP(pel,a),i,j,mm,mm)
      for (k=0; k < NVERT; k++){
         MMSET1(COEFF_NE_NMP(pel,h,k),COEFF_NE_NMP(pel,a,k),i,j,kk,mm)
         MMSET1(COEFF_FE_NMP(pel,h,k),COEFF_FE_NMP(pel,a,k),i,j,nn,mm)
         MMSET1(COEFF_EF_MNP(pel,h,k),COEFF_EF_MNP(pel,a,k),i,j,mm,nn)
         MMSET1(COEFF_EN_MNP(pel,h,k),COEFF_EN_MNP(pel,a,k),i,j,mm,kk)
      }
   }
}

void ILU_general(theGrid,a,h)
GRID *theGrid;
INT a,h;
{
   ELEMENT *pel;
   NODE *pnr, *pni, *pnj;
   FACE *pfr, *pfi, *pfj;
   LINK *plri, *plir, *plrj, *plij;
   NFLINK *pnfri,*pnfrj, *pnfij;
   FNLINK *pfnir, *pfnij;
   FLINK *pflri, *pflir, *pflrj, *pflij;
   NELINK *pneri, *pnerj;
   FLOAT pivot[50], piv;
   INT i, j, k, k1, p, q, r, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC,
                                                mm=N_OF_ELEM_FUNC;
   
   copy_mat_to_mat_general(theGrid,a,h);

   for (k = 0; k < kk; k++){
      k1 = k + 1;
      for (pnr = FIRSTNODE(theGrid); pnr; pnr = pnr->succ){
         if (k1 < kk){
            for (p = k1; p < kk; p++){
               pivot[p] = COEFF_NN(pnr,h,p,k) /= COEFF_NN(pnr,h,k,k);
               for (q = k1; q < kk; q++)
                  COEFF_NN(pnr,h,p,q) -= pivot[p]*COEFF_NN(pnr,h,k,q);
            }
            for (plrj = START(pnr); plrj; plrj=plrj->next){
               if ( INDEX(NBNODE(plrj)) > INDEX(pnr) )
                  for (p = k1; p < kk; p++)
                     COEFF_NN(plrj,h,p,k) -= pivot[p]*COEFF_NN(plrj,h,k,k);
               for (p = k1; p < kk; p++)
                  for (q = k1; q < kk; q++)
                     COEFF_NN(plrj,h,p,q) -= pivot[p]*COEFF_NN(plrj,h,k,q);
            }
            if (nn)
               for (pnfrj = NFSTART(pnr); pnfrj; pnfrj = pnfrj->next)
                  for (p = k1; p < kk; p++)
                     for (q = 0; q < nn; q++)
                        COEFF_NN(pnfrj,h,p,q) -= pivot[p]*COEFF_NN(pnfrj,h,k,q);
            if (mm)
               for (pnerj = NESTART(pnr); pnerj; pnerj = pnerj->next){
                  pel = NBELEM(pnerj);
                  NODE_INDEX(pel,pnr,r)
                  for (p = k1; p < kk; p++)
                     for (q = 0; q < mm; q++)
                        COEFF_NE_NM(pel,h,r,p,q) -= 
                                              pivot[p]*COEFF_NE_NM(pel,h,r,k,q);
               }
         }
         for (plri = START(pnr); plri; plri = plri->next){
            for (plir = START(pni=NBNODE(plri)); NBNODE(plir) != pnr; 
                                                             plir = plir->next);
            for (p = k1; p < kk; p++)
               pivot[p] = COEFF_NN(plir,h,p,k) /= COEFF_NN(pnr,h,k,k);
            if ( INDEX(pni) > INDEX(pnr) ){
               pivot[k] = COEFF_NN(plir,h,k,k) /= COEFF_NN(pnr,h,k,k);
               COEFF_NN(pni,h,k,k) -= pivot[k]*COEFF_NN(plri,h,k,k);
               for (p = k1; p < kk; p++){
                  COEFF_NN(pni,h,p,k)  -= pivot[p]*COEFF_NN(plri,h,k,k);
                  COEFF_NN(pni,h,k,p)  -= pivot[k]*COEFF_NN(plri,h,k,p);
                  COEFF_NN(plir,h,k,p) -= pivot[k]*COEFF_NN(pnr,h,k,p);
               }
            }
            for (p = k1; p < kk; p++)
               for (q = k1; q < kk; q++){
                  COEFF_NN(pni,h,p,q)  -= pivot[p]*COEFF_NN(plri,h,k,q);
                  COEFF_NN(plir,h,p,q) -= pivot[p]*COEFF_NN(pnr,h,k,q);
            }
            for (plrj = START(pnr); plrj; plrj=plrj->next)   
               if ((pnj=NBNODE(plrj)) != pni){
                  for (plij = START(pni); plij && NBNODE(plij) != pnj;
                                                             plij = plij->next);
                  if (plij){
                     if ( INDEX(pni) > INDEX(pnr) ){
                        if ( INDEX(pnj) > INDEX(pnr) )
                           COEFF_NN(plij,h,k,k) -= pivot[k]*COEFF_NN(plrj,h,k,k);
                        for (q = k1; q < kk; q++)
                           COEFF_NN(plij,h,k,q) -= pivot[k]*COEFF_NN(plrj,h,k,q);
                     }
                     if ( INDEX(pnj) > INDEX(pnr) )
                        for (p = k1; p < kk; p++)
                           COEFF_NN(plij,h,p,k) -= pivot[p]*COEFF_NN(plrj,h,k,k);
                     for (p = k1; p < kk; p++)
                        for (q = k1; q < kk; q++)
                           COEFF_NN(plij,h,p,q) -= pivot[p]*COEFF_NN(plrj,h,k,q);
                  }
               }
            if (nn)
               for (pnfrj = NFSTART(pnr); pnfrj; pnfrj = pnfrj->next){
                  pfj = NBFACE(pnfrj);
                  for (pnfij = NFSTART(pni); pnfij && NBFACE(pnfij) != pfj;
                                                           pnfij = pnfij->next);
                  if (pnfij){
                     if ( INDEX(pni) > INDEX(pnr) )
                        for (q = 0; q < nn; q++)
                           COEFF_NN(pnfij,h,k,q) -= 
                                                 pivot[k]*COEFF_NN(pnfrj,h,k,q);
                     for (p = k1; p < kk; p++)
                        for (q = 0; q < nn; q++)
                           COEFF_NN(pnfij,h,p,q) -= 
                                                 pivot[p]*COEFF_NN(pnfrj,h,k,q);
                  }
               }
            if (mm)
               for (pnerj = NESTART(pnr); pnerj; pnerj = pnerj->next){
                  pel = NBELEM(pnerj);
                  NODE_INDEX(pel,pni,i)
                  if (i > -1){
                     NODE_INDEX(pel,pnr,r)
                     if ( INDEX(pni) > INDEX(pnr) )
                        for (q = 0; q < mm; q++)
                           COEFF_NE_NM(pel,h,i,k,q) -= 
                                              pivot[k]*COEFF_NE_NM(pel,h,r,k,q);
                     for (p = k1; p < kk; p++)
                        for (q = 0; q < mm; q++)
                           COEFF_NE_NM(pel,h,i,p,q) -= 
                                              pivot[p]*COEFF_NE_NM(pel,h,r,k,q);
                  }
               }
         }
         if (nn)
            for (pnfri = NFSTART(pnr); pnfri; pnfri=pnfri->next){
               for (pfnir = FNSTART(pfi=NBFACE(pnfri)); NBNODE(pfnir) != pnr; 
                                                           pfnir = pfnir->next);
               for (p = 0; p < nn; p++){
                  pivot[p] = COEFF_NN(pfnir,h,p,k) /= COEFF_NN(pnr,h,k,k);
                  for (q = k1; q < kk; q++)
                     COEFF_NN(pfnir,h,p,q) -= pivot[p]*COEFF_NN(pnr,h,k,q); 
                  for (q = 0; q < nn; q++)
                     COEFF_NN(pfi,h,p,q) -= pivot[p]*COEFF_NN(pnfri,h,k,q); 
               }
               for (plrj = START(pnr); plrj; plrj=plrj->next){
                  pnj=NBNODE(plrj);
                  for (pfnij = FNSTART(pfi); pfnij && NBNODE(pfnij) != pnj; 
                                                           pfnij = pfnij->next);
                  if (pfnij){
                     if ( INDEX(pnj) > INDEX(pnr) )
                        for (p = 0; p < nn; p++)
                           COEFF_NN(pfnij,h,p,k) -= 
                                                  pivot[p]*COEFF_NN(plrj,h,k,k);
                     for (p = 0; p < nn; p++)
                        for (q = k1; q < kk; q++)
                           COEFF_NN(pfnij,h,p,q) -= 
                                                  pivot[p]*COEFF_NN(plrj,h,k,q);
                  }
               }
               for (pnfrj = NFSTART(pnr); pnfrj; pnfrj = pnfrj->next)
                  if ( (pfj = NBFACE(pnfrj)) != pfi){
                     for (pflij = FSTART(pfi); pflij && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                     if (pflij)
                        for (p = 0; p < nn; p++)
                           for (q = 0; q < nn; q++)
                              COEFF_NN(pflij,h,p,q) -= 
                                                 pivot[p]*COEFF_NN(pnfrj,h,k,q);
                  }
               if (mm)
                  for (pnerj = NESTART(pnr); pnerj; pnerj = pnerj->next){
                     pel = NBELEM(pnerj);
                     NODE_INDEX(pel,pnr,r)
                     FACE_INDEX(pel,pfi,i)
                     if (i > -1)
                        for (p = 0; p < nn; p++)
                           for (q = 0; q < mm; q++)
                              COEFF_FE_NM(pel,h,i,p,q) -= 
                                              pivot[p]*COEFF_NE_NM(pel,h,r,k,q);
                  }
            }
         if (mm)
            for (pneri = NESTART(pnr); pneri; pneri = pneri->next){
               pel = NBELEM(pneri);
               NODE_INDEX(pel,pnr,r)
               for (p = 0; p < mm; p++){
                  pivot[p] = COEFF_EN_MN(pel,h,r,p,k) /= COEFF_NN(pnr,h,k,k);
                  for (q = k1; q < kk; q++)
                     COEFF_EN_MN(pel,h,r,p,q) -= pivot[p]*COEFF_NN(pnr,h,k,q);
               }
               for (j = 0; j < NVERT; j++)
                  if (pel->n[j] != pnr && IS_FN(pel->n[j])){
                     for (plrj = START(pnr); NBNODE(plrj) != pel->n[j]; 
                                                             plrj = plrj->next);
                     if ( INDEX(pel->n[j]) > INDEX(pnr) )
                        for (p = 0; p < mm; p++)
                           COEFF_EN_MN(pel,h,j,p,k) -= 
                                                  pivot[p]*COEFF_NN(plrj,h,k,k);
                     for (p = 0; p < mm; p++)
                        for (q = k1; q < kk; q++)
                           COEFF_EN_MN(pel,h,j,p,q) -= 
                                                  pivot[p]*COEFF_NN(plrj,h,k,q);
                  }
               for (j = 0; j < SIDES; j++)
                  if (IS_FF(pel->f[j])){
                     for (pnfrj = NFSTART(pnr); NBFACE(pnfrj) != pel->f[j]; 
                                                           pnfrj = pnfrj->next);
                     for (p = 0; p < mm; p++)
                        for (q = 0; q < nn; q++)
                           COEFF_EF_MN(pel,h,j,p,q) -= 
                                                 pivot[p]*COEFF_NN(pnfrj,h,k,q);
                  }
               for (p = 0; p < mm; p++)
                  for (q = 0; q < mm; q++)
                     COEFF_EE_MM(pel,h,p,q) -= pivot[p]*COEFF_NE_NM(pel,h,r,k,q);
            }
      }
   }
   for (k = 0; k < nn; k++){
      k1 = k + 1;
      for (pfr = FIRSTFACE(theGrid); pfr; pfr = pfr->succ){
         if (mm){
            i = 1;
            for (pnerj = NESTART(NBNODE(TFNSTART(pfr))); i; pnerj =pnerj->next){
               pel = NBELEM(pnerj);
               FACE_INDEX(pel,pfr,r)
               if (r > -1) i = 0;
            }   /* pel possesses the face pfr */
            pnr = pel->n[(r+1) % NVERT];   /* vertex of pfr */
         }
         if (k1 < nn){
            for (p = k1; p < nn; p++){
               pivot[p] = COEFF_NN(pfr,h,p,k) /= COEFF_NN(pfr,h,k,k);
               for (q = k1; q < nn; q++)
                  COEFF_NN(pfr,h,p,q) -= pivot[p]*COEFF_NN(pfr,h,k,q);
            }
            for (pflrj = FSTART(pfr); pflrj; pflrj=pflrj->next){
               if ( INDEX(NBFACE(pflrj)) > INDEX(pfr) )
                  for (p = k1; p < nn; p++)
                     COEFF_NN(pflrj,h,p,k) -= pivot[p]*COEFF_NN(pflrj,h,k,k);
               for (p = k1; p < nn; p++)
                  for (q = k1; q < nn; q++)
                     COEFF_NN(pflrj,h,p,q) -= pivot[p]*COEFF_NN(pflrj,h,k,q);
            }
            if (mm)
               for (pnerj = NESTART(pnr); pnerj; pnerj = pnerj->next){
                  pel = NBELEM(pnerj);
                  FACE_INDEX(pel,pfr,r)
                  if (r > -1)
                     for (p = k1; p < nn; p++)
                        for (q = 0; q < mm; q++)
                           COEFF_FE_NM(pel,h,r,p,q) -= 
                                              pivot[p]*COEFF_FE_NM(pel,h,r,k,q);
               }
         }
         for (pflri = FSTART(pfr); pflri; pflri = pflri->next){
            for (pflir = FSTART(pfi=NBFACE(pflri)); NBFACE(pflir) != pfr; 
                                                           pflir = pflir->next);
            for (p = k1; p < nn; p++)
               pivot[p] = COEFF_NN(pflir,h,p,k) /= COEFF_NN(pfr,h,k,k);
            if ( INDEX(pfi) > INDEX(pfr) ){
               pivot[k] = COEFF_NN(pflir,h,k,k) /= COEFF_NN(pfr,h,k,k);
               COEFF_NN(pfi,h,k,k) -= pivot[k]*COEFF_NN(pflri,h,k,k);
               for (p = k1; p < nn; p++){
                  COEFF_NN(pfi,h,p,k)  -= pivot[p]*COEFF_NN(pflri,h,k,k);
                  COEFF_NN(pfi,h,k,p)  -= pivot[k]*COEFF_NN(pflri,h,k,p);
                  COEFF_NN(pflir,h,k,p) -= pivot[k]*COEFF_NN(pfr,h,k,p);
               }
            }
            for (p = k1; p < nn; p++)
               for (q = k1; q < nn; q++){
                  COEFF_NN(pfi,h,p,q)  -= pivot[p]*COEFF_NN(pflri,h,k,q);
                  COEFF_NN(pflir,h,p,q) -= pivot[p]*COEFF_NN(pfr,h,k,q);
            }
            for (pflrj = FSTART(pfr); pflrj; pflrj=pflrj->next)   
               if ((pfj=NBFACE(pflrj)) != pfi){
                  for (pflij = FSTART(pfi); pflij && NBFACE(pflij) != pfj;
                                                           pflij = pflij->next);
                  if (pflij){
                     if ( INDEX(pfi) > INDEX(pfr) ){
                        if ( INDEX(pfj) > INDEX(pfr) )
                           COEFF_NN(pflij,h,k,k) -= 
                                                 pivot[k]*COEFF_NN(pflrj,h,k,k);
                        for (q = k1; q < nn; q++)
                           COEFF_NN(pflij,h,k,q) -= 
                                                 pivot[k]*COEFF_NN(pflrj,h,k,q);
                     }
                     if ( INDEX(pfj) > INDEX(pfr) )
                        for (p = k1; p < nn; p++)
                           COEFF_NN(pflij,h,p,k) -= 
                                                 pivot[p]*COEFF_NN(pflrj,h,k,k);
                     for (p = k1; p < nn; p++)
                        for (q = k1; q < nn; q++)
                           COEFF_NN(pflij,h,p,q) -= 
                                                 pivot[p]*COEFF_NN(pflrj,h,k,q);
                  }
               }
            if (mm)
               for (pnerj = NESTART(pnr); pnerj; pnerj = pnerj->next){
                  pel = NBELEM(pnerj);
                  FACE_INDEX(pel,pfr,r)
                  FACE_INDEX(pel,pfi,i)
                  if (r > -1 && i > -1){
                     if ( INDEX(pfi) > INDEX(pfr) )
                        for (q = 0; q < mm; q++)
                           COEFF_FE_NM(pel,h,i,k,q) -= 
                                              pivot[k]*COEFF_FE_NM(pel,h,r,k,q);
                     for (p = k1; p < nn; p++)
                        for (q = 0; q < mm; q++)
                           COEFF_FE_NM(pel,h,i,p,q) -= 
                                              pivot[p]*COEFF_FE_NM(pel,h,r,k,q);
                  }
               }
         }
         if (mm)
            for (pneri = NESTART(pnr); pneri; pneri = pneri->next){
               pel = NBELEM(pneri);
               FACE_INDEX(pel,pfr,r)
               if (r > -1){
                  for (p = 0; p < mm; p++){
                     pivot[p] = COEFF_EF_MN(pel,h,r,p,k) /= COEFF_NN(pfr,h,k,k);
                     for (q = k1; q < nn; q++)
                        COEFF_EF_MN(pel,h,r,p,q) -= pivot[p]*COEFF_NN(pfr,h,k,q);
                  }
                  for (j = 0; j < SIDES; j++)
                     if (pel->f[j] != pfr && IS_FF(pel->f[j])){
                        for (pflrj = FSTART(pfr); NBFACE(pflrj) != pel->f[j]; 
                                                           pflrj = pflrj->next);
                        if ( INDEX(pel->f[j]) > INDEX(pfr) )
                           for (p = 0; p < mm; p++)
                              COEFF_EF_MN(pel,h,j,p,k) -= 
                                                 pivot[p]*COEFF_NN(pflrj,h,k,k);
                        for (p = 0; p < mm; p++)
                           for (q = k1; q < nn; q++)
                              COEFF_EF_MN(pel,h,j,p,q) -= 
                                                 pivot[p]*COEFF_NN(pflrj,h,k,q);
                     }
                  for (p = 0; p < mm; p++)
                     for (q = 0; q < mm; q++)
                        COEFF_EE_MM(pel,h,p,q) -= 
                                              pivot[p]*COEFF_FE_NM(pel,h,r,k,q);
               }
            }
      }
   }
   for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
      for (r = 0; r < mm; r++)
         for (i = r+1; i < mm; i++){
            piv = COEFF_EE_MM(pel,h,i,r) /= COEFF_EE_MM(pel,h,r,r);
            for (j = r+1; j < mm; j++)
               COEFF_EE_MM(pel,h,i,j) -= piv*COEFF_EE_MM(pel,h,r,j);
         }
}

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void solveILU_general(theGrid,h,x,y) 
GRID *theGrid;
INT h, x, y;
{
   ELEMENT *pel;
   NODE *pni, *pnj;
   FACE *pfi, *pfj;
   LINK *plij;
   NFLINK *pnfij;
   FNLINK *pfnij;
   FLINK *pflij;
   FLOAT sum;
   INT i, j, k, q, ind, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC, mm=N_OF_ELEM_FUNC;
 
   /* find x: L x = y */
   for (k = 0; k < kk; k++)
      for (pni = FIRSTNODE(theGrid); pni; pni = pni->succ){
         sum = 0.;
         for (q = 0; q < k; q++)
            sum += COEFF_NN(pni,h,k,q)*NDMV(pni,x,q);
         for (plij = START(pni); plij; plij = plij->next){
            if ( INDEX(pnj = NBNODE(plij)) < INDEX(pni) )
               sum += COEFF_NN(plij,h,k,k)*NDMV(pnj,x,k);
            for (q = 0; q < k; q++)
               sum += COEFF_NN(plij,h,k,q)*NDMV(pni,x,q);
         }
         NDMV(pni,x,k) = NDMV(pni,y,k) - sum;
      }
   if (kk && nn)
      for (pfi = FIRSTFACE(theGrid); pfi; pfi = pfi->succ)
         for (pfnij = FNSTART(pfi); pfnij; pfnij = pfnij->next)
            MMSET2(FDMVP(pfi,x),COEFF_NNP(pfnij,h),NDMVP(NBNODE(pfnij),x),
                                                                      i,j,nn,kk)
   for (k = 0; k < nn; k++)
      for (pfi = FIRSTFACE(theGrid); pfi; pfi = pfi->succ){
         sum = 0.;
         for (q = 0; q < k; q++)
            sum += COEFF_NN(pfi,h,k,q)*FDMV(pfi,x,q);
         for (pflij = FSTART(pfi); pflij; pflij = pflij->next){
            if ( INDEX(pfj = NBFACE(pflij)) < INDEX(pfi) )
               sum += COEFF_NN(pflij,h,k,k)*FDMV(pfj,x,k);
            for (q = 0; q < k; q++)
               sum += COEFF_NN(pflij,h,k,q)*FDMV(pfi,x,q);
         }
         FDMV(pfi,x,k) = FDMV(pfi,y,k) - FDMV(pfi,x,k) - sum;
      }
   if (kk)
      for (pni = FDBN(theGrid); pni != FIRSTNODE(theGrid); pni = pni->succ)
         GSET7(NDMVP(pni,x),0.,i,kk)
   if (nn)
      for (pfi = FDBF(theGrid); pfi != FIRSTFACE(theGrid); pfi = pfi->succ)
         GSET7(FDMVP(pfi,x),0.,i,nn)
   if (mm)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ){
         for (k=0; k < NVERT; k++){
            MMSET2(EDMVP(pel,x),COEFF_EN_MNP(pel,h,k),NDMVP(pel->n[k],x),
                                                                      i,j,mm,kk)
            MMSET4(EDMVP(pel,x),COEFF_EF_MNP(pel,h,k),FDMVP(pel->f[k],x),
                                                                      i,j,mm,nn)
         }
         for (i = 0; i < mm; i++){
            sum = 0.;
            for (j = 0; j < i; j++)
               sum += COEFF_EE_MM(pel,h,i,j)*EDMV(pel,x,j);
            EDMV(pel,x,i) = EDMV(pel,y,i) - EDMV(pel,x,i) - sum;
         }
      }
   
   /* find solution of U . = x and store it in x */
   for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
      for (i = mm-1; i >= 0; i--){
         sum = 0.;
         for (j = i+1; j < mm; j++)
            sum += COEFF_EE_MM(pel,h,i,j)*EDMV(pel,x,j);
         EDMV(pel,x,i) = (EDMV(pel,x,i) - sum)/COEFF_EE_MM(pel,h,i,i);
      }
   if (nn && mm)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         for (k = 0; k < NVERT; k++)
            MMSET9(FDMVP(pel->f[k],x),COEFF_FE_NMP(pel,h,k),EDMVP(pel,x),
                                                                      i,j,nn,mm)
   for (k = nn-1; k >= 0; k--){
      ind = INDEX(FIRSTFACE(theGrid)->prev);
      for (pfi = theGrid->lastFace; INDEX(pfi) != ind; pfi = pfi->prev){
         sum = 0.;
         for (q = k+1; q < nn; q++)
            sum += COEFF_NN(pfi,h,k,q)*FDMV(pfi,x,q);
         for (pflij = FSTART(pfi); pflij; pflij = pflij->next){
            if ( INDEX(pfj = NBFACE(pflij)) > INDEX(pfi) )
               sum += COEFF_NN(pflij,h,k,k)*FDMV(pfj,x,k);
            for (q = k+1; q < nn; q++)
               sum += COEFF_NN(pflij,h,k,q)*FDMV(pfj,x,q);
         }
         FDMV(pfi,x,k) = (FDMV(pfi,x,k)-sum)/COEFF_NN(pfi,h,k,k);
      }
   }
   if (kk && mm)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         for (k = 0; k < NVERT; k++)
            MMSET9(NDMVP(pel->n[k],x),COEFF_NE_NMP(pel,h,k),EDMVP(pel,x),
                                                                      i,j,kk,mm)
   if (kk && nn)
      for (pni = FIRSTNODE(theGrid); pni; pni = pni->succ)
         for (pnfij = NFSTART(pni); pnfij; pnfij = pnfij->next)
            MMSET9(NDMVP(pni,x),COEFF_NNP(pnfij,h),FDMVP(NBFACE(pnfij),x),
                                                                      i,j,kk,nn)
   for (k = kk-1; k >= 0; k--){
      ind = INDEX(FIRSTNODE(theGrid)->prev);
      for (pni = theGrid->lastNode; INDEX(pni) != ind; pni = pni->prev){
         sum = 0.;
         for (q = k+1; q < kk; q++)
            sum += COEFF_NN(pni,h,k,q)*NDMV(pni,x,q);
         for (plij = START(pni); plij; plij = plij->next){
            if ( INDEX(pnj = NBNODE(plij)) > INDEX(pni) )
               sum += COEFF_NN(plij,h,k,k)*NDMV(pnj,x,k);
            for (q = k+1; q < kk; q++)
               sum += COEFF_NN(plij,h,k,q)*NDMV(pnj,x,q);
         }
         NDMV(pni,x,k) = (NDMV(pni,x,k)-sum)/COEFF_NN(pni,h,k,k);
      }
   }
}

#else

void ILU_general(theGrid,a,h)
GRID *theGrid; INT a,h;
{  eprintf("Error: ILU_general not available.\n");  }

void solveILU_general(theGrid,h,x,y) 
GRID *theGrid; INT h, x, y;
{  eprintf("Error: solveILU_general not available.\n");  }

#endif

void make_ILU(theGrid,a,h,t,type,structure)
GRID *theGrid;
INT a, h, t, type, structure;
{
   switch(type){
   case Q_SN: sILU(theGrid,a,h,t);
        break;
   case Q_VN: 
           switch(structure){
           case Q_NEBDIAG: nn_copy_mat_to_mat(theGrid,a,h);
                           nn_ILU(theGrid,a,h);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_SF: ILU_f(theGrid,a,h);
        break;
   case Q_VF: 
           switch(structure){
           case Q_FULL: vILU_vf(theGrid,a,h);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_SNSF: sILU_sn_sf(theGrid,a,h);
        break;
   case Q_VNSF: 
           switch(structure){
           case Q_FULL: vILU_vn_sf_3d(theGrid,a,h);
                break; 
           case Q_NEBDIAG: sILU_vn_sf(theGrid,a,h);
                break; 
           case Q_NEBDIAG_NRED: sILU_vn_sf_nred(theGrid,a,h);
                break;
           case Q_NEBDIAG_NFRED: sILU_vn_sf_nfred(theGrid,a,h);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_VNVF: 
           switch(structure){
           case Q_FULL: vILU_vn_vf_2d(theGrid,a,h);
                break; 
           case Q_EBDIAG: sILU_sn_sf(theGrid,a,h);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_VNVE: 
           switch(structure){
           case Q_NEBDIAG_EDIAG: nn_copy_mat_to_mat(theGrid,a,h);
                                 nn_ILU(theGrid,a,h);
                                 copy(theGrid,a,h,0,Q_SE);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_VNVFVE: 
           switch(structure){
           case Q_EBDIAG: sILU_sn_sf_se(theGrid,a,h);
                break; 
           default:
                eprintf("Error: make_ILU not available.\n");
                break;
           }
        break;
   case Q_GENERAL: ILU_general(theGrid,a,h);
        break;
   default:
        eprintf("Error: make_ILU not available.\n");
        break;
   }
}

void make_ILU_beta_full(theGrid,a,h,q,beta,t,type,structure)
GRID *theGrid;
INT a, h, q, t, type, structure;
FLOAT beta;
{
   switch(type){
   case Q_SF: ILU_f_beta_full(theGrid,a,h,q,beta);
        break;
   default:
        eprintf("Error: make_ILU_beta_full not available.\n");
        break;
   }
}

void make_ILU_beta_sparse(theGrid,a,h,q,beta,t,type,structure)
GRID *theGrid;
INT a, h, q, t, type, structure;
FLOAT beta;
{
   switch(type){
   case Q_SN: ILU_n_beta_sparse(theGrid,a,h,q,beta);
        break;
   case Q_SF: ILU_f_beta_sparse(theGrid,a,h,q,beta);
        break;
   case Q_SNSF: eprintf("Error: make_ILU_beta_sparse not available.\n");
        break;
   default:
        eprintf("Error: make_ILU_beta_sparse not available.\n");
        break;
   }
}

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */
void solve_ILU(theGrid,h,x,y,t,type,structure)
GRID *theGrid;
INT h, x, y, t, type, structure;
{
   switch(type){
   case Q_SN: ssolveILU(theGrid,h,x,y,t);
        break;
   case Q_VN: 
           switch(structure){
           case Q_NEBDIAG: nn_solveILU(theGrid,h,x,y);
                break; 
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_SF: solveILU_f(theGrid,h,x,y);
        break;
   case Q_VF: 
           switch(structure){
           case Q_FULL: vsolveILU_vf_2d(theGrid,h,x,y);
                break; 
           case Q_FEBDIAG: ssolveILU_vf(theGrid,h,x,y);
                break; 
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_DVF: 
           switch(structure){
           case Q_DVFEBDIAG: dvsolveILU_f(theGrid,h,x,y);
                break; 
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_SNSF: ssolveILU_sn_sf(theGrid,h,x,y);
        break;
   case Q_VNSF: 
           switch(structure){
           case Q_FULL: vsolveILU_vn_sf_3d(theGrid,h,x,y);
                break; 
           case Q_NEBDIAG: ssolveILU_vn_sf(theGrid,h,x,y);
                break; 
           case Q_NEBDIAG_NRED: ssolveILU_vn_sf_nred(theGrid,h,x,y);
                break;
           case Q_NEBDIAG_NFRED: ssolveILU_vn_sf_nfred(theGrid,h,x,y);
                break;
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_VNVF: 
           switch(structure){
           case Q_FULL: vsolveILU_vn_vf_2d(theGrid,h,x,y);
                break; 
           case Q_EBDIAG: ssolveIL_vn_vf(theGrid,h,x,y);
                          ssolveIU_vn_vf(theGrid,h,x);
                break; 
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_VNVE: nn_solveILU(theGrid,h,x,y);
                vdivide_1e(theGrid,y,h,x);
        break;
   case Q_VNVFVE: 
           switch(structure){
           case Q_EBDIAG: ssolveILU_vn_vf_ve(theGrid,h,x,y);
                break; 
           default:
                eprintf("Error: solve_ILU not available.\n");
                break;
           }
        break;
   case Q_GENERAL: solveILU_general(theGrid,h,x,y);
        break;
   default:
        eprintf("Error: solve_ILU not available.\n");
        break;
   }
}

void ILU_step(tGrid,Z,h,f,u,d,q,t,type,structure)
GRID *tGrid;                                /*  u := u - (ILU)^{-1}(Z u - f)  */
INT Z, h, f, u, d, q, t, type, structure;
{
   defect(tGrid,Z,f,u,d,q,t,type,structure);
   solve_ILU(tGrid,h,q,d,t,type,structure);
   subtr(tGrid,u,q,u,t,type);
}

void multB_ILU_BT(tGrid,matrix_b,x,y,q,t_row,t_column,row_type,column_type,
         b_struct,ilu_matrix,u1,u2,t_for_u,u_type,a_struct,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_p; row_type = column_type = p_type */
INT x, y, q;  /*  p_type  */
INT u1, u2;   /*  u_type  */
INT matrix_b,t_row,t_column,row_type,column_type,b_struct,
    ilu_matrix,t_for_u,u_type,a_struct,p6,p7,p8,p9;
FLOAT r1, r2, r3;
{  
   mult_BT(tGrid,matrix_b,x,u1,u1,t_row,t_for_u,row_type,u_type,b_struct);
   solve_ILU(tGrid,ilu_matrix,u2,u1,t_for_u,u_type,a_struct);
   mult_B(tGrid,matrix_b,u2,y,y,t_row,t_for_u,row_type,u_type,b_struct);
   if (t_row & ZERO_MEAN)
      add_value(tGrid,BSC*sum(tGrid,x,t_row,row_type),y,t_row,row_type);
}

void multB_SSOR_BT(tGrid,matrix_b,x,y,q,t_row,t_column,row_type,column_type,
     b_struct,precond_matrix,u1,u2,t_for_u,u_type,a_struct,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_p; row_type = column_type = p_type */
INT x, y, q;  /*  p_type  */
INT u1, u2;   /*  u_type  */
INT matrix_b,t_row,t_column,row_type,column_type,b_struct,
    precond_matrix,t_for_u,u_type,a_struct,p6,p7,p8,p9;
FLOAT r1, r2, r3;
{  
   mult_BT(tGrid,matrix_b,x,u1,u1,t_row,t_for_u,row_type,u_type,b_struct);
   set_value(tGrid,0.,u2,t_for_u,u_type);
   SSOR_step(tGrid,precond_matrix,u2,u2,u2,t_for_u,t_for_u,u_type,u_type,
             a_struct,u1,1,2,3,4,5,6,7,8,9,1.,0.,0.);
   mult_B(tGrid,matrix_b,u2,y,y,t_row,t_for_u,row_type,u_type,b_struct);
   if (t_row & ZERO_MEAN)
      add_value(tGrid,BSC*sum(tGrid,x,t_row,row_type),y,t_row,row_type);
}

void multB_diag_BT(tGrid,matrix_b,x,y,q,t_row,t_column,row_type,column_type,
     b_struct,precond_matrix,u1,u2,t_for_u,u_type,a_struct,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_p; row_type = column_type = p_type */
INT x, y, q;  /*  p_type  */
INT u1, u2;   /*  u_type  */
INT matrix_b,t_row,t_column,row_type,column_type,b_struct,
    precond_matrix,t_for_u,u_type,a_struct,p6,p7,p8,p9;
FLOAT r1, r2, r3;
{  
   mult_BT(tGrid,matrix_b,x,u1,u1,t_row,t_for_u,row_type,u_type,b_struct);
   mult_inv_diag(tGrid,precond_matrix,u1,u2,u2,
                 t_for_u,t_for_u,u_type,u_type,a_struct);
   mult_B(tGrid,matrix_b,u2,y,y,t_row,t_for_u,row_type,u_type,b_struct);
   if (t_row & ZERO_MEAN)
      add_value(tGrid,BSC*sum(tGrid,x,t_row,row_type),y,t_row,row_type);
}

void multBBT(tGrid,matrix_b,x,y,q,t_row,t_column,row_type,column_type,
         b_struct,p0,u1,p2,t_for_u,u_type,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_p; row_type = column_type = p_type */
INT x, y, q;  /*  p_type  */
INT u1;   /*  u_type  */
INT matrix_b,t_row,t_column,row_type,column_type,b_struct,
    t_for_u,u_type,p5,p6,p7,p8,p9;
FLOAT r1, r2, r3;
{  
   mult_BT(tGrid,matrix_b,x,u1,u1,t_row,t_for_u,row_type,u_type,b_struct);
   mult_B(tGrid,matrix_b,u1,y,y,t_row,t_for_u,row_type,u_type,b_struct);
   if (t_row & ZERO_MEAN)
      add_value(tGrid,BSC*sum(tGrid,x,t_row,row_type),y,t_row,row_type);
}

/*  y = (matrix derived from precond_matrix) x  */
void preconditioning(theGrid,precond_matrix,x,y,t,type,q0,q1,q2,precond_type)
GRID *theGrid;
INT precond_matrix, x, y, t, type, q0, q1, q2, precond_type;
{
   switch(precond_type){
   case IDENT_PR: copy(theGrid,x,y,t,type);
        break;
   case DIAG_PR: divide(theGrid,x,precond_matrix,y,t,type);
        break;
   case MDIAG_PR: mult_inv_diag(theGrid,precond_matrix,x,y,-1,t,t,type,type,q0);
        break;
   case SSOR_PR:
           set_value(theGrid,0.,y,t,type);
           SSOR_step(theGrid,precond_matrix,y,y,y,t,t,type,type,q0,
                     x,1,2,3,4,5,6,7,8,9,1.,0.,0.);
        break;
   case ILU_PR: solve_ILU(theGrid,precond_matrix,y,x,t,type,q0);
        break;
   case UMF_PR:
           fill_Ax(theGrid,precond_matrix,Ap,Ai,Ax,type,q0); /* q0 = U_SPACE */
           make_vector_from_grid_data(theGrid,x,Rhs,type,q0);
           solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
           make_grid_data_from_vector(theGrid,y,Sol,type,q0);
        break;
   default:
        eprintf("Error: preconditioning not available.\n");
        break;
   }
}
