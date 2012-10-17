/******************************************************************************/
/*                                                                            */
/*                        restriction and prolongation                        */
/*                                                                            */
/******************************************************************************/

#if DIM == 3

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/* saves restriction of f in g; theGrid is the coarser grid */
void u_restriction(theGrid,f,g)
GRID *theGrid;
INT f, g;
{
   ELEMENT *pelem;
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLOAT sum0, sum1, sum2;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      sum0 = 0.0;
      sum1 = 0.0;
      sum2 = 0.0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         sum0 += ND(NBNODE(plink),f,0);
         sum1 += ND(NBNODE(plink),f,1);
         sum2 += ND(NBNODE(plink),f,2);
      }
      ND(pnode,g,0) = ND(pnode->son,f,0) + sum0/2.0;
      ND(pnode,g,1) = ND(pnode->son,f,1) + sum1/2.0;
      ND(pnode,g,2) = ND(pnode->son,f,2) + sum2/2.0;
   }
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      FD(pface,g) =0.0;
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      treat_cubic_functions(pelem,restrict_from_faces_of_son,f,g);
}

/* saves restriction of u in v; theGrid is the coarser grid */
void U_restriction(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   NODE *pnode;
   FACE *pface;
   LINK *plink;
   FLOAT sum0, sum1, sum2;
   INT i;
   
   for (pnode = FIRSTNODE(theGrid->finer); pnode != NULL; pnode = pnode->succ)
      if (pnode->father){
         ND(pnode,v,0) = ND(pnode->father,u,0);
         ND(pnode,v,1) = ND(pnode->father,u,1);
         ND(pnode,v,2) = ND(pnode->father,u,2);
      }
      else{
         ND(pnode,v,0) = 0.0;
         ND(pnode,v,1) = 0.0;
         ND(pnode,v,2) = 0.0;
      }
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      sum0 = sum1 = sum2 = 0.;
      i = 0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         sum0 += ND(NBNODE(plink),u,0);
         sum1 += ND(NBNODE(plink),u,1);
         sum2 += ND(NBNODE(plink),u,2);
         i++;
      }
      ND(pnode,v,0) = (ND(pnode->son,u,0) + sum0/i)/2.;
      ND(pnode,v,1) = (ND(pnode->son,u,1) + sum1/i)/2.;
      ND(pnode,v,2) = (ND(pnode->son,u,2) + sum2/i)/2.;
   }
      
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      FD(pface,v) = (FD(pface->sons[0],u) + FD(pface->sons[1],u) +
                     FD(pface->sons[2],u) + FD(pface->sons[3],u))/4.;
}

#else

void u_restriction(theGrid,f,g)
GRID *theGrid; INT f, g;
{  eprintf("Error: u_restriction not available.\n");  }

void U_restriction(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: U_restriction not available.\n");  }

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void p_restriction(theGrid,f,g)
GRID *theGrid;
INT f, g;
{
   ELEMENT *pel;
   
   for (pel = SUCC(FIRSTELEMENT(theGrid)); pel != NULL; pel = pel->succ)
      ED(pel,g) = ( ED(pel->sons[0],f) + ED(pel->sons[1],f) + ED(pel->sons[2],f)
                  + ED(pel->sons[3],f) + ED(pel->sons[4],f) + ED(pel->sons[5],f)
                  + ED(pel->sons[6],f) + ED(pel->sons[7],f) )/8.0;
}

void p_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
      ED(pelem->sons[0],v) = ED(pelem,u);
      ED(pelem->sons[1],v) = ED(pelem,u);
      ED(pelem->sons[2],v) = ED(pelem,u);
      ED(pelem->sons[3],v) = ED(pelem,u);
      ED(pelem->sons[4],v) = ED(pelem,u);
      ED(pelem->sons[5],v) = ED(pelem,u);
      ED(pelem->sons[6],v) = ED(pelem,u);
      ED(pelem->sons[7],v) = ED(pelem,u);
   }
}

#else

void p_restriction(theGrid,f,g)
GRID *theGrid; INT f, g;
{  eprintf("Error: p_restriction not available.\n");  }

void p_prolongation(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: p_prolongation not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && !(DATA_STR & LG_DATA)

/* saves prolongation of u in v; theGrid is the coarser grid */
void u_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   ELEMENT *pelem;
   NODE *pnode;
   LINK *plink;
   FLOAT val0, val1, val2;
   
   for (pnode = FIRSTNODE(theGrid->finer); pnode != NULL; pnode = pnode->succ)
      if (pnode->father){
         ND(pnode,v,0) = ND(pnode->father,u,0);
         ND(pnode,v,1) = ND(pnode->father,u,1);
         ND(pnode,v,2) = ND(pnode->father,u,2);
      }
      else{
         ND(pnode,v,0) = 0.0;
         ND(pnode,v,1) = 0.0;
         ND(pnode,v,2) = 0.0;
      }
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      val0 = ND(pnode,u,0)/2.0;
      val1 = ND(pnode,u,1)/2.0;
      val2 = ND(pnode,u,2)/2.0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         ND(NBNODE(plink),v,0) += val0;
         ND(NBNODE(plink),v,1) += val1;
         ND(NBNODE(plink),v,2) += val2;
      }
   }
      
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      treat_cubic_functions(pelem,prolong_to_faces_of_son,u,v);
}

#else

void u_prolongation(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: u_prolongation not available.\n");  }

#endif

#if DATA_STR & LG_DATA

/* saves prolongation of u in v; theGrid is the coarser grid */
void u_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   ELEMENT *pelem;
   NODE *pnode;
   LINK *plink;
   FLOAT val0, val1, val2;
   
   for (pnode = FIRSTNODE(theGrid->finer); pnode != NULL; pnode = pnode->succ)
      if (pnode->father){
         ND(pnode,v,0) = ND(pnode->father,u,0);
         ND(pnode,v,1) = ND(pnode->father,u,1);
         ND(pnode,v,2) = ND(pnode->father,u,2);
      }
      else{
         ND(pnode,v,0) = 0.0;
         ND(pnode,v,1) = 0.0;
         ND(pnode,v,2) = 0.0;
      }
   for (pnode = FDBN(theGrid); pnode != FIRSTNODE(theGrid); pnode =pnode->succ){
      val0 = ND(pnode,u,0)/2.0;
      val1 = ND(pnode,u,1)/2.0;
      val2 = ND(pnode,u,2)/2.0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
            ND(NBNODE(plink),v,0) += val0;
            ND(NBNODE(plink),v,1) += val1;
            ND(NBNODE(plink),v,2) += val2;
      }
   }
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      val0 = ND(pnode,u,0)/2.0;
      val1 = ND(pnode,u,1)/2.0;
      val2 = ND(pnode,u,2)/2.0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         ND(NBNODE(plink),v,0) += val0;
         ND(NBNODE(plink),v,1) += val1;
         ND(NBNODE(plink),v,2) += val2;
      }
   }
   for (pnode = FDBN(theGrid->finer); pnode != FIRSTNODE(theGrid->finer); 
                                                            pnode = pnode->succ)
      if (pnode->lgd)
         if(pnode->father){
            NDLG(pnode,v,0) = NDLG(pnode->father,u,0);
            NDLG(pnode,v,1) = NDLG(pnode->father,u,1);
         }
         else{
            NDLG(pnode,v,0) = 0.0;
            NDLG(pnode,v,1) = 0.0;
         }   
      
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      treat_cubic_functions(pelem,prolong_to_faces_of_son,u,v);
}

#endif

#else  /* if DIM == 2 */

/* saves restriction of f in g; theGrid is the coarser grid */
void u_restriction(theGrid,f,g)
GRID *theGrid;
INT f, g;
{
   eprintf("Error: u_restriction not available.\n");
}

/* saves restriction of u in v; theGrid is the coarser grid */
void U_restriction(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   eprintf("Error: U_restriction not available.\n");
}

void p_restriction(theGrid,f,g)
GRID *theGrid;
INT f, g;
{
   eprintf("Error: p_restriction not available.\n");
}

#if !(DATA_STR & LG_DATA)

/* saves prolongation of u in v; theGrid is the coarser grid */
void u_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   eprintf("Error: u_prolongation not available.\n");
}

#endif

#if DATA_STR & LG_DATA

/* saves prolongation of u in v; theGrid is the coarser grid */
void u_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   eprintf("Error: u_prolongation not available.\n");
}

#endif

void p_prolongation(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   eprintf("Error: p_prolongation not available.\n");
}

#endif  /* DIM == 2 */

