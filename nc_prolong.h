#if DIM == 2

void get_sons_of_face(f3,f31,f32,n11,n22)
FACE *f3, **f31, **f32;
NODE *n11, *n22;
{
   if (n11->index < n22->index){
      *f31 = f3->sons[0];
      *f32 = f3->sons[1];
   }
   else{
      *f31 = f3->sons[1];
      *f32 = f3->sons[0];
   }
}

NODE *unif_middle(n11,n22)
NODE *n11, *n22;
{
   LINK *pli, *plin;

   for (pli = n11->tstart; ; pli = pli->next){
      for(plin = n22->tstart; plin && plin->nbnode != pli->nbnode;
                                                               plin=plin->next);
      if (plin)
         return(plin->nbnode);
   }
}

void get_son_faces(pel,f11,f12,f13,f21,f22,f23,f31,f32,f33)
ELEMENT *pel;
FACE **f11, **f12, **f13, **f21, **f22, **f23, **f31, **f32, **f33;
{
   NODE *n1,  *n2,  *n3,  *n11, *n22, *n33;

   NODES_OF_ELEMENT(n1,n2,n3,pel);
   n11 = n1->son;
   n22 = n2->son;
   n33 = n3->son;
   get_sons_of_face(pel->f[0],f12,f13,n22,n33);
   get_sons_of_face(pel->f[1],f21,f23,n11,n33);
   get_sons_of_face(pel->f[2],f31,f32,n11,n22);
   *f11 = pel->sons[0]->f[0];
   *f22 = pel->sons[1]->f[0];
   *f33 = pel->sons[2]->f[0];
}
   
void get_sons(pel,n1,n2,n3,n12,n13,n23,n11,n22,n33,
              f11,f12,f13,f21,f22,f23,f31,f32,f33)
ELEMENT *pel;
NODE **n1,  **n2,  **n3,  **n12, **n13, **n23, **n11, **n22, **n33;
FACE **f11, **f12, **f13, **f21, **f22, **f23, **f31, **f32, **f33;
{
   NODES_OF_ELEMENT(*n1,*n2,*n3,pel);
   *n11 = (*n1)->son;
   *n22 = (*n2)->son;
   *n33 = (*n3)->son;
   *n12 = unif_middle(*n11,*n22);
   *n13 = unif_middle(*n11,*n33);
   *n23 = unif_middle(*n22,*n33);
   get_sons_of_face(pel->f[0],f12,f13,*n22,*n33);
   get_sons_of_face(pel->f[1],f21,f23,*n11,*n33);
   get_sons_of_face(pel->f[2],f31,f32,*n11,*n22);
   *f11 = pel->sons[0]->f[0];
   *f22 = pel->sons[1]->f[0];
   *f33 = pel->sons[2]->f[0];
}

INT compare_elements(pel,n1,n2,n3,f1,f2,f3)
ELEMENT *pel;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
{
   Order3(&n1,&n2,&n3,&f1,&f2,&f3);
   if (pel->n[0] == n1 && pel->n[1] == n2 && pel->n[2] == n3 &&
       pel->f[0] == f1 && pel->f[1] == f2 && pel->f[2] == f3)
      return(1);
   else
      return(0);
}

INT comp_elem_sons(pel,n1,n2,n3,f1,f2,f3)
ELEMENT *pel;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
{
   INT i=0, j;

   if (compare_elements(pel->sons[0],n1,n2,n3,f1,f2,f3)){
      i++; j = 1;
   }
   if (compare_elements(pel->sons[1],n1,n2,n3,f1,f2,f3)){
      i++; j = 10;
   }
   if (compare_elements(pel->sons[2],n1,n2,n3,f1,f2,f3)){
      i++; j = 100;
   }
   if (compare_elements(pel->sons[3],n1,n2,n3,f1,f2,f3)){
      i++; j = 1000;
   }
   if (i != 1)
      j = 0;
   return(j);
}

void p1c_coefficients(a1,a2,a3,u1,u2,u3,u12,u13,u23)
FLOAT a1, a2, a3, *u1, *u2, *u3, *u12, *u13, *u23;
{
   *u1 = a1;
   *u2 = a2;
   *u3 = a3;
   if (fabs(*u12) < 1.e-16)
      *u12 = 0.5*(a1+a2);
   if (fabs(*u13) < 1.e-16)
      *u13 = 0.5*(a1+a3);
   if (fabs(*u23) < 1.e-16)
      *u23 = 0.5*(a2+a3);
}

#if N_DATA & SCALAR_NODE_DATA

void scalar_p1c_prolongation(tGrid,u,v)
GRID *tGrid;     /*  coarse grid  */
INT u, v;        /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      n11 = n1->son;
      n22 = n2->son;
      n33 = n3->son;
      n12 = unif_middle(n11,n22);
      n13 = unif_middle(n11,n33);
      n23 = unif_middle(n22,n33);
      p1c_coefficients(NDS(n1,u),   NDS(n2,u),   NDS(n3,u),
                      &NDS(n11,v), &NDS(n22,v), &NDS(n33,v),
                      &NDS(n12,v), &NDS(n13,v), &NDS(n23,v));
   }
}

void scalar_p1c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid;      /*  coarse grid  */
INT f, g, t;      /*  f ... fine grid; g ... coarse grid  */
{
   NODE *theNode;
   LINK *plink;

   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode)){
      NDS(theNode,g) = 0.;
      for (plink = T_START(theNode->son,t); plink != NULL; plink = plink->next)
         NDS(theNode,g) += NDS(NBNODE(plink),f);
      NDS(theNode,g) = NDS(theNode->son,f) + 0.5*NDS(theNode,g);
   }
}

void scalar_p2c_rhs_restrict_BD();

void p_scalar_p1c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid;      /*  coarse grid  */
INT f, g, t;      /*  f ... fine grid; g ... coarse grid  */
{
   NODE *theNode, *first_n;

   if (t & WITHOUT_FIRST){      
      first_n = FIRSTN(tGrid->finer);
      NDS(first_n,f) = 0.;
      for (theNode = first_n->succ; theNode; theNode=SUCC(theNode))
         NDS(first_n,f) -= NDS(theNode,f);
   }
   scalar_p1c_rhs_restrict(tGrid,f,g,0);
   if (t & ADD_ZN_CONDITION)
      scalar_p2c_rhs_restrict_BD(tGrid,f,g);
}

#else

void scalar_p1c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p1c_prolongation not available.\n");  }

void scalar_p1c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid; INT f, g, t;
{  eprintf("Error: scalar_p1c_rhs_restrict not available.\n");  }

void p_scalar_p1c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid; INT f, g, t;
{  eprintf("Error: p_scalar_p1c_rhs_restrict not available.\n");  }

#endif

#if E_DATA & VECTOR_ELEMENT_DATA

void vect_prolong_bubble1(pel,u,v)
ELEMENT *pel;
INT u, v;
{
   DOUBLE a[DIM];

   SET2(a,EDVP(pel,u),0.5)
   SET1(EDVP(pel->sons[0],v),a)
   SET1(EDVP(pel->sons[1],v),a)
   SET1(EDVP(pel->sons[2],v),a)
   SET1(EDVP(pel->sons[3],v),EDVP(pel,u))
}

void vect_cubic_bubble1_rhs_restrict(tGrid,f,g)
GRID *tGrid;      /*  coarse grid  */
INT f, g;         /*  f ... fine grid; g ... coarse grid  */
{
   ELEMENT *pel, *e0, *e1, *e2, *e3;;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      e0 = pel->sons[0];
      e1 = pel->sons[1];
      e2 = pel->sons[2];
      e3 = pel->sons[3];
      EDV(pel,g,0) =0.5*(EDV(e0,f,0) + EDV(e1,f,0) + EDV(e2,f,0)) + EDV(e3,f,0);
      EDV(pel,g,1) =0.5*(EDV(e0,f,1) + EDV(e1,f,1) + EDV(e2,f,1)) + EDV(e3,f,1);
   }
}

#else

void vect_prolong_bubble1(pel,u,v)
ELEMENT *pel; INT u, v;
{  eprintf("Error: vect_prolong_bubble1 not available.\n");  }

void vect_cubic_bubble1_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_cubic_bubble1_rhs_restrict not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void vect_p1c_prolongation(tGrid,u,v,bubble)
GRID *tGrid;       /*  coarse grid  */
INT u, v, bubble;  /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      n11 = n1->son;
      n22 = n2->son;
      n33 = n3->son;
      n12 = unif_middle(n11,n22);
      n13 = unif_middle(n11,n33);
      n23 = unif_middle(n22,n33);
      p1c_coefficients(ND(n1,u,0),ND(n2,u,0),ND(n3,u,0),
                       &ND(n11,v,0), &ND(n22,v,0), &ND(n33,v,0),
                       &ND(n12,v,0), &ND(n13,v,0), &ND(n23,v,0));
      p1c_coefficients(ND(n1,u,1),ND(n2,u,1),ND(n3,u,1),
                       &ND(n11,v,1), &ND(n22,v,1), &ND(n33,v,1),
                       &ND(n12,v,1), &ND(n13,v,1), &ND(n23,v,1));
      if (bubble)
         vect_prolong_bubble1(pel,u,v);
   }
}

void vect_p1c_rhs_restrict(tGrid,f,g)
GRID *tGrid;      /*  coarse grid  */
INT f, g;         /*  f ... fine grid; g ... coarse grid  */
{
   NODE *theNode;
   LINK *plink;

   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      SET7(NDD(theNode,g),0.)
      for (plink = theNode->son->start; plink != NULL; plink = plink->next)
         SET5(NDD(theNode,g),NDD(NBNODE(plink),f))
      SET23(NDD(theNode,g),NDD(theNode->son,f),NDD(theNode,g),0.5)
   }
}

#else

void vect_p1c_prolongation(tGrid,u,v,bubble)
GRID *tGrid; INT u, v, bubble;
{  eprintf("Error: vect_p1c_prolongation not available.\n");  }

void vect_p1c_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_p1c_rhs_restrict not available.\n");  }

#endif

void p1c_values_at_midpoints(a1,a2,a3,u12,u13,u23)
FLOAT a1, a2, a3, *u12, *u13, *u23;
{
   if (fabs(*u12) < 1.e-16)
      *u12 = 0.5*(a1+a2);
   if (fabs(*u13) < 1.e-16)
      *u13 = 0.5*(a1+a3);
   if (fabs(*u23) < 1.e-16)
      *u23 = 0.5*(a2+a3);
}

void p2c_values_at_midpoints(a1,a2,a3,a12,a13,a23,u12,u13,u23)
FLOAT a1, a2, a3, a12, a13, a23, *u12, *u13, *u23;
{
   if (fabs(*u12) < 1.e-16)
      *u12 = 0.5*(a1+a2) + 0.25*a12;
   if (fabs(*u13) < 1.e-16)
      *u13 = 0.5*(a1+a3) + 0.25*a13;
   if (fabs(*u23) < 1.e-16)
      *u23 = 0.5*(a2+a3) + 0.25*a23;
}

void p2c_coefficients(a1,a2,a3,a12,a13,a23,u1,u2,u3,u12,u13,u23,
                      u112,u221,u113,u331,u223,u332,u11,u22,u33)
FLOAT a1, a2, a3, a12, a13, a23, *u1, *u2, *u3, *u12, *u13, *u23,
      *u112, *u221, *u113, *u331, *u223, *u332, *u11, *u22, *u33;
{
   *u1 = a1;
   *u2 = a2;
   *u3 = a3;
   *u11 = 0.25*a23;
   *u22 = 0.25*a13;
   *u33 = 0.25*a12;
   if (fabs(*u12) < 1.e-16){
      *u112 = *u221 = *u33;
      *u12 = 0.5*(a1+a2) + *u33;
   } 
   if (fabs(*u13) < 1.e-16){
      *u113 = *u331 = *u22;
      *u13 = 0.5*(a1+a3) + *u22;
   } 
   if (fabs(*u23) < 1.e-16){
      *u223 = *u332 = *u11;
      *u23 = 0.5*(a2+a3) + *u11;
   } 
}

void p2c_restr_coeff(a1,a2,a3,a12,a13,a23,u1,u2,u3,u12,u13,u23)
FLOAT a1, a2, a3, a12, a13, a23, *u1, *u2, *u3, *u12, *u13, *u23;
{
   *u1 = a1;
   *u2 = a2;
   *u3 = a3;
   if (fabs(*u12) < 1.e-16)
      *u12 = 2.*(a12 + a12 - a1 - a2);
   if (fabs(*u13) < 1.e-16)
      *u13 = 2.*(a13 + a13 - a1 - a3);
   if (fabs(*u23) < 1.e-16)
      *u23 = 2.*(a23 + a23 - a2 - a3);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void vect_prolong_bubble(pel,f11,f22,f33,u,v)
ELEMENT *pel;
FACE *f11, *f22, *f33;
INT u, v;
{
   DOUBLE a[DIM];

   SET2(a,EDVP(pel,u),0.125)
   SET1(EDVP(pel->sons[0],v),a)
   SET1(EDVP(pel->sons[1],v),a)
   SET1(EDVP(pel->sons[2],v),a)
   SET8(EDVP(pel->sons[3],v),a)
   SET5(FDVP(f11,v),a)
   SET5(FDVP(f22,v),a)
   SET5(FDVP(f33,v),a)
}

void vect_bubble_restr(pel,n11,n22,n33,n12,n13,n23,v,u)
ELEMENT *pel;
NODE *n11, *n22, *n33, *n12, *n13, *n23;
INT v, u;
{
   FACE *f11, *f22, *f33;

   f11 = pel->sons[0]->f[0];
   f22 = pel->sons[1]->f[0];
   f33 = pel->sons[2]->f[0];
   EDV(pel,u,0) = EDV(pel->sons[3],v,0) + 
                  3.*(   ND(n11,v,0) +  ND(n22,v,0) +  ND(n33,v,0)  
                      - (ND(n12,v,0) +  ND(n13,v,0) +  ND(n23,v,0))
                      + FDV(f11,v,0) + FDV(f22,v,0) + FDV(f33,v,0) );
   EDV(pel,u,1) = EDV(pel->sons[3],v,1) + 
                  3.*(   ND(n11,v,1) +  ND(n22,v,1) +  ND(n33,v,1)  
                      - (ND(n12,v,1) +  ND(n13,v,1) +  ND(n23,v,1))
                      + FDV(f11,v,1) + FDV(f22,v,1) + FDV(f33,v,1) );
}

void vect_cubic_bubble_rhs_restrict(tGrid,f,g)
GRID *tGrid;      /*  coarse grid  */
INT f, g;         /*  f ... fine grid; g ... coarse grid  */
{
   FACE *f11, *f22, *f33;
   ELEMENT *pel, *e0, *e1, *e2, *e3;;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      f11 = pel->sons[0]->f[0];
      f22 = pel->sons[1]->f[0];
      f33 = pel->sons[2]->f[0];
      e0 = pel->sons[0];
      e1 = pel->sons[1];
      e2 = pel->sons[2];
      e3 = pel->sons[3];
      EDV(pel,g,0) = 0.125*(FDV(f11,f,0) + FDV(f22,f,0) + FDV(f33,f,0) + 
                        EDV(e0,f,0) + EDV(e1,f,0) + EDV(e2,f,0) - EDV(e3,f,0));
      EDV(pel,g,1) = 0.125*(FDV(f11,f,1) + FDV(f22,f,1) + FDV(f33,f,1) + 
                        EDV(e0,f,1) + EDV(e1,f,1) + EDV(e2,f,1) - EDV(e3,f,1));
   }
}

#else

void vect_prolong_bubble(pel,f11,f22,f33,u,v)
ELEMENT *pel; FACE *f11, *f22, *f33; INT u, v;
{  eprintf("Error: vect_prolong_bubble not available.\n");  }

void vect_bubble_restr(pel,n11,n22,n33,n12,n13,n23,v,u)
ELEMENT *pel; NODE *n11, *n22, *n33, *n12, *n13, *n23; INT v, u;
{  eprintf("Error: vect_bubble_restr not available.\n");  }

void vect_cubic_bubble_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_cubic_bubble_rhs_restrict not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void scalar_p2c_prolongation(tGrid,u,v)
GRID *tGrid;                /*  coarse grid  */
INT u, v;                   /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p2c_coefficients(NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u),
                       FD(pel->f[2],u),FD(pel->f[1],u),FD(pel->f[0],u),
                       &NDS(n11,v), &NDS(n22,v), &NDS(n33,v),
                       &NDS(n12,v), &NDS(n13,v), &NDS(n23,v),
                        &FD(f31,v),  &FD(f32,v),  &FD(f21,v),
                        &FD(f23,v),  &FD(f12,v),  &FD(f13,v),
                        &FD(f11,v),  &FD(f22,v),  &FD(f33,v));
   }
}

void scalar_p2c_restriction(tGrid,v,u)
GRID *tGrid;        /*  coarse grid  */
INT v, u;           /*  v ... fine grid; u ... coarse grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      n11 = n1->son;
      n22 = n2->son;
      n33 = n3->son;
      n12 = unif_middle(n11,n22);
      n13 = unif_middle(n11,n33);
      n23 = unif_middle(n22,n33);
      p2c_restr_coeff(NDS(n11,v), NDS(n22,v), NDS(n33,v),
                      NDS(n12,v), NDS(n13,v), NDS(n23,v),
                     &NDS(n1,u), &NDS(n2,u), &NDS(n3,u), 
                      &FD(f3,u),  &FD(f2,u),  &FD(f1,u));
   }
}

void scalar_p2c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid;          /*  coarse grid  */
INT f, g, t;          /*  f ... fine grid; g ... coarse grid  */
{
   NODE *theNode, *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f1, *f2, *f3, *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   LINK *plink;
   ELEMENT *pel;

   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      NDS(theNode,g) = 0.;
      for (plink = theNode->son->start; plink != NULL; plink = plink->next)
         NDS(theNode,g) += NDS(NBNODE(plink),f);
      NDS(theNode,g) = NDS(theNode->son,f) + NDS(theNode,g)*0.5;
   }

   set_value(tGrid,0.,g,0,Q_SF);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      if (fabs(FD(f1,g)) < 1.e-16)
         FD(f1,g) = FD(f12,f) + FD(f13,f) + NDS(n23,f);
      if (fabs(FD(f2,g)) < 1.e-16)
         FD(f2,g) = FD(f21,f) + FD(f23,f) + NDS(n13,f);
      if (fabs(FD(f3,g)) < 1.e-16)
         FD(f3,g) = FD(f31,f) + FD(f32,f) + NDS(n12,f);
   }
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      f11 = pel->sons[0]->f[0];
      f22 = pel->sons[1]->f[0];
      f33 = pel->sons[2]->f[0];
      FD(f1,g) += FD(f11,f);
      FD(f2,g) += FD(f22,f);
      FD(f3,g) += FD(f33,f);
   }
   mult(tGrid,0.25,g,g,t,Q_SF);
}

void scalar_p1c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid;
INT u, v;                   /*  u ... p1c; v ... p1nc  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      p1c_values_at_midpoints(
                          NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u),
                          &FD(pel->f[2],v),&FD(pel->f[1],v),&FD(pel->f[0],v));
}

void scalar_p2c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid;
INT u, v;                   /*  u ... p2c; v ... p1nc  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      p2c_values_at_midpoints(
                          NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u),
                           FD(pel->f[2],u), FD(pel->f[1],u), FD(pel->f[0],u),
                          &FD(pel->f[2],v),&FD(pel->f[1],v),&FD(pel->f[0],v));
}

#else

void scalar_p2c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p2c_prolongation not available.\n");  }

void scalar_p2c_restriction(tGrid,v,u)
GRID *tGrid; INT v, u;
{  eprintf("Error: scalar_p2c_restriction not available.\n");  }

void scalar_p2c_rhs_restrict(tGrid,f,g,t)
GRID *tGrid; INT f, g, t;
{  eprintf("Error: scalar_p2c_rhs_restrict not available.\n");  }

void scalar_p1c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p1c_to_p1nc_restriction not available.\n");  }

void scalar_p2c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p2c_to_p1nc_restriction not available.\n");  }

#endif

void q2c_restr_coeff(a1,a2,a3,a4,a12,a23,a34,a41,a1234,
                     u1,u2,u3,u4,u12,u23,u34,u41,u1234)
FLOAT a1, a2, a3, a4, a12, a23, a34, a41, a1234, 
      *u1, *u2, *u3, *u4, *u12, *u23, *u34, *u41, *u1234;
{
   *u1 = a1;
   *u2 = a2;
   *u3 = a3;
   *u4 = a4;
   if (fabs(*u12) < 1.e-16)
      *u12 = a12  - 0.5*(a1 + a2);
   if (fabs(*u23) < 1.e-16)
      *u23 = a23  - 0.5*(a2 + a3);
   if (fabs(*u34) < 1.e-16)
      *u34 = a34  - 0.5*(a3 + a4);
   if (fabs(*u41) < 1.e-16)
      *u41 = a41  - 0.5*(a4 + a1);
   *u1234 = a1234 + 0.25*(a1+a2+a3+a4) - 0.5*(a12+a23+a34+a41);
}

#if (DIM == 2) && (ELEMENT_TYPE == CUBE) && (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)

void scalar_gq2c_restriction(tGrid,v,u)
GRID *tGrid;        /*  coarse grid  */
INT v, u;           /*  v ... fine grid; u ... coarse grid  */
{
   NODE *n1, *n2, *n3, *n4, *n11, *n22, *n33, *n44, 
        *n12, *n23, *n34, *n41, *n1234;
   FACE *f1, *f2, *f3, *f4;
   LINK *pl1, *pl2;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      FACES_OF_4ELEMENT(f1,f2,f3,f4,pel);
      n11 = n1->son;
      n22 = n2->son;
      n33 = n3->son;
      n44 = n4->son;
      n12 = unif_middle(n11,n22);
      n23 = unif_middle(n22,n33);
      n34 = unif_middle(n33,n44);
      n41 = unif_middle(n44,n11);
      n1234 = NULL;
      for (pl1 = START(n11); n1234 == NULL; pl1 = pl1->next)
         for (pl2 = START(n33); pl2 && n1234 == NULL; pl2 = pl2->next)
            if (NBNODE(pl1) == NBNODE(pl2))
               n1234 = NBNODE(pl1);
      q2c_restr_coeff(NDMV(n11,v,0), NDMV(n22,v,0), NDMV(n33,v,0), NDMV(n44,v,0),
                      NDMV(n12,v,0), NDMV(n23,v,0), NDMV(n34,v,0), NDMV(n41,v,0), 
                      NDMV(n1234,v,0),
                     &NDMV(n1,u,0), &NDMV(n2,u,0), &NDMV(n3,u,0), &NDMV(n4,u,0),
                     &FDMV(f1,u,0), &FDMV(f2,u,0), &FDMV(f3,u,0), &FDMV(f4,u,0),
                     &EDMV(pel,u,0));
   }
}

void scalar_gq2x4c_to_gq2c_restriction(tGrid,v,u)
GRID *tGrid;
INT v, u;           /*  v ... gq2x4c; u ... gq2c  */
{
   NODE *n1, *n2, *n3, *n4;
   FACE *f1, *f2, *f3, *f4;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      FACES_OF_4ELEMENT(f1,f2,f3,f4,pel);
      q2c_restr_coeff(NDMV(n1,v,0), NDMV(n2,v,0), NDMV(n3,v,0), NDMV(n4,v,0),
                      NDMV(f1,v,1), NDMV(f2,v,1), NDMV(f3,v,1), NDMV(f4,v,1), 
                      EDMV(pel,v,8),
                     &NDMV(n1,u,0), &NDMV(n2,u,0), &NDMV(n3,u,0), &NDMV(n4,u,0),
                     &FDMV(f1,u,0), &FDMV(f2,u,0), &FDMV(f3,u,0), &FDMV(f4,u,0),
                     &EDMV(pel,u,0));
   }
}

#else

void scalar_gq2c_restriction(tGrid,v,u)
GRID *tGrid; INT v, u;
{  eprintf("Error: scalar_gq2c_restriction not available.\n");  }

void scalar_gq2x4c_to_gq2c_restriction(tGrid,v,u)
GRID *tGrid; INT v, u;
{  eprintf("Error: scalar_gq2x4c_to_gq2c_restriction not available.\n");  }

#endif

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)

void transfer_data_from_3el_to_p2c(tGrid,u)
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   FACE *f0i, *f1i, *f2i;
   FLOAT r, s;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NDMV(pel->n[0]->son,u,0) = NDMV(pel->n[0],u,0);
      NDMV(pel->n[1]->son,u,0) = NDMV(pel->n[1],u,0);
      NDMV(pel->n[2]->son,u,0) = NDMV(pel->n[2],u,0);
      FDMV(pel->f[0]->sons[0],u,0) = FDMV(pel->f[0],u,0);
      FDMV(pel->f[1]->sons[0],u,0) = FDMV(pel->f[1],u,0);
      FDMV(pel->f[2]->sons[0],u,0) = FDMV(pel->f[2],u,0);
      r = NDMV(pel->n[0],u,0) + NDMV(pel->n[1],u,0) + NDMV(pel->n[2],u,0);
      s = FDMV(pel->f[0],u,0) + FDMV(pel->f[1],u,0) + FDMV(pel->f[2],u,0);
      NDMV(pel->sons[0]->n[2],u,0) = (3.*r + s)/9. + EDMV(pel,u,3);
      s *= 2./9;
      if (pel->sons[2]->n[0]->index == pel->n[0]->son->index){
         f0i = pel->sons[2]->f[1];   
         f1i = pel->sons[2]->f[0];   
      }
      else{
         f0i = pel->sons[2]->f[0];   
         f1i = pel->sons[2]->f[1];   
      }
      if (pel->sons[0]->n[0]->index == pel->n[1]->son->index)
         f2i = pel->sons[0]->f[0];
      else
         f2i = pel->sons[0]->f[1];
      FDMV(f0i,u,0) = EDMV(pel,u,0) - FDMV(pel->f[0],u,0)/3. + s;
      FDMV(f1i,u,0) = EDMV(pel,u,1) - FDMV(pel->f[1],u,0)/3. + s;
      FDMV(f2i,u,0) = EDMV(pel,u,2) - FDMV(pel->f[2],u,0)/3. + s;
   }
}

void transfer_data_from_4el_to_q1c(tGrid,u)
GRID *tGrid;
INT u;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NDMV(pel->n[0]->son,u,0) = NDMV(pel->n[0],u,0);
      NDMV(pel->n[1]->son,u,0) = NDMV(pel->n[1],u,0);
      NDMV(pel->n[2]->son,u,0) = NDMV(pel->n[2],u,0);
      NDMV(pel->n[3]->son,u,0) = NDMV(pel->n[3],u,0);
      NDMV(pel->sons[0]->n[0],u,0) = FDMV(pel->sons[0]->f[0]->father,u,0);
      NDMV(pel->sons[1]->n[0],u,0) = FDMV(pel->sons[1]->f[0]->father,u,0);
      NDMV(pel->sons[2]->n[0],u,0) = FDMV(pel->sons[2]->f[0]->father,u,0);
      NDMV(pel->sons[3]->n[0],u,0) = FDMV(pel->sons[3]->f[0]->father,u,0);
      NDMV(pel->sons[0]->n[3],u,0) = EDMV(pel,u,0);
   }
}

void transfer_data_from_4el_to_q2c(tGrid,u)
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   FLOAT f[4][2];
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NDMV(pel->n[0]->son,u,0) = NDMV(pel->n[0],u,0);
      NDMV(pel->n[1]->son,u,0) = NDMV(pel->n[1],u,0);
      NDMV(pel->n[2]->son,u,0) = NDMV(pel->n[2],u,0);
      NDMV(pel->n[3]->son,u,0) = NDMV(pel->n[3],u,0);
      NDMV(pel->sons[0]->n[0],u,0) = FDMV(pel->sons[0]->f[0]->father,u,1);
      NDMV(pel->sons[1]->n[0],u,0) = FDMV(pel->sons[1]->f[0]->father,u,1);
      NDMV(pel->sons[2]->n[0],u,0) = FDMV(pel->sons[2]->f[0]->father,u,1);
      NDMV(pel->sons[3]->n[0],u,0) = FDMV(pel->sons[3]->f[0]->father,u,1);
      NDMV(pel->sons[0]->n[3],u,0) = EDMV(pel,u,8);
      for (i = 0; i < 4; i++){
            j = i + 1;
            if (j == 4)
               j = 0;
            if (pel->n[i]->index < pel->n[j]->index){
               f[i][0] = FDMV(pel->f[i],u,0);
               f[i][1] = FDMV(pel->f[i],u,2);
            }
            else{
               f[i][0] = FDMV(pel->f[i],u,2);
               f[i][1] = FDMV(pel->f[i],u,0);
            }
         }
      FDMV(pel->sons[0]->f[0],u,0) = f[3][1];
      FDMV(pel->sons[0]->f[1],u,0) = f[0][0];
      FDMV(pel->sons[1]->f[0],u,0) = f[0][1];
      FDMV(pel->sons[1]->f[1],u,0) = f[1][0];
      FDMV(pel->sons[2]->f[0],u,0) = f[1][1];
      FDMV(pel->sons[2]->f[1],u,0) = f[2][0];
      FDMV(pel->sons[3]->f[0],u,0) = f[2][1];
      FDMV(pel->sons[3]->f[1],u,0) = f[3][0];
      FDMV(pel->sons[0]->f[2],u,0) = EDMV(pel,u,0);
      FDMV(pel->sons[1]->f[2],u,0) = EDMV(pel,u,1);
      FDMV(pel->sons[2]->f[2],u,0) = EDMV(pel,u,2);
      FDMV(pel->sons[3]->f[2],u,0) = EDMV(pel,u,3);
      EDMV(pel->sons[0],u,0) = EDMV(pel,u,4);
      EDMV(pel->sons[1],u,0) = EDMV(pel,u,5);
      EDMV(pel->sons[2],u,0) = EDMV(pel,u,6);
      EDMV(pel->sons[3],u,0) = EDMV(pel,u,7);
   }
}

#else

void transfer_data_from_3el_to_p2c(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: transfer_data_from_3el_to_p2c not available.\n");  }

void transfer_data_from_4el_to_q1c(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: transfer_data_from_4el_to_q1c not available.\n");  }

void transfer_data_from_4el_to_q2c(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: transfer_data_from_4el_to_q2c not available.\n");  }

#endif

void transfer_data_from_mult_el(tGrid,u,space)
GRID *tGrid;
INT u, space;
{
   if (space == GP2X3C)
      transfer_data_from_3el_to_p2c(tGrid,u);
   else if (space == GQ1X4C)
      transfer_data_from_4el_to_q1c(tGrid,u);
   else if (space == GQ2X4C)
      transfer_data_from_4el_to_q2c(tGrid,u);
   else
      eprintf("Error: transfer_data_from_mult_el not available.\n");
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vect_p2c_prolongation(tGrid,u,v,t,type,bubble)
GRID *tGrid;                /*  coarse grid  */
INT u, v, t, type, bubble;  /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p2c_coefficients(ND(pel->n[0],u,0),ND(pel->n[1],u,0),ND(pel->n[2],u,0),
                       FDV(pel->f[2],u,0),FDV(pel->f[1],u,0),FDV(pel->f[0],u,0),
                       &ND(n11,v,0), &ND(n22,v,0), &ND(n33,v,0),
                       &ND(n12,v,0), &ND(n13,v,0), &ND(n23,v,0),
                       &FDV(f31,v,0), &FDV(f32,v,0), &FDV(f21,v,0),
                       &FDV(f23,v,0), &FDV(f12,v,0), &FDV(f13,v,0),
                       &FDV(f11,v,0), &FDV(f22,v,0), &FDV(f33,v,0));
      p2c_coefficients(ND(pel->n[0],u,1),ND(pel->n[1],u,1),ND(pel->n[2],u,1),
                       FDV(pel->f[2],u,1),FDV(pel->f[1],u,1),FDV(pel->f[0],u,1),
                       &ND(n11,v,1), &ND(n22,v,1), &ND(n33,v,1),
                       &ND(n12,v,1), &ND(n13,v,1), &ND(n23,v,1),
                       &FDV(f31,v,1), &FDV(f32,v,1), &FDV(f21,v,1),
                       &FDV(f23,v,1), &FDV(f12,v,1), &FDV(f13,v,1),
                       &FDV(f11,v,1), &FDV(f22,v,1), &FDV(f33,v,1));
      if (bubble)
         vect_prolong_bubble(pel,f11,f22,f33,u,v);
   }
}

void vect_p2c_restriction(tGrid,v,u,bubble)
GRID *tGrid;        /*  coarse grid  */
INT v, u, bubble;   /*  v ... fine grid; u ... coarse grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      n11 = n1->son;
      n22 = n2->son;
      n33 = n3->son;
      n12 = unif_middle(n11,n22);
      n13 = unif_middle(n11,n33);
      n23 = unif_middle(n22,n33);
      p2c_restr_coeff(ND(n11,v,0), ND(n22,v,0), ND(n33,v,0),
                      ND(n12,v,0), ND(n13,v,0), ND(n23,v,0),
                     &ND(n1,u,0), &ND(n2,u,0), &ND(n3,u,0), 
                    &FDV(f3,u,0),&FDV(f2,u,0),&FDV(f1,u,0));
      p2c_restr_coeff(ND(n11,v,1), ND(n22,v,1), ND(n33,v,1),
                      ND(n12,v,1), ND(n13,v,1), ND(n23,v,1),
                     &ND(n1,u,1), &ND(n2,u,1), &ND(n3,u,1), 
                    &FDV(f3,u,1),&FDV(f2,u,1),&FDV(f1,u,1));
      if (bubble)
         vect_bubble_restr(pel,n11,n22,n33,n12,n13,n23,v,u);
   }
}

void vect_p2c_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid;          /*  coarse grid  */
INT f, g, t, type;    /*  f ... fine grid; g ... coarse grid  */
{
   NODE *theNode, *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f1, *f2, *f3, *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   LINK *plink;
   ELEMENT *pel;

   for (theNode=FIRSTNODE(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      SET7(NDD(theNode,g),0.)
      for (plink = theNode->son->start; plink != NULL; plink = plink->next)
         SET5(NDD(theNode,g),NDD(NBNODE(plink),f))
      SET23(NDD(theNode,g),NDD(theNode->son,f),NDD(theNode,g),0.5)
   }

   set_value(tGrid,0.,g,0,Q_VF);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      if (fabs(FDV(f1,g,0)) < 1.e-16)
         FDV(f1,g,0) = FDV(f12,f,0) + FDV(f13,f,0) + ND(n23,f,0);
      if (fabs(FDV(f2,g,0)) < 1.e-16)
         FDV(f2,g,0) = FDV(f21,f,0) + FDV(f23,f,0) + ND(n13,f,0);
      if (fabs(FDV(f3,g,0)) < 1.e-16)
         FDV(f3,g,0) = FDV(f31,f,0) + FDV(f32,f,0) + ND(n12,f,0);
      if (fabs(FDV(f1,g,1)) < 1.e-16)
         FDV(f1,g,1) = FDV(f12,f,1) + FDV(f13,f,1) + ND(n23,f,1);
      if (fabs(FDV(f2,g,1)) < 1.e-16)
         FDV(f2,g,1) = FDV(f21,f,1) + FDV(f23,f,1) + ND(n13,f,1);
      if (fabs(FDV(f3,g,1)) < 1.e-16)
         FDV(f3,g,1) = FDV(f31,f,1) + FDV(f32,f,1) + ND(n12,f,1);
   }
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      f11 = pel->sons[0]->f[0];
      f22 = pel->sons[1]->f[0];
      f33 = pel->sons[2]->f[0];
      FDV(f1,g,0) += FDV(f11,f,0);
      FDV(f2,g,0) += FDV(f22,f,0);
      FDV(f3,g,0) += FDV(f33,f,0);
      FDV(f1,g,1) += FDV(f11,f,1);
      FDV(f2,g,1) += FDV(f22,f,1);
      FDV(f3,g,1) += FDV(f33,f,1);
   }
   mult(tGrid,0.25,g,g,t,Q_VF);
}

void vect_p1c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid;
INT u, v;                   /*  u ... p1c; v ... p1nc  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      p1c_values_at_midpoints(
                     ND(pel->n[0],u,0),  ND(pel->n[1],u,0),  ND(pel->n[2],u,0),
                   &FDV(pel->f[2],v,0),&FDV(pel->f[1],v,0),&FDV(pel->f[0],v,0));
      p1c_values_at_midpoints(
                     ND(pel->n[0],u,1),  ND(pel->n[1],u,1),  ND(pel->n[2],u,1),
                   &FDV(pel->f[2],v,1),&FDV(pel->f[1],v,1),&FDV(pel->f[0],v,1));
   }
}

void vect_p2c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid;
INT u, v;                   /*  u ... p2c; v ... p1nc  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      p2c_values_at_midpoints(
                     ND(pel->n[0],u,0),  ND(pel->n[1],u,0),  ND(pel->n[2],u,0),
                    FDV(pel->f[2],u,0), FDV(pel->f[1],u,0), FDV(pel->f[0],u,0),
                   &FDV(pel->f[2],v,0),&FDV(pel->f[1],v,0),&FDV(pel->f[0],v,0));
      p2c_values_at_midpoints(
                     ND(pel->n[0],u,1),  ND(pel->n[1],u,1),  ND(pel->n[2],u,1),
                    FDV(pel->f[2],u,1), FDV(pel->f[1],u,1), FDV(pel->f[0],u,1),
                   &FDV(pel->f[2],v,1),&FDV(pel->f[1],v,1),&FDV(pel->f[0],v,1));
   }
}

#else

void vect_p2c_prolongation(tGrid,u,v,t,type,bubble)
GRID *tGrid; INT u, v, t, type, bubble;
{  eprintf("Error: vect_p2c_prolongation not available.\n");  }

void vect_p2c_restriction(tGrid,v,u,bubble)
GRID *tGrid; INT v, u, bubble;
{  eprintf("Error: vect_p2c_restriction not available.\n");  }

void vect_p2c_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid; INT f, g, t, type; 
{  eprintf("Error: vect_p2c_rhs_restrict not available.\n");  }

void vect_p1c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vect_p1c_to_p1nc_restriction not available.\n");  }

void vect_p2c_to_p1nc_restriction(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vect_p2c_to_p1nc_restriction not available.\n");  }

#endif

#if (N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA)

void p2c_coefficients_BD(f0,n1,n2,u,v)
FACE *f0;
NODE *n1, *n2;
INT u, v;        /*  u ... coarse grid; v ... fine grid  */
{
	NODE *n11, *n22, *n12;
	FACE *f01, *f02;

	n11 = n1->son;
	n22 = n2->son;
	n12 = unif_middle(n11,n22);
	get_sons_of_face(f0,&f01,&f02,n11,n22);
	FDBD(f01,v) = FDBD(f02,v) = FDBD(f0,u)/4.;
	if (IS_FN(n1)) {
		NDSBD(n11,v) = NDSBD(n1,u);
		if (IS_FN(n2)) {
			NDSBD(n22,v) = NDSBD(n2,u);
			NDSBD(n12,v) = FDBD(f01,v) + (NDSBD(n1,u) + NDSBD(n2,u))/2.;
		} else
			NDSBD(n12,v) = FDBD(f01,v) + NDSBD(n1,u)/2.;
	} else {
		NDSBD(n22,v) = NDSBD(n2,u);
		NDSBD(n12,v) = FDBD(f01,v) + NDSBD(n2,u)/2.;
	}
}

void scalar_p2c_prolongation_BD(tGrid,u,v)
GRID *tGrid;     /*  coarse grid  */
INT u, v;        /*  u ... coarse grid; v ... fine grid  */
{
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	ELEMENT *pel;

	for (pel = FIRSTELEMENT(tGrid); pel != FDBE(tGrid); pel = pel->succ) {
		NODES_OF_ELEMENT(n0,n1,n2,pel);
		FACES_OF_ELEMENT(f0,f1,f2,pel);
		if (IS_ZNF(f0))
			p2c_coefficients_BD(f0,n1,n2,u,v);
		else if (IS_ZNF(f1))
			p2c_coefficients_BD(f1,n0,n2,u,v);
		else if (IS_ZNF(f2))
			p2c_coefficients_BD(f2,n0,n1,u,v);
	}
}

void p2c_coeff_restrict_BD(f0,n1,n2,f,g)
FACE *f0;
NODE *n1, *n2;
INT f, g;	      /*  f ... fine grid; g ... coarse grid  */
{
	NODE *n11, *n22, *n12;

	n11 = n1->son;
	n22 = n2->son;
	n12 = unif_middle(n11,n22);
	if (IS_FN(n1)) {
		NDSBD(n1,g) = NDSBD(n11,f);
		if (IS_FN(n2)) {
			NDSBD(n2,g) = NDSBD(n22,f);
			FDBD(f0,g) = 4.*NDSBD(n12,f) - 2.*(NDSBD(n11,f) + NDSBD(n22,f));
		} else
			FDBD(f0,g) = 4.*NDSBD(n12,f) - 2.*NDSBD(n11,f);
	} else {
		NDSBD(n2,g) = NDSBD(n22,f);
		FDBD(f0,g) = 4.*NDSBD(n12,f) - 2.*NDSBD(n22,f);
	}
}

void scalar_p2c_rhs_restrict_BD(tGrid,f,g)
GRID *tGrid;      /*  coarse grid  */
INT f, g;	      /*  f ... fine grid; g ... coarse grid  */
{
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	ELEMENT *pel;

	for (pel = FIRSTELEMENT(tGrid); pel != FDBE(tGrid); pel = pel->succ) {
		NODES_OF_ELEMENT(n0,n1,n2,pel);
		FACES_OF_ELEMENT(f0,f1,f2,pel);
		if (IS_ZNF(f0))
			p2c_coeff_restrict_BD(f0,n1,n2,f,g);
		else if (IS_ZNF(f1))
			p2c_coeff_restrict_BD(f1,n0,n2,f,g);
		else if (IS_ZNF(f2))
			p2c_coeff_restrict_BD(f2,n0,n1,f,g);
	}
}

#else

void scalar_p2c_prolongation_BD(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p2c_prolongation_BD not available.\n");  }

void scalar_p2c_rhs_restrict_BD(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: scalar_p2c_rhs_restrict_BD not available.\n");  }

#endif

void p1nc_coefficients(a1,a2,a3,u11,u12,u13,u22,u21,u23,u33,u31,u32)
FLOAT a1, a2, a3, *u11, *u12, *u13, *u22, *u21, *u23, *u33, *u31, *u32;
{
   *u12 += 0.5*a1 - 0.25*(a2 - a3);
   *u13 += 0.5*a1 - 0.25*(a3 - a2);
   *u21 += 0.5*a2 - 0.25*(a1 - a3);
   *u23 += 0.5*a2 - 0.25*(a3 - a1);
   *u31 += 0.5*a3 - 0.25*(a1 - a2);
   *u32 += 0.5*a3 - 0.25*(a2 - a1);
   *u11  = 0.5*(a2 + a3);
   *u22  = 0.5*(a1 + a3);
   *u33  = 0.5*(a1 + a2);
}

FLOAT add_restrict_p1nc(fa1,fa2,fa3,f12,f13,f23,f32,f21,f31,f22,f33)
FACE *fa1, *fa2, *fa3;
FLOAT f12,f13,f23,f32,f21,f31,f22,f33;
{
   FLOAT s1=f12+f13, s2=f23-f21, s3=f32-f31;

   if (IS_BF(fa1))
      s1 += s1;
   if (IS_BF(fa2))
      s2 += s2;
   if (IS_BF(fa3))
      s3 += s3;
   return( ( 2.*(s1 + f22 + f33) + s2 + s3 )*0.25 );
}

void add3_restrict_p1nc(fa1,fa2,fa3,f11,f12,f13,f22,f21,f23,f33,f31,f32,f1,f2,f3)
FACE *fa1, *fa2, *fa3;
FLOAT f11, f12, f13, f22, f21, f23, f33, f31, f32, *f1, *f2, *f3;
{
   *f1 += add_restrict_p1nc(fa1,fa2,fa3,f12,f13,f23,f32,f21,f31,f22,f33);
   *f2 += add_restrict_p1nc(fa2,fa1,fa3,f21,f23,f13,f31,f12,f32,f11,f33);
   *f3 += add_restrict_p1nc(fa3,fa1,fa2,f31,f32,f12,f21,f13,f23,f11,f22);
}

#if F_DATA & SCALAR_FACE_DATA

void scalar_p1nc_prolongation(tGrid,u,v,t,type)  /* tGrid is the coarse grid */
GRID *tGrid;        /*  coarse grid  */
INT u, v, t, type;  /*  u ... coarse grid; v ... fine grid  */
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33, *pface;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p1nc_coefficients( FD(pel->f[0],u),FD(pel->f[1],u),FD(pel->f[2],u),
                        &FD(f11,v), &FD(f12,v), &FD(f13,v),
                        &FD(f22,v), &FD(f21,v), &FD(f23,v),
                        &FD(f33,v), &FD(f31,v), &FD(f32,v));
   }
   for (pface=FIRSTF(tGrid->finer); pface != NULL; pface=SUCC(pface))
      if (IS_BF(pface))
         FD(pface,v) += FD(pface,v);
}

void scalar_p1nc_restriction(tGrid,v,u,q,type)   /*  L2 projection  */
GRID *tGrid;        /*  coarse grid  */
INT v, u, q, type;  /*  v ... fine grid; u ... coarse grid; q ... aux. var.  */
{
   FACE *f1, *f2, *f3, *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33, *pface;
   ELEMENT *pel;
   FLOAT s, z;

   set_value(tGrid,0.,u,0,type);
   set_value(tGrid,0.,q,0,type);

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      s = VOLUME(pel);
      FD(f1,q) += s;
      FD(f2,q) += s;
      FD(f3,q) += s;
      s *= .25;
      z = s*.5;
      FD(f1,u) += s*(FD(f12,v)+FD(f13,v)+FD(f22,v)+FD(f33,v)) +
                  z*(FD(f23,v)+FD(f32,v)-FD(f21,v)-FD(f31,v));
      FD(f2,u) += s*(FD(f21,v)+FD(f23,v)+FD(f11,v)+FD(f33,v)) +
                  z*(FD(f13,v)+FD(f31,v)-FD(f12,v)-FD(f32,v));
      FD(f3,u) += s*(FD(f31,v)+FD(f32,v)+FD(f11,v)+FD(f22,v)) +
                  z*(FD(f12,v)+FD(f21,v)-FD(f13,v)-FD(f23,v));
   }
   for (pface=FIRSTF(tGrid); pface != NULL; pface=SUCC(pface))
      FD(pface,u) /= FD(pface,q);
}

void scalar_p1nc_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid;          /*  coarse grid  */
INT f, g, t, type;    /*  f ... fine grid; g ... coarse grid   */
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      add3_restrict_p1nc( pel->f[0],pel->f[1],pel->f[2],
                          FD(f11,f), FD(f12,f), FD(f13,f),
                          FD(f22,f), FD(f21,f), FD(f23,f),
                          FD(f33,f), FD(f31,f), FD(f32,f),
                         &FD(pel->f[0],g), &FD(pel->f[1],g), &FD(pel->f[2],g));
   }
}

#else  /*  if !(F_DATA & SCALAR_FACE_DATA)  */

void scalar_p1nc_prolongation(tGrid,u,v,t,type)  
GRID *tGrid; INT u, v, t, type;  
{  eprintf("Error: scalar_p1nc_prolongation not available.\n");  }

void scalar_p1nc_restriction(tGrid,v,u,q,type)
GRID *tGrid; INT v, u, q, type;
{  eprintf("Error: scalar_p1nc_restriction not available.\n");  }
 
void scalar_p1nc_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid; INT f, g, t, type;
{  eprintf("Error: scalar_p1nc_rhs_restrict not available.\n");  }
 
#endif

#if F_DATA & VECTOR_FACE_DATA

void vect_p1nc_prolongation(tGrid,u,v,t,type)  /* tGrid is the coarse grid */
GRID *tGrid;        /*  coarse grid  */
INT u, v, t, type;  /*  u ... coarse grid; v ... fine grid  */
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33, *pface;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p1nc_coefficients(
              FDV(pel->f[0],u,0),FDV(pel->f[1],u,0),FDV(pel->f[2],u,0),
              &FDV(f11,v,0), &FDV(f12,v,0), &FDV(f13,v,0),
              &FDV(f22,v,0), &FDV(f21,v,0), &FDV(f23,v,0),
              &FDV(f33,v,0), &FDV(f31,v,0), &FDV(f32,v,0));
      p1nc_coefficients(
              FDV(pel->f[0],u,1),FDV(pel->f[1],u,1),FDV(pel->f[2],u,1),
              &FDV(f11,v,1), &FDV(f12,v,1), &FDV(f13,v,1),
              &FDV(f22,v,1), &FDV(f21,v,1), &FDV(f23,v,1),
              &FDV(f33,v,1), &FDV(f31,v,1), &FDV(f32,v,1));
   }
   for (pface=FIRSTF(tGrid->finer); pface != NULL; pface=SUCC(pface))
      if (IS_BF(pface))
         SET5(FDVP(pface,v),FDVP(pface,v))
}

void vect_p1nc_restriction(tGrid,v,u,q,type)   /*  L2 projection  */
GRID *tGrid;        /*  coarse grid  */
INT v, u, q, type;  /*  v ... fine grid; u ... coarse grid; q ... aux. var.  */
{
   FACE *f1, *f2, *f3, *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33, *pface;
   ELEMENT *pel;
   FLOAT s, z;

   set_value(tGrid,0.,u,0,type);
   set_value(tGrid,0.,q,0,type);

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      s = VOLUME(pel);
      FDV(f1,q,0) += s;
      FDV(f2,q,0) += s;
      FDV(f3,q,0) += s;
      s *= .25;
      z = s*.5;
      FDV(f1,u,0) += s*(FDV(f12,v,0)+FDV(f13,v,0)+FDV(f22,v,0)+FDV(f33,v,0)) +
                     z*(FDV(f23,v,0)+FDV(f32,v,0)-FDV(f21,v,0)-FDV(f31,v,0));
      FDV(f2,u,0) += s*(FDV(f21,v,0)+FDV(f23,v,0)+FDV(f11,v,0)+FDV(f33,v,0)) +
                     z*(FDV(f13,v,0)+FDV(f31,v,0)-FDV(f12,v,0)-FDV(f32,v,0));
      FDV(f3,u,0) += s*(FDV(f31,v,0)+FDV(f32,v,0)+FDV(f11,v,0)+FDV(f22,v,0)) +
                     z*(FDV(f12,v,0)+FDV(f21,v,0)-FDV(f13,v,0)-FDV(f23,v,0));
      FDV(f1,u,1) += s*(FDV(f12,v,1)+FDV(f13,v,1)+FDV(f22,v,1)+FDV(f33,v,1)) +
                     z*(FDV(f23,v,1)+FDV(f32,v,1)-FDV(f21,v,1)-FDV(f31,v,1));
      FDV(f2,u,1) += s*(FDV(f21,v,1)+FDV(f23,v,1)+FDV(f11,v,1)+FDV(f33,v,1)) +
                     z*(FDV(f13,v,1)+FDV(f31,v,1)-FDV(f12,v,1)-FDV(f32,v,1));
      FDV(f3,u,1) += s*(FDV(f31,v,1)+FDV(f32,v,1)+FDV(f11,v,1)+FDV(f22,v,1)) +
                     z*(FDV(f12,v,1)+FDV(f21,v,1)-FDV(f13,v,1)-FDV(f23,v,1));
   }
   for (pface=FIRSTF(tGrid); pface != NULL; pface=SUCC(pface))
      SET13(FDVP(pface,u),FDV(pface,q,0))
}

void vect_p1nc_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid;          /*  coarse grid  */
INT f, g, t, type;    /*  f ... fine grid; g ... coarse grid   */
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      add3_restrict_p1nc(pel->f[0],pel->f[1],pel->f[2],
         FDV(f11,f,0), FDV(f12,f,0), FDV(f13,f,0),
         FDV(f22,f,0), FDV(f21,f,0), FDV(f23,f,0),
         FDV(f33,f,0), FDV(f31,f,0), FDV(f32,f,0),
        &FDV(pel->f[0],g,0), &FDV(pel->f[1],g,0), &FDV(pel->f[2],g,0));
      add3_restrict_p1nc(pel->f[0],pel->f[1],pel->f[2],
         FDV(f11,f,1), FDV(f12,f,1), FDV(f13,f,1),
         FDV(f22,f,1), FDV(f21,f,1), FDV(f23,f,1),
         FDV(f33,f,1), FDV(f31,f,1), FDV(f32,f,1),
        &FDV(pel->f[0],g,1), &FDV(pel->f[1],g,1), &FDV(pel->f[2],g,1));
   }
}

#else  /*  if !(F_DATA & VECTOR_FACE_DATA)  */

void vect_p1nc_prolongation(tGrid,u,v,t,type)  
GRID *tGrid; INT u, v, t, type;  
{  eprintf("Error: vect_p1nc_prolongation not available.\n");  }

void vect_p1nc_restriction(tGrid,v,u,q,type)
GRID *tGrid; INT v, u, q, type;
{  eprintf("Error: vect_p1nc_restriction not available.\n");  }
 
void vect_p1nc_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid; INT f, g, t, type;
{  eprintf("Error: vect_p1nc_rhs_restrict not available.\n");  }
 
#endif

#if (N_DATA & NUMBER_OF_N_NEL) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void scalar_p1nc_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid; 
INT u, v;  /*  u ... p1nc; v ... p1c  */
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      NDS(n1,v) += -FD(f1,u) + FD(f2,u) + FD(f3,u);
      NDS(n2,v) +=  FD(f1,u) - FD(f2,u) + FD(f3,u);
      NDS(n3,v) +=  FD(f1,u) + FD(f2,u) - FD(f3,u);
   }
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDS(pnode,v) /= pnode->nel;
   set_value(tGrid,0.,v,STOP_IS_FIRST_INNER,Q_SN);
}

#else  /*  !((N_DATA & NUMBER_OF_N_NEL) && (N_DATA & SCALAR_NODE_DATA) && 
             (F_DATA & SCALAR_FACE_DATA))  */

void scalar_p1nc_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p1nc_to_p1c_prolongation not available.\n");  }

#endif

#if (N_DATA & NUMBER_OF_N_NEL) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vect_p1nc_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid; 
INT u, v;  /*  u ... p1nc; v ... p1c  */
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      ND(n1,v,0) += -FDV(f1,u,0) + FDV(f2,u,0) + FDV(f3,u,0);
      ND(n2,v,0) +=  FDV(f1,u,0) - FDV(f2,u,0) + FDV(f3,u,0);
      ND(n3,v,0) +=  FDV(f1,u,0) + FDV(f2,u,0) - FDV(f3,u,0);
      ND(n1,v,1) += -FDV(f1,u,1) + FDV(f2,u,1) + FDV(f3,u,1);
      ND(n2,v,1) +=  FDV(f1,u,1) - FDV(f2,u,1) + FDV(f3,u,1);
      ND(n3,v,1) +=  FDV(f1,u,1) + FDV(f2,u,1) - FDV(f3,u,1);
   }
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      SET13(NDD(pnode,v),pnode->nel)
   set_value(tGrid,0.,v,STOP_IS_FIRST_INNER,Q_VN);
}

void vect_p1nc_to_p2c_prolongation(tGrid,u,v)
GRID *tGrid; 
INT u, v;  /*  u ... p1nc; v ... p2c  */
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      ND(n1,v,0) += -FDV(f1,u,0) + FDV(f2,u,0) + FDV(f3,u,0);
      ND(n2,v,0) +=  FDV(f1,u,0) - FDV(f2,u,0) + FDV(f3,u,0);
      ND(n3,v,0) +=  FDV(f1,u,0) + FDV(f2,u,0) - FDV(f3,u,0);
      ND(n1,v,1) += -FDV(f1,u,1) + FDV(f2,u,1) + FDV(f3,u,1);
      ND(n2,v,1) +=  FDV(f1,u,1) - FDV(f2,u,1) + FDV(f3,u,1);
      ND(n3,v,1) +=  FDV(f1,u,1) + FDV(f2,u,1) - FDV(f3,u,1);
   }
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      SET13(NDD(pnode,v),pnode->nel)
   set_value(tGrid,0.,v,STOP_IS_FIRST_INNER,Q_VN);
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      if (fabs(FDV(f1,v,0)) < 1.e-16)
         SET24(FDVP(f1,v),FDVP(f1,u),NDD(n2,v),NDD(n3,v))
      if (fabs(FDV(f2,v,0)) < 1.e-16)
         SET24(FDVP(f2,v),FDVP(f2,u),NDD(n1,v),NDD(n3,v))
      if (fabs(FDV(f3,v,0)) < 1.e-16)
         SET24(FDVP(f3,v),FDVP(f3,u),NDD(n1,v),NDD(n2,v))
   }
}

void vect_p2c_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid;
INT f, g; /*  f ... p2c; g ... p1nc  */
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      SET6(NDD(n2,f),FDVP(f1,f))   
      SET6(NDD(n3,f),FDVP(f1,f))   
      SET6(NDD(n1,f),FDVP(f2,f))   
      SET6(NDD(n3,f),FDVP(f2,f))
      SET6(NDD(n1,f),FDVP(f3,f))
      SET6(NDD(n2,f),FDVP(f3,f))
      if (IS_BF(f1)){
         SET6(NDD(n2,f),FDVP(f1,f))  
         SET6(NDD(n3,f),FDVP(f1,f))   
      }
      if (IS_BF(f2)){
         SET6(NDD(n1,f),FDVP(f2,f))   
         SET6(NDD(n3,f),FDVP(f2,f))
      }
      if (IS_BF(f3)){
         SET6(NDD(n1,f),FDVP(f3,f))
         SET6(NDD(n2,f),FDVP(f3,f))
      }
   }
   set_value(tGrid,0.,f,STOP_IS_FIRST_INNER,Q_VN);
   mult(tGrid,4.,f,g,0,Q_VF);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      FDV(f1,g,0) += ND(n2,f,0)/n2->nel+ND(n3,f,0)/n3->nel - ND(n1,f,0)/n1->nel;
      FDV(f2,g,0) += ND(n1,f,0)/n1->nel+ND(n3,f,0)/n3->nel - ND(n2,f,0)/n2->nel;
      FDV(f3,g,0) += ND(n1,f,0)/n1->nel+ND(n2,f,0)/n2->nel - ND(n3,f,0)/n3->nel;
      FDV(f1,g,1) += ND(n2,f,1)/n2->nel+ND(n3,f,1)/n3->nel - ND(n1,f,1)/n1->nel;
      FDV(f2,g,1) += ND(n1,f,1)/n1->nel+ND(n3,f,1)/n3->nel - ND(n2,f,1)/n2->nel;
      FDV(f3,g,1) += ND(n1,f,1)/n1->nel+ND(n2,f,1)/n2->nel - ND(n3,f,1)/n3->nel;
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      SET5(NDD(n2,f),FDVP(f1,f))   
      SET5(NDD(n3,f),FDVP(f1,f))   
      SET5(NDD(n1,f),FDVP(f2,f))   
      SET5(NDD(n3,f),FDVP(f2,f))
      SET5(NDD(n1,f),FDVP(f3,f))
      SET5(NDD(n2,f),FDVP(f3,f))
      if (IS_BF(f1)){
         SET5(NDD(n2,f),FDVP(f1,f))  
         SET5(NDD(n3,f),FDVP(f1,f))   
      }
      if (IS_BF(f2)){
         SET5(NDD(n1,f),FDVP(f2,f))   
         SET5(NDD(n3,f),FDVP(f2,f))
      }
      if (IS_BF(f3)){
         SET5(NDD(n1,f),FDVP(f3,f))
         SET5(NDD(n2,f),FDVP(f3,f))
      }
   }
}

#else  /*  !((N_DATA & NUMBER_OF_N_NEL) && (N_DATA & VECTOR_NODE_DATA) && 
             (F_DATA & VECTOR_FACE_DATA))  */

void vect_p1nc_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vect_p1nc_to_p1c_prolongation not available.\n");  }

void vect_p1nc_to_p2c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vect_p1nc_to_p2c_prolongation not available.\n");  }

void vect_p2c_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_p2c_to_p1nc_rhs_restrict not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vect_p2c_to_p1nc_rhs_restriction(tGrid,v,u,q)   /*  L2 projection  */
GRID *tGrid;
INT v, u, q;  /*  v ... p2c; u ... p1nc; q ... aux. var.  */
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3, *pface;
   ELEMENT *pel;
   FLOAT r, s;

   set_value(tGrid,0.,u,0,Q_VF);
   set_value(tGrid,0.,q,0,Q_VF);

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      s = VOLUME(pel);
      FDV(f1,q,0) += s;
      FDV(f2,q,0) += s;
      FDV(f3,q,0) += s;
      s *= .05;
      r = 10.*( ND(n1,v,0) +  ND(n2,v,0) +  ND(n3,v,0)) +
               FDV(f1,v,0) + FDV(f2,v,0) + FDV(f3,v,0);
      FDV(f1,u,0) += s*(r - 10.*ND(n1,v,0) + FDV(f1,v,0) + FDV(f1,v,0));
      FDV(f2,u,0) += s*(r - 10.*ND(n2,v,0) + FDV(f2,v,0) + FDV(f2,v,0));
      FDV(f3,u,0) += s*(r - 10.*ND(n3,v,0) + FDV(f3,v,0) + FDV(f3,v,0));
      r = 10.*( ND(n1,v,1) +  ND(n2,v,1) +  ND(n3,v,1)) +
               FDV(f1,v,1) + FDV(f2,v,1) + FDV(f3,v,1);
      FDV(f1,u,1) += s*(r - 10.*ND(n1,v,1) + FDV(f1,v,1) + FDV(f1,v,1));
      FDV(f2,u,1) += s*(r - 10.*ND(n2,v,1) + FDV(f2,v,1) + FDV(f2,v,1));
      FDV(f3,u,1) += s*(r - 10.*ND(n3,v,1) + FDV(f3,v,1) + FDV(f3,v,1));
   }
   for (pface=FIRSTF(tGrid); pface != NULL; pface=SUCC(pface))
      SET13(FDVP(pface,u),FDV(pface,q,0))
}

#else  /*  !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

void vect_p2c_to_p1nc_rhs_restriction(tGrid,v,u,q) 
GRID *tGrid; INT v, u, q;
{  eprintf("Error: vect_p2c_to_p1nc_rhs_restriction not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void set_bubbles_for_vect_p1nc_to_p2cb_prolongation(tGrid,u,v)
GRID *tGrid; 
INT u, v;  /*  u ... p1nc; v ... p2c  */
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   set_value(tGrid,0.,v,STOP_IS_FIRST_INNER,Q_VF);
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      EDV(pel,v,0) = 9.*( (FDV(f1,u,0) + FDV(f2,u,0) + FDV(f3,u,0))
                        - ( ND(n1,v,0) +  ND(n2,v,0) +  ND(n3,v,0)) )
                      -3.*(FDV(f1,v,0) + FDV(f2,v,0) + FDV(f3,v,0));
      EDV(pel,v,1) = 9.*( (FDV(f1,u,1) + FDV(f2,u,1) + FDV(f3,u,1))
                        - ( ND(n1,v,1) +  ND(n2,v,1) +  ND(n3,v,1)) )
                      -3.*(FDV(f1,v,1) + FDV(f2,v,1) + FDV(f3,v,1));
   }
}

#else

void set_bubbles_for_vect_p1nc_to_p2cb_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: set_bubbles_for_vect_p1nc_to_p2cb_prolongation not available.\n");  }

#endif

#if (N_DATA & NUMBER_OF_N_NEL) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void vect_p2cb_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid;
INT f, g; /*  f ... p2c; g ... p1nc  */
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      SET6(NDD(n2,f),FDVP(f1,f))   
      SET6(NDD(n3,f),FDVP(f1,f))   
      SET6(NDD(n1,f),FDVP(f2,f))   
      SET6(NDD(n3,f),FDVP(f2,f))
      SET6(NDD(n1,f),FDVP(f3,f))
      SET6(NDD(n2,f),FDVP(f3,f))
      if (IS_BF(f1)){
         SET6(NDD(n2,f),FDVP(f1,f))  
         SET6(NDD(n3,f),FDVP(f1,f))   
      }
      if (IS_BF(f2)){
         SET6(NDD(n1,f),FDVP(f2,f))   
         SET6(NDD(n3,f),FDVP(f2,f))
      }
      if (IS_BF(f3)){
         SET6(NDD(n1,f),FDVP(f3,f))
         SET6(NDD(n2,f),FDVP(f3,f))
      }
      SET2(EDVP(pel,f),EDVP(pel,f),3.)
      SET5(NDD(n1,f),EDVP(pel,f))
      SET5(NDD(n2,f),EDVP(pel,f))
      SET5(NDD(n3,f),EDVP(pel,f))
   }
   set_value(tGrid,0.,f,STOP_IS_FIRST_INNER,Q_VN);
   mult(tGrid,4.,f,g,0,Q_VF);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      FDV(f1,g,0) += ND(n2,f,0)/n2->nel+ND(n3,f,0)/n3->nel - ND(n1,f,0)/n1->nel
                     - EDV(pel,f,0);
      FDV(f2,g,0) += ND(n1,f,0)/n1->nel+ND(n3,f,0)/n3->nel - ND(n2,f,0)/n2->nel 
                     - EDV(pel,f,0);
      FDV(f3,g,0) += ND(n1,f,0)/n1->nel+ND(n2,f,0)/n2->nel - ND(n3,f,0)/n3->nel 
                     - EDV(pel,f,0);
      FDV(f1,g,1) += ND(n2,f,1)/n2->nel+ND(n3,f,1)/n3->nel - ND(n1,f,1)/n1->nel 
                     - EDV(pel,f,1);
      FDV(f2,g,1) += ND(n1,f,1)/n1->nel+ND(n3,f,1)/n3->nel - ND(n2,f,1)/n2->nel 
                     - EDV(pel,f,1);
      FDV(f3,g,1) += ND(n1,f,1)/n1->nel+ND(n2,f,1)/n2->nel - ND(n3,f,1)/n3->nel 
                     - EDV(pel,f,1);
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      SET5(NDD(n2,f),FDVP(f1,f))   
      SET5(NDD(n3,f),FDVP(f1,f))   
      SET5(NDD(n1,f),FDVP(f2,f))   
      SET5(NDD(n3,f),FDVP(f2,f))
      SET5(NDD(n1,f),FDVP(f3,f))
      SET5(NDD(n2,f),FDVP(f3,f))
      if (IS_BF(f1)){
         SET5(NDD(n2,f),FDVP(f1,f))  
         SET5(NDD(n3,f),FDVP(f1,f))   
      }
      if (IS_BF(f2)){
         SET5(NDD(n1,f),FDVP(f2,f))   
         SET5(NDD(n3,f),FDVP(f2,f))
      }
      if (IS_BF(f3)){
         SET5(NDD(n1,f),FDVP(f3,f))
         SET5(NDD(n2,f),FDVP(f3,f))
      }
      SET13(EDVP(pel,f),3.)
      SET6(NDD(n1,f),EDVP(pel,f))
      SET6(NDD(n2,f),EDVP(pel,f))
      SET6(NDD(n3,f),EDVP(pel,f))
   }
}

#else

void vect_p2cb_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_p2cb_to_p1nc_rhs_restrict not available.\n");  }

#endif

FLOAT average_for_son_face(a1,a2,a3,c3)
FLOAT a1, a2, a3, c3;
{
   return( (a1-a2)*0.25 + a3*0.5 + c3*0.03125 );
}

FLOAT average_for_inner_face(a1,a2,c1,c2)
FLOAT a1, a2, c1, c2;
{
   return( (a1+a2)*0.5 + (c1 + c2)/48. );
}

FLOAT moment_for_son_face(a1,a2,a3,c1,s12)
FLOAT a1, a2, a3, c1, s12;
{
   return( a1*15. + (a3-a2)*10. + c1 - s12*30. );
}

FLOAT moment_for_inner_face(a1,a2,c1,c2,c3,s33,n13,n23)
FLOAT a1, a2, c1, c2, c3, s33;
NODE *n13, *n23;
{
   return( (a1*160. + a2*80. + (c1+c2)*5. + c3 - s33*240.)*0.125*NINDI(n13,n23) );
}

void p1mod_coefficients(a1,a2,a3,b1,b2,b3,n1,n2,n3,n12,n13,n23,
                        u11,u12,u13,u22,u21,u23,u33,u31,u32,
                        v11,v12,v13,v22,v21,v23,v33,v31,v32)
FLOAT a1, a2, a3, b1, b2, b3,
      *u11, *u12, *u13, *u22, *u21, *u23, *u33, *u31, *u32,
      *v11, *v12, *v13, *v22, *v21, *v23, *v33, *v31, *v32;
NODE *n1, *n2, *n3, *n12, *n13, *n23;
{
   FLOAT c1, c2, c3, s12, s13, s21, s23, s31, s32;

   c1 = 10.*(a2 - a3) + b1*NINDI(n3,n2);
   c2 = 10.*(a1 - a3) + b2*NINDI(n3,n1);
   c3 = 10.*(a1 - a2) + b3*NINDI(n2,n1);
   
   *u12 += s12 = average_for_son_face(a3,a2,a1, c1);
   *u13 += s13 = average_for_son_face(a2,a3,a1,-c1);
   *u21 += s21 = average_for_son_face(a3,a1,a2, c2);
   *u23 += s23 = average_for_son_face(a1,a3,a2,-c2);
   *u31 += s31 = average_for_son_face(a2,a1,a3, c3);
   *u32 += s32 = average_for_son_face(a1,a2,a3,-c3);
   *u11  = average_for_inner_face(a2,a3, c2, c3);
   *u22  = average_for_inner_face(a1,a3, c1,-c3);
   *u33  = average_for_inner_face(a1,a2,-c1,-c2);
   *v12 += moment_for_son_face(a1,a2,a3, c1,s12);
   *v13 += moment_for_son_face(a1,a3,a2,-c1,s13);
   *v21 += moment_for_son_face(a2,a1,a3, c2,s21);
   *v23 += moment_for_son_face(a2,a3,a1,-c2,s23);
   *v31 += moment_for_son_face(a3,a1,a2, c3,s31);
   *v32 += moment_for_son_face(a3,a2,a1,-c3,s32);
   *v11 = moment_for_inner_face(a2,a3, c2, c3,-c1,*u11,n12,n13);
   *v22 = moment_for_inner_face(a1,a3, c1,-c3,-c2,*u22,n12,n23);
   *v33 = moment_for_inner_face(a1,a2,-c1,-c2,-c3,*u33,n13,n23);
}

FLOAT p1mod_int3(ne1,ne2,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) /* Simpson's rule */
FLOAT a1, a2, a3, b1, b2, b3, b[DIM2][DIM2];
NODE *ne1, *ne2, *n1, *n2, *n3;
{
   FLOAT x12[DIM];

   AVERAGE(ne1->myvertex->x,ne2->myvertex->x,x12);
   return( (p1mod_function(ne1->myvertex->x,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) +
            p1mod_function(ne2->myvertex->x,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) +
         4.*p1mod_function(x12,a1,a2,a3,b1,b2,b3,n1,n2,n3,b))/6. );
}

/* \int f = (7 f1 + 32 f112 + 12 f12 + 32 f122 + 7 f2)/90  exact for P4 */
/* \int p1mod * (l - 1/2) */
FLOAT p1mod_mint4(ne1,ne2,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)
FLOAT a1, a2, a3, b1, b2, b3, b[DIM2][DIM2];
NODE *ne1, *ne2, *n1, *n2, *n3;
{
   NODE *n;
   FLOAT x12[DIM], x112[DIM], x122[DIM];

   if ( ne1->index > ne2->index ){
      n = ne1;
      ne1 = ne2;
      ne2 = n;
   }
   
   AVERAGE(ne1->myvertex->x,ne2->myvertex->x,x12);
   AVERAGE(ne1->myvertex->x,x12,x112);
   AVERAGE(x12,ne2->myvertex->x,x122);
   return( 
    (7.*p1mod_function(ne1->myvertex->x,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) +
    16.*p1mod_function(x112,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) -
    16.*p1mod_function(x122,a1,a2,a3,b1,b2,b3,n1,n2,n3,b) -
     7.*p1mod_function(ne2->myvertex->x,a1,a2,a3,b1,b2,b3,n1,n2,n3,b))/180.);
}


void p1mod_coefficients_num(a1,a2,a3,b1,b2,b3,n1,n2,n3,n12,n13,n23,
                        u11,u12,u13,u22,u21,u23,u33,u31,u32,
                        v11,v12,v13,v22,v21,v23,v33,v31,v32)
FLOAT a1, a2, a3, b1, b2, b3,
      *u11, *u12, *u13, *u22, *u21, *u23, *u33, *u31, *u32,
      *v11, *v12, *v13, *v22, *v21, *v23, *v33, *v31, *v32;
NODE *n1, *n2, *n3, *n12, *n13, *n23;
{
   FLOAT b[DIM2][DIM2];

   barycentric_coordinates(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,b);

   *u12 += p1mod_int3(n2,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u13 += p1mod_int3(n3,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u21 += p1mod_int3(n1,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u23 += p1mod_int3(n3,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u31 += p1mod_int3(n1,n12,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u32 += p1mod_int3(n2,n12,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*0.5;
   *u11 = p1mod_int3(n12,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b);
   *u22 = p1mod_int3(n12,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b);
   *u33 = p1mod_int3(n13,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b);
   *v12 += p1mod_mint4(n2,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v13 += p1mod_mint4(n3,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v21 += p1mod_mint4(n1,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v23 += p1mod_mint4(n3,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v31 += p1mod_mint4(n1,n12,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v32 += p1mod_mint4(n2,n12,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*30.;
   *v11 = p1mod_mint4(n12,n13,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*60.;
   *v22 = p1mod_mint4(n12,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*60.;
   *v33 = p1mod_mint4(n13,n23,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)*60.;
}

FLOAT add_restrict_p1mod0(f12,f13,f23,f32,f21,f31,f11,f22,f33,
                          g23,g32,g21,g31,g22,g33,n12,n13,n23)
FLOAT f12,f13,f23,f32,f21,f31,f11,f22,f33,g23,g32,g21,g31,g22,g33;
NODE *n12, *n13, *n23;
{
   return( (f12 + f13)*.5 + (f21 + f31 - f23 - f32)/16. + 
           (10.*f11 + 7.*(f22 + f33))/24. + (g23 + g32 - g21 - g31)*1.875 +
           (g22*NINDI(n12,n23) + g33*NINDI(n13,n23))*3.75 );
}

FLOAT add_restrict_p1mod1(f12,f13,f22,f33,g12,g13,g11,n2,n3,n12,n13)
FLOAT f12,f13,f22,f33,g12,g13,g11;
NODE *n2, *n3, *n12, *n13;
{
   return( ( f22 - f33 + (f12 - f13)*1.5 + (g12 - g13)*3. 
                                   - 6.*g11*NINDI(n12,n13) )*NINDI(n3,n2)/48. );
}

void add3_restrict_p1mod(f11,f12,f13,f22,f21,f23,f33,f31,f32,
                         g11,g12,g13,g22,g21,g23,g33,g31,g32,
                         f1,f2,f3,g1,g2,g3,n1,n2,n3,n12,n13,n23)
FLOAT f11,f12,f13,f22,f21,f23,f33,f31,f32,g11,g12,g13,g22,g21,g23,g33,g31,g32,
      *f1, *f2, *f3, *g1, *g2, *g3;
NODE *n1, *n2, *n3, *n12, *n13, *n23;
{
   *f1 += add_restrict_p1mod0(f12,f13,f23,f32,f21,f31,f11,f22,f33,
                              g23,g32,g21,g31,g22,g33,n12,n13,n23);
   *f2 += add_restrict_p1mod0(f21,f23,f13,f31,f12,f32,f22,f11,f33,
                              g13,g31,g12,g32,g11,g33,n12,n23,n13);
   *f3 += add_restrict_p1mod0(f31,f32,f21,f12,f23,f13,f33,f11,f22,
                              g21,g12,g23,g13,g11,g22,n13,n23,n12);
   *g1 += add_restrict_p1mod1(f12,f13,f22,f33,g12,g13,g11,n2,n3,n12,n13);
   *g2 += add_restrict_p1mod1(f21,f23,f11,f33,g21,g23,g22,n1,n3,n12,n23);
   *g3 += add_restrict_p1mod1(f31,f32,f11,f22,g31,g32,g33,n1,n2,n13,n23);
}

#if F_DATA & VECTOR_FACE_DATA

void scalar_p1mod_prolongation(tGrid,u,v,t,type)  /* tGrid is the coarse grid */
GRID *tGrid;        /*  coarse grid  */
INT u, v, t, type;  /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p1mod_coefficients(
              FDV(pel->f[0],u,0),FDV(pel->f[1],u,0),FDV(pel->f[2],u,0),
              FDV(pel->f[0],u,1),FDV(pel->f[1],u,1),FDV(pel->f[2],u,1),
              n1,n2,n3,n12,n13,n23,
              &FDV(f11,v,0), &FDV(f12,v,0), &FDV(f13,v,0),
              &FDV(f22,v,0), &FDV(f21,v,0), &FDV(f23,v,0),
              &FDV(f33,v,0), &FDV(f31,v,0), &FDV(f32,v,0),
              &FDV(f11,v,1), &FDV(f12,v,1), &FDV(f13,v,1),
              &FDV(f22,v,1), &FDV(f21,v,1), &FDV(f23,v,1),
              &FDV(f33,v,1), &FDV(f31,v,1), &FDV(f32,v,1));
   }
}

void scalar_p1mod_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid;         /*  coarse grid  */
INT f, g, t, type;   /*  f ... fine grid; g ... coarse grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   set_value(tGrid,0.,g,0,type);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      add3_restrict_p1mod(
         FDV(f11,f,0), FDV(f12,f,0), FDV(f13,f,0),
         FDV(f22,f,0), FDV(f21,f,0), FDV(f23,f,0),
         FDV(f33,f,0), FDV(f31,f,0), FDV(f32,f,0),
         FDV(f11,f,1), FDV(f12,f,1), FDV(f13,f,1),
         FDV(f22,f,1), FDV(f21,f,1), FDV(f23,f,1),
         FDV(f33,f,1), FDV(f31,f,1), FDV(f32,f,1),
        &FDV(pel->f[0],g,0), &FDV(pel->f[1],g,0), &FDV(pel->f[2],g,0),
        &FDV(pel->f[0],g,1), &FDV(pel->f[1],g,1), &FDV(pel->f[2],g,1),
         n1,n2,n3,n12,n13,n23);
   }
}

#else  /*  if !(F_DATA & VECTOR_FACE_DATA)  */

void scalar_p1mod_prolongation(tGrid,u,v,t,type)  
GRID *tGrid; INT u, v, t;  
{  eprintf("Error: scalar_p1mod_prolongation not available.\n");  }

void scalar_p1mod_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid; INT f, g, t, type;
{  eprintf("Error: scalar_p1mod_rhs_restrict not available.\n");  }
 
#endif

#if F_DATA & DVECTOR_FACE_DATA

void vect_p1mod_prolongation(tGrid,u,v,t,type)  /* tGrid is the coarse grid */
GRID *tGrid;        /*  coarse grid  */
INT u, v, t, type;  /*  u ... coarse grid; v ... fine grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
/*      if(comp_elem_sons(pel,n11,n12,n13,f11,f21,f31) +
         comp_elem_sons(pel,n12,n22,n23,f12,f22,f32) +
         comp_elem_sons(pel,n13,n23,n33,f13,f23,f33) +
         comp_elem_sons(pel,n12,n23,n13,f33,f11,f22) != 1111)
         eprintf("Error***************************\n");
*/
      p1mod_coefficients(
              FDDV(pel->f[0],u,0,0),FDDV(pel->f[1],u,0,0),FDDV(pel->f[2],u,0,0),
              FDDV(pel->f[0],u,1,0),FDDV(pel->f[1],u,1,0),FDDV(pel->f[2],u,1,0),
              n1,n2,n3,n12,n13,n23,
              &FDDV(f11,v,0,0), &FDDV(f12,v,0,0), &FDDV(f13,v,0,0),
              &FDDV(f22,v,0,0), &FDDV(f21,v,0,0), &FDDV(f23,v,0,0),
              &FDDV(f33,v,0,0), &FDDV(f31,v,0,0), &FDDV(f32,v,0,0),
              &FDDV(f11,v,1,0), &FDDV(f12,v,1,0), &FDDV(f13,v,1,0),
              &FDDV(f22,v,1,0), &FDDV(f21,v,1,0), &FDDV(f23,v,1,0),
              &FDDV(f33,v,1,0), &FDDV(f31,v,1,0), &FDDV(f32,v,1,0));
      p1mod_coefficients(
              FDDV(pel->f[0],u,0,1),FDDV(pel->f[1],u,0,1),FDDV(pel->f[2],u,0,1),
              FDDV(pel->f[0],u,1,1),FDDV(pel->f[1],u,1,1),FDDV(pel->f[2],u,1,1),
              n1,n2,n3,n12,n13,n23,
              &FDDV(f11,v,0,1), &FDDV(f12,v,0,1), &FDDV(f13,v,0,1),
              &FDDV(f22,v,0,1), &FDDV(f21,v,0,1), &FDDV(f23,v,0,1),
              &FDDV(f33,v,0,1), &FDDV(f31,v,0,1), &FDDV(f32,v,0,1),
              &FDDV(f11,v,1,1), &FDDV(f12,v,1,1), &FDDV(f13,v,1,1),
              &FDDV(f22,v,1,1), &FDDV(f21,v,1,1), &FDDV(f23,v,1,1),
              &FDDV(f33,v,1,1), &FDDV(f31,v,1,1), &FDDV(f32,v,1,1));
   }
}

void vect_p1mod_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid;         /*  coarse grid  */
INT f, g, t, type;   /*  f ... fine grid; g ... coarse grid  */
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   ELEMENT *pel;

   set_value(tGrid,0.,g,0,type);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
               &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      add3_restrict_p1mod(
         FDDV(f11,f,0,0), FDDV(f12,f,0,0), FDDV(f13,f,0,0),
         FDDV(f22,f,0,0), FDDV(f21,f,0,0), FDDV(f23,f,0,0),
         FDDV(f33,f,0,0), FDDV(f31,f,0,0), FDDV(f32,f,0,0),
         FDDV(f11,f,1,0), FDDV(f12,f,1,0), FDDV(f13,f,1,0),
         FDDV(f22,f,1,0), FDDV(f21,f,1,0), FDDV(f23,f,1,0),
         FDDV(f33,f,1,0), FDDV(f31,f,1,0), FDDV(f32,f,1,0),
        &FDDV(pel->f[0],g,0,0), &FDDV(pel->f[1],g,0,0), &FDDV(pel->f[2],g,0,0),
        &FDDV(pel->f[0],g,1,0), &FDDV(pel->f[1],g,1,0), &FDDV(pel->f[2],g,1,0),
         n1,n2,n3,n12,n13,n23);
      add3_restrict_p1mod(
         FDDV(f11,f,0,1), FDDV(f12,f,0,1), FDDV(f13,f,0,1),
         FDDV(f22,f,0,1), FDDV(f21,f,0,1), FDDV(f23,f,0,1),
         FDDV(f33,f,0,1), FDDV(f31,f,0,1), FDDV(f32,f,0,1),
         FDDV(f11,f,1,1), FDDV(f12,f,1,1), FDDV(f13,f,1,1),
         FDDV(f22,f,1,1), FDDV(f21,f,1,1), FDDV(f23,f,1,1),
         FDDV(f33,f,1,1), FDDV(f31,f,1,1), FDDV(f32,f,1,1),
        &FDDV(pel->f[0],g,0,1), &FDDV(pel->f[1],g,0,1), &FDDV(pel->f[2],g,0,1),
        &FDDV(pel->f[0],g,1,1), &FDDV(pel->f[1],g,1,1), &FDDV(pel->f[2],g,1,1),
         n1,n2,n3,n12,n13,n23);
   }
}

void vect_p1mod_rhs_restrict_old(tGrid,f,g,q,t,type)
GRID *tGrid;          /*  coarse grid  */
INT f, g, q, t, type; /*  f ... fine grid; g ... coarse grid; q ... aux. var. */
{
   FACE *pface;

   set_value(tGrid,0.,q,0,type);
   for (pface=FIRSTFACE(tGrid); pface != NULL; pface=SUCC(pface)){
      FDDV(pface,q,0,0) = 1.;
      set_value(tGrid->finer,0.,q,0,type);
      vect_p1mod_prolongation(tGrid,q,q,t,type);
      FDDV(pface,g,0,0) = dot(tGrid->finer,f,q,t,type);

      FDDV(pface,q,0,0) = 0.;
      FDDV(pface,q,1,0) = 1.;
      set_value(tGrid->finer,0.,q,0,type);
      vect_p1mod_prolongation(tGrid,q,q,t,type);
      FDDV(pface,g,1,0) = dot(tGrid->finer,f,q,t,type);

      FDDV(pface,q,1,0) = 0.;
      FDDV(pface,q,0,1) = 1.;
      set_value(tGrid->finer,0.,q,0,type);
      vect_p1mod_prolongation(tGrid,q,q,t,type);
      FDDV(pface,g,0,1) = dot(tGrid->finer,f,q,t,type);

      FDDV(pface,q,0,1) = 0.;
      FDDV(pface,q,1,1) = 1.;
      set_value(tGrid->finer,0.,q,0,type);
      vect_p1mod_prolongation(tGrid,q,q,t,type);
      FDDV(pface,g,1,1) = dot(tGrid->finer,f,q,t,type);
      FDDV(pface,q,1,1) = 0.;
   }
}

#else  /*  if !(F_DATA & DVECTOR_FACE_DATA)  */

void vect_p1mod_prolongation(tGrid,u,v,t,type)  
GRID *tGrid; INT u, v, t;  
{  eprintf("Error: vect_p1mod_prolongation not available.\n");  }

void vect_p1mod_rhs_restrict(tGrid,f,g,t,type)
GRID *tGrid; INT f, g, t, type;
{  eprintf("Error: vect_p1mod_rhs_restrict not available.\n");  }
 
#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void p0_prolongation(tGrid,u,v,t)
GRID *tGrid;  /*  coarse grid  */
INT u, v, t;  /*  u ... coarse grid; v ... fine grid  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      ED(pel->sons[0],v) = ED(pel->sons[1],v) = 
      ED(pel->sons[2],v) = ED(pel->sons[3],v) = ED(pel,u);
}

void p0_rhs_restrict(tGrid,f,g,t)
GRID *tGrid;  /*  coarse grid  */
INT f, g, t;  /*  f ... fine grid; g ... coarse grid  */
{
   ELEMENT *pel, *first_el;

   if (t & WITHOUT_FIRST){
      first_el = FIRSTELEMENT(tGrid->finer);
      ED(first_el,f) = 0.;
      for (pel = first_el->succ; pel!=NULL; pel=pel->succ)
         ED(first_el,f) -= ED(pel,f);
   }
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      ED(pel,g) = ED(pel->sons[0],f) + ED(pel->sons[1],f) + 
                  ED(pel->sons[2],f) + ED(pel->sons[3],f);
}

#else  /*  if !(E_DATA & SCALAR_ELEMENT_DATA)  */

void p0_prolongation(tGrid,u,v,t)
GRID *tGrid; INT u, v, t;  
{  eprintf("Error: p0_prolongation not available.\n");  }

void p0_rhs_restrict(tGrid,f,g,t)
GRID *tGrid; INT f, g, t; 
{  eprintf("Error: p0_rhs_restrict not available.\n");  }

#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES

void p1_disc_prolong_on_son_el(pel,son3,soni,ai,aj,ak,k,v)
ELEMENT *pel, *son3, *soni;
FLOAT ai, aj, ak;
INT k, v;
{
   EDSN(soni,v,0) = (aj + ak)*0.5;
   if (IS_SON_FACE(soni->f[1],pel->f[k])){
      EDSN(soni,v,1) = ak + (aj - ai)*0.5;
      EDSN(soni,v,2) = aj + (ak - ai)*0.5;
   }
   else{
      EDSN(soni,v,1) = aj + (ak - ai)*0.5;
      EDSN(soni,v,2) = ak + (aj - ai)*0.5;
   }
   if (soni->f[0] == son3->f[0])
      EDSN(son3,v,0) = EDSN(soni,v,0);
   else if (soni->f[0] == son3->f[1])
      EDSN(son3,v,1) = EDSN(soni,v,0);
   else
      EDSN(son3,v,2) = EDSN(soni,v,0);
}
      
void p1_disc_prolongation(tGrid,u,v,t)
GRID *tGrid;  /*  coarse grid  */
INT u, v, t;  /*  u ... coarse grid; v ... fine grid  */
{
   ELEMENT *pel;
   FLOAT a0, a1, a2;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      a0 = EDSN(pel,u,0);
      a1 = EDSN(pel,u,1);
      a2 = EDSN(pel,u,2);
      p1_disc_prolong_on_son_el(pel,pel->sons[3],pel->sons[0],a0,a1,a2,2,v);
      p1_disc_prolong_on_son_el(pel,pel->sons[3],pel->sons[1],a1,a2,a0,0,v);
      p1_disc_prolong_on_son_el(pel,pel->sons[3],pel->sons[2],a2,a0,a1,1,v);
   }
}

void add_restrict_p1disc(pel,s,s0,k,e0,e1,e2,e3,f,g)
ELEMENT *pel, *e0, *e1, *e2, *e3;
FLOAT s, s0;
INT k, f, g;
{
   FLOAT c; 

   if (e0->f[0] == e3->f[0])
      c = EDSN(e0,f,0) - EDSN(e3,f,0);
   else if (e0->f[0] == e3->f[1])
      c = EDSN(e0,f,0) - EDSN(e3,f,1);
   else
      c = EDSN(e0,f,0) - EDSN(e3,f,2);
   if (IS_SON_FACE(e1->f[1],pel->f[k]))
      c += EDSN(e1,f,1);
   else
      c += EDSN(e1,f,2);
   if (IS_SON_FACE(e2->f[1],pel->f[k]))
      c += EDSN(e2,f,1);
   else
      c += EDSN(e2,f,2);
   EDSN(pel,g,k) = 0.5*(c + s) - s0;
}

void scalar_p1disc_rhs_restrict(tGrid,f,g,t)
GRID *tGrid;      /*  coarse grid  */
INT f, g, t;      /*  f ... fine grid; g ... coarse grid  */
{
   ELEMENT *pel, *e0, *e1, *e2, *e3;
   FLOAT s0, s1, s2, s3, s;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      e0 = pel->sons[0];
      e1 = pel->sons[1];
      e2 = pel->sons[2];
      e3 = pel->sons[3];
      s0 = EDSN(e0,f,0) + EDSN(e0,f,1) + EDSN(e0,f,2);
      s1 = EDSN(e1,f,0) + EDSN(e1,f,1) + EDSN(e1,f,2);
      s2 = EDSN(e2,f,0) + EDSN(e2,f,1) + EDSN(e2,f,2);
      s3 = EDSN(e3,f,0) + EDSN(e3,f,1) + EDSN(e3,f,2);
      s = s0 + s1 + s2 + s3;
      add_restrict_p1disc(pel,s,s0,0,e0,e1,e2,e3,f,g);
      add_restrict_p1disc(pel,s,s1,1,e1,e2,e0,e3,f,g);
      add_restrict_p1disc(pel,s,s2,2,e2,e0,e1,e3,f,g);
   }
}

#else  /* if  !(E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)  */

void p1_disc_prolongation(tGrid,u,v,t)
GRID *tGrid; INT u, v, t;  
{  eprintf("Error: p1_disc_prolongation not available.\n");  }

void scalar_p1disc_rhs_restrict(tGrid,f,g,t)
GRID *tGrid; INT f, g, t;
{  eprintf("Error: scalar_p1disc_rhs_restrict not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void vect_p1nc_to_p1mod_prolongation(tGrid,u,v)
GRID *tGrid; 
INT u, v;  /*  u ... p1nc; v ... p1mod  */
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3, *pface;
   ELEMENT *pel;

   for (pface=FIRSTFACE(tGrid); pface != NULL; pface=SUCC(pface))
      SET1(FDDVP(pface,v,0),FDVP(pface,u))
   for (pface=FIRSTF(tGrid); pface != FIRSTFACE(tGrid); pface=SUCC(pface))
      FDV(pface,u,0) = FDV(pface,u,1) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      FDDV(f1,v,1,0) += 5.*(FDV(f2,u,0)-FDV(f3,u,0))*NINDI(n2,n3);
      FDDV(f2,v,1,0) += 5.*(FDV(f3,u,0)-FDV(f1,u,0))*NINDI(n3,n1);
      FDDV(f3,v,1,0) += 5.*(FDV(f1,u,0)-FDV(f2,u,0))*NINDI(n1,n2);
      FDDV(f1,v,1,1) += 5.*(FDV(f2,u,1)-FDV(f3,u,1))*NINDI(n2,n3);
      FDDV(f2,v,1,1) += 5.*(FDV(f3,u,1)-FDV(f1,u,1))*NINDI(n3,n1);
      FDDV(f3,v,1,1) += 5.*(FDV(f1,u,1)-FDV(f2,u,1))*NINDI(n1,n2);
   }
}

void vect_p1mod_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid;
INT f, g; /*  f ... p1mod; g ... p1nc  */
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3, *pface;
   ELEMENT *pel;

   dvector_to_vector_f(tGrid,f,g,0);
   for (pface=FIRSTF(tGrid); pface != FIRSTFACE(tGrid); pface=SUCC(pface))
      FDDV(pface,f,1,0) = FDDV(pface,f,1,1) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      FDV(f1,g,0) +=
                 5.*(FDDV(f2,f,1,0)*NINDI(n1,n3) + FDDV(f3,f,1,0)*NINDI(n1,n2));      FDV(f2,g,0) +=
                 5.*(FDDV(f3,f,1,0)*NINDI(n2,n1) + FDDV(f1,f,1,0)*NINDI(n2,n3));      FDV(f3,g,0) +=
                 5.*(FDDV(f1,f,1,0)*NINDI(n3,n2) + FDDV(f2,f,1,0)*NINDI(n3,n1));      FDV(f1,g,1) +=
                 5.*(FDDV(f2,f,1,1)*NINDI(n1,n3) + FDDV(f3,f,1,1)*NINDI(n1,n2));      FDV(f2,g,1) +=
                 5.*(FDDV(f3,f,1,1)*NINDI(n2,n1) + FDDV(f1,f,1,1)*NINDI(n2,n3));      FDV(f3,g,1) +=
                 5.*(FDDV(f1,f,1,1)*NINDI(n3,n2) + FDDV(f2,f,1,1)*NINDI(n3,n1));   }
}

#else  /*  !((F_DATA & VECTOR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void vect_p1nc_to_p1mod_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vect_p1nc_to_p1mod_prolongation not available.\n");  }

void vect_p1mod_to_p1nc_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: vect_p1mod_to_p1nc_rhs_restrict not available.\n");  }

#endif

#if (N_DATA & NUMBER_OF_N_NEL) && (N_DATA & SCALAR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void scalar_p0_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid;
INT u, v;  /*  u ... p0; v ... p1c  */
{
   ELEMENT *pel;
   NODE *pnode;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NDS(pel->n[0],v) += ED(pel,u);
      NDS(pel->n[1],v) += ED(pel,u);
      NDS(pel->n[2],v) += ED(pel,u);
   }
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDS(pnode,v) /= pnode->nel;
}

void scalar_p1c_to_p0_rhs_restrict(tGrid,f,g)
GRID *tGrid;
INT f, g;         /*  f ... p1c; g ... p0  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      ED(pel,g) = NDS(pel->n[0],f)/pel->n[0]->nel +
                  NDS(pel->n[1],f)/pel->n[1]->nel +
                  NDS(pel->n[2],f)/pel->n[2]->nel;
}

#else

void scalar_p0_to_p1c_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_p0_to_p1c_prolongation not available.\n");  }

void scalar_p1c_to_p0_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: scalar_p1c_to_p0_rhs_restrict not available.\n");  }

#endif

#if (E_DATA & SCALAR_ELEMENT_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)

void scalar_p0_to_p1disc_prolongation(tGrid,u,v)
GRID *tGrid;
INT u, v;  /*  u ... p0; v ... p1disc  */
{  
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      EDSN(pel,v,0) = EDSN(pel,v,1) = EDSN(pel,v,2) = ED(pel,u);
}

void scalar_p1disc_to_p0_rhs_restrict(tGrid,f,g)
GRID *tGrid; 
INT f, g;         /*  f ... p1disc; g ... p0  */
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      ED(pel,g) = EDSN(pel,f,0) + EDSN(pel,f,1) + EDSN(pel,f,2);
}

#else  /*  !(E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)  */

void scalar_p0_to_p1disc_prolongation(tGrid,u,v)
GRID *tGrid; INT u, v;  
{  eprintf("Error: scalar_p0_to_p1disc_prolongation not available.\n");  }

void scalar_p1disc_to_p0_rhs_restrict(tGrid,f,g)
GRID *tGrid; INT f, g;
{  eprintf("Error: scalar_p1disc_to_p0_rhs_restrict not available.\n");  }

#endif

#define MAX_N_OF_NB_EL  20

#if (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_ELEMENTS) && (N_DATA & SCALAR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void ggem();

void scalar_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u,v)
GRID *tGrid;
INT u, v;     /*  u ... p0, v ... p1  */
{
   NODE *pn;
   LINK *pli;
   NELINK *pnel;
   ELEMENT *pel;
   DOUBLE a[MAX_N_OF_NB_EL][N_GGEM], b[MAX_N_OF_NB_EL], x[MAX_N_OF_NB_EL], 
          s, vol;
   INT ii[NVERT], i, j, n, m=MIN(MAX_N_OF_NB_EL,N_GGEM), ind[MAX_N_OF_NB_EL];

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ){
      ind[0] = INDEX(pn);
      n = 1;
      for (pli = TSTART(pn); pli && n < m; pli = pli->next)
         ind[n++] = INDEX(pli->nbnode);
      if (pli)
         eprintf("Error: two many neighbours in scalar_patchwise_p1_L2_prolong_of_p0_to_p1.\n");
      for (i = 0; i < n; i++){
         b[i] = 0.;
         for (j = 0; j < n; j++)
            a[i][j] = 0.;
      }
      for (pnel = NESTART(pn); pnel; pnel = pnel->next){
         pel = pnel->nbel;
         vol = VOLUME(pel);
         s = 4.*ED(pel,u)*vol;
         for (i = 0; i < NVERT; i++){
            for (j = 0; INDEX(pel->n[i]) != ind[j]; j++);
            ii[i] = j;
         }
         for (i = 0; i < NVERT; i++){
            b[ii[i]] += s;
            for (j = 0; j < NVERT; j++)
               a[ii[i]][ii[j]] += vol;
         }
      }
      for (i = 0; i < n; i++)
         a[i][i] += a[i][i];
      ggem(a,b,x,n);
      NDS(pn,v) = x[0];
   }
}

#else

void scalar_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: scalar_patchwise_p1_L2_prolong_of_p0_to_p1c not available.\n");  }

#endif

#if (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_ELEMENTS) && (N_DATA & VECTOR_NODE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void sggem();

void vector_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u,v)
GRID *tGrid;
INT u, v;     /*  u ... p0, v ... p1  */
{
   NODE *pn;
   LINK *pli;
   NELINK *pnel;
   ELEMENT *pel;
   DOUBLE a[MAX_N_OF_NB_EL][N_SGGEM], b[2][N_SGGEM], x[2][N_SGGEM], s[2], 
          z, vol;
   INT ii[NVERT], i, j, n, m=MIN(MAX_N_OF_NB_EL,N_SGGEM), ind[MAX_N_OF_NB_EL];

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ){
      ind[0] = INDEX(pn);
      n = 1;
      for (pli = TSTART(pn); pli && n < m; pli = pli->next)
         ind[n++] = INDEX(pli->nbnode);
      if (pli)
         eprintf("Error: two many neighbours in vector_patchwise_p1_L2_prolong_of_p0_to_p1.\n");
      for (i = 0; i < n; i++){
         b[0][i] = b[1][i] = 0.;
         for (j = 0; j < n; j++)
            a[i][j] = 0.;
      }
      for (pnel = NESTART(pn); pnel; pnel = pnel->next){
         pel = pnel->nbel;
         vol = VOLUME(pel);
         z = 4.*vol;
         SET2(s,EDVP(pel,u),z)
         for (i = 0; i < NVERT; i++){
            for (j = 0; INDEX(pel->n[i]) != ind[j]; j++);
            ii[i] = j;
         }
         for (i = 0; i < NVERT; i++){
            b[0][ii[i]] += s[0];
            b[1][ii[i]] += s[1];
            for (j = 0; j < NVERT; j++)
               a[ii[i]][ii[j]] += vol;
         }
      }
      for (i = 0; i < n; i++)
         a[i][i] += a[i][i];
      sggem(a,b,x,n,2);
      ND(pn,v,0) = x[0][0];
      ND(pn,v,1) = x[1][0];
   }
}

#else

void vector_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: vector_patchwise_p1_L2_prolong_of_p0_to_p1c not available.\n");  }

#endif

#if F_DATA & VECTOR_FACE_DATA

void prolong_low_to_high();

void vect_high_to_p1nc_rhs_restrict_test(tGrid,f,g,q,r)
GRID *tGrid;
INT f, g, q, r;
{
   FACE *pface;

   set_value(tGrid,0.,q,0,U_TYPE);
   for (pface=FIRSTFACE(tGrid); pface != NULL; pface=SUCC(pface)){
      FDV(pface,q,0) = 1.;
      prolong_low_to_high(tGrid,q,r,
                                T_FOR_U,U_TYPE_LOW,P1_NC,U_TYPE,U_SPACE,VECTOR);
      FDV(pface,g,0) = dot(tGrid,f,r,T_FOR_U,U_TYPE);
      FDV(pface,q,0) = 0.;

      FDV(pface,q,1) = 1.;
      prolong_low_to_high(tGrid,q,r,
                                T_FOR_U,U_TYPE_LOW,P1_NC,U_TYPE,U_SPACE,VECTOR);
      FDV(pface,g,1) = dot(tGrid,f,r,T_FOR_U,U_TYPE);
      FDV(pface,q,1) = 0.;
   }
}

#endif

void prolongation(tGrid,u,v,t,type,space,structure)  /* tGrid is coarse grid */
GRID *tGrid;                          /*  coarse grid  */
INT u, v, t, type, space, structure;  /*  u ... coarse grid; v ... fine grid  */
{
   if (structure == P_SCALAR){
      if (t & WITHOUT_FIRST)
         set_first(tGrid,0.,u,t,type);
   }
   else
      set_value(tGrid,0.,u,STOP_IS_FIRST_INNER,type);
   set_value(tGrid->finer,0.,v,0,type);
   if (t & ADD_ZN_CONDITION)
      set_value(tGrid->finer,0.,v,ADD_ZN_CONDITION,type);
   switch(space){
   case P1_NC:  if (structure == SCALAR)
                   scalar_p1nc_prolongation(tGrid,u,v,t,type);
                else if (structure == VECTOR)
                   vect_p1nc_prolongation(tGrid,u,v,t,type);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P1_MOD: if (structure == SCALAR)
                   scalar_p1mod_prolongation(tGrid,u,v,t,type);
                else if (structure == VECTOR)
                   vect_p1mod_prolongation(tGrid,u,v,t,type);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P0:     if (structure == P_SCALAR)
                   p0_prolongation(tGrid,u,v,t);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P1_DISC: if (structure == P_SCALAR)
                   p1_disc_prolongation(tGrid,u,v,t);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case IP1_DISC: if (structure == P_SCALAR)
                   p1_disc_prolongation(tGrid,u,v,t);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P1C:    if (structure == VECTOR)
                   vect_p1c_prolongation(tGrid,u,v,0); 
                else if (structure == SCALAR || structure == P_SCALAR){
                   scalar_p1c_prolongation(tGrid,u,v);
                   if (t & ADD_ZN_CONDITION)
                      scalar_p2c_prolongation_BD(tGrid,u,v);
                }
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case IP1C:   if (structure == VECTOR)
                   vect_p1c_prolongation(tGrid,u,v,0); 
                else if (structure == SCALAR || structure == P_SCALAR){
                   scalar_p1c_prolongation(tGrid,u,v);
                   if (t & ADD_ZN_CONDITION)
                      scalar_p2c_prolongation_BD(tGrid,u,v);
                }
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P1C_ELBUB: if (structure == VECTOR)
                   vect_p1c_prolongation(tGrid,u,v,1); 
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P2C:    if (structure == VECTOR)
                   vect_p2c_prolongation(tGrid,u,v,t,type,0);
                else if (structure == SCALAR || structure == P_SCALAR)
                   scalar_p2c_prolongation(tGrid,u,v);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case IP2C:    if (structure == VECTOR)
                   vect_p2c_prolongation(tGrid,u,v,t,type,0);
                else if (structure == SCALAR || structure == P_SCALAR)
                   scalar_p2c_prolongation(tGrid,u,v);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case P2C_ELBUB: if (structure == VECTOR)
                   vect_p2c_prolongation(tGrid,u,v,t,type,1);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   case IP2C_ELBUB: if (structure == VECTOR)
                   vect_p2c_prolongation(tGrid,u,v,t,type,1);
                else
                   eprintf("Error: prolongation not available.\n");
        break;
   default:
        eprintf("Error: prolongation not available.\n");
        break;
   }
   if (structure == P_SCALAR){
      if (t & WITHOUT_FIRST)
         subtr_first(tGrid->finer,v,t,type);
      else if (t & ZERO_MEAN)
         set_zero_mean(tGrid->finer,v,t,type);
   }
}

/* RESTRICTION OF A FUNCTION, NOT OF A RIGHT-HAND SIDE !!! */
void restriction(tGrid,v,u,q,t,type,space,structure)  /* tGrid is coarse grid */
GRID *tGrid;                          /*  coarse grid  */
INT v, u, q, t, type, space, structure;/* v ... fine grid; u ... coarse grid  */
{
   set_value(tGrid,0.,u,0,type);
   switch(space){
   case P1_NC:  if (structure == SCALAR)
                   scalar_p1nc_restriction(tGrid,v,u,q,type);
                else if (structure == VECTOR)
                   vect_p1nc_restriction(tGrid,v,u,q,type);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   case P2C:    if (structure == VECTOR)
                   vect_p2c_restriction(tGrid,v,u,0);
                else if (structure == SCALAR)
                   scalar_p2c_restriction(tGrid,v,u);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   case IP2C:   if (structure == VECTOR)
                   vect_p2c_restriction(tGrid,v,u,0);
                else if (structure == SCALAR)
                   scalar_p2c_restriction(tGrid,v,u);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   case P2C_ELBUB: if (structure == VECTOR)
                   vect_p2c_restriction(tGrid,v,u,1);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   case IP2C_ELBUB: if (structure == VECTOR)
                   vect_p2c_restriction(tGrid,v,u,1);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   case GQ2C:   if (structure == SCALAR)
                   scalar_gq2c_restriction(tGrid,v,u);
                else
                   eprintf("Error: restriction not available.\n");
        break;
   default:
        eprintf("Error: restriction not available.\n");
        break;
   }
}

void rhs_restrict(tGrid,f,g,t,type,space,structure)
GRID *tGrid;                     /*  coarse grid  */
INT f, g, t, type, space, structure; /*  f ... fine grid; g ... coarse grid  */
{
   if (structure != P_SCALAR)
      set_value(tGrid->finer,0.,f,STOP_IS_FIRST_INNER,type);
   set_value(tGrid,0.,g,0,type);
   if (t & ADD_ZN_CONDITION)
      set_value(tGrid,0.,g,ADD_ZN_CONDITION,type);
   switch(space){
   case P1_NC:  if (structure == SCALAR)
                   scalar_p1nc_rhs_restrict(tGrid,f,g,t,type);
                else if (structure == VECTOR)
                   vect_p1nc_rhs_restrict(tGrid,f,g,t,type);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P1_MOD: if (structure == SCALAR)
                   scalar_p1mod_rhs_restrict(tGrid,f,g,t,type);
                else if (structure == VECTOR)
                   vect_p1mod_rhs_restrict(tGrid,f,g,t,type);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P0:     if (structure == P_SCALAR)
                   p0_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P1_DISC: if (structure == P_SCALAR)
                   scalar_p1disc_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case IP1_DISC: if (structure == P_SCALAR)
                   scalar_p1disc_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P1C:    if (structure == VECTOR)
                   vect_p1c_rhs_restrict(tGrid,f,g);
                else if (structure == SCALAR)
                   scalar_p1c_rhs_restrict(tGrid,f,g,
                                              ONLY_INNER | NSTART_FROM_INNER);
                else if (structure == P_SCALAR)
                   p_scalar_p1c_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case IP1C:    if (structure == VECTOR)
                   vect_p1c_rhs_restrict(tGrid,f,g);
                else if (structure == SCALAR)
                   scalar_p1c_rhs_restrict(tGrid,f,g,
                                              ONLY_INNER | NSTART_FROM_INNER);
                else if (structure == P_SCALAR)
                   p_scalar_p1c_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P1C_ELBUB: if (structure == VECTOR){
                   vect_p1c_rhs_restrict(tGrid,f,g);
                   vect_cubic_bubble1_rhs_restrict(tGrid,f,g);
                }
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P2C:    if (structure == VECTOR)
                   vect_p2c_rhs_restrict(tGrid,f,g,t,type);
                else if (structure == SCALAR)
                   scalar_p2c_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case IP2C:   if (structure == VECTOR)
                   vect_p2c_rhs_restrict(tGrid,f,g,t,type);
                else if (structure == SCALAR)
                   scalar_p2c_rhs_restrict(tGrid,f,g,t);
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case P2C_ELBUB: if (structure == VECTOR){
                   vect_p2c_rhs_restrict(tGrid,f,g,t,type);
                   vect_cubic_bubble_rhs_restrict(tGrid,f,g);
                }
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   case IP2C_ELBUB: if (structure == VECTOR){
                   vect_p2c_rhs_restrict(tGrid,f,g,t,type);
                   vect_cubic_bubble_rhs_restrict(tGrid,f,g);
                }
                else
                   eprintf("Error: rhs_restrict not available.\n");
        break;
   default:
        eprintf("Error: rhs_restrict not available.\n");
        break;
   }
}

void prolong_low_to_high(tGrid,u_low,u_high,
                        low_t,low_type,low_space,high_type,high_space,structure)
GRID *tGrid;
INT u_low, u_high, low_t, low_type, low_space, high_type, high_space, structure;
{
   if (low_space == high_space)
      copy(tGrid,u_low,u_high,low_t,low_type);
   else{
      if (structure != P_SCALAR)
         set_value(tGrid,0.,u_low,STOP_IS_FIRST_INNER,low_type);
      set_value(tGrid,0.,u_high,0,high_type);
      switch(low_space){
         case P1_NC: if (high_space == P1_MOD && structure == VECTOR)
                        vect_p1nc_to_p1mod_prolongation(tGrid,u_low,u_high);
                     else if (high_space == P1C && structure == SCALAR)
                        scalar_p1nc_to_p1c_prolongation(tGrid,u_low,u_high);
                     else if (high_space == P1C && structure == VECTOR)
                        vect_p1nc_to_p1c_prolongation(tGrid,u_low,u_high);
                     else if ((high_space == P2C || high_space == IP2C) && 
                               structure == VECTOR)
                        vect_p1nc_to_p2c_prolongation(tGrid,u_low,u_high);
                     else if ((high_space == P2C_ELBUB || 
                             high_space == IP2C_ELBUB) && structure == VECTOR){
                        vect_p1nc_to_p2c_prolongation(tGrid,u_low,u_high);
/*
                        set_bubbles_for_vect_p1nc_to_p2cb_prolongation(tGrid,
                                                                 u_low,u_high);
*/
                     }
                     else
                        eprintf("Error: prolong_low_to_high not available.\n");
              break;
         case P0:    if ((high_space == P1C || high_space == IP1C) && 
                        (structure == SCALAR || structure == P_SCALAR))
                        scalar_p0_to_p1c_prolongation(tGrid,u_low,u_high);
          /* scalar_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u_low,u_high); */
                     else if ((high_space == P1_DISC || high_space == IP1_DISC)
                        && (structure == SCALAR || structure == P_SCALAR))
                        scalar_p0_to_p1disc_prolongation(tGrid,u_low,u_high);
                     else
                        eprintf("Error: prolong_low_to_high not available.\n");
          /* vector_patchwise_p1_L2_prolong_of_p0_to_p1c(tGrid,u_low,u_high); */
              break;
         case P1C:   if (high_space == P2C || high_space == IP2C ||
                         high_space == P2C_ELBUB || high_space == IP2C_ELBUB)
                        copy(tGrid,u_low,u_high,low_t,low_type);
                     else
                        eprintf("Error: prolong_low_to_high not available.\n");
              break;
         case P2C:   if (high_space == P2C_ELBUB || high_space == IP2C_ELBUB)
                        copy(tGrid,u_low,u_high,low_t,low_type);
                     else
                        eprintf("Error: prolong_low_to_high not available.\n");
              break;
         default:
              eprintf("Error: prolong_low_to_high not available.\n");
              break;
      }
   }
}

/* RESTRICTION OF A FUNCTION, NOT OF A RIGHT-HAND SIDE !!! */
void restriction_high_to_low(tGrid,u_high,u_low,
                             high_type,high_space,low_type,low_space,structure)
GRID *tGrid;
INT u_high, u_low, high_type, high_space, low_type, low_space, structure;
{
   if (low_space == high_space)
      copy(tGrid,u_high,u_low,0,low_type);
   else{
      set_value(tGrid,0.,u_low,0,low_type);
      switch(high_space){
         case P1C:    if (low_space == P1_NC && structure == SCALAR)
                         scalar_p1c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p1c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case IP1C:   if (low_space == P1_NC && structure == SCALAR)
                         scalar_p1c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p1c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case P2C:    if (low_space == P1_NC && structure == SCALAR)
                         scalar_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1C)
                         copy(tGrid,u_high,u_low,0,low_type);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case IP2C:   if (low_space == P1_NC && structure == SCALAR)
                         scalar_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1C)
                         copy(tGrid,u_high,u_low,0,low_type);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case P2C_ELBUB: if (low_space == P1_NC && structure == SCALAR)
                         scalar_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1C)
                         copy(tGrid,u_high,u_low,0,low_type);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case IP2C_ELBUB: if (low_space == P1_NC && structure == SCALAR)
                         scalar_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_restriction(tGrid,u_high,u_low);
                      else if (low_space == P1C)
                         copy(tGrid,u_high,u_low,0,low_type);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         case GQ2X4C: if (low_space == GQ2C && structure == SCALAR)
                         scalar_gq2x4c_to_gq2c_restriction(tGrid,u_high,u_low);
                      else
                         eprintf("Error: restriction_high_to_low not available.\n");
              break;
         default:
              eprintf("Error: restriction_high_to_low not available.\n");
              break;
      }
   }
}

void rhs_restrict_high_to_low(tGrid,f_high,f_low,
                        high_type,high_space,low_t,low_type,low_space,structure)
GRID *tGrid;
INT f_high, f_low, high_type, high_space, low_t, low_type, low_space, structure;
{
   if (low_space == high_space)
      copy(tGrid,f_high,f_low,low_t,low_type);
   else{
      if (structure != P_SCALAR)
         set_value(tGrid,0.,f_high,STOP_IS_FIRST_INNER,high_type);
      switch(high_space){
         case P1_MOD: if (low_space == P1_NC && structure == VECTOR)
                         vect_p1mod_to_p1nc_rhs_restrict(tGrid,f_high,f_low);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case P1_DISC: if (structure == P_SCALAR)
                         scalar_p1disc_to_p0_rhs_restrict(tGrid,f_high,f_low);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case IP1_DISC: if (structure == P_SCALAR)
                         scalar_p1disc_to_p0_rhs_restrict(tGrid,f_high,f_low);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case P1C:    if (structure == P_SCALAR)
                         scalar_p1c_to_p0_rhs_restrict(tGrid,f_high,f_low);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case IP1C:   if (structure == P_SCALAR)
                         scalar_p1c_to_p0_rhs_restrict(tGrid,f_high,f_low);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case P2C:    if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_rhs_restrict(tGrid,f_high,f_low);
                      else if (low_space == P1C && structure != P_SCALAR)
                         copy(tGrid,f_high,f_low,low_t,low_type);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case IP2C:   if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_rhs_restrict(tGrid,f_high,f_low);
                      else if (low_space == P1C && structure != P_SCALAR)
                         copy(tGrid,f_high,f_low,low_t,low_type);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case P2C_ELBUB: if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_rhs_restrict(tGrid,f_high,f_low);
/*                       vect_p2cb_to_p1nc_rhs_restrict(tGrid,f_high,f_low); */ 
                      else if (low_space == P2C || low_space == P1C)
                         copy(tGrid,f_high,f_low,low_t,low_type);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         case IP2C_ELBUB: if (low_space == P1_NC && structure == VECTOR)
                         vect_p2c_to_p1nc_rhs_restrict(tGrid,f_high,f_low);
/*                       vect_p2cb_to_p1nc_rhs_restrict(tGrid,f_high,f_low); */
                      else if (low_space == P2C || low_space == P1C)
                         copy(tGrid,f_high,f_low,low_t,low_type);
                      else
                         eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
         default:
              eprintf("Error: rhs_restrict_high_to_low not available.\n");
              break;
      }
   }
}

#else  /*  DIM != 2  */

void prolongation(tGrid,u,v,t,type,space,structure)  /* tGrid is coarse grid */
GRID *tGrid; INT u, v, t, type, space, structure;
{  eprintf("Error: prolongation not available.\n");  }

/* RESTRICTION OF A FUNCTION, NOT OF A RIGHT-HAND SIDE !!! */
void restriction(tGrid,v,u,q,t,type,space,structure)  /* tGrid is coarse grid */
GRID *tGrid; INT v, u, q, t, type, space, structure;
{  eprintf("Error: restriction not available.\n");  }

void rhs_restrict(tGrid,f,g,t,type,space,structure)
GRID *tGrid; INT f, g, t, type, space, structure; 
{  eprintf("Error: rhs_restrict not available.\n");  }

void prolong_low_to_high(tGrid,u_low,u_high,
                        low_t,low_type,low_space,high_type,high_space,structure)
GRID *tGrid; INT u_low, u_high, low_t, low_type, low_space, high_type, high_space, structure;
{  eprintf("Error: prolong_low_to_high not available.\n");  }

void rhs_restrict_high_to_low(tGrid,f_high,f_low,
                        high_type,high_space,low_t,low_type,low_space,structure)
GRID *tGrid; INT f_high, f_low, high_type, high_space, low_t, low_type, low_space, structure;
{  eprintf("Error: rhs_restrict_high_to_low not available.\n");  }

#endif

/* ========================================================================== */

/* u and ue are assumed to contain 0. on triangulation corresponding to tGrid.*/
void nc_smoothing_coarse_grid_solver(tGrid,level,m,n1,n2,print,r,alpha,
                   ZA,ZB,h,f,fe,u,ue,b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
                   imin,imax,r1,r2,t_u,t_p,u_type,p_type,u_space,p_space,
                   A_struct,B_struct,C_struct,
                   smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
                   coarse_grid_solver)
GRID *tGrid;
FLOAT r, alpha, r1, r2;
INT level, m, n1, n2, print;
INT ZA, ZB, h, f, fe, u, ue, b, c, d, g, p, q, w, de, ge, pe, qe, we, xe, ye, 
    imin, imax, t_u, t_p, u_type, p_type, u_space, p_space, A_struct, B_struct,
    C_struct, precond_matrix, q0, q1, q2, precond_type;
INT (*smoother)();
void (*multB_X_BT)();
void (*coarse_grid_solver)();
{
   FLOAT def1, def2, eps1, eps2, rhs1, rhs2;
   INT i=0;

   inv(tGrid,f,d,t_u,u_type);
   inv(tGrid,fe,de,t_p,p_type);
   rhs1 = def1 = dot(tGrid,d,d,t_u,u_type);
   rhs2 = def2 = dot(tGrid,de,de,t_p,p_type);
   if (rhs1 < 1.e-15)
      rhs1 = 1.;
   if (rhs2 < 1.e-15)
      rhs2 = 1.;
   eps1 = rhs1*r1*r1;
   eps2 = rhs2*r2*r2;
   if (eps1 < EPST || eps2 < EPST){
      imin = 2;
      eps1 = MAX(eps1,EPST);
      eps2 = MAX(eps2,EPST);
   }

   while (i < imin || (i < imax && (def1 > eps1 || def2 > eps2))){

      smoother(tGrid,0,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
               d,p,q,de,pe,qe,xe,ye,
               t_u,t_p,u_type,p_type,A_struct,B_struct,
               q0,q1,q2,precond_type,multB_X_BT);
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,
                    t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
      def1 = dot(tGrid,d,d,t_u,u_type);
      def2 = dot(tGrid,de,de,t_p,p_type);
      i++;
   }
   if (print == YES)
      printf("     %i iterations.   u_def/rhs: %e  p_def/rhs: %e\n",i,
                                    sqrt(dot(tGrid,d,d,t_u,u_type)/rhs1),
                                    sqrt(dot(tGrid,de,de,t_p,p_type)/rhs2));
}

void Stokes_mgm();

/* u and ue are assumed to contain 0. on triangulation corresponding to tGrid.*/
void Stokes_mg_coarse_grid_solver(tGrid,level,m,n1,n2,print,r,alpha,
                   ZA,ZB,h,f,fe,u,ue,b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
                   imin,imax,r1,r2,t_u,t_p,u_type,p_type,u_space,p_space,
                   A_struct,B_struct,C_struct,
                   smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
                   coarse_grid_solver)
GRID *tGrid;
FLOAT r, alpha, r1, r2;
INT level, m, n1, n2, print;
INT ZA, ZB, h, f, fe, u, ue, b, c, d, g, p, q, w, de, ge, pe, qe, we, xe, ye, 
    imin, imax, t_u, t_p, u_type, p_type, u_space, p_space, A_struct, B_struct,
    C_struct, precond_matrix, q0, q1, q2, precond_type;
INT (*smoother)();
void (*multB_X_BT)();
void (*coarse_grid_solver)();
{
   FLOAT def1, def2, eps1, eps2, rhs1, rhs2;
   INT i=0;
  
   inv(tGrid,f,d,t_u,u_type);
   inv(tGrid,fe,de,t_p,p_type);
   rhs1 = def1 = dot(tGrid,d,d,t_u,u_type);
   rhs2 = def2 = dot(tGrid,de,de,t_p,p_type);
   if (rhs2 < 1.e-15)
      rhs2 = 1.;
   eps1 = rhs1*r1*r1;
   eps2 = rhs2*r2*r2;
   if (eps1 < EPST || eps2 < EPST){
      imin = 2;
      eps1 = MAX(eps1,EPST);
      eps2 = MAX(eps2,EPST);
   }
   
   while (i < imin || (i < imax && (def1 > eps1 || def2 > eps2))){

      Stokes_mgm(tGrid,level,m,6,2,0,r,alpha,1.,1.,ZA,ZB,h,f,fe,u,ue, 
             b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
             t_u,t_p,u_type,p_type,u_space,p_space,A_struct,B_struct,C_struct,
             smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
             coarse_grid_solver,alpha);
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,
                    t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
      def1 = dot(tGrid,d,d,t_u,u_type);
      def2 = dot(tGrid,de,de,t_p,p_type);
/*     printf("(def1/rhs1)^2 = %e, (def2/rhs2)^2 = %e\n",def1/rhs1,def2/rhs2);*/
      i++;
   }
   printf("     %i iterations.   u_def/rhs: %e  p_def/rhs: %e\n",i,
                                    sqrt(dot(tGrid,d,d,t_u,u_type)/rhs1),
                                    sqrt(dot(tGrid,de,de,t_p,p_type)/rhs2));
}

void Stokes_mgm(tGrid,level,m,n1,n2,print,r,alpha,u_damp,p_damp,
    ZA,ZB,h,f,fe,u,ue,b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
    t_u,t_p,u_type,p_type,u_space,p_space,A_struct,B_struct,C_struct,
    smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
    coarse_grid_solver,coarse_grid_alpha)
GRID *tGrid;
FLOAT r, alpha, coarse_grid_alpha, u_damp, p_damp;
INT level, m, n1, n2, print; 
INT ZA, ZB, h, f, fe, u, ue, b, c, d, g, p, q, w, de, ge, pe, qe, we, xe, ye,
    t_u, t_p, u_type, p_type, u_space, p_space, A_struct, B_struct, C_struct,
    precond_matrix, q0, q1, q2, precond_type;
INT (*smoother)();
void (*multB_X_BT)();
void (*coarse_grid_solver)();
{
   INT i;

   if (tGrid->level == level){
      if (print == YES) printf("Level %1i:   solving.\n",tGrid->level);
      coarse_grid_solver(tGrid,1,1,n1,n2,print,r,coarse_grid_alpha,
              ZA,ZB,h,f,fe,u,ue,b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
              0,1000,7.e-8,7.e-8,
              t_u,t_p,u_type,p_type,u_space,p_space,A_struct,B_struct,C_struct,
              smoother,h,A_struct,q1,q2,ILU_PR,multB_ILU_BT,
              nc_smoothing_coarse_grid_solver);
   }
   else{
      if (print == YES) printf("Level %1i:   %i smoothings.\n",tGrid->level,n1);
if (tGrid->level == TOPLEVEL) AUX1 = dot(tGrid,u,u,t_u,u_type);
      for (i = 0; i < n1; i++)
         smoother(tGrid,2,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
                  d,p,q,de,pe,qe,xe,ye,
                  t_u,t_p,u_type,p_type,A_struct,B_struct,
                  q0,q1,q2,precond_type,multB_X_BT);
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,
                              t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
if (tGrid->level == TOPLEVEL) AUX1 = sqrt(dot(tGrid,d,d,t_u,u_type)/AUX1);
      rhs_restrict(tGrid->coarser,d,f,t_u,u_type,u_space,VECTOR);
/*    set_value(tGrid->coarser,0.,fe,t_p,p_type);  */
      rhs_restrict(tGrid->coarser,de,fe,t_p,p_type,p_space,P_SCALAR);
      set_value(tGrid->coarser,0.,u,t_u,u_type);
      set_value(tGrid->coarser,0.,ue,t_p,p_type);
      if (tGrid->coarser->level == level)
         m = 1;
      for (i = 0; i < m; i++)
         Stokes_mgm(tGrid->coarser,level,m,n1,n2,print,r,alpha,u_damp,p_damp,
               ZA,ZB,h,f,fe,u,ue,
               b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
               t_u,t_p,u_type,p_type,u_space,p_space,A_struct,B_struct,C_struct,
               smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
               coarse_grid_solver,coarse_grid_alpha);
      prolongation(tGrid->coarser,u,w,t_u,u_type,u_space,VECTOR);
      prolongation(tGrid->coarser,ue,we,t_p,p_type,p_space,P_SCALAR);
/*
      subtr(tGrid,u,w,u,t_u,u_type); 
      subtr(tGrid,ue,we,ue,t_p,p_type); 
*/
      if (fabs(u_damp-1.) < -1.e-15 && fabs(p_damp-1.) < 1.e-15){
         mult_by_Stokes_matrix(tGrid,ZA,ZB,w,we,p,pe,q,
                                       t_u,t_p,u_type,p_type,A_struct,B_struct);
         u_damp = p_damp = 
                       (dot(tGrid,d,p,t_u,u_type)+dot(tGrid,de,pe,t_p,p_type))/
                       (dot(tGrid,p,p,t_u,u_type)+dot(tGrid,pe,pe,t_p,p_type));
         u_damp = p_damp = dot(tGrid,d,p,t_u,u_type)/dot(tGrid,p,p,t_u,u_type);
         printf("u_damp = %e, p_damp = %e\n",u_damp,p_damp);
         u_damp = p_damp = 1.;
      }
      mult_and_add(tGrid,-u_damp,w,u,u,t_u,u_type);
      mult_and_add(tGrid,-p_damp,we,ue,ue,t_p,p_type);
      if (print == YES) printf("Level %1i:   %i smoothings.\n",tGrid->level,n2);
      for (i = 0; i < n2; i++)
         smoother(tGrid,2,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
                  d,p,q,de,pe,qe,xe,ye,
                  t_u,t_p,u_type,p_type,A_struct,B_struct,
                  q0,q1,q2,precond_type,multB_X_BT);
   }
}

void Stokes_mgm_low_and_high(tGrid,level,m,n1,n2,m_low,n1_low,n2_low,
    print,r,alpha,alpha_low,
    ZA,ZB,ZA_low,ZB_low,h,f,fe,f_low,fe_low,u,ue,u_low,ue_low,
    b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
    t_u,t_p,u_type,p_type,u_space,p_space,A_struct,B_struct,C_struct,
    low_t_u,low_t_p,low_u_type,low_p_type,low_u_space,low_p_space,
    low_A_struct,low_B_struct,low_C_struct,
    smoother,precond_matrix,q0,q1,q2,precond_type,multB_X_BT,
    low_smoother,low_precond_matrix,low_q0,low_q1,low_q2,low_precond_type,
    low_multB_X_BT,coarse_grid_solver,coarse_grid_alpha)
GRID *tGrid;
FLOAT r, alpha, alpha_low, coarse_grid_alpha;
INT level, m, m_low, n1, n2, n1_low, n2_low, print; 
INT ZA, ZB, ZA_low, ZB_low, h, f, fe, f_low, fe_low, u, ue, u_low, ue_low, 
    b, c, d, g, p, q, w, de, ge, pe, qe, we, xe, ye,
    t_u, t_p, u_type, p_type, u_space, p_space, A_struct, B_struct, C_struct,
    low_t_u, low_t_p, low_u_type, low_p_type, low_u_space, low_p_space, 
    low_A_struct, low_B_struct, low_C_struct,
    precond_matrix, q0, q1, q2, precond_type,
    low_precond_matrix, low_q0, low_q1, low_q2, low_precond_type;
INT (*smoother)(), (*low_smoother)();
void (*multB_X_BT)(), (*low_multB_X_BT)();
void (*coarse_grid_solver)();
{
   INT i;
   
   if (print == YES) 
      printf("Level %1i:   %i smoothings for higher order discretization.\n",
                                                               tGrid->level,n1);
   for (i = 0; i < n1; i++)
      smoother(tGrid,2,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
               d,p,q,de,pe,qe,xe,ye,
               t_u,t_p,u_type,p_type,A_struct,B_struct,
               q0,q1,q2,precond_type,multB_X_BT);
   Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,
                            t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
   rhs_restrict_high_to_low(tGrid,d,f_low,
                        u_type,u_space,low_t_u,low_u_type,low_u_space,VECTOR);
   rhs_restrict_high_to_low(tGrid,de,fe_low,
                        p_type,p_space,low_t_p,low_p_type,low_p_space,P_SCALAR);
/* set_value(tGrid,0.,fe_low,low_t_p,low_p_type); */
   set_value(tGrid,0.,u_low,low_t_u,low_u_type);
   set_value(tGrid,0.,ue_low,low_t_p,low_p_type);
   if (tGrid->level == level)
      m = 1;
   for (i = 0; i < m; i++)
      Stokes_mgm(tGrid,level,m_low,n1_low,n2_low,print,r,alpha_low,1.,1.,
             ZA_low,ZB_low,h,f_low,fe_low,u_low,ue_low,
             b,c,d,g,p,q,w,de,ge,pe,qe,we,xe,ye,
             low_t_u,low_t_p,low_u_type,low_p_type,low_u_space,low_p_space,
             low_A_struct,low_B_struct,low_C_struct,
             low_smoother,low_precond_matrix,low_q0,low_q1,low_q2,
             low_precond_type,low_multB_X_BT,
             coarse_grid_solver,coarse_grid_alpha);
   prolong_low_to_high(tGrid,u_low,w,
                        low_t_u,low_u_type,low_u_space,u_type,u_space,VECTOR);
   prolong_low_to_high(tGrid,ue_low,we,
                        low_t_p,low_p_type,low_p_space,p_type,p_space,P_SCALAR);
   subtr(tGrid,u,w,u,t_u,u_type); 
   subtr(tGrid,ue,we,ue,t_p,p_type); 
   if (print == YES) 
      printf("Level %1i:   %i smoothings for higher order discretization.\n",
                                                               tGrid->level,n2);
   for (i = 0; i < n2; i++)
      smoother(tGrid,2,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
               d,p,q,de,pe,qe,xe,ye,
               t_u,t_p,u_type,p_type,A_struct,B_struct,
               q0,q1,q2,precond_type,multB_X_BT);
}

void a_smoothing_coarse_grid_solver(tGrid,Z,h,f,u,d,q,imin,imax,r,
                t,type,A_struct,smoother_type)
GRID *tGrid;
FLOAT r;
INT Z, h, f, u, d, q, imin, imax, t, type, A_struct, smoother_type;
{
   INT i;
   FLOAT eps, def;
   
   set_value(tGrid,0.,u,t,type);
   def = dot(tGrid,f,f,t,type);
   eps = def*r*r;
   if (eps < 1.e-18){
      imin = 2;
      eps = 1.e-18;
   }
   i = 0;
   while (i < imin || (i < imax && def > eps)){
      a_smoother(tGrid,Z,h,f,u,d,q,t,type,A_struct,smoother_type);
      if (smoother_type != ILU_IT)
         defect(tGrid,Z,f,u,d,d,t,type,A_struct);
      def = dot(tGrid,d,d,t,type);
      i++;
   }
}

void mgm(tGrid,level,m,n1,n2,print,damp,Z,h,f,u,d,q,
                t,type,space,A_struct,structure,smoother_type)
GRID *tGrid;
FLOAT damp;
INT level, m, n1, n2, print; 
INT Z, h, f, u, d, q, t, type, space, A_struct, structure, smoother_type;
{
   INT i;
   
   if (tGrid->level == level){
      if (print == YES) printf("Level %1i:   solving.\n",tGrid->level);
      a_smoothing_coarse_grid_solver(tGrid,Z,h,f,u,d,q,1,1000,1.e-8,
                                     t,type,A_struct,smoother_type);
   }
   else{
      if (print == YES) printf("Level %1i:   %i smoothings.\n",tGrid->level,n1);
      for (i = 0; i < n1; i++)
         a_smoother(tGrid,Z,h,f,u,d,q,t,type,A_struct,smoother_type);
      defect(tGrid,Z,f,u,d,d,t,type,A_struct);
      rhs_restrict(tGrid->coarser,d,f,t,type,space,structure);
      set_value(tGrid->coarser,0.,u,t,type);
      if (tGrid->coarser->level == level)
         m = 1;
      for (i = 0; i < m; i++)
         mgm(tGrid->coarser,level,m,n1,n2,print,damp,Z,h,f,u,d,q,
             t,type,space,A_struct,structure,smoother_type);
      prolongation(tGrid->coarser,u,d,t,type,space,structure);
/*
      subtr(tGrid,u,d,u,t,type); 
*/
      mult_and_add(tGrid,-damp,d,u,u,t,type);
      if (print == YES) printf("Level %1i:   %i smoothings.\n",tGrid->level,n2);
      for (i = 0; i < n2; i++)
         a_smoother(tGrid,Z,h,f,u,d,q,t,type,A_struct,smoother_type);
   }
}

#if N_DATA & SCALAR_NODE_DATA

void smgm(tGrid,level,m,n1,n2,Z,h,f,u,d,p,q,t)
GRID *tGrid;
INT level, m, n1, n2, Z, h, f, u, d, p, q, t;
{
   INT i;
   
   printf("Level %1i:\n",tGrid->level);
   if (tGrid->level == level)
/*       PCG(tGrid,f,u,d,p,q,0,10000,1.e50,EPS,0,mult_A,Z,
         d,t,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
         h,0,0,0,ILU_PR);*/
      GMRES(tGrid,U,F,10,0,20,D,10,0.01,1.,NONZERO_INIT,mult_A,Z,
                                      U,t,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   else{
      for (i = 0; i < n1; i++)
         sGauss_Seidel_step(tGrid,Z,f,u,t);
      smultA(tGrid,Z,u,d,t);
      ssubtr(tGrid,d,f,d,t);
      scalar_p1c_rhs_restrict(tGrid->coarser,d,f,t);
      sset_value(tGrid->coarser,0.,u,t);
      for (i = 0; i < m; i++)
         smgm(tGrid->coarser,level,m,n1,n2,Z,h,f,u,d,p,q,t);
      printf("Level %1i:\n",tGrid->level);
      scalar_p1c_prolongation(tGrid->coarser,u,p);
      ssubtr(tGrid,u,p,u,t);
      for (i = 0; i < n2; i++)
         sGauss_Seidel_step(tGrid,Z,f,u,t);
   }
}

#else

void smgm(tGrid,level,m,n1,n2,Z,h,f,u,d,p,q,t)
GRID *tGrid; INT level, m, n1, n2, Z, h, f, u, d, p, q, t;
{  eprintf("Error: smgm not available.\n");  }

#endif
