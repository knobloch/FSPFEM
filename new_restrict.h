/******************************************************************************/
/*                                                                            */
/*                        restriction and prolongation                        */
/*                                                                            */
/******************************************************************************/

#if DIM == 3

#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR)

/* pi, qi, ri are values of linear functions p, q, r in vertices of 
   a triangle T. The value (60/|T|)\int_T pqr\dx is computed.                 */
FLOAT integr4(p1,p2,p3,q1,q2,q3,r1,r2,r3) 
FLOAT p1, p2, p3, q1, q2, q3, r1, r2, r3;
{
   FLOAT u, v, w;
   
   u = (p1+p2+p3)*(q1+q2+q3)*(r1+r2+r3);
   v = p1*q1*r1 + p2*q2*r2 + p3*q3*r3;
   w = p1*(q2*r3+q3*r2) + p2*(q1*r3+q3*r1) + p3*(q1*r2+q2*r1);
   return(u+u + v+v+v+v - w);
}

#if !(DATA_STR & LG_DATA)

void norm_vect(pel,i,n)
ELEMENT *pel;
INT i;
FLOAT *n;
{
   FLOAT q = FMULT/COEFF_BF(pel,0,i)/20.0;
   
   n[0] = -COEFF_BN(pel,0,i,0)*q;
   n[1] = -COEFF_BN(pel,0,i,1)*q;
   n[2] = -COEFF_BN(pel,0,i,2)*q;
}

void bar_coord(pel,b)
ELEMENT *pel;
FLOAT b[4][4];
{
   FLOAT *x0, *x1, vol;
   
   x0 = pel->n[0]->myvertex->x;
   x1 = pel->n[1]->myvertex->x;
   vol = COEFF_BN(pel,0,0,0)*(x1[0]-x0[0]) + 
         COEFF_BN(pel,0,0,1)*(x1[1]-x0[1]) + 
         COEFF_BN(pel,0,0,2)*(x1[2]-x0[2]);
   b[0][0] = -COEFF_BN(pel,0,0,0)/vol;
   b[0][1] = -COEFF_BN(pel,0,0,1)/vol;
   b[0][2] = -COEFF_BN(pel,0,0,2)/vol;
   b[1][0] = -COEFF_BN(pel,0,1,0)/vol;
   b[1][1] = -COEFF_BN(pel,0,1,1)/vol;
   b[1][2] = -COEFF_BN(pel,0,1,2)/vol;
   b[2][0] = -COEFF_BN(pel,0,2,0)/vol;
   b[2][1] = -COEFF_BN(pel,0,2,1)/vol;
   b[2][2] = -COEFF_BN(pel,0,2,2)/vol;
   b[0][3] = -DOT(b[0],pel->n[1]->myvertex->x);
   b[1][3] = -DOT(b[1],pel->n[2]->myvertex->x);
   b[2][3] = -DOT(b[2],pel->n[3]->myvertex->x);
}

#endif

#if DATA_STR & LG_DATA

void norm_vect(pel,i,n)
ELEMENT *pel;
INT i;
FLOAT *n;
{
   FLOAT q = FMULT/COEFF_BF(pel,0,i)/20.0;
   
   n[0] = -COEFF_BN(pel,0,i,0)*q;
   n[1] = -COEFF_BN(pel,0,i,1)*q;
   n[2] = -COEFF_BN(pel,0,i,2)*q;
   if (pel->n[i]->lgd){
      if (i == 0) normal_vector(MYVERTEX(pel->n[1])->x,MYVERTEX(pel->n[2])->x,
                                            MYVERTEX(pel->n[3])->x,pel->f[0],n);
      if (i == 1) normal_vector(MYVERTEX(pel->n[0])->x,MYVERTEX(pel->n[2])->x,
                                            MYVERTEX(pel->n[3])->x,pel->f[1],n);
      if (i == 2) normal_vector(MYVERTEX(pel->n[0])->x,MYVERTEX(pel->n[1])->x,
                                            MYVERTEX(pel->n[3])->x,pel->f[2],n);
      if (i == 3) normal_vector(MYVERTEX(pel->n[0])->x,MYVERTEX(pel->n[1])->x,
                                            MYVERTEX(pel->n[2])->x,pel->f[3],n);
   }
}

void bar_coord(pel,b)
ELEMENT *pel;
FLOAT b[4][4];
{
   FLOAT *x0, *x1, vol;
   
   if (pel->n[0]->lgd || pel->n[1]->lgd || 
       pel->n[2]->lgd || pel->n[3]->lgd)
      barycentric_coordinates(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                           pel->n[2]->myvertex->x,pel->n[3]->myvertex->x,b);
   else{
      x0 = pel->n[0]->myvertex->x;
      x1 = pel->n[1]->myvertex->x;
      vol = COEFF_BN(pel,0,0,0)*(x1[0]-x0[0]) + 
            COEFF_BN(pel,0,0,1)*(x1[1]-x0[1]) + 
            COEFF_BN(pel,0,0,2)*(x1[2]-x0[2]);
      b[0][0] = -COEFF_BN(pel,0,0,0)/vol;
      b[0][1] = -COEFF_BN(pel,0,0,1)/vol;
      b[0][2] = -COEFF_BN(pel,0,0,2)/vol;
      b[1][0] = -COEFF_BN(pel,0,1,0)/vol;
      b[1][1] = -COEFF_BN(pel,0,1,1)/vol;
      b[1][2] = -COEFF_BN(pel,0,1,2)/vol;
      b[2][0] = -COEFF_BN(pel,0,2,0)/vol;
      b[2][1] = -COEFF_BN(pel,0,2,1)/vol;
      b[2][2] = -COEFF_BN(pel,0,2,2)/vol;
      b[0][3] = -DOT(b[0],pel->n[1]->myvertex->x);
      b[1][3] = -DOT(b[1],pel->n[2]->myvertex->x);
      b[2][3] = -DOT(b[2],pel->n[3]->myvertex->x);
   }
}

#endif

void restrict_from_faces_of_son(pel,faces,b,n,f,g)
ELEMENT *pel;
FACE **faces;                         /*  faces: faces of father of pel       */
FLOAT b[4][4], n[4][3];               /*  n[i] ... normal vector to faces[i]  */
INT f, g;
{
   FLOAT *x0, *x1, *x2, *x3, l00, l01, l02, l03, l10, l11, l12, l13, 
         l20, l21, l22, l23, l30, l31, l32, l33,  /* lij = \lambda_i(x_j) */
         nn[3];
   
   if (COEFF_BF(pel,0,0) < 0 || COEFF_BF(pel,0,1) < 0 || 
       COEFF_BF(pel,0,2) < 0 || COEFF_BF(pel,0,3) < 0){
      x0 = pel->n[0]->myvertex->x;
      x1 = pel->n[1]->myvertex->x;
      x2 = pel->n[2]->myvertex->x;
      x3 = pel->n[3]->myvertex->x;
      l00 = DOT(b[0],x0) + b[0][3];
      l01 = DOT(b[0],x1) + b[0][3];
      l02 = DOT(b[0],x2) + b[0][3];
      l03 = DOT(b[0],x3) + b[0][3];
      l10 = DOT(b[1],x0) + b[1][3];
      l11 = DOT(b[1],x1) + b[1][3];
      l12 = DOT(b[1],x2) + b[1][3];
      l13 = DOT(b[1],x3) + b[1][3];
      l20 = DOT(b[2],x0) + b[2][3];
      l21 = DOT(b[2],x1) + b[2][3];
      l22 = DOT(b[2],x2) + b[2][3];
      l23 = DOT(b[2],x3) + b[2][3];
      l30 = 1.0 - l00 - l10 - l20;
      l31 = 1.0 - l01 - l11 - l21;
      l32 = 1.0 - l02 - l12 - l22;
      l33 = 1.0 - l03 - l13 - l23;
      if (COEFF_BF(pel,0,0) < 0){
         norm_vect(pel,0,nn);
         if (IS_FF(faces[0]))
            FD(faces[0],g) += DOT(nn,n[0])*
                   integr4(l11,l12,l13,l21,l22,l23,l31,l32,l33)*FD(pel->f[0],f);
         if (IS_FF(faces[1]))
            FD(faces[1],g) += DOT(nn,n[1])*
                   integr4(l01,l02,l03,l21,l22,l23,l31,l32,l33)*FD(pel->f[0],f);
         if (IS_FF(faces[2]))
            FD(faces[2],g) += DOT(nn,n[2])*
                   integr4(l01,l02,l03,l11,l12,l13,l31,l32,l33)*FD(pel->f[0],f);
         if (IS_FF(faces[3]))
            FD(faces[3],g) += DOT(nn,n[3])*
                   integr4(l01,l02,l03,l11,l12,l13,l21,l22,l23)*FD(pel->f[0],f);
      }
      if (COEFF_BF(pel,0,1) < 0){
         norm_vect(pel,1,nn);
         if (IS_FF(faces[0]))
            FD(faces[0],g) += DOT(nn,n[0])*
                   integr4(l10,l12,l13,l20,l22,l23,l30,l32,l33)*FD(pel->f[1],f);
         if (IS_FF(faces[1]))
            FD(faces[1],g) += DOT(nn,n[1])*
                   integr4(l00,l02,l03,l20,l22,l23,l30,l32,l33)*FD(pel->f[1],f);
         if (IS_FF(faces[2]))
            FD(faces[2],g) += DOT(nn,n[2])*
                   integr4(l00,l02,l03,l10,l12,l13,l30,l32,l33)*FD(pel->f[1],f);
         if (IS_FF(faces[3]))
            FD(faces[3],g) += DOT(nn,n[3])*
                   integr4(l00,l02,l03,l10,l12,l13,l20,l22,l23)*FD(pel->f[1],f);
      }
      if (COEFF_BF(pel,0,2) < 0){
         norm_vect(pel,2,nn);
         if (IS_FF(faces[0]))
            FD(faces[0],g) += DOT(nn,n[0])*
                   integr4(l10,l11,l13,l20,l21,l23,l30,l31,l33)*FD(pel->f[2],f);
         if (IS_FF(faces[1]))
            FD(faces[1],g) += DOT(nn,n[1])*
                   integr4(l00,l01,l03,l20,l21,l23,l30,l31,l33)*FD(pel->f[2],f);
         if (IS_FF(faces[2]))
            FD(faces[2],g) += DOT(nn,n[2])*
                   integr4(l00,l01,l03,l10,l11,l13,l30,l31,l33)*FD(pel->f[2],f);
         if (IS_FF(faces[3]))
            FD(faces[3],g) += DOT(nn,n[3])*
                   integr4(l00,l01,l03,l10,l11,l13,l20,l21,l23)*FD(pel->f[2],f);
      }
      if (COEFF_BF(pel,0,3) < 0){
         norm_vect(pel,3,nn);
         if (IS_FF(faces[0]))
            FD(faces[0],g) += DOT(nn,n[0])*
                   integr4(l10,l11,l12,l20,l21,l22,l30,l31,l32)*FD(pel->f[3],f);
         if (IS_FF(faces[1]))
            FD(faces[1],g) += DOT(nn,n[1])*
                   integr4(l00,l01,l02,l20,l21,l22,l30,l31,l32)*FD(pel->f[3],f);
         if (IS_FF(faces[2]))
            FD(faces[2],g) += DOT(nn,n[2])*
                   integr4(l00,l01,l02,l10,l11,l12,l30,l31,l32)*FD(pel->f[3],f);
         if (IS_FF(faces[3]))
            FD(faces[3],g) += DOT(nn,n[3])*
                   integr4(l00,l01,l02,l10,l11,l12,l20,l21,l22)*FD(pel->f[3],f);
      }
   }
}

void prolong_to_faces_of_son(pel,faces,b,n,u,v)
ELEMENT *pel;
FACE **faces;                          /*  faces: faces of father of pel       */
FLOAT b[4][4], n[4][3];                /*  n[i] ... normal vector to faces[i]  */
INT u, v;
{
   FLOAT *x0, *x1, *x2, *x3, l00, l01, l02, l03, l10, l11, l12, l13, 
         l20, l21, l22, l23, l30, l31, l32, l33,  /* lij = \lambda_i(x_j) */
         nn[3], sum;
   
   if (COEFF_BF(pel,0,0) < 0 || COEFF_BF(pel,0,1) < 0 || 
       COEFF_BF(pel,0,2) < 0 || COEFF_BF(pel,0,3) < 0){
      x0 = pel->n[0]->myvertex->x;
      x1 = pel->n[1]->myvertex->x;
      x2 = pel->n[2]->myvertex->x;
      x3 = pel->n[3]->myvertex->x;
      l00 = DOT(b[0],x0) + b[0][3];
      l01 = DOT(b[0],x1) + b[0][3];
      l02 = DOT(b[0],x2) + b[0][3];
      l03 = DOT(b[0],x3) + b[0][3];
      l10 = DOT(b[1],x0) + b[1][3];
      l11 = DOT(b[1],x1) + b[1][3];
      l12 = DOT(b[1],x2) + b[1][3];
      l13 = DOT(b[1],x3) + b[1][3];
      l20 = DOT(b[2],x0) + b[2][3];
      l21 = DOT(b[2],x1) + b[2][3];
      l22 = DOT(b[2],x2) + b[2][3];
      l23 = DOT(b[2],x3) + b[2][3];
      l30 = 1.0 - l00 - l10 - l20;
      l31 = 1.0 - l01 - l11 - l21;
      l32 = 1.0 - l02 - l12 - l22;
      l33 = 1.0 - l03 - l13 - l23;
      if (COEFF_BF(pel,0,0) < 0){
         norm_vect(pel,0,nn);
         sum = 0.0;
         if (IS_FF(faces[0]))
            sum += DOT(nn,n[0])*
                    integr4(l11,l12,l13,l21,l22,l23,l31,l32,l33)*FD(faces[0],u);
         if (IS_FF(faces[1]))
            sum += DOT(nn,n[1])*
                    integr4(l01,l02,l03,l21,l22,l23,l31,l32,l33)*FD(faces[1],u);
         if (IS_FF(faces[2]))
            sum += DOT(nn,n[2])*
                    integr4(l01,l02,l03,l11,l12,l13,l31,l32,l33)*FD(faces[2],u);
         if (IS_FF(faces[3]))
            sum += DOT(nn,n[3])*
                    integr4(l01,l02,l03,l11,l12,l13,l21,l22,l23)*FD(faces[3],u);
         FD(pel->f[0],v) = sum;
      }
      if (COEFF_BF(pel,0,1) < 0){
         norm_vect(pel,1,nn);
         sum = 0.0;
         if (IS_FF(faces[0]))
            sum += DOT(nn,n[0])*
                    integr4(l10,l12,l13,l20,l22,l23,l30,l32,l33)*FD(faces[0],u);
         if (IS_FF(faces[1]))
            sum += DOT(nn,n[1])*
                    integr4(l00,l02,l03,l20,l22,l23,l30,l32,l33)*FD(faces[1],u);
         if (IS_FF(faces[2]))
            sum += DOT(nn,n[2])*
                    integr4(l00,l02,l03,l10,l12,l13,l30,l32,l33)*FD(faces[2],u);
         if (IS_FF(faces[3]))
            sum += DOT(nn,n[3])*
                    integr4(l00,l02,l03,l10,l12,l13,l20,l22,l23)*FD(faces[3],u);
         FD(pel->f[1],v) = sum;
      }
      if (COEFF_BF(pel,0,2) < 0){
         norm_vect(pel,2,nn);
         sum = 0.0;
         if (IS_FF(faces[0]))
            sum += DOT(nn,n[0])*
                    integr4(l10,l11,l13,l20,l21,l23,l30,l31,l33)*FD(faces[0],u);
         if (IS_FF(faces[1]))
            sum += DOT(nn,n[1])*
                    integr4(l00,l01,l03,l20,l21,l23,l30,l31,l33)*FD(faces[1],u);
         if (IS_FF(faces[2]))
            sum += DOT(nn,n[2])*
                    integr4(l00,l01,l03,l10,l11,l13,l30,l31,l33)*FD(faces[2],u);
         if (IS_FF(faces[3]))
            sum += DOT(nn,n[3])*
                    integr4(l00,l01,l03,l10,l11,l13,l20,l21,l23)*FD(faces[3],u);
         FD(pel->f[2],v) = sum;
      }
      if (COEFF_BF(pel,0,3) < 0){
         norm_vect(pel,3,nn);
         sum = 0.0;
         if (IS_FF(faces[0]))
            sum += DOT(nn,n[0])*
                    integr4(l10,l11,l12,l20,l21,l22,l30,l31,l32)*FD(faces[0],u);
         if (IS_FF(faces[1]))
            sum += DOT(nn,n[1])*
                    integr4(l00,l01,l02,l20,l21,l22,l30,l31,l32)*FD(faces[1],u);
         if (IS_FF(faces[2]))
            sum += DOT(nn,n[2])*
                    integr4(l00,l01,l02,l10,l11,l12,l30,l31,l32)*FD(faces[2],u);
         if (IS_FF(faces[3]))
            sum += DOT(nn,n[3])*
                    integr4(l00,l01,l02,l10,l11,l12,l20,l21,l22)*FD(faces[3],u);
         FD(pel->f[3],v) = sum;
      }
   }
}

void treat_cubic_functions(pelem,treat_faces_of_son,f,g)
ELEMENT *pelem;
void (*treat_faces_of_son)();
INT f, g;
{
   FLOAT b[4][4], n[4][3];
   INT i;
   
   bar_coord(pelem,b);
   norm_vect(pelem,0,n[0]);
   norm_vect(pelem,1,n[1]);
   norm_vect(pelem,2,n[2]);
   norm_vect(pelem,3,n[3]);
   for (i = 0; i < 8; i++)
      (*treat_faces_of_son)(pelem->sons[i],pelem->f,b,n,f,g);
}

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void node_restriction(theGrid,f,g)/* (g,theGrid) := restrict(f,theGrid->finer)*/
GRID *theGrid;
INT f, g;
{
   NODE *pnode;
   LINK *plink;
   FLOAT sum0, sum1, sum2;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      sum0 = sum1 = sum2 = 0.0;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         sum0 += ND(NBNODE(plink),f,0);
         sum1 += ND(NBNODE(plink),f,1);
         sum2 += ND(NBNODE(plink),f,2);
      }
      ND(pnode,g,0) = ND(pnode->son,f,0) + sum0*0.5;
      ND(pnode,g,1) = ND(pnode->son,f,1) + sum1*0.5;
      ND(pnode,g,2) = ND(pnode->son,f,2) + sum2*0.5;
   }
}

#else

void node_restriction(theGrid,f,g)/* (g,theGrid) := restrict(f,theGrid->finer)*/
GRID *theGrid; INT f, g;
{  eprintf("Error: node_restriction not available.\n");  }

#endif

#if F_DATA & SCALAR_FACE_DATA

void face_restriction(theGrid,f,g)/* (g,theGrid) := restrict(f,theGrid->finer)*/
GRID *theGrid;
INT f, g;
{
   ELEMENT *pelem;
   FACE *pface;
   
   for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
      FD(pface,g) =0.0;
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      treat_cubic_functions(pelem,restrict_from_faces_of_son,f,g);
}

void face_prolongation(theGrid,u,v)/* (v,theGrid->finer) := prolong(u,theGrid)*/
GRID *theGrid;
INT u, v;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      treat_cubic_functions(pelem,prolong_to_faces_of_son,u,v);
}

#else

void face_restriction(theGrid,f,g)
GRID *theGrid; INT f, g;
{  eprintf("Error: face_restriction not available.\n");  }

void face_prolongation(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: face_prolongation not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

/* saves restriction of u in v; theGrid is the coarser grid */
void restriction_of_u(theGrid,u,v)
GRID *theGrid;
INT u, v;
{
   NODE *pnode;
   
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      ND(pnode,v,0) = ND(pnode->son,u,0);
      ND(pnode,v,1) = ND(pnode->son,u,1);
      ND(pnode,v,2) = ND(pnode->son,u,2);
   }
}

void node_prolongation(theGrid,u,v)  /*  (v,theGrid->finer) := prolong(u,theGrid)  */
GRID *theGrid;
INT u, v;
{
   NODE *pnode;
   LINK *plink;
   FLOAT val0, val1, val2;
   
   for (pnode = FIRSTNODE(theGrid->finer); pnode != NULL; pnode = pnode->succ)
      if (pnode->father){
         ND(pnode,v,0) = ND(pnode->father,u,0);
         ND(pnode,v,1) = ND(pnode->father,u,1);
         ND(pnode,v,2) = ND(pnode->father,u,2);
      }
      else
         ND(pnode,v,0) = ND(pnode,v,1) = ND(pnode,v,2) = 0.0;
   for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ){
      val0 = ND(pnode,u,0)*0.5;
      val1 = ND(pnode,u,1)*0.5;
      val2 = ND(pnode,u,2)*0.5;
      for (plink = START(pnode->son); plink != NULL; plink = plink->next){
         ND(NBNODE(plink),v,0) += val0;
         ND(NBNODE(plink),v,1) += val1;
         ND(NBNODE(plink),v,2) += val2;
      }
   }
}

#else

void restriction_of_u(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: restriction_of_u not available.\n");  }

void node_prolongation(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: node_prolongation not available.\n");  }

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void elem_restriction(theGrid,f,g)/* (g,theGrid) := restrict(f,theGrid->finer)*/
GRID *theGrid;
INT f, g;
{
   ELEMENT *pel;
   FLOAT sum=0;
   
   for (pel = SUCC(FIRSTELEMENT(theGrid->finer)); pel != NULL; pel = pel->succ)
      sum += ED(pel,f);
   ED(FIRSTELEMENT(theGrid->finer),f) = -sum;
   for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
      ED(pel,g) =   ED(pel->sons[0],f) + ED(pel->sons[1],f) + ED(pel->sons[2],f)
                  + ED(pel->sons[3],f) + ED(pel->sons[4],f) + ED(pel->sons[5],f)
                  + ED(pel->sons[6],f) + ED(pel->sons[7],f);
}

void elem_prolongation(theGrid,u,v)/* (v,theGrid->finer) := prolong(u,theGrid)*/
GRID *theGrid;
INT u, v;
{
   ELEMENT *pelem;

   ED(FIRSTELEMENT(theGrid),u) = 0.;
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
      ED(pelem->sons[0],v) = ED(pelem->sons[1],v) = ED(pelem->sons[2],v) = 
      ED(pelem->sons[3],v) = ED(pelem->sons[4],v) = ED(pelem->sons[5],v) = 
      ED(pelem->sons[6],v) = ED(pelem->sons[7],v) = ED(pelem,u);
   if (T_FOR_P & WITHOUT_FIRST)
      subtr_first_element(theGrid->finer,v);
}

#else

void elem_restriction(theGrid,f,g)/* (g,theGrid) := restrict(f,theGrid->finer)*/
GRID *theGrid; INT f, g;
{  eprintf("Error: elem_restriction not available.\n");  }

void elem_prolongation(theGrid,u,v)
GRID *theGrid; INT u, v;
{  eprintf("Error: elem_prolongation not available.\n");  }

#endif

#if !(U_SPACE == P1C_NEW_FBUB)

/* theGrid is the coarse grid */
void restrict_rhs(theGrid,Z,f_fine,g_fine,f_coarse,g_coarse,xe)
GRID *theGrid;
INT f_fine, g_fine, f_coarse, g_coarse, xe;
{
   node_restriction(theGrid,f_fine,f_coarse);
   face_restriction(theGrid,f_fine,f_coarse);
   elem_restriction(theGrid,g_fine,g_coarse);
}

/* theGrid is the coarse grid */
void prolong_solution(theGrid,u_coarse,p_coarse,u_fine,p_fine)
GRID *theGrid;
INT u_coarse, p_coarse, u_fine, p_fine;
{
   node_prolongation(theGrid,u_coarse,u_fine);
   face_prolongation(theGrid,u_coarse,u_fine);
   elem_prolongation(theGrid,p_coarse,p_fine);
}

void correction(tGrid,u,ue,d,de,Z,f)
GRID *tGrid;
INT u, ue, d, de, Z, f;
{
   vs_subtr_nf(tGrid,u,d,u,1);
   subtr_e(tGrid,ue,de,ue);
}

#endif

void exact_face_solution(tGrid,Z,f,pe,u,q)
GRID *tGrid;
INT Z, f, pe, u, q;  /* only boundary components of q */
{
   mult_BT(tGrid,0,pe,u,q,WITHOUT_FIRST,0,Q_SE,Q_SF,0);
   subtr_f(tGrid,f,u,u,1);
   mult_by_inverse_face_diagonal_sf(tGrid,Z,u,u);
}

#if U_SPACE == P1C_NEW_FBUB

/* theGrid is the coarse grid */
void restrict_rhs(theGrid,Z,f_fine,g_fine,f_coarse,g_coarse,xe)
GRID *theGrid;
INT f_fine, g_fine, f_coarse, g_coarse, xe;
{
/*
   mult_by_inverse_face_diagonal_sf(theGrid->finer,Z,f_fine,f_fine);
   multB_f(theGrid->finer,f_fine,xe);
   subtr_e(theGrid->finer,g_fine,xe,xe);
   node_restriction(theGrid,f_fine,f_coarse);
   set_value_f(theGrid,0.,f_coarse,1);
   elem_restriction(theGrid,xe,g_coarse);
*/
}

/* theGrid is the coarse grid */
void prolong_solution(theGrid,u_coarse,p_coarse,u_fine,p_fine)
GRID *theGrid;
INT u_coarse, p_coarse, u_fine, p_fine;
{
   node_prolongation(theGrid,u_coarse,u_fine);
   elem_prolongation(theGrid,p_coarse,p_fine);
}

void correction(tGrid,u,ue,d,de,Z,f)
GRID *tGrid;
INT u, ue, d, de, Z, f;
{
   v_subtr_n(tGrid,u,d,u,1);
   subtr_e(tGrid,ue,de,ue);
   exact_face_solution(tGrid,Z,f,ue,u,d);
}

#endif

#endif  /* DIM == 3 */
