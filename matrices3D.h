/******************************************************************************/
/*                                                                            */
/*                           3D stiffness matrices                            */
/*                                                                            */
/******************************************************************************/

/*  stiffness matrix corresponding to nu*\int_\om \nabla u\cdot\nabla v \dx   */

#if N_DATA & ONE_NODE_MATR

void putai3(n0,n1,n2,a0,a1,a2,Z)
NODE *n0, *n1, *n2;
FLOAT a0, a1, a2;
INT Z;
{
   COEFFN(n0,Z) += a0;
   COEFFN(n1,Z) += a1;
   COEFFN(n2,Z) += a2;
}

void putaij(plin,n1,n2,an1,an2,Z)
LINK *plin;
NODE *n1, *n2;
FLOAT an1, an2;
INT Z;
{
   LINK *pli;
   
   for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2; pli=pli->next);
   if (pli->nbnode == n1){
      COEFFL(pli,Z) += an1;
      for (pli=pli->next; pli->nbnode != n2; pli=pli->next); 
      COEFFL(pli,Z) += an2;
   }
   else{
      COEFFL(pli,Z) += an2;
      for (pli=pli->next; pli->nbnode != n1; pli=pli->next); 
      COEFFL(pli,Z) += an1;
   }
}

void putaij3(plin,n1,n2,n3,an1,an2,an3,Z)
LINK *plin;
NODE *n1, *n2, *n3;
FLOAT an1, an2, an3;
INT Z;
{
   LINK *pli;
   
   for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2 && pli->nbnode != n3; 
                                                                 pli=pli->next);
   if (pli->nbnode == n1){
      COEFFL(pli,Z) += an1;
      putaij(pli->next,n2,n3,an2,an3,Z);
   }
   else if (pli->nbnode == n2){
      COEFFL(pli,Z) += an2;
      putaij(pli->next,n1,n3,an1,an3,Z);
   }
   else{
      COEFFL(pli,Z) += an3;
      putaij(pli->next,n1,n2,an1,an2,Z);
   }
}

#else

void putai3(n0,n1,n2,a0,a1,a2,Z)
NODE *n0, *n1, *n2; FLOAT a0, a1, a2; INT Z;
{  eprintf("Error: putai3 not available.\n");  }

void putaij(plin,n1,n2,an1,an2,Z)
LINK *plin; NODE *n1, *n2; FLOAT an1, an2; INT Z;
{  eprintf("Error: putaij not available.\n");  }

void putaij3(plin,n1,n2,n3,an1,an2,an3,Z)
LINK *plin; NODE *n1, *n2, *n3; FLOAT an1, an2, an3; INT Z;
{  eprintf("Error: putaij3 not available.\n");  }

#endif
 
#if N_DATA & NxN_NODE_MATR

void putai4g(n0,n1,n2,n3,a0,a1,a2,a3,Z)
NODE *n0, *n1, *n2, *n3;
FLOAT a0, a1, a2, a3;
INT Z;
{
   COEFF_NN(n0,Z,0,0) += a0;
   COEFF_NN(n1,Z,0,0) += a1;
   COEFF_NN(n2,Z,0,0) += a2;
   COEFF_NN(n3,Z,0,0) += a3;
}


void putaijg(plin,n1,n2,an1,an2,Z)
LINK *plin;
NODE *n1, *n2;
FLOAT an1, an2;
INT Z;
{
   LINK *pli;
   
   for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2; pli=pli->next);
   if (pli->nbnode == n1){
      COEFF_NN(pli,Z,0,0) += an1;
      for (pli=pli->next; pli->nbnode != n2; pli=pli->next); 
      COEFF_NN(pli,Z,0,0) += an2;
   }
   else{
      COEFF_NN(pli,Z,0,0) += an2;
      for (pli=pli->next; pli->nbnode != n1; pli=pli->next); 
      COEFF_NN(pli,Z,0,0) += an1;
   }
}

#else

void putai4g(n0,n1,n2,n3,a0,a1,a2,a3,Z)
NODE *n0, *n1, *n2, *n3; FLOAT a0, a1, a2, a3; INT Z;
{  eprintf("Error: putai4g not available.\n");  }

void putaijg(plin,n1,n2,an1,an2,Z)
LINK *plin; NODE *n1, *n2; FLOAT an1, an2; INT Z;
{  eprintf("Error: putaijg not available.\n");  }

#endif
 
#if (N_DATA & DxD_NODE_MATR) && (DIM == 2)

 void vputaij(plin,n1,n2,an1,an2,Z)
 LINK *plin;
 NODE *n1, *n2;
 FLOAT an1, an2;
 INT Z;
 {
   LINK *pli;
   
   for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2; pli=pli->next);
   if (pli->nbnode == n1){
      SET_COEFFNN(pli,Z,an1,0.,0.,an1);
      for (pli=pli->next; pli->nbnode != n2; pli=pli->next); 
      SET_COEFFNN(pli,Z,an2,0.,0.,an2);
   }
   else {
      SET_COEFFNN(pli,Z,an2,0.,0.,an2);
      for (pli=pli->next; pli->nbnode != n1; pli=pli->next); 
      SET_COEFFNN(pli,Z,an1,0.,0.,an1);
   }
 }

#else

 void vputaij(plin,n1,n2,an1,an2,Z)
 LINK *plin; NODE *n1, *n2; FLOAT an1, an2; INT Z;
 {  eprintf("Error: vputaij not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_FACES) && (DIM == 2)

void vputasij(pnfl,fa1,fa2,pn1,pn2,Z)
NFLINK *pnfl; 
FACE *fa1, *fa2;
INT Z;
FLOAT pn1,pn2;
{
  NFLINK *pnf;
  
  for (pnf = pnfl; pnf->nbface != fa1 && pnf->nbface != fa2; pnf = pnf->next);
  if (pnf->nbface==fa1){
     SET_COEFFNN(pnf,Z,pn1,0.,0.,pn1);
     for (pnf = pnf->next; pnf->nbface != fa2; pnf = pnf->next); 
     SET_COEFFNN(pnf,Z,pn2,0.,0.,pn2);
  }
  else {
     SET_COEFFNN(pnf,Z,pn2,0.,0.,pn2);
     for (pnf = pnf->next; pnf->nbface != fa1; pnf = pnf->next); 
     SET_COEFFNN(pnf,Z,pn1,0.,0.,pn1);
  }
}

#else
 
void vputasij(pnfl,fa1,fa2,pn1,pn2,Z)
NFLINK *pnfl; FACE *fa1, *fa2; INT Z; FLOAT pn1,pn2;
{  eprintf("Error: vputasij not available.\n");  }

#endif

#if (DATA_S & F_LINK_TO_FACES) && (F_DATA & DxD_FACE_MATR) && (DIM == 2)

void vputppij(pflink,fa1,fa2,pp1,pp2,Z)
FLINK *pflink;
FACE *fa1, *fa2;
FLOAT pp1, pp2;
INT Z;
{
   FLINK *pfli;
  
   for (pfli=pflink; pfli->nbface!=fa1 && pfli->nbface!=fa2; pfli=pfli->next);
   if (pfli->nbface==fa1){
      SET_COEFFNN(pfli,Z,pp1,0.,0.,pp1);
      for (pfli=pfli->next; pfli->nbface!=fa2; pfli=pfli->next); 
      SET_COEFFNN(pfli,Z,pp2,0.,0.,pp2);
   }
   else {
      SET_COEFFNN(pfli,Z,pp2,0.,0.,pp2);
      for (pfli=pfli->next; pfli->nbface!=fa1; pfli=pfli->next); 
      SET_COEFFNN(pfli,Z,pp1,0.,0.,pp1);
   }
}
 
#else

void vputppij(pflink,fa1,fa2,pp1,pp2,Z)
FLINK *pflink; FACE *fa1, *fa2; FLOAT pp1, pp2; INT Z;
{  eprintf("Error: vputppij not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_FACES)

void putasij(pnfl,fa1,fa2,pn1,pn2,Z)
NFLINK *pnfl; 
FACE *fa1, *fa2;
INT Z;
FLOAT pn1,pn2;
{
  NFLINK *pnf;
  
  for (pnf = pnfl; pnf->nbface != fa1 && pnf->nbface != fa2; pnf = pnf->next);
  if (pnf->nbface==fa1){
     COEFFL(pnf,Z) += pn1;
     for (pnf = pnf->next; pnf->nbface != fa2; pnf = pnf->next); 
     COEFFL(pnf,Z) += pn2;
  }
  else {
     COEFFL(pnf,Z) += pn2;
     for (pnf = pnf->next; pnf->nbface != fa1; pnf = pnf->next); 
     COEFFL(pnf,Z) += pn1;
  }
}
 
#else

void putasij(pnfl,fa1,fa2,pn1,pn2,Z)
NFLINK *pnfl; FACE *fa1, *fa2; INT Z; FLOAT pn1,pn2;
{  eprintf("Error: putasij not available.\n");  }

#endif
 
#if !(DATA_STR & KORN_MATRIX)
 
 void Putaij(plin,n1,n2,an1,an2,Z)
 LINK *plin;
 NODE *n1, *n2;
 FLOAT an1, an2;
 INT Z;
 {
    Putaij(plin,n1,n2,an1,an2,Z);
 }

#if (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_FACES)

void putapij(pnfl,fa1,fa2,pn1,pn2,Z,nn1,nn2)
NFLINK *pnfl; 
FACE *fa1, *fa2;
INT Z;
FLOAT pn1,pn2,nn1[DIM],nn2[DIM];
{
  NFLINK *pnf;
  
  for (pnf = pnfl; pnf->nbface != fa1 && pnf->nbface != fa2; pnf = pnf->next);
  if (pnf->nbface==fa1){
     SET4(COEFF_NFP(pnf,Z),nn1,pn1)
     for (pnf = pnf->next; pnf->nbface != fa2; pnf = pnf->next); 
     SET4(COEFF_NFP(pnf,Z),nn2,pn2)
  }
  else {
     SET4(COEFF_NFP(pnf,Z),nn2,pn2)
     for (pnf = pnf->next; pnf->nbface != fa1; pnf = pnf->next); 
     SET4(COEFF_NFP(pnf,Z),nn1,pn1)
  }
}
 
#endif

#if (N_DATA & VECTOR_NODE_DATA) && (DATA_S & F_LINK_TO_NODES) && (F_DATA & IxD_FACE_NODE_MATR)

void putpaij(fan,n1,p1n,nn,Z,f,u)
FACE *fan;
NODE *n1;
FLOAT p1n, nn[DIM];
INT Z, f, u;
{
   FNLINK *pfn;
   
   if (IS_FN(n1)){
      for (pfn=fan->fnstart; pfn->nbnode!=n1; pfn=pfn->next);
      SET4(COEFF_FNP(pfn,Z),nn,p1n)
   }
   else if (IS_DN(n1))
      FD(fan,f) -= p1n*DOT(nn,NDD(n1,u));
}

#endif

#if N_DATA & VECTOR_NODE_DATA

void u0_to_f(pnode,f,n1,u,q)
NODE *pnode, *n1;
INT f,u;
FLOAT q;
{
   if (IS_DN(n1))
      SET9(NDD(pnode,f),NDD(n1,u),q)
}

#else

void u0_to_f(pnode,f,n1,u,q)
NODE *pnode, *n1; INT f,u; FLOAT q;
{  eprintf("Error: u0_to_f not available.\n");  }

#endif

#endif  /* !(DATA_STR & KORN_MATRIX) */

#if (DATA_S & F_LINK_TO_FACES) && (F_DATA & ONE_FACE_MATR)

 void putppij(pflink,fa1,fa2,pp1,pp2,Z)
 FLINK *pflink;
 FACE *fa1, *fa2;
 FLOAT pp1, pp2;
 INT Z;
 {
   FLINK *pfli;
   
   for (pfli=pflink; pfli->nbface!=fa1 && pfli->nbface!=fa2; pfli=pfli->next);
   if (pfli->nbface==fa1){
      COEFF_FL(pfli,Z) += pp1;
      for (pfli=pfli->next; pfli->nbface!=fa2; pfli=pfli->next); 
      COEFF_FL(pfli,Z) += pp2;
   }
   else {
      COEFF_FL(pfli,Z) += pp2;
      for (pfli=pfli->next; pfli->nbface!=fa1; pfli=pfli->next); 
      COEFF_FL(pfli,Z) += pp1;
   }
 }
 
#endif
 
#if DIM == 3

/* computes approximately \int_K f0 l_1\dx (exactly for linear f0) */
FLOAT integr1(x1,x2,x3,x4,f0,detB)         /*  detB = volume of K  */
FLOAT (*f0)();
FLOAT x1[DIM],x2[DIM],x3[DIM],x4[DIM],detB;
{
   return( (2.0* (*f0)(x1) + (*f0)(x2) + (*f0)(x3) + (*f0)(x4))*detB/20.0 );
}
 
 /* computes approximately \int_K f0\cdot p^0\dx (exactly for linear f0), where
    p^0=l_1*l_2*l_3*nn */
 FLOAT integr3(n0,n1,n2,n3,nn,detB)
 NODE *n0,*n1,*n2,*n3;
 FLOAT nn[DIM], detB;
 {
   FLOAT *x0, *x1, *x2, *x3;
   
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   
   return(( ( f01(x0)/2.0 + f01(x1) + f01(x2) + f01(x3) )*nn[0] +
            ( f02(x0)/2.0 + f02(x1) + f02(x2) + f02(x3) )*nn[1] +
            ( f03(x0)/2.0 + f03(x1) + f03(x2) + f03(x3) )*nn[2] )*detB/420.0);
 }
 
 /* analogous computation as integr3 on the tetrahedron (n1,n2,n3,xc)     */
 FLOAT integr3c(n0,n1,n2,n3,nn,detB)  /* bar. coord. w.r.t. (n1,n2,n3,xc) */
 NODE *n0,*n1,*n2,*n3;                /* detB = volume of (n1,n2,n3,n0)   */
 FLOAT nn[DIM], detB;
 {
   FLOAT *x0, *x1, *x2, *x3, xc[DIM];
   
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   POINT4(x0,x1,x2,x3,xc);
   
   return(( ( f01(xc)/2.0 + f01(x1) + f01(x2) + f01(x3) )*nn[0] +
            ( f02(xc)/2.0 + f02(x1) + f02(x2) + f02(x3) )*nn[1] +
            ( f03(xc)/2.0 + f03(x1) + f03(x2) + f03(x3) )*nn[2] )*detB/1680.0);
 }

/* 
 FLOAT integr3p(n0,n1,n2,n3,nn,detB)
 NODE *n0,*n1,*n2,*n3;
 FLOAT nn[DIM], detB;
 {
   FLOAT *x0, *x1, *x2, *x3, xc[DIM];
   
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   POINT4(x0,x1,x2,x3,xc);
   
   return(( ( fp1(xc)/2.0 + fp1(x1) + fp1(x2) + fp1(x3) )*nn[0] +
            ( fp2(xc)/2.0 + fp2(x1) + fp2(x2) + fp2(x3) )*nn[1] +
            ( fp3(xc)/2.0 + fp3(x1) + fp3(x2) + fp3(x3) )*nn[2] )*detB/1680.0);
 }
*/ 

/* computes approximately \int_K f0 l_1 l_2 l_3 l_4\dx (exactly for linear f0)*/
void integr3B(n0,n1,n2,n3,a,detB)
NODE *n0,*n1,*n2,*n3;
FLOAT a[DIM], detB;
{
   FLOAT *x0, *x1, *x2, *x3;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   
   a[0] = ( f01(x0) + f01(x1) + f01(x2) + f01(x3) )*detB/3360.;
   a[1] = ( f02(x0) + f02(x1) + f02(x2) + f02(x3) )*detB/3360.;
   a[2] = ( f03(x0) + f03(x1) + f03(x2) + f03(x3) )*detB/3360.;
}
 
#if !(DATA_STR & KORN_MATRIX)
 
#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA)

void nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB) /* grad u * grad v    */
NODE *n0, *n1, *n2, *n3;                               /* conforming pw. lin. */
INT Z, f, u;                                           /* b.c. in u           */
FLOAT ann, an1, an2, an3, rdetB;
{
   LINK *plin, *pli;
  
   ND(n0,f,0) += integr1(n0->myvertex->x,n1->myvertex->x,
                         n2->myvertex->x,n3->myvertex->x,f01,rdetB);
   ND(n0,f,1) += integr1(n0->myvertex->x,n1->myvertex->x,
                         n2->myvertex->x,n3->myvertex->x,f02,rdetB);
   ND(n0,f,2) += integr1(n0->myvertex->x,n1->myvertex->x,
                         n2->myvertex->x,n3->myvertex->x,f03,rdetB);
  
   COEFFN(n0,Z) += ann;
  
   plin=n0->start;
   if (IS_FN(n1) && IS_FN(n2) && IS_FN(n3)){
      for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2 &&
                     pli->nbnode != n3; pli=pli->next);
      if (pli->nbnode == n1){
         COEFFL(pli,Z) += an1;
         Putaij(pli->next,n2,n3,an2,an3,Z);
      }            
      else if(pli->nbnode == n2){
         COEFFL(pli,Z) += an2;
         Putaij(pli->next,n1,n3,an1,an3,Z);
      }
      else {    
         COEFFL(pli,Z) += an3;
         Putaij(pli->next,n1,n2,an1,an2,Z);
      }
   }
   else if (IS_FN(n1) && IS_FN(n2)){
      Putaij(plin,n1,n2,an1,an2,Z);
      u0_to_f(n0,f,n3,u,an3);
   }
   else if (IS_FN(n1) && IS_FN(n3)){
      Putaij(plin,n1,n3,an1,an3,Z);
      u0_to_f(n0,f,n2,u,an2);
   }
   else if (IS_FN(n2) && IS_FN(n3)){ 
      Putaij(plin,n2,n3,an2,an3,Z);
      u0_to_f(n0,f,n1,u,an1);
   }
   else if (IS_FN(n1)){
      for (pli=plin; pli->nbnode != n1; pli=pli->next);
      COEFFL(pli,Z) += an1;
      u0_to_f(n0,f,n2,u,an2);
      u0_to_f(n0,f,n3,u,an3);
   }                                         
   else if (IS_FN(n2)){
      for (pli=plin; pli->nbnode != n2; pli=pli->next);
      COEFFL(pli,Z) += an2; 
      u0_to_f(n0,f,n1,u,an1);
      u0_to_f(n0,f,n3,u,an3);                                    
   }                                       
   else if (IS_FN(n3)){
      for (pli=plin; pli->nbnode != n3; pli=pli->next);
      COEFFL(pli,Z) += an3;
      u0_to_f(n0,f,n1,u,an1);
      u0_to_f(n0,f,n2,u,an2);  
   }
   else{
      u0_to_f(n0,f,n1,u,an1);
      u0_to_f(n0,f,n2,u,an2);
      u0_to_f(n0,f,n3,u,an3);      
   }
}

#else

void nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB)
NODE *n0, *n1, *n2, *n3; INT Z, f, u; FLOAT ann, an1, an2, an3, rdetB;
{  eprintf("Error: nn_matrix not available.\n");  }

#endif

void laijb(n0,n1,n2,n3,Z,f,u,b0,b1,b2,b3,detB,rdetB)
NODE *n0, *n1, *n2, *n3;          /* detB = volume * nu;  rdetB = volume      */
INT Z, f, u;
FLOAT b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB;
{
   FLOAT ann, an1, an2, an3;
   
   if (IS_FN(n0)){  
      an1 = DOT(b0,b1)*detB;
      an2 = DOT(b0,b2)*detB;
      an3 = DOT(b0,b3)*detB;
      ann = DOT(b0,b0)*detB;
      nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB);
   }
}

#if !(DATA_STR & REDUCED) && !(U_SPACE == P1C_NEW_FBUB)

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) 

void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                   rdetB,fmult2)
NODE *n0, *n1, *n2, *n3;          /* detB = volume * nu;  rdetB = volume      */
FACE *fa0, *fa1, *fa2, *fa3;      /* NF matrix, FN matrix multiplied by FMULT */
INT Z, f, u;                      /* FF matrix multiplied by fmult2           */
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{
  NFLINK *pnf;
  FLINK *pflink, *pfli;
  FLOAT ann, an1, an2, an3, a11, a12, a13, a22, a23, a33, 
        pnn, pn1, pn2, pn3, p1n, p2n, p3n, 
        ppn, pp1, pp2, pp3;
  
  ann = DOT(b0,b0)*detB;
  an1 = DOT(b0,b1)*detB;
  an2 = DOT(b0,b2)*detB;
  an3 = DOT(b0,b3)*detB;
  
 if (IS_FN(n0)){  
  
  nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB);
  
  pnn = -ann/20.0*FMULT;
  pn1 = -an1/20.0*FMULT;
  pn2 = -an2/20.0*FMULT;
  pn3 = -an3/20.0*FMULT;
  
  for (pnf=n0->nfstart; pnf->nbface!=fa1 && pnf->nbface!=fa2 
                     && pnf->nbface!=fa3; pnf=pnf->next);
  if (pnf->nbface==fa1){
    SET4(COEFF_NFP(pnf,Z),nn1,pn1)
    putapij(pnf->next,fa2,fa3,pn2,pn3,Z,nn2,nn3);
  }
  else if (pnf->nbface==fa2){
    SET4(COEFF_NFP(pnf,Z),nn2,pn2)
    putapij(pnf->next,fa1,fa3,pn1,pn3,Z,nn1,nn3);
  }
  else {
    SET4(COEFF_NFP(pnf,Z),nn3,pn3)
    putapij(pnf->next,fa1,fa2,pn1,pn2,Z,nn1,nn2);
  } 
  
  if (NOT_FF(fa0))
    SET9(NDD(n0,f),nn0,FD(fa0,u)*pnn)
  else {
    for  (pnf=n0->nfstart; pnf->nbface!=fa0; pnf=pnf->next);
    SET4(COEFF_NFP(pnf,Z),nn0,pnn)
  }
 }
                                          
 if (IS_FF(fa0)){ 
 
  RHS_FOR_FACES 
  
  a11 = DOT(b1,b1)*detB;
  a12 = DOT(b1,b2)*detB;
  a13 = DOT(b1,b3)*detB;
  a22 = DOT(b2,b2)*detB;
  a23 = DOT(b2,b3)*detB;
  a33 = DOT(b3,b3)*detB;
  
  ppn = ( a11 + a22 + a33 + a12 + a13 + a23 )/210.0*fmult2;
  pp1 = ( an1 + an1 - a23 )*DOT(nn0,nn1)/420.0*fmult2;
  pp2 = ( an2 + an2 - a13 )*DOT(nn0,nn2)/420.0*fmult2; 
  pp3 = ( an3 + an3 - a12 )*DOT(nn0,nn3)/420.0*fmult2; 
  
  COEFF_FF(fa0,Z) += ppn;
  
  pflink = fa0->fstart;
  if (IS_FF(fa1) && IS_FF(fa2) && IS_FF(fa3)){
    for (pfli=pflink; pfli->nbface!=fa1 && pfli->nbface!=fa2 
                && pfli->nbface!=fa3 ; pfli=pfli->next);
    if (pfli->nbface==fa1){
       COEFF_FL(pfli,Z) += pp1;
       putppij(pfli->next,fa2,fa3,pp2,pp3,Z);
    }            
    else if(pfli->nbface==fa2){
       COEFF_FL(pfli,Z) += pp2;
       putppij(pfli->next,fa1,fa3,pp1,pp3,Z);
    }
    else {    
       COEFF_FL(pfli,Z) += pp3;
       putppij(pfli->next,fa1,fa2,pp1,pp2,Z);
    }
  }
  else if (IS_FF(fa1) && IS_FF(fa2)){
    putppij(pflink,fa1,fa2,pp1,pp2,Z);
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }
  else if (IS_FF(fa1) && IS_FF(fa3)){
    putppij(pflink,fa1,fa3,pp1,pp3,Z);
    FD(fa0,f) -= FD(fa2,u)*pp2;
  }
  else if (IS_FF(fa2) && IS_FF(fa3)){ 
    putppij(pflink,fa2,fa3,pp2,pp3,Z);
    FD(fa0,f) -= FD(fa1,u)*pp1;
  }
  else if (IS_FF(fa1)){
    for (pfli=pflink; pfli->nbface!=fa1; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }                                         
  else if (IS_FF(fa2)){
    for (pfli=pflink; pfli->nbface!=fa2; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp2; 
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }                                       
  else if (IS_FF(fa3)){
    for (pfli=pflink; pfli->nbface!=fa3; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp3;
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
  }
  else {
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }
   
  pnn = -ann/20.0*FMULT;
  p1n = -an1/20.0*FMULT;
  p2n = -an2/20.0*FMULT;
  p3n = -an3/20.0*FMULT;
  
  putpaij(fa0,n0,pnn,nn0,Z,f,u);
  putpaij(fa0,n1,p1n,nn0,Z,f,u);
  putpaij(fa0,n2,p2n,nn0,Z,f,u);
  putpaij(fa0,n3,p3n,nn0,Z,f,u);  
 }
}

#else

void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,rdetB,fmult2)
NODE *n0, *n1, *n2, *n3; FACE *fa0, *fa1, *fa2, *fa3; INT Z, f, u; FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{  eprintf("Error: aijb not available.\n");  }

#endif

#endif

#if !(DATA_STR & REDUCED) && (U_SPACE == P1C_NEW_FBUB)

void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,rdetB,fmult2)
NODE *n0, *n1, *n2, *n3; FACE *fa0, *fa1, *fa2, *fa3; INT Z, f, u; FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{  eprintf("Error: aijb not available.\n");  }

#endif

#if (DATA_STR & REDUCED) && !(U_SPACE == P1C_NEW_FBUB)

/* boundary condition has to be defined be means of pw. linear functions only */

void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                   rdetB,fmult2)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3;
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{
   FLINK *pflink, *pfli;
   FLOAT ann, an1, an2, an3, a11, a12, a13, a22, a23, a33, ppn, pp1, pp2, pp3;
   
   ann = DOT(b0,b0)*detB;
   an1 = DOT(b0,b1)*detB;
   an2 = DOT(b0,b2)*detB;
   an3 = DOT(b0,b3)*detB;
  
   if (IS_FN(n0))  
      nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB);
                                          
   if (IS_FF(fa0)){ 
      
      RHS_FOR_FACES 
      
      a11 = DOT(b1,b1)*detB;
      a12 = DOT(b1,b2)*detB;
      a13 = DOT(b1,b3)*detB;
      a22 = DOT(b2,b2)*detB;
      a23 = DOT(b2,b3)*detB;
      a33 = DOT(b3,b3)*detB;
     
      ppn = ( a11 + a22 + a33 + a12 + a13 + a23 )/210.0*fmult2;
      pp1 = ( an1 + an1 - a23 )*DOT(nn0,nn1)/420.0*fmult2;
      pp2 = ( an2 + an2 - a13 )*DOT(nn0,nn2)/420.0*fmult2; 
      pp3 = ( an3 + an3 - a12 )*DOT(nn0,nn3)/420.0*fmult2; 
    
      COEFF_FF(fa0,Z) += ppn;
     
      pflink = fa0->fstart;
      if (IS_FF(fa1) && IS_FF(fa2) && IS_FF(fa3)){
         for (pfli=pflink; pfli->nbface!=fa1 && pfli->nbface!=fa2 
                    && pfli->nbface!=fa3 ; pfli=pfli->next);
         if (pfli->nbface==fa1){
            COEFF_FL(pfli,Z) += pp1;
            putppij(pfli->next,fa2,fa3,pp2,pp3,Z);
         }            
         else if(pfli->nbface==fa2){
            COEFF_FL(pfli,Z) += pp2;
            putppij(pfli->next,fa1,fa3,pp1,pp3,Z);
         }
         else{    
            COEFF_FL(pfli,Z) += pp3;
            putppij(pfli->next,fa1,fa2,pp1,pp2,Z);
         }
      }
      else if (IS_FF(fa1) && IS_FF(fa2))
         putppij(pflink,fa1,fa2,pp1,pp2,Z);
      else if (IS_FF(fa1) && IS_FF(fa3))
         putppij(pflink,fa1,fa3,pp1,pp3,Z);
      else if (IS_FF(fa2) && IS_FF(fa3))
         putppij(pflink,fa2,fa3,pp2,pp3,Z);
      else if (IS_FF(fa1)){
         for (pfli=pflink; pfli->nbface!=fa1; pfli=pfli->next);
         COEFF_FL(pfli,Z) += pp1;
      }                                         
      else if (IS_FF(fa2)){
         for (pfli=pflink; pfli->nbface!=fa2; pfli=pfli->next);
         COEFF_FL(pfli,Z) += pp2; 
      }                                       
      else if (IS_FF(fa3)){
         for (pfli=pflink; pfli->nbface!=fa3; pfli=pfli->next);
         COEFF_FL(pfli,Z) += pp3;
      }
   }
}

#endif

#if (DATA_STR & REDUCED) && (U_SPACE == P1C_NEW_FBUB)

/* boundary condition has to be defined be means of pw. linear functions only */

void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                   rdetB,fmult2)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3;
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{
   FLOAT ann, an1, an2, an3, qnn;
   
   qnn = DOT(b0,b0);
   ann = qnn*detB;
  
   if (IS_FN(n0)){  
      an1 = DOT(b0,b1)*detB;
      an2 = DOT(b0,b2)*detB;
      an3 = DOT(b0,b3)*detB;
      nn_matrix(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,rdetB);
   }
                                          
   if (IS_FF(fa0)){ 
      
      RHS_FOR_FACESN 
/*      FD(fa0,f) += integr3p(n0,n1,n2,n3,nn0,rdetB)*FMULT;*/
      
      COEFF_FF(fa0,Z) += (DOT(b1,b1) + DOT(b2,b2) + DOT(b3,b3) + 
                                       21.*qnn)/1680.0*detB*fmult2;
   }
}

#endif

#endif  /* !(DATA_STR & KORN_MATRIX) */
  
#if (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & SCALAR_FACE_DATA)

#if !(DATA_STR & LG_DATA)

 void set_Bi(pelement,pnode,n1,fan,nn,arn,i,u,br,b0,detB)
 ELEMENT *pelement;
 NODE *pnode, *n1;
 FACE *fan; 
 INT br, u, i;
 FLOAT nn[DIM], arn, b0[DIM2], detB;
 {
    test_normal_vector_direction(nn,pnode->myvertex->x,n1->myvertex->x,&arn);
    SET2(COEFF_BNP(pelement,0,i),b0,(-detB))
    COEFF_BF(pelement,0,i) = -arn/60.0*FMULT;
    if (IS_DN(pnode))
       ED(pelement,br) += detB*DOT(b0,NDD(pnode,u));
    if (NOT_FF(fan))
       ED(pelement,br) += FD(fan,u)*arn/60.0*FMULT;
 }

#else

 void set_Bi(pelement,pnode,n1,fan,nn,arn,i,u,br,b0,detB)
 ELEMENT *pelement;
 NODE *pnode, *n1;
 FACE *fan; 
 INT br, u, i;
 FLOAT nn[DIM], arn, b0[DIM2], detB;
 {
    test_normal_vector_direction(nn,pnode->myvertex->x,n1->myvertex->x,&arn);
    if (IS_FN(pnode))
       SET2(COEFF_BNP(pelement,0,i),b0,(-detB))
    else if (pnode->lgd){
       COEFF_BN(pelement,0,i,0) = -detB*DOT(b0,pnode->lgd->t1);
       COEFF_BN(pelement,0,i,1) = -detB*DOT(b0,pnode->lgd->t2);
    }
    COEFF_BF(pelement,0,i) = -arn/60.0*FMULT;
    if (IS_DN(pnode))
       ED(pelement,br) += detB*DOT(b0,NDD(pnode,u));
    if (NOT_FF(fan))
       ED(pelement,br) += FD(fan,u)*arn/60.0*FMULT;
 }

#endif
 
void set_B(pelement,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,
                                          ar0,ar1,ar2,ar3,u,br,b0,b1,b2,b3,detB)
ELEMENT *pelement;
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3; 
INT br, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], ar0, ar1, ar2, ar3, 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   ED(pelement,br) = 0.0;
   set_Bi(pelement,n0,n1,fa0,nn0,ar0,0,u,br,b0,detB);
   set_Bi(pelement,n1,n2,fa1,nn1,ar1,1,u,br,b1,detB);
   set_Bi(pelement,n2,n3,fa2,nn2,ar2,2,u,br,b2,detB);
   set_Bi(pelement,n3,n0,fa3,nn3,ar3,3,u,br,b3,detB);
}

#else

void set_B(pelement,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,
                                          ar0,ar1,ar2,ar3,u,br,b0,b1,b2,b3,detB)
ELEMENT *pelement; NODE *n0, *n1, *n2, *n3; FACE *fa0, *fa1, *fa2, *fa3; 
INT br, u; FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], ar0, ar1, ar2, ar3, 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{  eprintf("Error: set_B not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (E_DATA & ExDN_MATR) && (N_DATA & SCALAR_NODE_DATA)

void set_bijb(n0,n1,Z,br,u,b1,tdetB)
NODE *n0, *n1;
INT Z, br, u;
FLOAT b1[DIM2], tdetB;
{
   LINK *pli;

   if (IS_FN(n1)){
      for (pli=n0->start; pli->nbnode != n1; pli=pli->next);
      SET9(COEFFBLP(pli,Z),b1,tdetB)
   }
   else if (IS_DN(n1))
      NDS(n0,br) += DOT(b1,NDD(n1,u))*tdetB;
}

void bijb(pelem,i,n0,n1,n2,n3,Z,br,u,b0,b1,b2,b3,tdetB) /* n0 = pelem->n[i] */
ELEMENT *pelem;                                         /* tdetB = |K|/4    */
NODE *n0, *n1, *n2, *n3;
INT i, Z, br, u;
FLOAT b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], tdetB;
{
   if (IS_DN(n0))
      NDS(n0,br) += DOT(b0,NDD(n0,u))*tdetB;
   set_bijb(n0,n1,Z,br,u,b1,tdetB);
   set_bijb(n0,n2,Z,br,u,b2,tdetB);
   set_bijb(n0,n3,Z,br,u,b3,tdetB);
   SET2(COEFF_BNP(pelem,0,i),b0,FMULT*tdetB/210.)
}

#else

void bijb(pelem,i,n0,n1,n2,n3,Z,br,u,b0,b1,b2,b3,tdetB)
ELEMENT *pelem; NODE *n0, *n1, *n2, *n3; INT i, Z, br, u; FLOAT b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], tdetB;
{  eprintf("Error: bijb not available.\n");  }

#endif

/******************************************************************************/

#if !(U_SPACE == P1C_NEW_FBUB) && (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA)

 void ff_copy(mg,fan,f,pp,Z,ff,pu)
 MULTIGRID *mg;
 FACE *fan, *f;
 FLOAT pp;
 INT Z, ff, pu;
 {
    GRID *pg;
    FLINK *pfl1, *pfl2;
    
    if (IS_FF(f)){
       for (pfl1 = fan->fstart; pfl1->nbface != f; pfl1 = pfl1->next);
       pg = TOP_GRID(mg)->coarser;
       while (fan->father->index < pg->minFaceIndex) pg = pg->coarser;
       while (f->index > pg->maxFaceIndex) f = f->father;
       for (pfl2 = fan->father->fstart; pfl2->nbface != f; pfl2 = pfl2->next);
       COEFF_FL(pfl1,Z) = COEFF_FL(pfl2,Z);
    }
    else
       FD(fan,ff) -= FD(f,pu)*pp;
 }

#else

void ff_copy(mg,fan,f,pp,Z,ff,pu)
MULTIGRID *mg; FACE *fan, *f; FLOAT pp; INT Z, ff, pu;
{  eprintf("Error: ff_copy not available.\n");  }

#endif

#if !(DATA_STR & REDUCED) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)

 void fn_copy(mg,f,n,Z)
 MULTIGRID *mg;
 FACE *f;
 NODE *n;
 INT Z;
 {
    GRID *pg;
    FNLINK *pfn1, *pfn2;
    
    for (pfn1 = f->fnstart; pfn1->nbnode != n; pfn1 = pfn1->next);
    pg = TOP_GRID(mg)->coarser;
    while (f->father->index < pg->minFaceIndex) pg = pg->coarser;
    while (n->index > pg->maxNodeIndex) n = n->father;
    for (pfn2 = f->father->fnstart; pfn2->nbnode != n; pfn2 = pfn2->next);
    COEFF_FN(pfn1,Z,0) = COEFF_FN(pfn2,Z,0);
    COEFF_FN(pfn1,Z,1) = COEFF_FN(pfn2,Z,1);
    COEFF_FN(pfn1,Z,2) = COEFF_FN(pfn2,Z,2);
 }
 
 void nf_copy(mg,n,f,Z) 
 MULTIGRID *mg;
 NODE *n; 
 FACE *f;
 INT Z;
 {
   GRID *pg;
   NFLINK *pnf1, *pnf2;
   
   for (pnf1 = n->nfstart; pnf1->nbface != f; pnf1 = pnf1->next);
   pg = TOP_GRID(mg)->coarser;
   while (n->father->index < pg->minNodeIndex) pg = pg->coarser;
   while (f->index > pg->maxFaceIndex) f = f->father;
   for (pnf2 = n->father->nfstart; pnf2->nbface != f; pnf2 = pnf2->next);
   COEFF_NF(pnf1,Z,0) = COEFF_NF(pnf2,Z,0);
   COEFF_NF(pnf1,Z,1) = COEFF_NF(pnf2,Z,1);
   COEFF_NF(pnf1,Z,2) = COEFF_NF(pnf2,Z,2);
 }

#endif

#if !(DATA_STR & KORN_MATRIX)

 void nn_copy(mg,pnode,n,pli1,Z)
 MULTIGRID *mg;
 NODE *pnode, *n;
 LINK *pli1;
 INT Z;
 {
   GRID *pg;
   LINK *pli2;
   
   pg = TOP_GRID(mg)->coarser;
   while (pnode->father->index < pg->minNodeIndex) pg = pg->coarser;
   while (n->index > pg->maxNodeIndex) n = n->father;
   for (pli2 = pnode->father->start; pli2->nbnode != n; pli2 = pli2->next); 
   COEFFL(pli1,Z) = COEFFL(pli2,Z);
 }

 void gputaij(mg,pnode,n1,an1,Z,f,u)
 MULTIGRID *mg;
 NODE *pnode, *n1;
 INT Z, f, u;
 FLOAT an1;
 {
   LINK *pli;
   
   if (NOT_FN(n1))
      u0_to_f(pnode,f,n1,u,an1);
   else
      if (ON_FINEST_GRID(n1,mg)){ 
         for (pli = pnode->start; pli->nbnode != n1; pli = pli->next);
         if (!IS_EDGE_TO_OLD(pli))
            COEFFL(pli,Z) += an1;
      }
 }

#if !(DATA_STR & REDUCED) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) 

 void gputapij(mg,pnode,fa1,pn1,nn1,Z)
 MULTIGRID *mg;
 NODE *pnode; 
 FACE *fa1;
 INT Z;
 FLOAT pn1, nn1[DIM];
 {
   NFLINK *pnf;
   
   if (F_ON_FINEST_GRID(fa1,mg)){
      for (pnf = pnode->nfstart; pnf->nbface != fa1; pnf = pnf->next);
      COEFF_NF(pnf,Z,0) += pn1*nn1[0];
      COEFF_NF(pnf,Z,1) += pn1*nn1[1];
      COEFF_NF(pnf,Z,2) += pn1*nn1[2];
   }
   else
      nf_copy(mg,pnode,fa1,Z);
 }

 INT gaijb(mg,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,
                                                              detB,rdetB,fmult2)
 MULTIGRID *mg;
 NODE *n0, *n1, *n2, *n3;
 FACE *fa0, *fa1, *fa2, *fa3; 
 INT Z, f, u;
 FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
       b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
 {
  FLOAT ann, an1, an2, an3, a11, a12, a13, a22, a23, a33, 
        pnn,ppn, pp1, pp2, pp3;
  
  if ( (NOT_FN(n0)  || !ON_FINEST_GRID(n0,mg)) &&
       (NOT_FF(fa0) || !F_ON_FINEST_GRID(fa0,mg)) ) return(0);
       
  ann = DOT(b0,b0)*detB;
  an1 = DOT(b0,b1)*detB;
  an2 = DOT(b0,b2)*detB;
  an3 = DOT(b0,b3)*detB;
  
 if (IS_FN(n0) && ON_FINEST_GRID(n0,mg)){  
  
  ND(n0,f,0)+=integr1(n0->myvertex->x,n1->myvertex->x,
                      n2->myvertex->x,n3->myvertex->x,f01,rdetB);
  ND(n0,f,1)+=integr1(n0->myvertex->x,n1->myvertex->x,
                      n2->myvertex->x,n3->myvertex->x,f02,rdetB);
  ND(n0,f,2)+=integr1(n0->myvertex->x,n1->myvertex->x,
                      n2->myvertex->x,n3->myvertex->x,f03,rdetB);
  COEFFN(n0,Z) += ann;
  
  gputaij(mg,n0,n1,an1,Z,f,u);
  gputaij(mg,n0,n2,an2,Z,f,u);
  gputaij(mg,n0,n3,an3,Z,f,u);
  
  gputapij(mg,n0,fa1,-an1/20.0*FMULT,nn1,Z);
  gputapij(mg,n0,fa2,-an2/20.0*FMULT,nn2,Z);
  gputapij(mg,n0,fa3,-an3/20.0*FMULT,nn3,Z);
  
  pnn = -ann/20.0*FMULT;
  
  if (NOT_FF(fa0)){
     ND(n0,f,0) -= FD(fa0,u)*pnn*nn0[0];
     ND(n0,f,1) -= FD(fa0,u)*pnn*nn0[1];
     ND(n0,f,2) -= FD(fa0,u)*pnn*nn0[2];
  }
  else 
     nf_copy(mg,n0,fa0,Z);
 }
                                          
 if (IS_FF(fa0) && F_ON_FINEST_GRID(fa0,mg)){ 
 
  RHS_FOR_FACES 
  
  a11 = DOT(b1,b1)*detB;
  a12 = DOT(b1,b2)*detB;
  a13 = DOT(b1,b3)*detB;
  a22 = DOT(b2,b2)*detB;
  a23 = DOT(b2,b3)*detB;
  a33 = DOT(b3,b3)*detB;
  
  ppn = ( a11 + a22 + a33 + a12 + a13 + a23 )/210.0*fmult2;
  pp1 = ( an1 + an1 - a23 )*DOT(nn0,nn1)/420.0*fmult2;
  pp2 = ( an2 + an2 - a13 )*DOT(nn0,nn2)/420.0*fmult2; 
  pp3 = ( an3 + an3 - a12 )*DOT(nn0,nn3)/420.0*fmult2; 
  
  COEFF_FF(fa0,Z) += ppn;
  
  ff_copy(mg,fa0,fa1,pp1,Z,f,u);
  ff_copy(mg,fa0,fa2,pp2,Z,f,u);
  ff_copy(mg,fa0,fa3,pp3,Z,f,u);
   
  if (IS_FN(n0)) 
     fn_copy(mg,fa0,n0,Z);
  else
     putpaij(fa0,n0,-ann/20.0*FMULT,nn0,Z,f,u); 
  putpaij(fa0,n1,-an1/20.0*FMULT,nn0,Z,f,u);
  putpaij(fa0,n2,-an2/20.0*FMULT,nn0,Z,f,u);
  putpaij(fa0,n3,-an3/20.0*FMULT,nn0,Z,f,u);  
 }
 return(0);
 }

#else

INT gaijb(mg,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,rdetB,fmult2)
MULTIGRID *mg; NODE *n0, *n1, *n2, *n3; FACE *fa0, *fa1, *fa2, *fa3; INT Z, f, u; FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
{  eprintf("gaijb not available.\n"); return(0);  }

#endif
  
#if DATA_STR & LG_DATA

void LGij(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                          rdetB)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3; 
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB;
{
   eprintf("LGij not available.\n");
}

#endif

#endif

/******************************************************************************/
/*                                                                            */
/*                             stiffness matrix 2                             */
/*                                                                            */
/******************************************************************************/

/*  stiffness matrix corresponding to 
         0.5*nu*\int_\om (\nabla u+\nabla u^T)\cdot(\nabla v+\nabla v^T) \dx  */

#if DATA_STR & KORN_MATRIX
 
 void fill(pli,an1,Z,b0,b1,detB)
 LINK *pli;
 INT Z;
 FLOAT an1, b0[DIM2], b1[DIM2], detB;
 {
    COEFFLL(pli,Z,0,0) += b0[0]*b1[0]*detB+an1;
    COEFFLL(pli,Z,1,1) += b0[1]*b1[1]*detB+an1;
    COEFFLL(pli,Z,2,2) += b0[2]*b1[2]*detB+an1;
    COEFFLL(pli,Z,0,1) += b0[1]*b1[0]*detB;
    COEFFLL(pli,Z,1,0) += b0[0]*b1[1]*detB;
    COEFFLL(pli,Z,0,2) += b0[2]*b1[0]*detB;
    COEFFLL(pli,Z,2,0) += b0[0]*b1[2]*detB;
    COEFFLL(pli,Z,1,2) += b0[2]*b1[1]*detB;
    COEFFLL(pli,Z,2,1) += b0[1]*b1[2]*detB;    
 } 
 
 void Putaij(plin,n1,n2,an1,an2,Z,b0,b1,b2,detB)
 LINK *plin;
 NODE *n1, *n2;
 FLOAT an1,an2, b0[DIM2], b1[DIM2], b2[DIM2], detB;
 INT Z;
 {
   LINK *pli;
   
   for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2; pli=pli->next);
   if (pli->nbnode == n1){
      fill(pli,an1,Z,b0,b1,detB);
      for (pli=pli->next; pli->nbnode != n2; pli=pli->next); 
      fill(pli,an2,Z,b0,b2,detB);
   }
   else {
      fill(pli,an2,Z,b0,b2,detB);
      for (pli=pli->next; pli->nbnode != n1; pli=pli->next); 
      fill(pli,an1,Z,b0,b1,detB);
   }
 }

 void putapij(pnfl,fa1,fa2,pn1,pn2,Z,nn1,nn2,dn1,dn2,b1,b2)
 NFLINK *pnfl; 
 FACE *fa1, *fa2;
 INT Z;
 FLOAT pn1, pn2, nn1[DIM], nn2[DIM], dn1, dn2, b1[DIM2], b2[DIM2];
 {
   NFLINK *pnf;
   
   for (pnf = pnfl; pnf->nbface != fa1 && pnf->nbface != fa2; pnf = pnf->next);
   if (pnf->nbface==fa1){
      COEFF_NF(pnf,Z,0) += pn1*nn1[0] + dn1*b1[0];
      COEFF_NF(pnf,Z,1) += pn1*nn1[1] + dn1*b1[1];
      COEFF_NF(pnf,Z,2) += pn1*nn1[2] + dn1*b1[2];
      for (pnf = pnf->next; pnf->nbface != fa2; pnf = pnf->next); 
      COEFF_NF(pnf,Z,0) += pn2*nn2[0] + dn2*b2[0];
      COEFF_NF(pnf,Z,1) += pn2*nn2[1] + dn2*b2[1];
      COEFF_NF(pnf,Z,2) += pn2*nn2[2] + dn2*b2[2];
   }
   else {
      COEFF_NF(pnf,Z,0) += pn2*nn2[0] + dn2*b2[0];
      COEFF_NF(pnf,Z,1) += pn2*nn2[1] + dn2*b2[1];
      COEFF_NF(pnf,Z,2) += pn2*nn2[2] + dn2*b2[2];
      for (pnf = pnf->next; pnf->nbface != fa1; pnf = pnf->next); 
      COEFF_NF(pnf,Z,0) += pn1*nn1[0] + dn1*b1[0];
      COEFF_NF(pnf,Z,1) += pn1*nn1[1] + dn1*b1[1];
      COEFF_NF(pnf,Z,2) += pn1*nn1[2] + dn1*b1[2];
   }
 }

 void putpaij(fan,n1,p1n,d1n,nn,b0,Z,f,u)
 FACE *fan;
 NODE *n1;
 FLOAT p1n, d1n, nn[DIM], b0[DIM2];
 INT Z, f, u;
 {
    FNLINK *pfn;
    
    if (IS_FN(n1)){
       for (pfn=fan->fnstart; pfn->nbnode!=n1; pfn=pfn->next);
       COEFF_FN(pfn,Z,0) += p1n*nn[0] + d1n*b0[0];
       COEFF_FN(pfn,Z,1) += p1n*nn[1] + d1n*b0[1];
       COEFF_FN(pfn,Z,2) += p1n*nn[2] + d1n*b0[2];
    }
    else if (IS_DN(n1))
       FD(fan,f) -= (p1n*nn[0] + d1n*b0[0])*ND(n1,u,0) +
                    (p1n*nn[1] + d1n*b0[1])*ND(n1,u,1) +
                    (p1n*nn[2] + d1n*b0[2])*ND(n1,u,2);
 }
 
 void u0_to_f(pnode,f,n1,u,q,b0,b1,detB)
 NODE *pnode, *n1;
 INT f, u;
 FLOAT q, b0[DIM2], b1[DIM2], detB;
 {
   FLOAT c;
   
   if (IS_DN(n1)){
      c = detB*DOT(NDD(n1,u),b0);
      ND(pnode,f,0) -= ND(n1,u,0)*q + c*b1[0];
      ND(pnode,f,1) -= ND(n1,u,1)*q + c*b1[1];
      ND(pnode,f,2) -= ND(n1,u,2)*q + c*b1[2];
   }
 }

 void aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                   rdetB,fmult2)
 NODE *n0, *n1, *n2, *n3;
 FACE *fa0, *fa1, *fa2, *fa3; 
 INT Z, f, u;
 FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
       b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
 {
  LINK *plin, *pli;
  NFLINK *pnf;
  FLINK *pflink, *pfli;
  FLOAT ann, an1, an2, an3, a11, a12, a13, a22, a23, a33, 
        pnn, pn1, pn2, pn3, p1n, p2n, p3n, ppn, pp1, pp2, pp3, 
        cnn, cn1, cn2, cn3, c1n, c11, c12, c13, 
        c2n, c21, c22, c23, c3n, c31, c32, c33,
        dnn, dn1, dn2, dn3, d1n, d2n, d3n, c;
         
  ann = DOT(b0,b0)*detB;
  an1 = DOT(b0,b1)*detB;
  an2 = DOT(b0,b2)*detB;
  an3 = DOT(b0,b3)*detB;
  cnn = DOT(b0,nn0);
  cn1 = DOT(b0,nn1);
  cn2 = DOT(b0,nn2);
  cn3 = DOT(b0,nn3);
  
 if (IS_FN(n0)){  
  
  ND(n0,f,0)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f01,rdetB);
  ND(n0,f,1)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f02,rdetB);
  ND(n0,f,2)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f03,rdetB);
  COEFFNN(n0,Z,0,0) += b0[0]*b0[0]*detB + ann;
  COEFFNN(n0,Z,1,1) += b0[1]*b0[1]*detB + ann;
  COEFFNN(n0,Z,2,2) += b0[2]*b0[2]*detB + ann;
  COEFFNN(n0,Z,0,1) += (c=b0[0]*b0[1]*detB);
  COEFFNN(n0,Z,1,0) += c;
  COEFFNN(n0,Z,0,2) += (c=b0[0]*b0[2]*detB);
  COEFFNN(n0,Z,2,0) += c;
  COEFFNN(n0,Z,1,2) += (c=b0[1]*b0[2]*detB);
  COEFFNN(n0,Z,2,1) += c;
  
  plin=n0->start;
  if (IS_FN(n1) && IS_FN(n2) && IS_FN(n3)){
    for (pli=plin; pli->nbnode != n1 && pli->nbnode != n2 
                && pli->nbnode != n3; pli=pli->next);
    if (pli->nbnode == n1){
       fill(pli,an1,Z,b0,b1,detB);
       Putaij(pli->next,n2,n3,an2,an3,Z,b0,b2,b3,detB);
    }            
    else if(pli->nbnode == n2){
       fill(pli,an2,Z,b0,b2,detB);
       Putaij(pli->next,n1,n3,an1,an3,Z,b0,b1,b3,detB);
    }
    else {    
       fill(pli,an3,Z,b0,b3,detB);
       Putaij(pli->next,n1,n2,an1,an2,Z,b0,b1,b2,detB);
    }
  }
  else if (IS_FN(n1) && IS_FN(n2)){
    Putaij(plin,n1,n2,an1,an2,Z,b0,b1,b2,detB);
    u0_to_f(n0,f,n3,u,an3,b0,b3,detB);
  }
  else if (IS_FN(n1) && IS_FN(n3)){
    Putaij(plin,n1,n3,an1,an3,Z,b0,b1,b3,detB);
    u0_to_f(n0,f,n2,u,an2,b0,b2,detB);
  }
  else if (IS_FN(n2) && IS_FN(n3)){ 
    Putaij(plin,n2,n3,an2,an3,Z,b0,b2,b3,detB);
    u0_to_f(n0,f,n1,u,an1,b0,b1,detB);
  }
  else if (IS_FN(n1)){
    for (pli=plin; pli->nbnode != n1; pli=pli->next);
    fill(pli,an1,Z,b0,b1,detB);
    u0_to_f(n0,f,n2,u,an2,b0,b2,detB);
    u0_to_f(n0,f,n3,u,an3,b0,b3,detB);
  }                                         
  else if (IS_FN(n2)){
    for (pli=plin; pli->nbnode != n2; pli=pli->next);
    fill(pli,an2,Z,b0,b2,detB); 
    u0_to_f(n0,f,n1,u,an1,b0,b1,detB);
    u0_to_f(n0,f,n3,u,an3,b0,b3,detB);                                    
  }                                       
  else if (IS_FN(n3)){
    for (pli=plin; pli->nbnode != n3; pli=pli->next);
    fill(pli,an3,Z,b0,b3,detB);
    u0_to_f(n0,f,n1,u,an1,b0,b1,detB);
    u0_to_f(n0,f,n2,u,an2,b0,b2,detB);  
  }
  else {
      u0_to_f(n0,f,n1,u,an1,b0,b1,detB);
      u0_to_f(n0,f,n2,u,an2,b0,b2,detB);
      u0_to_f(n0,f,n3,u,an3,b0,b3,detB);      
  }
  
  pnn = -ann/20.0*FMULT;
  pn1 = -an1/20.0*FMULT;
  pn2 = -an2/20.0*FMULT;
  pn3 = -an3/20.0*FMULT;
  dnn = -cnn*detB/20.0*FMULT;
  dn1 = -cn1*detB/20.0*FMULT;
  dn2 = -cn2*detB/20.0*FMULT;
  dn3 = -cn3*detB/20.0*FMULT;
  
  for (pnf=n0->nfstart; pnf->nbface!=fa1 && pnf->nbface!=fa2 
                        && pnf->nbface!=fa3; pnf=pnf->next);
  if (pnf->nbface==fa1){
    COEFF_NF(pnf,Z,0) += pn1*nn1[0] + dn1*b1[0];
    COEFF_NF(pnf,Z,1) += pn1*nn1[1] + dn1*b1[1];
    COEFF_NF(pnf,Z,2) += pn1*nn1[2] + dn1*b1[2];
    putapij(pnf->next,fa2,fa3,pn2,pn3,Z,nn2,nn3,dn2,dn3,b2,b3);
  }
  else if (pnf->nbface==fa2){
    COEFF_NF(pnf,Z,0) += pn2*nn2[0] + dn2*b2[0];
    COEFF_NF(pnf,Z,1) += pn2*nn2[1] + dn2*b2[1];
    COEFF_NF(pnf,Z,2) += pn2*nn2[2] + dn2*b2[2];
    putapij(pnf->next,fa1,fa3,pn1,pn3,Z,nn1,nn3,dn1,dn3,b1,b3);
  }
  else {
    COEFF_NF(pnf,Z,0) += pn3*nn3[0] + dn3*b3[0];
    COEFF_NF(pnf,Z,1) += pn3*nn3[1] + dn3*b3[1];
    COEFF_NF(pnf,Z,2) += pn3*nn3[2] + dn3*b3[2];
    putapij(pnf->next,fa1,fa2,pn1,pn2,Z,nn1,nn2,dn1,dn2,b1,b2);
  } 
  
  if (NOT_FF(fa0)){
    ND(n0,f,0) -= FD(fa0,u)*(pnn*nn0[0] + dnn*b0[0]);
    ND(n0,f,1) -= FD(fa0,u)*(pnn*nn0[1] + dnn*b0[1]);
    ND(n0,f,2) -= FD(fa0,u)*(pnn*nn0[2] + dnn*b0[2]);
  }
  else {
    for  (pnf=n0->nfstart; pnf->nbface!=fa0; pnf=pnf->next);
    COEFF_NF(pnf,Z,0) += pnn*nn0[0] + dnn*b0[0];
    COEFF_NF(pnf,Z,1) += pnn*nn0[1] + dnn*b0[1];
    COEFF_NF(pnf,Z,2) += pnn*nn0[2] + dnn*b0[2];  
  }
 }
                                          
 if (IS_FF(fa0)){ 
 
  RHS_FOR_FACES 
  
  a11 = DOT(b1,b1)*detB;
  a12 = DOT(b1,b2)*detB;
  a13 = DOT(b1,b3)*detB;
  a22 = DOT(b2,b2)*detB;
  a23 = DOT(b2,b3)*detB;
  a33 = DOT(b3,b3)*detB;
  
  c1n = DOT(b1,nn0);
  c11 = DOT(b1,nn1);
  c12 = DOT(b1,nn2);
  c13 = DOT(b1,nn3);
  c2n = DOT(b2,nn0);
  c21 = DOT(b2,nn1);
  c22 = DOT(b2,nn2);
  c23 = DOT(b2,nn3);
  c3n = -cnn - c1n - c2n;
  c31 = -cn1 - c11 - c21;  
  c32 = -cn2 - c12 - c22;  
  c33 = -cn3 - c13 - c23;   
  
  ppn = ( 2.0*(a11 + a22 + a33 + a12 + a13 + a23) + 
          (cnn*cnn + c1n*c1n + c2n*c2n + c3n*c3n)*detB )/420.0*fmult2;
  pp1 = ( (an1 + an1 - a23)*DOT(nn0,nn1) +
          (cn1*c1n + c11*cnn - (c21*c3n + c31*c2n)/2.0)*detB )/420.0*fmult2;
  pp2 = ( (an2 + an2 - a13)*DOT(nn0,nn2) +
          (cn2*c2n + c22*cnn - (c12*c3n + c32*c1n)/2.0)*detB )/420.0*fmult2; 
  pp3 = ( (an3 + an3 - a12)*DOT(nn0,nn3) +
          (cn3*c3n + c33*cnn - (c13*c2n + c23*c1n)/2.0)*detB )/420.0*fmult2; 
  
  COEFF_FF(fa0,Z) += ppn;
  
  pflink = fa0->fstart;
  if (IS_FF(fa1) && IS_FF(fa2) && IS_FF(fa3)){
    for (pfli=pflink; pfli->nbface!=fa1 && pfli->nbface!=fa2 
                && pfli->nbface!=fa3 ; pfli=pfli->next);
    if (pfli->nbface==fa1){
       COEFF_FL(pfli,Z) += pp1;
       putppij(pfli->next,fa2,fa3,pp2,pp3,Z);
    }            
    else if(pfli->nbface==fa2){
       COEFF_FL(pfli,Z) += pp2;
       putppij(pfli->next,fa1,fa3,pp1,pp3,Z);
    }
    else {    
       COEFF_FL(pfli,Z) += pp3;
       putppij(pfli->next,fa1,fa2,pp1,pp2,Z);
    }
  }
  else if (IS_FF(fa1) && IS_FF(fa2)){
    putppij(pflink,fa1,fa2,pp1,pp2,Z);
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }
  else if (IS_FF(fa1) && IS_FF(fa3)){
    putppij(pflink,fa1,fa3,pp1,pp3,Z);
    FD(fa0,f) -= FD(fa2,u)*pp2;
  }
  else if (IS_FF(fa2) && IS_FF(fa3)){ 
    putppij(pflink,fa2,fa3,pp2,pp3,Z);
    FD(fa0,f) -= FD(fa1,u)*pp1;
  }
  else if (IS_FF(fa1)){
    for (pfli=pflink; pfli->nbface!=fa1; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }                                         
  else if (IS_FF(fa2)){
    for (pfli=pflink; pfli->nbface!=fa2; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp2; 
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }                                       
  else if (IS_FF(fa3)){
    for (pfli=pflink; pfli->nbface!=fa3; pfli=pfli->next);
    COEFF_FL(pfli,Z) += pp3;
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
  }
  else {
    FD(fa0,f) -= FD(fa1,u)*pp1;
    FD(fa0,f) -= FD(fa2,u)*pp2;
    FD(fa0,f) -= FD(fa3,u)*pp3;
  }
   
  pnn = -ann/20.0*FMULT;
  p1n = -an1/20.0*FMULT;
  p2n = -an2/20.0*FMULT;
  p3n = -an3/20.0*FMULT;
  dnn = -cnn*detB/20.0*FMULT;
  d1n = -c1n*detB/20.0*FMULT;
  d2n = -c2n*detB/20.0*FMULT;
  d3n = -c3n*detB/20.0*FMULT;
  
  putpaij(fa0,n0,pnn,dnn,nn0,b0,Z,f,u);
  putpaij(fa0,n1,p1n,d1n,nn0,b0,Z,f,u);
  putpaij(fa0,n2,p2n,d2n,nn0,b0,Z,f,u);
  putpaij(fa0,n3,p3n,d3n,nn0,b0,Z,f,u);  
 }
 
 }

/******************************************************************************/

 void nn_copy(mg,pnode,n,pli1,Z)
 MULTIGRID *mg;
 NODE *pnode, *n;
 LINK *pli1;
 INT Z;
 {
   GRID *pg;
   LINK *pli2;
   
   pg = TOP_GRID(mg)->coarser;
   while (pnode->father->index < pg->minNodeIndex) pg = pg->coarser;
   while (n->index > pg->maxNodeIndex) n = n->father;
   for (pli2 = pnode->father->start; pli2->nbnode != n; pli2 = pli2->next); 
   COEFFLL(pli1,Z,0,0) = COEFFLL(pli2,Z,0,0);
   COEFFLL(pli1,Z,0,1) = COEFFLL(pli2,Z,0,1);
   COEFFLL(pli1,Z,0,2) = COEFFLL(pli2,Z,0,2);
   COEFFLL(pli1,Z,1,0) = COEFFLL(pli2,Z,1,0);
   COEFFLL(pli1,Z,1,1) = COEFFLL(pli2,Z,1,1);
   COEFFLL(pli1,Z,1,2) = COEFFLL(pli2,Z,1,2);
   COEFFLL(pli1,Z,2,0) = COEFFLL(pli2,Z,2,0);
   COEFFLL(pli1,Z,2,1) = COEFFLL(pli2,Z,2,1);
   COEFFLL(pli1,Z,2,2) = COEFFLL(pli2,Z,2,2);   
 } 
 
 void gputaij(mg,pnode,n1,an1,Z,f,u,b0,b1,detB)
 MULTIGRID *mg;
 NODE *pnode, *n1;
 INT Z, f, u;
 FLOAT an1, b0[DIM2], b1[DIM2], detB;
 {
   LINK *pli;
   
   if (NOT_FN(n1))
      u0_to_f(pnode,f,n1,u,an1,b0,b1,detB);   
   else
      if (ON_FINEST_GRID(n1,mg)){ 
         for (pli = pnode->start; pli->nbnode != n1; pli = pli->next);
         if (!IS_EDGE_TO_OLD(pli))
            fill(pli,an1,Z,b0,b1,detB);
      }
 }

 void gputapij(mg,pnode,fa1,pn1,nn1,dn1,b1,Z)
 MULTIGRID *mg;
 NODE *pnode; 
 FACE *fa1;
 INT Z;
 FLOAT pn1, nn1[DIM], dn1, b1[DIM2];
 {
   NFLINK *pnf;
   
   if (F_ON_FINEST_GRID(fa1,mg)){
      for (pnf = pnode->nfstart; pnf->nbface != fa1; pnf = pnf->next);
      COEFF_NF(pnf,Z,0) += pn1*nn1[0] + dn1*b1[0];
      COEFF_NF(pnf,Z,1) += pn1*nn1[1] + dn1*b1[1];
      COEFF_NF(pnf,Z,2) += pn1*nn1[2] + dn1*b1[2];
   }
   else
      nf_copy(mg,pnode,fa1,Z);
 }

 INT gaijb(mg,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,
                                                              detB,rdetB,fmult2)
 MULTIGRID *mg;
 NODE *n0, *n1, *n2, *n3;
 FACE *fa0, *fa1, *fa2, *fa3; 
 INT Z, f, u;
 FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
       b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, fmult2;
 {
  FLOAT ann, an1, an2, an3, a11, a12, a13, a22, a23, a33, 
        pnn, ppn, pp1, pp2, pp3, 
        cnn, cn1, cn2, cn3, c1n, c11, c12, c13, 
        c2n, c21, c22, c23, c3n, c31, c32, c33, dnn, c;
        
  if ( (NOT_FN(n0)  || !ON_FINEST_GRID(n0,mg)) &&
       (NOT_FF(fa0) || !F_ON_FINEST_GRID(fa0,mg)) ) return(0);
       
  ann = DOT(b0,b0)*detB;
  an1 = DOT(b0,b1)*detB;
  an2 = DOT(b0,b2)*detB;
  an3 = DOT(b0,b3)*detB;
  cnn = DOT(b0,nn0);
  cn1 = DOT(b0,nn1);
  cn2 = DOT(b0,nn2);
  cn3 = DOT(b0,nn3);
  
 if (IS_FN(n0) && ON_FINEST_GRID(n0,mg)){  
  
  ND(n0,f,0)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f01,rdetB);
  ND(n0,f,1)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f02,rdetB);
  ND(n0,f,2)+=integr1(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,
                                                      n3->myvertex->x,f03,rdetB);
  COEFFNN(n0,Z,0,0) += b0[0]*b0[0]*detB + ann;
  COEFFNN(n0,Z,1,1) += b0[1]*b0[1]*detB + ann;
  COEFFNN(n0,Z,2,2) += b0[2]*b0[2]*detB + ann;
  COEFFNN(n0,Z,0,1) += (c=b0[0]*b0[1]*detB);
  COEFFNN(n0,Z,1,0) += c;
  COEFFNN(n0,Z,0,2) += (c=b0[0]*b0[2]*detB);
  COEFFNN(n0,Z,2,0) += c;
  COEFFNN(n0,Z,1,2) += (c=b0[1]*b0[2]*detB);
  COEFFNN(n0,Z,2,1) += c;
  
  gputaij(mg,n0,n1,an1,Z,f,u,b0,b1,detB);
  gputaij(mg,n0,n2,an2,Z,f,u,b0,b2,detB);
  gputaij(mg,n0,n3,an3,Z,f,u,b0,b3,detB);
  
  gputapij(mg,n0,fa1,-an1/20.0*FMULT,nn1,-cn1*detB/20.0*FMULT,b1,Z);
  gputapij(mg,n0,fa2,-an2/20.0*FMULT,nn2,-cn2*detB/20.0*FMULT,b2,Z);
  gputapij(mg,n0,fa3,-an3/20.0*FMULT,nn3,-cn3*detB/20.0*FMULT,b3,Z);
  
  pnn = -ann/20.0*FMULT;
  dnn = -cnn*detB/20.0*FMULT;
  
  if (NOT_FF(fa0)){
    ND(n0,f,0) -= FD(fa0,u)*(pnn*nn0[0] + dnn*b0[0]);
    ND(n0,f,1) -= FD(fa0,u)*(pnn*nn0[1] + dnn*b0[1]);
    ND(n0,f,2) -= FD(fa0,u)*(pnn*nn0[2] + dnn*b0[2]);
  }
  else 
     nf_copy(mg,n0,fa0,Z);
 }
                                          
 if (IS_FF(fa0) && F_ON_FINEST_GRID(fa0,mg)){ 
 
  RHS_FOR_FACES 
  
  a11 = DOT(b1,b1)*detB;
  a12 = DOT(b1,b2)*detB;
  a13 = DOT(b1,b3)*detB;
  a22 = DOT(b2,b2)*detB;
  a23 = DOT(b2,b3)*detB;
  a33 = DOT(b3,b3)*detB;
  
  c1n = DOT(b1,nn0);
  c11 = DOT(b1,nn1);
  c12 = DOT(b1,nn2);
  c13 = DOT(b1,nn3);
  c2n = DOT(b2,nn0);
  c21 = DOT(b2,nn1);
  c22 = DOT(b2,nn2);
  c23 = DOT(b2,nn3);
  c3n = -cnn - c1n - c2n;
  c31 = -cn1 - c11 - c21;  
  c32 = -cn2 - c12 - c22;  
  c33 = -cn3 - c13 - c23;   
  
  ppn = ( 2.0*(a11 + a22 + a33 + a12 + a13 + a23) + 
          (cnn*cnn + c1n*c1n + c2n*c2n + c3n*c3n)*detB )/420.0*fmult2;
  pp1 = ( (an1 + an1 - a23)*DOT(nn0,nn1) +
          (cn1*c1n + c11*cnn - (c21*c3n + c31*c2n)/2.0)*detB )/420.0*fmult2;
  pp2 = ( (an2 + an2 - a13)*DOT(nn0,nn2) +
          (cn2*c2n + c22*cnn - (c12*c3n + c32*c1n)/2.0)*detB )/420.0*fmult2; 
  pp3 = ( (an3 + an3 - a12)*DOT(nn0,nn3) +
          (cn3*c3n + c33*cnn - (c13*c2n + c23*c1n)/2.0)*detB )/420.0*fmult2; 
  
  COEFF_FF(fa0,Z) += ppn;
  
  ff_copy(mg,fa0,fa1,pp1,Z,f,u);
  ff_copy(mg,fa0,fa2,pp2,Z,f,u);
  ff_copy(mg,fa0,fa3,pp3,Z,f,u);
   
  if (IS_FN(n0)) 
     fn_copy(mg,fa0,n0,Z);
  else
     putpaij(fa0,n0,-ann/20.0*FMULT,-cnn*detB/20.0*FMULT,nn0,b0,Z,f,u); 
  putpaij(fa0,n1,-an1/20.0*FMULT,-c1n*detB/20.0*FMULT,nn0,b0,Z,f,u);
  putpaij(fa0,n2,-an2/20.0*FMULT,-c2n*detB/20.0*FMULT,nn0,b0,Z,f,u);
  putpaij(fa0,n3,-an3/20.0*FMULT,-c3n*detB/20.0*FMULT,nn0,b0,Z,f,u);  
 }
 return(0);
 }
   
#if DATA_STR & LG_DATA

void put_nlg(pnode,n1,b0,b1,Z,detB)
NODE *pnode, *n1;
FLOAT b0[DIM2], b1[DIM2], detB;
INT Z;
{
   NLGLINK *pnlg;
   FLOAT p0, p1, p2;
   
   if (n1->lgd){
      p0 = DOT(b0,b1);
      p1 = DOT(n1->lgd->t1,b0);
      p2 = DOT(n1->lgd->t2,b0);
      for (pnlg = pnode->lgstart; pnlg->nbnode != n1; pnlg = pnlg->next);
      COEFF_NLG(pnlg,Z,0,0) += ( n1->lgd->t1[0]*p0 + b1[0]*p1 )*detB;
      COEFF_NLG(pnlg,Z,1,0) += ( n1->lgd->t1[1]*p0 + b1[1]*p1 )*detB;
      COEFF_NLG(pnlg,Z,2,0) += ( n1->lgd->t1[2]*p0 + b1[2]*p1 )*detB;
      COEFF_NLG(pnlg,Z,0,1) += ( n1->lgd->t2[0]*p0 + b1[0]*p2 )*detB;
      COEFF_NLG(pnlg,Z,1,1) += ( n1->lgd->t2[1]*p0 + b1[1]*p2 )*detB;
      COEFF_NLG(pnlg,Z,2,1) += ( n1->lgd->t2[2]*p0 + b1[2]*p2 )*detB;
   }
}

void put_flg(fan,n1,b0,b1,nn,Z,detB)
FACE *fan;
NODE *n1;
FLOAT b0[DIM2], b1[DIM2], nn[DIM], detB;
INT Z;
{
   FLGLINK *pflg;
   FLOAT p0, q0;
   
   if (n1->lgd){
      p0 = DOT(b0,b1)*FMULT*detB/20.;
      q0 = DOT(b1,nn)*FMULT*detB/20.;
      for (pflg = fan->lgstart; pflg->nbnode != n1; pflg = pflg->next);
      COEFF_FLG(pflg,Z,0) -= p0*DOT(n1->lgd->t1,nn) + q0*DOT(n1->lgd->t1,b0);
      COEFF_FLG(pflg,Z,1) -= p0*DOT(n1->lgd->t2,nn) + q0*DOT(n1->lgd->t2,b0);
   }
}

void put_lglg(pnode,n1,b0,b1,Z,detB)
NODE *pnode, *n1;
FLOAT b0[DIM2], b1[DIM2], detB;
INT Z;
{
   LGLGLINK *plglg;
   FLOAT p11, p12, p21, p22, p0, p1, p2, q1, q2;
   
   if (n1->lgd){
     p11 = DOT(pnode->lgd->t1,n1->lgd->t1);
     p12 = DOT(pnode->lgd->t1,n1->lgd->t2);
     p21 = DOT(pnode->lgd->t2,n1->lgd->t1);
     p22 = DOT(pnode->lgd->t2,n1->lgd->t2);
     p0 = DOT(b0,b1);
     p1 = DOT(pnode->lgd->t1,b1);
     p2 = DOT(pnode->lgd->t2,b1);
     q1 = DOT(n1->lgd->t1,b0);
     q2 = DOT(n1->lgd->t2,b0);
     for (plglg = pnode->lgd->lgstart; plglg->nbnode != n1; plglg = plglg->next);
     COEFF_LGL(plglg,Z,0,0) += (p11*p0 + p1*q1)*detB;
     COEFF_LGL(plglg,Z,0,1) += (p12*p0 + p1*q2)*detB;
     COEFF_LGL(plglg,Z,1,0) += (p21*p0 + p2*q1)*detB;
     COEFF_LGL(plglg,Z,1,1) += (p22*p0 + p2*q2)*detB;
   }
}

void put_lgn(pnode,n1,b0,b1,Z,f,u,detB)
NODE *pnode, *n1;
FLOAT b0[DIM2], b1[DIM2], detB;
INT Z, f, u;
{
   LGNLINK *plgn;
   FLOAT p0, p1, p2;
   
   p0 = DOT(b0,b1);
   p1 = DOT(pnode->lgd->t1,b1);
   p2 = DOT(pnode->lgd->t2,b1);
   if (IS_FN(n1)){
      for (plgn = pnode->lgd->nstart; plgn->nbnode != n1; plgn = plgn->next);
      COEFF_LGN(plgn,Z,0,0) += (pnode->lgd->t1[0]*p0 + b0[0]*p1)*detB;
      COEFF_LGN(plgn,Z,0,1) += (pnode->lgd->t1[1]*p0 + b0[1]*p1)*detB;
      COEFF_LGN(plgn,Z,0,2) += (pnode->lgd->t1[2]*p0 + b0[2]*p1)*detB;
      COEFF_LGN(plgn,Z,1,0) += (pnode->lgd->t2[0]*p0 + b0[0]*p2)*detB;
      COEFF_LGN(plgn,Z,1,1) += (pnode->lgd->t2[1]*p0 + b0[1]*p2)*detB;
      COEFF_LGN(plgn,Z,1,2) += (pnode->lgd->t2[2]*p0 + b0[2]*p2)*detB;
   }
   else if (n1->lgd == NULL && IS_DN(n1)){
      NDLG(pnode,f,0) -= ( ND(n1,u,0)*(pnode->lgd->t1[0]*p0 + b0[0]*p1) + 
                           ND(n1,u,1)*(pnode->lgd->t1[1]*p0 + b0[1]*p1) + 
                           ND(n1,u,2)*(pnode->lgd->t1[2]*p0 + b0[2]*p1) )*detB;
      NDLG(pnode,f,1) -= ( ND(n1,u,0)*(pnode->lgd->t2[0]*p0 + b0[0]*p2) + 
                           ND(n1,u,1)*(pnode->lgd->t2[1]*p0 + b0[1]*p2) + 
                           ND(n1,u,2)*(pnode->lgd->t2[2]*p0 + b0[2]*p2) )*detB;
   }
}
 
INT put_lgf(pnode,fa1,b0,b1,nn1,Z,f,u,detB)
NODE *pnode;
FACE *fa1;
FLOAT b0[DIM2], b1[DIM2], nn1[DIM], detB;
INT Z, f, u;
{
   LGFLINK *plgf;
   FLOAT p0, p1, p2, q0, q1, q2;
   
   p0 = DOT(b0,b1);
   p1 = DOT(pnode->lgd->t1,nn1);
   p2 = DOT(pnode->lgd->t2,nn1);
   q0 = DOT(b0,nn1);
   q1 = DOT(pnode->lgd->t1,b1);
   q2 = DOT(pnode->lgd->t2,b1);
   if (fa1->type == 0){
      for (plgf = pnode->lgd->fstart; plgf->nbface != fa1; plgf = plgf->next);
      COEFF_LGF(plgf,Z,0) -= (p1*p0 + q1*q0)*FMULT*detB/20.;
      COEFF_LGF(plgf,Z,1) -= (p2*p0 + q2*q0)*FMULT*detB/20.;
   }
   else{
      NDLG(pnode,f,0) += FD(fa1,u)*(p1*p0 + q1*q0)*FMULT*detB/20.;
      NDLG(pnode,f,1) += FD(fa1,u)*(p2*p0 + q2*q0)*FMULT*detB/20.;
   }
   return(0);
}

void LGij(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b0,b1,b2,b3,detB,
                                                                          rdetB)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3; 
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB;
{
   FLOAT p1, p2, p3;
         
   if (IS_FN(n0)){
      put_nlg(n0,n1,b0,b1,Z,detB);
      put_nlg(n0,n2,b0,b2,Z,detB);
      put_nlg(n0,n3,b0,b3,Z,detB);
   }
   if (IS_FF(fa0)){
      put_flg(fa0,n0,b0,b0,nn0,Z,detB);
      put_flg(fa0,n1,b0,b1,nn0,Z,detB);
      put_flg(fa0,n2,b0,b2,nn0,Z,detB);   
      put_flg(fa0,n3,b0,b3,nn0,Z,detB);
   }
   if (n0->lgd){
      p1 = integr1(n0->myvertex->x,n1->myvertex->x,
                   n2->myvertex->x,n3->myvertex->x,f01,rdetB);
      p2 = integr1(n0->myvertex->x,n1->myvertex->x,
                   n2->myvertex->x,n3->myvertex->x,f02,rdetB);
      p3 = integr1(n0->myvertex->x,n1->myvertex->x,
                   n2->myvertex->x,n3->myvertex->x,f03,rdetB);
      NDLG(n0,f,0) += p1*n0->lgd->t1[0] + p2*n0->lgd->t1[1] + p3*n0->lgd->t1[2];
      NDLG(n0,f,1) += p1*n0->lgd->t2[0] + p2*n0->lgd->t2[1] + p3*n0->lgd->t2[2];
      p1 = DOT(n0->lgd->t1,b0);
      p2 = DOT(n0->lgd->t2,b0);
      p3 = DOT(b0,b0);
      COEFF_LG(n0,Z,0,0) += (p3 + p1*p1)*detB;
      COEFF_LG(n0,Z,0,1) += p1*p2*detB;
      COEFF_LG(n0,Z,1,0) += p1*p2*detB;
      COEFF_LG(n0,Z,1,1) += (p3 + p2*p2)*detB;
      put_lglg(n0,n1,b0,b1,Z,detB);
      put_lglg(n0,n2,b0,b2,Z,detB);
      put_lglg(n0,n3,b0,b3,Z,detB);
      put_lgn(n0,n1,b0,b1,Z,f,u,detB);
      put_lgn(n0,n2,b0,b2,Z,f,u,detB);
      put_lgn(n0,n3,b0,b3,Z,f,u,detB);
      put_lgf(n0,fa0,b0,b0,nn0,Z,f,u,detB);
      put_lgf(n0,fa1,b0,b1,nn1,Z,f,u,detB);
      put_lgf(n0,fa2,b0,b2,nn2,Z,f,u,detB);
      put_lgf(n0,fa3,b0,b3,nn3,Z,f,u,detB);
   }
}

#endif

#endif

/******************************************************************************/
/*                                                                            */
/*                              stiffness matrix                              */
/*                                                                            */
/******************************************************************************/

#if !(DATA_STR & LG_DATA)

 void stiff_matr(mg,tGrid,nu,mult,Z,f,u,br) /*  in case of nonuniform         */
 MULTIGRID *mg;                          /*  refinement, tGrid = TOP_GRID(mg) */
 GRID *tGrid;
 FLOAT nu, mult;
 INT Z, f, u, br;
 {
    GRID *theGrid;
    NODE *pnode, *n0, *n1, *n2, *n3;
    FACE *fa0, *fa1, *fa2, *fa3;
    LINK *plink;
    ELEMENT *pelem;
    FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], ar0, ar1, ar2, ar3, 
          rdetB, ndetB, b[DIM2][DIM2], fmult2;
    fmult2 = FMULT2*mult;
    
    set_mat_value(tGrid,Z,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
    set_value(tGrid,0.,f,T_FOR_U,U_TYPE);
    
    for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      rdetB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                      n2->myvertex->x,n3->myvertex->x,b);
      ndetB = nu*rdetB;
      ar0 = normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      ar1 = normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      ar2 = normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      ar3 = normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b[0],b[1],b[2],b[3],
                                                            ndetB,rdetB,fmult2);
      aijb(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,Z,f,u,b[1],b[2],b[3],b[0],
                                                            ndetB,rdetB,fmult2);
      aijb(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,Z,f,u,b[2],b[3],b[0],b[1],
                                                            ndetB,rdetB,fmult2);
      aijb(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,Z,f,u,b[3],b[0],b[1],b[2],
                                                            ndetB,rdetB,fmult2);
      if (br > -1)
         set_B(pelem,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,
                                ar0,ar1,ar2,ar3,u,br,b[0],b[1],b[2],b[3],rdetB);
    } 
   
    for (theGrid = FIRSTGRID(mg); theGrid->finer != NULL; theGrid = 
                                                                 theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            TOPNODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
            if (ON_FINEST_GRID(n0,mg) || ON_FINEST_GRID(n1,mg) || 
                ON_FINEST_GRID(n2,mg) || ON_FINEST_GRID(n3,mg)){
               TOPFACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);      
               rdetB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                             n2->myvertex->x,n3->myvertex->x,b);
               ndetB = nu*rdetB;
               ar0 = normal_vector(n1->myvertex->x,n2->myvertex->x,
                                   n3->myvertex->x,pelem->f[0],nn0);
               ar1 = normal_vector(n0->myvertex->x,n2->myvertex->x,
                                   n3->myvertex->x,pelem->f[1],nn1);
               ar2 = normal_vector(n0->myvertex->x,n1->myvertex->x,
                                   n3->myvertex->x,pelem->f[2],nn2);
               ar3 = normal_vector(n0->myvertex->x,n1->myvertex->x,
                                   n2->myvertex->x,pelem->f[3],nn3);
               gaijb(mg,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,
                                        b[0],b[1],b[2],b[3],ndetB,rdetB,fmult2);
               gaijb(mg,n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,Z,f,u,
                                        b[1],b[2],b[3],b[0],ndetB,rdetB,fmult2);
               gaijb(mg,n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,Z,f,u,
                                        b[2],b[3],b[0],b[1],ndetB,rdetB,fmult2);
               gaijb(mg,n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,Z,f,u,
                                        b[3],b[0],b[1],b[2],ndetB,rdetB,fmult2);
               if (br > -1)
                  set_B(pelem,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,
                                ar0,ar1,ar2,ar3,u,br,b[0],b[1],b[2],b[3],rdetB);
            }
         }  
    if(TOP_GRID(mg)->coarser)
    for (pnode = FIRSTNODE(TOP_GRID(mg)); pnode != NULL; pnode = pnode->succ)
       for (plink = pnode->start; plink != NULL; plink = plink->next)
          if (IS_EDGE_TO_OLD(plink))
             nn_copy(mg,pnode,plink->nbnode,plink,Z);
 }
 
#endif

#if DATA_STR & LG_DATA

void stiff_matr(mg,tGrid,nu,mult,Z,f,u,br) /*  in case of nonuniform          */
MULTIGRID *mg;                           /*  refinement, tGrid = TOP_GRID(mg) */
GRID *tGrid;
FLOAT nu, mult;
INT Z, f, u, br;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], ar0, ar1, ar2, ar3, 
         rdetB, ndetB, b[DIM2][DIM2], fmult2;
   fmult2 = FMULT2*mult;
   
   set_mat_value(tGrid,Z,0.,-1,-1,U_TYPE,U_TYPE,A_STRUCT);
   set_value(tGrid,0.,f,T_FOR_U,U_TYPE);
   
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      rdetB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                      n2->myvertex->x,n3->myvertex->x,b);
      ndetB = nu*rdetB;
      ar0 = normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      ar1 = normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      ar2 = normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      ar3 = normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      aijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b[0],b[1],b[2],b[3],
                                                            ndetB,rdetB,fmult2);
      aijb(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,Z,f,u,b[1],b[2],b[3],b[0],
                                                            ndetB,rdetB,fmult2);
      aijb(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,Z,f,u,b[2],b[3],b[0],b[1],
                                                            ndetB,rdetB,fmult2);
      aijb(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,Z,f,u,b[3],b[0],b[1],b[2],
                                                            ndetB,rdetB,fmult2);
      if (br > -1)
         set_B(pelem,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,ar0,ar1,ar2,ar3,
                                                u,br,b[0],b[1],b[2],b[3],rdetB);
      if (pelem->n[0]->lgd || pelem->n[1]->lgd || pelem->n[2]->lgd 
                                               || pelem->n[3]->lgd){
         LGij(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,Z,f,u,b[0],b[1],b[2],
                                                              b[3],ndetB,rdetB);
         LGij(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,Z,f,u,b[1],b[2],b[3],
                                                              b[0],ndetB,rdetB);
         LGij(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,Z,f,u,b[2],b[3],b[0],
                                                              b[1],ndetB,rdetB);
         LGij(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,Z,f,u,b[3],b[0],b[1],
                                                              b[2],ndetB,rdetB);
      } 
   }
}
 
#endif

#if (E_DATA & SCALAR_ELEMENT_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

 void mini_stiff_matr(mg,tGrid,nu,mult,Z,Zbn,f,br,u) 
 MULTIGRID *mg;                     /* Z   ... node-node velocity matrix      */
 GRID *tGrid;                       /* Zbn ... node-node pressure matrix      */
 FLOAT nu, mult;                    /* br  ... rhs for continuity equation    */
 INT Z, Zbn, f, br, u;
 {
    NODE *n0, *n1, *n2, *n3;
    ELEMENT *pelem;
    FLOAT rdetB, ndetB, b[DIM2][DIM2];
    
    set_mat_value(tGrid,Z,0.,-1,-1,U_TYPE,U_TYPE,A_STRUCT);
    set_value(tGrid,0.,f,T_FOR_U,U_TYPE);
    set_value(tGrid,0.,br,T_FOR_P,P_TYPE);
    set_mat_value_snXvn(tGrid,Zbn,0.);
    for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      rdetB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                      n2->myvertex->x,n3->myvertex->x,b);
      ndetB = nu*rdetB;
      laijb(n0,n1,n2,n3,Z,f,u,b[0],b[1],b[2],b[3],ndetB,rdetB);
      laijb(n1,n2,n3,n0,Z,f,u,b[1],b[2],b[3],b[0],ndetB,rdetB);
      laijb(n2,n3,n0,n1,Z,f,u,b[2],b[3],b[0],b[1],ndetB,rdetB);
      laijb(n3,n0,n1,n2,Z,f,u,b[3],b[0],b[1],b[2],ndetB,rdetB);
      ndetB = rdetB/4.;
      bijb(pelem,0,n0,n1,n2,n3,Zbn,br,u,b[0],b[1],b[2],b[3],ndetB);
      bijb(pelem,1,n1,n2,n3,n0,Zbn,br,u,b[1],b[2],b[3],b[0],ndetB);
      bijb(pelem,2,n2,n3,n0,n1,Zbn,br,u,b[2],b[3],b[0],b[1],ndetB);
      bijb(pelem,3,n3,n0,n1,n2,Zbn,br,u,b[3],b[0],b[1],b[2],ndetB);
      ndetB = nu*mult*FMULT2*rdetB;
      ED(pelem,Z) = (DOT(b[0],b[0])+DOT(b[1],b[1])+
		       DOT(b[2],b[2])+DOT(b[3],b[3]))*ndetB/15120.;
      integr3B(n0,n1,n2,n3,EDVP(pelem,f),rdetB*FMULT);
    } 
 }
 
#else

void mini_stiff_matr(mg,tGrid,nu,mult,Z,Zbn,f,br,u) 
MULTIGRID *mg; GRID *tGrid; FLOAT nu, mult; INT Z, Zbn, f, br, u;
{  eprintf("Error: mini_stiff_matr not available.\n");  }
 
#endif

#endif  /* DIM == 3 */

/******************************************************************************/
/*                                                                            */
/*                                     NS                                     */
/*                                                                            */
/******************************************************************************/

#if !(DATA_STR & KORN_MATRIX)

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & VECTOR_NODE_DATA) && (DIM == 3)

void put_diag(n0,p,Z)
NODE *n0;
FLOAT p;
INT Z;
{
   COEFFN(n0,Z) += p;
}

void put_nn1(n0,n1,p,u,f,Z)
NODE *n0, *n1;
FLOAT p;
INT u, f, Z;
{
   LINK *pli;   
   
   if (IS_FN(n1)){
      for (pli=n0->start; pli->nbnode != n1; pli=pli->next); 
      COEFFL(pli,Z) += p;
   }
   else if (IS_DN(n1))
      SET9(NDD(n0,f),NDD(n1,u),p)
}  

#else

void put_diag(n0,p,Z)
NODE *n0; FLOAT p; INT Z;
{  eprintf("Error: put_diag not available.\n");  }

void put_nn1(n0,n1,p,u,f,Z)
NODE *n0, *n1; FLOAT p; INT u, f, Z;
{  eprintf("Error: put_nn1 not available.\n");  }

#endif

#endif

#if (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & IxD_FACE_NODE_MATR)

void put_fn1(fan,n1,nn,b1,z,u,f,Z)
FACE *fan;
NODE *n1;
INT u, f, Z;
FLOAT nn[DIM], b1[DIM2], z[DIM];
{
   FNLINK *pfn;
   FLOAT p;
   
   p = DOT(b1,z)*FMULT;
   if (IS_FN(n1)){
      for (pfn=fan->fnstart; pfn->nbnode != n1; pfn=pfn->next);
      SET4(COEFF_FNP(pfn,Z),nn,p)
   }
   else
      FD(fan,f) -= p*DOT(NDD(n1,u),nn);
}

#else

void put_fn1(fan,n1,nn,b1,z,u,f,Z)
FACE *fan; NODE *n1; INT u, f, Z; FLOAT nn[DIM], b1[DIM2], z[DIM];
{  eprintf("Error: put_fn1 not available.\n");  }

#endif

#if DIM == 3

#if DATA_STR & KORN_MATRIX

#if N_DATA & DxD_NODE_MATR

void put_diag(n0,p,Z)
NODE *n0;
FLOAT p;
INT Z;
{
   COEFFNN(n0,Z,0,0) += p;
   COEFFNN(n0,Z,1,1) += p;
   COEFFNN(n0,Z,2,2) += p;
}

#else

void put_diag(n0,p,Z)
NODE *n0; FLOAT p; INT Z;
{  eprintf("Error: put_diag not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES)  && (N_DATA & VECTOR_NODE_DATA)

void put_nn1(n0,n1,p,u,f,Z)
NODE *n0, *n1;
FLOAT p;
INT u, f, Z;
{
   LINK *pli;   

   if (IS_FN(n1)){
      for (pli=n0->start; pli->nbnode != n1; pli=pli->next); 
      COEFFLL(pli,Z,0,0) += p;
      COEFFLL(pli,Z,1,1) += p;
      COEFFLL(pli,Z,2,2) += p;   
   }
   else if (IS_DN(n1))
      SET9(NDD(n0,f),NDD(n1,u),p)
}

#else

void put_nn1(n0,n1,p,u,f,Z)
NODE *n0, *n1; FLOAT p; INT u, f, Z;
{  eprintf("Error: put_nn1 not available.\n");  }

#endif

#endif

#if !(U_SPACE == P1C_NEW_FBUB) && (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA) 

void put_ff1(fan,fa1,nn,nn1,sumv,summ,m0,m1,m2,m3,v1,v2,v3,b0,b1,b2,b3,u,f,Z,q)
FACE *fan, *fa1;
INT u, f, Z;
FLOAT nn[DIM], nn1[DIM], sumv[DIM], summ[DIM], 
      m0[DIM], m1[DIM], m2[DIM], m3[DIM], v1[DIM], v2[DIM], v3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], q;
{
   FLINK *pfli;
   FLOAT p, y[DIM], z[DIM];
   
   SET16(z,v1,v2,v2,v3,v3)
   SET21(y,sumv,15120.,summ,277200.)
   p = DOT(nn,nn1)*(-DOT(b1,y) + 
        (DOT(v2,b3)+DOT(v3,b2))/30240.+DOT(b0,z)/15120. - 
        (DOT(m2,b3)+DOT(m3,b2))/831600. + DOT2(b0,m1,m0)/554400.)*q*FMULT2;
   if (IS_FF(fa1)){
      for (pfli=fan->fstart; pfli->nbface!=fa1; pfli=pfli->next);
      COEFF_FL(pfli,Z) += p;
   }
   else
      FD(fan,f) -= FD(fa1,u)*p; 
}

#else

void put_ff1(fan,fa1,nn,nn1,sumv,summ,m0,m1,m2,m3,v1,v2,v3,b0,b1,b2,b3,u,f,Z,q)
FACE *fan, *fa1; INT u, f, Z; FLOAT nn[DIM], nn1[DIM], sumv[DIM], summ[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM], v1[DIM], v2[DIM], v3[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], q;
{  eprintf("Error: put_ff1 not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) 

#if !(DATA_STR & REDUCED)

/*  x = 3*v0 + v1 + v2 + v3, summ = m0 + m1 + m2 + m3 */
void put_nf1(pnode,fa1,nn1,x,summ,m1,m2,m3,v0,v1,v2,v3,b0,b1,b2,b3,Z,q)
NODE *pnode;
FACE *fa1;
INT Z;
FLOAT nn1[DIM], x[DIM], summ[DIM], m1[DIM], m2[DIM], m3[DIM], 
      v0[DIM], v1[DIM], v2[DIM], v3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], q;
{
   NFLINK *pnf;
   FLOAT y[DIM], p;

   for (pnf=pnode->nfstart; pnf->nbface!=fa1; pnf=pnf->next);
   SET17(y,m1,summ)
   SET21(y,x,420.,y,15120.)
   p = (-DOT(b1,y) + (DOT(v2,b3)+DOT(v3,b2))/420. - DOT2(b0,v1,v0)/840. + 
       (DOT2(b3,m2,m3)+DOT2(b2,m3,m2)-2.*DOT(b0,m1))/30240.)*q*FMULT;
   SET4(COEFF_NFP(pnf,Z),nn1,p)
}

void nijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,m0,m1,m2,m3,
                                               sumv,summ,Z,f,u,b0,b1,b2,b3,detB)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3;
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], 
      m2[DIM], m3[DIM], sumv[DIM], summ[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   NFLINK *pnf;
   FLOAT x[DIM], y[DIM], z[DIM], p;
  
   if (IS_FN(n0)){  
      SET10(x,v0,sumv)
      SET18(y,summ,m0)
      SET20(y,x,detB/20.,y,detB/840.)
      put_diag(n0,DOT(b0,y),Z);
      put_nn1(n0,n1,DOT(b1,y),u,f,Z);
      put_nn1(n0,n2,DOT(b2,y),u,f,Z);
      put_nn1(n0,n3,DOT(b3,y),u,f,Z);
      SET17(x,v0,sumv)
      put_nf1(n0,fa1,nn1,x,summ,m1,m2,m3,v0,v1,v2,v3,b0,b1,b2,b3,Z,detB);
      put_nf1(n0,fa2,nn2,x,summ,m2,m3,m1,v0,v2,v3,v1,b0,b2,b3,b1,Z,detB);
      put_nf1(n0,fa3,nn3,x,summ,m3,m1,m2,v0,v3,v1,v2,b0,b3,b1,b2,Z,detB);
      SET21(x,sumv,420.,summ,15120.)
      p = (-DOT(b0,x) - (DOT(v1,b1)+DOT(v2,b2)+DOT(v3,b3))/840. + 
           (DOT(m1,b1)+DOT(m2,b2)+DOT(m3,b3))/15120.)*detB*FMULT;
      if (NOT_FF(fa0)){
         p *= FD(fa0,u);
         SET9(NDD(n0,f),nn0,p)
      }
      else{
         for (pnf=n0->nfstart; pnf->nbface!=fa0; pnf=pnf->next);
         SET4(COEFF_NFP(pnf,Z),nn0,p)
      }
   }
   
   if (IS_FF(fa0)){
      SET17(y,m0,summ)
      SET19(x,sumv,v0)
      SET21(z,x,15120.,y,277200.)
      COEFF_FF(fa0,Z) += (-DOT(b0,z) - 
                       (DOT(v1,b1)+DOT(v2,b2)+DOT(v3,b3))/15120. + 
                       (DOT(m1,b1)+DOT(m2,b2)+DOT(m3,b3))/554400.)*detB*FMULT2;
      put_ff1(fa0,fa1,nn0,nn1,sumv,summ,m0,m1,m2,m3,v1,v2,v3,b0,b1,b2,b3,u,f,Z,detB);
      put_ff1(fa0,fa2,nn0,nn2,sumv,summ,m0,m2,m3,m1,v2,v3,v1,b0,b2,b3,b1,u,f,Z,detB);
      put_ff1(fa0,fa3,nn0,nn3,sumv,summ,m0,m3,m1,m2,v3,v1,v2,b0,b3,b1,b2,u,f,Z,detB);
      SET18(x,sumv,v0)
      SET10(y,summ,m0)
      SET20(z,x,detB/840.,y,detB/15120.)
      put_fn1(fa0,n0,nn0,b0,z,u,f,Z);
      put_fn1(fa0,n1,nn0,b1,z,u,f,Z);
      put_fn1(fa0,n2,nn0,b2,z,u,f,Z);
      put_fn1(fa0,n3,nn0,b3,z,u,f,Z);
   }
}

#else

void nijb(n0,n1,n2,n3,v0,sumv,Z,f,u,b0,b1,b2,b3,detB)
NODE *n0, *n1, *n2, *n3;
INT Z, f, u;
FLOAT v0[DIM], sumv[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   FLOAT x[DIM], p;
  
   if (IS_FN(n0)){  
      SET10(x,v0,sumv)
      p = detB/20.;
      SET2(x,x,p)
      put_diag(n0,DOT(b0,x),Z);
      put_nn1(n0,n1,DOT(b1,x),u,f,Z);
      put_nn1(n0,n2,DOT(b2,x),u,f,Z);
      put_nn1(n0,n3,DOT(b3,x),u,f,Z);
   }
}

#endif

#else

void nijb(n0,n1,n2,n3,v0,sumv,Z,f,u,b0,b1,b2,b3,detB)
NODE *n0, *n1, *n2, *n3; INT Z, f, u; FLOAT v0[DIM], sumv[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{  eprintf("Error: nijb not available.\n");  }

#endif

#if !(DATA_STR & LG_DATA) && !(DATA_STR & REDUCED) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/* matrix and rhs corresp. to \int_\om w_i\dot(\grad w_j)v\dx; u contains    */
void n_matr(tGrid,Z,f,u,v)                             /*  the Dirichlet b.c. */
GRID *tGrid;
INT Z, f, u, v;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
         v0[DIM], v1[DIM], v2[DIM], v3[3], m0[DIM], m1[DIM], m2[DIM], m3[DIM], 
         sumv[DIM], summ[DIM], detB, b[DIM2][DIM2], p;
  
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      SET1(v0,NDD(n0,v))
      SET1(v1,NDD(n1,v))
      SET1(v2,NDD(n2,v))
      SET1(v3,NDD(n3,v))
      SET3(m0,nn0,FD(fa0,v),FMULT,p)
      SET3(m1,nn1,FD(fa1,v),FMULT,p)
      SET3(m2,nn2,FD(fa2,v),FMULT,p)
      SET3(m3,nn3,FD(fa3,v),FMULT,p)
      SET15(sumv,v0,v1,v2,v3)
      SET15(summ,m0,m1,m2,m3)
      nijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,m0,m1,m2,m3,
                                      sumv,summ,Z,f,u,b[0],b[1],b[2],b[3],detB);
      nijb(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,v1,v2,v3,v0,m1,m2,m3,m0,
                                      sumv,summ,Z,f,u,b[1],b[2],b[3],b[0],detB);
      nijb(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,v2,v3,v0,v1,m2,m3,m0,m1,
                                      sumv,summ,Z,f,u,b[2],b[3],b[0],b[1],detB);
      nijb(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,v3,v0,v1,v2,m3,m0,m1,m2,
                                      sumv,summ,Z,f,u,b[3],b[0],b[1],b[2],detB);
   } 
}

#else

void n_matr(tGrid,Z,f,u,v)
GRID *tGrid; INT Z, f, u, v;
{  eprintf("Error: n_matr not available.\n");  }

#endif

#if !(DATA_STR & LG_DATA) && (DATA_STR & REDUCED)

/* matrix and rhs corresp. to \int_\om w_i\dot(\grad w_j)v\dx; u contains    */
void n_matr(tGrid,Z,f,u,v)                             /*  the Dirichlet b.c. */
GRID *tGrid;
INT Z, f, u, v;
{
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *n3;
   FLOAT v0[DIM], v1[DIM], v2[DIM], v3[DIM], sumv[DIM], detB, b[DIM2][DIM2];
  
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      SET1(v0,NDD(n0,v))
      SET1(v1,NDD(n1,v))
      SET1(v2,NDD(n2,v))
      SET1(v3,NDD(n3,v))
      SET15(sumv,v0,v1,v2,v3)
      nijb(n0,n1,n2,n3,v0,sumv,Z,f,u,b[0],b[1],b[2],b[3],detB);
      nijb(n1,n2,n3,n0,v1,sumv,Z,f,u,b[1],b[2],b[3],b[0],detB);
      nijb(n2,n3,n0,n1,v2,sumv,Z,f,u,b[2],b[3],b[0],b[1],detB);
      nijb(n3,n0,n1,n2,v3,sumv,Z,f,u,b[3],b[0],b[1],b[2],detB);
   } 
}

#endif

#if (DATA_STR & LG_DATA) && !(DATA_STR & REDUCED)

void Put_nlg(pnode,n1,b1,y,Z)
NODE *pnode, *n1;
FLOAT b1[DIM2], y[DIM];
INT Z;
{
   NLGLINK *pnlg;
   FLOAT q;
   
   if (n1->lgd){
      q = DOT(b1,y);
      for (pnlg = pnode->lgstart; pnlg->nbnode != n1; pnlg = pnlg->next);
      COEFF_NLG(pnlg,Z,0,0) += n1->lgd->t1[0]*q;
      COEFF_NLG(pnlg,Z,1,0) += n1->lgd->t1[1]*q;
      COEFF_NLG(pnlg,Z,2,0) += n1->lgd->t1[2]*q;
      COEFF_NLG(pnlg,Z,0,1) += n1->lgd->t2[0]*q;
      COEFF_NLG(pnlg,Z,1,1) += n1->lgd->t2[1]*q;
      COEFF_NLG(pnlg,Z,2,1) += n1->lgd->t2[2]*q;
   }
}

void Put_flg(fan,n1,nn,b1,z,Z)
FACE *fan;
NODE *n1;
FLOAT nn[DIM], b1[DIM2], z[DIM];
INT Z;
{
   FLGLINK *pflg;
   FLOAT q;
   
   if (n1->lgd){
      q = DOT(b1,z)*FMULT;
      for (pflg = fan->lgstart; pflg->nbnode != n1; pflg = pflg->next);
      COEFF_FLG(pflg,Z,0) += q*DOT(n1->lgd->t1,nn);
      COEFF_FLG(pflg,Z,1) += q*DOT(n1->lgd->t2,nn);
   }
}

void Put_lglg(pnode,n1,b1,y,Z)
NODE *pnode, *n1;
FLOAT b1[DIM2], y[DIM];
INT Z;
{
   LGLGLINK *plglg;
   FLOAT q;
   
   if (n1->lgd){
     q = DOT(b1,y);
     for (plglg = pnode->lgd->lgstart; plglg->nbnode != n1; plglg = plglg->next);
     COEFF_LGL(plglg,Z,0,0) += q*DOT(pnode->lgd->t1,n1->lgd->t1);
     COEFF_LGL(plglg,Z,0,1) += q*DOT(pnode->lgd->t1,n1->lgd->t2);
     COEFF_LGL(plglg,Z,1,0) += q*DOT(pnode->lgd->t2,n1->lgd->t1);
     COEFF_LGL(plglg,Z,1,1) += q*DOT(pnode->lgd->t2,n1->lgd->t2);
   }
}

void Put_lgn(pnode,n1,b1,y,Z,f,u)
NODE *pnode, *n1;
FLOAT b1[DIM2], y[DIM];
INT Z, f, u;
{
   LGNLINK *plgn;
   FLOAT q;
   
   q = DOT(b1,y);
   if (IS_FN(n1)){
      for (plgn = pnode->lgd->nstart; plgn->nbnode != n1; plgn = plgn->next);
      COEFF_LGN(plgn,Z,0,0) += pnode->lgd->t1[0]*q;
      COEFF_LGN(plgn,Z,0,1) += pnode->lgd->t1[1]*q;
      COEFF_LGN(plgn,Z,0,2) += pnode->lgd->t1[2]*q;
      COEFF_LGN(plgn,Z,1,0) += pnode->lgd->t2[0]*q;
      COEFF_LGN(plgn,Z,1,1) += pnode->lgd->t2[1]*q;
      COEFF_LGN(plgn,Z,1,2) += pnode->lgd->t2[2]*q;
   }
   else if (n1->lgd == NULL && IS_DN(n1)){
      NDLG(pnode,f,0) -= q*DOT(NDD(n1,u),pnode->lgd->t1);
      NDLG(pnode,f,1) -= q*DOT(NDD(n1,u),pnode->lgd->t2);
   }
}

void Put_lgfan(pnode,fan,nn,summ,m1,m2,m3,sumv,v1,v2,v3,b0,b1,b2,b3,Z,f,u,detB)
NODE *pnode;
FACE *fan;
INT Z, f, u;
FLOAT nn[DIM], summ[DIM], m1[DIM], m2[DIM], m3[DIM], sumv[DIM], 
      v1[DIM], v2[DIM], v3[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   LGFLINK *plgf;
   FLOAT p, p1, p2, x[DIM];
   
   SET21(x,sumv,420.,summ,15120.)
   p = ( -DOT(b0,x) - (DOT(v1,b1)+DOT(v2,b2)+DOT(v3,b3))/840. + 
         (DOT(m1,b1)+DOT(m2,b2)+DOT(m3,b3))/15120. )*detB*FMULT;
   p1 = DOT(pnode->lgd->t1,nn);
   p2 = DOT(pnode->lgd->t2,nn);
   if (IS_FF(fan)){
      for (plgf = pnode->lgd->fstart; plgf->nbface != fan; plgf = plgf->next);
      COEFF_LGF(plgf,Z,0) += p*p1;
      COEFF_LGF(plgf,Z,1) += p*p2;
   }
   else{
      NDLG(pnode,f,0) -= FD(fan,u)*p*p1;
      NDLG(pnode,f,1) -= FD(fan,u)*p*p2;
   }
}

void Put_lgf(pnode,fa1,nn1,x,summ,m1,m2,m3,v0,v1,v2,v3,b0,b1,b2,b3,Z,f,u,detB)
NODE *pnode;
FACE *fa1;
INT Z, f, u;
FLOAT nn1[DIM], x[DIM], summ[DIM], m1[DIM], m2[DIM], m3[DIM], 
      v0[DIM], v1[DIM], v2[DIM], v3[DIM], 
      b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   LGFLINK *plgf;
   FLOAT p, p1, p2, y[DIM];
   
   SET17(y,m1,summ)
   SET21(y,x,420.,y,15120.)
   p = ( -DOT(b1,y) + (DOT(v2,b3)+DOT(v3,b2))/420. - DOT2(b0,v1,v0)/840. + 
         (DOT2(b3,m2,m3)+DOT2(b2,m3,m2)-2.*DOT(b0,m1))/30240. )*detB*FMULT;
   p1 = DOT(pnode->lgd->t1,nn1);
   p2 = DOT(pnode->lgd->t2,nn1);
   if (IS_FF(fa1)){
      for (plgf = pnode->lgd->fstart; plgf->nbface != fa1; plgf = plgf->next);
      COEFF_LGF(plgf,Z,0) += p*p1;
      COEFF_LGF(plgf,Z,1) += p*p2;
   }
   else{
      NDLG(pnode,f,0) -= FD(fa1,u)*p*p1;
      NDLG(pnode,f,1) -= FD(fa1,u)*p*p2;
   }
}

void nLGij(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,m0,m1,m2,m3,
                                               sumv,summ,Z,f,u,b0,b1,b2,b3,detB)
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3; 
INT Z, f, u;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM], 
      sumv[DIM], summ[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   FLOAT x[DIM], y[DIM], z[DIM], q;
         
   if (IS_FN(n0)){
      SET10(x,v0,sumv)
      SET18(y,summ,m0)
      SET20(y,x,detB/20.,y,detB/840.)
      Put_nlg(n0,n1,b1,y,Z);
      Put_nlg(n0,n2,b2,y,Z);
      Put_nlg(n0,n3,b3,y,Z);
   }
   if (IS_FF(fa0)){
      SET18(x,sumv,v0)
      SET10(y,summ,m0)
      SET20(z,x,detB/840.,y,detB/15120.)
      Put_flg(fa0,n0,nn0,b0,z,Z);
      Put_flg(fa0,n1,nn0,b1,z,Z);
      Put_flg(fa0,n2,nn0,b2,z,Z);
      Put_flg(fa0,n3,nn0,b3,z,Z);
   }
   if (n0->lgd){
      SET10(x,v0,sumv)
      SET18(y,summ,m0)
      SET20(y,x,detB/20.,y,detB/840.)
      q = DOT(b0,y);
      COEFF_LG(n0,Z,0,0) += q;
      COEFF_LG(n0,Z,1,1) += q;
      Put_lglg(n0,n1,b1,y,Z);
      Put_lglg(n0,n2,b2,y,Z);
      Put_lglg(n0,n3,b3,y,Z);
      Put_lgn(n0,n1,b1,y,Z,f,u);
      Put_lgn(n0,n2,b2,y,Z,f,u);
      Put_lgn(n0,n3,b3,y,Z,f,u);
      SET17(x,v0,sumv)
      Put_lgfan(n0,fa0,nn0,summ,m1,m2,m3,sumv,v1,v2,v3,b0,b1,b2,b3,Z,f,u,detB);
      Put_lgf(n0,fa1,nn1,x,summ,m1,m2,m3,v0,v1,v2,v3,b0,b1,b2,b3,Z,f,u,detB);
      Put_lgf(n0,fa2,nn2,x,summ,m2,m3,m1,v0,v2,v3,v1,b0,b2,b3,b1,Z,f,u,detB);
      Put_lgf(n0,fa3,nn3,x,summ,m3,m1,m2,v0,v3,v1,v2,b0,b3,b1,b2,Z,f,u,detB);
   }
}
 
/* matrix and rhs corresp. to \int_\om w_i\dot(\grad w_j)v\dx; u contains the */
void n_matr(tGrid,Z,f,u,v)                                /*  Dirichlet b.c.  */
GRID *tGrid;
INT Z, f, u, v;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
         v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM],
         sumv[DIM], summ[DIM], detB, b[DIM2][DIM2], p;
  
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      SET1(v0,NDD(n0,v))
      SET1(v1,NDD(n1,v))
      SET1(v2,NDD(n2,v))
      SET1(v3,NDD(n3,v))
      SET3(m0,nn0,FD(fa0,v),FMULT,p)
      SET3(m1,nn1,FD(fa1,v),FMULT,p)
      SET3(m2,nn2,FD(fa2,v),FMULT,p)
      SET3(m3,nn3,FD(fa3,v),FMULT,p)
      SET15(sumv,v0,v1,v2,v3)
      SET15(summ,m0,m1,m2,m3)
      nijb(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,m0,m1,m2,m3,
                                      sumv,summ,Z,f,u,b[0],b[1],b[2],b[3],detB);
      nijb(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,v1,v2,v3,v0,m1,m2,m3,m0,
                                      sumv,summ,Z,f,u,b[1],b[2],b[3],b[0],detB);
      nijb(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,v2,v3,v0,v1,m2,m3,m0,m1,
                                      sumv,summ,Z,f,u,b[2],b[3],b[0],b[1],detB);
      nijb(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,v3,v0,v1,v2,m3,m0,m1,m2,
                                      sumv,summ,Z,f,u,b[3],b[0],b[1],b[2],detB);
      if (pelem->n[0]->lgd || pelem->n[1]->lgd || 
          pelem->n[2]->lgd || pelem->n[3]->lgd){
         nLGij(n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,
                          m0,m1,m2,m3,sumv,summ,Z,f,u,b[0],b[1],b[2],b[3],detB);
         nLGij(n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,v1,v2,v3,v0,
                          m1,m2,m3,m0,sumv,summ,Z,f,u,b[1],b[2],b[3],b[0],detB);
         nLGij(n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,v2,v3,v0,v1,
                          m2,m3,m0,m1,sumv,summ,Z,f,u,b[2],b[3],b[0],b[1],detB);
         nLGij(n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,v3,v0,v1,v2,
                          m3,m0,m1,m2,sumv,summ,Z,f,u,b[3],b[0],b[1],b[2],detB);
     }
   } 
}

#endif

#if (DATA_STR & LG_DATA) && (DATA_STR & REDUCED)

void n_matr_up(tGrid,Z,f,u,v,q1,q2,q3)
GRID *tGrid; INT Z, f, u, v, q1, q2, q3;
{  eprintf("Error: n_matr_up not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && !(DATA_STR & KORN_MATRIX) && !(DATA_STR & LG_DATA) && !(DATA_STR & REDUCED)

void set_q1(pel,i,fa0,n0,n1,n2,n3,v,q1,q2,q3)
ELEMENT *pel;
FACE *fa0;                /*  fa0 = pel->f[i]              */
NODE *n0, *n1, *n2, *n3;  /*  n1, n2, n3 are nodes of fa0  */
INT i, v, q1, q2, q3;
{
   FLOAT flux, a1, a2, a3, a4;  
   
   if (COEFF_BF(pel,0,i) < 0.0){
     flux = FMULT/60.0*(FD(fa0,v) -       /* flux = \int_fa0 v\dot n\ds/|fa0| */
     ((ND(n1,v,0) + ND(n2,v,0) + ND(n3,v,0))*COEFF_BN(pel,0,i,0) +
      (ND(n1,v,1) + ND(n2,v,1) + ND(n3,v,1))*COEFF_BN(pel,0,i,1) +
      (ND(n1,v,2) + ND(n2,v,2) + ND(n3,v,2))*COEFF_BN(pel,0,i,2))/COEFF_BF(pel,0,i));
     a1 = ND(n1,v,0)*ND(n1,v,0) + ND(n1,v,1)*ND(n1,v,1) + ND(n1,v,2)*ND(n1,v,2);
     a2 = ND(n2,v,0)*ND(n2,v,0) + ND(n2,v,1)*ND(n2,v,1) + ND(n2,v,2)*ND(n2,v,2);
     a3 = ND(n3,v,0)*ND(n3,v,0) + ND(n3,v,1)*ND(n3,v,1) + ND(n3,v,2)*ND(n3,v,2);
     a4 = FMULT*FD(fa0,v)/27.0;
     a1 = 2.0*( MAX(MAX(a1,a2),a3) + a4*a4 );          
     flux *= fabs(flux);
     if (fabs(flux) < a1)
        flux /= a1;
     else
        printf("flux = %e, a1 = %e\n",flux,a1);
     FD(fa0,q1) = flux;
     FD(fa0,q2) = 0.0;
     FD(fa0,q3) = 0.0;
  }
  ND(n0,q2,0) = 0.0;
  ND(n0,q3,0) = 0.0;
}

void Set_q1(pel,i,fa0,n0,n1,n2,n3,v,q1,q2,q3)
ELEMENT *pel;
FACE *fa0;                /*  fa0 = pel->f[i]              */
NODE *n0, *n1, *n2, *n3;  /*  n1, n2, n3 are nodes of fa0  */
INT i, v, q1, q2, q3;
{
   FLOAT flux, a1, a2, a3, a4;  
   
   if (COEFF_BF(pel,0,i) < 0.0){
     flux = FMULT/60.0*(FD(fa0,v) -       /* flux = \int_fa0 v\dot n\ds/|fa0| */
     ((ND(n1,v,0) + ND(n2,v,0) + ND(n3,v,0))*COEFF_BN(pel,0,i,0) +
      (ND(n1,v,1) + ND(n2,v,1) + ND(n3,v,1))*COEFF_BN(pel,0,i,1) +
      (ND(n1,v,2) + ND(n2,v,2) + ND(n3,v,2))*COEFF_BN(pel,0,i,2))/COEFF_BF(pel,0,i));
     a1 = ND(n1,v,0)*ND(n1,v,0) + ND(n1,v,1)*ND(n1,v,1) + ND(n1,v,2)*ND(n1,v,2);
     a2 = ND(n2,v,0)*ND(n2,v,0) + ND(n2,v,1)*ND(n2,v,1) + ND(n2,v,2)*ND(n2,v,2);
     a3 = ND(n3,v,0)*ND(n3,v,0) + ND(n3,v,1)*ND(n3,v,1) + ND(n3,v,2)*ND(n3,v,2);
     a4 = FMULT*FD(fa0,v)/27.0;
     a1 = sqrt(MAX(MAX(a1,a2),a3)) + fabs(a4);
     if (fabs(flux) < a1)
        flux /= a1;
     else
        printf("flux = %e, a1 = %e\n",flux,a1);
     FD(fa0,q1) = flux;
     FD(fa0,q2) = 0.0;
     FD(fa0,q3) = 0.0;
  }
  ND(n0,q2,0) = 0.0;
  ND(n0,q3,0) = 0.0;
}

void SSet_q1(pel,i,fa0,n0,n1,n2,n3,v,q1,q2,q3)
ELEMENT *pel;
FACE *fa0;                /*  fa0 = pel->f[i]              */
NODE *n0, *n1, *n2, *n3;  /*  n1, n2, n3 are nodes of fa0  */
INT i, v, q1, q2, q3;
{
   FLOAT flux, fl, nn[DIM], f[DIM];  
   
   if (COEFF_BF(pel,0,i) < 0.0){
      nn[0] = -COEFF_BN(pel,0,i,0)/COEFF_BF(pel,0,i)/20.0*FMULT;
      nn[1] = -COEFF_BN(pel,0,i,1)/COEFF_BF(pel,0,i)/20.0*FMULT;
      nn[2] = -COEFF_BN(pel,0,i,2)/COEFF_BF(pel,0,i)/20.0*FMULT;
      f[0] = 
        (FD(fa0,v)*nn[0]*FMULT/20.0 + ND(n1,v,0) + ND(n2,v,0) + ND(n3,v,0))/3.0;
      f[1] = 
        (FD(fa0,v)*nn[1]*FMULT/20.0 + ND(n1,v,1) + ND(n2,v,1) + ND(n3,v,1))/3.0;
      f[2] = 
        (FD(fa0,v)*nn[2]*FMULT/20.0 + ND(n1,v,2) + ND(n2,v,2) + ND(n3,v,2))/3.0;
      flux = DOT(f,nn);
      fl = sqrt(DOT(f,f));
      if (fl > 0.0)
         FD(fa0,q1) = flux/fl;
      else
         FD(fa0,q1) = 0.0;
      FD(fa0,q2) = 0.0;
      FD(fa0,q3) = 0.0;
   }
   ND(n0,q2,0) = 0.0;
   ND(n0,q3,0) = 0.0;
}

void SSSet_q1(pel,i,fa0,n0,n1,n2,n3,v,q1,q2,q3)
ELEMENT *pel;
FACE *fa0;                /*  fa0 = pel->f[i]              */
NODE *n0, *n1, *n2, *n3;  /*  n1, n2, n3 are nodes of fa0  */
INT i, v, q1, q2, q3;
{
   FLOAT c0, c1, c2, c3, b[DIM2][DIM2], vv[DIM];  
   
   barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
   SET1(vv,NDD(n0,v))
   c0 = DOT(vv,b[0]);
   if (c0 > 0.0 && 0 <= (c1=-DOT(vv,b[1])) && c1 <= c0
                && 0 <= (c2=-DOT(vv,b[2])) && c2 <= c0
                && 0 <= (c3=-DOT(vv,b[3])) && c3 <= c0)
      ND(n0,q1,0) = (FLOAT)(fa0->index);
   FD(fa0,q2) = 0.0;
   FD(fa0,q3) = 0.0;
   ND(n0,q2,0) = 0.0;
   ND(n0,q3,0) = 0.0;
}

void set_q2_and_q3(pel,i,fa0,n0,q1,q2,q3,detB)
ELEMENT *pel;
FACE *fa0;
NODE *n0;
INT i, q1, q2, q3;
FLOAT detB;
{
   FLOAT flux;
   
   FD(fa0,q2) += detB;
   ND(n0,q2,0) += detB;
   if (COEFF_BF(pel,0,i) < 0.0)
      flux =  FD(fa0,q1);
   else
      flux = -FD(fa0,q1);
   if (flux >= -EPSUP) FD(fa0,q3) += detB;
   if (flux <=  EPSUP) ND(n0,q3,0) += detB;
}

void Set_q2_and_q3(pel,i,fa0,n0,q1,q2,q3,detB)
ELEMENT *pel;
FACE *fa0;
NODE *n0;
INT i, q1, q2, q3;
FLOAT detB;
{
   FLOAT flux;
   
   FD(fa0,q2) += detB;
   ND(n0,q2,0) += detB;
   if (COEFF_BF(pel,0,i) < 0.0)
      flux =  FD(fa0,q1);
   else
      flux = -FD(fa0,q1);
   FD(fa0,q3) += detB*(KKK - flux);
   ND(n0,q3,0) += detB*(KKK + flux);
}

void set_auxiliary_fields(tGrid,q1,q2,q3,v)
GRID *tGrid;
INT q1, q2, q3, v;
{
   ELEMENT *pelem;
   NODE *pnode, *n0, *n1, *n2, *n3;
   FACE *pface, *fa0, *fa1, *fa2, *fa3;
   FLOAT detB, b[DIM2][DIM2];
    
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      SSet_q1(pelem,0,fa0,n0,n1,n2,n3,v,q1,q2,q3);
      SSet_q1(pelem,1,fa1,n1,n2,n3,n0,v,q1,q2,q3);
      SSet_q1(pelem,2,fa2,n2,n3,n0,n1,v,q1,q2,q3);
      SSet_q1(pelem,3,fa3,n3,n0,n1,n2,v,q1,q2,q3);
   } 
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      Set_q2_and_q3(pelem,0,fa0,n0,q1,q2,q3,detB);
      Set_q2_and_q3(pelem,1,fa1,n1,q1,q2,q3,detB);
      Set_q2_and_q3(pelem,2,fa2,n2,q1,q2,q3,detB);
      Set_q2_and_q3(pelem,3,fa3,n3,q1,q2,q3,detB);
   } 
   for (pnode=FIRSTNODE(tGrid);pnode!=NULL;pnode=pnode->succ)
      ND(pnode,q2,0) /= ND(pnode,q3,0);
   for (pface=FIRSTFACE(tGrid);pface!=NULL;pface=pface->succ)
      FD(pface,q2) /= FD(pface,q3);
}

void nijb_up(pelem,i,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,
                             m0,m1,m2,m3,sumv,summ,Z,f,u,q1,q2,b0,b1,b2,b3,detB)
ELEMENT *pelem;
NODE *n0, *n1, *n2, *n3;
FACE *fa0, *fa1, *fa2, *fa3;
INT i, Z, f, u, q1, q2;
FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
      v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM], 
      sumv[DIM], summ[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB;
{
   NFLINK *pnf;
   FLOAT x[DIM], y[DIM], z[DIM], flux, p, q;
  
   if (COEFF_BF(pelem,0,i) < 0.0)
      flux =  FD(fa0,q1);
   else
      flux = -FD(fa0,q1);
   
/*   if (IS_FN(n0) && flux <= EPSUP){
      q = ND(n0,q2,0)*detB;*/
   if (IS_FN(n0)){
      q = ND(n0,q2,0)*(KKK - flux)*detB;
      SET10(x,v0,sumv)
      SET18(y,summ,m0)
      SET20(y,x,q/20.,y,q/840.);
      COEFFN(n0,Z) += DOT(b0,y);
      put_nn1(n0,n1,DOT(b1,y),u,f,Z);
      put_nn1(n0,n2,DOT(b2,y),u,f,Z);
      put_nn1(n0,n3,DOT(b3,y),u,f,Z);
      SET17(x,v0,sumv)
      put_nf1(n0,fa1,nn1,x,summ,m1,m2,m3,v0,v1,v2,v3,b0,b1,b2,b3,Z,q);
      put_nf1(n0,fa2,nn2,x,summ,m2,m3,m1,v0,v2,v3,v1,b0,b2,b3,b1,Z,q);
      put_nf1(n0,fa3,nn3,x,summ,m3,m1,m2,v0,v3,v1,v2,b0,b3,b1,b2,Z,q);
      SET21(x,sumv,420.,summ,15120.)
      p = (-DOT(b0,x) - (DOT(v1,b1)+DOT(v2,b2)+DOT(v3,b3))/840. + 
           (DOT(m1,b1)+DOT(m2,b2)+DOT(m3,b3))/15120.)*q*FMULT;
      if (NOT_FF(fa0)){
         ND(n0,f,0) -= FD(fa0,u)*p*nn0[0];
         ND(n0,f,1) -= FD(fa0,u)*p*nn0[1];
         ND(n0,f,2) -= FD(fa0,u)*p*nn0[2];
      }
      else{
         for (pnf=n0->nfstart; pnf->nbface!=fa0; pnf=pnf->next);
         COEFF_NF(pnf,Z,0) += p*nn0[0];
         COEFF_NF(pnf,Z,1) += p*nn0[1];
         COEFF_NF(pnf,Z,2) += p*nn0[2];  
      }
   }
   
/*   if (IS_FF(fa0) && flux >= -EPSUP){
      q = FD(fa0,q2)*detB;*/
   if (IS_FF(fa0)){
     q = FD(fa0,q2)*(KKK + flux)*detB;
     SET17(y,m0,summ)
     SET19(x,sumv,v0)
     SET21(z,x,15120.,y,277200.)
     COEFF_FF(fa0,Z) += (-DOT(b0,z) - 
                      (DOT(v1,b1)+DOT(v2,b2)+DOT(v3,b3))/15120. + 
                      (DOT(m1,b1)+DOT(m2,b2)+DOT(m3,b3))/554400.)*q*FMULT2;
     put_ff1(fa0,fa1,nn0,nn1,sumv,summ,m0,m1,m2,m3,v1,v2,v3,b0,b1,b2,b3,u,f,Z,q);
     put_ff1(fa0,fa2,nn0,nn2,sumv,summ,m0,m2,m3,m1,v2,v3,v1,b0,b2,b3,b1,u,f,Z,q);
     put_ff1(fa0,fa3,nn0,nn3,sumv,summ,m0,m3,m1,m2,v3,v1,v2,b0,b3,b1,b2,u,f,Z,q);
     SET18(x,sumv,v0)
     SET10(y,summ,m0)
     SET20(z,x,q/840.,y,q/15120.)
     put_fn1(fa0,n0,nn0,b0,z,u,f,Z);
     put_fn1(fa0,n1,nn0,b1,z,u,f,Z);
     put_fn1(fa0,n2,nn0,b2,z,u,f,Z);
     put_fn1(fa0,n3,nn0,b3,z,u,f,Z);
   }
}

void n_matr_up(tGrid,Z,f,u,v,q1,q2,q3)  /*  n_matr with upwinding  */
GRID *tGrid;
INT Z, f, u, v, q1, q2, q3; /* q1, q2, q3 are auxiliary node-face fields */
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
         v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM],
         sumv[DIM], summ[DIM], detB, b[DIM2][DIM2], p;
  
   set_auxiliary_fields(tGrid,q1,q2,q3,v);
   
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      SET1(v0,NDD(n0,v))
      SET1(v1,NDD(n1,v))
      SET1(v2,NDD(n2,v))
      SET1(v3,NDD(n3,v))
      SET3(m0,nn0,FD(fa0,v),FMULT,p)
      SET3(m1,nn1,FD(fa1,v),FMULT,p)
      SET3(m2,nn2,FD(fa2,v),FMULT,p)
      SET3(m3,nn3,FD(fa3,v),FMULT,p)
      SET15(sumv,v0,v1,v2,v3)
      SET15(summ,m0,m1,m2,m3)
      nijb_up(pelem,0,n0,n1,n2,n3,fa0,fa1,fa2,fa3,nn0,nn1,nn2,nn3,v0,v1,v2,v3,
                    m0,m1,m2,m3,sumv,summ,Z,f,u,q1,q2,b[0],b[1],b[2],b[3],detB);
      nijb_up(pelem,1,n1,n2,n3,n0,fa1,fa2,fa3,fa0,nn1,nn2,nn3,nn0,v1,v2,v3,v0,
                    m1,m2,m3,m0,sumv,summ,Z,f,u,q1,q2,b[1],b[2],b[3],b[0],detB);
      nijb_up(pelem,2,n2,n3,n0,n1,fa2,fa3,fa0,fa1,nn2,nn3,nn0,nn1,v2,v3,v0,v1,
                    m2,m3,m0,m1,sumv,summ,Z,f,u,q1,q2,b[2],b[3],b[0],b[1],detB);
      nijb_up(pelem,3,n3,n0,n1,n2,fa3,fa0,fa1,fa2,nn3,nn0,nn1,nn2,v3,v0,v1,v2,
                    m3,m0,m1,m2,sumv,summ,Z,f,u,q1,q2,b[3],b[0],b[1],b[2],detB);
   } 
}

#else

void n_matr_up(tGrid,Z,f,u,v,q1,q2,q3)
GRID *tGrid;
INT Z, f, u, v, q1, q2, q3;
{  eprintf("Error: n_matr_up not available.\n");  }

#endif

#if !(DATA_STR & KORN_MATRIX) && !(DATA_STR & LG_DATA) && (N_DATA & ONE_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & IxD_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) 

void add_flux2(a,v,n3,n4,x3,x4,x1,x2)
NODE *n3, *n4;
FLOAT *a, *x3, *x4, *x1, *x2;
INT v;
{
   FLOAT m34[DIM], x34[DIM], flux34;
   LINK *pli;
   
   AVERAGE(x3,x4,x34);
   m34[0] = (x34[1]-x1[1])*(x34[2]-x2[2]) - (x34[2]-x1[2])*(x34[1]-x2[1]);
   m34[1] = (x34[2]-x1[2])*(x34[0]-x2[0]) - (x34[0]-x1[0])*(x34[2]-x2[2]);
   m34[2] = (x34[0]-x1[0])*(x34[1]-x2[1]) - (x34[1]-x1[1])*(x34[0]-x2[0]);
   flux34 = ( (a[0] + ND(n3,v,0) + ND(n4,v,0))*m34[0] + 
              (a[1] + ND(n3,v,1) + ND(n4,v,1))*m34[1] + 
              (a[2] + ND(n3,v,2) + ND(n4,v,2))*m34[2] )/54.;
   if (m34[0]*(x4[0]-x3[0]) + m34[1]*(x4[1]-x3[1]) + m34[2]*(x4[2]-x3[2]) < 0.) 
                                                               flux34 = -flux34;
   if (IS_FN(n4)){
      for (pli = n3->start; pli->nbnode != n4; pli = pli->next);
      pli->flux += flux34;
   }
   if (IS_FN(n3)){
      for (pli = n4->start; pli->nbnode != n3; pli = pli->next);
      pli->flux -= flux34;
   }
}

void add_flux3(a,v,a12,a34,n3,n4,x3,x4,x1,x2)
NODE *n3, *n4;
FLOAT *a, *a12, *a34, *x3, *x4, *x1, *x2;
INT v;
{
   FLOAT m34[DIM], x34[DIM], flux34;
   LINK *pli;
   
   AVERAGE(x3,x4,x34);
   m34[0] = (x34[1]-x1[1])*(x34[2]-x2[2]) - (x34[2]-x1[2])*(x34[1]-x2[1]);
   m34[1] = (x34[2]-x1[2])*(x34[0]-x2[0]) - (x34[0]-x1[0])*(x34[2]-x2[2]);
   m34[2] = (x34[0]-x1[0])*(x34[1]-x2[1]) - (x34[1]-x1[1])*(x34[0]-x2[0]);
   flux34 = ( (a[0] + ND(n3,v,0) + ND(n4,v,0))*m34[0] + 
              (a[1] + ND(n3,v,1) + ND(n4,v,1))*m34[1] + 
              (a[2] + ND(n3,v,2) + ND(n4,v,2))*m34[2] )/54. +
            ( 285.*DOT(a12,m34) + 97.*DOT(a34,m34) )/207360.;
   if (m34[0]*(x4[0]-x3[0]) + m34[1]*(x4[1]-x3[1]) + m34[2]*(x4[2]-x3[2]) < 0.) 
                                                               flux34 = -flux34;
   if (IS_FN(n4)){
      for (pli = n3->start; pli->nbnode != n4; pli = pli->next);
      pli->flux += flux34;
   }
   if (IS_FN(n3)){
      for (pli = n4->start; pli->nbnode != n3; pli = pli->next);
      pli->flux -= flux34;
   }
}

FLOAT llambda(x)
FLOAT x;
{
   if (x < 0.)
      return(0.);
   else
      return(1.);
}  

FLOAT lambda(x)
FLOAT x;
{
   if (x < 0.)
      return(0.5*NU/(NU - x));
   else
      return((0.5*NU + x)/(NU + x));
} 

/* matrix and rhs corresp. to \int_\om w_i\dot(\grad w_j)v\dx; u contains     */
void n_matr_up2(tGrid,Z,f,u,v)                         /*  the Dirichlet b.c. */
GRID *tGrid;
INT Z, f, u, v;
{
   NODE *n0, *n1, *n2, *n3, *pnode, *nbnode;
   ELEMENT *pelem;
   LINK *pli;
   FLOAT a[DIM], *x0, *x1, *x2, *x3;
  
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next)
         pli->flux = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      x0 = n0->myvertex->x;
      x1 = n1->myvertex->x;
      x2 = n2->myvertex->x;
      x3 = n3->myvertex->x;
      a[0] = (ND(n0,v,0) + ND(n1,v,0) + ND(n2,v,0) + ND(n3,v,0))*0.625;
      a[1] = (ND(n0,v,1) + ND(n1,v,1) + ND(n2,v,1) + ND(n3,v,1))*0.625;
      a[2] = (ND(n0,v,2) + ND(n1,v,2) + ND(n2,v,2) + ND(n3,v,2))*0.625;
      add_flux2(a,v,n0,n1,x0,x1,x2,x3);
      add_flux2(a,v,n0,n2,x0,x2,x1,x3);
      add_flux2(a,v,n0,n3,x0,x3,x1,x2);
      add_flux2(a,v,n1,n2,x1,x2,x0,x3);
      add_flux2(a,v,n1,n3,x1,x3,x0,x2);
      add_flux2(a,v,n2,n3,x2,x3,x0,x1);
   }
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next)
         if (IS_FN(nbnode=pli->nbnode)){
            COEFFN(nbnode,Z) += pli->flux*(0.5 - lambda(-pli->flux));
            ND(nbnode,f,0) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,0);
            ND(nbnode,f,1) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,1);
            ND(nbnode,f,2) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,2);
         }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         COEFFN(pli->nbnode,Z) += pli->flux*(0.5 - lambda(-pli->flux));
         COEFFL(pli,Z) += pli->flux*(1. - lambda(pli->flux));
      }
}

/* matrix and rhs corresp. to \int_\om w_i\dot(\grad w_j)v\dx; u contains     */
void n_matr_up3(tGrid,Z,f,u,v)                         /*  the Dirichlet b.c. */
GRID *tGrid;
INT Z, f, u, v;
{
   NODE *n0, *n1, *n2, *n3, *pnode, *nbnode;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   LINK *pli;
   FLOAT a[DIM],  a01[DIM], a02[DIM], a03[DIM], a12[DIM], a13[DIM], a23[DIM], 
         nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], *x0, *x1, *x2, *x3;
  
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next)
         pli->flux = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      x0 = n0->myvertex->x;
      x1 = n1->myvertex->x;
      x2 = n2->myvertex->x;
      x3 = n3->myvertex->x;
      normal_vector(x1,x2,x3,fa0,nn0);
      normal_vector(x0,x2,x3,fa1,nn1);
      normal_vector(x0,x1,x3,fa2,nn2);
      normal_vector(x0,x1,x2,fa3,nn3);
      a[0] = (ND(n0,v,0) + ND(n1,v,0) + ND(n2,v,0) + ND(n3,v,0))*0.625;
      a[1] = (ND(n0,v,1) + ND(n1,v,1) + ND(n2,v,1) + ND(n3,v,1))*0.625;
      a[2] = (ND(n0,v,2) + ND(n1,v,2) + ND(n2,v,2) + ND(n3,v,2))*0.625;
      ADDMULT(FD(fa0,v),nn0,FD(fa1,v),nn1,a01);
      ADDMULT(FD(fa0,v),nn0,FD(fa2,v),nn2,a02);
      ADDMULT(FD(fa0,v),nn0,FD(fa3,v),nn3,a03);
      ADDMULT(FD(fa1,v),nn1,FD(fa2,v),nn2,a12);
      ADDMULT(FD(fa1,v),nn1,FD(fa3,v),nn3,a13);
      ADDMULT(FD(fa2,v),nn2,FD(fa3,v),nn3,a23);
      add_flux3(a,v,a23,a01,n0,n1,x0,x1,x2,x3);
      add_flux3(a,v,a13,a02,n0,n2,x0,x2,x1,x3);
      add_flux3(a,v,a12,a03,n0,n3,x0,x3,x1,x2);
      add_flux3(a,v,a03,a12,n1,n2,x1,x2,x0,x3);
      add_flux3(a,v,a02,a13,n1,n3,x1,x3,x0,x2);
      add_flux3(a,v,a01,a23,n2,n3,x2,x3,x0,x1);
   }
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next)
         if (IS_FN(nbnode=pli->nbnode)){
            COEFFN(nbnode,Z) += pli->flux*(0.5 - lambda(-pli->flux));
            ND(nbnode,f,0) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,0);
            ND(nbnode,f,1) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,1);
            ND(nbnode,f,2) -= pli->flux*(lambda(-pli->flux) - 1.)*ND(pnode,u,2);
         }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      for (pli = pnode->start; pli != NULL; pli = pli->next){
         COEFFN(pli->nbnode,Z) += pli->flux*(0.5 - lambda(-pli->flux));
         COEFFL(pli,Z) += pli->flux*(1. - lambda(pli->flux));
      }
}

#endif

#if DATA_STR & KORN_MATRIX || DATA_STR & LG_DATA

void n_matr_up2(tGrid,Z,f,u,v,q1,q2,q3)
GRID *tGrid;
INT Z, f, u, v, q1, q2, q3;
{  eprintf("Error: n_matr_up2 not available.\n");  }

#endif

#endif  /* DIM == 3 */

/******************************************************************************/
/*                                                                            */
/*          stiffness matrices for scalar piecewise linear functions          */
/*                                                                            */
/******************************************************************************/

#if N_DATA & SCALAR_NODE_DATA

void sbound_val(mg,u,u0,t)
MULTIGRID *mg;
FLOAT (*u0)();
INT u, t;
{
   GRID *theGrid;
   NODE *pnode;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode != NULL; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode) && IS_DIR(pnode,t))
            NDS(pnode,u) = u0(pnode->myvertex->x);
}

#else

void sbound_val(mg,u,u0,t)
MULTIGRID *mg; FLOAT (*u0)(); INT u, t;
{  eprintf("Error: sbound_val not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (N_DATA & ONE_NODE_MATR)

void sputaij(plin,n1,n2,an1,an2,Z)
LINK *plin;
NODE *n1, *n2;
INT Z;
FLOAT an1, an2;
{
  LINK *pli;

  for (pli=plin; NBNODE(pli) != n1  && NBNODE(pli) != n2 ; pli=pli->next);
  if (NBNODE(pli) == n1){
     COEFFLS(pli,Z) += an1;
     for (pli=pli->next; NBNODE(pli) != n2 ; pli=pli->next); 
     COEFFLS(pli,Z) += an2;
  }
  else{
     COEFFLS(pli,Z) += an2;
     for (pli=pli->next; NBNODE(pli) != n1 ; pli=pli->next); 
     COEFFLS(pli,Z) += an1;
  }
}
 
#endif  /*  (N_DATA & SCALAR_NODE_DATA) && (N_DATA & ONE_NODE_MATR)  */

#if DIM == 3

#if (N_DATA & SCALAR_NODE_DATA) && (N_DATA & ONE_NODE_MATR)

void saij(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,t)
NODE *n0, *n1, *n2, *n3;
INT Z, f, u, t;
FLOAT ann, an1, an2, an3;
{
   LINK *plin, *pli;

   COEFFS(n0,Z) += ann;
   plin = TSTART(n0);
   if (!IS_DIR(n1,t) && !IS_DIR(n2,t) && !IS_DIR(n3,t)){
      for (pli=plin; NBNODE(pli) != n1 && NBNODE(pli) != n2
                                && NBNODE(pli) != n3 ; pli=pli->next);
      if (NBNODE(pli) == n1){
         COEFFLS(pli,Z) += an1;
         sputaij(pli->next,n2,n3,an2,an3,Z);
      }            
      else if(NBNODE(pli) == n2){
         COEFFLS(pli,Z) += an2;
         sputaij(pli->next,n1,n3,an1,an3,Z);
      }
      else{    
         COEFFLS(pli,Z) += an3;
         sputaij(pli->next,n1,n2,an1,an2,Z);
      }
   }
   else if (!IS_DIR(n1,t) && !IS_DIR(n2,t)){
      sputaij(plin,n1,n2,an1,an2,Z);
      NDS(n0,f) -= an3*NDS(n3,u);
   }
   else if (!IS_DIR(n1,t) && !IS_DIR(n3,t)){
      sputaij(plin,n1,n3,an1,an3,Z);
      NDS(n0,f) -= an2*NDS(n2,u);
   }
   else if (!IS_DIR(n2,t) && !IS_DIR(n3,t)){ 
      sputaij(plin,n2,n3,an2,an3,Z);
      NDS(n0,f) -= an1*NDS(n1,u);
   }
   else if (!IS_DIR(n1,t)){
      for (pli=plin; NBNODE(pli) != n1 ; pli=pli->next);
      COEFFLS(pli,Z) += an1;
      NDS(n0,f) -= an2*NDS(n2,u) + an3*NDS(n3,u);
   }                                         
   else if (!IS_DIR(n2,t)){
      for (pli=plin; NBNODE(pli) != n2 ; pli=pli->next);
      COEFFLS(pli,Z) += an2; 
      NDS(n0,f) -= an1*NDS(n1,u) + an3*NDS(n3,u);
   }                                       
   else if (!IS_DIR(n3,t)){
      for (pli=plin; NBNODE(pli) != n3 ; pli=pli->next);
      COEFFLS(pli,Z) += an3;
      NDS(n0,f) -= an1*NDS(n1,u) + an2*NDS(n2,u);
   }
   else
      NDS(n0,f) -= an1*NDS(n1,u) + an2*NDS(n2,u) + an3*NDS(n3,u);
}

void saijb(n0,n1,n2,n3,Z,f,u,b0,b1,b2,b3,detB,rdetB,t,t0)
NODE *n0, *n1, *n2, *n3;
INT Z, f, u, t;
FLOAT b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, rdetB, (*t0)();
{
   FLOAT ann, an1, an2, an3;
  
   if (!IS_DIR(n0,t)){  
      ann = DOT(b0,b0)*detB;
      an1 = DOT(b0,b1)*detB;
      an2 = DOT(b0,b2)*detB;
      an3 = DOT(b0,b3)*detB;
    
      NDS(n0,f) += integr1(n0->myvertex->x,n1->myvertex->x,
                           n2->myvertex->x,n3->myvertex->x,t0,rdetB);
      saij(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,t);
   }
}
  
/*  stiffness matrix on one grid corresponding to 
                                       nu*\int_\om \nabla u\cdot\nabla v \dx  */
void sstiff_matr(tGrid,nu,Z,f,u,t,t0)
GRID *tGrid;
FLOAT nu, (*t0)();
INT Z, f, u, t;
{
   NODE *n0, *n1, *n2, *n3;
   ELEMENT *pelem;
   FLOAT rdetB, ndetB, b[DIM2][DIM2];

   set_mat_value(tGrid,Z,0.,t,t,Q_SN,Q_SN,-1);
   set_value(tGrid,0.,f,T_FOR_U,U_TYPE);
   
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      rdetB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                      n2->myvertex->x,n3->myvertex->x,b);
      ndetB = nu*rdetB;
      saijb(n0,n1,n2,n3,Z,f,u,b[0],b[1],b[2],b[3],ndetB,rdetB,t,t0);
      saijb(n1,n2,n3,n0,Z,f,u,b[1],b[2],b[3],b[0],ndetB,rdetB,t,t0);
      saijb(n2,n3,n0,n1,Z,f,u,b[2],b[3],b[0],b[1],ndetB,rdetB,t,t0);
      saijb(n3,n0,n1,n2,Z,f,u,b[3],b[0],b[1],b[2],ndetB,rdetB,t,t0);
   }
}

void snijb(n0,n1,n2,n3,vn,sumv,mn,summ,Z,f,u,b0,b1,b2,b3,detB,t)
NODE *n0, *n1, *n2, *n3;
INT Z, f, u, t;
FLOAT b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2], detB, 
      vn[DIM], sumv[DIM], mn[DIM], summ[DIM];
{
   FLOAT ann, an1, an2, an3, sum[DIM];
  
   if (!IS_DIR(n0,t)){
      sum[0] = (vn[0]+sumv[0] + (summ[0]+summ[0]-mn[0])*FMULT/42.)*detB/20.;
      sum[1] = (vn[1]+sumv[1] + (summ[1]+summ[1]-mn[1])*FMULT/42.)*detB/20.;
      sum[2] = (vn[2]+sumv[2] + (summ[2]+summ[2]-mn[2])*FMULT/42.)*detB/20.;
      ann = DOT(sum,b0);
      an1 = DOT(sum,b1);
      an2 = DOT(sum,b2);
      an3 = DOT(sum,b3);
      saij(n0,n1,n2,n3,ann,an1,an2,an3,Z,f,u,t);
   }
}
  
#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/*  stiffness matrix on one grid corresponding to 
                                          \int_\om u_i(v\cdot\nabla u_j) \dx  */
void sn_matr(tGrid,Z,f,u,v,t)
GRID *tGrid;
INT Z,f,u,t;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT b[DIM2][DIM2], nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
         v0[DIM], v1[DIM], v2[DIM], v3[DIM], m0[DIM], m1[DIM], m2[DIM], m3[DIM],
         sumv[DIM], summ[DIM], detB;

   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      SET1(v0,NDD(n0,v))
      SET1(v1,NDD(n1,v))
      SET1(v2,NDD(n2,v))
      SET1(v3,NDD(n3,v))
      SET2(m0,nn0,FD(fa0,v))
      SET2(m1,nn1,FD(fa1,v))
      SET2(m2,nn2,FD(fa2,v))
      SET2(m3,nn3,FD(fa3,v))
      SET15(sumv,v0,v1,v2,v3)
      SET15(summ,m0,m1,m2,m3)
      snijb(n0,n1,n2,n3,v0,sumv,m0,summ,Z,f,u,b[0],b[1],b[2],b[3],detB,t);
      snijb(n1,n2,n3,n0,v1,sumv,m1,summ,Z,f,u,b[1],b[2],b[3],b[0],detB,t);
      snijb(n2,n3,n0,n1,v2,sumv,m2,summ,Z,f,u,b[2],b[3],b[0],b[1],detB,t);
      snijb(n3,n0,n1,n2,v3,sumv,m3,summ,Z,f,u,b[3],b[0],b[1],b[2],detB,t);
   }
}

#endif

FLOAT tarea(n1,n2,n3)
NODE *n1, *n2, *n3;
{
   FLOAT *x1, *x2, *x3, a[DIM], b[DIM], n[DIM]; 
                   
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   
   a[0] = x2[0]-x1[0];
   a[1] = x2[1]-x1[1];
   a[2] = x2[2]-x1[2];
   b[0] = x3[0]-x2[0];
   b[1] = x3[1]-x2[1];
   b[2] = x3[2]-x2[2];
   n[0] = a[1]*b[2]-a[2]*b[1];
   n[1] = a[2]*b[0]-a[0]*b[2];
   n[2] = a[0]*b[1]-a[1]*b[0];

   return(sqrt(DOT(n,n))/2.);
}

FLOAT integr6(n1,n2,n3,g) /* \int_K' g \lambda_1 \dx exact for g\in P_3 */
NODE *n1, *n2, *n3;
FLOAT (*g)();
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x123[DIM]; 
                   
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   POINT2(x1,x2,x112,2.0,3.0);
   POINT2(x1,x3,x113,2.0,3.0);
   POINT3(x1,x2,x3,x123);
   
   return((2.*(*g)(x1) + (*g)(x2) + (*g)(x3) + 
           9.*((*g)(x112) + (*g)(x113)) + 18.*(*g)(x123))*tarea(n1,n2,n3)/120.);
} 

FLOAT integr7(n1,n2,n3,g) /* \int_K' g \lambda_1 \lambda_2 \dx exact          */
NODE *n1, *n2, *n3;                                       /*    for g\in P_3  */
FLOAT (*g)();
{
  FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], 
                                  x331[DIM], x332[DIM], x123[DIM]; 
 
  x1 = n1->myvertex->x;
  x2 = n2->myvertex->x;
  x3 = n3->myvertex->x;
  POINT2(x1,x2,x112,2.0,3.0);
  POINT2(x1,x3,x113,2.0,3.0);
  POINT2(x2,x1,x221,2.0,3.0);
  POINT2(x2,x3,x223,2.0,3.0);
  POINT2(x3,x1,x331,2.0,3.0);
  POINT2(x3,x2,x332,2.0,3.0);
  POINT3(x1,x2,x3,x123);
  
  return(((*g)(x1) + (*g)(x2) + 2.*(*g)(x3) + 
         12.*((*g)(x112) + (*g)(x221)) + 6.*((*g)(x113) + (*g)(x223)) 
         - 3.*((*g)(x331) + (*g)(x332)) + 36.*(*g)(x123))*tarea(n1,n2,n3)/840.);
} 

FLOAT integr8(n1,n2,n3,g) /* \int_K' g (\lambda_1)^2 \dx exact for g\in P_3 */
NODE *n1, *n2, *n3;
FLOAT (*g)();
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], 
                                   x331[DIM], x332[DIM], x123[DIM]; 
                   
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   POINT2(x1,x2,x112,2.0,3.0);
   POINT2(x1,x3,x113,2.0,3.0);
   POINT2(x2,x1,x221,2.0,3.0);
   POINT2(x2,x3,x223,2.0,3.0);
   POINT2(x3,x1,x331,2.0,3.0);
   POINT2(x3,x2,x332,2.0,3.0);
   POINT3(x1,x2,x3,x123);
   
   return((( 3.*(*g)(x1) + (*g)(x2) + (*g)(x3) )/210. + 
          ( 15.*((*g)(x112) + (*g)(x113)) - 3.*((*g)(x221) + (*g)(x331)) 
          - ((*g)(x223) + (*g)(x332)) + 18.*(*g)(x123) )/280.)*tarea(n1,n2,n3));
} 

FLOAT integr9(n1,n2,n3,u1,u2,u3,g) /* \int_K' g(u) \lambda_1 \lambda_2 \dx    */
NODE *n1, *n2, *n3;                /*                      exact for g\in P_3 */
FLOAT u1, u2, u3, (*g)();          /*  u is linear on K'; u(n_i) = u_i        */
{
   FLOAT u112, u113, u221, u223, u331, u332, u123; 
                   
   u112 = (2.*u1 + u2)/3.;
   u221 = (2.*u2 + u1)/3.;
   u113 = (2.*u1 + u3)/3.;
   u331 = (2.*u3 + u1)/3.;
   u223 = (2.*u2 + u3)/3.;
   u332 = (2.*u3 + u2)/3.;
   u123 = (u1 + u2 + u3)/3.;
   
   return(((*g)(u1) + (*g)(u2) + 2.*(*g)(u3) + 
          12.*((*g)(u112) + (*g)(u221)) + 6.*((*g)(u113) + (*g)(u223)) 
         - 3.*((*g)(u331) + (*g)(u332)) + 36.*(*g)(u123))*tarea(n1,n2,n3)/840.);
} 

FLOAT integr10(n1,n2,n3,u1,u2,u3,g) /* \int_K' g(u) (\lambda_1)^2 \dx         */
NODE *n1, *n2, *n3;                 /*                     exact for g\in P_3 */
FLOAT u1, u2, u3, (*g)();           /*  u is linear on K'; u(n_i) = u_i       */
{
   FLOAT u112, u113, u221, u223, u331, u332, u123; 
                   
   u112 = (2.*u1 + u2)/3.;
   u221 = (2.*u2 + u1)/3.;
   u113 = (2.*u1 + u3)/3.;
   u331 = (2.*u3 + u1)/3.;
   u223 = (2.*u2 + u3)/3.;
   u332 = (2.*u3 + u2)/3.;
   u123 = (u1 + u2 + u3)/3.;
   
   return((( 3.*(*g)(u1) + (*g)(u2) + (*g)(u3) )/210. + 
          ( 15.*((*g)(u112) + (*g)(u113)) - 3.*((*g)(u221) + (*g)(u331)) 
          - ((*g)(u223) + (*g)(u332)) + 18.*(*g)(u123) )/280.)*tarea(n1,n2,n3));
}

void sNeumann_bc_on_face(n1,n2,n3,g,f,t)
NODE *n1, *n2, *n3;
FLOAT (*g)();
INT f, t;
{
   if(!IS_DIR(n1,t))
      NDS(n1,f) -= integr6(n1,n2,n3,g);
   if(!IS_DIR(n2,t))
      NDS(n2,f) -= integr6(n2,n3,n1,g);
   if(!IS_DIR(n3,t))
      NDS(n3,f) -= integr6(n3,n1,n2,g);
}

void sNeumann_bc(tGrid,g,f,t)
GRID *tGrid;
FLOAT (*g)();
INT f, t;
{
   ELEMENT *pel;
   
   for (pel = FDBE(tGrid); pel != NULL; pel = pel->succ){
      if (NOT_FF(pel->f[0]))
         sNeumann_bc_on_face(pel->n[1],pel->n[2],pel->n[3],g,f,t);
      if (NOT_FF(pel->f[1]))
         sNeumann_bc_on_face(pel->n[0],pel->n[2],pel->n[3],g,f,t);
      if (NOT_FF(pel->f[2]))
         sNeumann_bc_on_face(pel->n[0],pel->n[1],pel->n[3],g,f,t);
      if (NOT_FF(pel->f[3]))
         sNeumann_bc_on_face(pel->n[0],pel->n[1],pel->n[2],g,f,t);
   }
}

void sNeu_ij(n1,n2,n3,u1,u2,u3,Z,f,g,t)
NODE *n1, *n2, *n3;
FLOAT u1, u2, u3, (*g)();
INT Z, f, t;
{
   LINK *pli;
   
   COEFFS(n1,Z) += integr10(n1,n2,n3,u1,u2,u3,g);
   if (IS_DIR(n2,t))
      NDS(n1,f) -= integr9(n1,n2,n3,u1,u2,u3,g)*u2;
   else{
      for (pli = TSTART(n1); NBNODE(pli) != n2; pli = pli->next);
      COEFFLS(pli,Z) += integr9(n1,n2,n3,u1,u2,u3,g);
   }
   if (IS_DIR(n3,t))
      NDS(n1,f) -= integr9(n1,n3,n2,u1,u3,u2,g)*u3;
   else{
      for (pli = TSTART(n1); NBNODE(pli) != n3; pli = pli->next);
      COEFFLS(pli,Z) += integr9(n1,n3,n2,u1,u3,u2,g);
   }
}

void streat_Neumann_bc_on_face(n1,n2,n3,u,Z,f,g,t)
NODE *n1, *n2, *n3;
FLOAT (*g)();
INT u, Z, f, t;
{
   if(!IS_DIR(n1,t)) sNeu_ij(n1,n2,n3,NDS(n1,u),NDS(n2,u),NDS(n3,u),Z,f,g,t);
   if(!IS_DIR(n2,t)) sNeu_ij(n2,n3,n1,NDS(n2,u),NDS(n3,u),NDS(n1,u),Z,f,g,t);
   if(!IS_DIR(n3,t)) sNeu_ij(n3,n1,n2,NDS(n3,u),NDS(n1,u),NDS(n2,u),Z,f,g,t);
}

void streat_Neumann_bc(tGrid,u,Z,f,g,t)
GRID *tGrid;
FLOAT (*g)();
INT u, Z, f, t;
{
   ELEMENT *pel;
   
   for (pel = FDBE(tGrid); pel != NULL; pel = pel->succ){
      if (NOT_FF(pel->f[0]))
         streat_Neumann_bc_on_face(pel->n[1],pel->n[2],pel->n[3],u,Z,f,g,t);
      if (NOT_FF(pel->f[1]))
         streat_Neumann_bc_on_face(pel->n[0],pel->n[2],pel->n[3],u,Z,f,g,t);
      if (NOT_FF(pel->f[2]))
         streat_Neumann_bc_on_face(pel->n[0],pel->n[1],pel->n[3],u,Z,f,g,t);
      if (NOT_FF(pel->f[3]))
         streat_Neumann_bc_on_face(pel->n[0],pel->n[1],pel->n[2],u,Z,f,g,t);
   }
}

#else

void sstiff_matr(tGrid,nu,Z,f,u,t,t0)
GRID *tGrid; FLOAT nu, (*t0)(); INT Z, f, u, t;
{  eprintf("Error: sstiff_matr not available.\n");  }

#endif  /*  (N_DATA & SCALAR_NODE_DATA) && (N_DATA & ONE_NODE_MATR)  */

#endif  /* DIM == 3 */

/******************************************************************************/
/*                                                                            */
/*                                  buoyancy                                  */
/*                                                                            */
/******************************************************************************/

#if DIM == 3

#if (N_DATA & SCALAR_NODE_DATA) && (DATA_STR & LG_DATA)

void set_buo(n0,fa0,sum,dsum,n3,q1,q2,f,u)
NODE *n0;
FACE *fa0;
FLOAT sum, dsum, n3, q1, q2;
INT f, u;
{
   if (IS_FN(n0))
      ND(n0,f,2) += q1*(sum + NDS(n0,u));
   else if (n0->lgd){
      NDLG(n0,f,0) += q1*(sum + NDS(n0,u))*n0->lgd->t1[2];
      NDLG(n0,f,1) += q1*(sum + NDS(n0,u))*n0->lgd->t2[2];
   }
   if (IS_FF(fa0))
      FD(fa0,f) += q2*n3*(dsum - NDS(n0,u));
}

void buoyancy(tGrid,a,f,u) /*  a*\int_om \w_i\cdot\e_3 u\dx is added  */
GRID *tGrid;              /*                    to the velocity rhs  */
FLOAT a;
INT f, u;
{
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   FLOAT a1, a2, q1, q2, sum, dsum, nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM],detB;
   
   a1 = a/20.;
   a2 = a*FMULT/840.;
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = volume(n0->myvertex->x,n1->myvertex->x,
                    n2->myvertex->x,n3->myvertex->x);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      sum = NDS(n0,u) + NDS(n1,u) + NDS(n2,u) + NDS(n3,u);
      dsum = sum + sum;
      q1 = a1*detB;
      q2 = a2*detB;
      set_buo(n0,fa0,sum,dsum,nn0[2],q1,q2,f,u);
      set_buo(n1,fa1,sum,dsum,nn1[2],q1,q2,f,u);
      set_buo(n2,fa2,sum,dsum,nn2[2],q1,q2,f,u);
      set_buo(n3,fa3,sum,dsum,nn3[2],q1,q2,f,u);
    } 
}

#endif

/******************************************************************************/
/*                                                                            */
/*                              Marangoni effect                              */
/*                                                                            */
/******************************************************************************/

#if (N_DATA & SCALAR_NODE_DATA) && (DATA_STR & LG_DATA)

/* value of (\grad u)\cdot[w div m - m div w + (\grad w)m - (\grad m)w]      */
/* for w = h a, where h is a scalar function and and a is a constant vector  */
FLOAT Mar_val(gu,h,gh,a,m,divm,gm)
FLOAT gu[DIM], h, gh[DIM], a[DIM], m[DIM], divm, gm[DIM][DIM];
{
   FLOAT hdivm, agh, mgh;
   
   hdivm = h*divm;
   agh = DOT(a,gh);
   mgh = DOT(m,gh);
   return( gu[0]*( a[0]*hdivm - m[0]*agh + a[0]*mgh - h*DOT(gm[0],a) ) +
           gu[1]*( a[1]*hdivm - m[1]*agh + a[1]*mgh - h*DOT(gm[1],a) ) +
           gu[2]*( a[2]*hdivm - m[2]*agh + a[2]*mgh - h*DOT(gm[2],a) ) );
}

void Marij(n0,fan,nn,gu,f,b0,b1,b2,b3,m0,divm0,gm0,m1,divm1,gm1,m2,divm2,gm2,
                                                      m3,divm3,gm3,mc,divmc,gmc)
NODE *n0;
FACE *fan;
FLOAT m0[DIM], m1[DIM], m2[DIM], m3[DIM], mc[DIM], divm0, divm1, divm2, divm3, 
      divmc, gm0[DIM][DIM], gm1[DIM][DIM], gm2[DIM][DIM], gm3[DIM][DIM], 
      gmc[DIM][DIM], nn[DIM], gu[DIM], b0[DIM2], b1[DIM2], b2[DIM2], b3[DIM2];
INT f;
{
   FLOAT e1[DIM], e2[DIM], e3[DIM], gc[DIM];

   if (IS_FN(n0)){
      e1[0] = 1.;  e1[1] = 0.;  e1[2] = 0.;  
      e2[0] = 0.;  e2[1] = 1.;  e2[2] = 0.;  
      e3[0] = 0.;  e3[1] = 0.;  e3[2] = 1.;  
      ND(n0,f,0) -= Mar_val(gu,1.,b0,e1,m0,divm0,gm0) + 
                    Mar_val(gu,0.,b0,e1,m1,divm1,gm1) + 
                    Mar_val(gu,0.,b0,e1,m2,divm2,gm2) + 
                    Mar_val(gu,0.,b0,e1,m3,divm3,gm3) + 
                16.*Mar_val(gu,0.25,b0,e1,mc,divmc,gmc);
      ND(n0,f,1) -= Mar_val(gu,1.,b0,e2,m0,divm0,gm0) + 
                    Mar_val(gu,0.,b0,e2,m1,divm1,gm1) + 
                    Mar_val(gu,0.,b0,e2,m2,divm2,gm2) + 
                    Mar_val(gu,0.,b0,e2,m3,divm3,gm3) + 
                16.*Mar_val(gu,0.25,b0,e2,mc,divmc,gmc);
      ND(n0,f,2) -= Mar_val(gu,1.,b0,e3,m0,divm0,gm0) + 
                    Mar_val(gu,0.,b0,e3,m1,divm1,gm1) + 
                    Mar_val(gu,0.,b0,e3,m2,divm2,gm2) + 
                    Mar_val(gu,0.,b0,e3,m3,divm3,gm3) + 
                16.*Mar_val(gu,0.25,b0,e3,mc,divmc,gmc);
   }
   else if (n0->lgd){
      NDLG(n0,f,0) -= Mar_val(gu,1.,b0,n0->lgd->t1,m0,divm0,gm0) + 
                      Mar_val(gu,0.,b0,n0->lgd->t1,m1,divm1,gm1) + 
                      Mar_val(gu,0.,b0,n0->lgd->t1,m2,divm2,gm2) + 
                      Mar_val(gu,0.,b0,n0->lgd->t1,m3,divm3,gm3) + 
                  16.*Mar_val(gu,0.25,b0,n0->lgd->t1,mc,divmc,gmc);
      NDLG(n0,f,1) -= Mar_val(gu,1.,b0,n0->lgd->t2,m0,divm0,gm0) + 
                      Mar_val(gu,0.,b0,n0->lgd->t2,m1,divm1,gm1) + 
                      Mar_val(gu,0.,b0,n0->lgd->t2,m2,divm2,gm2) + 
                      Mar_val(gu,0.,b0,n0->lgd->t2,m3,divm3,gm3) + 
                  16.*Mar_val(gu,0.25,b0,n0->lgd->t2,mc,divmc,gmc);
/*      NDLG(n0,f,0) -= 20.*Mar_val(gu,0.25,b0,n0->lgd->t1,mc,divmc,gmc);
      NDLG(n0,f,1) -= 20.*Mar_val(gu,0.25,b0,n0->lgd->t2,mc,divmc,gmc);*/
   }
/*   if (IS_FF(fan)){
      gc[0] = b1[0] + b2[0] + b3[0];
      gc[1] = b1[1] + b2[1] + b3[1];
      gc[2] = b1[2] + b2[2] + b3[2];
      FD(fan,f) -= FMULT*Mar_val(gu,0.25,gc,nn,mc,divmc,gmc);
   }*/
}

void Marangoni_effect(tGrid,a,f,u)
GRID *tGrid;
FLOAT a;
INT f, u;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   ELEMENT *pelem;
   FLOAT m0[DIM], m1[DIM], m2[DIM], m3[DIM], mc[DIM], 
      divm0, divm1, divm2, divm3, divmc, 
      gm0[DIM][DIM], gm1[DIM][DIM], gm2[DIM][DIM], gm3[DIM][DIM], gmc[DIM][DIM],
      nn0[DIM], nn1[DIM], nn2[DIM], nn3[DIM], gu[DIM], xc[DIM], 
      b[DIM2][DIM2], detB;
  
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      FACES_OF_ELEMENT(fa0,fa1,fa2,fa3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      normal_vector(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa0,nn0);
      normal_vector(n0->myvertex->x,n2->myvertex->x,n3->myvertex->x,fa1,nn1);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n3->myvertex->x,fa2,nn2);
      normal_vector(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,fa3,nn3);
      POINT4(MYVERTEX(n0)->x,MYVERTEX(n1)->x,MYVERTEX(n2)->x,MYVERTEX(n3)->x,xc);
      values_of_m(MYVERTEX(n0)->x,m0,&divm0,gm0);
      values_of_m(MYVERTEX(n1)->x,m1,&divm1,gm1);
      values_of_m(MYVERTEX(n2)->x,m2,&divm2,gm2);
      values_of_m(MYVERTEX(n3)->x,m3,&divm3,gm3);
      values_of_m(xc,mc,&divmc,gmc);
      a *= detB/20.;
      gu[0] = a*(NDS(n0,u)*b[0][0] + NDS(n1,u)*b[1][0] + NDS(n2,u)*b[2][0]
                                                       + NDS(n3,u)*b[3][0]);
      gu[1] = a*(NDS(n0,u)*b[0][1] + NDS(n1,u)*b[1][1] + NDS(n2,u)*b[2][1]
                                                       + NDS(n3,u)*b[3][1]);
      gu[2] = a*(NDS(n0,u)*b[0][2] + NDS(n1,u)*b[1][2] + NDS(n2,u)*b[2][2]
                                                       + NDS(n3,u)*b[3][2]);
      Marij(n0,fa0,nn0,gu,f,b[0],b[1],b[2],b[3],m0,divm0,gm0,m1,divm1,gm1,m2,
                                           divm2,gm2,m3,divm3,gm3,mc,divmc,gmc);
      Marij(n1,fa1,nn1,gu,f,b[1],b[2],b[3],b[0],m1,divm1,gm1,m2,divm2,gm2,m3,
                                           divm3,gm3,m0,divm0,gm0,mc,divmc,gmc);
      Marij(n2,fa2,nn2,gu,f,b[2],b[3],b[0],b[1],m2,divm2,gm2,m3,divm3,gm3,m0,
                                           divm0,gm0,m1,divm1,gm1,mc,divmc,gmc);
      Marij(n3,fa3,nn3,gu,f,b[3],b[0],b[1],b[2],m3,divm3,gm3,m0,divm0,gm0,m1,
                                           divm1,gm1,m2,divm2,gm2,mc,divmc,gmc);
   } 
}

void marij(n0,gu,f,b0,mc,divmc,gmc)
NODE *n0;
FLOAT gu[DIM], b0[DIM2], mc[DIM], divmc, gmc[DIM][DIM];
INT f;
{
   if (n0->lgd){
      NDLG(n0,f,0) -= Mar_val(gu,0.25,b0,n0->lgd->t1,mc,divmc,gmc);
      NDLG(n0,f,1) -= Mar_val(gu,0.25,b0,n0->lgd->t2,mc,divmc,gmc);
   }
}

void marangoni_effect(tGrid,a,f,u)
GRID *tGrid;
FLOAT a;
INT f, u;
{
   NODE *n0, *n1, *n2, *n3;
   ELEMENT *pelem;
   FLOAT m0[DIM], m1[DIM], m2[DIM], m3[DIM], mc[DIM], divm, 
         gm[DIM][DIM], gu[DIM], b[DIM2][DIM2], detB;
  
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pelem);
      detB = barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,
                                     n2->myvertex->x,n3->myvertex->x,b);
      values_of_m(MYVERTEX(n0)->x,m0,&divm,gm);
      values_of_m(MYVERTEX(n1)->x,m1,&divm,gm);
      values_of_m(MYVERTEX(n2)->x,m2,&divm,gm);
      values_of_m(MYVERTEX(n3)->x,m3,&divm,gm);
      gm[0][0] = m0[0]*b[0][0] + m1[0]*b[1][0] + m2[0]*b[2][0] + m3[0]*b[3][0];
      gm[0][1] = m0[0]*b[0][1] + m1[0]*b[1][1] + m2[0]*b[2][1] + m3[0]*b[3][1];
      gm[0][2] = m0[0]*b[0][2] + m1[0]*b[1][2] + m2[0]*b[2][2] + m3[0]*b[3][2];
      gm[1][0] = m0[1]*b[0][0] + m1[1]*b[1][0] + m2[1]*b[2][0] + m3[1]*b[3][0];
      gm[1][1] = m0[1]*b[0][1] + m1[1]*b[1][1] + m2[1]*b[2][1] + m3[1]*b[3][1];
      gm[1][2] = m0[1]*b[0][2] + m1[1]*b[1][2] + m2[1]*b[2][2] + m3[1]*b[3][2];
      gm[2][0] = m0[2]*b[0][0] + m1[2]*b[1][0] + m2[2]*b[2][0] + m3[2]*b[3][0];
      gm[2][1] = m0[2]*b[0][1] + m1[2]*b[1][1] + m2[2]*b[2][1] + m3[2]*b[3][1];
      gm[2][2] = m0[2]*b[0][2] + m1[2]*b[1][2] + m2[2]*b[2][2] + m3[2]*b[3][2];
      divm = gm[0][0] + gm[1][1] + gm[2][2];
      POINT4(m0,m1,m2,m3,mc);
      a *= detB;
      gu[0] = a*(NDS(n0,u)*b[0][0] + NDS(n1,u)*b[1][0] + NDS(n2,u)*b[2][0]
                                                       + NDS(n3,u)*b[3][0]);
      gu[1] = a*(NDS(n0,u)*b[0][1] + NDS(n1,u)*b[1][1] + NDS(n2,u)*b[2][1]
                                                       + NDS(n3,u)*b[3][1]);
      gu[2] = a*(NDS(n0,u)*b[0][2] + NDS(n1,u)*b[1][2] + NDS(n2,u)*b[2][2]
                                                       + NDS(n3,u)*b[3][2]);
      marij(n0,gu,f,b[0],mc,divm,gm);
      marij(n1,gu,f,b[1],mc,divm,gm); 
      marij(n2,gu,f,b[2],mc,divm,gm);
      marij(n3,gu,f,b[3],mc,divm,gm); 
   }  
}

#endif

#endif  /* DIM == 3 */

