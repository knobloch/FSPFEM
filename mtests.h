/******************************************************************************/
/*                                                                            */
/*                                   tests                                    */
/*                                                                            */
/******************************************************************************/

#define N1 1262
#define N2 1262 
#define N3   40 
#define N4    5 

void get_son_faces();
void get_sons();
NODE *unif_middle();

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA)

void p1nc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,
                        uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid;
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], 
      pp[N2], gg[N2];
INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], 
    *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   INT i, j, k;
   ELEMENT *pel;
   FACE *pface;
   FLINK *pfl;

   i = 0;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ)
      pface->index = i++;
   *n = i;
   i = 0;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      pel->status = i++;
   *m = i;
   printf("N1 = %i,  N2 = %i,  n (velocity) = %i, m (pressure) = %i\n",N1,N2,*n,*m);
   for (i = 0; i < *n; i++)
      for(j = 0; j < *n; j++){
         aa[i][j] = 0.;
         pattern[i][j] = 0;
      }
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ){
      i = pface->index;
      aa[i][i] = COEFF_FF(pface,ZA);
      pattern[i][i] = 1;
      uu0[i] = FDV(pface,u,0);
      uu1[i] = FDV(pface,u,1);
      ff0[i] = FDV(pface,f,0);
      ff1[i] = FDV(pface,f,1);
      for (pfl = pface->fstart; pfl != NULL; pfl = pfl->next){
	  j = pfl->nbface->index;
	  aa[i][j] = COEFF_FL(pfl,ZA);
          pattern[i][j] = 1;
      }
   }
   for (i = 0; i < *m; i++)
      for (j = 0; j < *n; j++)
         bb0[i][j] = bb1[i][j] = 0.;
   *vanka_np = 1;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
     i = pel->status;
     gg[i] = ED(pel,g);
     pp[i] = ED(pel,p);
     vanka_p[i] = i;
     vanka_nu[i] = 0;
     vanka_u[i][1] = vanka_u[i][2] = -1;
     for (j=0; j < 3; j++)
        if (IS_FF(pel->f[j])){
           bb0[i][pel->f[j]->index] = COEFF_BDF(pel,ZB,j,0);
           bb1[i][pel->f[j]->index] = COEFF_BDF(pel,ZB,j,1);
           vanka_u[i][(vanka_nu[i])++] = pel->f[j]->index;
        }
   }
}

#else

void p1nc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid; FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1], pp[N2], gg[N2]; INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{  eprintf("Error: p1nc_structures_to_matrices not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (N_DATA & NODE_IFLAG) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) 

void p2c_p1c_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,
                        uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid;
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], 
      pp[N2], gg[N2];
INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], 
    *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   INT i, j, k;
   NODE *pnode;
   FACE *pface;
   LINK *pl;
   NFLINK *pnf;
   FLINK *pfl;
   FNLINK *pfnl;

   i = 0;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
      pnode->index = i++;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ)
      pface->index = i++;
   *n = i;
   i = 0;
   for (pnode = FIRSTN(tGrid); pnode; pnode=pnode->succ)
      pnode->iflag = i++;
   *m = i;
   printf("N1 = %i,  N2 = %i,  n (velocity) = %i, m (pressure) = %i\n",N1,N2,*n,*m);
   for (i = 0; i < *n; i++)
      for(j = 0; j < *n; j++){
         aa[i][j] = 0.;
         pattern[i][j] = 0;
      }
   for (pnode=FIRSTNODE(tGrid); pnode; pnode=SUCC(pnode)){
      i = pnode->index;
      aa[i][i] = COEFFN(pnode,ZA);
      pattern[i][i] = 1;
      uu0[i] = ND(pnode,u,0);
      uu1[i] = ND(pnode,u,1);
      ff0[i] = ND(pnode,f,0);
      ff1[i] = ND(pnode,f,1);
      for (pl=START(pnode); pl; pl=NEXT(pl)){
         j = NBNODE(pl)->index;
         aa[i][j] = COEFFL(pl,ZA);
         pattern[i][j] = 1;
      }
      for (pnf=NFSTART(pnode); pnf; pnf=NEXT(pnf)){
         j = NBFACE(pnf)->index;
         aa[i][j] = COEFFL(pnf,ZA);
         pattern[i][j] = 1;
      }
   }
   for (pface=FIRSTFACE(tGrid); pface; pface=SUCC(pface)){
      i = pface->index;
      aa[i][i] = COEFF_FF(pface,ZA);
      pattern[i][i] = 1;
      uu0[i] = FDV(pface,u,0);
      uu1[i] = FDV(pface,u,1);
      ff0[i] = FDV(pface,f,0);
      ff1[i] = FDV(pface,f,1);
      for (pfl=FSTART(pface); pfl; pfl=NEXT(pfl)){
         j = NBFACE(pfl)->index;
         aa[i][j] = COEFF_FL(pfl,ZA);
         pattern[i][j] = 1;
      }
      for (pfnl=FNSTART(pface); pfnl; pfnl=NEXT(pfnl)){
         j = NBNODE(pfnl)->index;
         aa[i][j] = COEFF_FL(pfnl,ZA);
         pattern[i][j] = 1;
      }
   }
   for (i = 0; i < *m; i++)
      for (j = 0; j < *n; j++)
         bb0[i][j] = bb1[i][j] = 0.;
   *vanka_np = 1;
   for (pnode = FDBN(tGrid); pnode; pnode = pnode->succ){
      i = pnode->iflag;
      gg[i] = NDS(pnode,g);
      pp[i] = NDS(pnode,p);
      vanka_p[i] = i;
      vanka_nu[i] = 0;
      if (IS_FN(pnode)){
         j = pnode->index;
         bb0[i][j] = COEFFB(pnode,ZB,0);
         bb1[i][j] = COEFFB(pnode,ZB,1);
         vanka_u[i][(vanka_nu[i])++] = j;
      }
      for (pl = pnode->start; pl; pl = pl->next){
         j = NBNODE(pl)->index;
         bb0[i][j] = COEFFBL(pl,ZB,0);
         bb1[i][j] = COEFFBL(pl,ZB,1);
         vanka_u[i][(vanka_nu[i])++] = j;
      }
      for (pnf = pnode->nfstart; pnf; pnf = pnf->next){
         j = NBFACE(pnf)->index;
         bb0[i][j] = COEFF_NF(pnf,ZB,0);
         bb1[i][j] = COEFF_NF(pnf,ZB,1);
         vanka_u[i][(vanka_nu[i])++] = j;
      }
   }
}

#else

void p2c_p1c_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid; FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1], pp[N2], gg[N2]; INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{  eprintf("Error: p2c_p1c_structures_to_matrices not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (E_DATA & VECTOR_ELEMENT_DATA)

void p2c_elbub_p1_disc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,
pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid;
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1],
      pp[N2], gg[N2];
INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], 
    *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   LINK *pl;
   NFLINK *pnf;
   FLINK *pfl;
   FNLINK *pfnl;
   INT i, j, k, l;

   i = 0;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
      pnode->index = i++;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ)
      pface->index = i++;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      pel->eflag = i++;
   *n = i;
   i = 0;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      pel->status = i;
      i += 3;
   }
   *m = i;
   printf("N1 = %i,  N2 = %i,  n (velocity) = %i, m (pressure) = %i\n",N1,N2,*n,*m);
   for (i = 0; i < *n; i++)
      for(j = 0; j < *n; j++){
         aa[i][j] = 0.;
         pattern[i][j] = 0;
      }
   for (pnode=FIRSTNODE(tGrid); pnode; pnode=SUCC(pnode)){
      i = pnode->index;
      aa[i][i] = COEFFN(pnode,ZA);
      pattern[i][i] = 1;
      uu0[i] = ND(pnode,u,0);
      uu1[i] = ND(pnode,u,1);
      ff0[i] = ND(pnode,f,0);
      ff1[i] = ND(pnode,f,1);
      for (pl=START(pnode); pl; pl=NEXT(pl)){
         j = NBNODE(pl)->index;
         aa[i][j] = COEFFL(pl,ZA);
         pattern[i][j] = 1;
      }
      for (pnf=NFSTART(pnode); pnf; pnf=NEXT(pnf)){
         j = NBFACE(pnf)->index;
         aa[i][j] = COEFFL(pnf,ZA);
         pattern[i][j] = 1;
      }
   }
   for (pface=FIRSTFACE(tGrid); pface; pface=SUCC(pface)){
      i = pface->index;
      aa[i][i] = COEFF_FF(pface,ZA);
      pattern[i][i] = 1;
      uu0[i] = FDV(pface,u,0);
      uu1[i] = FDV(pface,u,1);
      ff0[i] = FDV(pface,f,0);
      ff1[i] = FDV(pface,f,1);
      for (pfl=FSTART(pface); pfl; pfl=NEXT(pfl)){
         j = NBFACE(pfl)->index;
         aa[i][j] = COEFF_FL(pfl,ZA);
         pattern[i][j] = 1;
      }
      for (pfnl=FNSTART(pface); pfnl; pfnl=NEXT(pfnl)){
         j = NBNODE(pfnl)->index;
         aa[i][j] = COEFF_FL(pfnl,ZA);
         pattern[i][j] = 1;
      }
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      i = pel->eflag;
      aa[i][i] = COEFF_EE(pel,ZA);
      pattern[i][i] = 1;
      uu0[i] = EDV(pel,u,0);
      uu1[i] = EDV(pel,u,1);
      ff0[i] = EDV(pel,f,0);
      ff1[i] = EDV(pel,f,1);
      if (IS_FN(pel->n[0])){
         aa[i][pel->n[0]->index] = COEFF_EN(pel,ZA,0); 
         aa[pel->n[0]->index][i] = COEFF_NE(pel,ZA,0);
         pattern[i][pel->n[0]->index] = 1;
         pattern[pel->n[0]->index][i] = 1;
      }
      if (IS_FN(pel->n[1])){
         aa[i][pel->n[1]->index] = COEFF_EN(pel,ZA,1); 
         aa[pel->n[1]->index][i] = COEFF_NE(pel,ZA,1);
         pattern[i][pel->n[1]->index] = 1;
         pattern[pel->n[1]->index][i] = 1;
      }
      if (IS_FN(pel->n[2])){
         aa[i][pel->n[2]->index] = COEFF_EN(pel,ZA,2); 
         aa[pel->n[2]->index][i] = COEFF_NE(pel,ZA,2);
         pattern[i][pel->n[2]->index] = 1;
         pattern[pel->n[2]->index][i] = 1;
      }
      if (IS_FF(pel->f[0])){
         aa[i][pel->f[0]->index] = COEFF_EF(pel,ZA,0); 
         aa[pel->f[0]->index][i] = COEFF_FE(pel,ZA,0);
         pattern[i][pel->f[0]->index] = 1;
         pattern[pel->f[0]->index][i] = 1;
      }
      if (IS_FF(pel->f[1])){
         aa[i][pel->f[1]->index] = COEFF_EF(pel,ZA,1); 
         aa[pel->f[1]->index][i] = COEFF_FE(pel,ZA,1);
         pattern[i][pel->f[1]->index] = 1;
         pattern[pel->f[1]->index][i] = 1;
      }
      if (IS_FF(pel->f[2])){
         aa[i][pel->f[2]->index] = COEFF_EF(pel,ZA,2); 
         aa[pel->f[2]->index][i] = COEFF_FE(pel,ZA,2);
         pattern[i][pel->f[2]->index] = 1;
         pattern[pel->f[2]->index][i] = 1;
      }
   }
   for (i = 0; i < *m; i++)
      for (j = 0; j < *n; j++)
         bb0[i][j] = bb1[i][j] = 0.;
   *vanka_np = 3;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      i = pel->status;
      gg[i]   = EDSN(pel,g,0);
      gg[i+1] = EDSN(pel,g,1);
      gg[i+2] = EDSN(pel,g,2);
      pp[i]   = EDSN(pel,p,0);
      pp[i+1] = EDSN(pel,p,1);
      pp[i+2] = EDSN(pel,p,2);
      for (k = 0; k < NVERT; k++){
         vanka_p[i+k] = i+k;
         bb0[i+k][pel->eflag] = COEFF_BDF(pel,ZB,k,0);
         bb1[i+k][pel->eflag] = COEFF_BDF(pel,ZB,k,1);
         vanka_u[i+k][0]= pel->eflag;
         vanka_nu[i+k] = 1;
      }
      for (l = 0; l < NVERT; l++)
         if (IS_FN(pel->n[l])){
            j = pel->n[l]->index;
            for (k = 0; k < NVERT; k++){
               bb0[i+k][j] = COEFF_BFDN(pel,ZB,k,l,0);
               bb1[i+k][j] = COEFF_BFDN(pel,ZB,k,l,1);
               vanka_u[i+k][vanka_nu[i]] = j;
            }
            (vanka_nu[i])++;
         }
      for (l = 0; l < SIDES; l++)
         if (IS_FF(pel->f[l])){
            j = pel->f[l]->index;
            for (k = 0; k < NVERT; k++){
               bb0[i+k][j] = COEFF_BFDF(pel,ZB,k,l,0);
               bb1[i+k][j] = COEFF_BFDF(pel,ZB,k,l,1);
               vanka_u[i+k][vanka_nu[i]] = j;
            }
            (vanka_nu[i])++;
         }
   }
}

#else

void p2c_elbub_p1_disc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid; FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1], pp[N2], gg[N2]; INT ZA, ZB, u, p, f, g, *n, *m, pattern[N1][N1], *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{  eprintf("Error: p2c_elbub_p1_disc_structures_to_matrices not available.\n");  }

#endif

void structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,u_space,p_space,n,m,aa,
pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u)
GRID *tGrid;
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1],
      pp[N2], gg[N2];
INT ZA, ZB, u, p, f, g, u_space, p_space, *n, *m, pattern[N1][N1], 
    *vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   if (u_space == P1_NC && p_space == P0)
      p1nc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,
                       uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u);
   else if (u_space == P2C && p_space == P1C)
      p2c_p1c_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,pattern,bb0,bb1,
                       uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u);
   else if (u_space == P2C_ELBUB && p_space == P1_DISC)
      p2c_elbub_p1_disc_structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,n,m,aa,
       pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,vanka_np,vanka_p,vanka_nu,vanka_u);
   else
      eprintf("Error: structures_to_matrices not available.\n");
}

void p_p0_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      pr[pel->sons[0]->status][pel->status] =
      pr[pel->sons[1]->status][pel->status] =
      pr[pel->sons[2]->status][pel->status] =
      pr[pel->sons[3]->status][pel->status] = 1.;
}

#if N_DATA & NODE_IFLAG

void p_p1c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   NODE *n1,  *n2,  *n3,  *n12, *n13, *n23, *n11, *n22, *n33;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      n11 = (n1)->son;
      n22 = (n2)->son;
      n33 = (n3)->son;
      n12 = unif_middle(n11,n22);
      n13 = unif_middle(n11,n33);
      n23 = unif_middle(n22,n33);
      pr[n11->iflag][n1->iflag] = 1.;
      pr[n12->iflag][n1->iflag] = 0.5;
      pr[n13->iflag][n1->iflag] = 0.5;
      pr[n22->iflag][n2->iflag] = 1.;
      pr[n12->iflag][n2->iflag] = 0.5;
      pr[n23->iflag][n2->iflag] = 0.5;
      pr[n33->iflag][n3->iflag] = 1.;
      pr[n13->iflag][n3->iflag] = 0.5;
      pr[n23->iflag][n3->iflag] = 0.5;
   }
}

#else

void p_p1c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid; FLOAT pr[N1][N1]; INT nc, nf;
{  eprintf("Error: p_p1c_prolong_matr not available.\n");  }

#endif

INT face_index(pel,f)
ELEMENT *pel;
FACE *f;
{
   if (f == pel->f[0])
      return(0);
   else if (f == pel->f[1])
      return(1);
   else if (f == pel->f[2])
      return(2);
   else
      return(-1);
}

void p1disc_prolong_coeff(pel,f1,f12,f13,f11,pr)
ELEMENT *pel;
FACE *f1, *f12, *f13, *f11;
FLOAT pr[N1][N1];
{
   INT j, i1=face_index(pel,f1), i2, i3, i00, i11, i22, i33;

   j   = pel->status + i1;
   i00 =  face_index(pel->sons[ 3],f11);
   i11 =  face_index(pel->sons[i1],f11);
   i2  = (i1+1)%3;
   if ((i22=face_index(pel->sons[i2],f12)) < 0){
      i3 = i2;
      i2 = (i1+2)%3;
      i22 = face_index(pel->sons[i2],f12);
   }
   else
      i3 = (i1+2)%3;
   i33 = face_index(pel->sons[i3],f13);
   pr[pel->sons[ 3]->status][j] = pr[pel->sons[ 3]->status+1][j] 
                                = pr[pel->sons[ 3]->status+2][j] =  0.5;
   pr[pel->sons[i1]->status][j] = pr[pel->sons[i1]->status+1][j] 
                                = pr[pel->sons[i1]->status+2][j] = -0.5;
   pr[pel->sons[i2]->status][j] = pr[pel->sons[i2]->status+1][j] 
                                = pr[pel->sons[i2]->status+2][j] =  0.5;
   pr[pel->sons[i3]->status][j] = pr[pel->sons[i3]->status+1][j] 
                                = pr[pel->sons[i3]->status+2][j] =  0.5;
   pr[pel->sons[i2]->status+i22][j] = pr[pel->sons[i3]->status+i33][j] = 1.;
   pr[pel->sons[i1]->status+i11][j] = pr[pel->sons[ 3]->status+i00][j] = 0.;
}

void p_p1disc_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p1disc_prolong_coeff(pel,pel->f[0],f12,f13,f11,pr);
      p1disc_prolong_coeff(pel,pel->f[1],f21,f23,f22,pr);
      p1disc_prolong_coeff(pel,pel->f[2],f31,f32,f33,pr);
   }
}

INT number_of_neigh_el(pnode)
NODE *pnode;
{
   LINK *pli;
   INT n=0;

   for (pli=TSTART(pnode); pli; pli=pli->next, n++);
   if (IS_BN(pnode))
      n--;
   return(n);
}

#if N_DATA & NODE_IFLAG

void p_p0_to_p1c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = pel->status;
      pr[pel->n[0]->iflag][j] = 1./number_of_neigh_el(pel->n[0]);
      pr[pel->n[1]->iflag][j] = 1./number_of_neigh_el(pel->n[1]);
      pr[pel->n[2]->iflag][j] = 1./number_of_neigh_el(pel->n[2]);
   }
}

#else

void p_p0_to_p1c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid; FLOAT pr[N1][N1]; INT nc, nf;
{  eprintf("Error: p_p0_to_p1c_prolong_matr not available.\n");  }

#endif

void p_p0_to_p1disc_prolong_matr(pr,nf,nc)
FLOAT pr[N1][N1];
INT nc, nf;
{
   INT i=0, j, k;

   for (j=0; j < nc; j++, i+=3){
      for (k=0; k < nf; k++)
         pr[k][j] = 0.;
      pr[i][j] = pr[i+1][j] = pr[i+2][j] = 1.;
   }
}

void p1nc_prolong_coeff(f1,f12,f13,f22,f33,f23,f32,f21,f31,pr)
FACE *f1, *f12, *f13, *f22, *f33, *f23, *f32, *f21, *f31;
FLOAT pr[N1][N1];
{
   if (IS_FF(f1)){
      pr[f12->index][f1->index] = 1.;
      pr[f13->index][f1->index] = 1.;
      pr[f22->index][f1->index] = 0.5;
      pr[f33->index][f1->index] = 0.5;
      if (IS_FF(f23)){
         if (IS_BF(f23)){
            pr[f23->index][f1->index] =  0.5;
            pr[f21->index][f1->index] = -0.5;
         }
         else{
            pr[f23->index][f1->index] =  0.25;
            pr[f21->index][f1->index] = -0.25;
         }
      }
      if (IS_FF(f32)){
         if (IS_BF(f32)){
            pr[f32->index][f1->index] =  0.5;
            pr[f31->index][f1->index] = -0.5;
         }
         else{
            pr[f32->index][f1->index] =  0.25;
            pr[f31->index][f1->index] = -0.25;
         }
      }
   }
}

void u_p1nc_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      get_son_faces(pel,&f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p1nc_prolong_coeff(pel->f[0],f12,f13,f22,f33,f23,f32,f21,f31,pr);
      p1nc_prolong_coeff(pel->f[1],f21,f23,f11,f33,f31,f13,f32,f12,pr);
      p1nc_prolong_coeff(pel->f[2],f31,f32,f11,f22,f21,f12,f23,f13,pr);
   }
}

void p2c_prolong_coeff(n1,n11,n12,n13,n23,f1,f12,f13,f11,pr)
NODE *n1, *n11, *n12, *n13, *n23;
FACE *f1, *f12, *f13, *f11;
FLOAT pr[N1][N1];
{
   if (IS_FN(n1)){
      pr[n11->index][n1->index] = 1.;
      pr[n12->index][n1->index] = 0.5;
      pr[n13->index][n1->index] = 0.5;
   }
   if (IS_FF(f1)){
      pr[n23->index][f1->index] = 0.25;
      pr[f12->index][f1->index] = 0.25;
      pr[f13->index][f1->index] = 0.25;
      pr[f11->index][f1->index] = 0.25;
   }
}

void u_p2c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   NODE *n1,  *n2,  *n3,  *n12, *n13, *n23, *n11, *n22, *n33;
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   INT i, j;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      get_sons(pel,&n1,&n2,&n3,&n12,&n13,&n23,&n11,&n22,&n33,
                   &f11,&f12,&f13,&f21,&f22,&f23,&f31,&f32,&f33);
      p2c_prolong_coeff(n1,n11,n12,n13,n23,pel->f[0],f12,f13,f11,pr);
      p2c_prolong_coeff(n2,n22,n12,n23,n13,pel->f[1],f21,f23,f22,pr);
      p2c_prolong_coeff(n3,n33,n13,n23,n12,pel->f[2],f31,f32,f33,pr);
   }
}

void u_p2c_elbub_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   INT j;

   u_p2c_prolong_matr(tGrid,pr,nf,nc);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = pel->eflag;
      pr[pel->sons[3]->f[0]->index][j] = pr[pel->sons[3]->f[1]->index][j] = 
      pr[pel->sons[3]->f[2]->index][j] = pr[pel->sons[0]->eflag][j] = 
      pr[pel->sons[1]->eflag][j] = pr[pel->sons[2]->eflag][j] = 0.125;
      pr[pel->sons[3]->eflag][j] = -0.125;
   }
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void set_node_values(pface,n0,n1,n2)
FACE *pface;
NODE *n0, *n1, *n2;
{
   FLOAT a=2.;

   if (IS_BF(pface))
      a = 1.; 
   if (IS_FN(n0))
      ND(n0,D,0) = -1./number_of_neigh_el(n0);
   if (IS_FN(n1))
      ND(n1,D,0) = a/number_of_neigh_el(n1);
   if (IS_FN(n2))
      ND(n2,D,0) = a/number_of_neigh_el(n2);
}

void u_p1nc_to_p2c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface, *pf;
   INT i, j, nn, jj;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   i = 0;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
      pnode->index = i++;
   nn = i;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ)
      pface->index = i++;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ){
      jj = pface->index-nn;
      for (pnode = FIRSTN(tGrid); pnode; pnode=pnode->succ)
         ND(pnode,D,0) = 0.;
      for (pf = FIRSTF(tGrid); pf; pf=pf->succ)
         FDV(pf,D,0) = 0.;
      FDV(pface,D,0) = 1.;
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         if (pel->f[0] == pface)
            set_node_values(pface,pel->n[0],pel->n[1],pel->n[2]);
         if (pel->f[1] == pface)
            set_node_values(pface,pel->n[1],pel->n[2],pel->n[0]);
         if (pel->f[2] == pface)
            set_node_values(pface,pel->n[2],pel->n[0],pel->n[1]);
      }
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         if (IS_FN(pel->n[0]))
            pr[pel->n[0]->index][jj] = ND(pel->n[0],D,0);
         if (IS_FN(pel->n[1]))
            pr[pel->n[1]->index][jj] = ND(pel->n[1],D,0);
         if (IS_FN(pel->n[2]))
            pr[pel->n[2]->index][jj] = ND(pel->n[2],D,0);
         if (IS_FF(pel->f[0]))
            pr[pel->f[0]->index][jj] = 4.*FDV(pel->f[0],D,0)
                                     - 2.*(ND(pel->n[1],D,0)+ND(pel->n[2],D,0));
         if (IS_FF(pel->f[1]))
            pr[pel->f[1]->index][jj] = 4.*FDV(pel->f[1],D,0)
                                     - 2.*(ND(pel->n[0],D,0)+ND(pel->n[2],D,0));
         if (IS_FF(pel->f[2]))
            pr[pel->f[2]->index][jj] = 4.*FDV(pel->f[2],D,0)
                                     - 2.*(ND(pel->n[1],D,0)+ND(pel->n[0],D,0));
      }
   }
}

void u_p1nc_to_p2cb_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid;
FLOAT pr[N1][N1];
INT nc, nf;
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface, *pf;
   FLOAT a0, a1, a2, s0, s1, s2;
   INT i, j, nn, jj;

   for (i=0; i < nf; i++)
      for (j=0; j < nc; j++)
         pr[i][j] = 0.;
   i = 0;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
      pnode->index = i++;
   nn = i;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ)
      pface->index = i++;
   for (pface = FIRSTFACE(tGrid); pface; pface=pface->succ){
      jj = pface->index-nn;
      for (pnode = FIRSTN(tGrid); pnode; pnode=pnode->succ)
         ND(pnode,D,0) = 0.;
      for (pf = FIRSTF(tGrid); pf; pf=pf->succ)
         FDV(pf,D,0) = 0.;
      FDV(pface,D,0) = 1.;
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         if (pel->f[0] == pface)
            set_node_values(pface,pel->n[0],pel->n[1],pel->n[2]);
         if (pel->f[1] == pface)
            set_node_values(pface,pel->n[1],pel->n[2],pel->n[0]);
         if (pel->f[2] == pface)
            set_node_values(pface,pel->n[2],pel->n[0],pel->n[1]);
      }
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         s0 = s1 = s2 = a0 = a1 = a2 = 0.;
         if (IS_FN(pel->n[0]))
            s0 = pr[pel->n[0]->index][jj] = ND(pel->n[0],D,0);
         if (IS_FN(pel->n[1]))
            s1 = pr[pel->n[1]->index][jj] = ND(pel->n[1],D,0);
         if (IS_FN(pel->n[2]))
            s2 = pr[pel->n[2]->index][jj] = ND(pel->n[2],D,0);
         if (IS_FF(pel->f[0]))
            a0 = pr[pel->f[0]->index][jj] = 4.*FDV(pel->f[0],D,0)
                                     - 2.*(ND(pel->n[1],D,0)+ND(pel->n[2],D,0));
         if (IS_FF(pel->f[1]))
            a1 = pr[pel->f[1]->index][jj] = 4.*FDV(pel->f[1],D,0)
                                     - 2.*(ND(pel->n[0],D,0)+ND(pel->n[2],D,0));
         if (IS_FF(pel->f[2]))
            a2 = pr[pel->f[2]->index][jj] = 4.*FDV(pel->f[2],D,0)
                                     - 2.*(ND(pel->n[1],D,0)+ND(pel->n[0],D,0));
         if (pel->f[0] == pface || pel->f[1] == pface || pel->f[2] == pface)
            pr[pel->eflag][jj] = 9.-9.*(s0+s1+s2) - 3.*(a0+a1+a2);
         else
            pr[pel->eflag][jj] = -9.*(s0+s1+s2) - 3.*(a0+a1+a2);
      }
   }
}

#else

void u_p1nc_to_p2c_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid; FLOAT pr[N1][N1]; INT nc, nf;
{  eprintf("Error: u_p1nc_to_p2c_prolong_matr not available.\n");  }

void u_p1nc_to_p2cb_prolong_matr(tGrid,pr,nf,nc)
GRID *tGrid; FLOAT pr[N1][N1]; INT nc, nf;
{  eprintf("Error: u_p1nc_to_p2cb_prolong_matr not available.\n");  }

#endif

void u_prolong_matr(tGrid,pr,nf,nc,u_space)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf, u_space;
{
   if (u_space == P1_NC)
      u_p1nc_prolong_matr(tGrid,pr,nf,nc);
   else if (u_space == P2C)
      u_p2c_prolong_matr(tGrid,pr,nf,nc);
   else if (u_space == P2C_ELBUB)
      u_p2c_elbub_prolong_matr(tGrid,pr,nf,nc);
   else
      eprintf("Error: u_prolong_matr not available.\n");
}

void p_prolong_matr(tGrid,pr,nf,nc,p_space)
GRID *tGrid;  /*  coarse grid  */
FLOAT pr[N1][N1];
INT nc, nf, p_space;
{
   if (p_space == P0)
      p_p0_prolong_matr(tGrid,pr,nf,nc,p_space);
   else if (p_space == P1C)
      p_p1c_prolong_matr(tGrid,pr,nf,nc);
   else if (p_space == P1_DISC)
      p_p1disc_prolong_matr(tGrid,pr,nf,nc);
   else
      eprintf("Error: p_prolong_matr not available.\n");
}

void m_multA(aa,xx,bb,n1,n2) /*  aa  is   n1 x n2  */
FLOAT aa[N2][N2], *xx, *bb;
INT n1, n2;
{
   INT i, j;

   for (i=0; i < n1; i++){
      bb[i] = 0.;
      for (j=0; j < n2; j++)
         bb[i] += aa[i][j]*xx[j];
   }
}

void m_multAT(aa,xx,bb,n1,n2) /*  aa  is   n1 x n2  */
FLOAT aa[N2][N2], *xx, *bb;
INT n1, n2;
{
   INT i, j;

   for (i=0; i < n2; i++){
      bb[i] = 0.;
      for (j=0; j < n1; j++)
         bb[i] += aa[j][i]*xx[j];
   }
}

void m_set_value(r,zz,n)
FLOAT r, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = r;
}

FLOAT m_max_abs_value(zz,n)
FLOAT *zz;
INT n;
{
    INT i;
    FLOAT max=0.;

    for (i=0; i < n; i++)
       max = MAX(max,fabs(zz[i]));
    return(max);
}

FLOAT m_dot(xx,yy,n)
FLOAT *xx, *yy;
INT n;
{
    FLOAT sum=0.;
    INT i;

    for (i=0; i < n; i++)
       sum += xx[i]*yy[i];
    return(sum);
}

void m_copy(xx,zz,n)
FLOAT *xx, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = xx[i];
}

void m_inv(xx,zz,n)
FLOAT *xx, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = -xx[i];
}

void m_add(xx,yy,zz,n)
FLOAT *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = xx[i] + yy[i];
}

void m_subtr(xx,yy,zz,n)
FLOAT *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = xx[i] - yy[i];
}

void m_mult_and_add(r,xx,yy,zz,n)
FLOAT r, *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = r*xx[i] + yy[i];
}

void m_mult_and_subtr(r,xx,yy,zz,n)
FLOAT r, *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = r*xx[i] - yy[i];
}

void m_damp(r,xx,yy,zz,n)
FLOAT r, *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = xx[i] + r*(yy[i] - xx[i]);
}

void m_multiply(xx,yy,zz,n)
FLOAT *xx, *yy, *zz;
INT n;
{
    INT i;

    for (i=0; i < n; i++)
       zz[i] = xx[i]*yy[i];
}

void m_copy_mat_to_mat(aa,hh,n)
FLOAT aa[N2][N2], hh[N2][N2];
INT n;
{
   INT i, j;

   for (i=0; i < n; i++)
      for (j=0; j < n; j++)
         hh[i][j] = aa[i][j];
}

void m_ILU(aa,hh,pattern,n)
FLOAT aa[N2][N2], hh[N2][N2];
INT pattern[N2][N2];
INT n;
{
   INT r, i, j;
   FLOAT pivot;
   
   m_copy_mat_to_mat(aa,hh,n);

   for (r = 0; r < n; r++)
      for (i = 0; i < n; i++)
         if ( pattern[r][i] && i > r ){
           pivot = hh[i][r] /= hh[r][r];
           for (j = 0; j < n; j++)
             if ( pattern[r][j] && j > r && pattern[i][j])
                hh[i][j] -= pivot*hh[r][j];
         }         
}   

/* Matrix h contains an ILU-decomposition; LU x = y is solved.  */           
void m_solveILU(hh,xx,yy,n)
FLOAT hh[N2][N2], xx[N1], yy[N1];
INT n;
{
   INT i, j;
   FLOAT sumf;
   
   /* find x: L x = y */
   for (i = 0; i < n; i++){
      sumf = 0.0;
      for (j = 0; j < n; j++)
         if ( j < i )
            sumf += hh[i][j]*xx[j];
      xx[i] = yy[i] - sumf;
   }
   
   /* find solution of U . = x and store it in x */
   for (i = n-1; i >= 0; i--){
         sumf = 0.0;
         for (j = 0; j < n; j++)
            if ( j > i )
               sumf += hh[i][j]*xx[j];
         xx[i] = (xx[i] - sumf)/hh[i][i];
   }
}

INT tm_CG(aa,bb,pp,n,n2)
FLOAT aa[N1][N1], bb[N1], pp[N1];
INT n;
{
   FLOAT ss[N1], gg[N1], ww[N1], qq[N1], uu[DIM][N1], uue[DIM][N1], 
         abs_g2=0., gamma, ro, sum, eps1;
   INT i=0, j;

   for (j=0; j<n; j++){
      pp[j] = 0.;
      gg[j] = -bb[j];
      ww[j] = gg[j];
      abs_g2 += gg[j]*gg[j];
   }
   eps1 = abs_g2*EPS;
   while (abs_g2>eps1){
      printf("%e\n",abs_g2);
   /* printf("%i:  %e\n",i,abs_g2); */
      m_multAT(aa,ww,ss,n,n2);
      m_multA(aa,ss,qq,n,n2);
      ro = abs_g2/m_dot(qq,ww,n);
      sum = 0.;
      for (j = 0; j<n; j++){
         pp[j] -= ro*ww[j];
         gg[j] -= ro*qq[j];
         sum += gg[j]*gg[j];
      }
      gamma = sum/abs_g2;
      abs_g2 = sum;
      m_mult_and_add(gamma,ww,gg,ww,n);
      i++;
   }
   return(i);
}

INT m_CG(aa,bb,pp,n)
FLOAT aa[N1][N1], bb[N1], pp[N1];
INT n;
{
   FLOAT gg[N1], ww[N1], qq[N1], uu[DIM][N1], uue[DIM][N1], 
         abs_g2=0., gamma, ro, sum, eps1;
   INT i=0, j;

   for (j=0; j<n; j++){
      pp[j] = 0.;
      gg[j] = -bb[j];
      ww[j] = gg[j];
      abs_g2 += gg[j]*gg[j];
   }
   eps1 = abs_g2*EPS;
   while (abs_g2>eps1){
      printf("%f\n",abs_g2);
   /* printf("%i:  %e\n",i,abs_g2); */
      m_multA(aa,ww,qq,n,n);
      ro = abs_g2/m_dot(qq,ww,n);
      sum = 0.;
      for (j = 0; j<n; j++){
         pp[j] -= ro*ww[j];
         gg[j] -= ro*qq[j];
         sum += gg[j]*gg[j];
      }
      gamma = sum/abs_g2;
      abs_g2 = sum;
      m_mult_and_add(gamma,ww,gg,ww,n);
      i++;
   }
   return(i);
}

INT m_PCG(aa,hh,pattern,bb,pp,n)
FLOAT aa[N1][N1], hh[N1][N1], bb[N1], pp[N1];
INT pattern[N1][N1];
INT n;
{
   FLOAT gg[N1], ww[N1], qq[N1], uu[DIM][N1], uue[DIM][N1];
   DOUBLE dot_qg, gamma, ro, sum, eps;
   INT i=0, j;

   for (j=0; j<n; j++){
      pp[j] = 0.;
      gg[j] = -bb[j];
   }

   m_solveILU(hh,ww,gg,n);
   dot_qg = m_dot(ww,gg,n);
   eps = EPS*m_dot(gg,gg,n);

   while ( (ro=m_dot(gg,gg,n))>eps ){
      printf("%f\n",ro);
      m_multA(aa,ww,qq,n,n);
      ro = dot_qg/m_dot(qq,ww,n);
      for (j = 0; j<n; j++){
         pp[j] -= ro*ww[j];
         gg[j] -= ro*qq[j];
      }
      m_solveILU(hh,qq,gg,n);
      sum = m_dot(qq,gg,n);
      gamma = sum/dot_qg;
      dot_qg = sum;
      m_mult_and_add(gamma,ww,qq,ww,n);
      i++;
   }
   return(i);
}

double fcn_for_minimization(n,x)
double *x;
int n;
{
   return(x[0]*x[0] + 3.*(x[1]+17.)*(x[1]+17.) + 2.*(x[2]-2.)*(x[2]-2.) + 5.);
}

void grad_fcn_for_minimization(n,x,y)
double *x, *y;
int n;
{
   y[0] = 2.*x[0];
   y[1] = 6.*(x[1]+17.);
   y[2] = 4.*(x[2]-2.);
}

double m_line_phi(a,n,x,p,z)
double a, *x, *p, *z;
int n;
{
   m_mult_and_add(a,p,x,z,n);
   return(fcn_for_minimization(n,z));
}

double m_der_line_phi(a,n,x,p,y,z)
double a, *x, *p, *y, *z;
int n;
{
   m_mult_and_add(a,p,x,z,n);
   grad_fcn_for_minimization(n,z,y);
   return(m_dot(y,p,n));
}

double m_zoom(a1,a2,phi1,phi2,dphi1,phi0,dphi0,c1,c2,n,x,p,y,z)
double a1, a2, phi1, phi2, dphi1, phi0, dphi0, c1, c2, *x, *p, *y, *z;
int n;
{
   double c, a3, phi3, dphi3;
   int i=0, k=1, max_it=100;

   while (k && i < max_it){
//    a3 = argmin q(x) where q(x) = phi1 + dphi1*(x-a1) + c*(x-a1)*(x-a1) 
//    with q(a2)=phi2 
      c = (phi2 - phi1 - dphi1*(a2-a1))/((a2-a1)*(a2-a1));
      a3 = a1 - 0.5*dphi1/c;
      if (!( (a1 <= a3 && a3 <= a2) || (a2 <= a3 && a3 <= a1) ))
         printf("Error: a3 not between a1 and a2. (i = %i)\n",i);
      phi3 =  m_line_phi(a3,n,x,p,z);
      if (phi3 > phi0 + c1*a3*dphi0 || phi3 >= phi1){
         a2 = a3;
         phi2 = phi3;
      }
      else{
         dphi3 = m_der_line_phi(a3,n,x,p,y,z);
         if (fabs(dphi3) <= -c2*dphi0)
            k = 0;
         if (dphi3*(a2-a1) >= 0.){
            a2 = a1;
            phi2 = phi1;
         }
         a1 = a3;
         phi1 = phi3;
         dphi1 = dphi3;
      }
      i++;
   }
   if (k)
      printf("Error in m_zoom: maximum number of iterations exceeded.\n");
   return(a3);
}

double m_line_search(phi0,dphi0,dphi0_old,a_old,c1,c2,n,x,p,y,z)
double phi0, dphi0, dphi0_old, a_old, c1, c2, *x, *p, *y, *z;
int n;
{
   double a1=0., a2, phi1=phi0, phi2, dphi1=dphi0, dphi2;
   int i=0, k=1, max_it=1000;

   if (dphi0 >= 0.)
      printf("Error in line_search: not a descent direction.\n");
   a2 = a_old*dphi0_old/dphi0;
   while (k && i < max_it){
      phi2 = m_line_phi(a2,n,x,p,z);
      if (phi2 > phi0 + c1*a2*dphi0 || (i && phi2 >= phi1)){
         a2 = m_zoom(a1,a2,phi1,phi2,dphi1,phi0,dphi0,c1,c2,n,x,p,y,z);
         k = 0;
      }
      else{
         dphi2 = m_der_line_phi(a2,n,x,p,y,z);
         if (fabs(dphi2) <= -c2*dphi0)
            k = 0;
         else if (dphi2 >= 0){
               a2 = m_zoom(a2,a1,phi2,phi1,dphi2,phi0,dphi0,c1,c2,n,x,p,y,z);
               k = 0;
         }
         else{
            phi1 = phi2;
            a1 = a2;
            a2 *= 2.;
            i++;
         }
      }
   }
   if (k)
      printf("Error in m_line_search: maximum number of iterations exceeded.\n");
   return(a2);
}

void m_minimize(fcn,grad_fcn)
double (*fcn)(), (*grad_fcn)();
{
   double a=1., c1=1.e-4, c2=0.1, eps=1.e-5, q, r, sum, phi0, dphi0, dphi0_old, 
          x[100], p[100], y[100], z[100], der_f[100], der_f_old[100];
   int k=0, n=3;

   m_set_value(1.,x,n);
   phi0 = fcn(n,x);
   grad_fcn(n,x,der_f);
   m_copy(der_f,der_f_old,n);
   m_inv(der_f,p,n);
   r = m_dot(der_f,der_f,n);
   dphi0_old = -r;
   while (m_max_abs_value(der_f,n) > eps*(1.+fabs(phi0))){
      dphi0 = m_dot(der_f,p,n);
      a = m_line_search(phi0,dphi0,dphi0_old,a,c1,c2,n,x,p,y,z);
      m_mult_and_add(a,p,x,x,n);
      grad_fcn(n,x,der_f);
      sum = m_dot(der_f,der_f,n);
      q = m_dot(der_f,der_f_old,n);
      if (fabs(q) > 0.1*r || q > sum)
         m_inv(der_f,p,n);
      else
         m_mult_and_subtr((sum-q)/r,p,der_f,p,n);
      r = sum;
      phi0 = fcn(n,x);
      dphi0_old = dphi0;
      m_copy(der_f,der_f_old,n);
   }
   printf("x = (%e,%e,%e)\n",x[0],x[1],x[2]);
}

void m_Stokes_defect(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,dd0,dd1,dd,n,m)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], 
      pp[N2], gg[N2], dd0[N1], dd1[N1], dd[N2];
INT n, m;
{
   FLOAT cc0[N1], cc1[N1];

   m_multA(aa,uu0,cc0,n,n);
   m_multA(aa,uu1,cc1,n,n);
   m_subtr(ff0,cc0,dd0,n);
   m_subtr(ff1,cc1,dd1,n);
   m_multAT(bb0,pp,cc0,m,n);
   m_multAT(bb1,pp,cc1,m,n);
   m_subtr(dd0,cc0,dd0,n);
   m_subtr(dd1,cc1,dd1,n);
   m_multA(bb0,uu0,cc0,m,n);
   m_multA(bb1,uu1,cc1,m,n);
   m_subtr(gg,cc0,dd,m);
   m_subtr(dd,cc1,dd,m);
}

void m_Vanka_step_1p(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,n_p,i_p,n_u,u_dofs,n,m,om)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1],
      pp[N2], gg[N2], om;
INT n_p, i_p, n_u, *u_dofs, n, m;
{
    FLOAT a[N3][N_SGGEM], b[N3][2], y[N3][N_SGGEM], u_rhs[N3][2],
          p_rhs=gg[i_p], q, s, pressure;
    INT i, j, k;

    for (i=0; i < n_u; i++){
       for (j=0; j < n_u; j++)
          a[i][j] = aa[u_dofs[i]][u_dofs[j]];          
       b[i][0] = bb0[i_p][u_dofs[i]];
       b[i][1] = bb1[i_p][u_dofs[i]];
       u_rhs[i][0] = ff0[u_dofs[i]];
       u_rhs[i][1] = ff1[u_dofs[i]];
       for (j=0; j < n; j++){
          for (k=0; k < n_u && j != u_dofs[k]; k++);
          if (k == n_u){
             u_rhs[i][0] -= aa[u_dofs[i]][j]*uu0[j];
             u_rhs[i][1] -= aa[u_dofs[i]][j]*uu1[j];
          }
       }
       for (j=0; j < m; j++)
          if (j != i_p){
             u_rhs[i][0] -= bb0[j][u_dofs[i]]*pp[j];
             u_rhs[i][1] -= bb1[j][u_dofs[i]]*pp[j];
          }
   }
/*
   for (i=0; i < n; i++){
       for (k=0; k < n_u && i != u_dofs[k]; k++);
       if (k == n_u)
          p_rhs -= bb0[i_p][i]*uu0[i] + bb1[i_p][i]*uu1[i];
   }
*/
   solve_local_Stokes(n_u,a,b,u_rhs,p_rhs,y,&pressure);
   pp[i_p] += om*(pressure - pp[i_p]);
   for (i=0; i < n_u; i++){
      uu0[u_dofs[i]] += om*(y[0][i] - uu0[u_dofs[i]]);
      uu1[u_dofs[i]] += om*(y[1][i] - uu1[u_dofs[i]]);
   }
}

void m_Vanka_step_mp(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,n_p,i_p,n_u,u_dofs,n,m,om)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1],
      pp[N2], gg[N2], om;
INT n_p, i_p, n_u, *u_dofs, n, m;
{
   FLOAT a[N3][N_SGGEM], b0[N4][N3], b1[N4][N3], c[N3][N_SGGEM], 
         y[N3][N_SGGEM], z[2][N_SGGEM], u_rhs[N3][2], r[N4], x[N4], q, s;
   INT i, j, k, i0=2*n_p, i1=2*n_p+1;

   for (i=0; i < n_u; i++){
      for (j=0; j < n_u; j++)
         a[i][j] = aa[u_dofs[i]][u_dofs[j]];          
      for (k = 0; k < n_p; k++){
         b0[k][i] = bb0[i_p+k][u_dofs[i]];
         b1[k][i] = bb1[i_p+k][u_dofs[i]];
      }
      u_rhs[i][0] = ff0[u_dofs[i]];
      u_rhs[i][1] = ff1[u_dofs[i]];
      for (j=0; j < n; j++){
         for (k=0; k < n_u && j != u_dofs[k]; k++);
         if (k == n_u){
            u_rhs[i][0] -= aa[u_dofs[i]][j]*uu0[j];
            u_rhs[i][1] -= aa[u_dofs[i]][j]*uu1[j];
         }
      }
      for (j=0; j < m; j++)
         if (j < i_p || j >= i_p+n_p){
            u_rhs[i][0] -= bb0[j][u_dofs[i]]*pp[j];
            u_rhs[i][1] -= bb1[j][u_dofs[i]]*pp[j];
         }
   }
   for (i = 0; i < n_u; i++){
      for (k = 0; k < n_p; k++){
         c[k    ][i] = b0[k][i];
         c[k+n_p][i] = b1[k][i];
      }
      c[i0][i] = u_rhs[i][0];
      c[i1][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n_u,i1+1);
   for (i = 0; i < n_p; i++)
      for (j = 0; j < n_p; j++)
         a[i][j] = 0.;
   for (i = 0; i < n_p; i++)
      r[i] = -gg[i_p+i];
   for (k = 0; k < n_u; k++)
      for (i = 0; i < n_p; i++){
         for (j = 0; j < n_p; j++)
            a[i][j] += b0[i][k]*y[j][k] + b1[i][k]*y[j+n_p][k];
         r[i] += b0[i][k]*y[i0][k] + b1[i][k]*y[i1][k];
      }
   ggem(a,r,x,n_p);
   for (i = 0; i < n_u; i++){
      z[0][i] = y[i0][i];
      z[1][i] = y[i1][i];
      for (k=0; k < n_p; k++){
         z[0][i] -= y[k][i]*x[k];
         z[1][i] -= y[k+n_p][i]*x[k];
      }
   }
   for (k = 0; k < n_p; k++)
      pp[i_p+k] += om*(x[k] - pp[i_p+k]);
   for (i=0; i < n_u; i++){
      uu0[u_dofs[i]] += om*(z[0][i] - uu0[u_dofs[i]]);
      uu1[u_dofs[i]] += om*(z[1][i] - uu1[u_dofs[i]]);
   }
}

void m_Vanka_step(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,n_p,i_p,n_u,u_dofs,n,m,om)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1], ff0[N1], ff1[N1],
      pp[N2], gg[N2], om;
INT n_p, i_p, n_u, *u_dofs, n, m;
{
   if (n_p == 1)
    m_Vanka_step_1p(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,n_p,i_p,n_u,u_dofs,n,m,om);
   else
    m_Vanka_step_mp(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,n_p,i_p,n_u,u_dofs,n,m,om);
}

void m_Vanka_smoother(j,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,vanka_np,vanka_p,vanka_nu,vanka_u,n,m)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], 
      pp[N2], gg[N2], omega;
INT j, n, m, vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   INT i;

   if (j == 1 || j == 3)
      printf(" M Smoothing.\n");
   for (i=0; i < m; i+=vanka_np)
      m_Vanka_step(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,
                   vanka_np,vanka_p[i],vanka_nu[i],vanka_u[i],n,m,omega);
}

void m_Stokes_smoothing_solver(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,tol1,tol2,omega,vanka_np,vanka_p,vanka_nu,vanka_u,n,m)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], pp[N2], gg[N2], tol1, tol2, omega;
INT n, m, vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3];
{
   FLOAT dd0[N1], dd1[N1], dd[N2], u_rhs, p_rhs, err44, err55, eps1, eps2;
   INT i=0, imin=0, imax=1000;
INT k;
   u_rhs = m_dot(ff0,ff0,n)+m_dot(ff1,ff1,n);
   p_rhs = m_dot(gg,gg,m);
   if (u_rhs < 1.e-15)
      u_rhs = 1.;
   if (p_rhs < 1.e-15)
      p_rhs = 1.;
   eps1 = u_rhs*tol1*tol1;
   eps2 = p_rhs*tol2*tol2;
   if (eps1 < EPST || eps2 < EPST){
      imin = 2;
      eps1 = MAX(eps1,EPST);
      eps2 = MAX(eps2,EPST);
   }

   m_set_value(0.,uu0,n);
   m_set_value(0.,uu1,n);
   m_set_value(0.,pp,m);
   err44 = err55 = 1.;
   while (i < imin || (i < imax && (err44 > eps1 || err55 > eps2))){
      m_Vanka_smoother(0,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,
                       vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
      m_Stokes_defect(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,dd0,dd1,dd,n,m);
      err44=m_dot(dd0,dd0,n)+m_dot(dd1,dd1,n);
      err55=m_dot(dd,dd,m);
      i++;
   }
   printf("     %i iterations   u_def/rhs: %e  p_def/rhs: %e\n",i,sqrt(err44/u_rhs),sqrt(err55/p_rhs));
}

void two_level_method(aa,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,aa_c,bb0_c,bb1_c,u_pr,p_pr,
                     tol1,tol2,omega,n,m,n_c,m_c,
    vanka_np,vanka_p,vanka_nu,vanka_u,vanka_np_c,vanka_p_c,vanka_nu_c,vanka_u_c)
FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  ff0[N1], ff1[N1], 
      pp[N2], gg[N2], aa_c[N1][N1], bb0_c[N2][N1], bb1_c[N2][N1], 
      u_pr[N1][N1], p_pr[N1][N1], tol1, tol2, omega;
INT n, m, n_c, m_c, vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3], 
    vanka_np_c, vanka_p_c[N2], vanka_nu_c[N2], vanka_u_c[N2][N3];
{
   FLOAT uu0_c[N1], uu1_c[N1],  ff0_c[N1], ff1_c[N1], pp_c[N2], gg_c[N2], 
         dd0[N1], dd1[N1], dd[N2];

   printf("Level %1i:   %i smoothings.\n",2,2);
   m_Vanka_smoother(0,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
   m_Vanka_smoother(0,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
   m_Stokes_defect(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,dd0,dd1,dd,n,m);
   m_multAT(u_pr,dd0,ff0_c,n,n_c);
   m_multAT(u_pr,dd1,ff1_c,n,n_c);
   m_multAT(p_pr,dd,gg_c,m,m_c);
   printf("Level %1i:   solving.\n",1);
   m_Stokes_smoothing_solver(aa_c,bb0_c,bb1_c,uu0_c,uu1_c,pp_c,ff0_c,ff1_c,gg_c,
                         tol1,tol2,omega,vanka_np_c,vanka_p_c,vanka_nu_c,vanka_u_c,n_c,m_c);
   m_multA(u_pr,uu0_c,dd0,n,n_c);
   m_multA(u_pr,uu1_c,dd1,n,n_c);
   m_multA(p_pr,pp_c,dd,m,m_c);
   m_add(dd0,uu0,uu0,n);
   m_add(dd1,uu1,uu1,n);
   m_add(dd,pp,pp,m);
   printf("Level %1i:   %i smoothings.\n",2,2);
   m_Vanka_smoother(0,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
   m_Vanka_smoother(0,omega,aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
}

void tlm_mtest(tGrid,ZA,ZB,u,p,f,g)
GRID *tGrid;      
INT ZA, ZB, u, p, f, g;
{
   FLOAT aa[N1][N1], bb0[N2][N1], bb1[N2][N1], uu0[N1], uu1[N1],  
         ff0[N1], ff1[N1], pp[N2], gg[N2], 
         aa_c[N1][N1], bb0_c[N2][N1], bb1_c[N2][N1], uu0_c[N1], uu1_c[N1],  
         ff0_c[N1], ff1_c[N1], pp_c[N2], gg_c[N2], 
         u_pr[N1][N1], p_pr[N1][N1], dd0[N1], dd1[N1], dd[N2], 
         u_rhs, p_rhs, err4, err44, err5, err55, 
         tol1=1.e-8, tol2=1.e-8, omega=1.;
   INT i, n, m, n_c, m_c, pattern[N1][N1], pattern_c[N1][N1], 
       vanka_np, vanka_p[N2], vanka_nu[N2], vanka_u[N2][N3], 
       vanka_np_c, vanka_p_c[N2], vanka_nu_c[N2], vanka_u_c[N2][N3],
       u_space=U_SPACE, p_space=P_SPACE;

   structures_to_matrices(tGrid,ZA,ZB,u,p,f,g,u_space,p_space,&n,&m,
                  aa,pattern,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,
                  &vanka_np,vanka_p,vanka_nu,vanka_u);
   if (LOW_AND_HIGH == YES){
      structures_to_matrices(tGrid,A_LOW,B_LOW,u,p,f,g,P1_NC,P0,&n_c,&m_c,
                  aa_c,pattern_c,bb0_c,bb1_c,uu0_c,uu1_c,ff0_c,ff1_c,pp_c,gg_c,
                  &vanka_np_c,vanka_p_c,vanka_nu_c,vanka_u_c);
      if (u_space == P2C)
         u_p1nc_to_p2c_prolong_matr(tGrid,u_pr,n,n_c);
      else
         u_p1nc_to_p2cb_prolong_matr(tGrid,u_pr,n,n_c);
      if (p_space == P1C)
         p_p0_to_p1c_prolong_matr(tGrid,p_pr,m,m_c);
      else
         p_p0_to_p1disc_prolong_matr(p_pr,m,m_c);
   }
   else{
      structures_to_matrices(tGrid->coarser,ZA,ZB,u,p,f,g,u_space,p_space,
               &n_c,&m_c,aa_c,pattern_c,bb0_c,bb1_c,uu0_c,uu1_c,ff0_c,ff1_c,
               pp_c,gg_c,&vanka_np_c,vanka_p_c,vanka_nu_c,vanka_u_c);
      u_prolong_matr(tGrid->coarser,u_pr,n,n_c,u_space);
      p_prolong_matr(tGrid->coarser,p_pr,m,m_c,p_space);
   }
   u_rhs = m_dot(ff0,ff0,n)+m_dot(ff1,ff1,n);
   if (u_rhs < 1.e-15)
      u_rhs = 1.;
   p_rhs = m_dot(gg,gg,m);
   if (p_rhs < 1.e-15)
      p_rhs = 1.;
   m_set_value(0.,uu0,n);
   m_set_value(0.,uu1,n);
   m_set_value(0.,pp,m);
   m_Stokes_defect(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,dd0,dd1,dd,n,m);
   err4=sqrt((m_dot(dd0,dd0,n)+m_dot(dd1,dd1,n))/u_rhs);
   err5=sqrt(m_dot(dd,dd,m)/p_rhs);
   printf("MM: u_def: %e  p_def: %e\n",sqrt((m_dot(dd0,dd0,n)+m_dot(dd1,dd1,n))),sqrt(m_dot(dd,dd,m)));
   printf("MM: u_def/rhs: %e  p_def/rhs: %e\n",err4,err5);
   err44 = err55 = 1.;
   for (i=0; i<200000  && (err44 > tol1 || err55 > tol2); i++){
      printf("%i\n",i);
      two_level_method(aa,bb0,bb1,uu0,uu1,ff0,ff1,pp,gg,aa_c,bb0_c,bb1_c,u_pr,p_pr,
                      tol1,tol2,omega,n,m,n_c,m_c,
                      vanka_np,vanka_p,vanka_nu,vanka_u,vanka_np_c,vanka_p_c,vanka_nu_c,vanka_u_c);
      m_Stokes_defect(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,dd0,dd1,dd,n,m);
      err44=sqrt((m_dot(dd0,dd0,n)+m_dot(dd1,dd1,n))/u_rhs);
      err55=sqrt(m_dot(dd,dd,m)/p_rhs);
      printf("MM: u_def/rhs: %e  p_def/rhs: %e\n",err44,err55);
   }
   printf("%i iterations\n",i);
//   m_Stokes_smoothing_solver(aa,bb0,bb1,uu0,uu1,pp,ff0,ff1,gg,tol1,tol2,omega,vanka_np,vanka_p,vanka_nu,vanka_u,n,m);
}
