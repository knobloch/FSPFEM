/******************************************************************************/
/*                                                                            */
/*                                   tests                                    */
/*                                                                            */
/******************************************************************************/

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA)

void nc_structure_to_matrix(tGrid,Z,u,f,n,m,aa,pattern,uu,ff)
GRID *tGrid;      
INT Z,u,f;
INT *n, *m;
FLOAT aa[N2][N2], uu[N1], ff[N1];
INT pattern[N2][N2];
{
   INT i, j, k, l;
   ELEMENT *pel;
   FACE *pface;
   FLINK *pfl;

   for (i = 0; i < N2; i++)
      for(j = 0; j < N2; j++){
         aa[i][j] = 0.;
         pattern[i][j] = 0;
      }
   i = 0;
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface=pface->succ)
      pface->index = i++;
   *n = i;
   i = 0;
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ)
      pel->status = i++;
   *m = i;
   printf("N1 = %i,  N2 = %i,  n = %i, m = %i\n",N1,N2,*n,*m);
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface=pface->succ){
      i = pface->index;
      aa[i][i] = COEFF_FF(pface,Z);
      pattern[i][i] = 1;
      uu[i] = FD(pface,u);
      ff[i] = FD(pface,f);
      for (pfl = pface->fstart; pfl != NULL; pfl = pfl->next){
	  j = pfl->nbface->index;
	  aa[i][j] = COEFF_FL(pfl,Z);
          pattern[i][j] = 1;
      }
   }
}

void nc_structure_to_matrix_bub(tGrid,Z1,Z2,Z3,Z4,u1,u2,f1,f2,n,me,aa,pattern,uu,ff)
GRID *tGrid;      
INT Z1,Z2,Z3,Z4,u1,u2,f1,f2;
INT *n, *me;
FLOAT aa[N2][N2], uu[N1], ff[N1];
INT pattern[N2][N2];
{
   INT i, j, k, l, m;
   ELEMENT *pel;
   FACE *pface;
   FLINK *pfl;

   for (i = 0; i < N2; i++)
      for(j = 0; j < N2; j++){
         aa[i][j] = 0.;
         pattern[i][j] = 0;
      }
   i = 0;
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface=pface->succ)
      pface->index = i++;
   m = i;
   *n = 2*i;
   i = 0;
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ)
      pel->status = i++;
   *me = i;
   printf("N1 = %i,  N2 = %i,  n = %i, m = %i\n",N1,N2,*n,*me);
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface=pface->succ){
      i = pface->index;
      aa[i][i] = COEFF_FF(pface,Z1);
      aa[i][i+m] = COEFF_FF(pface,Z2);
      aa[i+m][i] = COEFF_FF(pface,Z3);
      aa[i+m][i+m] = COEFF_FF(pface,Z4);
      pattern[i][i] = 1;
      pattern[i][i+m] = 1;
      pattern[i+m][i] = 1;
      pattern[i+m][i+m] = 1;
      uu[i] = FD(pface,u1);
      uu[i+m] = FD(pface,u2);
      ff[i] = FD(pface,f1);
      ff[i+m] = FD(pface,f2);
      for (pfl = pface->fstart; pfl != NULL; pfl = pfl->next){
	  j = pfl->nbface->index;
	  aa[i][j] = COEFF_FL(pfl,Z1);
	  aa[i][j+m] = COEFF_FL(pfl,Z2);
	  aa[i+m][j] = COEFF_FL(pfl,Z3);
	  aa[i+m][j+m] = COEFF_FL(pfl,Z4);
          pattern[i][j] = 1;
          pattern[i][j+m] = 1;
          pattern[i+m][j] = 1;
          pattern[i+m][j+m] = 1;
      }
   }
}

#else

void nc_structure_to_matrix(tGrid,Z,u,f,n,m,aa,pattern,uu,ff)
GRID *tGrid; INT Z,u,f; INT *n, *m; FLOAT aa[N2][N2], uu[N1], ff[N1]; INT pattern[N2][N2];
{  eprintf("Error: nc_structure_to_matrix not available.\n");  }

void nc_structure_to_matrix_bub(tGrid,Z1,Z2,Z3,Z4,u1,u2,f1,f2,n,me,aa,pattern,uu,ff)
GRID *tGrid; INT Z1,Z2,Z3,Z4,u1,u2,f1,f2; INT *n, *me; FLOAT aa[N2][N2], uu[N1], ff[N1]; INT pattern[N2][N2];
{  eprintf("Error: nc_structure_to_matrix_bub not available.\n");  }

#endif

#if (E_DATA & ExDF_MATR) && (F_DATA & VECTOR_FACE_DATA)

void nc_b_to_matr(tGrid,bb,uu,ff,u,f,n,m)
GRID *tGrid;      
FLOAT bb[N2][N2], uu[N1], ff[N1];
INT u, f, n, m;
{
   GRID *theGrid;
   ELEMENT *pel;
   FACE *pface;
   INT i, j;

   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         bb[i][j] = 0.;
   for (pel = SUCC(FIRSTELEMENT(tGrid)); pel!=NULL; pel=pel->succ){
     i = pel->status;
     for (j=0; j < 3; j++)
        if (IS_FF(pel->f[j])){
           bb[i][pel->f[j]->index]   = COEFF_BDF(pel,0,j,0);
           bb[i][pel->f[j]->index+n] = COEFF_BDF(pel,0,j,1);
        }
     ff[i] = ED(pel,f);
   }
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface=pface->succ){
      i = pface->index;
      uu[i]   = FDV(pface,u,0);
      uu[i+n] = FDV(pface,u,1);
   }
}

#else

void nc_b_to_matr(tGrid,bb,uu,ff,u,f,n,m)
GRID *tGrid;      
FLOAT bb[N2][N2], uu[N1], ff[N1];
INT u, f, n, m;
{  eprintf("Error: nc_b_to_matr not available.\n");  }

#endif

void nc_matrix_test(tGrid,Z,u,f)
GRID *tGrid;      
INT Z,u,f;
{
   FLOAT aa[N2][N2], hh[N2][N2], bb[N2][N2], uu[N1], ff[N1], gg[N1];
   INT pattern[N2][N2];
   INT n, m, i;

/* nc_structure_to_matrix(tGrid,1,u,f,&n,&m,hh,pattern,uu,ff); */
   nc_structure_to_matrix(tGrid,Z,u,f,&n,&m,aa,pattern,uu,ff);
   nc_b_to_matr(tGrid,bb,uu,ff,U,F,n,m);
   tm_CG(bb,ff,uu,m,2*n); 
/*   m_ILU(aa,hh,pattern,n); */
/*   m_PCG(aa,hh,pattern,ff,uu,n); */
/* m_CG(aa,ff,uu,n); */
}

void nc_matrix_test_bub(tGrid,Z1,Z2,Z3,Z4,u1,u2,f1,f2)
GRID *tGrid;      
INT Z1,Z2,Z3,Z4,u1,u2,f1,f2;
{
   FLOAT aa[N2][N2], hh[N2][N2], uu[N1], ff[N1];
   INT pattern[N2][N2];
   INT n, m, i;

   nc_structure_to_matrix_bub(tGrid,Z1,Z2,Z3,Z4,u1,u2,f1,f2,&n,&m,aa,pattern,uu,ff);
   m_ILU(aa,hh,pattern,n); 
   m_PCG(aa,hh,pattern,ff,uu,n); 
/*   m_CG(aa,ff,uu,n); */
}
