/******************************************************************************/
/*                                                                            */
/*                           comparison of matrices                           */
/*                                                                            */
/******************************************************************************/

#define F_NAME        "Aug"

void nc_print_matr_mbub(tGrid,Z,f)
GRID *tGrid;
INT Z, f;
{
   NODE *n0, *n1, *n2;
   FACE *pface, *fa0, *fa1, *fa2;
   FLINK *pflink;
   ELEMENT *pelem;
   FLOAT rdetB, ndetB, b[DIM2][DIM2];

   for (pface=FIRSTFACE(tGrid);pface!=NULL;pface=pface->succ) {
      printf("diag: %e,  %e,  %e,  %e\n",COEFFNN(pface,Z,0,0),
      COEFFNN(pface,Z,0,1),COEFFNN(pface,Z,1,0),COEFFNN(pface,Z,1,1));
      printf("rhs:  %e, %e\n",FDV(pface,f,0),FDV(pface,f,1));
      for (pflink=pface->fstart;pflink!=NULL;pflink=pflink->next)
         printf("      %e,  %e,  %e,  %e\n",COEFFNN(pflink,Z,0,0),
         COEFFNN(pflink,Z,0,1),COEFFNN(pflink,Z,1,0),COEFFNN(pflink,Z,1,1));
   }
   printf("\n");
}

void nc_print_diff_matr_mbub(tGrid,Z,s,f,g)
GRID *tGrid;
INT Z, s, f, g;
{
   NODE *n0, *n1, *n2;
   FACE *pface, *fa0, *fa1, *fa2;
   FLINK *pflink;
   ELEMENT *pelem;
   FLOAT rdetB, ndetB, b[DIM2][DIM2];

   for (pface=FIRSTFACE(tGrid);pface!=NULL;pface=pface->succ) {
      printf("diag: %e,  %e,  %e,  %e\n",
      COEFFNN(pface,Z,0,0)-COEFFNN(pface,s,0,0),
      COEFFNN(pface,Z,0,1)-COEFFNN(pface,s,0,1),
      COEFFNN(pface,Z,1,0)-COEFFNN(pface,s,1,0),
      COEFFNN(pface,Z,1,1)-COEFFNN(pface,s,1,1));
      printf("rhs:  %e, %e\n",FDV(pface,f,0)-FDV(pface,g,0),
                              FDV(pface,f,1)-FDV(pface,g,1));
      for (pflink=pface->fstart;pflink!=NULL;pflink=pflink->next)
         printf("      %e,  %e,  %e,  %e\n",
         COEFFNN(pflink,Z,0,0)-COEFFNN(pflink,s,0,0),
         COEFFNN(pflink,Z,0,1)-COEFFNN(pflink,s,0,1),
         COEFFNN(pflink,Z,1,0)-COEFFNN(pflink,s,1,0),
         COEFFNN(pflink,Z,1,1)-COEFFNN(pflink,s,1,1));
   }
   printf("\n");
}

void stest6(mg,plink,u0,u01,u02,u03,t0,t)  /* u0 ... solution; t0 ... rhs */
MULTIGRID *mg;
LINK **plink;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*u03)(), (*t0)();
INT t;
{
   NODE *pnode;
   FLOAT err, sum1=0.0, sum2=0.0, max=0.0;
   FLOAT err1, err1lin, err1cube, err2, err2lin, err2cube, err3, err4, err5,
      err44, err55;
   FILE *fp;

   fp = fopen(F_NAME,"w");
   fprintf(fp,"--------------------\n");
   fclose(fp);

   boundary_values(mg,u0,u0,u0,U,D,G,U_SPACE,SCALAR);
   add_Laplace_matr(TOP_GRID(mg),TNU,A,U_SPACE,A_STRUCT,SCALAR,KORN_LAPLACE);
   integrate_rhs(TOP_GRID(mg),F,ft0,ft0,ft0,
                                      T_FOR_U,U_TYPE,U_SPACE,SCALAR,RHS_INTEGR);
   add_Laplace_matr(TOP_GRID(mg),TNU,1,U_SPACE,A_STRUCT,SCALAR,KORN_LAPLACE);
   integrate_rhs(TOP_GRID(mg),G,ft0,ft0,ft0,
                                      T_FOR_U,U_TYPE,U_SPACE,SCALAR,RHS_INTEGR);
   nc_conv_term_to_stiff_matr_mbub(TOP_GRID(mg),A,F,U,bb0,bb1,CONV);
   nc_conv_term_to_stiff_matr_mbub(TOP_GRID(mg),1,G,U,bb0,bb1,SKEW);
   nc_print_diff_matr_mbub(TOP_GRID(mg),A,1,F,G);
   printf("convective:\n");
   nc_print_matr_mbub(TOP_GRID(mg),A,F);
   printf("skew:\n");
   nc_print_matr_mbub(TOP_GRID(mg),1,G);
   mult_A(TOP_GRID(mg),A,F,U,U,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   printf("convective: %e\n",vdot_f(TOP_GRID(mg),F,U,1));
   mult_A(TOP_GRID(mg),1,F,U,U,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   printf("skew: %e\n",vdot_f(TOP_GRID(mg),F,U,1));
}
