/******************************************************************************/
/*                                                                            */
/*                                   solver                                   */
/*                                                                            */
/******************************************************************************/

#if INCLUDE_UMFPACK == YES

void solve_system_using_umfpack(n,Ap,Ai,Ax,b,x)
INT n, *Ap, *Ai;
DOUBLE *Ax, *b, *x;
{
    DOUBLE *null = (double *) NULL ;

    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
}

#else

void solve_system_using_umfpack(n,Ap,Ai,Ax,b,x)
INT n, *Ap, *Ai; DOUBLE *Ax, *b, *x;
{  eprintf("Error: solve_system_using_umfpack not available.\n");  }

#endif

void augmented_lagrange(tGrid,r,Z,f,fe,u,ue,b,g,w,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, r1, r2;
INT Z;      /* velocity matrix */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q, d, de;  /* auxiliary fields */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, rhs_f, rhs_fe;
   INT i=0;
  
   printf(" Augmented Lagrange.\n");
   
   rhs_f = vs_dot_nf(tGrid,f,f,1);
   rhs_fe = dot_e(tGrid,fe,fe);
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   vs_set_value_nf(tGrid,0.0,u,1);
   set_value_e(tGrid,0.0,ue);
   mult_e(tGrid,r,fe,de);
   mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
   vs_add_nf(tGrid,f,b,b,1);
   CG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
      d,T_FOR_U,U_TYPE,A_STRUCT,
      Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
   Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                 U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
   udef = vs_dot_nf(tGrid,d,d,1);
   pdef = dot_e(tGrid,de,de);   
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,r,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      CG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
         d,T_FOR_U,U_TYPE,A_STRUCT,
         Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
      Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
/*  printf("  (u_def/rhs)^2: %e  (p_def/rhs)^2: %e\n",udef/rhs_f,pdef/rhs_fe);*/
   }
   printf("     %3i outer cykles.   u_def/rhs: %e  p_def/rhs: %e\n",
                                          i,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
}

INT raugmented_lagrange2(tGrid,r,rho,Z,f,fe,u,ue,b,g,w,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, rho, r1, r2;
INT Z;      /* velocity matrix */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q;  /* auxiliary fields */
INT d, de;       /* defect; then auxiliary fields */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, rhs_f, rhs_fe;
   INT i=0, j=0;
  
   printf(" Augmented Lagrange.\n");
   
   udef = rhs_f = vs_dot_nf(tGrid,d,d,1);
   pdef = rhs_fe = dot_e(tGrid,de,de);  
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,rho,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      j += CG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
              d,T_FOR_U,U_TYPE,A_STRUCT,
              Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
/*    j += PCG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
               d,T_FOR_U,U_TYPE,A_STRUCT,
               Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.,
               X,0,0,0,DIAG_PR); */
      Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
   }
   printf("     %3i outer cykles.     %3i inner cykles.   u_def/in_def: %e  p_def/in_def: %e\n",
                                        i,j,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
   return(j);
}

INT augmented_lagrange2(tGrid,r,Z,f,fe,u,ue,b,g,w,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, r1, r2;
INT Z;      /* velocity matrix */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q;  /* auxiliary fields */
INT d, de;       /* defect; then auxiliary fields */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, rhs_f, rhs_fe;
   INT i=0, j=0;;
  
   printf(" Augmented Lagrange.\n");
   
   udef = rhs_f = vs_dot_nf(tGrid,d,d,1);
   pdef = rhs_fe = dot_e(tGrid,de,de);  
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,r,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      j += CG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
              d,T_FOR_U,U_TYPE,A_STRUCT,
              Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
/*
      j += PCG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
               d,T_FOR_U,U_TYPE,A_STRUCT,
               Z,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.,
               X,0,0,0,DIAG_PR);
*/
      Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
   }
   printf("     %3i outer cykles.     %3i inner cykles.   u_def/in_def: %e  p_def/in_def: %e\n",
                                        i,j,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
   return(j);
}

INT augmented_lagrange3(tGrid,r,Z1,Z2,f,fe,u,ue,b,g,h,w,p,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, r1, r2;
INT Z1, Z2; /* velocity matrices */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, h, w, p, q;  /* auxiliary fields */
INT d, de;             /* defect; then auxiliary fields */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, eps3, rhs_f, rhs_fe;
   INT i=0, j=0, k;
  
   printf(" Augmented Lagrange.\n");
   
   udef = rhs_f = vs_dot_nf(tGrid,d,d,1);
   pdef = rhs_fe = dot_e(tGrid,de,de);  
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,r,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      mult_Ar(tGrid,Z1,u,h,d,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
              0,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
      vs_subtr_nf(tGrid,b,h,h,1);
      eps3 = EPS_CGAR*vs_dot_nf(tGrid,h,h,1);
      k = 0.;
      while (vs_dot_nf(tGrid,h,h,1) > eps3 && k < 30){
         vs_set_value_nf(tGrid,0.,p,1);
         j += PCG(tGrid,h,p,g,w,q,0,300,.7,1.e50,NONZERO_INIT,mult_Ar,Z2,
                  d,T_FOR_U,U_TYPE,A_STRUCT,
                  0,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.,
                  X,0,0,0,DIAG_PR);
         vs_mult_and_add_nf(tGrid,1.,p,u,u,1);
         mult_Ar(tGrid,Z1,u,h,d,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
                 0,B_STRUCT,d,de,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
         vs_subtr_nf(tGrid,b,h,h,1);
         k++;
      }
      printf("   %i\n",k);
      Stokes_defect(tGrid,Z1,Z1,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
   }
   printf("     %3i outer cykles.     %3i inner cykles.   u_def/in_def: %e  p_def/in_def: %e\n",
                                        i,j,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
   return(j);
}

void Augmented_lagrange(tGrid,r,Z,f,fe,u,pe,b,g,w,q,h,he,imin,imax,r1,eps)
GRID *tGrid;
FLOAT r, r1, eps;
INT Z;      /* velocity matrix */
INT u, pe;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q, h, he;  /* auxiliary fields */
INT imin, imax;
{  
   FLOAT res, eps1;
   INT i=0;
  
   printf(" Augmented Lagrange.\n");
   
   vs_set_value_nf(tGrid,0.0,u,1);
   set_value_e(tGrid,0.0,pe);
   res = vs_dot_nf(tGrid,f,f,1) + dot_e(tGrid,fe,fe);
   eps1 = MIN(eps,res*r1);
   
   while (i < imin || (i < imax && res > eps1)){
      i++;
      mult_and_add_e(tGrid,-r,fe,pe,he);
      mult_BT(tGrid,0,he,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_subtr_nf(tGrid,f,b,b,1);
      CG(tGrid,b,u,g,w,q,0,300,EPS_CGAR,1.e50,NONZERO_INIT,mult_Ar,Z,
         h,T_FOR_U,U_TYPE,A_STRUCT,
         Z,B_STRUCT,h,he,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
      mult_B(tGrid,0,u,he,0,0,0,Q_SE,Q_VNSF,0);
      subtr_e(tGrid,he,fe,he);
      res = dot_e(tGrid,he,he);
      mult_and_add_e(tGrid,r,he,pe,pe);
      mult_BT(tGrid,0,pe,h,h,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      mult_A(tGrid,Z,u,q,q,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
         0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
      vs_add_nf(tGrid,q,h,h,1);
      vs_subtr_nf(tGrid,f,h,h,1);
      res += vs_dot_nf(tGrid,h,h,1);
   }
   printf("     %3i outer cykles.   u_def: %e  p_def: %e\n",
                                     i,vs_dot_nf(tGrid,h,h,1),dot_e(tGrid,he,he));
}

void augmented_lagrange_GMRES(tGrid,r,Z,f,fe,u,ue,b,g,w,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, r1, r2;
INT Z;      /* velocity matrix */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q, d, de;  /* auxiliary fields; g, w not used */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, rhs_f, rhs_fe;
   INT i=0;
  
   printf(" Augmented Lagrange.\n");
   
   rhs_f = vs_dot_nf(tGrid,f,f,1);
   rhs_fe = dot_e(tGrid,fe,fe);
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   vs_set_value_nf(tGrid,0.0,u,1);
   set_value_e(tGrid,0.0,ue);
   mult_e(tGrid,r,fe,de);
   mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
   vs_add_nf(tGrid,f,b,b,1);
   GMRES(tGrid,u,b,0,0,10,10,10,0.1,1.e50,NONZERO_INIT,mult_Ar,Z,
         d,T_FOR_U,U_TYPE,A_STRUCT,Z,B_STRUCT,d,de,
         T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
   Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                 U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
   udef = vs_dot_nf(tGrid,d,d,1);
   pdef = dot_e(tGrid,de,de);
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,r,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      GMRES(tGrid,u,b,0,0,10,10,10,0.1,1.e50,NONZERO_INIT,mult_Ar,Z,
            d,T_FOR_U,U_TYPE,A_STRUCT,Z,B_STRUCT,d,de,
            T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
      Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
/*  printf("  (u_def/rhs)^2: %e  (p_def/rhs)^2: %e\n",udef/rhs_f,pdef/rhs_fe);*/
   }
   printf("     %3i outer cykles.   u_def/rhs: %e  p_def/rhs: %e\n",
                                          i,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
}

INT augmented_lagrange2_GMRES(tGrid,r,Z,f,fe,u,ue,b,g,w,q,d,de,imin,imax,r1,r2)
GRID *tGrid;
FLOAT r, r1, r2;
INT Z;      /* velocity matrix */
INT u, ue;  /* solution */
INT f, fe;  /* rhs      */
INT b, g, w, q;  /* auxiliary fields; g, w not used */
INT d, de;       /* defect; then auxiliary fields */
INT imin, imax;
{  
   FLOAT udef, pdef, eps1, eps2, rhs_f, rhs_fe;
   INT i=0, j=0;
  
   printf(" Augmented Lagrange.\n");
   
   udef = rhs_f = vs_dot_nf(tGrid,d,d,1);
   pdef = rhs_fe = dot_e(tGrid,de,de);  
   eps1 = rhs_f*r1*r1;
   eps2 = rhs_fe*r2*r2;
   
   while (i < imin || (i < imax && (udef > eps1 || pdef > eps2))){
      i++;
      mult_and_add_e(tGrid,r,de,ue,ue);
      mult_e(tGrid,r,fe,de);
      subtr_e(tGrid,de,ue,de);
      mult_BT(tGrid,0,de,b,b,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
      vs_add_nf(tGrid,f,b,b,1);
      j += GMRES(tGrid,u,b,0,0,50,10,10,0.224,1.e50,NONZERO_INIT,mult_Ar,Z,
                 d,T_FOR_U,U_TYPE,A_STRUCT,Z,B_STRUCT,d,de,
                 T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
/*
      j += PGMRES_Ar(tGrid,Z,u,b,de,d,0,0,50,NONZERO_INIT,10,10,r,0.224,1.e50,X);
*/
      Stokes_defect(tGrid,Z,Z,f,fe,u,ue,d,de,q,T_FOR_U,T_FOR_P,
                    U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
      udef = vs_dot_nf(tGrid,d,d,1);
      pdef = dot_e(tGrid,de,de);
   }
   printf("     %3i outer cykles.     %3i inner cykles.   u_def/in_def: %e  p_def/in_def: %e\n",
                                        i,j,sqrt(udef/rhs_f),sqrt(pdef/rhs_fe));
   return(j);
}

