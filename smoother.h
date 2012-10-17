/******************************************************************************/
/*                                                                            */
/*                                 smoothing                                  */
/*                                                                            */
/******************************************************************************/

#if (SOLVER_IN_SMOOTHER == BASIC_CG)
#define SOLVER_IN_BS_SMOOTHER(multB_X_BT)                                      \
   j = CG(tGrid,de,pe,qe,xe,ye,0,5000,0.01,EPSE,0,multB_X_BT,ZB,               \
          0,t_p,p_type,B_struct,precond_matrix,p,q,t_u,u_type,A_struct,        \
          6,7,8,9,0.,0.,0.);
#define BS_PRECOND(d,q)                                                        \
   ii = PCG(tGrid,d,q,J,K,L,0,500,0.01,1.,0,mult_A,A,                          \
       G,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,                 \
       1,A_STRUCT,0,0,ILU_PR);

#elif (SOLVER_IN_SMOOTHER == BASIC_GMRES)
#define SOLVER_IN_BS_SMOOTHER(multB_X_BT)                                      \
   j = GMRES(tGrid,pe,de,3,0,5000,qe,m,0.1,sqrt(EPSE),0,multB_X_BT,ZB,         \
          0,t_p,p_type,B_struct,precond_matrix,p,q,t_u,u_type,A_struct,        \
          6,7,8,9,0.,0.,0.);
#define BS_PRECOND(d,q)                                                        \
   ii = PGMRES(tGrid,q,d,J,K,0,0,500,L,15,0.1,1.,0,mult_A,A,                      \
       G,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,                 \
       1,A_STRUCT,0,0,ILU_PR);

#else
#define SOLVER_IN_BS_SMOOTHER(multB_X_BT)                                      \
   eprintf("Error: SOLVER_IN_BS_SMOOTHER not available.\n");
#define BS_PRECOND(d,q)                                                        \
   eprintf("Error: BS_PRECOND not available.\n");
#endif

INT BS_smoothing_step(tGrid,i,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
                      d,p,q,de,pe,qe,xe,ye,
                      t_u,t_p,u_type,p_type,A_struct,B_struct,
                      q0,q1,q2,precond_type,multB_X_BT)
GRID *tGrid;
FLOAT alpha;
INT i, ZA, ZB, precond_matrix, f, fe, u, ue, d, p, q, de, pe, qe, xe, ye;
INT t_u, t_p, u_type, p_type, A_struct, B_struct, q0, q1, q2, precond_type;
void (*multB_X_BT)();
{
   INT j, m=15, n; /*  GMRES used  =>  qe,...,qe+m+1 are auxiliary variables  */
   FLOAT s;
   INT C_struct=C_STRUCT;

   if (i == 1){ 
      printf(" Smoothing.\n");
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
   }
   else if (i == 2)
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
   else if (i == 3)
      printf(" Smoothing.\n");
   preconditioning(tGrid,precond_matrix,d,q,t_u,u_type,q0,q1,q2,precond_type); 
   mult_B(tGrid,ZB,q,qe,qe,t_p,t_u,p_type,u_type,B_struct);
   mult_and_add(tGrid,-alpha,de,qe,de,t_p,p_type);
   SOLVER_IN_BS_SMOOTHER(multB_X_BT)
/* For MINI also the following is possible:
   j = PCG(tGrid,de,pe,qe,xe,ye,0,500,0.09,EPSE,0,multB_X_BT,ZB,
           0,t_p,p_type,B_struct,precond_matrix,p,q,t_u,u_type,A_struct,
           6,7,8,9,0.,0.,0.,X,0,0,0,DIAG_PR); */
   mult_BT(tGrid,ZB,pe,p,p,t_p,t_u,p_type,u_type,B_struct);
   subtr(tGrid,d,p,d,t_u,u_type);
   preconditioning(tGrid,precond_matrix,d,q,t_u,u_type,q0,q1,q2,precond_type); 
   mult_and_add(tGrid,-1./alpha,q,u,u,t_u,u_type);
   subtr(tGrid,ue,pe,ue,t_p,p_type);
   if (t_p & ZERO_MEAN){
      s = sum_n(tGrid,ue,&n,t_p,p_type);
      add_value(tGrid,-s/n,ue,t_p,p_type);
   }
   if (i == 1 || i == 3)
      printf("     %i iterations.\n",j);
   return(j);
}

void multB_ILU_BT_ex(tGrid,matrix_b,x,y,q,t_row,t_column,row_type,column_type,
         b_struct,ilu_matrix,u1,u2,t_for_u,u_type,a_struct,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;  /*  t_row = t_column = t_for_p; row_type = column_type = p_type */
INT x, y, q;  /*  p_type  */
INT u1, u2;   /*  u_type  */
INT matrix_b,t_row,t_column,row_type,column_type,b_struct,
    ilu_matrix,t_for_u,u_type,a_struct,p6,p7,p8,p9;
FLOAT r1, r2, r3;
{  
   INT ii;

   mult_BT(tGrid,matrix_b,x,u1,u1,t_row,t_for_u,row_type,u_type,b_struct);
   BS_PRECOND(u1,u2)
   printf("A^{-1} it.: %i\n",ii);
   mult_B(tGrid,matrix_b,u2,y,y,t_row,t_for_u,row_type,u_type,b_struct);
   if (t_row & ZERO_MEAN)
      add_value(tGrid,BSC*sum(tGrid,x,t_row,row_type),y,t_row,row_type);
}

INT BS_smoothing_step_ex(tGrid,i,alpha,ZA,ZB,precond_matrix,f,fe,u,ue,
                      d,p,q,de,pe,qe,xe,ye,
                      t_u,t_p,u_type,p_type,A_struct,B_struct,
                      q0,q1,q2,precond_type,multB_X_BT)
GRID *tGrid;
FLOAT alpha;
INT i, ZA, ZB, precond_matrix, f, fe, u, ue, d, p, q, de, pe, qe, xe, ye;
INT t_u, t_p, u_type, p_type, A_struct, B_struct, q0, q1, q2, precond_type;
void (*multB_X_BT)();
{
   INT ii, j, m=15, n;/* GMRES used  =>  qe,...,qe+m+1 are auxiliary variables*/
   FLOAT s;
   INT C_struct=C_STRUCT;

   if (i == 1){ 
      printf(" Smoothing.\n");
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
   }
   else if (i == 2)
      Stokes_defect(tGrid,ZA,ZB,f,fe,u,ue,d,de,q,t_u,t_p,u_type,p_type,A_struct,B_struct,C_struct);
   else if (i == 3)
      printf(" Smoothing.\n");
   BS_PRECOND(d,q)
   mult_B(tGrid,ZB,q,qe,qe,t_p,t_u,p_type,u_type,B_struct);
   mult_and_add(tGrid,-alpha,de,qe,de,t_p,p_type);
   SOLVER_IN_BS_SMOOTHER(multB_ILU_BT_ex)
   mult_BT(tGrid,ZB,pe,p,p,t_p,t_u,p_type,u_type,B_struct);
   subtr(tGrid,d,p,d,t_u,u_type);
   BS_PRECOND(d,q)
   mult_and_add(tGrid,-1./alpha,q,u,u,t_u,u_type);
   subtr(tGrid,ue,pe,ue,t_p,p_type);
   if (t_p & ZERO_MEAN){
      s = sum_n(tGrid,ue,&n,t_p,p_type);
      add_value(tGrid,-s/n,ue,t_p,p_type);
   }
   if (i == 1 || i == 3)
      printf("     %i iterations.\n",j);
   return(j);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

FLOAT first_def(tGrid,alpha,Z,h,f,fe,u,ue,d,p,q,de,pe,qe)
GRID *tGrid;
FLOAT alpha;
INT Z, h, f, fe, u, ue, d, p ,q, de, pe, qe;
{
   vs_set_value_nf(tGrid,0.0,u,1);
   set_value_e(tGrid,0.0,ue);
/*   smoothing_step(tGrid,2,alpha,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);  */
   mult_A(tGrid,Z,u,d,d,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   mult_BT(tGrid,0,ue,q,q,WITHOUT_FIRST,0,Q_SE,Q_VNSF,0);
   vs_add_nf(tGrid,d,q,q,1);
   vs_subtr_nf(tGrid,f,q,q,1);
   return(vs_dot_nf(tGrid,q,q,1));
}    

void quarter_xy(tGrid,x,y,z,fx,fy,fz,Z,h,f,fe,u,ue,d,p,q,de,pe,qe)
GRID *tGrid;
FLOAT *x, *y, *z, *fx, *fy, *fz;  /* x<y, z = (x+y)/2; fz<fx || fz<fy */
INT Z, h, f, fe, u, ue, d, p ,q, de, pe, qe;
{
   FLOAT v, w, fv, fw;
   
   if (*fz > *fx && *fz > *fy)
      eprintf("Error in quarter_xy.\n");
   else if (*fz > *fx || *fz > *fy){
      if(*fz > *fx){     /* minimum on (x,z) */
         *y = *z;
         *fy = *fz;
      }
      else{              /* minimum on (z,y) */
         *x = *z;
         *fx = *fz;
      }
      *z = (*x + *y)/2.0;
      *fz = first_def(tGrid,*z,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
   }
   else{
      v = (*x + *z)/2.0;
      fv = first_def(tGrid,v,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
      if (fv < *fz){     /* minimum on (x,z) */
         *y = *z;
         *fy = *fz;
         *z = v;
         *fz = fv;
      }
      else{
         w = (*z + *y)/2.0;
         fw = first_def(tGrid,w,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
         if (fw < *fz){  /* minimum on (z,y) */
            *x = *z;
            *fx = *fz;
            *z = w;
            *fz = fw;
         }
         else{           /* minimum on (v,w) */
            *x = v;
            *fx = fv;
            *y = w;
            *fy = fw;
         }
      }
   }
}

FLOAT find_alpha(tGrid,Z,h,f,fe,u,ue,d,p,q,de,pe,qe,a_old,eps)
GRID *tGrid;
INT Z, h, f, fe, u, ue, d, p ,q, de, pe, qe;
FLOAT eps, a_old;
{
   FLOAT x, y, z, fx, fy, fz;
   
   if (a_old > 0.){
      z = a_old;
      x = z - eps;
      y = z + eps;
   }
   else{
/*      z = norm_of_A(tGrid,Z,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT); */
      z = 1.;
      x = z*0.1;
   }   
   fz = first_def(tGrid,z,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
   fx = first_def(tGrid,x,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
   if (fx > fz){
      y = z+eps;
      while ( (fy=first_def(tGrid,y,Z,h,f,fe,u,ue,d,p,q,de,pe,qe)) < fz){
         z = y;
         fz = fy;
         y *= 1.5;
      }
   }
   else{
      y = z;
      fy = fz;
      z = x;
      fz = fx;
      x *= 0.5;
      while ( (fx=first_def(tGrid,x,Z,h,f,fe,u,ue,d,p,q,de,pe,qe)) < fz){
         z = x;
         fz = fx;
         x *= 0.5;
      }
   }
   z = (x + y)/2.0;
   fz = first_def(tGrid,z,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
   while ( y-x > eps)
      quarter_xy(tGrid,&x,&y,&z,&fx,&fy,&fz,Z,h,f,fe,u,ue,d,p,q,de,pe,qe);
   if (fx < fz) z = x;
   if (fy < fz) z = y;
   return(z);
} 
   
void set_alpha(mg,Z,h,f,fe,u,ue,d,p,q,de,pe,qe)
MULTIGRID *mg;
INT Z, h, f, fe, u, ue, d, p ,q, de, pe, qe;
{
   GRID *theGrid;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
     theGrid->alpha = find_alpha(theGrid,Z,h,f,fe,u,ue,d,p,q,de,pe,qe,-1.,0.05);
}
   
#endif

/* CAUTION! The rhs b is changed!!! */
void ggem(a,b,x,n)
FLOAT a[][N_GGEM], *b, *x;
INT n;
{
   FLOAT z, max=1.;
   INT p[N_GGEM], i, ii, j, jj, k, n1=n-1;

   if (n > N_GGEM)
     eprintf("Error: N_GGEM too small.\n");
   else{
   for (i = 0; i < n; i++)
      p[i] = i;
   for (k = 0; k < n1 && max > EPSA; k++){
      max = 0.;
      for (i = k; i < n; i++)
         for (j = k; j < n; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            for (j = k; j < n; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(b[k],b[ii],z)
         }
         if (jj != k){
            EXCHANGE(p[k],p[jj],i)
            for (i = 0; i < n; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
         }
         for (i = k+1; i < n; i++){
            z = a[i][k]/a[k][k];
            for (j = k+1; j < n; j++)
               a[i][j] -= z*a[k][j];
            b[i] -= b[k]*z;
         }
      }
   }
   if (max < EPSA || fabs(a[n1][n1]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      for (i = n1; i >= 0; i--){
         z = b[i];
         for (j = i+1; j < n; j++)
            z -= x[p[j]]*a[i][j];
         x[p[i]] = z/a[i][i];
      }
   }
   }
}

/* CAUTION! The rhs c is changed!!! */
void sggem(a,c,x,n,m)
FLOAT a[][N_SGGEM], c[][N_SGGEM], x[][N_SGGEM];
INT n, m;
{
   FLOAT z, max=1.;
   INT p[N_SGGEM], i, ii, j, jj, k, n1=n-1;

   if (n > N_SGGEM)
     eprintf("Error: N_SGGEM too small.\n");
   else{
   for (i = 0; i < n; i++)
      p[i] = i;
   for (k = 0; k < n1 && max > EPSA; k++){
      max = 0.;
      for (i = k; i < n; i++)
         for (j = k; j < n; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            for (j = k; j < n; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            for (j = 0; j < m; j++)
               EXCHANGE(c[j][k],c[j][ii],z)
         }
         if (jj != k){
            EXCHANGE(p[k],p[jj],i)
            for (i = 0; i < n; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
         }
         for (i = k+1; i < n; i++){
            z = a[i][k]/a[k][k];
            for (j = k+1; j < n; j++)
               a[i][j] -= z*a[k][j];
            for (j = 0; j < m; j++)
               c[j][i] -= c[j][k]*z;
         }
      }
   }
   if (max < EPSA || fabs(a[n1][n1]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      for (i = n1; i >= 0; i--)
         for (k = 0; k < m; k++){
            z = c[k][i];
            for (j = i+1; j < n; j++)
               z -= x[k][p[j]]*a[i][j];
            x[k][p[i]] = z/a[i][i];
         }
   }
   }
}

void solve_local_Stokes(n,a,b,u_rhs,g,y,p)
FLOAT a[][N_SGGEM], b[][2], u_rhs[][2], g, y[][N_SGGEM], *p;
INT n;
{
   FLOAT c[4][N_SGGEM], q=0., s=0.;
   INT i;
   	
   for (i = 0; i < n; i++){
      c[0][i] = b[i][0];   	
      c[1][i] = b[i][1];   	
      c[2][i] = u_rhs[i][0];
      c[3][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,4);
   for (i = 0; i < n; i++){
      s += b[i][0]*y[0][i] + b[i][1]*y[1][i];
      q += b[i][0]*y[2][i] + b[i][1]*y[3][i];
   }
   *p = (q - g)/s;
   for (i = 0; i < n; i++){
      y[0][i] = y[2][i] - y[0][i]*(*p);
      y[1][i] = y[3][i] - y[1][i]*(*p);
   }
}

void solve_local_Stokes_with_a1_a2(n,a1,a2,b,u_rhs,g,y,p)
FLOAT a1[][N_SGGEM], a2[][N_SGGEM], b[][2], u_rhs[][2], g, y[][N_SGGEM], *p;
INT n;
{
   FLOAT c[2][N_SGGEM], d[2][N_SGGEM], q=0., s=0.;
   INT i;
   	
   for (i = 0; i < n; i++){
      c[0][i] = b[i][0];   	
      c[1][i] = u_rhs[i][0];
      d[0][i] = b[i][1];   	
      d[1][i] = u_rhs[i][1];
   }
   sggem(a1,c,y,n,2);
   sggem(a2,d,c,n,2);
   for (i = 0; i < n; i++){
      s += b[i][0]*y[0][i] + b[i][1]*c[0][i];
      q += b[i][0]*y[1][i] + b[i][1]*c[1][i];
   }
   *p = (q - g)/s;
   for (i = 0; i < n; i++){
      y[0][i] = y[1][i] - y[0][i]*(*p);
      y[1][i] = c[1][i] - c[0][i]*(*p);
   }
}

void solve_local_Stokes_Korn(n,n2,a,b,u_rhs,g,y,p)
FLOAT a[][N_SGGEM], b[][2], u_rhs[][2], g, y[][N_SGGEM], *p;
INT n, n2;
{
   FLOAT c[2][N_SGGEM], q=0., s=0.;
   INT i, in;
   	
   for (i = 0; i < n; i++){
      in = i + n;
      c[0][i]  = b[i][0];   	
      c[0][in] = b[i][1];   	
      c[1][i]  = u_rhs[i][0];
      c[1][in] = u_rhs[i][1];
   }
   sggem(a,c,y,n2,2);
   for (i = 0; i < n; i++){
      in = i + n;
      s += b[i][0]*y[0][i] + b[i][1]*y[0][in];
      q += b[i][0]*y[1][i] + b[i][1]*y[1][in];
   }
   *p = (q - g)/s;
   for (i = 0; i < n; i++){
      in = i + n;
      y[0][i] = y[1][i]  - y[0][i ]*(*p);
      y[1][i] = y[1][in] - y[0][in]*(*p);
   }
}

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p1nc_p0_old(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   FACE *fa0, *fa1, *fa2;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[10][N_SGGEM], u_rhs[4][2], rhs[10], x[10], q, s;
   INT i, j, m, m1, n=1, nn1, ifa0=0, ifa1=0, ifa2=0;
   	
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   n--;
   m = 2*n;
   m1 = m + 1;
   for (i = 0; i < m1; i++)
      for (j = 0; j < m1; j++)
         a[i][j] = 0.;
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,0,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,0,1);
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,1,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,1,1);
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,2,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,2,1);
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < SIDES; i++){
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }
   rhs[0] = ED(pelem,g);

   nn1 = n + 1;
   for (i = 1; i < nn1; i++){
      for (j = 1; j < nn1; j++)
         a[i+n][j+n] = a[i][j];
      rhs[i]   = u_rhs[i][0];
      rhs[i+n] = u_rhs[i][1];
   }

   ggem(a,rhs,x,m1);

/*
   s = fabs(ED(pelem,p) - x[0]);
   if (IS_FF(fa0) && (q=MAX(fabs(FDV(fa0,u,0)-x[ifa0]),
                            fabs(FDV(fa0,u,1)-x[ifa0+n]))) > s)
      s = q;
   if (IS_FF(fa1) && (q=MAX(fabs(FDV(fa1,u,0)-x[ifa1]),
                            fabs(FDV(fa1,u,1)-x[ifa1+n]))) > s)
      s = q;
   if (IS_FF(fa2) && (q=MAX(fabs(FDV(fa2,u,0)-x[ifa2]),
                            fabs(FDV(fa2,u,1)-x[ifa2+n]))) > s)
      s = q;
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(x[0] - ED(pelem,p));
   if (IS_FF(fa0)){
      FDV(fa0,u_new,0) = FDV(fa0,u,0) + om*(x[ifa0]   - FDV(fa0,u,0));
      FDV(fa0,u_new,1) = FDV(fa0,u,1) + om*(x[ifa0+n] - FDV(fa0,u,1));
   }
   if (IS_FF(fa1)){
      FDV(fa1,u_new,0) = FDV(fa1,u,0) + om*(x[ifa1]   - FDV(fa1,u,0));
      FDV(fa1,u_new,1) = FDV(fa1,u,1) + om*(x[ifa1+n] - FDV(fa1,u,1));
   }
   if (IS_FF(fa2)){
      FDV(fa2,u_new,0) = FDV(fa2,u,0) + om*(x[ifa2]   - FDV(fa2,u,0));
      FDV(fa2,u_new,1) = FDV(fa2,u,1) + om*(x[ifa2+n] - FDV(fa2,u,1));
   }
}

void full_Vanka_step_p1nc_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   FACE *fa0, *fa1, *fa2;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[4][N_SGGEM], b[4][2], c[4][N_SGGEM], y[4][N_SGGEM], u_rhs[4][2], 
         q, s, pressure;
   INT i, j, n=0, ifa0, ifa1, ifa2;

   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         a[i][j] = 0.;
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,0))
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,1))
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,2))
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < SIDES; i++){
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }
   solve_local_Stokes(n,a,b,u_rhs,ED(pelem,g),y,&pressure);
/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
   UPDATE_FDV_VALUES
}

#else

void full_Vanka_step_p1nc_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p1nc_p0 not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p1nc_p0_Korn(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   FACE *fa0, *fa1, *fa2;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[6][N_SGGEM], b[3][2], y[2][N_SGGEM], u_rhs[3][2], 
         q, s, pressure;
   INT i, j, n=0, n2, ifa0=0, ifa1=0, ifa2=0, ifa0n, ifa1n, ifa2n;

   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
// LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   if (IS_FF(fa0)) ifa0 = n++;
   if (IS_FF(fa1)) ifa1 = n++;
   if (IS_FF(fa2)) ifa2 = n++;
   ifa0n = ifa0 + n;
   ifa1n = ifa1 + n;
   ifa2n = ifa2 + n;
   n2 = 2*n;
   for (i = 0; i < n2; i++)
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;
   if (IS_FF(fa0)){
      i  = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,0))
      MFILL_2x2(a,COEFFNNP(fa0,ZA),i,i,ifa0n,ifa0n)
      for (pfl = FSTART(fa0); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa1)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa1,ifa0n,ifa1n)
         else if (NBFACE(pfl) == fa2)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa2,ifa0n,ifa2n)
         else 
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   if (IS_FF(fa1)){
      i  = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,1))
      MFILL_2x2(a,COEFFNNP(fa1,ZA),i,i,ifa1n,ifa1n)
      for (pfl = FSTART(fa1); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa0)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa0,ifa1n,ifa0n)
         else if (NBFACE(pfl) == fa2)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa2,ifa1n,ifa2n)
         else
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   if (IS_FF(fa2)){
      i  = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,2))
      MFILL_2x2(a,COEFFNNP(fa2,ZA),i,i,ifa2n,ifa2n)
      for (pfl = FSTART(fa2); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa0)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa0,ifa2n,ifa0n)
         else if (NBFACE(pfl) == fa1)
            MFILL_2x2(a,COEFFNNP(pfl,ZA),i,ifa1,ifa2n,ifa1n)
         else
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < SIDES; i++){
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }

   solve_local_Stokes_Korn(n,n2,a,b,u_rhs,ED(pelem,g),y,&pressure);
/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
// UPDATE_FDV_VALUES
   if (IS_FF(fa0)){
      FDV(fa0,u_new,0) = FDV(fa0,u,0) + om*(y[0][ifa0] - FDV(fa0,u,0));
      FDV(fa0,u_new,1) = FDV(fa0,u,1) + om*(y[1][ifa0] - FDV(fa0,u,1));
   }
   if (IS_FF(fa1)){
      FDV(fa1,u_new,0) = FDV(fa1,u,0) + om*(y[0][ifa1] - FDV(fa1,u,0));
      FDV(fa1,u_new,1) = FDV(fa1,u,1) + om*(y[1][ifa1] - FDV(fa1,u,1));
   }
   if (IS_FF(fa2)){
      FDV(fa2,u_new,0) = FDV(fa2,u,0) + om*(y[0][ifa2] - FDV(fa2,u,0));
      FDV(fa2,u_new,1) = FDV(fa2,u,1) + om*(y[1][ifa2] - FDV(fa2,u,1));
   }
}

void half_full_Vanka_step_p1nc_p0_Korn(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   FACE *fa0, *fa1, *fa2;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a0[3][N_SGGEM], a1[3][N_SGGEM], b[3][2], y[2][N_SGGEM], u_rhs[3][2], 
         q, s, pressure;
   INT i, j, n=0, ifa0=0, ifa1=0, ifa2=0;

   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
// LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   if (IS_FF(fa0)) ifa0 = n++;
   if (IS_FF(fa1)) ifa1 = n++;
   if (IS_FF(fa2)) ifa2 = n++;
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         a0[i][j] = a1[i][j] = 0.;
   if (IS_FF(fa0)){
      i  = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,0))
      MFILL_SUBTR_2x2(a0,a1,COEFFNNP(fa0,ZA),i,i,u_rhs[i],FDVP(fa0,u))
      for (pfl = FSTART(fa0); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa1)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa1,u_rhs[i],FDVP(fa1,u))
         else if (NBFACE(pfl) == fa2)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa2,u_rhs[i],FDVP(fa2,u))
         else 
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   if (IS_FF(fa1)){
      i  = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,1))
      MFILL_SUBTR_2x2(a0,a1,COEFFNNP(fa1,ZA),i,i,u_rhs[i],FDVP(fa1,u))
      for (pfl = FSTART(fa1); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa0)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa0,u_rhs[i],FDVP(fa0,u))
         else if (NBFACE(pfl) == fa2)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa2,u_rhs[i],FDVP(fa2,u))
         else
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   if (IS_FF(fa2)){
      i  = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      SET1(b[i],COEFF_BDFP(pelem,ZB,2))
      MFILL_SUBTR_2x2(a0,a1,COEFFNNP(fa2,ZA),i,i,u_rhs[i],FDVP(fa2,u))
      for (pfl = FSTART(fa2); pfl; pfl = NEXT(pfl)){
         if (NBFACE(pfl) == fa0)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa0,u_rhs[i],FDVP(fa0,u))
         else if (NBFACE(pfl) == fa1)
            MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pfl,ZA),i,ifa1,u_rhs[i],FDVP(fa1,u))
         else
            MSET9(u_rhs[i],COEFFNNP(pfl,ZA),FDVP(NBFACE(pfl),u))
      }
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < SIDES; i++){
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }

   solve_local_Stokes_with_a1_a2(n,a0,a1,b,u_rhs,ED(pelem,g),y,&pressure);
/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
// UPDATE_FDV_VALUES
   if (IS_FF(fa0)){
      FDV(fa0,u_new,0) = FDV(fa0,u,0) + om*(y[0][ifa0] - FDV(fa0,u,0));
      FDV(fa0,u_new,1) = FDV(fa0,u,1) + om*(y[1][ifa0] - FDV(fa0,u,1));
   }
   if (IS_FF(fa1)){
      FDV(fa1,u_new,0) = FDV(fa1,u,0) + om*(y[0][ifa1] - FDV(fa1,u,0));
      FDV(fa1,u_new,1) = FDV(fa1,u,1) + om*(y[1][ifa1] - FDV(fa1,u,1));
   }
   if (IS_FF(fa2)){
      FDV(fa2,u_new,0) = FDV(fa2,u,0) + om*(y[0][ifa2] - FDV(fa2,u,0));
      FDV(fa2,u_new,1) = FDV(fa2,u,1) + om*(y[1][ifa2] - FDV(fa2,u,1));
   }
}

#else

void full_Vanka_step_p1nc_p0_Korn(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p1nc_p0_Korn not available.\n");  }

void half_full_Vanka_step_p1nc_p0_Korn(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: half_full_Vanka_step_p1nc_p0_Korn not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExE_MATR) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA)

void full_Vanka_step_p1c_p0_div_stab(pelem,ZA,ZB,ZC,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, ZC, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   LINK *pli;
   ELINK *peli;
   FLOAT a[3][N_SGGEM], b0[3], b1[3], u_rhs[3][2], c[4][N_SGGEM], y[4][N_SGGEM],
         q=0., s=0., pressure;
   INT i, j, n=0, in0, in1, in2;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n)
   for (i = 0; i < 3; i++){
      b0[i] = b1[i] = 0.;
      for (j = 0; j < 3; j++)
         a[i][j] = 0.;
   }
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      b0[i] = COEFF_BN(pelem,ZB,0,0);
      b1[i] = COEFF_BN(pelem,ZB,0,1);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      b0[i] = COEFF_BN(pelem,ZB,1,0);
      b1[i] = COEFF_BN(pelem,ZB,1,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      b0[i] = COEFF_BN(pelem,ZB,2,0);
      b1[i] = COEFF_BN(pelem,ZB,2,1);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
   }
   for (peli = pelem->estart; peli; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0)
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n1)
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n2)
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
         }
         if (IS_FF(pel->f[i]) && (pel->f[i] == pelem->f[0] || 
                                  pel->f[i] == pelem->f[1] || 
                                  pel->f[i] == pelem->f[2]))
            q -= COEFF_BF(pel,ZC,i)*ED(pel,p);
      }
   }

   for (i = 0; i < n; i++){
      c[0][i] = b0[i];
      c[1][i] = b1[i];
      c[2][i] = u_rhs[i][0];
      c[3][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,4);
   for (i = 0; i < n; i++){
      s += b0[i]*y[0][i] + b1[i]*y[1][i];
      q += b0[i]*y[2][i] + b1[i]*y[3][i];
   }
   pressure = (q  - ED(pelem,g))/(s + COEFF_EE(pelem,ZC));
   for (i = 0; i < n; i++){
      y[0][i] = y[2][i] - y[0][i]*pressure;
      y[1][i] = y[3][i] - y[1][i]*pressure;
   }

/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_ND_DIFF
   if (fabs(s) < 1.e-13)
      printf("Max. diff. = %e\n",0.);
   else
      printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
   UPDATE_ND_VALUES
}

#else

void full_Vanka_step_p1c_p0_div_stab(pelem,ZA,ZB,ZC,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, ZC, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p1c_p0_div_stab not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (F_DATA & ONE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExDN_MATR) && (E_DATA & ExF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void full_Vanka_step_p1c_new_fbub_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   ELINK *peli;
   FLOAT a[3][N_SGGEM], af[3], b0[3], b1[3], bf[3], u_rhs[3][2], uf_rhs[3],
         c[4][N_SGGEM], y[4][N_SGGEM],
         q, r, s, pressure;
   INT i, j, nn=0, nf=0, in0, in1, in2, ifa0, ifa1, ifa2;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   if (IS_FN(n0)) in0 = nn++;
   if (IS_FN(n1)) in1 = nn++;
   if (IS_FN(n2)) in2 = nn++;
   if (IS_FF(fa0)) ifa0 = nf++;
   if (IS_FF(fa1)) ifa1 = nf++;
   if (IS_FF(fa2)) ifa2 = nf++;
   for (i = 0; i < 3; i++){
      b0[i] = b1[i] = af[i] = 0.;
      for (j = 0; j < 3; j++)
         a[i][j] = 0.;
   }
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      b0[i] = COEFF_BN(pelem,ZB,0,0);
      b1[i] = COEFF_BN(pelem,ZB,0,1);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      b0[i] = COEFF_BN(pelem,ZB,1,0);
      b1[i] = COEFF_BN(pelem,ZB,1,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      b0[i] = COEFF_BN(pelem,ZB,2,0);
      b1[i] = COEFF_BN(pelem,ZB,2,1);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      uf_rhs[i] = FD(fa0,f);
      bf[i] = COEFF_BF(pelem,ZB,0);
      af[i] = COEFF_FF(fa0,ZA);
   }
   if (IS_FF(fa1)){
      i = ifa1;
      uf_rhs[i] = FD(fa1,f);
      bf[i] = COEFF_BF(pelem,ZB,1);
      af[i] = COEFF_FF(fa1,ZA);
   }
   if (IS_FF(fa2)){
      i = ifa2;
      uf_rhs[i] = FD(fa2,f);
      bf[i] = COEFF_BF(pelem,ZB,2);
      af[i] = COEFF_FF(fa2,ZA);
   }
   for (peli = pelem->estart; peli; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0)
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n1)
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n2)
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               uf_rhs[ifa0] -= COEFF_BF(pel,ZB,i)*ED(pel,p);
            else if (pel->f[i] == fa1)
               uf_rhs[ifa1] -= COEFF_BF(pel,ZB,i)*ED(pel,p);
            else if (pel->f[i] == fa2)
               uf_rhs[ifa2] -= COEFF_BF(pel,ZB,i)*ED(pel,p);
         }
      }
   }

   for (i = 0; i < nn; i++){
      c[0][i] = b0[i];
      c[1][i] = b1[i];
      c[2][i] = u_rhs[i][0];
      c[3][i] = u_rhs[i][1];
   }
   sggem(a,c,y,nn,4);
   s = q = 0.;
   for (i = 0; i < nn; i++){
      s += b0[i]*y[0][i] + b1[i]*y[1][i];
      q += b0[i]*y[2][i] + b1[i]*y[3][i];
   }
   for (i = 0; i < nf; i++){
      r  = bf[i]/af[i];
      s += r*bf[i];
      q += r*uf_rhs[i];
   }
   pressure = (q - ED(pelem,g))/s;
   for (i = 0; i < nn; i++){
      y[0][i] = y[2][i] - y[0][i]*pressure;
      y[1][i] = y[3][i] - y[1][i]*pressure;
   }
   for (i = 0; i < nf; i++)
      y[2][i] = (uf_rhs[i] - bf[i]*pressure)/af[i];

/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_ND_DIFF
   if (IS_FF(fa0) && (q=fabs(FD(fa0,u)-y[2][ifa0])) > s)
      s = q;
   if (IS_FF(fa1) && (q=fabs(FD(fa1,u)-y[2][ifa1])) > s)
      s = q;
   if (IS_FF(fa2) && (q=fabs(FD(fa2,u)-y[2][ifa2])) > s)
      s = q;
   if (fabs(s) < 1.e-13)
      printf("Max. diff. = %e\n",0.);
   else
      printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
   UPDATE_ND_VALUES
   if (IS_FF(fa0))
      FD(fa0,u_new) = FD(fa0,u) + om*(y[2][ifa0] - FD(fa0,u));
   if (IS_FF(fa1))
      FD(fa1,u_new) = FD(fa1,u) + om*(y[2][ifa1] - FD(fa1,u));
   if (IS_FF(fa2))
      FD(fa2,u_new) = FD(fa2,u) + om*(y[2][ifa2] - FD(fa2,u));
}

void stab_eq_rhs_to_div_eq_rhs(tGrid,ZA,ZB,f,g,p)
GRID *tGrid;
INT ZA, ZB, f, g, p;
{
   ELEMENT *pel;
   INT i;
   	
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         if (NOT_BF(pel->f[i]))
            ED(pel,g) -= FD(pel->f[i],f)*COEFF_BF(pel,ZB,i)/COEFF_FF(pel->f[i],ZA);
}

void compute_stab_velocity(tGrid,ZA,ZB,u,p,f)
GRID *tGrid;
INT ZA, ZB, u, p, f;
{
   ELEMENT *pel;
   FACE *pface;
   INT i;
   	
   copy(tGrid,f,u,0,Q_SF);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         FD(pel->f[i],u) -= COEFF_BF(pel,ZB,i)*ED(pel,p);
   for (pface = FIRSTF(tGrid); pface; pface = pface->succ)
      if (IS_BF(pface))
         FD(pface,u) = 0.;
      else
         FD(pface,u) /= COEFF_FF(pface,ZA);
}

#else

void full_Vanka_step_p1c_new_fbub_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p1c_new_fbub_p0 not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExDN_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p2_p0_old(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[15][N_SGGEM], u_rhs[7][2], rhs[15], x[15], q, s;
   INT i, j, m, m1, n=0, nn1, in0=0, in1=0, in2=0, ifa0=0, ifa1=0, ifa2=0;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   if (IS_FN(n0)) in0 = ++n;
   if (IS_FN(n1)) in1 = ++n;
   if (IS_FN(n2)) in2 = ++n;
   if (IS_FF(fa0)) ifa0 = ++n;
   if (IS_FF(fa1)) ifa1 = ++n;
   if (IS_FF(fa2)) ifa2 = ++n;
   m = 2*n;
   m1 = m + 1;
   for (i = 0; i < m1; i++)
      for (j = 0; j < m1; j++)
         a[i][j] = 0.;
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,0,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,0,1);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,1,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,1,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,2,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,2,1);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,0,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,0,1);
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,1,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,1,1);
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,2,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,2,1);
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0)
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n1)
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n2)
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }
   rhs[0] = ED(pelem,g);

   nn1 = n + 1;
   for (i = 1; i < nn1; i++){
      for (j = 1; j < nn1; j++)
         a[i+n][j+n] = a[i][j];
      rhs[i]   = u_rhs[i][0];
      rhs[i+n] = u_rhs[i][1];
   }

   ggem(a,rhs,x,m1);

/*
   s = fabs(ED(pelem,p) - x[0]);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-x[in0]),
                           fabs(ND(n0,u,1)-x[in0+n]))) > s)
      s = q;
   if (IS_FN(n1) && (q=MAX(fabs(ND(n1,u,0)-x[in1]),
                           fabs(ND(n1,u,1)-x[in1+n]))) > s)
      s = q;
   if (IS_FN(n2) && (q=MAX(fabs(ND(n2,u,0)-x[in2]),
                           fabs(ND(n2,u,1)-x[in2+n]))) > s)
      s = q;
   if (IS_FF(fa0) && (q=MAX(fabs(FDV(fa0,u,0)-x[ifa0]),
                            fabs(FDV(fa0,u,1)-x[ifa0+n]))) > s)
      s = q;
   if (IS_FF(fa1) && (q=MAX(fabs(FDV(fa1,u,0)-x[ifa1]),
                            fabs(FDV(fa1,u,1)-x[ifa1+n]))) > s)
      s = q;
   if (IS_FF(fa2) && (q=MAX(fabs(FDV(fa2,u,0)-x[ifa2]),
                            fabs(FDV(fa2,u,1)-x[ifa2+n]))) > s)
      s = q;
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(x[0] - ED(pelem,p));
   UPDATE_ND_VALUES
   UPDATE_FDV_VALUES
}

void full_Vanka_step_p2_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[7][N_SGGEM], b0[7], b1[7], c[4][N_SGGEM], y[4][N_SGGEM], u_rhs[7][2],
         q, s, pressure;
   INT i, j, n=0, in0, in1, in2, ifa0, ifa1, ifa2;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n)
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   for (i = 0; i < n; i++){
      b0[i] = b1[i] = 0.;
      for (j = 0; j < n; j++)
         a[i][j] = 0.;
   }
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      b0[i] = COEFF_BN(pelem,ZB,0,0);
      b1[i] = COEFF_BN(pelem,ZB,0,1);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      b0[i] = COEFF_BN(pelem,ZB,1,0);
      b1[i] = COEFF_BN(pelem,ZB,1,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      b0[i] = COEFF_BN(pelem,ZB,2,0);
      b1[i] = COEFF_BN(pelem,ZB,2,1);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      b0[i] = COEFF_BDF(pelem,ZB,0,0);
      b1[i] = COEFF_BDF(pelem,ZB,0,1);
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      b0[i] = COEFF_BDF(pelem,ZB,1,0);
      b1[i] = COEFF_BDF(pelem,ZB,1,1);
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      b0[i] = COEFF_BDF(pelem,ZB,2,0);
      b1[i] = COEFF_BDF(pelem,ZB,2,1);
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0)
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n1)
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
            else if (pel->n[i] == n2)
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa1)
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
            else if (pel->f[i] == fa2)
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
         }
      }
   }

   for (i = 0; i < n; i++){
      c[0][i] = b0[i];
      c[1][i] = b1[i];
      c[2][i] = u_rhs[i][0];
      c[3][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,4);
   s = q = 0.;
   for (i = 0; i < n; i++){
      s += b0[i]*y[0][i] + b1[i]*y[1][i];
      q += b0[i]*y[2][i] + b1[i]*y[3][i];
   }
   pressure = (q - ED(pelem,g))/s;
   for (i = 0; i < n; i++){
      y[0][i] = y[2][i] - y[0][i]*pressure;
      y[1][i] = y[3][i] - y[1][i]*pressure;
   }

/*
   s = fabs(ED(pelem,p) - pressure);
   COMPUTE_MAX_ND_DIFF
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
   UPDATE_ND_VALUES
   UPDATE_FDV_VALUES
}

#else

void full_Vanka_step_p2_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p0 not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExDN_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (E_DATA & VECTOR_ELEMENT_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p2b_p0_old(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[17][N_SGGEM], u_rhs[8][2], rhs[17], x[17], q, s;
   INT i, j, m, m1, n=1, nn1, in0=0, in1=0, in2=0, ifa0=0, ifa1=0, ifa2=0;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   if (IS_FN(n0)) in0 = ++n;
   if (IS_FN(n1)) in1 = ++n;
   if (IS_FN(n2)) in2 = ++n;
   if (IS_FF(fa0)) ifa0 = ++n;
   if (IS_FF(fa1)) ifa1 = ++n;
   if (IS_FF(fa2)) ifa2 = ++n;
   m = 2*n;
   m1 = m + 1;
   for (i = 0; i < m1; i++)
      for (j = 0; j < m1; j++)
         a[i][j] = 0.;
   a[1][1] = COEFF_EE(pelem,ZA);
   SET1(u_rhs[1],EDVP(pelem,f))
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,0,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,0,1);
      a[1][i] = COEFF_EN(pelem,ZA,0);
      a[i][1] = COEFF_NE(pelem,ZA,0);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,1,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,1,1);
      a[1][i] = COEFF_EN(pelem,ZA,1);
      a[i][1] = COEFF_NE(pelem,ZA,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      a[0][i]   = a[i][0]   = COEFF_BN(pelem,ZB,2,0);
      a[0][i+n] = a[i+n][0] = COEFF_BN(pelem,ZB,2,1);
      a[1][i] = COEFF_EN(pelem,ZA,2);
      a[i][1] = COEFF_NE(pelem,ZA,2);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,0,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,0,1);
      a[1][i] = COEFF_EF(pelem,ZA,0);
      a[i][1] = COEFF_FE(pelem,ZA,0);
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,1,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,1,1);
      a[1][i] = COEFF_EF(pelem,ZA,1);
      a[i][1] = COEFF_FE(pelem,ZA,1);
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      a[0][i]   = a[i][0]   = COEFF_BDF(pelem,ZB,2,0);
      a[0][i+n] = a[i+n][0] = COEFF_BDF(pelem,ZB,2,1);
      a[1][i] = COEFF_EF(pelem,ZA,2);
      a[i][1] = COEFF_FE(pelem,ZA,2);
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0){
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in0],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n1){
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in1],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n2){
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in2],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0){
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa0],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa1){
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa1],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa2){
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa2],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
         }
      }
   }
   rhs[0] = ED(pelem,g);

   nn1 = n + 1;
   for (i = 1; i < nn1; i++){
      for (j = 1; j < nn1; j++)
         a[i+n][j+n] = a[i][j];
      rhs[i]   = u_rhs[i][0];
      rhs[i+n] = u_rhs[i][1];
   }

   ggem(a,rhs,x,m1);

/*
   s = MAX(fabs(EDV(pelem,u,0) - x[1]),fabs(EDV(pelem,u,1) - x[n+1]));
   s = MAX(s,fabs(ED(pelem,p) - x[0]));
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-x[in0]),
                           fabs(ND(n0,u,1)-x[in0+n]))) > s)
      s = q;
   if (IS_FN(n1) && (q=MAX(fabs(ND(n1,u,0)-x[in1]),
                           fabs(ND(n1,u,1)-x[in1+n]))) > s)
      s = q;
   if (IS_FN(n2) && (q=MAX(fabs(ND(n2,u,0)-x[in2]),
                           fabs(ND(n2,u,1)-x[in2+n]))) > s)
      s = q;
   if (IS_FF(fa0) && (q=MAX(fabs(FDV(fa0,u,0)-x[ifa0]),
                            fabs(FDV(fa0,u,1)-x[ifa0+n]))) > s)
      s = q;
   if (IS_FF(fa1) && (q=MAX(fabs(FDV(fa1,u,0)-x[ifa1]),
                            fabs(FDV(fa1,u,1)-x[ifa1+n]))) > s)
      s = q;
   if (IS_FF(fa2) && (q=MAX(fabs(FDV(fa2,u,0)-x[ifa2]),
                            fabs(FDV(fa2,u,1)-x[ifa2+n]))) > s)
      s = q;
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(x[0] - ED(pelem,p));
   EDV(pelem,u_new,0) = EDV(pelem,u,0) + om*(x[1]   - EDV(pelem,u,0));
   EDV(pelem,u_new,1) = EDV(pelem,u,1) + om*(x[1+n] - EDV(pelem,u,1));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(x[in0]   - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(x[in0+n] - ND(n0,u,1));
   }
   if (IS_FN(n1)){
      ND(n1,u_new,0) = ND(n1,u,0) + om*(x[in1]   - ND(n1,u,0));
      ND(n1,u_new,1) = ND(n1,u,1) + om*(x[in1+n] - ND(n1,u,1));
   }
   if (IS_FN(n2)){
      ND(n2,u_new,0) = ND(n2,u,0) + om*(x[in2]   - ND(n2,u,0));
      ND(n2,u_new,1) = ND(n2,u,1) + om*(x[in2+n] - ND(n2,u,1));
   }
   if (IS_FF(fa0)){
      FDV(fa0,u_new,0) = FDV(fa0,u,0) + om*(x[ifa0]   - FDV(fa0,u,0));
      FDV(fa0,u_new,1) = FDV(fa0,u,1) + om*(x[ifa0+n] - FDV(fa0,u,1));
   }
   if (IS_FF(fa1)){
      FDV(fa1,u_new,0) = FDV(fa1,u,0) + om*(x[ifa1]   - FDV(fa1,u,0));
      FDV(fa1,u_new,1) = FDV(fa1,u,1) + om*(x[ifa1+n] - FDV(fa1,u,1));
   }
   if (IS_FF(fa2)){
      FDV(fa2,u_new,0) = FDV(fa2,u,0) + om*(x[ifa2]   - FDV(fa2,u,0));
      FDV(fa2,u_new,1) = FDV(fa2,u,1) + om*(x[ifa2+n] - FDV(fa2,u,1));
   }
}

void full_Vanka_step_p2b_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[8][N_SGGEM], b0[8], b1[8], c[4][N_SGGEM], y[4][N_SGGEM], u_rhs[8][2],
         q, s, pressure;
   INT i, j, n=1, in0, in1, in2, ifa0, ifa1, ifa2;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n)
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   for (i = 0; i < n; i++){
      b0[i] = b1[i] = 0.;
      for (j = 0; j < n; j++)
         a[i][j] = 0.;
   }
   a[0][0] = COEFF_EE(pelem,ZA);
   SET1(u_rhs[0],EDVP(pelem,f))
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      b0[i] = COEFF_BN(pelem,ZB,0,0);
      b1[i] = COEFF_BN(pelem,ZB,0,1);
      a[0][i] = COEFF_EN(pelem,ZA,0);
      a[i][0] = COEFF_NE(pelem,ZA,0);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      b0[i] = COEFF_BN(pelem,ZB,1,0);
      b1[i] = COEFF_BN(pelem,ZB,1,1);
      a[0][i] = COEFF_EN(pelem,ZA,1);
      a[i][0] = COEFF_NE(pelem,ZA,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      b0[i] = COEFF_BN(pelem,ZB,2,0);
      b1[i] = COEFF_BN(pelem,ZB,2,1);
      a[0][i] = COEFF_EN(pelem,ZA,2);
      a[i][0] = COEFF_NE(pelem,ZA,2);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      b0[i] = COEFF_BDF(pelem,ZB,0,0);
      b1[i] = COEFF_BDF(pelem,ZB,0,1);
      a[0][i] = COEFF_EF(pelem,ZA,0);
      a[i][0] = COEFF_FE(pelem,ZA,0);
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      b0[i] = COEFF_BDF(pelem,ZB,1,0);
      b1[i] = COEFF_BDF(pelem,ZB,1,1);
      a[0][i] = COEFF_EF(pelem,ZA,1);
      a[i][0] = COEFF_FE(pelem,ZA,1);
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      b0[i] = COEFF_BDF(pelem,ZB,2,0);
      b1[i] = COEFF_BDF(pelem,ZB,2,1);
      a[0][i] = COEFF_EF(pelem,ZA,2);
      a[i][0] = COEFF_FE(pelem,ZA,2);
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0){
               SET9(u_rhs[in0],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in0],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n1){
               SET9(u_rhs[in1],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in1],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n2){
               SET9(u_rhs[in2],COEFF_BNP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[in2],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0){
               SET9(u_rhs[ifa0],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa0],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa1){
               SET9(u_rhs[ifa1],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa1],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa2){
               SET9(u_rhs[ifa2],COEFF_BDFP(pel,ZB,i),ED(pel,p))
               SET9(u_rhs[ifa2],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
         }
      }
   }

   for (i = 0; i < n; i++){
      c[0][i] = b0[i];
      c[1][i] = b1[i];
      c[2][i] = u_rhs[i][0];
      c[3][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,4);
   s = q = 0.;
   for (i = 0; i < n; i++){
      s += b0[i]*y[0][i] + b1[i]*y[1][i];
      q += b0[i]*y[2][i] + b1[i]*y[3][i];
   }
   pressure = (q - ED(pelem,g))/s;
   for (i = 0; i < n; i++){
      y[0][i] = y[2][i] - y[0][i]*pressure;
      y[1][i] = y[3][i] - y[1][i]*pressure;
   }

/*
   s = MAX(fabs(EDV(pelem,u,0) - y[0][0]),fabs(EDV(pelem,u,1) - y[1][0]));
   s = MAX(s,fabs(ED(pelem,p) - pressure));
   COMPUTE_MAX_ND_DIFF
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   ED(pelem,p_new) = ED(pelem,p) + om*(pressure - ED(pelem,p));
   EDV(pelem,u_new,0) = EDV(pelem,u,0) + om*(y[0][0] - EDV(pelem,u,0));
   EDV(pelem,u_new,1) = EDV(pelem,u,1) + om*(y[1][0] - EDV(pelem,u,1));
   UPDATE_ND_VALUES
   UPDATE_FDV_VALUES
}

#else

void full_Vanka_step_p2b_p0(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2b_p0 not available.\n");  }

#endif

#if (N_DATA & IxD_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA)

void subtract_Bni_x_sn(rhs,ni,n0,Z,p)
FLOAT rhs[2];
NODE *ni, *n0;
INT Z, p;
{
   NODE *pn;
   LINK *pli, *pli2;

   SET9(rhs,COEFFBP(ni,Z),NDS(ni,p))
   for (pli = TSTART(ni); pli; pli = NEXT(pli)){
      pn=NBNODE(pli);
      if (pn != n0){
         for (pli2 = START(pn); ni != NBNODE(pli2); pli2 = NEXT(pli2));
         SET9(rhs,COEFFBLP(pli2,Z),NDS(pn,p))
      }
   }
}

#else

void subtract_Bni_x_sn(rhs,ni,n0,Z,p)
FLOAT rhs[2]; NODE *ni, *n0; INT Z, p;
{  eprintf("Error: subtract_Bni_x_sn not available.\n");  }

#endif

#if (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA)

void subtract_Bfi_x_sn(rhs,fi,n0,Z,p)
FLOAT rhs[2];
FACE *fi;
NODE *n0;
INT Z, p;
{
   NODE *pn;
   FNLINK *pfn;
   NFLINK *pnfl;

   for (pfn = TFNSTART(fi); pfn; pfn = NEXT(pfn)){
      pn = NBNODE(pfn);
      if (pn != n0){
         for (pnfl = NFSTART(pn); fi != NBFACE(pnfl); pnfl = NEXT(pnfl));
         SET9(rhs,COEFF_NFP(pnfl,Z),NDS(pn,p))
      }
   }
}

#else

void subtract_Bfi_x_sn(rhs,fi,n0,Z,p)
FLOAT rhs[2]; FACE *fi; NODE *n0; INT Z, p;
{  eprintf("Error: subtract_Bfi_x_sn not available.\n");  }

#endif

#if DATA_S & SPECIAL_NODES_AND_FACES

void subtract_Dni_x_ssn(rhs,ni,n0,Z,p,l)
FLOAT rhs[2];
NODE *ni, *n0;
INT Z, p, l;
{
   NODE *pn;

   if (IS_LN(ni)){
      if (l)
         SET9(rhs,ni->s_node->ann0[Z],SNDS(ni->s_node,p))
      pn = ni->s_node->n1;
      if (IS_LN(pn) && pn != n0){
         if (pn->s_node->n1 == ni)
            SET9(rhs,pn->s_node->ann1[Z],SNDS(pn->s_node,p))
         else
            SET9(rhs,pn->s_node->ann2[Z],SNDS(pn->s_node,p))
      }
      pn = ni->s_node->n2;
      if (IS_LN(pn) && pn != n0 && ni->s_node->pel2){
         if (pn->s_node->n1 == ni)
            SET9(rhs,pn->s_node->ann1[Z],SNDS(pn->s_node,p))
         else
            SET9(rhs,pn->s_node->ann2[Z],SNDS(pn->s_node,p))
      }
   }
}

void subtract_Dni_x_ssf(rhs,ni,f0,Z,p)
FLOAT rhs[2];
NODE *ni;
FACE *f0;
INT Z, p;
{
   FACE *pf;

   if (IS_LN(ni)){
      pf = ni->s_node->f1;
      if (IS_LF(pf) && pf != f0){
         if (pf->s_face->n1 == ni)
            SET9(rhs,pf->s_face->afn1[Z],SFDS(pf->s_face,p))
         else
            SET9(rhs,pf->s_face->afn2[Z],SFDS(pf->s_face,p))
      }
      pf = ni->s_node->f2;
      if (IS_LF(pf) && pf != f0 && ni->s_node->pel2){
         if (pf->s_face->n1 == ni)
            SET9(rhs,pf->s_face->afn1[Z],SFDS(pf->s_face,p))
         else
            SET9(rhs,pf->s_face->afn2[Z],SFDS(pf->s_face,p))
      }
   }
}

void subtract_Dfi_x_ssn(rhs,fi,n0,Z,p)
FLOAT rhs[2];
FACE *fi;
NODE *n0;
INT Z, p;
{
   NODE *pn;

   if (IS_LF(fi)){
      pn = fi->s_face->n1;
      if (IS_LN(pn) && pn != n0){
         if (pn->s_node->f1 == fi)
            SET9(rhs,pn->s_node->anf1[Z],SNDS(pn->s_node,p))
         else
            SET9(rhs,pn->s_node->anf2[Z],SNDS(pn->s_node,p))
      }
      pn = fi->s_face->n2;
      if (IS_LN(pn) && pn != n0){
         if (pn->s_node->f1 == fi)
            SET9(rhs,pn->s_node->anf1[Z],SNDS(pn->s_node,p))
         else
            SET9(rhs,pn->s_node->anf2[Z],SNDS(pn->s_node,p))
      }
   }
}

void subtract_Dfi_x_ssf(rhs,fi,Z,p)
FLOAT rhs[2];
FACE *fi;
INT Z, p;
{
   if (IS_LF(fi))
      SET9(rhs,fi->s_face->aff[Z],SFDS(fi->s_face,p))
}

#else

void subtract_Dni_x_ssn(rhs,ni,n0,Z,p,l)
FLOAT rhs[2]; NODE *ni, *n0; INT Z, p, l;
{  eprintf("Error: subtract_Dni_x_ssn not available.\n");  }

void subtract_Dni_x_ssf(rhs,ni,f0,Z,p)
FLOAT rhs[2]; NODE *ni; FACE *f0; INT Z, p;
{  eprintf("Error: subtract_Dni_x_ssf not available.\n");  }

void subtract_Dfi_x_ssn(rhs,fi,n0,Z,p)
FLOAT rhs[2]; FACE *fi; NODE *n0; INT Z, p;
{  eprintf("Error: subtract_Dfi_x_ssn not available.\n");  }

void subtract_Dfi_x_ssf(rhs,fi,Z,p)
FLOAT rhs[2]; FACE *fi; INT Z, p;
{  eprintf("Error: subtract_Dfi_x_ssf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void diag_Vanka_step_p2_p1(n0,ZA,ZB,f,g,u,p,u_new,p_new)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT uu[N_SGGEM][2], u_rhs[N_SGGEM][2], b[N_SGGEM][2], d[N_SGGEM], 
         pp, q0, q1, s;
   INT i, n=0, n1;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      IFLAG(pnode=NBNODE(pl)) = ++n;
      SET1(u_rhs[n],NDD(pnode,f))
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      IFLAG(pface=NBFACE(pnfl)) = ++n;
      SET1(u_rhs[n],FDVP(pface,f))
   }
   if (IS_FN(n0)){
      IFLAG(n0) = ++n;
      SET1(u_rhs[n],NDD(n0,f))
   }
   n1 = n + 1;

   if (IS_FN(n0)){
      SET1(b[n],COEFFBP(n0,ZB))
      d[n] = COEFFN(n0,ZA);
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[n],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode);
         if (IS_FN(n0))
            SET9(u_rhs[n],NDD(pnode,u),COEFFL(pl,ZA))
         d[i] = COEFFN(pnode,ZA);
         SET1(b[i],COEFFBLP(pl,ZB))
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn))
               SET9(u_rhs[i],NDD(pn,u),COEFFL(pli,ZA))
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IS_FF(pf=NBFACE(pnfl)))
               SET9(u_rhs[i],FDVP(pf,u),COEFFL(pnfl,ZA))
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf))
      if (IS_FF(pface=NBFACE(pnf))){
         i = IFLAG(pface);
         if (IS_FN(n0))
            SET9(u_rhs[n],FDVP(pface,u),COEFFL(pnf,ZA))
         d[i] = COEFF_FF(pface,ZA);
         SET1(b[i],COEFF_NFP(pnf,ZB))
         for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
            pn = NBNODE(pfn);
            if (pn != n0){
               for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
               SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
            }
            if (IS_FN(pn))
               SET9(u_rhs[i],NDD(pn,u),COEFFL(pfn,ZA))
         }
         for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
            if (IS_FF(pf=NBFACE(pfl)))
               SET9(u_rhs[i],FDVP(pf,u),COEFFL(pfl,ZA))
      } 
   pp = -NDS(n0,g);
   s = 0;
   for (i = 1; i < n1; i++){
      q0 = b[i][0]/d[i];
      q1 = b[i][1]/d[i];
      pp += u_rhs[i][0]*q0 + u_rhs[i][1]*q1;
      s += b[i][0]*q0 + b[i][1]*q1;
   }
   pp /= s;
   for (i = 1; i < n1; i++){
      uu[i][0] = (u_rhs[i][0] - pp*b[i][0])/d[i];
      uu[i][1] = (u_rhs[i][1] - pp*b[i][1])/d[i];
   }

/*
   s = fabs(NDS(n0,p) - pp);
   if (IS_FN(n0) && (q0=MAX_DIFF(NDD(n0,u),uu[n])) > s)
      s = q0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) &&
                               (q0=MAX_DIFF(NDD(pnode,u),uu[IFLAG(pnode)])) > s)
         s = q0;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) && 
                              (q0=MAX_DIFF(FDVP(pface,u),uu[IFLAG(pface)])) > s)
         s = q0;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = pp;
   if (IS_FN(n0))
      SET1(NDD(n0,u_new),uu[n])
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         SET1(NDD(pnode,u_new),uu[IFLAG(pnode)])
         IFLAG(pnode) = 0;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         SET1(FDVP(pface,u_new),uu[IFLAG(pface)])
         IFLAG(pface) = 0;
      }
}

void full_Vanka_step_p2_p1_old(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[N_SGGEM][N_SGGEM], u_rhs[N_SGGEM][2], rhs[N_SGGEM], x[N_SGGEM], q, s;
   INT i, j, m, m1, n=0, n1;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         IFLAG(pnode) = ++n;
         SET1(u_rhs[n],NDD(pnode,f))
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         IFLAG(pface) = ++n;
         SET1(u_rhs[n],FDVP(pface,f))
      }
   if (IS_FN(n0)){
      IFLAG(n0) = ++n;
      SET1(u_rhs[n],NDD(n0,f))
   }
   m = 2*n;
   m1 = m + 1;
   for (i = 0; i < m1; i++)
      for (j = 0; j < m1; j++)
         a[i][j] = 0.;

   rhs[0] = NDS(n0,g);
   if (IS_FN(n0)){
      a[0][n] = a[n][0] = COEFFB(n0,ZB,0);
      a[0][m] = a[m][0] = COEFFB(n0,ZB,1);
      a[n][n] = COEFFN(n0,ZA);
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[n],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode);
         if (IS_FN(n0))
            a[n][i] = COEFFL(pl,ZA); 
         a[i][i] = COEFFN(pnode,ZA);
         a[0][i] = a[i][0] = COEFFBL(pl,ZB,0);
         a[0][i+n] = a[i+n][0] = COEFFBL(pl,ZB,1);
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn)){
               if (IFLAG(pn))
                  a[i][IFLAG(pn)] = COEFFL(pli,ZA);
               else
                  SET9(u_rhs[i],NDD(pn,u),COEFFL(pli,ZA))
            }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IS_FF(pf=NBFACE(pnfl))){
               if (IFLAG(pf))
                  a[i][IFLAG(pf)] = COEFFL(pnfl,ZA);
               else
                  SET9(u_rhs[i],FDVP(pf,u),COEFFL(pnfl,ZA))
            } 
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf))
      if (IS_FF(pface=NBFACE(pnf))){
         i = IFLAG(pface);
         if (IS_FN(n0))
            a[n][i] = COEFFL(pnf,ZA); 
         a[i][i] = COEFF_FF(pface,ZA);
         a[0][i] = a[i][0] = COEFF_NF(pnf,ZB,0);
         a[0][i+n] = a[i+n][0] = COEFF_NF(pnf,ZB,1);
         for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
            pn = NBNODE(pfn);
            if (pn != n0){
               for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
               SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
            }
            if (IS_FN(pn)){
               if (IFLAG(pn))
                  a[i][IFLAG(pn)] = COEFFL(pfn,ZA);
               else
                  SET9(u_rhs[i],NDD(pn,u),COEFFL(pfn,ZA))
            }
         }
         for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
            if (IS_FF(pf=NBFACE(pfl))){
               if (IFLAG(pf))
                  a[i][IFLAG(pf)] = COEFFL(pfl,ZA);
               else
                  SET9(u_rhs[i],FDVP(pf,u),COEFFL(pfl,ZA))
            } 
      } 
   n1 = n + 1;
   for (i = 1; i < n1; i++){
      for (j = 1; j < n1; j++)
         a[i+n][j+n] = a[i][j];
      rhs[i] = u_rhs[i][0];
      rhs[i+n] = u_rhs[i][1];
   }

   ggem(a,rhs,x,m1);

/*
   s = fabs(NDS(n0,p) - x[0]);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-x[n]),fabs(ND(n0,u,1)-x[m]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-x[IFLAG(pnode)]),
                                    fabs(ND(pnode,u,1)-x[IFLAG(pnode)+n]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-x[IFLAG(pface)]),
                                   fabs(FDV(pface,u,1)-x[IFLAG(pface)+n]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(x[0] - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(x[n] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(x[m] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         ND(pnode,u_new,0) = ND(pnode,u,0) + om*(x[IFLAG(pnode)]-ND(pnode,u,0));
         ND(pnode,u_new,1) = ND(pnode,u,1) + om*(x[IFLAG(pnode)+n]-ND(pnode,u,1));
         IFLAG(pnode) = 0;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         FDV(pface,u_new,0) = FDV(pface,u,0) + om*(x[IFLAG(pface)]-FDV(pface,u,0));
         FDV(pface,u_new,1) = FDV(pface,u,1) + om*(x[IFLAG(pface)+n]-FDV(pface,u,1));
         IFLAG(pface) = 0;
      }
}

void full_Vanka_step_p2_p1(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new, lagr;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[N_SGGEM][N_SGGEM], b[N_SGGEM][2], y[4][N_SGGEM], u_rhs[N_SGGEM][2], 
         q, s, pressure;
   INT i, j, n=0, k;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      SET1(u_rhs[n],NDD(pnode,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],pnode,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],pnode,NULL,ZB,p);
      }
      IFLAG(pnode) = ++n;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      SET1(u_rhs[n],FDVP(pface,f))
      if (lagr == YES){
         subtract_Dfi_x_ssn(u_rhs[n],pface,NULL,ZB,p);
         subtract_Dfi_x_ssf(u_rhs[n],pface,ZB,p);
      }
      IFLAG(pface) = ++n;
   }
   if (IS_FN(n0)){
      SET1(u_rhs[n],NDD(n0,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],n0,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],n0,NULL,ZB,p);
      }
      IFLAG(n0) = ++n;
   }
   k = n - 1;
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         a[i][j] = 0.;

   if (IS_FN(n0)){
      SET1(b[k],COEFFBP(n0,ZB))
      a[k][k] = COEFFN(n0,ZA);
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode) - 1;
         if (IS_FN(n0))
            a[k][i] = COEFFL(pl,ZA); 
         a[i][i] = COEFFN(pnode,ZA);
         SET1(b[i],COEFFBLP(pl,ZB))
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn)){
               if (IFLAG(pn))
                  a[i][IFLAG(pn)-1] = COEFFL(pli,ZA);
               else
                  SET9(u_rhs[i],NDD(pn,u),COEFFL(pli,ZA))
            }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IFLAG(pf=NBFACE(pnfl)))
               a[i][IFLAG(pf)-1] = COEFFL(pnfl,ZA);
            else
               SET9(u_rhs[i],FDVP(pf,u),COEFFL(pnfl,ZA))
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf)){
      i = IFLAG(pface=NBFACE(pnf))-1;
      if (IS_FN(n0))
         a[k][i] = COEFFL(pnf,ZA); 
      a[i][i] = COEFF_FF(pface,ZA);
      SET1(b[i],COEFF_NFP(pnf,ZB))
      for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
         pn = NBNODE(pfn);
         if (pn != n0){
            for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
            SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
         }
         if (IS_FN(pn)){
            if (IFLAG(pn))
               a[i][IFLAG(pn)-1] = COEFFL(pfn,ZA);
            else
               SET9(u_rhs[i],NDD(pn,u),COEFFL(pfn,ZA))
         }
      }
      for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
         if (IFLAG(pf=NBFACE(pfl)))
            a[i][IFLAG(pf)-1] = COEFFL(pfl,ZA);
         else
            SET9(u_rhs[i],FDVP(pf,u),COEFFL(pfl,ZA))
   } 

   solve_local_Stokes(n,a,b,u_rhs,NDS(n0,g),y,&pressure);

/*
   s = fabs(NDS(n0,p) - pressure);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(pressure - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
      ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
      IFLAG(pnode) = 0;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
      FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
      IFLAG(pface) = 0;
   }
}

#if DATA_S & SPECIAL_NODES_AND_FACES

void A_contribution_for_p2_p2ln(sn,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn;
NODE *ni[3];
FACE *fi[2];
FLOAT a[][N_SGGEM], u_rhs[][2];
INT n_ni, n_fi, ZA, f, u;
{
   NODE *pn;
   FACE *pf;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   INT i, l;
   	
   for (l = 0; l < n_ni; l++){
      i = IFLAG(ni[l]) - 1;
      a[i][i] = COEFFN(ni[l],ZA);
      for (pli = START(ni[l]); pli; pli = NEXT(pli)){
         pn=NBNODE(pli);
         if (IS_LN(pn) && (pn == sn->n0 || pn == sn->n1 || pn == sn->n2))
            a[i][IFLAG(pn)-1] = COEFFL(pli,ZA);
         else
            SET9(u_rhs[i],NDD(pn,u),COEFFL(pli,ZA))
      }
      for (pnfl = NFSTART(ni[l]); pnfl; pnfl = NEXT(pnfl)){
         pf=NBFACE(pnfl);
         if (IS_LF(pf) && (pf == sn->f1 || pf == sn->f2))
            a[i][IFLAG(pf)-1] = COEFFL(pnfl,ZA);
         else
            SET9(u_rhs[i],FDVP(pf,u),COEFFL(pnfl,ZA))
      }
   }
   for (l = 0; l < n_fi; l++){
      i = IFLAG(fi[l])-1;
      a[i][i] = COEFF_FF(fi[l],ZA);
      for (pfn = FNSTART(fi[l]); pfn; pfn = NEXT(pfn)){
         pn = NBNODE(pfn);
         if (IS_LN(pn) && (pn == sn->n0 || pn == sn->n1 || pn == sn->n2))
            a[i][IFLAG(pn)-1] = COEFFL(pfn,ZA);
         else
            SET9(u_rhs[i],NDD(pn,u),COEFFL(pfn,ZA))
      }
      for (pfl = FSTART(fi[l]); pfl; pfl = NEXT(pfl)){
         pf=NBFACE(pfl);
         if (IS_LF(pf) && (pf == sn->f1 || pf == sn->f2))
            a[i][IFLAG(pf)-1] = COEFFL(pfl,ZA);
         else
            SET9(u_rhs[i],FDVP(pf,u),COEFFL(pfl,ZA))
      }
   } 
}

#else

void A_contribution_for_p2_p2ln(sn,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn; NODE *ni[3]; FACE *fi[2]; FLOAT a[][N_SGGEM], u_rhs[][2]; INT n_ni, n_fi, ZA, f, u;
{  eprintf("Error: A_contribution_for_p2_p2ln not available.\n");  }

#endif

void A_contribution_for_p2_p2lf(f0,n1,n2,a,u_rhs,ZA,f,u)
FACE *f0;
NODE *n1, *n2;
FLOAT a[][N_SGGEM], u_rhs[][2];
INT ZA, f, u;
{
   NODE *pn;
   FACE *pf;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   	
   a[0][0] = COEFF_FF(f0,ZA);
   for (pfn = FNSTART(f0); pfn; pfn = NEXT(pfn)){
      pn = NBNODE(pfn);
      if (pn == n1)
         a[0][1] = COEFFL(pfn,ZA);
      else if (pn == n2)
         a[0][2] = COEFFL(pfn,ZA);
      else
         SET9(u_rhs[0],NDD(pn,u),COEFFL(pfn,ZA))
   }
   for (pfl = FSTART(f0); pfl; pfl = NEXT(pfl))
      SET9(u_rhs[0],FDVP(NBFACE(pfl),u),COEFFL(pfl,ZA))

   a[1][1] = COEFFN(n1,ZA);
   for (pli = START(n1); pli; pli = NEXT(pli))
      if ((pn=NBNODE(pli)) == n2)
         a[1][2] = COEFFL(pli,ZA);
      else
         SET9(u_rhs[1],NDD(pn,u),COEFFL(pli,ZA))
   for (pnfl = NFSTART(n1); pnfl; pnfl = NEXT(pnfl))
      if ((pf=NBFACE(pnfl)) == f0)
         a[1][0] = COEFFL(pnfl,ZA);
      else
         SET9(u_rhs[1],FDVP(pf,u),COEFFL(pnfl,ZA))
   if (n2){
      a[2][2] = COEFFN(n2,ZA);
      for (pli = START(n2); pli; pli = NEXT(pli))
         if ((pn=NBNODE(pli)) == n1)
            a[2][1] = COEFFL(pli,ZA);
         else
            SET9(u_rhs[2],NDD(pn,u),COEFFL(pli,ZA))
      for (pnfl = NFSTART(n2); pnfl; pnfl = NEXT(pnfl))
         if ((pf=NBFACE(pnfl)) == f0)
            a[2][0] = COEFFL(pnfl,ZA);
         else
            SET9(u_rhs[2],FDVP(pf,u),COEFFL(pnfl,ZA))
   }
}

#else

void diag_Vanka_step_p2_p1(n0,ZA,ZB,f,g,u,p,u_new,p_new)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new;
{  eprintf("Error: diag_Vanka_step_p2_p1 not available.\n");  }

void full_Vanka_step_p2_p1(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new, lagr; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1 not available.\n");  }

void A_contribution_for_p2_p2ln(sn,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn; NODE *ni[3]; FACE *fi[2]; FLOAT a[][N_SGGEM], u_rhs[][2]; INT n_ni, n_fi, ZA, f, u;
{  eprintf("Error: A_contribution_for_p2_p2ln not available.\n");  }

void A_contribution_for_p2_p2lf(f0,n1,n2,a,u_rhs,ZA,f,u)
FACE *f0; NODE *n1, *n2; FLOAT a[][N_SGGEM], u_rhs[][2]; INT ZA, f, u;
{  eprintf("Error: A_contribution_for_p2_p2lf not available.\n");  }

#endif

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p2_p1_Korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new, lagr;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[N_SGGEM][N_SGGEM], b[N_SGGEM][2], y[2][N_SGGEM], u_rhs[N_SGGEM][2], 
         q, s, pressure;
   INT i, in, j, jn, n=0, n2, k, k2;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      SET1(u_rhs[n],NDD(pnode,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],pnode,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],pnode,NULL,ZB,p);
      }
      IFLAG(pnode) = ++n;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      SET1(u_rhs[n],FDVP(pface,f))
      if (lagr == YES){
         subtract_Dfi_x_ssn(u_rhs[n],pface,NULL,ZB,p);
         subtract_Dfi_x_ssf(u_rhs[n],pface,ZB,p);
      }
      IFLAG(pface) = ++n;
   }
   if (IS_FN(n0)){
      SET1(u_rhs[n],NDD(n0,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],n0,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],n0,NULL,ZB,p);
      }
      IFLAG(n0) = ++n;
   }
   n2 = 2*n;
   k = n - 1;
   k2 = n2 - 1;
   for (i = 0; i < n2; i++)
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;

   if (IS_FN(n0)){
      SET1(b[k],COEFFBP(n0,ZB))
      MFILL_2x2(a,COEFFNNP(n0,ZA),k,k,k2,k2)
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode) - 1;
         in = i + n;
         if (IS_FN(n0))
            MFILL_2x2(a,COEFFLLP(pl,ZA),k,i,k2,in)
         MFILL_2x2(a,COEFFNNP(pnode,ZA),i,i,in,in)
         SET1(b[i],COEFFBLP(pl,ZB))
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn)){
               if (IFLAG(pn)){
                  j = IFLAG(pn)-1;
                  jn = j + n;
                  MFILL_2x2(a,COEFFLLP(pli,ZA),i,j,in,jn)
               }
               else
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
            }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IFLAG(pf=NBFACE(pnfl))){
               j = IFLAG(pf)-1;
               jn = j + n;
               MFILL_2x2(a,COEFFLLP(pnfl,ZA),i,j,in,jn)
            }
            else
               MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf)){
      pface=NBFACE(pnf);
      i = IFLAG(pface)-1;
      in = i + n;
      if (IS_FN(n0))
         MFILL_2x2(a,COEFFLLP(pnf,ZA),k,i,k2,in)
      MFILL_2x2(a,COEFFNNP(pface,ZA),i,i,in,in)
      SET1(b[i],COEFF_NFP(pnf,ZB))
      for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
         pn = NBNODE(pfn);
         if (pn != n0){
            for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
            SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
         }
         if (IS_FN(pn)){
            if (IFLAG(pn)){
               j = IFLAG(pn)-1;
               jn = j + n;
               MFILL_2x2(a,COEFFLLP(pfn,ZA),i,j,in,jn)
            }
            else
               MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
         }
      }
      for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
         if (IFLAG(pf=NBFACE(pfl))){
            j = IFLAG(pf)-1;
            jn = j + n;
            MFILL_2x2(a,COEFFLLP(pfl,ZA),i,j,in,jn)
         }
         else
            MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
   } 

   solve_local_Stokes_Korn(n,n2,a,b,u_rhs,NDS(n0,g),y,&pressure);
/*
   s = fabs(NDS(n0,p) - pressure);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(pressure - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
      ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
      IFLAG(pnode) = 0;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
      FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
      IFLAG(pface) = 0;
   }
}

void half_full_Vanka_step_p2_p1_Korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new, lagr;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a0[N_SGGEM][N_SGGEM], a1[N_SGGEM][N_SGGEM], b[N_SGGEM][2], 
         y[2][N_SGGEM], u_rhs[N_SGGEM][2], q, s, pressure;
   INT i, j, n=0, k;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      SET1(u_rhs[n],NDD(pnode,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],pnode,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],pnode,NULL,ZB,p);
      }
      IFLAG(pnode) = ++n;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      SET1(u_rhs[n],FDVP(pface,f))
      if (lagr == YES){
         subtract_Dfi_x_ssn(u_rhs[n],pface,NULL,ZB,p);
         subtract_Dfi_x_ssf(u_rhs[n],pface,ZB,p);
      }
      IFLAG(pface) = ++n;
   }
   if (IS_FN(n0)){
      SET1(u_rhs[n],NDD(n0,f))
      if (lagr == YES){
         subtract_Dni_x_ssn(u_rhs[n],n0,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[n],n0,NULL,ZB,p);
      }
      IFLAG(n0) = ++n;
   }
   k = n - 1;
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         a0[i][j] = a1[i][j] = 0.;

   if (IS_FN(n0)){
      SET1(b[k],COEFFBP(n0,ZB))
      MFILL_SUBTR_2x2(a0,a1,COEFFNNP(n0,ZA),k,k,u_rhs[k],NDD(n0,u))
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode) - 1;
         if (IS_FN(n0))
            MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pl,ZA),k,i,u_rhs[k],NDD(pnode,u))
         MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pnode,ZA),i,i,u_rhs[i],NDD(pnode,u))
         SET1(b[i],COEFFBLP(pl,ZB))
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn)){
               if (IFLAG(pn)){
                  j = IFLAG(pn)-1;
                  MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pli,ZA),i,j,u_rhs[i],NDD(pn,u))
               }
               else
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
            }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IFLAG(pf=NBFACE(pnfl))){
               j = IFLAG(pf)-1;
               MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pnfl,ZA),i,j,u_rhs[i],FDVP(pf,u))
            }
            else
               MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf)){
      pface=NBFACE(pnf);
      i = IFLAG(pface)-1;
      if (IS_FN(n0))
         MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pnf,ZA),k,i,u_rhs[k],FDVP(pface,u))
      MFILL_SUBTR_2x2(a0,a1,COEFFNNP(pface,ZA),i,i,u_rhs[i],FDVP(pface,u))
      SET1(b[i],COEFF_NFP(pnf,ZB))
      for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
         pn = NBNODE(pfn);
         if (pn != n0){
            for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
            SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
         }
         if (IS_FN(pn)){
            if (IFLAG(pn)){
               j = IFLAG(pn)-1;
               MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pfn,ZA),i,j,u_rhs[i],NDD(pn,u))
            }
            else
               MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
         }
      }
      for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
         if (IFLAG(pf=NBFACE(pfl))){
            j = IFLAG(pf)-1;
            MFILL_SUBTR_2x2(a0,a1,COEFFLLP(pfl,ZA),i,j,u_rhs[i],FDVP(pf,u))
         }
         else
            MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
   } 

   solve_local_Stokes_with_a1_a2(n,a0,a1,b,u_rhs,NDS(n0,g),y,&pressure);
/*
   s = fabs(NDS(n0,p) - pressure);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(pressure - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl)){
      pnode=NBNODE(pl);
      ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
      ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
      IFLAG(pnode) = 0;
   }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl)){
      pface=NBFACE(pnfl);
      FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
      FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
      IFLAG(pface) = 0;
   }
}

#if DATA_S & SPECIAL_NODES_AND_FACES

void A_Korn_contribution_for_p2_p2ln(sn,n,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn;
NODE *ni[3];
FACE *fi[2];
FLOAT a[][N_SGGEM], u_rhs[][2];
INT n, n_ni, n_fi, ZA, f, u;
{
   NODE *pn;
   FACE *pf;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   INT i, in, j, jn, l;
   	
   for (l = 0; l < n_ni; l++){
      i = IFLAG(ni[l]) - 1;
      in = i + n;
      MFILL_2x2(a,COEFFNNP(ni[l],ZA),i,i,in,in)
      for (pli = START(ni[l]); pli; pli = NEXT(pli)){
         pn=NBNODE(pli);
         j  = IFLAG(pn)-1;
         jn = j + n;
         if (IS_LN(pn) && (pn == sn->n0 || pn == sn->n1 || pn == sn->n2))
            MFILL_2x2(a,COEFFLLP(pli,ZA),i,j,in,jn)
         else
            MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
      }
      for (pnfl = NFSTART(ni[l]); pnfl; pnfl = NEXT(pnfl)){
         pf=NBFACE(pnfl);
         j = IFLAG(pf)-1;
         jn = j + n;
         if (IS_LF(pf) && (pf == sn->f1 || pf == sn->f2))
            MFILL_2x2(a,COEFFLLP(pnfl,ZA),i,j,in,jn)
         else
            MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
      }
   }
   for (l = 0; l < n_fi; l++){
      i = IFLAG(fi[l])-1;
      in = i + n;
      MFILL_2x2(a,COEFFNNP(fi[l],ZA),i,i,in,in)
      for (pfn = FNSTART(fi[l]); pfn; pfn = NEXT(pfn)){
         pn = NBNODE(pfn);
         j  = IFLAG(pn)-1;
         jn = j + n;
         if (IS_LN(pn) && (pn == sn->n0 || pn == sn->n1 || pn == sn->n2))
            MFILL_2x2(a,COEFFLLP(pfn,ZA),i,j,in,jn)
         else
            MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
      }
      for (pfl = FSTART(fi[l]); pfl; pfl = NEXT(pfl)){
         pf=NBFACE(pfl);
         j = IFLAG(pf)-1;
         jn = j + n;
         if (IS_LF(pf) && (pf == sn->f1 || pf == sn->f2))
            MFILL_2x2(a,COEFFLLP(pfl,ZA),i,j,in,jn)
         else
            MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
      }
   } 
}

#else

void A_Korn_contribution_for_p2_p2ln(sn,n,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn; NODE *ni[3]; FACE *fi[2]; FLOAT a[][N_SGGEM], u_rhs[][2]; INT n, n_ni, n_fi, ZA, f, u;
{  eprintf("Error: A_Korn_contribution_for_p2_p2ln not available.\n");  }

#endif

void A_Korn_contribution_for_p2_p2lf(n,f0,n1,n2,a,u_rhs,ZA,f,u)
FACE *f0;
NODE *n1, *n2;
FLOAT a[][N_SGGEM], u_rhs[][2];
INT n, ZA, f, u;
{
   NODE *pn;
   FACE *pf;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   INT m1=n+1, m2=n+2;
   	
   MFILL_2x2(a,COEFFNNP(f0,ZA),0,0,n,n)
   for (pfn = FNSTART(f0); pfn; pfn = NEXT(pfn)){
      pn = NBNODE(pfn);
      if (pn == n1)
         MFILL_2x2(a,COEFFLLP(pfn,ZA),0,1,n,m1)
      else if (pn == n2)
         MFILL_2x2(a,COEFFLLP(pfn,ZA),0,2,n,m2)
      else
         MSET9(u_rhs[0],COEFFLLP(pfn,ZA),NDD(pn,u))
   }
   for (pfl = FSTART(f0); pfl; pfl = NEXT(pfl))
      MSET9(u_rhs[0],COEFFLLP(pfl,ZA),FDVP(NBFACE(pfl),u))

   MFILL_2x2(a,COEFFNNP(n1,ZA),1,1,m1,m1)
   for (pli = START(n1); pli; pli = NEXT(pli))
      if ((pn=NBNODE(pli)) == n2)
         MFILL_2x2(a,COEFFLLP(pli,ZA),1,2,m1,m2)
      else
         MSET9(u_rhs[1],COEFFLLP(pli,ZA),NDD(pn,u))
   for (pnfl = NFSTART(n1); pnfl; pnfl = NEXT(pnfl))
      if ((pf=NBFACE(pnfl)) == f0)
         MFILL_2x2(a,COEFFLLP(pnfl,ZA),1,0,m1,n)
      else
         MSET9(u_rhs[1],COEFFLLP(pnfl,ZA),FDVP(pf,u))
   if (n2){
      MFILL_2x2(a,COEFFNNP(n2,ZA),2,2,m2,m2)
      for (pli = START(n2); pli; pli = NEXT(pli))
         if ((pn=NBNODE(pli)) == n1)
            MFILL_2x2(a,COEFFLLP(pli,ZA),2,1,m2,m1)
         else
            MSET9(u_rhs[2],COEFFLLP(pli,ZA),NDD(pn,u))
      for (pnfl = NFSTART(n2); pnfl; pnfl = NEXT(pnfl))
         if ((pf=NBFACE(pnfl)) == f0)
            MFILL_2x2(a,COEFFLLP(pnfl,ZA),2,0,m2,n)
         else
            MSET9(u_rhs[2],COEFFLLP(pnfl,ZA),FDVP(pf,u))
   }
}

#else

void full_Vanka_step_p2_p1_Korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new, lagr; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1_Korn not available.\n");  }

void half_full_Vanka_step_p2_p1_Korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om,lagr)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new, lagr; FLOAT om;
{  eprintf("Error: half_full_Vanka_step_p2_p1_Korn not available.\n");  }

void A_Korn_contribution_for_p2_p2ln(sn,n,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u)
SNODE *sn; NODE *ni[3]; FACE *fi[2]; FLOAT a[][N_SGGEM], u_rhs[][2]; INT n, n_ni, n_fi, ZA, f, u;
{  eprintf("Error: A_Korn_contribution_for_p2_p2ln not available.\n");  }

void A_Korn_contribution_for_p2_p2lf(n,f0,n1,n2,a,u_rhs,ZA,f,u)
FACE *f0; NODE *n1, *n2; FLOAT a[][N_SGGEM], u_rhs[][2]; INT n, ZA, f, u;
{  eprintf("Error: A_Korn_contribution_for_p2_p2lf not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (N_DATA & NODE_IFLAG) && (F_DATA & FACE_IFLAG) && (DATA_S & SPECIAL_NODES_AND_FACES)

void full_Vanka_step_p2_p2ln(sn,ZA,ZB,f,g,u,p,u_new,p_new,om,A_struct)
SNODE *sn;
INT ZA, ZB, f, g, u, p, u_new, p_new, A_struct;
FLOAT om;
{
   NODE *ni[3];
   FACE *fi[2];
   FLOAT a[10][N_SGGEM], b[6][2], y[4][N_SGGEM], u_rhs[6][2], q, s, multiplier;
   INT i, j, n=1, l, n_ni=1, n_fi=0;
   	
   if (IS_LN(sn->n0)){
      ni[0] = sn->n0; 
      SET1(b[0],sn->ann0[ZB])
      if (IS_LN(sn->n1)){
         SET1(b[n],sn->ann1[ZB])
         n++;
         ni[n_ni++] = sn->n1;
      }
      if (IS_LN(sn->n2) && sn->pel2){
         SET1(b[n],sn->ann2[ZB])
         n++;
         ni[n_ni++] = sn->n2;
      }
      if (IS_LF(sn->f1)){
         SET1(b[n],sn->anf1[ZB])
         n++;
         fi[n_fi++] = sn->f1;
      }
      if (IS_LF(sn->f2) && sn->pel2){
         SET1(b[n],sn->anf2[ZB])
         n++;
         fi[n_fi++] = sn->f2;
      }
      n = 0;
      for (l = 0; l < n_ni; l++){
         SET1(u_rhs[n],NDD(ni[l],f))
         subtract_Dni_x_ssn(u_rhs[n],ni[l],ni[0],ZB,p,l);
         subtract_Dni_x_ssf(u_rhs[n],ni[l],NULL,ZB,p);
         subtract_Bni_x_sn(u_rhs[n],ni[l],NULL,ZB,p);
         IFLAG(ni[l]) = ++n;
      }
      for (l = 0; l < n_fi; l++){
         SET1(u_rhs[n],FDVP(fi[l],f))
         subtract_Dfi_x_ssn(u_rhs[n],fi[l],ni[0],ZB,p);
         subtract_Dfi_x_ssf(u_rhs[n],fi[l],ZB,p);
         subtract_Bfi_x_sn(u_rhs[n],fi[l],NULL,ZB,p);
         IFLAG(fi[l]) = ++n;
      }
      for (i = 0; i < 10; i++)
         for (j = 0; j < 10; j++)
            a[i][j] = 0.;

      if (A_struct & Q_FULL){
         A_Korn_contribution_for_p2_p2ln(sn,n,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u);
         solve_local_Stokes_Korn(n,2*n,a,b,u_rhs,SNDS(sn,g),y,&multiplier);
      }
      else if (A_struct & Q_EBDIAG){
         A_contribution_for_p2_p2ln(sn,ni,fi,n_ni,n_fi,a,u_rhs,ZA,f,u);
         solve_local_Stokes(n,a,b,u_rhs,SNDS(sn,g),y,&multiplier);
      }
      else
         eprintf("Error: wrong A_struct in full_Vanka_step_p2_p2ln.\n");
/*
   s = fabs(multiplier - SNDS(sn,p));
   for (l = 0; l < n_ni; l++)
      if ((q=MAX(fabs(y[0][IFLAG(ni[l])-1] - ND(ni[l],u,0)),
                 fabs(y[1][IFLAG(ni[l])-1] - ND(ni[l],u,1)) ))>s) s = q;
   for (l = 0; l < n_fi; l++)
      if ((q=MAX(fabs(y[0][IFLAG(fi[l])-1]-FDV(fi[l],u,0)),
                 fabs(y[1][IFLAG(fi[l])-1]-FDV(fi[l],u,1)) ))>s) s = q;
   printf("Max. diff. = %e\n",s);
*/

      SNDS(sn,p_new) = SNDS(sn,p) + om*(multiplier - SNDS(sn,p));
      for (l = 0; l < n_ni; l++){
         ND(ni[l],u_new,0) = 
                      ND(ni[l],u,0) + om*(y[0][IFLAG(ni[l])-1] - ND(ni[l],u,0));
         ND(ni[l],u_new,1) = 
                      ND(ni[l],u,1) + om*(y[1][IFLAG(ni[l])-1] - ND(ni[l],u,1));
      }
      for (l = 0; l < n_fi; l++){
         FDV(fi[l],u_new,0) = 
                      FDV(fi[l],u,0) + om*(y[0][IFLAG(fi[l])-1]-FDV(fi[l],u,0));
         FDV(fi[l],u_new,1) = 
                      FDV(fi[l],u,1) + om*(y[1][IFLAG(fi[l])-1]-FDV(fi[l],u,1));
      }
      IFLAG(sn->n0) = IFLAG(sn->n1) = IFLAG(sn->n2) 
                    = IFLAG(sn->f1) = IFLAG(sn->f2) = 0;
   }
}

void full_Vanka_step_p2_p2lf(sf,ZA,ZB,f,g,u,p,u_new,p_new,om,A_struct)
SFACE *sf;
INT ZA, ZB, f, g, u, p, u_new, p_new, A_struct;
FLOAT om;
{
   NODE *n1, *n2;
   FACE *f0;
   FLOAT a[6][N_SGGEM], b[4][2], y[4][N_SGGEM], u_rhs[4][2], q, s, multiplier;
   INT i, j, n=2;
   	
   if (IS_LF(sf->f)){
      f0 = sf->f;
      SET1(b[0],sf->aff[ZB])
      if (IS_LN(sf->n1)){
         n1 = sf->n1;
         SET1(b[1],sf->afn1[ZB])
         if (IS_LN(sf->n2)){
            n2 = sf->n2;
            SET1(b[2],sf->afn2[ZB])
         }
         else
            n2 = NULL;
      }
      else{
         n1 = sf->n2;
         SET1(b[1],sf->afn2[ZB])
         n2 = NULL;
      }
      SET1(u_rhs[0],FDVP(f0,f))
      SET1(u_rhs[1],NDD(n1,f))
      subtract_Dni_x_ssn(u_rhs[1],n1,NULL,ZB,p,1);
      subtract_Dni_x_ssf(u_rhs[1],n1,f0,ZB,p);
      subtract_Bni_x_sn(u_rhs[1],n1,NULL,ZB,p);
      if (n2){
         n = 3;
         SET1(u_rhs[2],NDD(n2,f))
         subtract_Dni_x_ssn(u_rhs[2],n2,NULL,ZB,p,1);
         subtract_Dni_x_ssf(u_rhs[2],n2,f0,ZB,p);
         subtract_Bni_x_sn(u_rhs[2],n2,NULL,ZB,p);
      }

      if (n1->s_node->f1 == f0)
         SET9(u_rhs[0],n1->s_node->anf1[ZB],SNDS(n1->s_node,p))
      else
         SET9(u_rhs[0],n1->s_node->anf2[ZB],SNDS(n1->s_node,p))
      if (n2 && n2->s_node->f1 == f0)
         SET9(u_rhs[0],n2->s_node->anf1[ZB],SNDS(n2->s_node,p))
      else if (n2 && n2->s_node->f2 == f0)
         SET9(u_rhs[0],n2->s_node->anf2[ZB],SNDS(n2->s_node,p))
      subtract_Bfi_x_sn(u_rhs[0],f0,NULL,ZB,p);

      for (i = 0; i < 6; i++)
         for (j = 0; j < 6; j++)
            a[i][j] = 0.;

      if (A_struct & Q_FULL){
         A_Korn_contribution_for_p2_p2lf(n,f0,n1,n2,a,u_rhs,ZA,f,u);
         solve_local_Stokes_Korn(n,2*n,a,b,u_rhs,SFDS(sf,g),y,&multiplier);
      }
      else if (A_struct & Q_EBDIAG){
         A_contribution_for_p2_p2lf(f0,n1,n2,a,u_rhs,ZA,f,u);
         solve_local_Stokes(n,a,b,u_rhs,SFDS(sf,g),y,&multiplier);
      }
      else
         eprintf("Error: wrong A_struct in full_Vanka_step_p2_p2lf.\n");
/*
   s = fabs(multiplier - SFDS(sf,p));
   if ((q=MAX(fabs(y[0][0]-FDV(f0,u,0)),fabs(y[1][0]-FDV(f0,u,1)))) > s) s = q;
   if ((q=MAX(fabs(y[0][1]-ND(n1,u,0)),fabs(y[1][1]-ND(n1,u,1)))) > s) s = q;
   if (n2 && (q=MAX(fabs(y[0][2]-ND(n2,u,0)),fabs(y[1][2]-ND(n2,u,1)))) > s) 
                                                                       s = q;
   printf("Max. diff. = %e\n",s);
*/
      SFDS(sf,p_new) = SFDS(sf,p) + om*(multiplier - SFDS(sf,p));
      FDV(f0,u_new,0) = FDV(f0,u,0) + om*(y[0][0]-FDV(f0,u,0));
      FDV(f0,u_new,1) = FDV(f0,u,1) + om*(y[1][0]-FDV(f0,u,1));
      ND(n1,u_new,0) = ND(n1,u,0) + om*(y[0][1] - ND(n1,u,0)); 
      ND(n1,u_new,1) = ND(n1,u,1) + om*(y[1][1] - ND(n1,u,1));
      if (n2){
         ND(n2,u_new,0) = ND(n2,u,0) + om*(y[0][2] - ND(n2,u,0)); 
         ND(n2,u_new,1) = ND(n2,u,1) + om*(y[1][2] - ND(n2,u,1)); 
      }
   }
}

#else

void full_Vanka_step_p2_p2ln(sn,ZA,ZB,f,g,u,p,u_new,p_new,om,A_struct)
SNODE *sn; INT ZA, ZB, f, g, u, p, u_new, p_new, A_struct; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p2ln not available.\n");  }

void full_Vanka_step_p2_p2lf(sf,ZA,ZB,f,g,u,p,u_new,p_new,om,A_struct)
SFACE *sf; INT ZA, ZB, f, g, u, p, u_new, p_new, A_struct; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p2lf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExFxDN_MATR) && (E_DATA & ExFxDF_MATR) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p2_p1disc(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[8][N_SGGEM], b0[3][8], b1[3][8], c[8][N_SGGEM], y[8][N_SGGEM],
         u_rhs[8][2], r[3], x[3], q, s;
   INT i, j, n=0, in0, in1, in2, ifa0, ifa1, ifa2;

   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n)
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   for (i = 0; i < n; i++){
      b0[0][i] = b0[1][i] = b0[2][i] = b1[0][i] = b1[1][i] = b1[2][i] = 0.;
      for (j = 0; j < n; j++)
         a[i][j] = 0.;
   }
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      FILL_P1DISC_B_MATR(pelem,i,0,COEFF_BFDN)
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      FILL_P1DISC_B_MATR(pelem,i,1,COEFF_BFDN)
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      FILL_P1DISC_B_MATR(pelem,i,2,COEFF_BFDN)
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      FILL_P1DISC_B_MATR(pelem,i,0,COEFF_BFDF)
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      FILL_P1DISC_B_MATR(pelem,i,1,COEFF_BFDF)
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      FILL_P1DISC_B_MATR(pelem,i,2,COEFF_BFDF)
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0)
               SUBTR_P1DISC_B_MATR(u_rhs[in0],pel,i,COEFF_BFDNP)
            else if (pel->n[i] == n1)
               SUBTR_P1DISC_B_MATR(u_rhs[in1],pel,i,COEFF_BFDNP)
            else if (pel->n[i] == n2)
               SUBTR_P1DISC_B_MATR(u_rhs[in2],pel,i,COEFF_BFDNP)
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0)
               SUBTR_P1DISC_B_MATR(u_rhs[ifa0],pel,i,COEFF_BFDFP)
            else if (pel->f[i] == fa1)
               SUBTR_P1DISC_B_MATR(u_rhs[ifa1],pel,i,COEFF_BFDFP)
            else if (pel->f[i] == fa2)
               SUBTR_P1DISC_B_MATR(u_rhs[ifa2],pel,i,COEFF_BFDFP)
         }
      }
   }

   for (i = 0; i < n; i++){
      c[0][i] = b0[0][i];
      c[1][i] = b0[1][i];
      c[2][i] = b0[2][i];
      c[3][i] = b1[0][i];
      c[4][i] = b1[1][i];
      c[5][i] = b1[2][i];
      c[6][i] = u_rhs[i][0];
      c[7][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,8);
   a[0][0] = a[0][1] = a[0][2] =
   a[1][0] = a[1][1] = a[1][2] =
   a[2][0] = a[2][1] = a[2][2] = 0.;
   r[0] = -EDSN(pelem,g,0);
   r[1] = -EDSN(pelem,g,1);
   r[2] = -EDSN(pelem,g,2);
   for (i = 0; i < n; i++){
      a[0][0] += b0[0][i]*y[0][i] + b1[0][i]*y[3][i];
      a[0][1] += b0[0][i]*y[1][i] + b1[0][i]*y[4][i];
      a[0][2] += b0[0][i]*y[2][i] + b1[0][i]*y[5][i];
      a[1][0] += b0[1][i]*y[0][i] + b1[1][i]*y[3][i];
      a[1][1] += b0[1][i]*y[1][i] + b1[1][i]*y[4][i];
      a[1][2] += b0[1][i]*y[2][i] + b1[1][i]*y[5][i];
      a[2][0] += b0[2][i]*y[0][i] + b1[2][i]*y[3][i];
      a[2][1] += b0[2][i]*y[1][i] + b1[2][i]*y[4][i];
      a[2][2] += b0[2][i]*y[2][i] + b1[2][i]*y[5][i];
      r[0] += b0[0][i]*y[6][i] + b1[0][i]*y[7][i];
      r[1] += b0[1][i]*y[6][i] + b1[1][i]*y[7][i];
      r[2] += b0[2][i]*y[6][i] + b1[2][i]*y[7][i];
   }
   ggem(a,r,x,3);
   for (i = 0; i < n; i++){
      y[0][i] = y[6][i] - y[0][i]*x[0] - y[1][i]*x[1] - y[2][i]*x[2];
      y[1][i] = y[7][i] - y[3][i]*x[0] - y[4][i]*x[1] - y[5][i]*x[2];
   }

/*
   s = 0.;
   if((q=fabs(EDSN(pelem,p,0)-x[0])) > s) s = q;
   if((q=fabs(EDSN(pelem,p,1)-x[1])) > s) s = q;
   if((q=fabs(EDSN(pelem,p,2)-x[2])) > s) s = q;
   COMPUTE_MAX_ND_DIFF
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   EDSN(pelem,p_new,0) = EDSN(pelem,p,0) + om*(x[0] - EDSN(pelem,p,0));
   EDSN(pelem,p_new,1) = EDSN(pelem,p,1) + om*(x[1] - EDSN(pelem,p,1));
   EDSN(pelem,p_new,2) = EDSN(pelem,p,2) + om*(x[2] - EDSN(pelem,p,2));
   UPDATE_ND_VALUES
   UPDATE_FDV_VALUES
}

#else

void full_Vanka_step_p2_p1disc(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1disc not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & E_E_NEIGHBOURS) && (E_DATA & ExFxDN_MATR) && (E_DATA & ExFxDF_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (E_DATA & VECTOR_ELEMENT_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void full_Vanka_step_p2b_p1disc(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   LINK *pli;
   NFLINK *pnfl;
   FNLINK *pfnl;
   FLINK *pfl;
   ELINK *peli;
   FLOAT a[8][N_SGGEM], b0[3][8], b1[3][8], c[8][N_SGGEM], y[8][N_SGGEM], 
         u_rhs[8][2], r[3], x[3], q, s;
   INT i, j, n=1, in0, in1, in2, ifa0, ifa1, ifa2;
   	
   NODES_OF_ELEMENT(n0,n1,n2,pelem);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
   LOCAL_NODE_INDICES(n0,n1,n2,in0,in1,in2,n)
   LOCAL_FACE_INDICES(fa0,fa1,fa2,ifa0,ifa1,ifa2,n)
   for (i = 0; i < n; i++){
      b0[0][i] = b0[1][i] = b0[2][i] = b1[0][i] = b1[1][i] = b1[2][i] = 0.;
      for (j = 0; j < n; j++)
         a[i][j] = 0.;
   }
   a[0][0] = COEFF_EE(pelem,ZA);
   b0[0][0] = COEFF_BDF(pelem,ZB,0,0);
   b0[1][0] = COEFF_BDF(pelem,ZB,1,0);
   b0[2][0] = COEFF_BDF(pelem,ZB,2,0);
   b1[0][0] = COEFF_BDF(pelem,ZB,0,1);
   b1[1][0] = COEFF_BDF(pelem,ZB,1,1);
   b1[2][0] = COEFF_BDF(pelem,ZB,2,1);
   SET1(u_rhs[0],EDVP(pelem,f))
   if (IS_FN(n0)){
      i = in0;
      SET1(u_rhs[i],NDD(n0,f))
      FILL_P1DISC_B_MATR(pelem,i,0,COEFF_BFDN)
      a[0][i] = COEFF_EN(pelem,ZA,0);
      a[i][0] = COEFF_NE(pelem,ZA,0);
      FILL_AND_SUBTR_NN_MATR(n0,n1,n2,in1,in2)
      FILL_AND_SUBTR_NF_MATR(n0)
   }
   if (IS_FN(n1)){
      i = in1;
      SET1(u_rhs[i],NDD(n1,f))
      FILL_P1DISC_B_MATR(pelem,i,1,COEFF_BFDN)
      a[0][i] = COEFF_EN(pelem,ZA,1);
      a[i][0] = COEFF_NE(pelem,ZA,1);
      FILL_AND_SUBTR_NN_MATR(n1,n0,n2,in0,in2)
      FILL_AND_SUBTR_NF_MATR(n1)
   }
   if (IS_FN(n2)){
      i = in2;
      SET1(u_rhs[i],NDD(n2,f))
      FILL_P1DISC_B_MATR(pelem,i,2,COEFF_BFDN)
      a[0][i] = COEFF_EN(pelem,ZA,2);
      a[i][0] = COEFF_NE(pelem,ZA,2);
      FILL_AND_SUBTR_NN_MATR(n2,n0,n1,in0,in1)
      FILL_AND_SUBTR_NF_MATR(n2)
   }
   if (IS_FF(fa0)){
      i = ifa0;
      SET1(u_rhs[i],FDVP(fa0,f))
      FILL_P1DISC_B_MATR(pelem,i,0,COEFF_BFDF)
      a[0][i] = COEFF_EF(pelem,ZA,0);
      a[i][0] = COEFF_FE(pelem,ZA,0);
      FILL_AND_SUBTR_FN_MATR(fa0)
      FILL_AND_SUBTR_FF_MATR(fa0,fa1,fa2,ifa1,ifa2)
   }
   if (IS_FF(fa1)){
      i = ifa1;
      SET1(u_rhs[i],FDVP(fa1,f))
      FILL_P1DISC_B_MATR(pelem,i,1,COEFF_BFDF)
      a[0][i] = COEFF_EF(pelem,ZA,1);
      a[i][0] = COEFF_FE(pelem,ZA,1);
      FILL_AND_SUBTR_FN_MATR(fa1)
      FILL_AND_SUBTR_FF_MATR(fa1,fa0,fa2,ifa0,ifa2)
   }
   if (IS_FF(fa2)){
      i = ifa2;
      SET1(u_rhs[i],FDVP(fa2,f))
      FILL_P1DISC_B_MATR(pelem,i,2,COEFF_BFDF)
      a[0][i] = COEFF_EF(pelem,ZA,2);
      a[i][0] = COEFF_FE(pelem,ZA,2);
      FILL_AND_SUBTR_FN_MATR(fa2)
      FILL_AND_SUBTR_FF_MATR(fa2,fa0,fa1,ifa0,ifa1)
   }
   for (peli = pelem->estart; peli != NULL; peli = peli->next){
      pel = peli->nbel;
      for (i = 0; i < NVERT; i++){
         if (IS_FN(pel->n[i])){
            if (pel->n[i] == n0){
               SUBTR_P1DISC_B_MATR(u_rhs[in0],pel,i,COEFF_BFDNP)
               SET9(u_rhs[in0],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n1){
               SUBTR_P1DISC_B_MATR(u_rhs[in1],pel,i,COEFF_BFDNP)
               SET9(u_rhs[in1],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
            else if (pel->n[i] == n2){
               SUBTR_P1DISC_B_MATR(u_rhs[in2],pel,i,COEFF_BFDNP)
               SET9(u_rhs[in2],EDVP(pel,u),COEFF_NE(pel,ZA,i))
            }
         }
         if (IS_FF(pel->f[i])){
            if (pel->f[i] == fa0){
               SUBTR_P1DISC_B_MATR(u_rhs[ifa0],pel,i,COEFF_BFDFP)
               SET9(u_rhs[ifa0],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa1){
               SUBTR_P1DISC_B_MATR(u_rhs[ifa1],pel,i,COEFF_BFDFP)
               SET9(u_rhs[ifa1],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
            else if (pel->f[i] == fa2){
               SUBTR_P1DISC_B_MATR(u_rhs[ifa2],pel,i,COEFF_BFDFP)
               SET9(u_rhs[ifa2],EDVP(pel,u),COEFF_FE(pel,ZA,i))
            }
         }
      }
   }

   for (i = 0; i < n; i++){
      c[0][i] = b0[0][i];
      c[1][i] = b0[1][i];
      c[2][i] = b0[2][i];
      c[3][i] = b1[0][i];
      c[4][i] = b1[1][i];
      c[5][i] = b1[2][i];
      c[6][i] = u_rhs[i][0];
      c[7][i] = u_rhs[i][1];
   }
   sggem(a,c,y,n,8);
   a[0][0] = a[0][1] = a[0][2] =
   a[1][0] = a[1][1] = a[1][2] =
   a[2][0] = a[2][1] = a[2][2] = 0.;
   r[0] = -EDSN(pelem,g,0);
   r[1] = -EDSN(pelem,g,1);
   r[2] = -EDSN(pelem,g,2);
   for (i = 0; i < n; i++){
      a[0][0] += b0[0][i]*y[0][i] + b1[0][i]*y[3][i]; 
      a[0][1] += b0[0][i]*y[1][i] + b1[0][i]*y[4][i];
      a[0][2] += b0[0][i]*y[2][i] + b1[0][i]*y[5][i];
      a[1][0] += b0[1][i]*y[0][i] + b1[1][i]*y[3][i];
      a[1][1] += b0[1][i]*y[1][i] + b1[1][i]*y[4][i];
      a[1][2] += b0[1][i]*y[2][i] + b1[1][i]*y[5][i];
      a[2][0] += b0[2][i]*y[0][i] + b1[2][i]*y[3][i];
      a[2][1] += b0[2][i]*y[1][i] + b1[2][i]*y[4][i];
      a[2][2] += b0[2][i]*y[2][i] + b1[2][i]*y[5][i];
      r[0] += b0[0][i]*y[6][i] + b1[0][i]*y[7][i];
      r[1] += b0[1][i]*y[6][i] + b1[1][i]*y[7][i];
      r[2] += b0[2][i]*y[6][i] + b1[2][i]*y[7][i];
   }
   ggem(a,r,x,3);
   for (i = 0; i < n; i++){
      y[0][i] = y[6][i] - y[0][i]*x[0] - y[1][i]*x[1] - y[2][i]*x[2];
      y[1][i] = y[7][i] - y[3][i]*x[0] - y[4][i]*x[1] - y[5][i]*x[2];
   }

/*
   s = MAX(fabs(EDV(pelem,u,0) - y[0][0]),fabs(EDV(pelem,u,1) - y[1][0]));
   if((q=fabs(EDSN(pelem,p,0)-x[0])) > s) s = q;
   if((q=fabs(EDSN(pelem,p,1)-x[1])) > s) s = q;
   if((q=fabs(EDSN(pelem,p,2)-x[2])) > s) s = q;
   COMPUTE_MAX_ND_DIFF
   COMPUTE_MAX_FDV_DIFF
   printf("Max. diff. = %e\n",s);
*/

   EDSN(pelem,p_new,0) = EDSN(pelem,p,0) + om*(x[0] - EDSN(pelem,p,0));
   EDSN(pelem,p_new,1) = EDSN(pelem,p,1) + om*(x[1] - EDSN(pelem,p,1));
   EDSN(pelem,p_new,2) = EDSN(pelem,p,2) + om*(x[2] - EDSN(pelem,p,2));
   EDV(pelem,u_new,0) = EDV(pelem,u,0) + om*(y[0][0] - EDV(pelem,u,0));
   EDV(pelem,u_new,1) = EDV(pelem,u,1) + om*(y[1][0] - EDV(pelem,u,1));
   UPDATE_ND_VALUES
   UPDATE_FDV_VALUES
}

#else

void full_Vanka_step_p2b_p1disc(pelem,ZA,ZB,f,g,u,p,u_new,p_new,om)
ELEMENT *pelem; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2b_p1disc not available.\n");  }

#endif

#if (E_DATA & E_E_NEIGHBOURS) && (DIM == 3)

void p_to_rhs_diag(pel,nrhs,frhs,xe)
ELEMENT *pel;
FLOAT nrhs[4][DIM], frhs[4];
INT xe;
{
   ELEMENT *pelem;
   ELINK *peli;
   INT i;
   	
   for (peli = pel->estart; peli != NULL; peli = peli->next){
      pelem = peli->nbel;
      for (i = 0; i < 4; i++){
         if (pelem->n[i] == pel->n[0])
            SET9(nrhs[0],COEFF_BNP(pelem,0,i),ED(pelem,xe))
         else if (pelem->n[i] == pel->n[1])
            SET9(nrhs[1],COEFF_BNP(pelem,0,i),ED(pelem,xe))
         else if (pelem->n[i] == pel->n[2])
            SET9(nrhs[2],COEFF_BNP(pelem,0,i),ED(pelem,xe))
         else if (pelem->n[i] == pel->n[3])
            SET9(nrhs[3],COEFF_BNP(pelem,0,i),ED(pelem,xe))
         if (pelem->f[i] == pel->f[0])
            frhs[0] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
         else if (pelem->f[i] == pel->f[1])
            frhs[1] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
         else if (pelem->f[i] == pel->f[2])
            frhs[2] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
         else if (pelem->f[i] == pel->f[3])
            frhs[3] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
      }
   }
}

void p_to_rhs_stab(pel,rhs,xe)
ELEMENT *pel;
FLOAT rhs[17];
INT xe;
{
   ELEMENT *pelem;
   ELINK *peli;
   INT i;
   	
   for (peli = pel->estart; peli != NULL; peli = peli->next){
      pelem = peli->nbel;
      for (i = 0; i < 4; i++){
         if (IS_FN(pelem->n[i])){
            if (pelem->n[i] == pel->n[0]){
               rhs[0] -= COEFF_BN(pelem,0,i,0)*ED(pelem,xe);
               rhs[4] -= COEFF_BN(pelem,0,i,1)*ED(pelem,xe);
               rhs[8] -= COEFF_BN(pelem,0,i,2)*ED(pelem,xe);
            }
            else if (pelem->n[i] == pel->n[1]){
               rhs[1] -= COEFF_BN(pelem,0,i,0)*ED(pelem,xe);
               rhs[5] -= COEFF_BN(pelem,0,i,1)*ED(pelem,xe);
               rhs[9] -= COEFF_BN(pelem,0,i,2)*ED(pelem,xe);
            }
            else if (pelem->n[i] == pel->n[2]){
               rhs[2]  -= COEFF_BN(pelem,0,i,0)*ED(pelem,xe);
               rhs[6]  -= COEFF_BN(pelem,0,i,1)*ED(pelem,xe);
               rhs[10] -= COEFF_BN(pelem,0,i,2)*ED(pelem,xe);
            }
            else if (pelem->n[i] == pel->n[3]){
               rhs[3]  -= COEFF_BN(pelem,0,i,0)*ED(pelem,xe);
               rhs[7]  -= COEFF_BN(pelem,0,i,1)*ED(pelem,xe);
               rhs[11] -= COEFF_BN(pelem,0,i,2)*ED(pelem,xe);
            }
         }
         if (IS_FF(pelem->f[i])){
            if (pelem->f[i] == pel->f[0])
               rhs[12] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
            else if (pelem->f[i] == pel->f[1])
               rhs[13] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
            else if (pelem->f[i] == pel->f[2])
               rhs[14] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
            else if (pelem->f[i] == pel->f[3])
               rhs[15] -= COEFF_BF(pelem,0,i)*ED(pelem,xe);
         }
      }
   }
}

#else

void p_to_rhs_diag(pel,nrhs,frhs,xe)
ELEMENT *pel; FLOAT nrhs[4][3], frhs[4]; INT xe;
{  eprintf("Error: p_to_rhs_diag not available.\n");  }

void p_to_rhs_stab(pel,rhs,xe)
ELEMENT *pel; FLOAT rhs[17]; INT xe;
{  eprintf("Error: p_to_rhs_stab not available.\n");  }

#endif

/******************************************************************************/
/*                                                                            */
/*                        smoothing by Martin Sodomka                         */
/*                                                                            */
/******************************************************************************/

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR) && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR)

void full_Vanka_step_p2_p1_korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[N_SGGEM][N_SGGEM], b0[N_SGGEM], b1[N_SGGEM],
	     c[2][N_SGGEM], y[2][N_SGGEM], u_rhs[N_SGGEM][2],
         q, s, pressure;
   INT i, in, j, jn, n=0, n2, k, k2;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         SET1(u_rhs[n],NDD(pnode,f))
         IFLAG(pnode) = ++n;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         SET1(u_rhs[n],FDVP(pface,f))
         IFLAG(pface) = ++n;
      }
   if (IS_FN(n0)){
      SET1(u_rhs[n],NDD(n0,f))
      IFLAG(n0) = ++n;
   }
   n2 = 2*n;
   k = n - 1;
   k2 = n2 - 1;
   for (i = 0; i < n2; i++)
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;

   if (IS_FN(n0)){
      b0[k] = COEFFB(n0,ZB,0);
      b1[k] = COEFFB(n0,ZB,1);
      a[k][k] = COEFFNN(n0,ZA,0,0);
	  a[k][k2] = COEFFNN(n0,ZA,0,1);
	  a[k2][k] = COEFFNN(n0,ZA,1,0);
	  a[k2][k2] = COEFFNN(n0,ZA,1,1);
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode) - 1;
		 in = i + n;
         if (IS_FN(n0)) {
            a[k][i] = COEFFLL(pl,ZA,0,0);
			a[k][in] = COEFFLL(pl,ZA,0,1);
			a[k2][i] = COEFFLL(pl,ZA,1,0);
			a[k2][in] = COEFFLL(pl,ZA,1,1);
		 }
         a[i][i] = COEFFNN(pnode,ZA,0,0);
		 a[i][in] = COEFFNN(pnode,ZA,0,1);
		 a[in][i] = COEFFNN(pnode,ZA,1,0);
		 a[in][in] = COEFFNN(pnode,ZA,1,1);
         b0[i] = COEFFBL(pl,ZB,0);
         b1[i] = COEFFBL(pl,ZB,1);
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pli,ZA,0,0);
				  a[i][jn] = COEFFLL(pli,ZA,0,1);
				  a[in][j] = COEFFLL(pli,ZA,1,0);
				  a[in][jn] = COEFFLL(pli,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IS_FF(pf=NBFACE(pnfl))){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pnfl,ZA,0,0);
				  a[i][jn] = COEFFLL(pnfl,ZA,0,1);
				  a[in][j] = COEFFLL(pnfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pnfl,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
            } 
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf))
      if (IS_FF(pface=NBFACE(pnf))){
         i = IFLAG(pface)-1;
		 in = i + n;
         if (IS_FN(n0)) {
            a[k][i] = COEFFLL(pnf,ZA,0,0);
			a[k][in] = COEFFLL(pnf,ZA,0,1);
			a[k2][i] = COEFFLL(pnf,ZA,1,0);
			a[k2][in] = COEFFLL(pnf,ZA,1,1);
		 }
         a[i][i] = COEFFNN(pface,ZA,0,0);
		 a[i][in] = COEFFNN(pface,ZA,0,1);
		 a[in][i] = COEFFNN(pface,ZA,1,0);
		 a[in][in] = COEFFNN(pface,ZA,1,1);
         b0[i] = COEFF_NF(pnf,ZB,0);
         b1[i] = COEFF_NF(pnf,ZB,1);
         for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
            pn = NBNODE(pfn);
            if (pn != n0){
               for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
               SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pfn,ZA,0,0);
				  a[i][jn] = COEFFLL(pfn,ZA,0,1);
				  a[in][j] = COEFFLL(pfn,ZA,1,0);
				  a[in][jn] = COEFFLL(pfn,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
         }
         for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
            if (IS_FF(pf=NBFACE(pfl))){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pfl,ZA,0,0);
				  a[i][jn] = COEFFLL(pfl,ZA,0,1);
				  a[in][j] = COEFFLL(pfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pfl,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
            } 
      } 

   for (i = 0; i < n; i++){
	  in = i + n;
      c[0][i] = b0[i];   	
      c[0][in] = b1[i];   	
      c[1][i] = u_rhs[i][0];
      c[1][in] = u_rhs[i][1];
   }
   sggem(a,c,y,n2,2);
   s = q = 0.;
   for (i = 0; i < n; i++){
	  in = i + n;
      s += b0[i]*y[0][i] + b1[i]*y[0][in];
      q += b0[i]*y[1][i] + b1[i]*y[1][in];
   }
   pressure = (q - NDS(n0,g))/s;
   for (i = 0; i < n; i++){
	  in = i + n;
      y[0][i] = y[1][i] - y[0][i]*pressure;
      y[1][i] = y[1][in] - y[0][in]*pressure;
   }

/*
   s = fabs(NDS(n0,p) - pressure);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(pressure - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
         ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
         IFLAG(pnode) = 0;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
         FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
         IFLAG(pface) = 0;
      }
}

#else /* #if !( (N_DATA & DxD_NODE_MATR) && (N_DATA & DxD_NODE_FACE_MATR)
                && (F_DATA & DxD_FACE_MATR) && (F_DATA & DxD_FACE_NODE_MATR) ) */

void full_Vanka_step_p2_p1_korn(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1_korn not available.\n");  }

#endif

#if (N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA)

void full_Vanka_step_p2_p1_korn_BD(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[N_SGGEM][N_SGGEM], b0[N_SGGEM], b1[N_SGGEM],
	     c[2][N_SGGEM], y[2][N_SGGEM], u_rhs[N_SGGEM][2],
         q, s, pressure;
   INT i, in, j, jn, n=0, n2, k, k2;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         SET1(u_rhs[n],NDD(pnode,f))
         IFLAG(pnode) = ++n;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         SET1(u_rhs[n],FDVP(pface,f))
         IFLAG(pface) = ++n;
      }
   if (IS_FN(n0)){
      SET1(u_rhs[n],NDD(n0,f))
      IFLAG(n0) = ++n;
   }
   n2 = 2*n;
   k = n - 1;
   k2 = n2 - 1;
   for (i = 0; i < n2; i++) {
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;
   }

   if (IS_FN(n0)){
      b0[k] = COEFFB(n0,ZB,0);
      b1[k] = COEFFB(n0,ZB,1);
      a[k][k] = COEFFNN(n0,ZA,0,0);
	  a[k][k2] = COEFFNN(n0,ZA,0,1);
	  a[k2][k] = COEFFNN(n0,ZA,1,0);
	  a[k2][k2] = COEFFNN(n0,ZA,1,1);
	  if (IS_ZNN(n0))
		  SET9(u_rhs[k],COEFFDP(n0,ZB),NDSBD(n0,p))
   }
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
      if (IS_FN(n0)){
         for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
         SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
		 if (IS_ZNN(n0) && IS_ZNN(pnode))
			 SET9(u_rhs[k],COEFFDP(pli,ZB),NDSBD(pnode,p))
      }
      if (IS_FN(pnode)){
         i = IFLAG(pnode) - 1;
		 in = i + n;
         if (IS_FN(n0)) {
            a[k][i] = COEFFLL(pl,ZA,0,0);
			a[k][in] = COEFFLL(pl,ZA,0,1);
			a[k2][i] = COEFFLL(pl,ZA,1,0);
			a[k2][in] = COEFFLL(pl,ZA,1,1);
		 }
         a[i][i] = COEFFNN(pnode,ZA,0,0);
		 a[i][in] = COEFFNN(pnode,ZA,0,1);
		 a[in][i] = COEFFNN(pnode,ZA,1,0);
		 a[in][in] = COEFFNN(pnode,ZA,1,1);
         b0[i] = COEFFBL(pl,ZB,0);
         b1[i] = COEFFBL(pl,ZB,1);
         SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
		 if (IS_ZNN(pnode)) {
			SET9(u_rhs[i],COEFFDP(pnode,ZB),NDSBD(pnode,p))
			if (IS_ZNN(n0))
				SET9(u_rhs[i],COEFFDP(pl,ZB),NDSBD(n0,p))
		 }
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
               SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
			   if (IS_ZNN(pnode) && IS_ZNN(pn))
				   SET9(u_rhs[i],COEFFDP(pli2,ZB),NDSBD(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pli,ZA,0,0);
				  a[i][jn] = COEFFLL(pli,ZA,0,1);
				  a[in][j] = COEFFLL(pli,ZA,1,0);
				  a[in][jn] = COEFFLL(pli,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
			   }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl)){
			pf=NBFACE(pnfl);
			if (IS_ZNN(pnode) && IS_ZNF(pf)){
				for (pfn = FNSTART(pf); pnode != NBNODE(pfn); pfn = NEXT(pfn));
				SET9(u_rhs[i],COEFFDP(pfn,ZB),FDBD(pf,p))
			}
            if (IS_FF(pf)){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pnfl,ZA,0,0);
				  a[i][jn] = COEFFLL(pnfl,ZA,0,1);
				  a[in][j] = COEFFLL(pnfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pnfl,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
			   }
            }
		 }
      } 
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf)){
	  pface=NBFACE(pnf);
	  if (IS_ZNN(n0) && IS_ZNF(pface)){
         for (pfn = FNSTART(pface); n0 != NBNODE(pfn); pfn = NEXT(pfn));
         SET9(u_rhs[k],COEFFDP(pfn,ZB),FDBD(pface,p))
      }
      if (IS_FF(pface)){
         i = IFLAG(pface)-1;
		 in = i + n;
         if (IS_FN(n0)) {
            a[k][i] = COEFFLL(pnf,ZA,0,0);
			a[k][in] = COEFFLL(pnf,ZA,0,1);
			a[k2][i] = COEFFLL(pnf,ZA,1,0);
			a[k2][in] = COEFFLL(pnf,ZA,1,1);
		 }
         a[i][i] = COEFFNN(pface,ZA,0,0);
		 a[i][in] = COEFFNN(pface,ZA,0,1);
		 a[in][i] = COEFFNN(pface,ZA,1,0);
		 a[in][in] = COEFFNN(pface,ZA,1,1);
         b0[i] = COEFF_NF(pnf,ZB,0);
         b1[i] = COEFF_NF(pnf,ZB,1);
		 if (IS_ZNF(pface)) {
			SET9(u_rhs[i],COEFFDP(pface,ZB),NDSBD(pface,p))
			if (IS_ZNN(n0))
				SET9(u_rhs[i],COEFFDP(pnf,ZB),NDSBD(n0,p))
		 }
         for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
            pn = NBNODE(pfn);
            if (pn != n0){
               for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
               SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
			   if (IS_ZNF(pface) && IS_ZNN(pn))
				   SET9(u_rhs[i],COEFFDP(pnfl,ZB),NDSBD(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pfn,ZA,0,0);
				  a[i][jn] = COEFFLL(pfn,ZA,0,1);
				  a[in][j] = COEFFLL(pfn,ZA,1,0);
				  a[in][jn] = COEFFLL(pfn,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
			   }
         }
         for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
            if (IS_FF(pf=NBFACE(pfl))){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j] = COEFFLL(pfl,ZA,0,0);
				  a[i][jn] = COEFFLL(pfl,ZA,0,1);
				  a[in][j] = COEFFLL(pfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pfl,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
			   }
            }
      }
   }

   for (i = 0; i < n; i++){
      in = i + n;
      c[0][i] = b0[i];   	
      c[0][in] = b1[i];   	
      c[1][i] = u_rhs[i][0];
      c[1][in] = u_rhs[i][1];
   }
   sggem(a,c,y,n2,2);
   s = q = 0.;
   for (i = 0; i < n; i++){
      in = i + n;
      s += b0[i]*y[0][i] + b1[i]*y[0][in];
      q += b0[i]*y[1][i] + b1[i]*y[1][in];
   }
   pressure = (q - NDS(n0,g))/s;
   for (i = 0; i < n; i++){
      in = i + n;
      y[0][i] = y[1][i] - y[0][i]*pressure;
      y[1][i] = y[1][in] - y[0][in]*pressure;
   }

/*
   s = fabs(NDS(n0,p) - pressure);
   if (IS_FN(n0) && (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s)
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. diff. = %e\n",s);
*/

   NDS(n0,p_new) = NDS(n0,p) + om*(pressure - NDS(n0,p));
   if (IS_FN(n0)){
      ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
      ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   }
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_FN(pnode=NBNODE(pl))){
         ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
         ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
         IFLAG(pnode) = 0;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_FF(pface=NBFACE(pnfl))){
         FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
         FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
         IFLAG(pface) = 0;
      }
}

void full_Vanka_step_p2_p1_korn_BD_node(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl, *pli, *pli2;
   NFLINK *pnf, *pnfl;
   FNLINK *pfn;
   FLINK *pfl;
   FLOAT a[10][N_SGGEM], d0[5], d1[5],
	     c[2][N_SGGEM], y[2][N_SGGEM], u_rhs[5][2],
         q, s, pressure;
   INT i, in, j, jn, n=0, n2, k, k2;
   	
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_ZNN(pnode=NBNODE(pl))){
         SET1(u_rhs[n],NDD(pnode,f))
         IFLAG(pnode) = ++n;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_ZNF(pface=NBFACE(pnfl))){
         SET1(u_rhs[n],FDVP(pface,f))
         IFLAG(pface) = ++n;
      }
   SET1(u_rhs[n],NDD(n0,f))
   IFLAG(n0) = ++n;
   
   n2 = 2*n;
   k = n - 1;
   k2 = n2 - 1;
   for (i = 0; i < n; i++) {
	  d0[i] = d1[i] = 0.;
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;
   }
   for (; i < n2; i++)
	   for (j = 0; j < n2; j++)
		   a[i][j] = 0.;

   d0[k] = COEFFD(n0,ZB,0);
   d1[k] = COEFFD(n0,ZB,1);
   a[k][k]   = COEFFNN(n0,ZA,0,0);
   a[k][k2]  = COEFFNN(n0,ZA,0,1);
   a[k2][k]  = COEFFNN(n0,ZA,1,0);
   a[k2][k2] = COEFFNN(n0,ZA,1,1);
   SET9(u_rhs[k],COEFFBP(n0,ZB),NDS(n0,p))
   for (pl = TSTART(n0); pl != NULL; pl = NEXT(pl)){
      pnode = NBNODE(pl);
	  for (pli = START(pnode); n0 != NBNODE(pli); pli = NEXT(pli));
	  SET9(u_rhs[k],COEFFBLP(pli,ZB),NDS(pnode,p))
      if (IS_ZNN(pnode)){
         SET9(u_rhs[k],COEFFDP(pli,ZB),NDSBD(pnode,p))
         i = IFLAG(pnode) - 1;
		 in = i + n;
         a[k][i]   = COEFFLL(pl,ZA,0,0);
		 a[k][in]  = COEFFLL(pl,ZA,0,1);
		 a[k2][i]  = COEFFLL(pl,ZA,1,0);
		 a[k2][in] = COEFFLL(pl,ZA,1,1);
         a[i][i]   = COEFFNN(pnode,ZA,0,0);
		 a[i][in]  = COEFFNN(pnode,ZA,0,1);
		 a[in][i]  = COEFFNN(pnode,ZA,1,0);
		 a[in][in] = COEFFNN(pnode,ZA,1,1);
	     d0[i] = COEFFD(pl,ZB,0);
		 d1[i] = COEFFD(pl,ZB,1);
		 SET9(u_rhs[i],COEFFDP(pnode,ZB),NDSBD(pnode,p))
		 SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
		 SET9(u_rhs[i],COEFFBLP(pl,ZB),NDS(n0,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            if (pn != n0){
               for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
			   SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
               if (IS_ZNN(pn))
				   SET9(u_rhs[i],COEFFDP(pli2,ZB),NDSBD(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pli,ZA,0,0);
				  a[i][jn]  = COEFFLL(pli,ZA,0,1);
				  a[in][j]  = COEFFLL(pli,ZA,1,0);
				  a[in][jn] = COEFFLL(pli,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl)){
			pf=NBFACE(pnfl);
			if (IS_ZNF(pf)){
				for (pfn = FNSTART(pf); pnode != NBNODE(pfn); pfn = NEXT(pfn));
				SET9(u_rhs[i],COEFFDP(pfn,ZB),FDBD(pf,p))
			}
            if (IS_FF(pf)){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pnfl,ZA,0,0);
				  a[i][jn]  = COEFFLL(pnfl,ZA,0,1);
				  a[in][j]  = COEFFLL(pnfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pnfl,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
            } 
		 }
      } else /* NOT_ZNN(pnode) */
		  if (IS_FN(pnode))
			  MSET9(u_rhs[k],COEFFLLP(pl,ZA),NDD(pnode,u))
   }
   for (pnf = NFSTART(n0); pnf != NULL; pnf = NEXT(pnf)){
	  pface=NBFACE(pnf);
	  if (IS_ZNF(pface)) {
     	 for (pfn = FNSTART(pface); n0 != NBNODE(pfn); pfn = NEXT(pfn));
		 SET9(u_rhs[k],COEFFDP(pfn,ZB),FDBD(pface,p))
         i = IFLAG(pface)-1;
		 in = i + n;
         a[k][i]   = COEFFLL(pnf,ZA,0,0);
		 a[k][in]  = COEFFLL(pnf,ZA,0,1);
		 a[k2][i]  = COEFFLL(pnf,ZA,1,0);
		 a[k2][in] = COEFFLL(pnf,ZA,1,1);
         a[i][i]   = COEFFNN(pface,ZA,0,0);
		 a[i][in]  = COEFFNN(pface,ZA,0,1);
		 a[in][i]  = COEFFNN(pface,ZA,1,0);
		 a[in][in] = COEFFNN(pface,ZA,1,1);
		 d0[i] = COEFFD(pnf,ZB,0);
		 d1[i] = COEFFD(pnf,ZB,1);
		 SET9(u_rhs[i],COEFFDP(pface,ZB),FDBD(pface,p))
		 SET9(u_rhs[i],COEFF_NFP(pnf,ZB),NDS(n0,p))
         for (pfn = TFNSTART(pface); pfn != NULL; pfn = NEXT(pfn)){
            pn = NBNODE(pfn);
            if (pn != n0){
               for (pnfl = NFSTART(pn); pface != NBFACE(pnfl); pnfl = NEXT(pnfl));
			   SET9(u_rhs[i],COEFF_NFP(pnfl,ZB),NDS(pn,p))
               if (IS_ZNN(pn))
				   SET9(u_rhs[i],COEFFDP(pnfl,ZB),NDSBD(pn,p))
            }
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pfn,ZA,0,0);
				  a[i][jn]  = COEFFLL(pfn,ZA,0,1);
				  a[in][j]  = COEFFLL(pfn,ZA,1,0);
				  a[in][jn] = COEFFLL(pfn,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pfn,ZA),NDD(pn,u))
         }
         for (pfl = FSTART(pface); pfl != NULL; pfl = NEXT(pfl))
            if (IS_FF(pf=NBFACE(pfl))){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pfl,ZA,0,0);
				  a[i][jn]  = COEFFLL(pfl,ZA,0,1);
				  a[in][j]  = COEFFLL(pfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pfl,ZA,1,1);
               } else
                  MSET9(u_rhs[i],COEFFLLP(pfl,ZA),FDVP(pf,u))
            } 
      } else
		  if (IS_FF(pface))
			  MSET9(u_rhs[k],COEFFLLP(pnf,ZA),FDVP(pface,u))
   }

   for (i = 0; i < n; i++){
	  in = i + n;
      c[0][i]  = d0[i];   	
      c[0][in] = d1[i];   	
      c[1][i]  = u_rhs[i][0];
      c[1][in] = u_rhs[i][1];
   }
   sggem(a,c,y,n2,2);
   s = q = 0.;
   for (i = 0; i < n; i++){
	  in = i + n;
      s += d0[i]*y[0][i] + d1[i]*y[0][in];
      q += d0[i]*y[1][i] + d1[i]*y[1][in];
   }
   pressure = (q - NDSBD(n0,g))/s;
   for (i = 0; i < n; i++){
	  in = i + n;
      y[0][i] = y[1][i]  - y[0][i]*pressure;
      y[1][i] = y[1][in] - y[0][in]*pressure;
   }

/*
   s = fabs(NDSBD(n0,p) - pressure);
   if ( (q=MAX(fabs(ND(n0,u,0)-y[0][k]),fabs(ND(n0,u,1)-y[1][k]))) > s )
      s = q; 
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_ZNN(pnode=NBNODE(pl)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_ZNF(pface=NBFACE(pnfl)) &&
                            (q=MAX(fabs(FDV(pface,u,0)-y[0][IFLAG(pface)-1]),
                                   fabs(FDV(pface,u,1)-y[1][IFLAG(pface)-1]))) > s)
         s = q;
   printf("Max. n diff. = %e\n",s);
*/

   NDSBD(n0,p_new) = NDSBD(n0,p) + om*(pressure - NDSBD(n0,p));
   ND(n0,u_new,0) = ND(n0,u,0) + om*(y[0][k] - ND(n0,u,0));
   ND(n0,u_new,1) = ND(n0,u,1) + om*(y[1][k] - ND(n0,u,1));
   IFLAG(n0) = 0;
   for (pl = START(n0); pl != NULL; pl = NEXT(pl))
      if (IS_ZNN(pnode=NBNODE(pl))){
         ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
         ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
         IFLAG(pnode) = 0;
      }
   for (pnfl = NFSTART(n0); pnfl != NULL; pnfl = NEXT(pnfl))
      if (IS_ZNF(pface=NBFACE(pnfl))){
         FDV(pface,u_new,0) = FDV(pface,u,0) + om*(y[0][IFLAG(pface)-1]-FDV(pface,u,0));
         FDV(pface,u_new,1) = FDV(pface,u,1) + om*(y[1][IFLAG(pface)-1]-FDV(pface,u,1));
         IFLAG(pface) = 0;
      }
}

void full_Vanka_step_p2_p1_korn_BD_face(f0,ZA,ZB,f,g,u,p,u_new,p_new,om)
FACE *f0;
INT ZA, ZB, f, g, u, p, u_new, p_new;
FLOAT om;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pli, *pli2;
   NFLINK *pnfl;
   FNLINK *pfn, *pfn2;
   FLINK *pfl;
   FLOAT a[6][N_SGGEM], d0[3], d1[3],
	     c[2][N_SGGEM], y[2][N_SGGEM], u_rhs[3][2],
         q, s, pressure;
   INT i, in, j, jn, n=0, n2, k, k2;
   	
   for (pfn = FNSTART(f0); pfn != NULL; pfn = NEXT(pfn))
      if (IS_ZNN(pnode=NBNODE(pfn))){
         SET1(u_rhs[n],NDD(pnode,f))
         IFLAG(pnode) = ++n;
      }
   SET1(u_rhs[n],FDVP(f0,f))
   IFLAG(f0) = ++n;
   
   n2 = 2*n;
   k = n - 1;
   k2 = n2 - 1;
   for (i = 0; i < n; i++) {
	  d0[i] = d1[i] = 0.;
      for (j = 0; j < n2; j++)
         a[i][j] = 0.;
   }
   for (; i < n2; i++)
	   for (j = 0; j < n2; j++)
		   a[i][j] = 0.;

   d0[k] = COEFFD(f0,ZB,0);
   d1[k] = COEFFD(f0,ZB,1);
   a[k][k]   = COEFFNN(f0,ZA,0,0);
   a[k][k2]  = COEFFNN(f0,ZA,0,1);
   a[k2][k]  = COEFFNN(f0,ZA,1,0);
   a[k2][k2] = COEFFNN(f0,ZA,1,1);
   for (pfn = TFNSTART(f0); pfn != NULL; pfn = NEXT(pfn)){
      pnode = NBNODE(pfn);
	  for (pnfl = NFSTART(pnode); f0 != NBFACE(pnfl); pnfl = NEXT(pnfl));
      SET9(u_rhs[k],COEFF_NFP(pnfl,ZB),NDS(pnode,p))
      if (IS_ZNN(pnode)){
         i = IFLAG(pnode) - 1;
		 in = i + n;
		 SET9(u_rhs[k],COEFFDP(pnfl,ZB),NDSBD(pnode,p))
         a[k][i]   = COEFFLL(pfn,ZA,0,0);
		 a[k][in]  = COEFFLL(pfn,ZA,0,1);
		 a[k2][i]  = COEFFLL(pfn,ZA,1,0);
		 a[k2][in] = COEFFLL(pfn,ZA,1,1);
         a[i][i]   = COEFFNN(pnode,ZA,0,0);
		 a[i][in]  = COEFFNN(pnode,ZA,0,1);
		 a[in][i]  = COEFFNN(pnode,ZA,1,0);
		 a[in][in] = COEFFNN(pnode,ZA,1,1);
		 d0[i] = COEFFD(pfn,ZB,0);
		 d1[i] = COEFFD(pfn,ZB,1);
		 SET9(u_rhs[i],COEFFDP(pnode,ZB),NDSBD(pnode,p))
		 SET9(u_rhs[i],COEFFBP(pnode,ZB),NDS(pnode,p))
         for (pli = TSTART(pnode); pli != NULL; pli = NEXT(pli)){
            pn=NBNODE(pli);
            for (pli2 = START(pn); pnode != NBNODE(pli2); pli2 = NEXT(pli2));
            SET9(u_rhs[i],COEFFBLP(pli2,ZB),NDS(pn,p))
            if (IS_FN(pn))
			   if (IFLAG(pn)) {
				  SET9(u_rhs[i],COEFFDP(pli2,ZB),NDSBD(pn,p))
				  j = IFLAG(pn)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pli,ZA,0,0);
				  a[i][jn]  = COEFFLL(pli,ZA,0,1);
				  a[in][j]  = COEFFLL(pli,ZA,1,0);
				  a[in][jn] = COEFFLL(pli,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pli,ZA),NDD(pn,u))
				  if (IS_ZNN(pn))
					  SET9(u_rhs[i],COEFFDP(pli2,ZB),NDSBD(pn,p))
			   }
         }
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if (IS_FF(pf=NBFACE(pnfl))){
			   if (IFLAG(pf)) {
				  j = IFLAG(pf)-1;
				  jn = j + n;
                  a[i][j]   = COEFFLL(pnfl,ZA,0,0);
				  a[i][jn]  = COEFFLL(pnfl,ZA,0,1);
				  a[in][j]  = COEFFLL(pnfl,ZA,1,0);
				  a[in][jn] = COEFFLL(pnfl,ZA,1,1);
               } else {
                  MSET9(u_rhs[i],COEFFLLP(pnfl,ZA),FDVP(pf,u))
				  if (IS_ZNF(pf)) {
					  for (pfn2 = FNSTART(pf); pnode != NBNODE(pfn2); pfn2 = NEXT(pfn2));
					  SET9(u_rhs[i],COEFFDP(pfn2,ZB),FDBD(pf,p))
				  }
			   }
            } 
      } else /* NOT_ZNN(pnode) */
		  if (IS_FN(pnode))
			  MSET9(u_rhs[k],COEFFLLP(pfn,ZA),NDD(pnode,u))
   }
   for (pfl = FSTART(f0); pfl != NULL; pfl = NEXT(pfl))
      if (IS_FF(pface=NBFACE(pfl))){
		  MSET9(u_rhs[k],COEFFLLP(pfl,ZA),FDVP(pface,u))
      } 

   for (i = 0; i < n; i++){
	  in = i + n;
      c[0][i]  = d0[i];   	
      c[0][in] = d1[i];   	
      c[1][i]  = u_rhs[i][0];
      c[1][in] = u_rhs[i][1];
   }
   sggem(a,c,y,n2,2);
   s = q = 0.;
   for (i = 0; i < n; i++){
	  in = i + n;
      s += d0[i]*y[0][i] + d1[i]*y[0][in];
      q += d0[i]*y[1][i] + d1[i]*y[1][in];
   }
   pressure = (q - FDBD(f0,g))/s;
   for (i = 0; i < n; i++){
	  in = i + n;
      y[0][i] = y[1][i]  - y[0][i]*pressure;
      y[1][i] = y[1][in] - y[0][in]*pressure;
   }

/*
   s = fabs(FDBD(f0,p) - pressure);
   if ( (q=MAX(fabs(FDV(f0,u,0)-y[0][k]),fabs(FDV(f0,u,1)-y[1][k]))) > s )
      s = q; 
   for (pfn = FNSTART(f0); pfn != NULL; pfn = NEXT(pfn))
      if (IS_ZNN(pnode=NBNODE(pfn)) && 
                             (q=MAX(fabs(ND(pnode,u,0)-y[0][IFLAG(pnode)-1]),
                                    fabs(ND(pnode,u,1)-y[1][IFLAG(pnode)-1]))) > s)
         s = q;
   printf("Max. f diff. = %e\n",s);
*/

   FDBD(f0,p_new) = FDBD(f0,p) + om*(pressure - FDBD(f0,p));
   FDV(f0,u_new,0) = FDV(f0,u,0) + om*(y[0][k] - FDV(f0,u,0));
   FDV(f0,u_new,1) = FDV(f0,u,1) + om*(y[1][k] - FDV(f0,u,1));
   IFLAG(f0) = 0;
   for (pfn = FNSTART(f0); pfn != NULL; pfn = NEXT(pfn))
      if (IS_ZNN(pnode=NBNODE(pfn))){
         ND(pnode,u_new,0) = ND(pnode,u,0) + om*(y[0][IFLAG(pnode)-1]-ND(pnode,u,0));
         ND(pnode,u_new,1) = ND(pnode,u,1) + om*(y[1][IFLAG(pnode)-1]-ND(pnode,u,1));
         IFLAG(pnode) = 0;
      }
}

#else

void full_Vanka_step_p2_p1_korn_BD(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1_korn_BD not available.\n");  }

void full_Vanka_step_p2_p1_korn_BD_node(n0,ZA,ZB,f,g,u,p,u_new,p_new,om)
NODE *n0; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1_korn_BD_node not available.\n");  }

void full_Vanka_step_p2_p1_korn_BD_face(f0,ZA,ZB,f,g,u,p,u_new,p_new,om)
FACE *f0; INT ZA, ZB, f, g, u, p, u_new, p_new; FLOAT om;
{  eprintf("Error: full_Vanka_step_p2_p1_korn_BD_face not available.\n");  }

#endif

/******************************************************************************/
/*                     end of smoothing by Martin Sodomka                     */
/******************************************************************************/

/*--------------------------------- DIM = 3 ----------------------------------*/

#if DIM == 3

#if !(DATA_STR & KORN_MATRIX) && !(DATA_STR & REDUCED) 

void gem(a,b,x)
FLOAT a[17][17], b[17], x[17];
{
   FLOAT z, z4, z8, max=1.;
   INT p[17], i, ii, j, jj, k, i4, i8, k4, k8, ii4, ii8, jj4, jj8;

   for (i = 0; i < 17; i++)
      p[i] = i;
   for (k = 0; k < 4 && max > EPSA; k++){
      k4 = k + 4;
      k8 = k + 8;
      max = 0.0;
      for (i = k; i < 4; i++)
         for (j = k; j < 4; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            ii4 = ii + 4;
            ii8 = ii + 8;
            for (j = k; j < 4; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            for (j = 12; j < 17; j++){
               EXCHANGE(a[k ][j],a[ii ][j],z)
               EXCHANGE(a[k4][j],a[ii4][j],z)
               EXCHANGE(a[k8][j],a[ii8][j],z)
            }
            EXCHANGE(b[k ],b[ii ],z)
            EXCHANGE(b[k4],b[ii4],z)
            EXCHANGE(b[k8],b[ii8],z)
         }
         if (jj != k){
            jj4 = jj + 4;
            jj8 = jj + 8;
            EXCHANGE(p[k ],p[jj ],z)
            EXCHANGE(p[k4],p[jj4],z)
            EXCHANGE(p[k8],p[jj8],z)
            for (i = 0; i < 4; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
            for (i = 12; i < 17; i++){
               EXCHANGE(a[i][k ],a[i][jj ],z)
               EXCHANGE(a[i][k4],a[i][jj4],z)
               EXCHANGE(a[i][k8],a[i][jj8],z)
            }
         }
         for (i = k+1; i < 4; i++){
            i4 = i + 4;
            i8 = i + 8;
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 4; j++)
               a[i][j] -= z*a[k][j];
            for (j = 12; j < 17; j++){
               a[i ][j] -= z*a[k ][j];
               a[i4][j] -= z*a[k4][j];
               a[i8][j] -= z*a[k8][j];
            }
            b[i ] -= b[k ]*z;
            b[i4] -= b[k4]*z;
            b[i8] -= b[k8]*z;
         }
         for (i = 12; i < 17; i++){
            z  = a[i][k ]/a[k][k];
            z4 = a[i][k4]/a[k][k];
            z8 = a[i][k8]/a[k][k];
            for (j = k+1; j < 4; j++){
               a[i][j]   -= z *a[k][j];
               a[i][j+4] -= z4*a[k][j];
               a[i][j+8] -= z8*a[k][j];
            }
            for (j = 12; j < 17; j++)
               a[i][j] -= z*a[k][j] + z4*a[k4][j] + z8*a[k8][j];
            b[i] -= b[k]*z + b[k4]*z4 + b[k8]*z8;
         }
      }
   }
   for (k = 12; k < 16 && max > EPSA; k++){
      max = 0.0;
      for (i = k; i < 17; i++)
         for (j = k; j < 17; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            for (j = k; j < 17; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(b[k],b[ii],z)
         }
         if (jj != k){
            EXCHANGE(p[k],p[jj],i)
            for (i = 0; i < 17; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
         }
         for (i = k+1; i < 17; i++){
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 17; j++)
               a[i][j] -= z*a[k][j];
            b[i] -= b[k]*z;
         }
      }
   }
   if (max < EPSA || fabs(a[16][16]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      x[p[16]] = b[16]/a[16][16];
      x[p[15]] = (b[15] - x[p[16]]*a[15][16])/a[15][15];
      x[p[14]] = (b[14] - x[p[15]]*a[14][15] - x[p[16]]*a[14][16])/a[14][14];
      x[p[13]] = (b[13] - x[p[14]]*a[13][14] - x[p[15]]*a[13][15] - 
                                               x[p[16]]*a[13][16])/a[13][13];
      x[p[12]] = (b[12] - x[p[13]]*a[12][13] - x[p[14]]*a[12][14] - 
                          x[p[15]]*a[12][15] - x[p[16]]*a[12][16])/a[12][12];
      for (i = 3; i >= 0; i--){
         i4 = i + 4;
         i8 = i + 8;
         z  =b[i ]-x[p[12]]*a[i ][12] - x[p[13]]*a[i ][13] - x[p[14]]*a[i ][14]
                                      - x[p[15]]*a[i ][15] - x[p[16]]*a[i ][16];
         z4 =b[i4]-x[p[12]]*a[i4][12] - x[p[13]]*a[i4][13] - x[p[14]]*a[i4][14]
                                      - x[p[15]]*a[i4][15] - x[p[16]]*a[i4][16];
         z8 =b[i8]-x[p[12]]*a[i8][12] - x[p[13]]*a[i8][13] - x[p[14]]*a[i8][14]
                                      - x[p[15]]*a[i8][15] - x[p[16]]*a[i8][16];
         for (j = i+1; j < 4; j++){
            z  -= x[p[j  ]]*a[i][j];
            z4 -= x[p[j+4]]*a[i][j];
            z8 -= x[p[j+8]]*a[i][j];
         }
         x[p[i ]] = z /a[i][i];
         x[p[i4]] = z4/a[i][i];
         x[p[i8]] = z8/a[i][i];
      }
   }
}

#endif

#if !(DATA_STR & KORN_MATRIX) && (DATA_STR & REDUCED) && !(U_SPACE == P1C_NEW_FBUB) 

void gem(a,b,x)
FLOAT a[17][17], b[17], x[17];
{
   FLOAT z, z4, z8, max=1.;
   INT p[17], i, ii, j, jj, k, i4, i8, k4, k8, ii4, ii8, jj4, jj8;

   for (i = 0; i < 17; i++)
      p[i] = i;
   for (k = 0; k < 4 && max > EPSA; k++){
      k4 = k + 4;
      k8 = k + 8;
      max = 0.0;
      for (i = k; i < 4; i++)
         for (j = k; j < 4; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            ii4 = ii + 4;
            ii8 = ii + 8;
            for (j = k; j < 4; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(a[k ][16],a[ii ][16],z)
            EXCHANGE(a[k4][16],a[ii4][16],z)
            EXCHANGE(a[k8][16],a[ii8][16],z)
            EXCHANGE(b[k ],b[ii ],z)
            EXCHANGE(b[k4],b[ii4],z)
            EXCHANGE(b[k8],b[ii8],z)
         }
         if (jj != k){
            jj4 = jj + 4;
            jj8 = jj + 8;
            EXCHANGE(p[k ],p[jj ],z)
            EXCHANGE(p[k4],p[jj4],z)
            EXCHANGE(p[k8],p[jj8],z)
            for (i = 0; i < 4; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
            EXCHANGE(a[16][k ],a[16][jj ],z)
            EXCHANGE(a[16][k4],a[16][jj4],z)
            EXCHANGE(a[16][k8],a[16][jj8],z)
         }
         for (i = k+1; i < 4; i++){
            i4 = i + 4;
            i8 = i + 8;
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 4; j++)
               a[i][j] -= z*a[k][j];
            a[i ][16] -= z*a[k ][16];
            a[i4][16] -= z*a[k4][16];
            a[i8][16] -= z*a[k8][16];
            b[i ] -= b[k ]*z;
            b[i4] -= b[k4]*z;
            b[i8] -= b[k8]*z;
         }
         z  = a[16][k ]/a[k][k];
         z4 = a[16][k4]/a[k][k];
         z8 = a[16][k8]/a[k][k];
         for (j = k+1; j < 4; j++){
            a[16][j]   -= z *a[k][j];
            a[16][j+4] -= z4*a[k][j];
            a[16][j+8] -= z8*a[k][j];
         }
         a[16][16] -= z*a[k][16] + z4*a[k4][16] + z8*a[k8][16];
         b[16] -= b[k]*z + b[k4]*z4 + b[k8]*z8;
      }
   }
   for (k = 12; k < 16 && max > EPSA; k++){
      max = 0.0;
      for (i = k; i < 17; i++)
         for (j = k; j < 17; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            for (j = k; j < 17; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(b[k],b[ii],z)
         }
         if (jj != k){
            EXCHANGE(p[k],p[jj],i)
            for (i = 0; i < 17; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
         }
         for (i = k+1; i < 17; i++){
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 17; j++)
               a[i][j] -= z*a[k][j];
            b[i] -= b[k]*z;
         }
      }
   }
   if (max < EPSA || fabs(a[16][16]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      x[p[16]] = b[16]/a[16][16];
      x[p[15]] = (b[15] - x[p[16]]*a[15][16])/a[15][15];
      x[p[14]] = (b[14] - x[p[15]]*a[14][15] - x[p[16]]*a[14][16])/a[14][14];
      x[p[13]] = (b[13] - x[p[14]]*a[13][14] - x[p[15]]*a[13][15] - 
                                               x[p[16]]*a[13][16])/a[13][13];
      x[p[12]] = (b[12] - x[p[13]]*a[12][13] - x[p[14]]*a[12][14] - 
                          x[p[15]]*a[12][15] - x[p[16]]*a[12][16])/a[12][12];
      for (i = 3; i >= 0; i--){
         i4 = i + 4;
         i8 = i + 8;
         z  = b[i ] - x[p[16]]*a[i ][16];
         z4 = b[i4] - x[p[16]]*a[i4][16];
         z8 = b[i8] - x[p[16]]*a[i8][16];
         for (j = i+1; j < 4; j++){
            z  -= x[p[j  ]]*a[i][j];
            z4 -= x[p[j+4]]*a[i][j];
            z8 -= x[p[j+8]]*a[i][j];
         }
         x[p[i ]] = z /a[i][i];
         x[p[i4]] = z4/a[i][i];
         x[p[i8]] = z8/a[i][i];
      }
   }
}

#endif

#if !(DATA_STR & KORN_MATRIX) && (DATA_STR & REDUCED) && (U_SPACE == P1C_NEW_FBUB) 

void gem(a,b,x)
FLOAT a[17][17], b[17], x[17];
{
   FLOAT z, z4, z8, max=1.;
   INT p[17], i, ii, j, jj, k, i4, i8, k4, k8, ii4, ii8, jj4, jj8;

   for (i = 0; i < 17; i++)
      p[i] = i;
   for (k = 0; k < 4 && max > EPSA; k++){
      k4 = k + 4;
      k8 = k + 8;
      max = 0.0;
      for (i = k; i < 4; i++)
         for (j = k; j < 4; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            ii4 = ii + 4;
            ii8 = ii + 8;
            for (j = k; j < 4; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(a[k ][16],a[ii ][16],z)
            EXCHANGE(a[k4][16],a[ii4][16],z)
            EXCHANGE(a[k8][16],a[ii8][16],z)
            EXCHANGE(b[k ],b[ii ],z)
            EXCHANGE(b[k4],b[ii4],z)
            EXCHANGE(b[k8],b[ii8],z)
         }
         if (jj != k){
            jj4 = jj + 4;
            jj8 = jj + 8;
            EXCHANGE(p[k ],p[jj ],z)
            EXCHANGE(p[k4],p[jj4],z)
            EXCHANGE(p[k8],p[jj8],z)
            for (i = 0; i < 4; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
            EXCHANGE(a[16][k ],a[16][jj ],z)
            EXCHANGE(a[16][k4],a[16][jj4],z)
            EXCHANGE(a[16][k8],a[16][jj8],z)
         }
         for (i = k+1; i < 4; i++){
            i4 = i + 4;
            i8 = i + 8;
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 4; j++)
               a[i][j] -= z*a[k][j];
            a[i ][16] -= z*a[k ][16];
            a[i4][16] -= z*a[k4][16];
            a[i8][16] -= z*a[k8][16];
            b[i ] -= b[k ]*z;
            b[i4] -= b[k4]*z;
            b[i8] -= b[k8]*z;
         }
         z  = a[16][k ]/a[k][k];
         z4 = a[16][k4]/a[k][k];
         z8 = a[16][k8]/a[k][k];
         for (j = k+1; j < 4; j++){
            a[16][j]   -= z *a[k][j];
            a[16][j+4] -= z4*a[k][j];
            a[16][j+8] -= z8*a[k][j];
         }
         a[16][16] -= z*a[k][16] + z4*a[k4][16] + z8*a[k8][16];
         b[16] -= b[k]*z + b[k4]*z4 + b[k8]*z8;
      }
   }
   for (k = 12; k < 16 && max > EPSA; k++){
      if ((max=fabs(a[k][k])) > EPSA){
            z = a[16][k]/a[k][k];
            a[16][16] -= z*a[k][16];
            b[16] -= b[k]*z;
      }
   }
   if (max < EPSA || fabs(a[16][16]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      x[p[16]] = b[16]/a[16][16];
      x[p[15]] = (b[15] - x[p[16]]*a[15][16])/a[15][15];
      x[p[14]] = (b[14] - x[p[16]]*a[14][16])/a[14][14];
      x[p[13]] = (b[13] - x[p[16]]*a[13][16])/a[13][13];
      x[p[12]] = (b[12] - x[p[16]]*a[12][16])/a[12][12];
      for (i = 3; i >= 0; i--){
         i4 = i + 4;
         i8 = i + 8;
         z  = b[i ] - x[p[16]]*a[i ][16];
         z4 = b[i4] - x[p[16]]*a[i4][16];
         z8 = b[i8] - x[p[16]]*a[i8][16];
         for (j = i+1; j < 4; j++){
            z  -= x[p[j  ]]*a[i][j];
            z4 -= x[p[j+4]]*a[i][j];
            z8 -= x[p[j+8]]*a[i][j];
         }
         x[p[i ]] = z /a[i][i];
         x[p[i4]] = z4/a[i][i];
         x[p[i8]] = z8/a[i][i];
      }
   }
}

#endif

#if DATA_STR & KORN_MATRIX

void gem(a,b,x)
FLOAT a[17][17], b[17], x[17];
{
   FLOAT z, max=1.;
   INT p[17], i, ii, j, jj, k;

   for (i = 0; i < 17; i++)
      p[i] = i;
   for (k = 0; k < 16 && max > EPSA; k++){
      max = 0.0;
      for (i = k; i < 17; i++)
         for (j = k; j < 17; j++)
            if (fabs(a[i][j]) > max){
               max = fabs(a[i][j]);
               ii = i;
               jj = j;
            }
      if (max > EPSA){
         if (ii != k){
            for (j = k; j < 17; j++)
               EXCHANGE(a[k][j],a[ii][j],z)
            EXCHANGE(b[k],b[ii],z)
         }
         if (jj != k){
            EXCHANGE(p[k],p[jj],i)
            for (i = 0; i < 17; i++)
               EXCHANGE(a[i][k],a[i][jj],z)
         }
         for (i = k+1; i < 17; i++){
            z = a[i][k]/a[k][k];
            for (j = k+1; j < 17; j++)
               a[i][j] -= z*a[k][j];
            b[i] -= b[k]*z;
         }
      }
   }
   if (max < EPSA || fabs(a[16][16]) < EPSA)
      eprintf("Matrix is singular.\n");
   else{
      for (i = 16; i >= 0; i--){
         z = b[i];
         for (j = i+1; j < 17; j++)
            z -= x[p[j]]*a[i][j];
         x[p[i]] = z/a[i][i];
      }
   }
}

#endif

#if !(DATA_STR & KORN_MATRIX)

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & ExDN_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

#if !(DATA_STR & REDUCED)

void Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   FLOAT sum1, sum2, sum3, sumf, q, nrhs[4][3], frhs[4];
   INT i, j;
   	
   for (i = 0; i < 4; i++){
      sum1 = sum2 = sum3 = sumf = 0.;
   	   
      pnode = pel->n[i];
      for (pl=START(pnode); pl!=NULL; pl=NEXT(pl)){
         pn = NBNODE(pl);
	 sum1 += COEFFL(pl,Z)*ND(pn,x,0);
         sum2 += COEFFL(pl,Z)*ND(pn,x,1);
         sum3 += COEFFL(pl,Z)*ND(pn,x,2);
      }
      for (pnfl=NFSTART(pnode); pnfl != NULL; pnfl=NEXT(pnfl)){
         sum1 += COEFF_NF(pnfl,Z,0)*FD(NBFACE(pnfl),x);
         sum2 += COEFF_NF(pnfl,Z,1)*FD(NBFACE(pnfl),x);
         sum3 += COEFF_NF(pnfl,Z,2)*FD(NBFACE(pnfl),x);
      } 
      nrhs[i][0] = ND(pnode,f,0) - sum1;
      nrhs[i][1] = ND(pnode,f,1) - sum2;
      nrhs[i][2] = ND(pnode,f,2) - sum3;
      
      pface = pel->f[i];
      for (pfl=FSTART(pface); pfl != NULL; pfl=NEXT(pfl))
	 sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      for (pfnl=FNSTART(pface); pfnl != NULL; pfnl=NEXT(pfnl))
	 sumf += COEFF_FN(pfnl,Z,0)*ND(NBNODE(pfnl),x,0) +
	      	 COEFF_FN(pfnl,Z,1)*ND(NBNODE(pfnl),x,1) +
	      	 COEFF_FN(pfnl,Z,2)*ND(NBNODE(pfnl),x,2); 
      frhs[i] = FD(pface,f) - sumf;
   }
   p_to_rhs_diag(pel,nrhs,frhs,xe);
   sum1 = 0.;
   sum2 = -ED(pel,fe);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++){
            q = COEFF_BN(pel,0,i,j)/COEFFN(pel->n[i],Z);
            sum1 += COEFF_BN(pel,0,i,j)*q;
            sum2 += nrhs[i][j]*q;
         }
      if (IS_FF(pel->f[i])){
         q = COEFF_BF(pel,0,i)/COEFF_FF(pel->f[i],Z);
         sum1 += COEFF_BF(pel,0,i)*q;
         sum2 += frhs[i]*q;
      }
   }
   q = ED(pel,xe) = sum2/sum1;
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++)
            ND(pel->n[i],x,j) = 
                        (nrhs[i][j] - COEFF_BN(pel,0,i,j)*q)/COEFFN(pel->n[i],Z);
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = (frhs[i] - COEFF_BF(pel,0,i)*q)/COEFF_FF(pel->f[i],Z);
   }
}

void Vanka_step_stab(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   FLOAT sum1, sum2, sum3, sumf, a[17][17], rhs[17], u[17];
   INT i, j, i4, i8, i12;
   	
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
         a[i][j] = 0.;
   for (i = 0; i < 12; i++)
      for (j = 12; j < 17; j++)
         a[i][j] = a[j][i] = 0.;
   for (i = 12; i < 17; i++)
      for (j = 12; j < 17; j++)
         a[i][j] = 0.;
   
   for (i = 0; i < 4; i++){
      i4  = i + 4;
      i8  = i + 8;
      i12 = i + 12;
      sum1 = sum2 = sum3 = sumf = 0.;
         
      if (IS_FN(pnode=pel->n[i])){
         a[i][i] = COEFFN(pnode,Z);
         for (pl = START(pnode); pl != NULL; pl = NEXT(pl))
            if ((pn=NBNODE(pl)) == pel->n[0])
               a[i][0] = COEFFL(pl,Z);
            else if (pn == pel->n[1])
               a[i][1] = COEFFL(pl,Z);
            else if (pn == pel->n[2])
               a[i][2] = COEFFL(pl,Z);
            else if (pn == pel->n[3])
               a[i][3] = COEFFL(pl,Z);
            else{
               sum1 += COEFFL(pl,Z)*ND(pn,x,0);
               sum2 += COEFFL(pl,Z)*ND(pn,x,1);
               sum3 += COEFFL(pl,Z)*ND(pn,x,2);
            }
         
         for (pnfl = NFSTART(pnode); pnfl != NULL; pnfl = NEXT(pnfl))
            if ((pf=NBFACE(pnfl)) == pel->f[0]){
               a[i ][12] = COEFF_NF(pnfl,Z,0);
               a[i4][12] = COEFF_NF(pnfl,Z,1);
               a[i8][12] = COEFF_NF(pnfl,Z,2);
            }
            else if (pf == pel->f[1]){
               a[i ][13] = COEFF_NF(pnfl,Z,0);
               a[i4][13] = COEFF_NF(pnfl,Z,1);
               a[i8][13] = COEFF_NF(pnfl,Z,2);
            }
            else if (pf == pel->f[2]){
               a[i ][14] = COEFF_NF(pnfl,Z,0);
               a[i4][14] = COEFF_NF(pnfl,Z,1);
               a[i8][14] = COEFF_NF(pnfl,Z,2);
            }
            else if (pf == pel->f[3]){
               a[i ][15] = COEFF_NF(pnfl,Z,0);
               a[i4][15] = COEFF_NF(pnfl,Z,1);
               a[i8][15] = COEFF_NF(pnfl,Z,2);
            }
            else{
               sum1 += COEFF_NF(pnfl,Z,0)*FD(pf,x);
               sum2 += COEFF_NF(pnfl,Z,1)*FD(pf,x);
               sum3 += COEFF_NF(pnfl,Z,2)*FD(pf,x);
            }
         a[i ][16] = a[16][i ] = COEFF_BN(pel,0,i,0);
         a[i4][16] = a[16][i4] = COEFF_BN(pel,0,i,1);
         a[i8][16] = a[16][i8] = COEFF_BN(pel,0,i,2);
         rhs[i ] = ND(pnode,f,0) - sum1;
         rhs[i4] = ND(pnode,f,1) - sum2;
         rhs[i8] = ND(pnode,f,2) - sum3;
      }
      else{
         a[i][i] = 1.;
         rhs[i ] = 0.;
         rhs[i4] = 0.;
         rhs[i8] = 0.;
      }
      
      if (IS_FF(pface=pel->f[i])){
         a[i12][i12] = COEFF_FF(pface,Z);
         for (pfl=FSTART(pface); pfl != NULL; pfl=NEXT(pfl))
            if ((pf=NBFACE(pfl)) == pel->f[0])
               a[i12][12] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[1])
               a[i12][13] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[2])
               a[i12][14] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[3])
               a[i12][15] = COEFF_FL(pfl,Z);
            else
	       sumf += COEFF_FL(pfl,Z)*FD(pf,x);
         for (pfnl=FNSTART(pface); pfnl != NULL; pfnl=NEXT(pfnl))
            if ((pn=NBNODE(pfnl)) == pel->n[0]){ 
               a[i12][0] = COEFF_FN(pfnl,Z,0);
               a[i12][4] = COEFF_FN(pfnl,Z,1);
               a[i12][8] = COEFF_FN(pfnl,Z,2);
            }
            else if(pn == pel->n[1]){
               a[i12][1] = COEFF_FN(pfnl,Z,0);
               a[i12][5] = COEFF_FN(pfnl,Z,1);
               a[i12][9] = COEFF_FN(pfnl,Z,2);
            }
            else if(pn == pel->n[2]){
               a[i12][2 ] = COEFF_FN(pfnl,Z,0);
               a[i12][6 ] = COEFF_FN(pfnl,Z,1);
               a[i12][10] = COEFF_FN(pfnl,Z,2);
            }
            else if(pn == pel->n[3]){
               a[i12][3 ] = COEFF_FN(pfnl,Z,0);
               a[i12][7 ] = COEFF_FN(pfnl,Z,1);
               a[i12][11] = COEFF_FN(pfnl,Z,2);
            }
            else
               sumf += COEFF_FN(pfnl,Z,0)*ND(pn,x,0) +
	               COEFF_FN(pfnl,Z,1)*ND(pn,x,1) +
	               COEFF_FN(pfnl,Z,2)*ND(pn,x,2);
	 a[i12][16] = a[16][i12] = COEFF_BF(pel,0,i);
         rhs[i12] = FD(pface,f) - sumf;
      }
      else{
         a[i12][i12] = 1.;
         rhs[i12] = 0.;
      }
   }
   p_to_rhs_stab(pel,rhs,xe);
   rhs[16] = ED(pel,fe);
   gem(a,rhs,u);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i])){
         ND(pel->n[i],x,0) = u[i];
         ND(pel->n[i],x,1) = u[i+4];
         ND(pel->n[i],x,2) = u[i+8];
      }
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = u[i+12];
   }
   ED(pel,xe) = u[16];   
}

#endif

#if (DATA_STR & REDUCED) && !(U_SPACE == P1C_NEW_FBUB)

void Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface;
   LINK *pl;
   FLINK *pfl;
   FLOAT sum1, sum2, sum3, sumf, q, nrhs[4][3], frhs[4];
   INT i, j;
   	
   for (i = 0; i < 4; i++){
      sum1 = sum2 = sum3 = sumf = 0.;
   	   
      pnode = pel->n[i];
      for (pl=START(pnode); pl!=NULL; pl=NEXT(pl)){
         pn = NBNODE(pl);
	 sum1 += COEFFL(pl,Z)*ND(pn,x,0);
         sum2 += COEFFL(pl,Z)*ND(pn,x,1);
         sum3 += COEFFL(pl,Z)*ND(pn,x,2);
      }
      nrhs[i][0] = ND(pnode,f,0) - sum1;
      nrhs[i][1] = ND(pnode,f,1) - sum2;
      nrhs[i][2] = ND(pnode,f,2) - sum3;
      
      pface = pel->f[i];
      for (pfl=FSTART(pface); pfl != NULL; pfl=NEXT(pfl))
	 sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      frhs[i] = FD(pface,f) - sumf;
   }
   p_to_rhs_diag(pel,nrhs,frhs,xe);
   sum1 = 0.;
   sum2 = -ED(pel,fe);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++){
            q = COEFF_BN(pel,0,i,j)/COEFFN(pel->n[i],Z);
            sum1 += COEFF_BN(pel,0,i,j)*q;
            sum2 += nrhs[i][j]*q;
         }
      if (IS_FF(pel->f[i])){
         q = COEFF_BF(pel,0,i)/COEFF_FF(pel->f[i],Z);
         sum1 += COEFF_BF(pel,0,i)*q;
         sum2 += frhs[i]*q;
      }
   }
   q = ED(pel,xe) = sum2/sum1;
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++)
            ND(pel->n[i],x,j) = 
                        (nrhs[i][j] - COEFF_BN(pel,0,i,j)*q)/COEFFN(pel->n[i],Z);
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = (frhs[i] - COEFF_BF(pel,0,i)*q)/COEFF_FF(pel->f[i],Z);
   }
}

void Vanka_step_stab(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl;
   FLINK *pfl;
   FLOAT sum1, sum2, sum3, sumf, a[17][17], rhs[17], u[17];
   INT i, j, i4, i8, i12;
   	
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
         a[i][j] = 0.;
   for (i = 0; i < 12; i++)
      a[i][16] = a[16][i] = 0.;
   for (i = 12; i < 17; i++)
      for (j = 12; j < 17; j++)
         a[i][j] = 0.;
   
   for (i = 0; i < 4; i++){
      i4  = i + 4;
      i8  = i + 8;
      i12 = i + 12;
      sum1 = sum2 = sum3 = sumf = 0.;
         
      if (IS_FN(pnode=pel->n[i])){
         a[i][i] = COEFFN(pnode,Z);
         for (pl = START(pnode); pl != NULL; pl = NEXT(pl))
            if ((pn=NBNODE(pl)) == pel->n[0])
               a[i][0] = COEFFL(pl,Z);
            else if (pn == pel->n[1])
               a[i][1] = COEFFL(pl,Z);
            else if (pn == pel->n[2])
               a[i][2] = COEFFL(pl,Z);
            else if (pn == pel->n[3])
               a[i][3] = COEFFL(pl,Z);
            else{
               sum1 += COEFFL(pl,Z)*ND(pn,x,0);
               sum2 += COEFFL(pl,Z)*ND(pn,x,1);
               sum3 += COEFFL(pl,Z)*ND(pn,x,2);
            }
         a[i ][16] = a[16][i ] = COEFF_BN(pel,0,i,0);
         a[i4][16] = a[16][i4] = COEFF_BN(pel,0,i,1);
         a[i8][16] = a[16][i8] = COEFF_BN(pel,0,i,2);
         rhs[i ] = ND(pnode,f,0) - sum1;
         rhs[i4] = ND(pnode,f,1) - sum2;
         rhs[i8] = ND(pnode,f,2) - sum3;
      }
      else{
         a[i][i] = 1.;
         rhs[i ] = 0.;
         rhs[i4] = 0.;
         rhs[i8] = 0.;
      }
      
      if (IS_FF(pface=pel->f[i])){
         a[i12][i12] = COEFF_FF(pface,Z);
         for (pfl=FSTART(pface); pfl != NULL; pfl=NEXT(pfl))
            if ((pf=NBFACE(pfl)) == pel->f[0])
               a[i12][12] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[1])
               a[i12][13] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[2])
               a[i12][14] = COEFF_FL(pfl,Z);
            else if (pf == pel->f[3])
               a[i12][15] = COEFF_FL(pfl,Z);
            else
	       sumf += COEFF_FL(pfl,Z)*FD(pf,x);
	 a[i12][16] = a[16][i12] = COEFF_BF(pel,0,i);
         rhs[i12] = FD(pface,f) - sumf;
      }
      else{
         a[i12][i12] = 1.;
         rhs[i12] = 0.;
      }
   }
   p_to_rhs_stab(pel,rhs,xe);
   rhs[16] = ED(pel,fe);
   gem(a,rhs,u);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i])){
         ND(pel->n[i],x,0) = u[i];
         ND(pel->n[i],x,1) = u[i+4];
         ND(pel->n[i],x,2) = u[i+8];
      }
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = u[i+12];
   }
   ED(pel,xe) = u[16];   
}

#endif

#if (DATA_STR & REDUCED) && (U_SPACE == P1C_NEW_FBUB)

void Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   LINK *pl;
   FLOAT sum1, sum2, sum3, q, nrhs[4][3], frhs[4];
   INT i, j;
   	
   for (i = 0; i < 4; i++){
      sum1 = sum2 = sum3 = 0.;
      pnode = pel->n[i];
      for (pl=START(pnode); pl!=NULL; pl=NEXT(pl)){
         pn = NBNODE(pl);
	 sum1 += COEFFL(pl,Z)*ND(pn,x,0);
         sum2 += COEFFL(pl,Z)*ND(pn,x,1);
         sum3 += COEFFL(pl,Z)*ND(pn,x,2);
      }
      nrhs[i][0] = ND(pnode,f,0) - sum1;
      nrhs[i][1] = ND(pnode,f,1) - sum2;
      nrhs[i][2] = ND(pnode,f,2) - sum3;
      frhs[i] = FD(pel->f[i],f);
   }
   p_to_rhs_diag(pel,nrhs,frhs,xe);
   sum1 = 0.;
   sum2 = -ED(pel,fe);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++){
            q = COEFF_BN(pel,0,i,j)/COEFFN(pel->n[i],Z);
            sum1 += COEFF_BN(pel,0,i,j)*q;
            sum2 += nrhs[i][j]*q;
         }
      if (IS_FF(pel->f[i])){
         q = COEFF_BF(pel,0,i)/COEFF_FF(pel->f[i],Z);
         sum1 += COEFF_BF(pel,0,i)*q;
         sum2 += frhs[i]*q;
      }
   }
   q = ED(pel,xe) = sum2/sum1;
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++)
            ND(pel->n[i],x,j) = 
                        (nrhs[i][j] - COEFF_BN(pel,0,i,j)*q)/COEFFN(pel->n[i],Z);
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = (frhs[i] - COEFF_BF(pel,0,i)*q)/COEFF_FF(pel->f[i],Z);
   }
}

void Vanka_step_stab(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface, *pf;
   LINK *pl;
   FLOAT sum1, sum2, sum3, a[17][17], rhs[17], u[17];
   INT i, j, i4, i8, i12;
   	
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
         a[i][j] = 0.;
   for (i = 0; i < 16; i++)
      a[i][16] = a[16][i] = 0.;
   a[12][12] = a[13][13] = a[14][14] = a[15][15] = a[16][16] = 0.;
   
   for (i = 0; i < 4; i++){
      i4  = i + 4;
      i8  = i + 8;
      i12 = i + 12;
      sum1 = sum2 = sum3 = 0.;
         
      if (IS_FN(pnode=pel->n[i])){
         a[i][i] = COEFFN(pnode,Z);
         for (pl = START(pnode); pl != NULL; pl = NEXT(pl))
            if ((pn=NBNODE(pl)) == pel->n[0])
               a[i][0] = COEFFL(pl,Z);
            else if (pn == pel->n[1])
               a[i][1] = COEFFL(pl,Z);
            else if (pn == pel->n[2])
               a[i][2] = COEFFL(pl,Z);
            else if (pn == pel->n[3])
               a[i][3] = COEFFL(pl,Z);
            else{
               sum1 += COEFFL(pl,Z)*ND(pn,x,0);
               sum2 += COEFFL(pl,Z)*ND(pn,x,1);
               sum3 += COEFFL(pl,Z)*ND(pn,x,2);
            }
         a[i ][16] = a[16][i ] = COEFF_BN(pel,0,i,0);
         a[i4][16] = a[16][i4] = COEFF_BN(pel,0,i,1);
         a[i8][16] = a[16][i8] = COEFF_BN(pel,0,i,2);
         rhs[i ] = ND(pnode,f,0) - sum1;
         rhs[i4] = ND(pnode,f,1) - sum2;
         rhs[i8] = ND(pnode,f,2) - sum3;
      }
      else{
         a[i][i] = 1.;
         rhs[i ] = 0.;
         rhs[i4] = 0.;
         rhs[i8] = 0.;
      }
      
      if (IS_FF(pface=pel->f[i])){
         a[i12][i12] = COEFF_FF(pface,Z);
	 a[i12][16] = a[16][i12] = COEFF_BF(pel,0,i);
         rhs[i12] = FD(pface,f);
      }
      else{
         a[i12][i12] = 1.;
         rhs[i12] = 0.;
      }
   }
   p_to_rhs_stab(pel,rhs,xe);
   rhs[16] = ED(pel,fe);
   gem(a,rhs,u);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i])){
         ND(pel->n[i],x,0) = u[i];
         ND(pel->n[i],x,1) = u[i+4];
         ND(pel->n[i],x,2) = u[i+8];
      }
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = u[i+12];
   }
   ED(pel,xe) = u[16];   
}

#endif

#else

void Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid; ELEMENT *pel; INT Z, f, fe, x, xe;
{  eprintf("Error: Vanka_step_diag not available.\n");  }

#endif

#endif

#if DATA_STR & KORN_MATRIX

#if (N_DATA & DxD_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & ExDN_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid;
ELEMENT *pel;
INT Z, f, fe, x, xe;
{
   NODE *pnode, *pn;
   FACE *pface;
   LINK *pl;
   NFLINK *pnfl;
   FLINK *pfl;
   FNLINK *pfnl;
   FLOAT sum1, sum2, sum3, sumf, q, nrhs[4][3], frhs[4];
   INT i, j;
   	
   for (i = 0; i < 4; i++){
      pnode = pel->n[i];
      sum1 = COEFFNN(pnode,Z,0,1)*ND(pnode,x,1) + 
             COEFFNN(pnode,Z,0,2)*ND(pnode,x,2);
      sum2 = COEFFNN(pnode,Z,1,0)*ND(pnode,x,0) + 
             COEFFNN(pnode,Z,1,2)*ND(pnode,x,2);
      sum3 = COEFFNN(pnode,Z,2,0)*ND(pnode,x,0) + 
             COEFFNN(pnode,Z,2,1)*ND(pnode,x,1);
   	   
      for (pl=START(pnode); pl!=NULL; pl=NEXT(pl)){
         pn = NBNODE(pl);
	 sum1 += COEFFLL(pl,Z,0,0)*ND(pn,x,0) + COEFFLL(pl,Z,0,1)*ND(pn,x,1) + 
	                                        COEFFLL(pl,Z,0,2)*ND(pn,x,2);
         sum2 += COEFFLL(pl,Z,1,0)*ND(pn,x,0) + COEFFLL(pl,Z,1,1)*ND(pn,x,1) + 
                                                COEFFLL(pl,Z,1,2)*ND(pn,x,2);
         sum3 += COEFFLL(pl,Z,2,0)*ND(pn,x,0) + COEFFLL(pl,Z,2,1)*ND(pn,x,1) + 
                                                COEFFLL(pl,Z,2,2)*ND(pn,x,2);
      }
      for (pnfl=NFSTART(pnode); pnfl != NULL; pnfl=NEXT(pnfl)){
         sum1 += COEFF_NF(pnfl,Z,0)*FD(NBFACE(pnfl),x);
         sum2 += COEFF_NF(pnfl,Z,1)*FD(NBFACE(pnfl),x);
         sum3 += COEFF_NF(pnfl,Z,2)*FD(NBFACE(pnfl),x);
      } 
      nrhs[i][0] = ND(pnode,f,0) - sum1;
      nrhs[i][1] = ND(pnode,f,1) - sum2;
      nrhs[i][2] = ND(pnode,f,2) - sum3;
      
      pface = pel->f[i];
      sumf = 0.;
      for (pfl=FSTART(pface); pfl != NULL; pfl=NEXT(pfl))
	 sumf += COEFF_FL(pfl,Z)*FD(NBFACE(pfl),x);
      for (pfnl=FNSTART(pface); pfnl != NULL; pfnl=NEXT(pfnl))
	 sumf += COEFF_FN(pfnl,Z,0)*ND(NBNODE(pfnl),x,0) +
	      	 COEFF_FN(pfnl,Z,1)*ND(NBNODE(pfnl),x,1) +
	      	 COEFF_FN(pfnl,Z,2)*ND(NBNODE(pfnl),x,2); 
      frhs[i] = FD(pface,f) - sumf;
   }
   p_to_rhs_diag(pel,nrhs,frhs,xe);
   sum1 = 0.;
   sum2 = -ED(pel,fe);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++){
            q = COEFF_BN(pel,0,i,j)/COEFFNN(pel->n[i],Z,j,j);
            sum1 += COEFF_BN(pel,0,i,j)*q;
            sum2 += nrhs[i][j]*q;
         }
      if (IS_FF(pel->f[i])){
         q = COEFF_BF(pel,0,i)/COEFF_FF(pel->f[i],Z);
         sum1 += COEFF_BF(pel,0,i)*q;
         sum2 += frhs[i]*q;
      }
   }
   q = ED(pel,xe) = sum2/sum1;
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         for (j = 0; j < 3; j++)
            ND(pel->n[i],x,j) = 
                   (nrhs[i][j] - COEFF_BN(pel,0,i,j)*q)/COEFFNN(pel->n[i],Z,j,j);
      if (IS_FF(pel->f[i]))
         FD(pel->f[i],x) = (frhs[i] - COEFF_BF(pel,0,i)*q)/COEFF_FF(pel->f[i],Z);
   }
}

void Vanka_step_stab(tGrid,pel,Z,f,fe,x,xe)
GRID *tGrid; ELEMENT *pel; INT Z, f, fe, x, xe;
{  eprintf("Error: Vanka_step_stab not available.\n");  }

#endif

#endif

void Vanka_smoother1(tGrid,Z,f,fe,u,ue,x,xe,om)
GRID *tGrid;
FLOAT om;
INT Z, f, fe, u, ue, x, xe;
{
   ELEMENT *pel;

/*   printf(" Smoothing.\n");*/
   vs_copy_nf(tGrid,u,x,1);
   copy_e(tGrid,ue,xe);
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      Vanka_step_diag(tGrid,pel,Z,f,fe,x,xe);
   vs_damp_nf(tGrid,om,u,x,u,1);
   damp_e(tGrid,om,ue,xe,ue);
   if (T_FOR_P & WITHOUT_FIRST)
      subtr_first_element(tGrid,ue);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

FLOAT kontrola(tGrid,pel,x,xe)  /* v multB a multBT se nesmi vynechavat radek */
GRID *tGrid;
ELEMENT *pel;
INT x, xe;
{
   INT i;
   FLOAT sum=0.;
   
   Stokes_defect(tGrid,A,A,F,EF,x,xe,R,ER,Q,T_FOR_U,T_FOR_P,
                 U_TYPE,P_TYPE,A_STRUCT,B_STRUCT,C_STRUCT);
   for (i = 0; i < 4; i++){
      if (IS_FN(pel->n[i]))
         sum += fabs(ND(pel->n[i],R,0)) + fabs(ND(pel->n[i],R,0)) + 
                                          fabs(ND(pel->n[i],R,0)); 
      if (IS_FF(pel->f[i]))
         sum += fabs(FD(pel->f[i],R));
   }
   sum += fabs(ED(pel,ER));
   return(sum);
}

void Vanka_smoother3(tGrid,Z,f,fe,u,ue,x,xe,om)
GRID *tGrid;
FLOAT om;
INT Z, f, fe, u, ue, x, xe;
{
   ELEMENT *pel;
   FLOAT max=0., sum;

   printf(" Smoothing.\n");
   vs_copy_nf(tGrid,u,x,1);
   copy_e(tGrid,ue,xe);
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      Vanka_step_stab(tGrid,pel,Z,f,fe,x,xe);
      sum = kontrola(tGrid,pel,x,xe);
      if (sum > max) max = sum;
   }
   printf("max = %e\n",max);
   vs_damp_nf(tGrid,om,u,x,u,1);
   damp_e(tGrid,om,ue,xe,ue);
   if (T_FOR_P & WITHOUT_FIRST)
      subtr_first_element(tGrid,ue);
}

#endif

#endif

/*------------------------------ end of DIM = 3 ------------------------------*/

#if DIM == 2

void Vanka_smoother1(tGrid,Z,f,fe,u,ue,x,xe,om)
GRID *tGrid; FLOAT om; INT Z, f, fe, u, ue, x, xe;
{  eprintf("Error: Vanka_smoother1 not available.\n");  }

#endif

INT Vanka_smoother(tGrid,i,omega,ZA,ZB,j0,f,g,u,p,      j1,j2,j3,j4,j5,j6,j7,j8,
                   t_u,t_p,u_type,p_type,A_struct,B_struct, j9,j10,j11,j12,proc)
GRID *tGrid;
FLOAT omega;
INT i, ZA, ZB, f, g, u, p, 
    j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12;
INT t_u, t_p, u_type, p_type, A_struct, B_struct;
void (*proc)();
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   SNODE *psnode;
   SFACE *psface;
   INT u_new, p_new, k;

   if (i == 1 || i == 3) 
      printf(" Smoothing.\n");
   u_new = u;
   p_new = p;
   switch(u_type){
   case Q_VN: if (p_type == Q_SE){
//                 copy(tGrid,g,R,t_p,p_type);
//                 stab_eq_rhs_to_div_eq_rhs(tGrid,ZA,ZB,f,R,p);
                 for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
                    full_Vanka_step_p1c_p0_div_stab(pel,ZA,ZB,ZA,f,g,u,p,
                                                             u_new,p_new,omega);
//                 compute_stab_velocity(tGrid,ZA,ZB,u,p,f);
              }
              else
                 eprintf("Error: Vanka_smoother not available.\n");
        break;
   case Q_VF: if (p_type == Q_SE){
                   if (A_struct & Q_FULL)
                      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                         full_Vanka_step_p1nc_p0_Korn(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
//                         half_full_Vanka_step_p1nc_p0_Korn(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                   else if (A_struct & Q_FEBDIAG)
                      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                         full_Vanka_step_p1nc_p0(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                   else
                      eprintf("Error: Vanka_smoother not available.\n");
                }
                else
                   eprintf("Error: Vanka_smoother not available.\n");
        break;
   case Q_VNSF: if (p_type == Q_SE)
                   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
                      full_Vanka_step_p1c_new_fbub_p0(pel,ZA,ZB,f,g,u,p,
                                                             u_new,p_new,omega);
                else
                   eprintf("Error: Vanka_smoother not available.\n");
        break;
   case Q_VNVF: if (p_type == Q_SE)
                   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                      full_Vanka_step_p2_p0(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                else if (p_type == Q_SN){
                   if (A_struct & Q_FULL){
                      if (!(B_struct & Q_BDAUGMENT))
                         for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
                            full_Vanka_step_p2_p1_Korn(pnode,ZA,ZB,f,g,u,p,
                                                         u_new,p_new,omega,NO);
/*
                            half_full_Vanka_step_p2_p1_Korn(pnode,ZA,ZB,f,g,u,p,
                                                         u_new,p_new,omega,NO);
*/
                      else{
                         for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
                            full_Vanka_step_p2_p1_korn_BD(pnode,ZA,ZB,f,g,u,p,
                                                             u_new,p_new,omega);
                         for (i = 0; i < j11; i++) {
                            for (pnode=FIRSTNODE(tGrid);pnode;pnode=SUCC(pnode))
                               if (IS_ZNN(pnode))
                                  full_Vanka_step_p2_p1_korn_BD_node(pnode,
                                               ZA,ZB,f,g,u,p,u_new,p_new,omega);
                            for (pface=FIRSTFACE(tGrid);pface;pface=SUCC(pface))
                               if (IS_ZNF(pface))
                                  full_Vanka_step_p2_p1_korn_BD_face(pface,
                                               ZA,ZB,f,g,u,p,u_new,p_new,omega);
                         }
                      }
                   }
                   else if (A_struct & Q_EBDIAG)
                      for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
                         full_Vanka_step_p2_p1(pnode,ZA,ZB,f,g,u,p,
                                                         u_new,p_new,omega,NO);
                   else
                      eprintf("Error: Vanka_smoother not available.\n");
                }
                else if (p_type == Q_SNSSNSSF){
                   if (A_struct & Q_FULL){
                      for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
                         full_Vanka_step_p2_p1_Korn(pnode,ZA,ZB,f,g,u,p,
                                                        u_new,p_new,omega,YES);
/*
                         half_full_Vanka_step_p2_p1_Korn(pnode,ZA,ZB,f,g,u,p,
                                                        u_new,p_new,omega,YES);
*/
                      for (k = 0; k < 10; k++){
                        for (psnode=FIRSTSN(tGrid); psnode; psnode=SUCC(psnode))
                           full_Vanka_step_p2_p2ln(psnode,ZA,ZB,f,g,u,p,
                                                    u_new,p_new,omega,A_struct);
                        for (psface=FIRSTSF(tGrid); psface; psface=SUCC(psface))
                           full_Vanka_step_p2_p2lf(psface,ZA,ZB,f,g,u,p,
                                                    u_new,p_new,omega,A_struct);
                      }
                   }
                   else if (A_struct & Q_EBDIAG){
                      for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
                         full_Vanka_step_p2_p1(pnode,ZA,ZB,f,g,u,p,
                                                        u_new,p_new,omega,YES);
                      for (k = 0; k < 10; k++){
                        for (psnode=FIRSTSN(tGrid); psnode; psnode=SUCC(psnode))
                           full_Vanka_step_p2_p2ln(psnode,ZA,ZB,f,g,u,p,
                                                    u_new,p_new,omega,A_struct);
                        for (psface=FIRSTSF(tGrid); psface; psface=SUCC(psface))
                           full_Vanka_step_p2_p2lf(psface,ZA,ZB,f,g,u,p,
                                                    u_new,p_new,omega,A_struct);
                      }
                   }
                   else
                      eprintf("Error: Vanka_smoother not available.\n");
                }
                else if (p_type == Q_SNE)
                   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                      full_Vanka_step_p2_p1disc(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                else
                   eprintf("Error: Vanka_smoother not available.\n");
        break;
   case Q_VNVFVE: if (p_type == Q_SE)
                   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                      full_Vanka_step_p2b_p0(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                else if (p_type == Q_SNE)
                   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
                      full_Vanka_step_p2b_p1disc(pel,ZA,ZB,f,g,u,p,u_new,p_new,omega);
                else
                   eprintf("Error: Vanka_smoother not available.\n");
        break;
   default:
        eprintf("Error: Vanka_smoother not available.\n");
        break;
   }
   if (t_p & ZERO_MEAN || t_p & WITHOUT_FIRST)
      set_zero_mean(tGrid,p_new,t_p,p_type);
   return(1);
}

void a_smoother(tGrid,Z,h,f,u,d,q,t,type,A_struct,smoother_type)
GRID *tGrid;
INT Z, h, f, u, d, q, t, type, A_struct, smoother_type;
{
   switch(smoother_type){
   case ILU_IT: ILU_step(tGrid,Z,h,f,u,d,q,t,type,A_struct);
        break;
   case SOR_F:  SOR_step_forward(tGrid,Z,u,u,d,t,t,type,type,A_struct,f,
                                 1,2,3,4,5,6,7,8,9,1.,0.,0.);
        break;
   default:
        eprintf("Error: a_smoother not available.\n");
        break;
   }
}

#if (N_DATA & SCALAR_NODE_DATA) && (N_DATA & ONE_NODE_MATR)

void sGauss_Seidel_step(tGrid,Z,f,x,t)
GRID *tGrid;
INT Z, f, x, t;
{
   GRID *theGrid;
   NODE *theNode, *pnode;
   LINK *pl;
   DOUBLE sum;
	
   if (!(t & USE_IS_DIR))
      eprintf("Error: t & USE_IS_DIR should be true\n");
   for (theNode=FIRSTN(tGrid); theNode!= NULL; theNode=SUCC(theNode))
      if (!IS_DIR(theNode,t)){
         sum = 0.0;
         for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
            pnode = NBNODE(pl);
            if (!IS_DIR(pnode,t))
	       sum += COEFFLS(pl,Z)*NDS(pnode,x);
         }
	 NDS(theNode,x) = (NDS(theNode,f) - sum)/COEFFS(theNode,Z);
      }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FIRSTNODE(theGrid); theNode!= NULL; theNode=SUCC(theNode)) 
         if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t)){
            sum = 0.0;
            for (pl=TSTART(theNode); pl!=NULL; pl=NEXT(pl)){
               pnode = ltop_node(NBNODE(pl),tGrid);
               if (!IS_DIR(pnode,t))
	          sum += COEFFLS(pl,Z)*NDS(pnode,x);
            }
            NDS(theNode,x) = (NDS(theNode,f) - sum)/COEFFS(theNode,Z);
         }
   }
}

#else

void sGauss_Seidel_step(tGrid,Z,f,x,t)
GRID *tGrid; INT Z, f, x, t;
{  eprintf("Error: sGauss_Seidel_step not available.\n");  }

#endif
