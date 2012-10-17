/******************************************************************************/
/*                                                                            */
/*                                     CG                                     */
/*                                                                            */
/******************************************************************************/

/* computes the solution p of the system X p = b */
INT CG(tGrid,b,p,g,w,v,imin,imax,r,eps,init,multX,matrix_X,
       q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;
INT b, p, g, w, v;
INT imin, imax, init, matrix_X, q, t, type, structure, 
    p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
FLOAT r, eps, r1, r2, r3;
void (*multX)();
{
   DOUBLE abs_g2, gamma, ro, sum, eps1;
   INT i=0;
  
   if (init & NONZERO_INIT){
      multX(tGrid,matrix_X,p,g,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      subtr(tGrid,g,b,g,t,type);
   }
   else{
      set_value(tGrid,0.,p,t,type);
      inv(tGrid,b,g,t,type);
   }
   copy(tGrid,g,w,t,type);
   abs_g2 = dot(tGrid,g,g,t,type);
   eps1 = MIN(eps,abs_g2*r);
   
   while (i < imin || (i < imax && abs_g2 > eps1)){
/*      printf("%e\n",abs_g2); */
      multX(tGrid,matrix_X,w,v,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      ro = abs_g2/dot(tGrid,v,w,t,type);
      mult_and_add(tGrid,-ro,w,p,p,t,type);
      mult_and_add(tGrid,-ro,v,g,g,t,type);
      sum = dot(tGrid,g,g,t,type);
      gamma = sum/abs_g2;
      abs_g2 = sum;
      mult_and_add(tGrid,gamma,w,g,w,t,type);
      i++;
   }
   return(i);
}

/* computes solution p of Z p = b */
INT PCG(tGrid,b,p,g,w,v,imin,imax,r,eps,init,multX,matrix_X,
        q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3,
        precond_matrix,q0,q1,q2,precond_type)
GRID *tGrid;
INT b, p, g, w, v,
    imin, imax, init, matrix_X, q, t, type, structure,
    p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, 
    precond_matrix, q0, q1, q2, precond_type;
FLOAT r, eps, r1, r2, r3;
void (*multX)();
{
   DOUBLE abs_g2, dot_vg, gamma, ro, sum, eps1;
   INT i=0;
  
   if (init & NONZERO_INIT){
      multX(tGrid,matrix_X,p,g,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      subtr(tGrid,g,b,g,t,type);
   }
   else{
      set_value(tGrid,0.,p,t,type);
      inv(tGrid,b,g,t,type);
   }
   preconditioning(tGrid,precond_matrix,g,w,t,type,q0,q1,q2,precond_type);
   dot_vg = dot(tGrid,w,g,t,type);
   abs_g2=dot(tGrid,g,g,t,type);
   eps1 = MIN(eps,abs_g2*r);
  
   while (i < imin || (i < imax && abs_g2 > eps1)){
/*    printf("%e\n",abs_g2);  */
      multX(tGrid,matrix_X,w,v,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      ro = dot_vg/dot(tGrid,v,w,t,type);
      mult_and_add(tGrid,-ro,w,p,p,t,type);
      mult_and_add(tGrid,-ro,v,g,g,t,type);
      preconditioning(tGrid,precond_matrix,g,v,t,type,q0,q1,q2,precond_type);
      sum = dot(tGrid,v,g,t,type);      
      gamma = sum/dot_vg;
      dot_vg = sum;
      mult_and_add(tGrid,gamma,w,v,w,t,type);
      abs_g2=dot(tGrid,g,g,t,type);
      i++;
   }
   return(i);
}

/******************************************************************************/
/*                                                                            */
/*                                   GMRES                                    */
/*                                                                            */
/******************************************************************************/

#define MACRO1_FOR_GMRES       for (i = 0; i < j; i++){                        \
                                  r[i][j] = c[i]*(d=r[i][j]) - s[i]*r[i+1][j]; \
                                  r[i+1][j] = s[i]*d + c[i]*r[i+1][j];         \
                               }                                               \
                               d = sqrt(r[j][j]*r[j][j] + h[j+1][j]*h[j+1][j]);\
                               c[j] = r[j][j]/d;                               \
                               s[j] = -h[j+1][j]/d;                            \
                               r[j][j] = c[j]*r[j][j] - s[j]*h[j+1][j];        \
                               g[j+1] = s[j]*g[j];                             \
                               g[j] *= c[j];

#define MACRO2_FOR_GMRES(t,type)  h[j+1][j] = sqrt(dot(tGrid,l,l,t,type));     \
                                  mult(tGrid,1./h[j+1][j],l,l,t,type);         \
                                  MACRO1_FOR_GMRES

#define MACRO3_FOR_GMRES(t,type)  d = 0;                                       \
                               for (i = 0; i <= j; i++){                       \
                                  r[i][j] = h[i][j] = dot(tGrid,l,k+i,t,type); \
                                  d += h[i][j]*h[i][j];                        \
                               }                                               \
                               h[j+1][j] = sqrt(dot(tGrid,l,l,t,type)-d);      \
                               MACRO1_FOR_GMRES

#define MACRO4_FOR_GMRES       for (i = j; i >= 0; i--){                       \
                                  for (l = i+1; l <= j; l++)                   \
                                     g[i] -= g[l]*r[i][l];                     \
                                  g[i] /= r[i][i];                             \
                               }

FLOAT GMRES_m(tGrid,u,f,jmin,k,m,eps,jj,multX,Z,
              q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;
INT u, f, Z;      /* estimate of sol., rhs, matrix       */
INT k, m;         /* k,...,k+m+1 auxiliary variables     */
INT jmin, *jj;
INT q, t, type, structure, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
FLOAT eps, r1, r2, r3;
void (*multX)();
{
   INT i, j, l;
   FLOAT h[GMN][GMN], r[GMN][GMN], g[GMN], c[GMN], s[GMN], d;
   
   multX(tGrid,Z,u,k,q,t,t,type,type,structure,
         p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
   subtr(tGrid,f,k,k,t,type);
   g[0] = sqrt(dot(tGrid,k,k,t,type));
   if (g[0] < EPS_GMRES)
      return(g[0]);
   mult(tGrid,1./g[0],k,k,t,type);
   
   for (j = 0; fabs(g[j]) > EPS_GMRES &&
              (j < jmin || (j < m  && fabs(g[j]) > eps)); j++){
      l = k + j + 1;
      multX(tGrid,Z,k+j,l,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      for (i = 0; i <= j; i++)
         r[i][j] = h[i][j] = dot(tGrid,l,k+i,t,type);
      mult_and_add1_for_GMRES(tGrid,h,l,j,k,t,type);
      MACRO2_FOR_GMRES(t,type)
   }
   if (fabs(g[j]) <= EPS_GMRES) *jj += (j--);
   else{
      *jj += j;
      l = k + j + 1;
      multX(tGrid,Z,k+j,l,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      MACRO3_FOR_GMRES(t,type)
   }
   MACRO4_FOR_GMRES
   mult_and_add2_for_GMRES(tGrid,g,u,j,k,t,type);
   return(fabs(g[j+1]));
}

FLOAT PGMRES_m(tGrid,u,f,x,y,jmin,k,m,eps,jj,multX,Z,
               q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3,
               precond_matrix,q0,q1,q2,precond_type)
GRID *tGrid;
INT u, f, Z;  /* estimate of sol., rhs, matrix, ILU  */
INT x, y;     /* auxiliary variables                 */
INT k, m;     /* k,...,k+m+1 auxiliary variables     */
INT jmin, *jj, q, t, type, structure, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9,
    precond_matrix, q0, q1, q2, precond_type;
FLOAT eps, r1, r2, r3;
void (*multX)();
{
   INT i, j, l;
   FLOAT h[GMN][GMN], r[GMN][GMN], g[GMN], c[GMN], s[GMN], d;
   
   multX(tGrid,Z,u,y,q,t,t,type,type,structure,
         p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
   subtr(tGrid,f,y,y,t,type);
   preconditioning(tGrid,precond_matrix,y,k,t,type,q0,q1,q2,precond_type);
   g[0] = sqrt(dot(tGrid,k,k,t,type));
   if (g[0] < EPS_GMRES)
      return(g[0]);
   mult(tGrid,1./g[0],k,k,t,type);
   
   for (j = 0; fabs(g[j]) > EPS_GMRES &&
              (j < jmin || (j < m  && fabs(g[j]) > eps)); j++){
      l = k + j + 1;
      multX(tGrid,Z,k+j,y,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      preconditioning(tGrid,precond_matrix,y,l,t,type,q0,q1,q2,precond_type);
      for (i = 0; i <= j; i++)
         r[i][j] = h[i][j] = dot(tGrid,l,k+i,t,type);
      mult_and_add1_for_GMRES(tGrid,h,l,j,k,t,type);
      MACRO2_FOR_GMRES(t,type)
   }
   if (fabs(g[j]) <= EPS_GMRES) *jj += (j--);
   else{
      *jj += j;
      l = k + j + 1;
      multX(tGrid,Z,k+j,y,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      preconditioning(tGrid,precond_matrix,y,l,t,type,q0,q1,q2,precond_type);
      MACRO3_FOR_GMRES(t,type)
   }
   MACRO4_FOR_GMRES
   mult_and_add2_for_GMRES(tGrid,g,u,j,k,t,type);
   return(fabs(g[j+1]));
}

INT GMRES(tGrid,u,f,jmin,imin,imax,k,m,r,eps,init,multX,Z,
          q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3)
GRID *tGrid;
INT u, f, Z;      /* (estimate of) sol., rhs, matrix     */
INT k, m;         /* k,...,k+m+1 auxiliary variables     */
INT jmin, imin, imax, init;
INT q, t, type, structure, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
FLOAT r, eps, r1, r2, r3;
void (*multX)();
{
   FLOAT eps1, res;
   INT i=1, j=0;
   
 if (jmin > m){ 
    eprintf("Error: wrong relation between jmin and m in GMRES.\n");
    printf("       (jmin=%i but m=%i)\n",jmin,m);
 }
 else if (GMN < m+2){
    eprintf("Error: wrong relation between GMN and m in GMRES.\n");
    printf("       (GMN=%i but m=%i)\n",GMN,m);
 }
 else{
   if (init & NONZERO_INIT){
      multX(tGrid,Z,u,k,q,t,t,type,type,structure,
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      subtr(tGrid,f,k,k,t,type);
   }
   else{
      set_value(tGrid,0.,u,t,type);
      copy(tGrid,f,k,t,type);
   }
   eps1 = MIN(eps,sqrt(dot(tGrid,k,k,t,type))*r);
   if (eps1 < EPST){
      imin = 2;
      eps1 = EPST;
   }
   res = GMRES_m(tGrid,u,f,jmin,k,m,eps1,&j,multX,Z,
                 q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
   
   while (i < imin || (i < imax && res > eps1)){
      printf("%e\n",res);
      res = GMRES_m(tGrid,u,f,jmin,k,m,eps1,&j,multX,Z,
                    q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3);
      i++;
   }
 /*  printf("   GMRES:  i = %i,   res = %lf\n",i,res);*/
 }
 return(j);
}

INT PGMRES(tGrid,u,f,x,y,jmin,imin,imax,k,m,r,eps,init,multX,Z,
           q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3,
           precond_matrix,q0,q1,q2,precond_type)
GRID *tGrid;
INT u, f, Z;  /* (estimate of) sol., rhs, matrix, ILU  */
INT x, y;     /* auxiliary variables                   */
INT k, m;     /* k,...,k+m+1 auxiliary variables       */
INT jmin, imin, imax, init,
    q, t, type, structure, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9,
    precond_matrix, q0, q1, q2, precond_type;
FLOAT r, eps, r1, r2, r3;
void (*multX)();
{
   FLOAT eps1, res;
   INT i=1, j=0;
   
 if (jmin > m){ 
    eprintf("Error: wrong relation between jmin and m in PGMRES.\n");
    printf("       (jmin=%i but m=%i)\n",jmin,m);
 }
 else if (GMN < m+2){
    eprintf("Error: wrong relation between GMN and m in PGMRES.\n");
    printf("       (GMN=%i but m=%i)\n",GMN,m);
 }
 else{
   if (!(init & NONZERO_INIT)){
      set_value(tGrid,0.,u,t,type);
      eps1 = MIN(eps,sqrt(dot(tGrid,f,f,t,type))*r);
   }
   else{
      defect(tGrid,Z,f,u,x,x,t,type,structure);
      eps1 = MIN(eps,sqrt(dot(tGrid,x,x,t,type))*r);
   }
   if (eps1 < EPST){
      imin = 2;
      eps1 = EPST;
   }
   res = PGMRES_m(tGrid,u,f,x,y,jmin,k,m,eps1,&j,multX,Z,
                  q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3,
                  precond_matrix,q0,q1,q2,precond_type);
   
   while (i < imin || (i < imax && res > eps1)){
      printf("%e\n",res);
      res = PGMRES_m(tGrid,u,f,x,y,jmin,k,m,eps1,&j,multX,Z,
                     q,t,type,structure,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,r1,r2,r3,
                     precond_matrix,q0,q1,q2,precond_type);
      i++;
   }
 /*  printf("   GMRES:  i = %i,   res = %lf\n",i,res);*/
 }
 return(j);
}

INT PGMRES_Ar(tGrid,Z,x,f,xe,y,jmin,imin,imax,init,k,m,r,r1,eps,d)
GRID *tGrid;
INT d;           /* d^2 ... diagonal preconditioner       */
INT x, f, Z;     /* (estimate of) sol., rhs, matrix       */
INT xe, y;       /* auxiliary variables                   */
INT k, m;        /* k,...,k+m+1 aux. node-face variables  */
INT jmin, imin, imax,init;
FLOAT r, r1, eps;
{   
   FLOAT eps1, res;
   INT i=1, j=0;
   
 if (jmin > m){ 
    eprintf("Error: wrong relation between jmin and m in PGMRES_Ar.\n");
    printf("       (jmin=%i but m=%i)\n",jmin,m);
 }
 else if (GMN < m+2){
    eprintf("Error: wrong relation between GMN and m in PGMRES_Ar.\n");
    printf("       (GMN=%i but m=%i)\n",GMN,m);
 }
 else{
   if (!(init & NONZERO_INIT))
      vs_set_value_nf(tGrid,0.,x,1);
   mult_Ar(tGrid,Z,x,y,k,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
           Z,B_STRUCT,k,xe,T_FOR_P,P_TYPE,6,7,8,9,r,0.,0.);
   vs_divide_nf(tGrid,y,d,y,1);
   vs_multiply_nf(tGrid,x,d,x,1);
   vs_divide_nf(tGrid,f,d,f,1);
   vs_subtr_nf(tGrid,f,y,y,1);
   eps1 = MIN(eps,sqrt(vs_dot_nf(tGrid,y,y,1))*r1);
   if (eps1 < EPST){
      imin = 2;
      eps1 = EPST;
   }
   res = GMRES_m(tGrid,x,f,jmin,k,m,eps1,&j,Pmult_Ar,Z,
                 y,T_FOR_U,U_TYPE,A_STRUCT,Z,B_STRUCT,y,xe,
                 T_FOR_P,P_TYPE,d,7,8,9,r,0.,0.);
   
   while (i < imin || (i < imax && res > eps1)){
      printf("%e\n",res);
      res = GMRES_m(tGrid,x,f,jmin,k,m,eps1,&j,Pmult_Ar,Z,
                    y,T_FOR_U,U_TYPE,A_STRUCT,Z,B_STRUCT,y,xe,
                    T_FOR_P,P_TYPE,d,7,8,9,r,0.,0.);
      i++;
   }
 /* printf("   GMRES:  i = %i,   res = %lf\n",i,res); */
   vs_divide_nf(tGrid,x,d,x,1);
 }
 return(j);
}

DOUBLE fcn_for_minim(tGrid,x,t,type)
GRID *tGrid;
INT x, t, type;
{
   mult_A(tGrid,A,x,D,F,t,t,type,type,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   return(0.5*dot(tGrid,D,x,t,type)-dot(tGrid,F,x,t,type));
}

void grad_fcn_for_minim(tGrid,x,y,t,type)
GRID *tGrid;
INT x, y, t, type;
{
   defect(tGrid,A,F,x,y,F,t,type,A_STRUCT);
}

DOUBLE line_phi(tGrid,a,x,p,z,t,type)
GRID *tGrid;
DOUBLE a;
INT x, p, z, t, type;
{
   mult_and_add(tGrid,a,p,x,z,t,type);
   return(fcn_for_minim(tGrid,z,t,type));
}

DOUBLE der_line_phi(tGrid,a,x,p,y,z,t,type)
GRID *tGrid;
DOUBLE a;
INT x, p, y, z, t, type;
{
   mult_and_add(tGrid,a,p,x,z,t,type);
   grad_fcn_for_minim(tGrid,z,y,t,type);
   return(dot(tGrid,y,p,t,type));
}

DOUBLE zoom(tGrid,zoom_it,a1,a2,phi1,phi2,dphi1,phi0,dphi0,phi3,dphi3,
            c1,c2,x,p,y,z,t,type)
GRID *tGrid;
DOUBLE a1, a2, phi1, phi2, dphi1, phi0, dphi0, *phi3, *dphi3, c1, c2;
INT zoom_it, x, p, y, z, t, type;
{
   DOUBLE c, a3;
   INT i=0, k=1;

   while (k && i < zoom_it){
//    a3 = argmin q(x) where q(x) = phi1 + dphi1*(x-a1) + c*(x-a1)*(x-a1) 
//    with q(a2)=phi2 
      c = (phi2 - phi1 - dphi1*(a2-a1))/((a2-a1)*(a2-a1));
      a3 = a1 - 0.5*dphi1/c;
      if (!( (a1 <= a3 && a3 <= a2) || (a2 <= a3 && a3 <= a1) ))
         a3 = 0.5*(a1 + a2);
      *phi3 =  line_phi(tGrid,a3,x,p,z,t,type);
      if (*phi3 > phi0 + c1*a3*dphi0 || *phi3 >= phi1){
         a2 = a3;
         phi2 = *phi3;
      }
      else{
         *dphi3 = der_line_phi(tGrid,a3,x,p,y,z,t,type);
         if (fabs(*dphi3) <= -c2*dphi0)
            k = 0;
         if ((*dphi3)*(a2-a1) >= 0.){
            a2 = a1;
            phi2 = phi1;
         }
         a1 = a3;
         phi1 = *phi3;
         dphi1 = *dphi3;
      }
      i++;
   }
// printf("%i iterations in zoom\n",i);
   return(a3);
}

INT line_search(tGrid,line_it,zoom_it,phi0,dphi0,phi0_new,dphi0_old,a_old,
                   c1,c2,x,p,y,z,t,type)
GRID *tGrid;
DOUBLE phi0, dphi0, *phi0_new, dphi0_old, *a_old, c1, c2;
INT line_it, zoom_it, x, p, y, z, t, type;
{
   DOUBLE a1=0., a2, phi1=phi0, phi2, dphi1=dphi0, dphi2, phi3, dphi3;
   INT i=0, k=1;

   if (dphi0 >= 0.)
      eprintf("Error in line_search: not a descent direction.\n");
   a2 = (*a_old)*dphi0_old/dphi0;
   while (k && i < line_it){
      phi2 = line_phi(tGrid,a2,x,p,z,t,type);
      if (phi2 > phi0 + c1*a2*dphi0 || (i && phi2 >= phi1)){
         a2 = zoom(tGrid,zoom_it,a1,a2,phi1,phi2,dphi1,phi0,dphi0,
                   &phi2,&dphi2,c1,c2,x,p,y,z,t,type);
         k = 0;
      }
      else{
         dphi2 = der_line_phi(tGrid,a2,x,p,y,z,t,type);
         if (fabs(dphi2) <= -c2*dphi0)
            k = 0;
         else if (dphi2 >= 0){
               a2 = zoom(tGrid,zoom_it,a2,a1,phi2,phi1,dphi2,phi0,dphi0,
                         &phi2,&dphi2,c1,c2,x,p,y,z,t,type);
               k = 0;
         }
         else{
            phi1 = phi2;
            dphi1 = dphi2;
            a1 = a2;
            a2 *= 2.;
            i++;
         }
      }
   }
// printf("%i iterations in line_search\n",i);
   if (phi2 <= phi0 + c1*a2*dphi0 && fabs(dphi2) <= -c2*dphi0){
      *a_old = a2;
      *phi0_new = phi2;
      return(1);
   }
   else{
      eprintf("Error in line_search: parameter not found.\n");
      return(0);
   }
}

void nonlinear_CG(tGrid,max_it,line_it,zoom_it,eps,x,p,y,z,v,w,t,type,
                  fcn,grad_fcn)
GRID *tGrid;
INT max_it, line_it, zoom_it, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE a=1., c1=1.e-4, c2=0.1, q, r, sum, phi0, dphi0, dphi0_old;
   INT i=0, k1=1, k2=1;

   set_value(tGrid,0.,x,t,type);
   phi0 = fcn(tGrid,x,t,type);
   grad_fcn(tGrid,x,v,t,type);
   copy(tGrid,v,w,t,type);
   inv(tGrid,v,p,t,type);
   r = dot(tGrid,v,v,t,type);
   dphi0_old = -r;
   while (k2 && 
          max_abs_value(tGrid,v,t,type) > eps*(1.+fabs(phi0)) && i < max_it){
      dphi0 = dot(tGrid,v,p,t,type);
      if (line_search(tGrid,line_it,zoom_it,phi0,dphi0,&phi0,dphi0_old,&a,c1,c2,
                      x,p,y,z,t,type)){
         mult_and_add(tGrid,a,p,x,x,t,type);
         grad_fcn(tGrid,x,v,t,type);
         sum = dot(tGrid,v,v,t,type);
         q = dot(tGrid,v,w,t,type);
         if (fabs(q) > 0.1*r || q > sum)
            inv(tGrid,v,p,t,type);
         else
            mult_and_subtr(tGrid,(sum-q)/r,p,v,p,t,type);
         r = sum;
         dphi0_old = dphi0;
         copy(tGrid,v,w,t,type);
         k1 = 1;
      }
      else if (k1){
         k1 = 0;
         inv(tGrid,v,p,t,type);
         a = 1.;
      }
      else
         k2 = 0;
      i++;
   }
   printf("%i iterations in nonlinear CG.\n",i);
}
