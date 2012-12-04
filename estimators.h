/******************************************************************************/
/*                                                                            */
/*                              error estimators                              */
/*                                                                            */
/******************************************************************************/

DOUBLE dot();
void coord_of_barycentre();
DOUBLE L1_norm_of_p1_fcn_on_triangle();

 void estimate1(mg)
 MULTIGRID *mg;
 {
    ELEMENT *pelem;
    GRID *theGrid;
    
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          if (IS_TOP_ELEMENT(pelem)) pelem->eflag = 2;
 }    
 
 void estimate2(mg)
 MULTIGRID *mg;
 {
    ELEMENT *pelem;
    GRID *theGrid;
    INT i=0;
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          if (IS_TOP_ELEMENT(pelem)){
           if(i++<10) pelem->eflag = 2;
           else pelem->eflag = 0;
          }
 } 
 
INT is_in_element();

void mark_element_with_point(tGrid,x,y)
GRID *tGrid;
FLOAT x, y;
{
   ELEMENT *pel;
   FLOAT x0[2]={x,y}, x_ref[2], eps=0.;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel && i; pel = pel->succ)
      if (is_in_element(x0,pel,eps,x_ref)){
         pel->eflag = 2;
         i = 0;
      }
   if (i)
      eprintf("Error in mark_element_with_point.\n"); 
}

#if E_DATA & SCALAR_ELEMENT_DATA

 /*  err ... error indicator on one element; err(E1)+err(E2)=err(E1+E2);
     q   ... relative part of elements to refine primarily (q\in(0,1)).       */
 void estimate(mg,q,e,f,err)
 MULTIGRID *mg;
 FLOAT q, (*err)();
 INT e,f;
 {
    GRID *theGrid;
    ELEMENT *pelem, *pel;
    FLOAT c, min=1.0e18, max=0.0, n;
    INT *a, i, ncat=1000, ncat1=1001, nelem=0, ntopel=0, nrefel=0, sum;
    
    a = (INT *) calloc(ncat1,sizeof(INT));
    
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          ED(pelem,f) = 0.0;
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          if (IS_TOP_ELEMENT(pelem)){
             if (pelem->status == 0){
                ED(pelem,e) = (*err)(TOP_GRID(mg),pelem,U);
                nelem++;
             }
             else if (ED(pelem,f) < 0.5){
                c = 0.0;
                i = 0;
                while (pel=pelem->father->sons[i++]) 
                   c += (*err)(TOP_GRID(mg),pel,U);
                i = 0;
                while (pel=pelem->father->sons[i++]){
                   ED(pel,e) = c;
                   ED(pel,f) = 1.0;
                }
                nelem++;
             }
             if (ED(pelem,e) > max) max = ED(pelem,e);
             else if (ED(pelem,e) < min) min = ED(pelem,e);
             ntopel++;
          }   
    c = ncat/(max - min);
    for (i = 0; i < ncat1; i++) a[i] = 0;
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          if ( IS_TOP_ELEMENT(pelem) && 
                                 (pelem->status == 0 || ED(pelem,f) > 0.5) ){
             a[ (int)(floor( (ED(pelem,e)-min)*c )) ]++;
             i = 0;
             if (pelem->status)
                while (pel=pelem->father->sons[i++])
                   ED(pel,f) = 0.0;
          }
    n = q*nelem;
    sum = a[ncat];
    i = ncat;
    while (sum < n) sum += a[--i];
    if ( sum - n > n - sum + a[i] ){
       sum -= a[i];
       i++;
    }
    c = i/c + min;
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
          if (IS_TOP_ELEMENT(pelem))
             if ( ED(pelem,e) >= c){
                pelem->eflag = 2;
                nrefel++;
             }
             else
                pelem->eflag = 0;
    printf("%4.2f %% of top (macro)elements ",100.0*sum/nelem);
    printf("have to be refined primarily (green elements building up a ");
    printf("forefather element are regarded as one element).\n");
    printf("Primary increase of top element number: %i elements (%5.1f %% ",
                                      8*sum-nrefel,100.0*(8*sum-nrefel)/ntopel);
    printf("of all top elements).\n");
    printf("Min. estimated error = %e, max estimated error = %e, ",min,max);
    printf("(max-min)/min = %5.1f;\n",(max-min)/min);
    printf("min. error for marking to refine: ref_err = %e ",c);
    printf("( (ref_err-min)/min = %5.1f)\n",(c-min)/min);
 } 

#else  /*  !(E_DATA & SCALAR_ELEMENT_DATA)  */

 void estimate(mg,q,e,f,err)
 MULTIGRID *mg; FLOAT q, (*err)(); INT e,f;
 {  eprintf("Error: estimate not available.\n");  }

#endif

#if DIM == 2

void mark_elements(tGrid)
GRID *tGrid;
{
   ELEMENT *pelem;
   FLOAT *x1, *x2, *x3;
   
   for (pelem = FIRSTELEMENT(tGrid);pelem != NULL;pelem = pelem->succ){
      VERTICES_OF_ELEMENT(x1,x2,x3,pelem);
      if ( fabs(x1[0]-1.)<EPSA || fabs(x2[0]-1.)<EPSA || fabs(x3[0]-1.)<EPSA ||
           fabs(x1[1]-1.)<EPSA || fabs(x2[1]-1.)<EPSA || fabs(x3[1]-1.)<EPSA)
         pelem->eflag = 2;
      else
         pelem->eflag = 0;
   }
}

#else

void mark_elements(tGrid)
GRID *tGrid;
{  eprintf("Error: mark_elements not available.\n");  }

#endif

#if PAR_OPT == YES

#if N_DATA & SCALAR_NODE_DATA

void mark_outflow_boundary_nodes(tGrid,bb0,bb1,q)
GRID *tGrid;
DOUBLE (*bb0)(), (*bb1)();
INT q;  /*  q ... scalar node variable  */
{
   ELEMENT *pelem;
   DOUBLE a[DIM], b[DIM], nn[DIM], bb[DIM], xjk[DIM], r;
   INT i, j, k;
   
   set_value(tGrid,0.,q,0,Q_SN); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      for (i = 0; i < SIDES; i++)
         if (NOT_FF(pelem->f[i])){
            j = (i+1)%SIDES;
            k = (i+2)%SIDES;
            SUBTR(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,a);
            SUBTR(pelem->n[j]->myvertex->x,pelem->n[i]->myvertex->x,b);
            ORT_VECT_DIR(nn,a,b)
            AVERAGE(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,xjk)
            bb[0] = bb0(xjk);
            bb[1] = bb1(xjk);
            if (fabs(bb[0]) + fabs(bb[1]) > 1.e-30)
               if ((r=DOT(bb,nn)/sqrt(DOT(bb,bb)*DOT(nn,nn))) > -0.01){
                  if (r < 0.5)  /* -> characteristic layer */
                     NDS(pelem->n[j],q) = NDS(pelem->n[k],q) = 20.;
                  else{
                     if (NDS(pelem->n[j],q) < 5.)
                        NDS(pelem->n[j],q) = 10.;
                     if (NDS(pelem->n[k],q) < 5.)
                        NDS(pelem->n[k],q) = 10.;
                  }
               }
         }
}

void mark_elements_for_error_ind(tGrid,bb0,bb1,q)
GRID *tGrid;
DOUBLE (*bb0)(), (*bb1)();
INT q;  /*  q ... scalar node variable  */
{
   ELEMENT *pelem;
   INT i;
   
   if (FCN_AT_INFLOW == YES){
      mark_outflow_boundary_nodes(tGrid,bb0,bb1,q);
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
         if (NDS(pelem->n[0],q) > 5. || NDS(pelem->n[1],q) > 5. 
                                     || NDS(pelem->n[2],q) > 5.){
            pelem->eflag = 1;
            if (NDS(pelem->n[0],q) > 15. || NDS(pelem->n[1],q) > 15. 
                                         || NDS(pelem->n[2],q) > 15.){
               for (i = 0; i < NVERT; i++)
                  if (NDS(pelem->n[i],q) < 6.)
                     NDS(pelem->n[i],q) = 4.;
            }
            else
               for (i = 0; i < NVERT; i++)
                  if (NDS(pelem->n[i],q) < 3.)
                     NDS(pelem->n[i],q) = 2.;
         }
         else
            pelem->eflag = 0;
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
         if (pelem->eflag == 0 && (NDS(pelem->n[0],q) > 1. 
              || NDS(pelem->n[1],q) > 1. || NDS(pelem->n[2],q) > 1.))
            if (NDS(pelem->n[0],q) > 3. || NDS(pelem->n[1],q) > 3. 
                                        || NDS(pelem->n[2],q) > 3.)
               pelem->eflag = 3;
            else
               pelem->eflag = 2;
   }
   else
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
         if (!(IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                                  && IS_FN(pelem->n[2])))
            pelem->eflag = 1;
         else
            pelem->eflag = 0;
}

void set_zero_along_outflow_boundary(tGrid,u)
GRID *tGrid;
INT u;  /*  u ... scalar node variable  */
{
   ELEMENT *pelem;
   NODE *pn;

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      if (IS_BN(pn))
         NDS(pn,u) = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag == 1)
         NDS(pelem->n[0],u) = NDS(pelem->n[1],u) = NDS(pelem->n[2],u) = 0.;
}

#else

void mark_elements_for_error_ind(tGrid,bb0,bb1,q)
GRID *tGrid; DOUBLE (*bb0)(), (*bb1)(); INT q;
{  eprintf("Error: mark_elements_for_error_ind not available.\n");  }

void set_zero_along_outflow_boundary(tGrid,u)
GRID *tGrid; INT u;  /*  u ... scalar node variable  */
{  eprintf("Error: set_zero_along_outflow_boundary not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

DOUBLE langevin(x)
DOUBLE x;
{
   if (x < 0.01)
      return(x/3.);
   else
      return(1./tanh(x)-1./x);
}

DOUBLE der_langevin(x)
DOUBLE x;
{
   double a;

   if (x < 0.01)
      return(1./3.);
   else{
      a = sinh(x);
      a *= a;
      return(-1./a+1./(x*x));
   }
}

DOUBLE fcn0(x)
DOUBLE x;
{
// return(x);
/*
   if (x > 1.)
      return(x);
   else
      return(2.*x*x-x*x*x);
*/
   if (x > 1.)
      return(sqrt(x));
   else
      return(2.5*x*x-1.5*x*x*x);
}

DOUBLE der_fcn0(x)
DOUBLE x;
{
// return(1.);
/*
   if (x > 1.)
      return(1.);
   else
      return(4.*x-3.*x*x);
*/
   if (x > 1.)
      return(0.5/sqrt(x));
   else
      return(5.*x-4.5*x*x);
}

DOUBLE fcn1(x)
DOUBLE x;
{
   if (x > 1.)
      return(2.);
   else
//    return(4.*x-2.*x*x);
      return(x*x*x*x-2.*x*x*x-x*x+4.*x);
//    return(2.*x*x*x*x-4.*x*x*x+4.*x);
//    return(-4.*x*x*x*x+8.*x*x*x-6.*x*x+4.*x);
}

DOUBLE der_fcn1(x)
DOUBLE x;
{
   if (x > 1.)
      return(0.);
   else
//    return(4.-4.*x);
      return(4.*x*x*x-6.*x*x-2.*x+4.);
//    return(8.*x*x*x-12.*x*x+4.);
//    return(-16.*x*x*x+24.*x*x-12.*x+4.);
}

DOUBLE fcn_in_ind(x)
DOUBLE x;
{
   double a=A_IN_FCN;

   return(sqrt(a)*fcn0(x/a));
}

DOUBLE der_of_fcn_in_ind(x)
DOUBLE x;
{
   double a=A_IN_FCN;

   return(der_fcn0(x/a)/sqrt(a));
}

DOUBLE fcn_in_res(x)
DOUBLE x;
{
   double a=RES_SCALING*RES_SCALING;

   return(fcn1(x*x/a));
}

DOUBLE der_of_fcn_in_res(x)
DOUBLE x;
{
   double a=RES_SCALING*RES_SCALING;

   return(der_fcn1(x*x/a)*2.*x/a);
}

#if E_DATA & VECTOR_ELEMENT_DATA

void compute_res(tGrid,u,v,eps,bb0,bb1,react,rhs)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u, v;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], rdetB;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag != 1){
         rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                         pelem->n[1]->myvertex->x,
                                         pelem->n[2]->myvertex->x,bar);
         coord_of_barycentre(pelem,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
         EDV(pelem,v,0) = fabs(DOT(bb,grad) + 
               react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                            +NDS(pelem->n[2],u))/3. - rhs(xc));
         EDV(pelem,v,1) = rdetB;
      }
      else
         EDV(pelem,v,0) = EDV(pelem,v,1) = 0.;
}

void compute_norm_of_res(tGrid,v,norm,p_in_norm)
GRID *tGrid;
DOUBLE p_in_norm;
INT v, norm;
{
   ELEMENT *pelem;
   FLOAT sum=0.;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (norm == L1_NORM)
         sum += EDV(pelem,v,0)*EDV(pelem,v,1);
      else if (norm == L2_NORM)
         sum += EDV(pelem,v,0)*EDV(pelem,v,0)*EDV(pelem,v,1);
      else if (norm == LP_NORM)
         sum += pow(EDV(pelem,v,0),p_in_norm)*EDV(pelem,v,1);
      else if (norm == RES_FCN)
         sum += fcn_in_res(EDV(pelem,v,0))*EDV(pelem,v,1);
   if (norm == L1_NORM)
      printf("L1 norm of res: %e\n",sum);
   else if (norm == L2_NORM)
      printf("L2 norm of res: %e\n",sqrt(sum));
   else if (norm == LP_NORM)
      printf("Lp norm of res with p=%2.0f: %e\n",p_in_norm,
                                                 pow(sum,1./p_in_norm));
   else if (norm == RES_FCN)
      printf("L1 norm of fcn(res): %e\n",sum);
}

void compute_mid_res(tGrid,v)
GRID *tGrid;
INT v;
{
   ELEMENT *pelem;
   FLOAT sum=0., ar=0., mid;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag != 1){
         sum += EDV(pelem,v,0)*EDV(pelem,v,1);
         ar  += EDV(pelem,v,1);
      }
   mid = sum/ar;
   printf("mid = %e\n",mid);
   sum = ar = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag != 1 && EDV(pelem,v,0) < mid){
         sum += EDV(pelem,v,0)*EDV(pelem,v,1);
         ar  += EDV(pelem,v,1);
      }
   mid = sum/ar;
   printf("mid = %e\n",mid);
}

double max_crosswind_derivative(tGrid,u,bb0,bb1)
GRID *tGrid;
FLOAT (*bb0)(), (*bb1)();
INT u;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bbo[DIM], q, max=0.;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag != 1){
         barycentric_coordinates(pelem->n[0]->myvertex->x,
                                 pelem->n[1]->myvertex->x,
                                 pelem->n[2]->myvertex->x,bar);
         coord_of_barycentre(pelem,xc);
         bbo[0] = -bb1(xc);
         bbo[1] =  bb0(xc);
         q = sqrt(DOT(bbo,bbo));
         if (q > 1.e-30){
            SET22(bbo,bbo,q)
            sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
            q = fabs(DOT(bbo,grad));
            max = MAX(max,q);
         }
      }
   if (max < 1.e-30)
      eprintf("Error: too small crosswind derivative.\n");
   return(max);
}

double compute_crosswind_derivative(tGrid,u,v,bb0,bb1)
GRID *tGrid;
FLOAT (*bb0)(), (*bb1)();
INT u, v;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bbo[DIM], q;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if (pelem->eflag != 1){
         barycentric_coordinates(pelem->n[0]->myvertex->x,
                                 pelem->n[1]->myvertex->x,
                                 pelem->n[2]->myvertex->x,bar);
         coord_of_barycentre(pelem,xc);
         bbo[0] = -bb1(xc);
         bbo[1] =  bb0(xc);
         q = sqrt(DOT(bbo,bbo));
         if (q > 1.e-30){
            SET22(bbo,bbo,q)
            sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
            EDV(pelem,v,0) = DOT(bbo,grad);
         }
         else
            EDV(pelem,v,0) = 0.;
      }
      else
         EDV(pelem,v,0) = 0.;
}

#else

void compute_res(tGrid,u,v,eps,bb0,bb1,react,rhs)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(); INT u, v;
{  eprintf("Error: compute_res not available.\n"); return(0.);  }

void compute_norm_of_res(tGrid,v,norm,p_in_norm)
GRID *tGrid; DOUBLE p_in_norm; INT v, norm;
{  eprintf("Error: compute_norm_of_res not available.\n"); return(0.);  }

void compute_mid_res(tGrid,v)
GRID *tGrid; INT v;
{  eprintf("Error: compute_mid_res not available.\n"); return(0.);  }

double max_crosswind_derivative(tGrid,u,bb0,bb1)
GRID *tGrid; FLOAT (*bb0)(), (*bb1)(); INT u;
{  eprintf("Error: max_crosswind_derivative not available.\n"); return(0.);  }

double compute_crosswind_derivative(tGrid,u,v,bb0,bb1)
GRID *tGrid; FLOAT (*bb0)(), (*bb1)(); INT u, v;
{  eprintf("Error: compute_crosswind_derivative not available.\n"); return(0.);  }

#endif

double res_and_der_based_error_indicator(tGrid,u,eps,bb0,bb1,react,rhs)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], bbo[DIM], 
         rdetB, q, res, err=0.;
FLOAT *x0, *x1, *x2;
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if ((RESTRICT_FCN == NO && pelem->eflag != 1) ||
          (RESTRICT_FCN == YES && (pelem->eflag == 2 || pelem->eflag == 3))){
         rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                         pelem->n[1]->myvertex->x,
                                         pelem->n[2]->myvertex->x,bar);
         coord_of_barycentre(pelem,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         ORT_VECT(bbo,bb)
         q = sqrt(DOT(bbo,bbo));
         if (q > 1.e-30)
            SET22(bbo,bbo,q)
         sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
         if (RES_NORM_IN_FCN == L1_NORM || RES_NORM_IN_FCN == L2_NORM ||
             RES_NORM_IN_FCN == LP_NORM || RES_NORM_IN_FCN == RES_FCN){
           if (RES_NORM_IN_FCN == L1_NORM){
           VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
           res = L1_norm_of_p1_fcn_on_triangle(x0,x1,x2,
             bb0(x0)*grad[0]+bb1(x0)*grad[1]+react(x0)*NDS(pelem->n[0],u)-rhs(x0),
             bb0(x1)*grad[0]+bb1(x1)*grad[1]+react(x1)*NDS(pelem->n[1],u)-rhs(x1),
             bb0(x2)*grad[0]+bb1(x2)*grad[1]+react(x2)*NDS(pelem->n[2],u)-rhs(x2));
           }
           else if (RES_NORM_IN_FCN == L2_NORM){
//           res = L2norm2_of_cd_res(pelem,bar,u,P1C,eps,bb0,bb1,react,rhs);
             res = DOT(bb,grad) + 
                react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                             +NDS(pelem->n[2],u))/3. - rhs(xc);
             res = res*res*rdetB;
           }
           else if (RES_NORM_IN_FCN == LP_NORM){
             res = DOT(bb,grad) + 
               react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                            +NDS(pelem->n[2],u))/3. - rhs(xc);
             res = pow(fabs(res),P_IN_LP_NORM)*rdetB;
           }
           else if (RES_NORM_IN_FCN == RES_FCN){
             res = DOT(bb,grad) + 
               react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                            +NDS(pelem->n[2],u))/3. - rhs(xc);
             res = fcn_in_res(fabs(res))*rdetB;
           }
           res *= RES_NORM_WEIGHT;

#if E_DATA & VECTOR_ELEMENT_DATA

           if (pelem->eflag != 1){
             EDV(pelem,ERR_VAR,0) = res;
             EDV(pelem,ERR_VAR,1) = EDV(pelem,CROSS_FCN_W,1)*CROSS_DER_WEIGHT*
               fcn_in_ind(fabs(DOT(bbo,grad))/MAX_CROSS_DER)*rdetB;
           }
           else{
             EDV(pelem,ERR_VAR,0) = res;
             EDV(pelem,ERR_VAR,1) = 0.;
           }
           err += EDV(pelem,ERR_VAR,0) + EDV(pelem,ERR_VAR,1);

#else

           err += res + fcn_in_ind(fabs(DOT(bbo,grad)))*rdetB;

#endif

        }
        else if (RES_NORM_IN_FCN == RES_X_CROSS){
          res = DOT(bb,grad) + 
               react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                            +NDS(pelem->n[2],u))/3. - rhs(xc);
//        if (pelem->eflag == 2 || pelem->eflag == 3)
//          res = res*res*rdetB;
//        else
          err += fcn_in_res(fabs(res))
                *fcn_in_ind(fabs(DOT(bbo,grad))/MAX_CROSS_DER)*rdetB;
        }
      }
   return(err);
}

double derivatives_of_res_and_der_based_error_indicator(tGrid,u,d,
                                                          eps,bb0,bb1,react,rhs)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u, d;  /*  u ... solution at nodes, d ... derivatives at nodes  */
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], bbo[DIM],
         rdetB, res, p, q, err=0.;
   INT i;

   set_value(tGrid,0.,d,0,U_TYPE); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if ((RESTRICT_FCN == NO && pelem->eflag != 1) ||
          (RESTRICT_FCN == YES && (pelem->eflag == 2 || pelem->eflag == 3))){
         rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                         pelem->n[1]->myvertex->x,
                                         pelem->n[2]->myvertex->x,bar);
         coord_of_barycentre(pelem,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         ORT_VECT(bbo,bb)
         q = sqrt(DOT(bbo,bbo));
         if (q > 1.e-30)
            SET22(bbo,bbo,q)
         sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
         q = react(xc)/3.;
         res = DOT(bb,grad)
               + q*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u) +NDS(pelem->n[2],u)) 
               - rhs(xc);
      if (RES_NORM_IN_FCN == L1_NORM || RES_NORM_IN_FCN == L2_NORM ||
          RES_NORM_IN_FCN == LP_NORM || RES_NORM_IN_FCN == RES_FCN){
         if (RES_NORM_IN_FCN == L1_NORM){
            if (fabs(res) > 1.e-30){
               res = SGN(res);
               for (i = 0; i < NVERT; i++)
                  NDS(pelem->n[i],d) +=
                                   res*(DOT(bb,bar[i])+q)*rdetB*RES_NORM_WEIGHT;
            }
         }
         else if (RES_NORM_IN_FCN == L2_NORM){
            for (i = 0; i < NVERT; i++)
               NDS(pelem->n[i],d) +=
//     2.*cd_res_times_test_cd_res(pelem,bar,u,Q,i,eps,bb0,bb1,react,rhs)*RES_NORM_WEIGHT;
             2.*res*(DOT(bb,bar[i])+q)*rdetB*RES_NORM_WEIGHT;
         }
         else if (RES_NORM_IN_FCN == LP_NORM){
            for (i = 0; i < NVERT; i++)
               NDS(pelem->n[i],d) +=
                  P_IN_LP_NORM*pow(fabs(res),P_IN_LP_NORM-1.)*SGN(res)
                              *(DOT(bb,bar[i])+q)*rdetB*RES_NORM_WEIGHT;
         }
         else if (RES_NORM_IN_FCN == RES_FCN){
            for (i = 0; i < NVERT; i++)
               NDS(pelem->n[i],d) +=
                  der_of_fcn_in_res(fabs(res))*SGN(res)
                              *(DOT(bb,bar[i])+q)*rdetB*RES_NORM_WEIGHT;
         }
         if (pelem->eflag != 1){
            q = DOT(bbo,grad)/MAX_CROSS_DER;
            if (fabs(q) > 1.e-30){

#if E_DATA & VECTOR_ELEMENT_DATA

               q = EDV(pelem,CROSS_FCN_W,1)*CROSS_DER_WEIGHT*
                   SGN(q)*der_of_fcn_in_ind(fabs(q))/MAX_CROSS_DER;

#else

               q = CROSS_DER_WEIGHT*
                   SGN(q)*der_of_fcn_in_ind(fabs(q))/MAX_CROSS_DER;

#endif

               for (i = 0; i < NVERT; i++)
                  NDS(pelem->n[i],d) += q*DOT(bbo,bar[i])*rdetB;
            }
         }
      }
      else if (RES_NORM_IN_FCN == RES_X_CROSS){
/*
         if (pelem->eflag == 2 || pelem->eflag == 3)
            for (i = 0; i < NVERT; i++)
               NDS(pelem->n[i],d) += 2.*res*(DOT(bb,bar[i])+q)*rdetB;
         else
*/
         {
         p = DOT(bbo,grad)/MAX_CROSS_DER;
         if (fabs(p) > 1.e-30)
            p = SGN(p)*der_of_fcn_in_ind(fabs(p))/MAX_CROSS_DER;
         else
            p = 0.;
         for (i = 0; i < NVERT; i++)
            NDS(pelem->n[i],d) += der_of_fcn_in_res(fabs(res))*SGN(res)
        *(DOT(bb,bar[i])+q)*fcn_in_ind(fabs(DOT(bbo,grad))/MAX_CROSS_DER)*rdetB
        +fcn_in_res(fabs(res))*p*DOT(bbo,bar[i])*rdetB;
         }
      }
   }
}

double gradient_based_error_indicator(tGrid,u,use_bel)
GRID *tGrid;
INT u, use_bel;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], rdetB, err=0.;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
   if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                           && IS_FN(pelem->n[2])) || use_bel == YES){
      rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                      pelem->n[1]->myvertex->x,
                                      pelem->n[2]->myvertex->x,bar);
      sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
      err += sqrt(DOT(grad,grad))*rdetB;
   }
   return(err);
}

void derivatives_of_gradient_based_error_indicator(tGrid,u,d,use_bel)
GRID *tGrid;
INT u, d, use_bel; /* u ... solution at nodes, d ... derivatives at nodes  */
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], rdetB, q;
   INT i;

   set_value(tGrid,0.,d,0,U_TYPE); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
   if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                           && IS_FN(pelem->n[2])) || use_bel == YES){
      rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                      pelem->n[1]->myvertex->x,
                                      pelem->n[2]->myvertex->x,bar);
      sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
      q = sqrt(DOT(grad,grad));
      if (q > 1.e-30){
         for (i = 0; i < NVERT; i++)
            NDS(pelem->n[i],d) += DOT(grad,bar[i])/q*rdetB;
      }
   }
}

void compute_jumps_and_edge_lengths(tGrid,u,v,w,eps,use_bel)
GRID *tGrid;
FLOAT eps;
INT u, v, w, use_bel;  /* u ... P1 solution, v ... edge lengths, w ... jumps  */
{
   ELEMENT *pelem;
   FACE *pface;
   FLOAT a[DIM], b[DIM], bar[DIM2][DIM2], nn[DIM], grad[DIM], p;
   INT i, j, k;

   for (pface = FIRSTF(tGrid); pface; pface = pface->succ)
      FD(pface,v) = FD(pface,w) = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      barycentric_coordinates(pelem->n[0]->myvertex->x,
                              pelem->n[1]->myvertex->x,
                              pelem->n[2]->myvertex->x,bar);
      sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
      for (i = 0; i < SIDES; i++)
         if (IS_FF(pelem->f[i])){
            j = (i+1)%SIDES;
            k = (i+2)%SIDES;
            SUBTR(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,a);
            SUBTR(pelem->n[j]->myvertex->x,pelem->n[i]->myvertex->x,b);
            ORT_VECT_DIR(nn,a,b)
            if (FD(pelem->f[i],v) < 1.e-12)
               FD(pelem->f[i],v) = sqrt(DOT(a,a));
            FD(pelem->f[i],w) -= eps*DOT(nn,grad)/FD(pelem->f[i],v);
         }
   }
   if (use_bel == NO)
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
         for (i = 0; i < SIDES; i++)
            FD(pelem->f[i],w) = 0.;
}

double residual_based_error_estimator(tGrid,u,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT u, v, w, use_bel;
{
   ELEMENT *pelem;
   FACE *pface;
   FLOAT a[DIM], b[DIM], bar[DIM2][DIM2], nn[DIM], grad[DIM], p, 
         err=0., errK, alphaE, alphaK, b12=sqrt(beta), bm12, 
         epsm12=1./sqrt(eps);
   INT i, j, k;

if (USE_GRAD_IND == YES)
   err = res_and_der_based_error_indicator(tGrid,u,eps,bb0,bb1,react,rhs);
else{
//printf("\nvalues of error est\n");
   if (b12 > 0.)
      bm12 = 1./b12;
   compute_jumps_and_edge_lengths(tGrid,u,v,w,eps,use_bel);
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
   if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                           && IS_FN(pelem->n[2])) || use_bel == YES){
      barycentric_coordinates(pelem->n[0]->myvertex->x,
                              pelem->n[1]->myvertex->x,
                              pelem->n[2]->myvertex->x,bar);
      p = 0.;
      for (i = 0; i < SIDES; i++)
         if (IS_FF(pelem->f[i])){
            alphaE = epsm12*FD(pelem->f[i],v);
            if (alphaE*b12 > 1.)
               alphaE = bm12;
            p += alphaE*FD(pelem->f[i],v)*FD(pelem->f[i],w)*FD(pelem->f[i],w);
         }
      alphaK = epsm12*diameter(pelem);
      if (alphaK*b12 > 1.)
         alphaK = bm12;
      errK = p*epsm12
       + alphaK*alphaK*L2norm2_of_cd_res(pelem,bar,u,P1C,eps,bb0,bb1,react,rhs);
/*
printf("(%f,%f) (%f,%f) (%f,%f) %e\n",
        pelem->n[0]->myvertex->x[0],pelem->n[0]->myvertex->x[1],
        pelem->n[1]->myvertex->x[0],pelem->n[1]->myvertex->x[1],
        pelem->n[2]->myvertex->x[0],pelem->n[2]->myvertex->x[1],errK);
*/
      err += errK;
   }
// printf("error estimate = %e\n",sqrt(err));
}
   return(err);
}

double jump_indicator(tGrid,u,v,w,eps)
GRID *tGrid;
FLOAT eps;
INT u, v, w;
{
   ELEMENT *pelem;
   FLOAT err=0.;
   INT i;

   compute_jumps_and_edge_lengths(tGrid,u,v,w,eps,YES);
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      for (i = 0; i < SIDES; i++)
         if (IS_FF(pelem->f[i]))
            err += FD(pelem->f[i],v)*FD(pelem->f[i],v)*FD(pelem->f[i],w)*FD(pelem->f[i],w);
   return(err);
}

#if E_DATA & VECTOR_ELEMENT_DATA

void res_err_to_p1c_prolongation(tGrid,u,v,q)
GRID *tGrid;
INT u, v;  /*  u ... p0; v ... p1c  */
{
   ELEMENT *pel;
   NODE *pnode;

   set_value(tGrid,0.,v,0,U_TYPE);
   set_value(tGrid,0.,q,0,U_TYPE);
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
      NDS(pel->n[0],v) = MAX(NDS(pel->n[0],v),EDV(pel,u,0));
      NDS(pel->n[1],v) = MAX(NDS(pel->n[1],v),EDV(pel,u,0));
      NDS(pel->n[2],v) = MAX(NDS(pel->n[2],v),EDV(pel,u,0));
   }
}

#else

void res_err_to_p1c_prolongation(tGrid,u,v,q)
GRID *tGrid; INT u, v;
{  eprintf("Error: res_err_to_p1c_prolongation not available.\n"); return(0.);  }

#endif

void derivatives_of_residual_based_error_estimator(tGrid,u,d,g,v,w,
                                                   eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT u, d, g, v, w, use_bel; /* u ... solution at nodes, d ... derivatives at nodes  */
                /* g ... auxiliary scalar node variable  */
                /* v, w ... auxiliary scalar face variables  */
{
   ELEMENT *pelem;
   FACE *pface;
   FLOAT a[DIM], b[DIM], bar[DIM2][DIM2], nn[SIDES][DIM], grad[DIM], m[DIM], q, 
         alphaE, alphaK, b12=sqrt(beta), bm12, 
         eps12=sqrt(eps), epsm12=1./sqrt(eps), rdetB;
   INT i, j, k, ii=0;

if (USE_GRAD_IND == YES)
   derivatives_of_res_and_der_based_error_indicator(tGrid,u,d,
                                                         eps,bb0,bb1,react,rhs);
else{
   if (b12 > 0.)
      bm12 = 1./b12;
   compute_jumps_and_edge_lengths(tGrid,u,v,w,eps,use_bel);
   set_value(tGrid,0.,d,0,U_TYPE); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
   if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                           && IS_FN(pelem->n[2])) || use_bel == YES){
      rdetB = barycentric_coordinates(pelem->n[0]->myvertex->x,
                                      pelem->n[1]->myvertex->x,
                                      pelem->n[2]->myvertex->x,bar);
      SET7(m,0.)
      for (i = 0; i < SIDES; i++){
         j = (i+1)%SIDES;
         k = (i+2)%SIDES;
         SUBTR(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,a);
         SUBTR(pelem->n[j]->myvertex->x,pelem->n[i]->myvertex->x,b);
         ORT_VECT_DIR(nn[i],a,b)
         if (IS_FF(pelem->f[i])){
            alphaE = epsm12*FD(pelem->f[i],v);
            if (alphaE*b12 > 1.)
               alphaE = bm12;
            q = alphaE*FD(pelem->f[i],w)*eps12/rdetB;
            if (NOT_BF(pelem->f[i]))
               q += q;
            SET4(m,nn[i],q)
         }
      }
      alphaK = epsm12*diameter(pelem);
      if (alphaK*b12 > 1.)
         alphaK = bm12;
      q = 2.*alphaK*alphaK;
      for (i = 0; i < NVERT; i++)
         NDS(pelem->n[i],d) += DOT(m,nn[i])
            + q*cd_res_times_test_cd_res(pelem,bar,u,g,i,eps,bb0,bb1,react,rhs);
ii++;
   }
//printf("%i triangles for error est.\n",ii);
/*
   printf("norm %e\n",sqrt(dot(tGrid,d,d,T_FOR_U,U_TYPE)));
NODE *pnode;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
   printf("(%f,%f) %e\n",pnode->myvertex->x[0],pnode->myvertex->x[1],
          NDS(pnode,d));
   for (pnode = FIRSTNODE(tGrid); pnode; pnode=pnode->succ)
   printf("(%f,%f) u = %e\n",pnode->myvertex->x[0],pnode->myvertex->x[1],
          NDS(pnode,u));
*/
}
}

void compute_b_and_sign_for_edge(pelem,i,bb0,bb1,bb,sign,a)
ELEMENT *pelem;
FLOAT (*bb0)(), (*bb1)(), *bb, *sign, *a;
INT i;
{
   FLOAT b[DIM], nn[DIM], xjk[DIM];
   INT j, k;

   j = (i+1)%SIDES;
   k = (i+2)%SIDES;
   AVERAGE(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,xjk)
   SUBTR(pelem->n[j]->myvertex->x,pelem->n[k]->myvertex->x,a);
   SUBTR(pelem->n[j]->myvertex->x,pelem->n[i]->myvertex->x,b);
   ORT_VECT_DIR(nn,a,b)
   bb[0] = bb0(xjk);
   bb[1] = bb1(xjk);
   *sign = DOT(bb,nn)/(fabs(nn[0])+fabs(nn[1]));
   if (fabs(*sign) < 1.e-30)
      *sign = 0.;
   else if ((*sign) > 0.)
      *sign = 1.;
   else
      *sign = -1.;
}

void compute_conv_jumps_and_edge_lengths(tGrid,u,v,w,bb0,bb1,use_bel)
GRID *tGrid;
FLOAT (*bb0)(), (*bb1)();
INT u, v, w, use_bel;  /* u ... P1 solution, v ... edge lengths, w ... jumps  */
{
   ELEMENT *pelem;
   FACE *pface;
   FLOAT a[DIM], bb[DIM], bar[DIM2][DIM2], grad[DIM], q;
   INT i, j, k;

   for (pface = FIRSTF(tGrid); pface; pface = pface->succ)
      FD(pface,v) = FD(pface,w) = 0.;
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      barycentric_coordinates(pelem->n[0]->myvertex->x,
                              pelem->n[1]->myvertex->x,
                              pelem->n[2]->myvertex->x,bar);
      sp1c_grad(pelem->n[0],pelem->n[1],pelem->n[2],bar,u,grad);
      for (i = 0; i < SIDES; i++)
         if (NOT_BF(pelem->f[i])){
            compute_b_and_sign_for_edge(pelem,i,bb0,bb1,bb,&q,a);
            FD(pelem->f[i],w) += q*DOT(bb,grad);
            if (FD(pelem->f[i],v) < 1.e-30)
               FD(pelem->f[i],v) = sqrt(DOT(a,a));
         }
   }
   if (use_bel == NO)
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
         if (!(IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                                  && IS_FN(pelem->n[2])))
            FD(pelem->f[0],w) = FD(pelem->f[1],w) = FD(pelem->f[2],w) = 0.;
}

DOUBLE conv_jumps_based_error_indicator(tGrid,u,v,w,bb0,bb1,use_bel)
GRID *tGrid;
FLOAT (*bb0)(), (*bb1)();
INT u, v, w, use_bel;  /* u ... P1 solution, v ... edge lengths, w ... jumps  */
{
   FACE *pface;
   DOUBLE sum=0.;

   compute_conv_jumps_and_edge_lengths(tGrid,u,v,w,bb0,bb1,use_bel);
   for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
      sum += FD(pface,w)*FD(pface,w)*FD(pface,v);
   return(sum);
}

void derivatives_of_conv_jumps_based_error_indicator(tGrid,u,d,v,w,bb0,bb1,
                                                                        use_bel)
GRID *tGrid;
FLOAT (*bb0)(), (*bb1)();
INT u, d, v, w, use_bel; /* u ... solution at nodes, d ... derivatives at nodes  */
                /* v, w ... auxiliary scalar face variables  */
{
   ELEMENT *pelem;
   FLOAT a[DIM], bb[DIM], bar[DIM2][DIM2], q;
   INT i, j, k;

   compute_conv_jumps_and_edge_lengths(tGrid,u,v,w,bb0,bb1,use_bel);
   set_value(tGrid,0.,d,0,U_TYPE); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
   if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                           && IS_FN(pelem->n[2])) || use_bel == YES){
      barycentric_coordinates(pelem->n[0]->myvertex->x,
                              pelem->n[1]->myvertex->x,
                              pelem->n[2]->myvertex->x,bar);
      for (i = 0; i < SIDES; i++){
         compute_b_and_sign_for_edge(pelem,i,bb0,bb1,bb,&q,a);
         q *= 2.*FD(pelem->f[i],v)*FD(pelem->f[i],w);
         for (j = 0; j < NVERT; j++)
            NDS(pelem->n[j],d) += q*DOT(bb,bar[j]);
      }
   }
}

#else

double residual_based_error_estimator(tGrid,u,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, v, w, use_bel;
{  eprintf("Error: residual_based_error_estimator not available.\n"); return(0.);  }

void derivatives_of_residual_based_error_estimator(tGrid,u,d,g,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, d, v, w, use_bel;
{  eprintf("Error: derivatives_of_residual_based_error_estimator not available.\n");  }

DOUBLE conv_jumps_based_error_indicator(tGrid,u,v,w,bb0,bb1,use_bel)
GRID *tGrid; FLOAT (*bb0)(), (*bb1)(); INT u, v, w, use_bel;
{  eprintf("Error: conv_jumps_based_error_indicator not available.\n"); return(0.);  }

void derivatives_of_conv_jumps_based_error_indicator(tGrid,u,d,v,w,bb0,bb1,use_bel)
GRID *tGrid; FLOAT (*bb0)(), (*bb1)(); INT u, d, v, w, use_bel;
{  eprintf("Error: derivatives_of_conv_jumps_based_error_indicatornot available.\n");  }

#endif

#else

double residual_based_error_estimator(tGrid,u,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, v, w, use_bel;
{  eprintf("Error: residual_based_error_estimator not available.\n"); return(0.);  }

void derivatives_of_residual_based_error_estimator(tGrid,u,d,g,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, d, v, w, use_bel;
{  eprintf("Error: derivatives_of_residual_based_error_estimator not available.\n");  }

#endif
