/******************************************************************************/
/*                                                                            */
/*                              error estimators                              */
/*                                                                            */
/******************************************************************************/

DOUBLE fcn_in_ind(x)
DOUBLE x;
{
   if (x > 1.)
      return(sqrt(x));
   else
      return(2.5*x*x-1.5*x*x*x);
}

DOUBLE der_of_fcn_in_ind(x)
DOUBLE x;
{
   if (x > 1.)
      return(0.5/sqrt(x));
   else
      return(5.*x-4.5*x*x);
}

double res_and_der_based_error_indicator(tGrid,u,eps,bb0,bb1,react,rhs)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u;
{
   ELEMENT *pelem;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], bbo[DIM], 
         rdetB, q, res, err=0.;

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1])
                              && IS_FN(pelem->n[2])) || USE_BEL == YES){
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
         res = 0.;
         res = DOT(bb,grad) + 
               react(xc)*(NDS(pelem->n[0],u)+NDS(pelem->n[1],u)
                                            +NDS(pelem->n[2],u))/3. - rhs(xc);
         res = res*res*rdetB;
         err += res + fcn_in_ind(fabs(DOT(bbo,grad)))*rdetB;
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
         rdetB, res, q, err=0.;
   INT i;

   set_value(tGrid,0.,d,0,U_TYPE); 
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1])
                              && IS_FN(pelem->n[2])) || USE_BEL == YES){
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
         for (i = 0; i < NVERT; i++)
            NDS(pelem->n[i],d) += 2.*res*(DOT(bb,bar[i])+q)*rdetB;
         q = DOT(bbo,grad);
         if (fabs(q) > 1.e-30){
            q = SGN(q)*der_of_fcn_in_ind(fabs(q));
            for (i = 0; i < NVERT; i++)
               NDS(pelem->n[i],d) += q*DOT(bbo,bar[i])*rdetB;
         }
   }
}










DOUBLE dot();

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

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

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

#if U_SPACE == P1C

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
         err += errK;
      }
   }
   return(err);
}

#elif U_SPACE == P2C

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
      if (b12 > 0.)
         bm12 = 1./b12;
      //compute_jumps_and_edge_lengths(tGrid,u,v,w,eps,use_bel); // u from P1C
      for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ)
      if ((IS_FN(pelem->n[0]) && IS_FN(pelem->n[1]) 
                              && IS_FN(pelem->n[2])) || use_bel == YES){
         barycentric_coordinates(pelem->n[0]->myvertex->x,
                                 pelem->n[1]->myvertex->x,
                                 pelem->n[2]->myvertex->x,bar);
         p = 0.;
         for (i = 0; i < SIDES; i++)
            if (IS_FF(pelem->f[i])){
               //alphaE = epsm12*FD(pelem->f[i],v);
               if (alphaE*b12 > 1.)
                  alphaE = bm12;
               //p += alphaE*FD(pelem->f[i],v)*FD(pelem->f[i],w)*FD(pelem->f[i],w);
            }
         alphaK = epsm12*diameter(pelem);
         if (alphaK*b12 > 1.)
            alphaK = bm12;
         errK = p*epsm12
          + alphaK*alphaK*L2norm2_of_cd_res(pelem,bar,u,P2C,eps,bb0,bb1,react,rhs);
         err += errK;
      }
   }
   return(err);
}

#else

#endif

#if U_SPACE == P1C

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

#elif U_SPACE == P2C

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
   INT i, j, k;

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
         for (i = 0; i < NVERT; i++){
            NDS(pelem->n[i],d) += DOT(m,nn[i])
               + q*cd_res_times_p2c_test_cd_res(pelem,bar,u,g,i,eps,bb0,bb1,react,rhs);
            FD(pelem->f[i],d) += DOT(m,nn[i])
               + q*cd_res_times_p2c_ftest_cd_res(pelem,bar,u,g,i,eps,bb0,bb1,react,rhs);
         }
      }
   }
}

#else

#endif

#else

double residual_based_error_estimator(tGrid,u,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, v, w, use_bel;
{  eprintf("Error: residual_based_error_estimator not available.\n"); return(0.)  }

void derivatives_of_residual_based_error_estimator(tGrid,u,d,g,v,w,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT u, d, v, w, use_bel;
{  eprintf("Error: derivatives_of_residual_based_error_estimator not available.\n");  }

#endif
