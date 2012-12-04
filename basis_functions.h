/*
      v gscalar_conv_diff_res chybi Laplace !!!!!!!!!
*/

/******************************************************************************/
/*                                                                            */
/*                              basis functions                               */
/*                                                                            */
/******************************************************************************/

FLOAT pq3();

FLOAT fcn_lin_value(x,a,c,f)
FLOAT *x, a[][DIM], *c, (*f)();
{
   FLOAT y[DIM];

   V_LIN_VALUE(y,a,c,x[0],x[1])
   return(f(y));
}

FLOAT fcn_bilin_value(x,a,c,alpha,f)
FLOAT *x, a[][DIM], *c, *alpha, (*f)();
{
   FLOAT y[DIM];

   V_BILIN_VALUE(y,a,c,alpha,x[0],x[1])
   return(f(y));
}

FLOAT fcn_ref_map_value(x,ref_map,f)
FLOAT *x, (*f)();
REF_MAPPING *ref_map;
{
   if (ref_map->type == Q1_REF_MAP)
      return(fcn_bilin_value(x,ref_map->q1_a,ref_map->q1_c,ref_map->q1_alpha,f));
   else if (ref_map->type == P1_REF_MAP)
      return(fcn_lin_value(x,ref_map->p1_a,ref_map->p1_c,f));
   else
      eprintf("Error: fcn_ref_map_value not available.\n");
}

void ref_map_point(x,y,ref_map)
FLOAT *x, *y;
REF_MAPPING *ref_map;
{
   if (ref_map->type == Q1_REF_MAP)
      V_BILIN_VALUE(y,ref_map->q1_a,ref_map->q1_c,ref_map->q1_alpha,x[0],x[1])
   else if (ref_map->type == P1_REF_MAP)
      V_LIN_VALUE(y,ref_map->p1_a,ref_map->p1_c,x[0],x[1])
   else
      eprintf("Error: ref_map_point not available.\n");
}

FLOAT no_rd(v,ref_map)
FLOAT (*v)(); REF_MAPPING *ref_map;
{ eprintf("Error: rd not available.\n"); return(0.);  }

#if DIM == 2

FLOAT p32();
void points();

FLOAT zero_fcn(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT one_fcn(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_l0(x)
FLOAT x[DIM];
{ return(1. - x[0] - x[1]); }

FLOAT r_l1(x)
FLOAT x[DIM];
{ return(x[0]); }

FLOAT r_l2(x)
FLOAT x[DIM];
{ return(x[1]); }

FLOAT r_l0_0(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_l0_1(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_l1_0(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_l1_1(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l2_0(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l2_1(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_l0l1(x)
FLOAT x[DIM];
{ return(x[0]*(1. - x[0] - x[1])); }

FLOAT r_l0l2(x)
FLOAT x[DIM];
{ return(x[1]*(1. - x[0] - x[1])); }

FLOAT r_l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]); }

FLOAT r_l0l1_0(x)
FLOAT x[DIM];
{ return(1. - 2.*x[0] - x[1]); }

FLOAT r_l0l1_1(x)
FLOAT x[DIM];
{ return(-x[0]); }

FLOAT r_l0l1_00(x)
FLOAT x[DIM];
{ return(-2.); }

FLOAT r_l0l1_01(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_l0l1_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l0l2_0(x)
FLOAT x[DIM];
{ return(-x[1]); }

FLOAT r_l0l2_1(x)
FLOAT x[DIM];
{ return(1. - x[0] - 2.*x[1]); }

FLOAT r_l0l2_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l0l2_01(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_l0l2_11(x)
FLOAT x[DIM];
{ return(-2.); }

FLOAT r_l1l2_0(x)
FLOAT x[DIM];
{ return(x[1]); }

FLOAT r_l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]); }

FLOAT r_l1l2_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l1l2_01(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_l1l2_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_l0l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*(1.-x[0]-x[1])); }

FLOAT r_l0l1l2_0(x)
FLOAT x[DIM];
{ return(x[1]*(1. - 2.*x[0] - x[1])); }

FLOAT r_l0l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]*(1. - x[0] - 2.*x[1])); }

FLOAT r_l0l1l2_00(x)
FLOAT x[DIM];
{ return(-2.*x[1]); }

FLOAT r_l0l1l2_01(x)
FLOAT x[DIM];
{ return(1. - 2.*x[0] - 2.*x[1]); }

FLOAT r_l0l1l2_11(x)
FLOAT x[DIM];
{ return(-2.*x[0]); }

FLOAT rd_l0l1l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ //FLOAT x[DIM]={0.3333333333,0.3333333333}; 
  //return(fcn_ref_map_value(x,ref_map,v)); 
  eprintf("Error: rd_l0l1l2 not available.\n"); return(0.);  }

FLOAT r_l0l0l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*(1.-x[0]-x[1])*(1.-x[0]-x[1])); }

FLOAT r_l0l0l1l2_0(x)
FLOAT x[DIM];
{ return(x[1]*(1.-x[0]-x[1])*(1.-3.*x[0]-x[1])); }

FLOAT r_l0l0l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]*(1.-x[0]-x[1])*(1.-x[0]-3.*x[1])); }

FLOAT r_l0l0l1l2_00(x)
FLOAT x[DIM];
{ return(x[1]*(6.*x[0]+4.*x[1]-4.)); }

FLOAT r_l0l0l1l2_01(x)
FLOAT x[DIM];
{ return(3.*(x[0]*x[0]+x[1]*x[1]) + 8.*x[0]*x[1] - 4.*(x[0]+x[1]) + 1.); }

FLOAT r_l0l0l1l2_11(x)
FLOAT x[DIM];
{ return(x[0]*(4.*x[0]+6.*x[1]-4.)); }

FLOAT rd_l0l0l1l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l0l1l2 not available.\n"); return(0.);  }

FLOAT r_l0l1l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*x[1]*(1.-x[0]-x[1])); }

FLOAT r_l0l1l1l2_0(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*(2.-3.*x[0]-2.*x[1])); }

FLOAT r_l0l1l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*(1.-x[0]-2.*x[1])); }

FLOAT r_l0l1l1l2_00(x)
FLOAT x[DIM];
{ return(x[1]*(2.-6.*x[0]-2.*x[1])); }

FLOAT r_l0l1l1l2_01(x)
FLOAT x[DIM];
{ return(x[0]*(2.-3.*x[0]-4.*x[1])); }

FLOAT r_l0l1l1l2_11(x)
FLOAT x[DIM];
{ return(-2.*x[0]*x[0]); }

FLOAT rd_l0l1l1l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l1l1l2 not available.\n"); return(0.);  }

FLOAT r_l0l1l2l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*x[1]*(1.-x[0]-x[1])); }

FLOAT r_l0l1l2l2_0(x)
FLOAT x[DIM];
{ return(x[1]*x[1]*(1.-2.*x[0]-x[1])); }

FLOAT r_l0l1l2l2_1(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*(2.-2.*x[0]-3.*x[1])); }

FLOAT r_l0l1l2l2_00(x)
FLOAT x[DIM];
{ return(-2.*x[1]*x[1]); }

FLOAT r_l0l1l2l2_01(x)
FLOAT x[DIM];
{ return(x[1]*(2.-4.*x[0]-3.*x[1])); }

FLOAT r_l0l1l2l2_11(x)
FLOAT x[DIM];
{ return(x[0]*(2.-2.*x[0]-6.*x[1])); }

FLOAT rd_l0l1l2l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l1l2l2 not available.\n"); return(0.);  }

FLOAT r_l0l0l0l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*(1.-x[0]-x[1])*(1.-x[0]-x[1])*(1.-x[0]-x[1])); }

FLOAT r_l0l0l0l1l2_0(x)
FLOAT x[DIM];
{ return(x[1]*(1.-x[0]-x[1])*(1.-x[0]-x[1])*(1.-4.*x[0]-x[1])); }

FLOAT r_l0l0l0l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]*(1.-x[0]-x[1])*(1.-x[0]-x[1])*(1.-x[0]-4.*x[1])); }

FLOAT r_l0l0l0l1l2_00(x)
FLOAT x[DIM];
{ return(6.*x[1]*(1.-x[0]-x[1])*(2.*x[0]+x[1]-1.)); }

FLOAT r_l0l0l0l1l2_01(x)
FLOAT x[DIM];
{ return(6.*x[0]*x[1]*(1.-x[0]-x[1])
         +(1.-4.*x[0]-4.*x[1])*(1.-x[0]-x[1])*(1.-x[0]-x[1])); }

FLOAT r_l0l0l0l1l2_11(x)
FLOAT x[DIM];
{ return(6.*x[0]*(1.-x[0]-x[1])*(x[0]+2.*x[1]-1.)); }

FLOAT rd_l0l0l0l1l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l0l0l1l2 not available.\n"); return(0.);  }

FLOAT r_l0l1l1l1l2(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*x[0]*x[1]*(1.-x[0]-x[1])); }

FLOAT r_l0l1l1l1l2_0(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*x[1]*(3.-4.*x[0]-3.*x[1])); }

FLOAT r_l0l1l1l1l2_1(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*x[0]*(1.-x[0]-2.*x[1])); }

FLOAT r_l0l1l1l1l2_00(x)
FLOAT x[DIM];
{ return(6.*x[0]*x[1]*(1.-2.*x[0]-x[1])); }

FLOAT r_l0l1l1l1l2_01(x)
FLOAT x[DIM];
{ return(x[0]*x[0]*(3.-4.*x[0]-6.*x[1])); }

FLOAT r_l0l1l1l1l2_11(x)
FLOAT x[DIM];
{ return(-2.*x[0]*x[0]*x[0]); }

FLOAT rd_l0l1l1l1l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l1l1l1l2 not available.\n"); return(0.);  }

FLOAT r_l0l1l2l2l2(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*x[1]*x[1]*(1.-x[0]-x[1])); }

FLOAT r_l0l1l2l2l2_0(x)
FLOAT x[DIM];
{ return(x[1]*x[1]*x[1]*(1.-2.*x[0]-x[1])); }

FLOAT r_l0l1l2l2l2_1(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*x[1]*(3.-3.*x[0]-4.*x[1])); }

FLOAT r_l0l1l2l2l2_00(x)
FLOAT x[DIM];
{ return(-2.*x[1]*x[1]*x[1]); }

FLOAT r_l0l1l2l2l2_01(x)
FLOAT x[DIM];
{ return(x[1]*x[1]*(3.-6.*x[0]-4.*x[1])); }

FLOAT r_l0l1l2l2l2_11(x)
FLOAT x[DIM];
{ return(6.*x[0]*x[1]*(1.-x[0]-2.*x[1])); }

FLOAT rd_l0l1l2l2l2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_l0l1l2l2l2 not available.\n"); return(0.);  }

FLOAT nc_ze(x,b0,b1,b2)  /*  nonconforming P1 basis function ze_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2];
{
   FLOAT l0;

   l0 = b0[0]*x[0] + b0[1]*x[1] + b0[2];
   return(1. - 2.*l0 );
}

FLOAT nc_gze(x,i,b0,b1,b2)  /*  x[i] derivative of ze_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2];
INT i;
{
   return(-2.*b0[i]);
}

FLOAT nc_bgze(x,b0,b1,b2,bb0,bb1)  /*  bb*grad ze_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
{
   return( bb0(x)*nc_gze(x,0,b0,b1,b2) + bb1(x)*nc_gze(x,1,b0,b1,b2) );
}

FLOAT nc_fi(x,b0,b1,b2)  /*  basis function fi_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2];
{
   FLOAT l0, l1, l2;

   l0 = b0[0]*x[0] + b0[1]*x[1] + b0[2];
   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return(1. - 2.*l0 - 10.*(l1*l1*l0 - l1*l0*l0)
                     - 10.*(l2*l2*l0 - l2*l0*l0) );
}

FLOAT nc_gfi(x,i,b0,b1,b2)  /*  x[i] derivative of fi_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2];
INT i;
{
   FLOAT l0, l1, l2;

   l0 = b0[0]*x[0] + b0[1]*x[1] + b0[2];
   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return(-2.*b0[i] - 20.*l1*l0*(b1[i]-b0[i]) - 10.*(l1*l1*b0[i]-l0*l0*b1[i])
                    - 20.*l2*l0*(b2[i]-b0[i]) - 10.*(l2*l2*b0[i]-l0*l0*b2[i]));
}

FLOAT nc_lfi(x,b0,b1,b2)  /*  laplacian of fi_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2];
{
   FLOAT l0, l1, l2;

   l0 = b0[0]*x[0] + b0[1]*x[1] + b0[2];
   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return(  -40.*(l1-l0)*DOT(b1,b0) + 20.*l1*DOT(b0,b0) - 20.*l0*DOT(b1,b1) - 
             40.*(l2-l0)*DOT(b2,b0) + 20.*l2*DOT(b0,b0) - 20.*l0*DOT(b2,b2) );
}

FLOAT nc_bgfi(x,b0,b1,b2,bb0,bb1)  /*  bb*grad fi_12  */
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
{
   return( bb0(x)*nc_gfi(x,0,b0,b1,b2) + bb1(x)*nc_gfi(x,1,b0,b1,b2) );
}

/*  -nu*lapl fi_12 + bb*grad fi_12 + c*fi_12                                  */
FLOAT nc_lbgfi(x,b0,b1,b2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
{
   return( -nu*nc_lfi(x,b0,b1,b2) + nc_bgfi(x,b0,b1,b2,bb0,bb1) + 
                                                       c(x)*nc_fi(x,b0,b1,b2) );
}

FLOAT nc_chi(x,b1,b2,n1,n2)  /*  basis function chi_12  */
FLOAT x[DIM], b1[DIM2], b2[DIM2];
NODE *n1, *n2;
{
   FLOAT l1, l2;

   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return( (l1*l1*l2 - l1*l2*l2)*NINDI(n2,n1) );
}

FLOAT nc_gchi(x,i,b1,b2,n1,n2)  /*  x[i] derivative of chi_12  */
FLOAT x[DIM], b1[DIM2], b2[DIM2];
NODE *n1, *n2;
INT i;
{
   FLOAT l1, l2;

   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return( ( 2.*l1*l2*(b1[i]-b2[i]) + l1*l1*b2[i]-l2*l2*b1[i] )*NINDI(n2,n1) );
}

FLOAT nc_lchi(x,b1,b2,n1,n2)  /*  laplacian of chi_12  */
FLOAT x[DIM], b1[DIM2], b2[DIM2];
NODE *n1, *n2;
{
   FLOAT l1, l2;

   l1 = b1[0]*x[0] + b1[1]*x[1] + b1[2];
   l2 = b2[0]*x[0] + b2[1]*x[1] + b2[2];
   return( ( 4.*(l1-l2)*DOT(b1,b2) - 2.*l1*DOT(b2,b2) + 
                                     2.*l2*DOT(b1,b1) )*NINDI(n2,n1) );
}

FLOAT nc_bgchi(x,b1,b2,n1,n2,bb0,bb1)  /*  bb*grad chi_12  */
FLOAT x[DIM], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( bb0(x)*nc_gchi(x,0,b1,b2,n1,n2) + bb1(x)*nc_gchi(x,1,b1,b2,n1,n2) );
}

/*  -nu*lapl chi_12 + bb*grad chi_12 + c*chi_12                               */
FLOAT nc_lbgchi(x,b1,b2,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n1, *n2;
{
   return( -nu*nc_lchi(x,b1,b2,n1,n2) + nc_bgchi(x,b1,b2,n1,n2,bb0,bb1) +
                                                   c(x)*nc_chi(x,b1,b2,n1,n2) );
}

FLOAT nc_ze12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_ze(x,b0,b1,b2) );
}

FLOAT nc_fi12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_fi(x,b0,b1,b2) );
}

FLOAT nc_chi12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_chi(x,b1,b2,n1,n2) );
}

FLOAT nc_bgze12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_bgze(x,b0,b1,b2,bb0,bb1) );
}

FLOAT nc_bgfi12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_bgfi(x,b0,b1,b2,bb0,bb1) );
}

FLOAT nc_bgchi12(x,b0,b1,b2,n1,n2,bb0,bb1)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)();
NODE *n1, *n2;
{
   return( nc_bgchi(x,b1,b2,n1,n2,bb0,bb1) );
}

FLOAT nc_lbgfi12(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgfi(x,b0,b1,b2,bb0,bb1,c,nu) );
}

FLOAT nc_lbgfi02(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgfi(x,b1,b0,b2,bb0,bb1,c,nu) );
}

FLOAT nc_lbgfi01(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgfi(x,b2,b0,b1,bb0,bb1,c,nu) );
}

FLOAT nc_lbgchi12(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgchi(x,b1,b2,n1,n2,bb0,bb1,c,nu) );
}

FLOAT nc_lbgchi02(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgchi(x,b0,b2,n0,n2,bb0,bb1,c,nu) );
}

FLOAT nc_lbgchi01(x,b0,b1,b2,n0,n1,n2,bb0,bb1,c,nu)
FLOAT x[DIM], b0[DIM2], b1[DIM2], b2[DIM2], (*bb0)(), (*bb1)(), (*c)(), nu;
NODE *n0, *n1, *n2;
{
   return( nc_lbgchi(x,b0,b1,n0,n1,bb0,bb1,c,nu) );
}

FLOAT p1mod_function(x,a1,a2,a3,b1,b2,b3,n1,n2,n3,b)
FLOAT x[DIM], a1, a2, a3, b1, b2, b3, b[DIM2][DIM2];
NODE *n1, *n2, *n3; 
{
   return( a1*nc_fi(x,b[0],b[1],b[2]) +
           a2*nc_fi(x,b[1],b[2],b[0]) +
           a3*nc_fi(x,b[2],b[0],b[1]) +
           b1*nc_chi(x,b[1],b[2],n2,n3) +
           b2*nc_chi(x,b[0],b[2],n1,n3) +
           b3*nc_chi(x,b[0],b[1],n1,n2));
}

/* p1 basis funtions */

FLOAT r_p1_0(x)
FLOAT x[DIM];
{ return(1. - x[0] - x[1]); }

FLOAT r_p1_1(x)
FLOAT x[DIM];
{ return(x[0]); }

FLOAT r_p1_2(x)
FLOAT x[DIM];
{ return(x[1]); }

FLOAT r_p1_0_0(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_p1_0_1(x)
FLOAT x[DIM];
{ return(-1.); }

FLOAT r_p1_1_0(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_p1_1_1(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_2_0(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_2_1(x)
FLOAT x[DIM];
{ return(1.); }

FLOAT r_p1_0_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_0_01(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_0_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_1_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_1_01(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_1_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_2_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_2_01(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p1_2_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT rd_p1_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.,0.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_p1_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,0.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_p1_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.,1.}; return(fcn_ref_map_value(x,ref_map,v)); }

/* p1x3 basis funtions */

FLOAT r_p1x3_b(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(3.*(1.-x[0]-x[1]));
  else if (x[0] < x[1]) return(3.*x[0]); else return(3.*x[1]); }

FLOAT r_p1x3_b_0(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-3.);
  else if (x[0] < x[1]) return(3.); else return(0.); }

FLOAT r_p1x3_b_1(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-3.);
  else if (x[0] < x[1]) return(0.); else return(3.); }

FLOAT rd_p1x3_b(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, 
        x[DIM]={0.3333333333333333333333,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)-(fcn_ref_map_value(x0,ref_map,v)+
         fcn_ref_map_value(x1,ref_map,v)+fcn_ref_map_value(x2,ref_map,v))/3.); }

/* p2 basis funtions */

FLOAT rd_p2_f01(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x[DIM]={0.5,0.}; 
  return(4.*fcn_ref_map_value(x,ref_map,v) 
   - 2.*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v))); }

FLOAT rd_p2_f02(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x2[DIM]={0.,1.}, x[DIM]={0.,0.5}; 
  return(4.*fcn_ref_map_value(x,ref_map,v) 
   - 2.*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x2,ref_map,v))); }

FLOAT rd_p2_f12(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, x[DIM]={0.5,0.5}; 
  return(4.*fcn_ref_map_value(x,ref_map,v) 
   - 2.*(fcn_ref_map_value(x1,ref_map,v) + fcn_ref_map_value(x2,ref_map,v))); }

/* p2x3 basis funtions */

FLOAT r_p2x3_0(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(0.);
  else if (x[0] < x[1]) return(3.*x[0]*(1.-2.*x[0]-x[1]));
  else return(3.*x[1]*(1.-x[0]-2.*x[1])); }

FLOAT r_p2x3_1(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.)
     return(3.*(1.-x[0]-x[1])*(2.*x[0]+x[1]-1.));
  else if (x[0] < x[1]) return(0.); else return(3.*x[1]*(x[0]-x[1])); }

FLOAT r_p2x3_2(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.)
     return(3.*(1.-x[0]-x[1])*(x[0]+2.*x[1]-1.));
  else if (x[0] < x[1]) return(3.*x[0]*(x[1]-x[0])); else return(0.); }

FLOAT r_p2x3_0_0(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(0.);
  else if (x[0] < x[1]) return(3.*(1.-4.*x[0]-x[1])); else return(-3.*x[1]); }

FLOAT r_p2x3_0_1(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(0.);
  else if (x[0] < x[1]) return(-3.*x[0]); else return(3.*(1.-x[0]-4.*x[1])); }

FLOAT r_p2x3_1_0(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(9.-12.*x[0]-9.*x[1]);
  else if (x[0] < x[1]) return(0.); else return(3.*x[1]); }

FLOAT r_p2x3_1_1(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(6.-9.*x[0]-6.*x[1]);
  else if (x[0] < x[1]) return(0.); else return(3.*x[0]-6.*x[1]); }

FLOAT r_p2x3_2_0(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(6.-6.*x[0]-9.*x[1]);
  else if (x[0] < x[1]) return(3.*x[1]-6.*x[0]); else return(0.); }

FLOAT r_p2x3_2_1(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(9.-9.*x[0]-12.*x[1]);
  else if (x[0] < x[1]) return(3.*x[0]); else return(0.); }

FLOAT r_p2x3_0_00(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] < 1. && x[0] < x[1]) return(-12.); else return(0.); }

FLOAT r_p2x3_0_01(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(0.); else return(-3.); }

FLOAT r_p2x3_0_11(x)
FLOAT x[DIM];
{ if (x[0]+2.*x[1] < 1. && x[0] > x[1]) return(-12.); else return(0.); }

FLOAT r_p2x3_1_00(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-12.);
  else return(0.); }

FLOAT r_p2x3_1_01(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-9.);
  else if (x[0] < x[1]) return(0.); else return(3.); }

FLOAT r_p2x3_1_11(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. || x[0] > x[1]) return(-6.); else return(0.); }

FLOAT r_p2x3_2_00(x)
FLOAT x[DIM];
{ if (x[0]+2.*x[1] > 1. || x[0] < x[1]) return(-6.); else return(0.); }

FLOAT r_p2x3_2_01(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-9.);
  else if (x[0] < x[1]) return(3.); else return(0.); }

FLOAT r_p2x3_2_11(x)
FLOAT x[DIM];
{ if (2.*x[0]+x[1] > 1. && x[0]+2.*x[1] > 1.) return(-12.); else return(0.); }

FLOAT rd_p2x3_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, 
        x01[DIM]={0.5,0.}, x12[DIM]={0.5,0.5}, x20[DIM]={0.,0.5}, 
        x012[DIM]={0.3333333333333333333333,0.3333333333333333333333},
        x[DIM]={0.1666666666666666666666,0.1666666666666666666666}, s;
  s =((fcn_ref_map_value(x01,ref_map,v)+fcn_ref_map_value(x12,ref_map,v)+
       fcn_ref_map_value(x20,ref_map,v))*4. -fcn_ref_map_value(x0,ref_map,v) -
       fcn_ref_map_value(x1,ref_map,v) -fcn_ref_map_value(x2,ref_map,v))/18.+
       fcn_ref_map_value(x012,ref_map,v)*0.5;
  return(4.*(fcn_ref_map_value(x,ref_map,v) - s +
      (fcn_ref_map_value(x12,ref_map,v)-fcn_ref_map_value(x0,ref_map,v))/3.)); }

FLOAT rd_p2x3_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, 
        x01[DIM]={0.5,0.}, x12[DIM]={0.5,0.5}, x20[DIM]={0.,0.5}, 
        x012[DIM]={0.3333333333333333333333,0.3333333333333333333333},
        x[DIM]={0.6666666666666666666666,0.1666666666666666666666}, s;
  s =((fcn_ref_map_value(x01,ref_map,v)+fcn_ref_map_value(x12,ref_map,v)+
       fcn_ref_map_value(x20,ref_map,v))*4. -fcn_ref_map_value(x0,ref_map,v) -
       fcn_ref_map_value(x1,ref_map,v) -fcn_ref_map_value(x2,ref_map,v))/18.+
       fcn_ref_map_value(x012,ref_map,v)*0.5;
  return(4.*(fcn_ref_map_value(x,ref_map,v) - s +
      (fcn_ref_map_value(x20,ref_map,v)-fcn_ref_map_value(x1,ref_map,v))/3.)); }

FLOAT rd_p2x3_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, 
        x01[DIM]={0.5,0.}, x12[DIM]={0.5,0.5}, x20[DIM]={0.,0.5}, 
        x012[DIM]={0.3333333333333333333333,0.3333333333333333333333},
        x[DIM]={0.1666666666666666666666,0.6666666666666666666666}, s;
  s =((fcn_ref_map_value(x01,ref_map,v)+fcn_ref_map_value(x12,ref_map,v)+
       fcn_ref_map_value(x20,ref_map,v))*4. -fcn_ref_map_value(x0,ref_map,v) -
       fcn_ref_map_value(x1,ref_map,v) -fcn_ref_map_value(x2,ref_map,v))/18.+
       fcn_ref_map_value(x012,ref_map,v)*0.5;
  return(4.*(fcn_ref_map_value(x,ref_map,v) - s +
      (fcn_ref_map_value(x01,ref_map,v)-fcn_ref_map_value(x2,ref_map,v))/3.)); }

FLOAT rd_p2x3_b(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={0.,0.}, x1[DIM]={1.,0.}, x2[DIM]={0.,1.}, 
        x01[DIM]={0.5,0.}, x12[DIM]={0.5,0.5}, x20[DIM]={0.,0.5}, 
        x[DIM]={0.3333333333333333333333,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v) + (fcn_ref_map_value(x0,ref_map,v) + 
       fcn_ref_map_value(x1,ref_map,v) +fcn_ref_map_value(x2,ref_map,v) -
      (fcn_ref_map_value(x01,ref_map,v)+fcn_ref_map_value(x12,ref_map,v)+
       fcn_ref_map_value(x20,ref_map,v))*4.)/9.); }

/* some p3 basis functions */

FLOAT r_p3_01(x)
FLOAT x[DIM];
{ return((1. - x[1])*x[0]*x[1]); }

FLOAT r_p3_02(x)
FLOAT x[DIM];
{ return(x[0]*x[1]*x[1]); }

FLOAT r_p3_12(x)
FLOAT x[DIM];
{ return((x[0] + x[1])*(1. - x[0] - x[1])*x[1]); }

FLOAT r_p3_10(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*(1. - x[0] - x[1])*x[1]); }

FLOAT r_p3_20(x)
FLOAT x[DIM];
{ return((1. - x[0])*(1. - x[0] - x[1])*x[0]); }

FLOAT r_p3_21(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*x[0]*x[0]); }

FLOAT r_p3_01_0(x)
FLOAT x[DIM];
{ return((1. - x[1])*x[1]); }

FLOAT r_p3_01_1(x)
FLOAT x[DIM];
{ return((1. - 2.*x[1])*x[0]); }

FLOAT r_p3_02_0(x)
FLOAT x[DIM];
{ return(x[1]*x[1]); }

FLOAT r_p3_02_1(x)
FLOAT x[DIM];
{ return(2.*x[0]*x[1]); }

FLOAT r_p3_12_0(x)
FLOAT x[DIM];
{ return(x[1] - 2.*(x[0] + x[1])*x[1]); }

FLOAT r_p3_12_1(x)
FLOAT x[DIM];
{ return(x[1] + (x[0] + x[1])*(1. - x[0] - 3.*x[1])); }

FLOAT r_p3_10_0(x)
FLOAT x[DIM];
{ return(-2.*(1. - x[0] - x[1])*x[1]); }

FLOAT r_p3_10_1(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*(1. - x[0] - 3.*x[1])); }

FLOAT r_p3_20_0(x)
FLOAT x[DIM];
{ return((1. - 2.*x[0])*(1. - 2.*x[0] - x[1]) - x[0]*x[0]); }

FLOAT r_p3_20_1(x)
FLOAT x[DIM];
{ return(-(1. - x[0])*x[0]); }

FLOAT r_p3_21_0(x)
FLOAT x[DIM];
{ return(2.*x[0]*(1. - x[0] - x[1]) - x[0]*x[0]); }

FLOAT r_p3_21_1(x)
FLOAT x[DIM];
{ return(-x[0]*x[0]); }

FLOAT r_p3_01_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p3_01_01(x)
FLOAT x[DIM];
{ return(1. - 2.*x[1]); }

FLOAT r_p3_01_11(x)
FLOAT x[DIM];
{ return(-2.*x[0]); }

FLOAT r_p3_02_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p3_02_01(x)
FLOAT x[DIM];
{ return(2.*x[1]); }

FLOAT r_p3_02_11(x)
FLOAT x[DIM];
{ return(2.*x[0]); }

FLOAT r_p3_12_00(x)
FLOAT x[DIM];
{ return(-2.*x[1]); }

FLOAT r_p3_12_01(x)
FLOAT x[DIM];
{ return(1. - 2.*(x[0] + 2.*x[1])); }

FLOAT r_p3_12_11(x)
FLOAT x[DIM];
{ return(2. - 4.*x[0] - 6.*x[1]); }

FLOAT r_p3_10_00(x)
FLOAT x[DIM];
{ return(2.*x[1]); }

FLOAT r_p3_10_01(x)
FLOAT x[DIM];
{ return(-2. + 2.*x[0] + 4.*x[1]); }

FLOAT r_p3_10_11(x)
FLOAT x[DIM];
{ return(-4. + 4.*x[0] + 6.*x[1]); }

FLOAT r_p3_20_00(x)
FLOAT x[DIM];
{ return(-4. + 6.*x[0] + 2.*x[1]); }

FLOAT r_p3_20_01(x)
FLOAT x[DIM];
{ return(-1. + 2.*x[0]); }

FLOAT r_p3_20_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_p3_21_00(x)
FLOAT x[DIM];
{ return(2. - 6.*x[0] - 2.*x[1]); }

FLOAT r_p3_21_01(x)
FLOAT x[DIM];
{ return(-2.*x[0]); }

FLOAT r_p3_21_11(x)
FLOAT x[DIM];
{ return(0.); }

/* some p4 basis functions (of the type li*l0*l1*l2) */

FLOAT r_0012(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*(1. - x[0] - x[1])*x[0]*x[1]); }

FLOAT r_1012(x)
FLOAT x[DIM];
{ return(x[0]*(1. - x[0] - x[1])*x[0]*x[1]); }

FLOAT r_2012(x)
FLOAT x[DIM];
{ return(x[1]*(1. - x[0] - x[1])*x[0]*x[1]); }

FLOAT r_0012_0(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*(1. - 3.*x[0] - x[1])*x[1]); }

FLOAT r_0012_1(x)
FLOAT x[DIM];
{ return((1. - x[0] - x[1])*(1. - x[0] - 3.*x[1])*x[0]); }

FLOAT r_1012_0(x)
FLOAT x[DIM];
{ return((2. - 3.*x[0] - 2.*x[1])*x[0]*x[1]); }

FLOAT r_1012_1(x)
FLOAT x[DIM];
{ return((1. - x[0] - 2.*x[1])*x[0]*x[0]); }

FLOAT r_2012_0(x)
FLOAT x[DIM];
{ return((1. - 2.*x[0] - x[1])*x[1]*x[1]); }

FLOAT r_2012_1(x)
FLOAT x[DIM];
{ return((2. - 2.*x[0] - 3.*x[1])*x[0]*x[1]); }

FLOAT r_0012_00(x)
FLOAT x[DIM];
{ return(-(4. - 6.*x[0] - 4.*x[1])*x[1]); }

FLOAT r_0012_01(x)
FLOAT x[DIM];
{ return((1. - 3.*x[0] - x[1])*(1. - x[0] - 3.*x[1]) - 2.*x[0]*x[1]); }

FLOAT r_0012_11(x)
FLOAT x[DIM];
{ return(-(4. - 4.*x[0] - 6.*x[1])*x[0]); }

FLOAT r_1012_00(x)
FLOAT x[DIM];
{ return((2. - 6.*x[0] - 2.*x[1])*x[1]); }

FLOAT r_1012_01(x)
FLOAT x[DIM];
{ return((2. - 3.*x[0] - 4.*x[1])*x[0]); }

FLOAT r_1012_11(x)
FLOAT x[DIM];
{ return(-2.*x[0]*x[0]); }

FLOAT r_2012_00(x)
FLOAT x[DIM];
{ return(-2.*x[1]*x[1]); }

FLOAT r_2012_01(x)
FLOAT x[DIM];
{ return((2. - 4.*x[0] - 3.*x[1])*x[1]); }

FLOAT r_2012_11(x)
FLOAT x[DIM];
{ return((2. - 2.*x[0] - 6.*x[1])*x[0]); }

/* q1 basis funtions */

FLOAT r_q1_0(x)
FLOAT x[DIM];
{ return(0.25*(1. - x[0])*(1. - x[1])); }
/* { return((1. - x[0])*(1. - x[1])); } */

FLOAT r_q1_1(x)
FLOAT x[DIM];
{ return(0.25*(1. + x[0])*(1. - x[1])); }
/* { return(x[0]*(1. - x[1])); } */

FLOAT r_q1_2(x)
FLOAT x[DIM];
{ return(0.25*(1. + x[0])*(1. + x[1])); }
/* { return(x[0]*x[1]); } */

FLOAT r_q1_3(x)
FLOAT x[DIM];
{ return(0.25*(1. - x[0])*(1. + x[1])); }
/* { return((1. - x[0])*x[1]); } */

FLOAT r_q1_0_0(x)
FLOAT x[DIM];
{ return(0.25*(x[1] - 1.)); }
/* { return(x[1] - 1.); } */

FLOAT r_q1_0_1(x)
FLOAT x[DIM];
{ return(0.25*(x[0] - 1.)); }
/* { return(x[0] - 1.); } */

FLOAT r_q1_1_0(x)
FLOAT x[DIM];
{ return(-0.25*(x[1] - 1.)); }
/* { return(1. - x[1]); } */

FLOAT r_q1_1_1(x)
FLOAT x[DIM];
{ return(-0.25*(x[0] + 1.)); }
/* { return(-x[0]); } */

FLOAT r_q1_2_0(x)
FLOAT x[DIM];
{ return(0.25*(x[1] + 1.)); }
/* { return(x[1]); } */

FLOAT r_q1_2_1(x)
FLOAT x[DIM];
{ return(0.25*(x[0] + 1.)); }
/* { return(x[0]); } */

FLOAT r_q1_3_0(x)
FLOAT x[DIM];
{ return(-0.25*(x[1] + 1.)); }
/* { return(-x[1]); } */

FLOAT r_q1_3_1(x)
FLOAT x[DIM];
{ return(-0.25*(x[0] - 1.)); }
/* { return(1. - x[0]); } */

FLOAT r_q1_0_01(x)
FLOAT x[DIM];
{ return(0.25); }

FLOAT r_q1_1_01(x)
FLOAT x[DIM];
{ return(-0.25); }

FLOAT r_q1_2_01(x)
FLOAT x[DIM];
{ return(0.25); }

FLOAT r_q1_3_01(x)
FLOAT x[DIM];
{ return(-0.25); }

FLOAT rd_q1_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
//{ FLOAT x[DIM]={-1.,-1.}; return(v(x)); }
{ FLOAT x[DIM]={-1.,-1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
//{ FLOAT x[DIM]={ 1.,-1.}; return(v(x)); }
{ FLOAT x[DIM]={ 1.,-1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
//{ FLOAT x[DIM]={ 1., 1.}; return(v(x)); }
{ FLOAT x[DIM]={ 1., 1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1_3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
//{ FLOAT x[DIM]={-1., 1.}; return(v(x)); }
{ FLOAT x[DIM]={-1., 1.}; return(fcn_ref_map_value(x,ref_map,v)); }

/* q1x4 basis funtions */

FLOAT r_q1x4_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(x[0]*x[1]); else return(0.); }

FLOAT r_q1x4_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-x[0]*x[1]); else return(0.); }

FLOAT r_q1x4_2(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(x[0]*x[1]); else return(0.); }

FLOAT r_q1x4_3(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-x[0]*x[1]); else return(0.); }

FLOAT r_q1x4_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-(1.+x[0])*x[1]); else return((x[0]-1.)*x[1]); }

FLOAT r_q1x4_12(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(x[0]*(1.+x[1])); else return(x[0]*(1.-x[1])); }

FLOAT r_q1x4_23(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return((1.+x[0])*x[1]); else return((1.-x[0])*x[1]); }

FLOAT r_q1x4_30(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-x[0]*(1.+x[1])); else return(x[0]*(x[1]-1.)); }

FLOAT r_q1x4_b(x)
FLOAT x[DIM];
{ if (x[1] < 0.){ if (x[0] < 0.) return((1.+x[0])*(1.+x[1]));
                  else           return((1.-x[0])*(1.+x[1])); }
  else{           if (x[0] > 0.) return((1.-x[0])*(1.-x[1]));
                  else           return((1.+x[0])*(1.-x[1])); } }

FLOAT r_q1x4_0_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(x[1]); else return(0.); }

FLOAT r_q1x4_0_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(x[0]); else return(0.); }

FLOAT r_q1x4_1_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-x[1]); else return(0.); }

FLOAT r_q1x4_1_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-x[0]); else return(0.); }

FLOAT r_q1x4_2_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(x[1]); else return(0.); }

FLOAT r_q1x4_2_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(x[0]); else return(0.); }

FLOAT r_q1x4_3_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-x[1]); else return(0.); }

FLOAT r_q1x4_3_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-x[0]); else return(0.); }

FLOAT r_q1x4_01_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-x[1]); else return(x[1]); }

FLOAT r_q1x4_01_1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-1.-x[0]); else return(x[0]-1.); }

FLOAT r_q1x4_12_0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(1.+x[1]); else return(1.-x[1]); }

FLOAT r_q1x4_12_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(x[0]); else return(-x[0]); }

FLOAT r_q1x4_23_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(x[1]); else return(-x[1]); }

FLOAT r_q1x4_23_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(1.+x[0]); else return(1.-x[0]); }

FLOAT r_q1x4_30_0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-1.-x[1]); else return(x[1]-1.); }

FLOAT r_q1x4_30_1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-x[0]); else return(x[0]); }

FLOAT r_q1x4_b_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.){ if (x[0] < 0.) return( 1.+x[1]);
                  else           return(-1.-x[1]); }
  else{           if (x[0] > 0.) return(-1.+x[1]);
                  else           return( 1.-x[1]); } }

FLOAT r_q1x4_b_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.){ if (x[0] < 0.) return( 1.+x[0]);
                  else           return( 1.-x[0]); }
  else{           if (x[0] > 0.) return(-1.+x[0]);
                  else           return(-1.-x[0]); } }

FLOAT r_q1x4_0_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(1.); else return(0.); }

FLOAT r_q1x4_1_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-1.); else return(0.); }

FLOAT r_q1x4_2_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(1.); else return(0.); }

FLOAT r_q1x4_3_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-1.); else return(0.); }

FLOAT r_q1x4_01_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-1.); else return(1.); }

FLOAT r_q1x4_12_01(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(1.); else return(-1.); }

FLOAT r_q1x4_23_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(1.); else return(-1.); }

FLOAT r_q1x4_30_01(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-1.); else return(1.); }

FLOAT r_q1x4_b_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.){ if (x[0] < 0.) return( 1.);
                  else           return(-1.); }
  else{           if (x[0] > 0.) return( 1.);
                  else           return(-1.); } }

FLOAT rd_q1x4_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,-1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={ 1.,-1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={ 1., 1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1., 1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_01(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.,-1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_12(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,0.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_23(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.,1.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_30(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,0.}; return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q1x4_b(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.,0.}; return(fcn_ref_map_value(x,ref_map,v)); }

/* q1cb basis funtions */

FLOAT rd_q1cb_b(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,-1.}, x1[DIM]={1.,-1.}, x2[DIM]={1.,1.}, x3[DIM]={-1.,1.},
      x[DIM] = {0.,0.};
  return(fcn_ref_map_value(x,ref_map,v)
  -0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v)));}

/* q2 basis funtions */

FLOAT r_q2_0(x)
FLOAT x[DIM];
{ return(0.5*(1. - x[0]*x[0])*(1. - x[1])); }

FLOAT r_q2_1(x)
FLOAT x[DIM];
{ return(0.5*(1. + x[0])*(1. - x[1]*x[1])); }

FLOAT r_q2_2(x)
FLOAT x[DIM];
{ return(0.5*(1. - x[0]*x[0])*(1. + x[1])); }

FLOAT r_q2_3(x)
FLOAT x[DIM];
{ return(0.5*(1. - x[0])*(1. - x[1]*x[1])); }

FLOAT r_q2_b(x)
FLOAT x[DIM];
{ return((1. - x[0]*x[0])*(1. - x[1]*x[1])); }

FLOAT r_q2_0_0(x)
FLOAT x[DIM];
{ return(-x[0]*(1. - x[1])); }

FLOAT r_q2_0_1(x)
FLOAT x[DIM];
{ return(-0.5*(1. - x[0]*x[0])); }

FLOAT r_q2_1_0(x)
FLOAT x[DIM];
{ return(0.5*(1. - x[1]*x[1])); }

FLOAT r_q2_1_1(x)
FLOAT x[DIM];
{ return(-(1. + x[0])*x[1]); }

FLOAT r_q2_2_0(x)
FLOAT x[DIM];
{ return(-x[0]*(1. + x[1])); }

FLOAT r_q2_2_1(x)
FLOAT x[DIM];
{ return(0.5*(1. - x[0]*x[0])); }

FLOAT r_q2_3_0(x)
FLOAT x[DIM];
{ return(-0.5*(1. - x[1]*x[1])); }

FLOAT r_q2_3_1(x)
FLOAT x[DIM];
{ return(-(1. - x[0])*x[1]); }

FLOAT r_q2_b_0(x)
FLOAT x[DIM];
{ return(-2.*x[0]*(1. - x[1]*x[1])); }

FLOAT r_q2_b_1(x)
FLOAT x[DIM];
{ return(-2.*(1. - x[0]*x[0])*x[1]); }

FLOAT r_q2_0_00(x)
FLOAT x[DIM];
{ return(-1. + x[1]); }

FLOAT r_q2_0_01(x)
FLOAT x[DIM];
{ return(x[0]); }

FLOAT r_q2_0_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2_1_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2_1_01(x)
FLOAT x[DIM];
{ return(-x[1]); }

FLOAT r_q2_1_11(x)
FLOAT x[DIM];
{ return(-1. - x[0]); }

FLOAT r_q2_2_00(x)
FLOAT x[DIM];
{ return(-1. - x[1]); }

FLOAT r_q2_2_01(x)
FLOAT x[DIM];
{ return(-x[0]); }

FLOAT r_q2_2_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2_3_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2_3_01(x)
FLOAT x[DIM];
{ return(x[1]); }

FLOAT r_q2_3_11(x)
FLOAT x[DIM];
{ return(-1. + x[0]); }

FLOAT r_q2_b_00(x)
FLOAT x[DIM];
{ return(-2.*(1. - x[1]*x[1])); }

FLOAT r_q2_b_01(x)
FLOAT x[DIM];
{ return(4.*x[0]*x[1]); }

FLOAT r_q2_b_11(x)
FLOAT x[DIM];
{ return(-2.*(1. - x[0]*x[0])); }

FLOAT rd_q2_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,-1.}, x1[DIM]={ 1.,-1.}, x[DIM]={ 0.,-1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v))); }

FLOAT rd_q2_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x1[DIM]={ 1.,-1.}, x2[DIM]={ 1., 1.}, x[DIM]={ 1., 0.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x1,ref_map,v) + fcn_ref_map_value(x2,ref_map,v))); }

FLOAT rd_q2_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x2[DIM]={ 1., 1.}, x3[DIM]={-1., 1.}, x[DIM]={ 0., 1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))); }

FLOAT rd_q2_3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x3[DIM]={-1., 1.}, x0[DIM]={-1.,-1.}, x[DIM]={-1., 0.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x3,ref_map,v) + fcn_ref_map_value(x0,ref_map,v))); }

FLOAT rd_q2_b(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,-1.}, x1[DIM]={1.,-1.}, x2[DIM]={1.,1.}, x3[DIM]={-1.,1.},
      x01[DIM]={0.,-1.}, x12[DIM]={1.,0.}, x23[DIM]={0.,1.}, x30[DIM]={-1.,0.},
      x[DIM] = {0.,0.};
  return(fcn_ref_map_value(x,ref_map,v)
  +0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x12,ref_map,v)
       + fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x30,ref_map,v)));}

/* q2x4 basis funtions */

FLOAT r_q2x4_00(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*x[0]*(1.+x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*x[0]*(1.-x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_11(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*x[0]*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_12(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*x[0]*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_22(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*x[0]*(1.-x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_23(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*x[0]*(1.+x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_33(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*x[0]*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_30(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*x[0]*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-4.*(1.+x[0])*x[1]*(1.+x[1])); 
  else                return(-4.*(1.-x[0])*x[1]*(1.+x[1])); }

FLOAT r_q2x4_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(4.*x[0]*(1.-x[0])*(1.+x[1])); 
  else                return(4.*x[0]*(1.-x[0])*(1.-x[1])); }

FLOAT r_q2x4_2(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(4.*(1.+x[0])*x[1]*(1.-x[1])); 
  else                return(4.*(1.-x[0])*x[1]*(1.-x[1])); }

FLOAT r_q2x4_3(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-4.*x[0]*(1.+x[0])*(1.+x[1])); 
  else                return(-4.*x[0]*(1.+x[0])*(1.-x[1])); }

FLOAT r_q2x4_b0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(16.*x[0]*x[1]*(1.+x[0])*(1.+x[1])); 
  else return(0.); }

FLOAT r_q2x4_b1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-16.*x[0]*x[1]*(1.-x[0])*(1.+x[1]));
  else return(0.); }

FLOAT r_q2x4_b2(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(16.*x[0]*x[1]*(1.-x[0])*(1.-x[1]));
  else return(0.); }

FLOAT r_q2x4_b3(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-16.*x[0]*x[1]*(1.+x[0])*(1.-x[1]));
  else return(0.); }

FLOAT r_q2x4_00_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*(1.+2.*x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_00_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*x[0]*(1.+x[0])); else return(0.); }

FLOAT r_q2x4_01_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*(1.-2.*x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_01_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*x[0]*(1.-x[0])); else return(0.); }

FLOAT r_q2x4_11_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_11_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*x[0]*(1.+2.*x[1])); else return(0.); }

FLOAT r_q2x4_12_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_12_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*x[0]*(1.-2.*x[1])); else return(0.); }

FLOAT r_q2x4_22_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*(1.-2.*x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_22_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*x[0]*(1.-x[0])); else return(0.); }

FLOAT r_q2x4_23_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*(1.+2.*x[0])*x[1]); else return(0.); }

FLOAT r_q2x4_23_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*x[0]*(1.+x[0])); else return(0.); }

FLOAT r_q2x4_33_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_33_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*x[0]*(1.-2.*x[1])); else return(0.); }

FLOAT r_q2x4_30_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_30_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*x[0]*(1.+2.*x[1])); else return(0.); }

FLOAT r_q2x4_0_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-4.*x[1]*(1.+x[1])); 
  else                return( 4.*x[1]*(1.+x[1])); }

FLOAT r_q2x4_0_1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-4.*(1.+x[0])*(1.+2.*x[1])); 
  else                return(-4.*(1.-x[0])*(1.+2.*x[1])); }

FLOAT r_q2x4_1_0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(4.*(1.-2.*x[0])*(1.+x[1])); 
  else                return(4.*(1.-2.*x[0])*(1.-x[1])); }

FLOAT r_q2x4_1_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return( 4.*x[0]*(1.-x[0])); 
  else                return(-4.*x[0]*(1.-x[0])); }

FLOAT r_q2x4_2_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return( 4.*x[1]*(1.-x[1])); 
  else                return(-4.*x[1]*(1.-x[1])); }

FLOAT r_q2x4_2_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(4.*(1.+x[0])*(1.-2.*x[1])); 
  else                return(4.*(1.-x[0])*(1.-2.*x[1])); }

FLOAT r_q2x4_3_0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-4.*(1.+2.*x[0])*(1.+x[1])); 
  else                return(-4.*(1.+2.*x[0])*(1.-x[1])); }

FLOAT r_q2x4_3_1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-4.*x[0]*(1.+x[0])); 
  else                return( 4.*x[0]*(1.+x[0])); }

FLOAT r_q2x4_b0_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(16.*x[1]*(1.+2.*x[0])*(1.+x[1])); 
  else return(0.); }

FLOAT r_q2x4_b0_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(16.*x[0]*(1.+x[0])*(1.+2.*x[1])); 
  else return(0.); }

FLOAT r_q2x4_b1_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-16.*x[1]*(1.-2.*x[0])*(1.+x[1]));
  else return(0.); }

FLOAT r_q2x4_b1_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-16.*x[0]*(1.-x[0])*(1.+2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_b2_0(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(16.*x[1]*(1.-2.*x[0])*(1.-x[1]));
  else return(0.); }

FLOAT r_q2x4_b2_1(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(16.*x[0]*(1.-x[0])*(1.-2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_b3_0(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-16.*x[1]*(1.+2.*x[0])*(1.-x[1]));
  else return(0.); }

FLOAT r_q2x4_b3_1(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-16.*x[0]*(1.+x[0])*(1.-2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_00_00(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(8.*x[1]); else return(0.); }

FLOAT r_q2x4_00_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*(1.+2.*x[0])); else return(0.); }

FLOAT r_q2x4_00_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_01_00(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(8.*x[1]); else return(0.); }

FLOAT r_q2x4_01_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*(1.-2.*x[0])); else return(0.); }

FLOAT r_q2x4_01_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_11_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_11_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-4.*(1.+2.*x[1])); else return(0.); }

FLOAT r_q2x4_11_11(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-8.*x[0]); else return(0.); }

FLOAT r_q2x4_12_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_12_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*(1.-2.*x[1])); else return(0.); }

FLOAT r_q2x4_12_11(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(-8.*x[0]); else return(0.); }

FLOAT r_q2x4_22_00(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(-8.*x[1]); else return(0.); }

FLOAT r_q2x4_22_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(4.*(1.-2.*x[0])); else return(0.); }

FLOAT r_q2x4_22_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_23_00(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-8.*x[1]); else return(0.); }

FLOAT r_q2x4_23_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*(1.+2.*x[0])); else return(0.); }

FLOAT r_q2x4_23_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_33_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_33_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-4.*(1.-2.*x[1])); else return(0.); }

FLOAT r_q2x4_33_11(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(8.*x[0]); else return(0.); }

FLOAT r_q2x4_30_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_30_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(4.*(1.+2.*x[1])); else return(0.); }

FLOAT r_q2x4_30_11(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(8.*x[0]); else return(0.); }

FLOAT r_q2x4_0_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_0_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-4.*(1.+2.*x[1])); 
  else                return( 4.*(1.+2.*x[1])); }

FLOAT r_q2x4_0_11(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(0.);
  else if (x[0] < 0.) return(-8.*(1.+x[0])); 
  else                return(-8.*(1.-x[0])); }

FLOAT r_q2x4_1_00(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return(-8.*(1.+x[1])); 
  else                return(-8.*(1.-x[1])); }

FLOAT r_q2x4_1_01(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(0.);
  else if (x[1] < 0.) return( 4.*(1.-2.*x[0])); 
  else                return(-4.*(1.-2.*x[0])); }

FLOAT r_q2x4_1_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_2_00(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_2_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return( 4.*(1.-2.*x[1])); 
  else                return(-4.*(1.-2.*x[1])); }

FLOAT r_q2x4_2_11(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(0.);
  else if (x[0] < 0.) return(-8.*(1.+x[0])); 
  else                return(-8.*(1.-x[0])); }

FLOAT r_q2x4_3_00(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-8.*(1.+x[1])); 
  else                return(-8.*(1.-x[1])); }

FLOAT r_q2x4_3_01(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(0.);
  else if (x[1] < 0.) return(-4.*(1.+2.*x[0])); 
  else                return( 4.*(1.+2.*x[0])); }

FLOAT r_q2x4_3_11(x)
FLOAT x[DIM];
{ return(0.); }

FLOAT r_q2x4_b0_00(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(32.*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_b0_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(16.*(1.+2.*x[0])*(1.+2.*x[1])); 
  else return(0.); }

FLOAT r_q2x4_b0_11(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] < 0.) return(32.*x[0]*(1.+x[0])); else return(0.); }

FLOAT r_q2x4_b1_00(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(32.*x[1]*(1.+x[1])); else return(0.); }

FLOAT r_q2x4_b1_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-16.*(1.-2.*x[0])*(1.+2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_b1_11(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] < 0.) return(-32.*x[0]*(1.-x[0])); else return(0.); }

FLOAT r_q2x4_b2_00(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(-32.*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_b2_01(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(16.*(1.-2.*x[0])*(1.-2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_b2_11(x)
FLOAT x[DIM];
{ if (x[0] > 0. && x[1] > 0.) return(-32.*x[0]*(1.-x[0])); else return(0.); }

FLOAT r_q2x4_b3_00(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-32.*x[1]*(1.-x[1])); else return(0.); }

FLOAT r_q2x4_b3_01(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(-16.*(1.+2.*x[0])*(1.-2.*x[1]));
  else return(0.); }

FLOAT r_q2x4_b3_11(x)
FLOAT x[DIM];
{ if (x[0] < 0. && x[1] > 0.) return(32.*x[0]*(1.+x[0]));
  else return(0.); }

FLOAT rd_q2x4_00(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,-1.}, x01[DIM]={ 0.,-1.}, x[DIM]={-0.5,-1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x01,ref_map,v))); }

FLOAT rd_q2x4_01(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x01[DIM]={ 0.,-1.}, x1[DIM]={ 1.,-1.}, x[DIM]={ 0.5,-1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x1,ref_map,v))); }

FLOAT rd_q2x4_11(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x1[DIM]={ 1.,-1.}, x12[DIM]={ 1., 0.}, x[DIM]={ 1.,-0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x1,ref_map,v) + fcn_ref_map_value(x12,ref_map,v))); }

FLOAT rd_q2x4_12(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x12[DIM]={ 1., 0.}, x2[DIM]={ 1., 1.}, x[DIM]={ 1., 0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x12,ref_map,v) + fcn_ref_map_value(x2,ref_map,v))); }

FLOAT rd_q2x4_22(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x2[DIM]={ 1., 1.}, x23[DIM]={ 0., 1.}, x[DIM]={ 0.5, 1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x23,ref_map,v))); }

FLOAT rd_q2x4_23(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x23[DIM]={ 0., 1.}, x3[DIM]={-1., 1.}, x[DIM]={-0.5, 1.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))); }

FLOAT rd_q2x4_33(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x3[DIM]={-1., 1.}, x30[DIM]={-1., 0.}, x[DIM]={-1., 0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x3,ref_map,v) + fcn_ref_map_value(x30,ref_map,v))); }

FLOAT rd_q2x4_30(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x30[DIM]={-1., 0.}, x0[DIM]={-1.,-1.}, x[DIM]={-1.,-0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x30,ref_map,v) + fcn_ref_map_value(x0,ref_map,v))); }

FLOAT rd_q2x4_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x01[DIM]={ 0.,-1.}, xc[DIM]={ 0., 0.}, x[DIM]={ 0.,-0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(xc,ref_map,v))); }

FLOAT rd_q2x4_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x12[DIM]={ 1., 0.}, xc[DIM]={ 0., 0.}, x[DIM]={ 0.5, 0.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x12,ref_map,v) + fcn_ref_map_value(xc,ref_map,v))); }

FLOAT rd_q2x4_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x23[DIM]={ 0., 1.}, xc[DIM]={ 0., 0.}, x[DIM]={ 0., 0.5};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(xc,ref_map,v))); }

FLOAT rd_q2x4_3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x30[DIM]={-1., 0.}, xc[DIM]={ 0., 0.}, x[DIM]={-0.5, 0.};
  return(fcn_ref_map_value(x,ref_map,v)
   -0.5*(fcn_ref_map_value(x30,ref_map,v) + fcn_ref_map_value(xc,ref_map,v))); }

FLOAT rd_q2x4_b0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,-1.}, x1[DIM]={0.,-1.}, x2[DIM]={0.,0.}, x3[DIM]={-1.,0.},
        x01[DIM]={-0.5,-1.}, x12[DIM]={ 0.,-0.5}, 
        x23[DIM]={-0.5, 0.}, x30[DIM]={-1.,-0.5}, x[DIM] = {-0.5,-0.5};
  return(fcn_ref_map_value(x,ref_map,v)
  +0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x12,ref_map,v)
       + fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x30,ref_map,v)));}

FLOAT rd_q2x4_b1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={1.,-1.}, x1[DIM]={1.,0.}, x2[DIM]={0.,0.}, x3[DIM]={0.,-1.},
        x01[DIM]={ 1.,-0.5}, x12[DIM]={ 0.5, 0.},
        x23[DIM]={ 0.,-0.5}, x30[DIM]={ 0.5,-1.}, x[DIM] = { 0.5,-0.5};
  return(fcn_ref_map_value(x,ref_map,v)
  +0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x12,ref_map,v)
       + fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x30,ref_map,v)));}

FLOAT rd_q2x4_b2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={1.,1.}, x1[DIM]={0.,1.}, x2[DIM]={0.,0.}, x3[DIM]={1.,0.},
        x01[DIM]={ 0.5, 1.}, x12[DIM]={ 0., 0.5}, 
        x23[DIM]={ 0.5, 0.}, x30[DIM]={ 1., 0.5}, x[DIM] = { 0.5, 0.5};
  return(fcn_ref_map_value(x,ref_map,v)
  +0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x12,ref_map,v)
       + fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x30,ref_map,v)));}

FLOAT rd_q2x4_b3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x0[DIM]={-1.,1.}, x1[DIM]={-1.,0.}, x2[DIM]={0.,0.}, x3[DIM]={0.,1.},
        x01[DIM]={-1., 0.5}, x12[DIM]={-0.5, 0.},
        x23[DIM]={ 0., 0.5}, x30[DIM]={-0.5, 1.}, x[DIM] = {-0.5, 0.5};
  return(fcn_ref_map_value(x,ref_map,v)
  +0.25*(fcn_ref_map_value(x0,ref_map,v) + fcn_ref_map_value(x1,ref_map,v)
       + fcn_ref_map_value(x2,ref_map,v) + fcn_ref_map_value(x3,ref_map,v))
   -0.5*(fcn_ref_map_value(x01,ref_map,v) + fcn_ref_map_value(x12,ref_map,v)
       + fcn_ref_map_value(x23,ref_map,v) + fcn_ref_map_value(x30,ref_map,v)));}

/* q2nb basis funtions */

FLOAT r_q2nb0(x)
FLOAT x[DIM];
{ return((x[0] - x[0]*x[0]*x[0])*(1. - x[1]*x[1])); }

FLOAT r_q2nb0_0(x)
FLOAT x[DIM];
{ return((1. - 3.*x[0]*x[0])*(1. - x[1]*x[1])); }

FLOAT r_q2nb0_1(x)
FLOAT x[DIM];
{ return(-2.*(x[0] - x[0]*x[0]*x[0])*x[1]); }

FLOAT r_q2nb0_00(x)
FLOAT x[DIM];
{ return(-6.*x[0]*(1. - x[1]*x[1])); }

FLOAT r_q2nb0_01(x)
FLOAT x[DIM];
{ return(-2.*(1. - 3.*x[0]*x[0])*x[1]); }

FLOAT r_q2nb0_11(x)
FLOAT x[DIM];
{ return(-2.*(x[0] - x[0]*x[0]*x[0])); }

FLOAT rd_q2nb0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2nb0 not available.\n"); return(0.);  }

FLOAT r_q2nb1(x)
FLOAT x[DIM];
{ return((1. - x[0]*x[0])*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb1_0(x)
FLOAT x[DIM];
{ return(-2.*x[0]*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb1_1(x)
FLOAT x[DIM];
{ return((1. - x[0]*x[0])*(1. - 3.*x[1]*x[1])); }

FLOAT r_q2nb1_00(x)
FLOAT x[DIM];
{ return(-2.*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb1_01(x)
FLOAT x[DIM];
{ return(-2.*x[0]*(1. - 3.*x[1]*x[1])); }

FLOAT r_q2nb1_11(x)
FLOAT x[DIM];
{ return(-6.*(1. - x[0]*x[0])*x[1]); }

FLOAT rd_q2nb1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2nb1 not available.\n"); return(0.);  }

FLOAT r_q2nb01(x)
FLOAT x[DIM];
{ return((x[0] - x[0]*x[0]*x[0])*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb01_0(x)
FLOAT x[DIM];
{ return((1. - 3.*x[0]*x[0])*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb01_1(x)
FLOAT x[DIM];
{ return((x[0] - x[0]*x[0]*x[0])*(1. - 3.*x[1]*x[1])); }

FLOAT r_q2nb01_00(x)
FLOAT x[DIM];
{ return(-6.*x[0]*(x[1] - x[1]*x[1]*x[1])); }

FLOAT r_q2nb01_01(x)
FLOAT x[DIM];
{ return((1. - 3.*x[0]*x[0])*(1. - 3.*x[1]*x[1])); }

FLOAT r_q2nb01_11(x)
FLOAT x[DIM];
{ return(-6.*(x[0] - x[0]*x[0]*x[0])*x[1]); }

FLOAT rd_q2nb01(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2nb01 not available.\n"); return(0.);  }

/* q2sbt basis funtions */

FLOAT r_q2sbt(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - x[0]*x[0])*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*x[0]*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt_1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - x[0]*x[0])*(1. - 2.*x[1]));  else return(0.); }

FLOAT r_q2sbt_00(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*x[0]*(1. - 2.*x[1]));  else return(0.); }

FLOAT r_q2sbt_11(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*(1. - x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbt(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbt not available.\n"); return(0.); }

FLOAT r_q2sbt0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((x[0] - x[0]*x[0]*x[0])*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt0_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - 3.*x[0]*x[0])*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt0_1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((x[0] - x[0]*x[0]*x[0])*(1. - 2.*x[1]));  else return(0.); }

FLOAT r_q2sbt0_00(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-6.*x[0]*(x[1] - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt0_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - 3.*x[0]*x[0])*(1. - 2.*x[1]));  else return(0.); }

FLOAT r_q2sbt0_11(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*(x[0] - x[0]*x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbt0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbt0 not available.\n"); return(0.); }

FLOAT r_q2sbt1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - x[0]*x[0])*(x[1]*x[1] - x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt1_0(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*x[0]*(x[1]*x[1] - x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt1_1(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - x[0]*x[0])*(2.*x[1] - 3.*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt1_00(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*(x[1]*x[1] - x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt1_01(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return(-2.*x[0]*(2.*x[1] - 3.*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbt1_11(x)
FLOAT x[DIM];
{ if (x[1] > 0.) return((1. - x[0]*x[0])*(2. - 6.*x[1]));  else return(0.); }

FLOAT rd_q2sbt1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbt1 not available.\n"); return(0.); }

/* q2sbb basis funtions */

FLOAT r_q2sbb(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - x[0]*x[0])*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*x[0]*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - x[0]*x[0])*(1. + 2.*x[1]));  else return(0.); }

FLOAT r_q2sbb_00(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*x[0]*(1. + 2.*x[1]));  else return(0.); }

FLOAT r_q2sbb_11(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(2.*(1. - x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbb(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbb not available.\n"); return(0.); }

FLOAT r_q2sbb0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((x[0] - x[0]*x[0]*x[0])*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb0_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - 3.*x[0]*x[0])*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb0_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((x[0] - x[0]*x[0]*x[0])*(1. + 2.*x[1]));  else return(0.); }

FLOAT r_q2sbb0_00(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-6.*x[0]*(x[1] + x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb0_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - 3.*x[0]*x[0])*(1. + 2.*x[1]));  else return(0.); }

FLOAT r_q2sbb0_11(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(2.*(x[0] - x[0]*x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbb0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbb0 not available.\n"); return(0.); }

FLOAT r_q2sbb1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - x[0]*x[0])*(x[1]*x[1] + x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb1_0(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*x[0]*(x[1]*x[1] + x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb1_1(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - x[0]*x[0])*(2.*x[1] + 3.*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb1_00(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*(x[1]*x[1] + x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb1_01(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return(-2.*x[0]*(2.*x[1] + 3.*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbb1_11(x)
FLOAT x[DIM];
{ if (x[1] < 0.) return((1. - x[0]*x[0])*(2. + 6.*x[1]));  else return(0.); }

FLOAT rd_q2sbb1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbb1 not available.\n"); return(0.); }

/* q2sbr basis funtions */

FLOAT r_q2sbr(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - x[1]*x[1])*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr_0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - x[1]*x[1])*(1. - 2.*x[0]));  else return(0.); }

FLOAT r_q2sbr_1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*x[1]*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr_00(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*(1. - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbr_01(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*x[1]*(1. - 2.*x[0]));  else return(0.); }

FLOAT r_q2sbr_11(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbr(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbr not available.\n"); return(0.); }

FLOAT r_q2sbr0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - x[1]*x[1])*(x[0]*x[0] - x[0]*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr0_0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - x[1]*x[1])*(2.*x[0] - 3.*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr0_1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*x[1]*(x[0]*x[0] - x[0]*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr0_00(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - x[1]*x[1])*(2. - 6.*x[0]));  else return(0.); }

FLOAT r_q2sbr0_01(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*x[1]*(2.*x[0] - 3.*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr0_11(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*(x[0]*x[0] - x[0]*x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbr0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbr0 not available.\n"); return(0.); }

FLOAT r_q2sbr1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((x[1] - x[1]*x[1]*x[1])*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr1_0(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((x[1] - x[1]*x[1]*x[1])*(1. - 2.*x[0]));  else return(0.); }

FLOAT r_q2sbr1_1(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - 3.*x[1]*x[1])*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbr1_00(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-2.*(x[1] - x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbr1_01(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return((1. - 3.*x[1]*x[1])*(1. - 2.*x[0]));  else return(0.); }

FLOAT r_q2sbr1_11(x)
FLOAT x[DIM];
{ if (x[0] > 0.) return(-6.*x[1]*(x[0] - x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbr1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbr1 not available.\n"); return(0.); }

/* q2sbl basis funtions */

FLOAT r_q2sbl(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - x[1]*x[1])*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl_0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - x[1]*x[1])*(1. + 2.*x[0]));  else return(0.); }

FLOAT r_q2sbl_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*x[1]*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl_00(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(2.*(1. - x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbl_01(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*x[1]*(1. + 2.*x[0]));  else return(0.); }

FLOAT r_q2sbl_11(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbl(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbl not available.\n"); return(0.); }

FLOAT r_q2sbl0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - x[1]*x[1])*(x[0]*x[0] + x[0]*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl0_0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - x[1]*x[1])*(2.*x[0] + 3.*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl0_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*x[1]*(x[0]*x[0] + x[0]*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl0_00(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - x[1]*x[1])*(2. + 6.*x[0]));  else return(0.); }

FLOAT r_q2sbl0_01(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*x[1]*(2.*x[0] + 3.*x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl0_11(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-2.*(x[0]*x[0] + x[0]*x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbl0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbl0 not available.\n"); return(0.); }

FLOAT r_q2sbl1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((x[1] - x[1]*x[1]*x[1])*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl1_0(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((x[1] - x[1]*x[1]*x[1])*(1. + 2.*x[0]));  else return(0.); }

FLOAT r_q2sbl1_1(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - 3.*x[1]*x[1])*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT r_q2sbl1_00(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(2.*(x[1] - x[1]*x[1]*x[1]));  else return(0.); }

FLOAT r_q2sbl1_01(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return((1. - 3.*x[1]*x[1])*(1. + 2.*x[0]));  else return(0.); }

FLOAT r_q2sbl1_11(x)
FLOAT x[DIM];
{ if (x[0] < 0.) return(-6.*x[1]*(x[0] + x[0]*x[0]));  else return(0.); }

FLOAT rd_q2sbl1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ eprintf("Error: rd_q2sbl1 not available.\n"); return(0.); }

/* q3 basis funtions */

FLOAT r_q3_0(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_1(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_2(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_3(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_4(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_5(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_6(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_7(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_8(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_9(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_10(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_12(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_13(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_14(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_15(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_0_0(x)
FLOAT x[DIM];
{ return((27.*x[0]*x[0]-18.*x[0]-1.)*
         (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_0_1(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_1_0(x)
FLOAT x[DIM];
{ return((-27.*x[0]*x[0]-18.*x[0]+1.)*
         (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_1_1(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_2_0(x)
FLOAT x[DIM];
{ return((-27.*x[0]*x[0]-18.*x[0]+1.)*
         (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_2_1(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_3_0(x)
FLOAT x[DIM];
{ return((27.*x[0]*x[0]-18.*x[0]-1.)*
         (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_3_1(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
         (-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_4_0(x)
FLOAT x[DIM];
{ return(-9.*(9.*x[0]*x[0]-2.*x[0]-3.)*
             (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_4_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_5_0(x)
FLOAT x[DIM];
{ return(-9.*(-27.*x[0]*x[0]-18.*x[0]+1.)*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_5_1(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_6_0(x)
FLOAT x[DIM];
{ return(-9.*(-9.*x[0]*x[0]-2.*x[0]+3.)*
             (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_6_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_7_0(x)
FLOAT x[DIM];
{ return(-9.*(27.*x[0]*x[0]-18.*x[0]-1.)*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_7_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_8_0(x)
FLOAT x[DIM];
{ return(-9.*(-9.*x[0]*x[0]-2.*x[0]+3.)*
             (1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_8_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_9_0(x)
FLOAT x[DIM];
{ return(-9.*(-27.*x[0]*x[0]-18.*x[0]+1.)*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_9_1(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_10_0(x)
FLOAT x[DIM];
{ return(-9.*(9.*x[0]*x[0]-2.*x[0]-3.)*
             (1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_10_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_11_0(x)
FLOAT x[DIM];
{ return(-9.*(27.*x[0]*x[0]-18.*x[0]-1.)*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_11_1(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*
             (9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_12_0(x)
FLOAT x[DIM];
{ return(81.*(9.*x[0]*x[0]-2.*x[0]-3.)*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_12_1(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_13_0(x)
FLOAT x[DIM];
{ return(81.*(-9.*x[0]*x[0]-2.*x[0]+3.)*
             (1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_13_1(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_14_0(x)
FLOAT x[DIM];
{ return(81.*(-9.*x[0]*x[0]-2.*x[0]+3.)*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_14_1(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*
             (-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_15_0(x)
FLOAT x[DIM];
{ return(81.*(9.*x[0]*x[0]-2.*x[0]-3.)*
             (1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_15_1(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*
             (-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_0_00(x)
FLOAT x[DIM];
{ return((54.*x[0]-18.)*(1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_0_01(x)
FLOAT x[DIM];
{ return((27.*x[0]*x[0]-18.*x[0]-1.)*(27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_0_11(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(54.*x[1]-18.)/256.); }

FLOAT r_q3_1_00(x)
FLOAT x[DIM];
{ return((-54.*x[0]-18.)*(1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_1_01(x)
FLOAT x[DIM];
{ return((-27.*x[0]*x[0]-18.*x[0]+1.)*(27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_1_11(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(54.*x[1]-18.)/256.); }

FLOAT r_q3_2_00(x)
FLOAT x[DIM];
{ return((-54.*x[0]-18.)*(1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_2_01(x)
FLOAT x[DIM];
{ return((-27.*x[0]*x[0]-18.*x[0]+1.)*(-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_2_11(x)
FLOAT x[DIM];
{ return((1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(-54.*x[1]-18.)/256.); }

FLOAT r_q3_3_00(x)
FLOAT x[DIM];
{ return((54.*x[0]-18.)*(1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_3_01(x)
FLOAT x[DIM];
{ return((27.*x[0]*x[0]-18.*x[0]-1.)*(-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_3_11(x)
FLOAT x[DIM];
{ return((1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(-54.*x[1]-18.)/256.); }

FLOAT r_q3_4_00(x)
FLOAT x[DIM];
{ return(-9.*(18.*x[0]-2.)*(1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_4_01(x)
FLOAT x[DIM];
{ return(-9.*(9.*x[0]*x[0]-2.*x[0]-3.)*(27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_4_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*(54.*x[1]-18.)/256.); }

FLOAT r_q3_5_00(x)
FLOAT x[DIM];
{ return(-9.*(-54.*x[0]-18.)*(1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_5_01(x)
FLOAT x[DIM];
{ return(-9.*(-27.*x[0]*x[0]-18.*x[0]+1.)*(9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_5_11(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(18.*x[1]-2.)/256.); }

FLOAT r_q3_6_00(x)
FLOAT x[DIM];
{ return(-9.*(-18.*x[0]-2.)*(1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_6_01(x)
FLOAT x[DIM];
{ return(-9.*(-9.*x[0]*x[0]-2.*x[0]+3.)*(-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_6_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*(-54.*x[1]-18.)/256.); }

FLOAT r_q3_7_00(x)
FLOAT x[DIM];
{ return(-9.*(54.*x[0]-18.)*(1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_7_01(x)
FLOAT x[DIM];
{ return(-9.*(27.*x[0]*x[0]-18.*x[0]-1.)*(-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_7_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(-18.*x[1]-2.)/256.); }

FLOAT r_q3_8_00(x)
FLOAT x[DIM];
{ return(-9.*(-18.*x[0]-2.)*(1.-x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_8_01(x)
FLOAT x[DIM];
{ return(-9.*(-9.*x[0]*x[0]-2.*x[0]+3.)*(27.*x[1]*x[1]-18.*x[1]-1.)/256.); }

FLOAT r_q3_8_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*(54.*x[1]-18.)/256.); }

FLOAT r_q3_9_00(x)
FLOAT x[DIM];
{ return(-9.*(-54.*x[0]-18.)*(1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_9_01(x)
FLOAT x[DIM];
{ return(-9.*(-27.*x[0]*x[0]-18.*x[0]+1.)*(-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_9_11(x)
FLOAT x[DIM];
{ return(-9.*(1.+x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(-18.*x[1]-2.)/256.); }

FLOAT r_q3_10_00(x)
FLOAT x[DIM];
{ return(-9.*(18.*x[0]-2.)*(1.+x[1])*(1.-3.*x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_10_01(x)
FLOAT x[DIM];
{ return(-9.*(9.*x[0]*x[0]-2.*x[0]-3.)*(-27.*x[1]*x[1]-18.*x[1]+1.)/256.); }

FLOAT r_q3_10_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*(-54.*x[1]-18.)/256.); }

FLOAT r_q3_11_00(x)
FLOAT x[DIM];
{ return(-9.*(54.*x[0]-18.)*(1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_11_01(x)
FLOAT x[DIM];
{ return(-9.*(27.*x[0]*x[0]-18.*x[0]-1.)*(9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_11_11(x)
FLOAT x[DIM];
{ return(-9.*(1.-x[0])*(1.-3.*x[0])*(1.+3.*x[0])*(18.*x[1]-2.)/256.); }

FLOAT r_q3_12_00(x)
FLOAT x[DIM];
{ return(81.*(18.*x[0]-2.)*(1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_12_01(x)
FLOAT x[DIM];
{ return(81.*(9.*x[0]*x[0]-2.*x[0]-3.)*(9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_12_11(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*(18.*x[1]-2.)/256.); }

FLOAT r_q3_13_00(x)
FLOAT x[DIM];
{ return(81.*(-18.*x[0]-2.)*(1.-x[1])*(1.+x[1])*(1.-3.*x[1])/256.); }

FLOAT r_q3_13_01(x)
FLOAT x[DIM];
{ return(81.*(-9.*x[0]*x[0]-2.*x[0]+3.)*(9.*x[1]*x[1]-2.*x[1]-3.)/256.); }

FLOAT r_q3_13_11(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*(18.*x[1]-2.)/256.); }

FLOAT r_q3_14_00(x)
FLOAT x[DIM];
{ return(81.*(-18.*x[0]-2.)*(1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_14_01(x)
FLOAT x[DIM];
{ return(81.*(-9.*x[0]*x[0]-2.*x[0]+3.)*(-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_14_11(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.+3.*x[0])*(-18.*x[1]-2.)/256.); }

FLOAT r_q3_15_00(x)
FLOAT x[DIM];
{ return(81.*(18.*x[0]-2.)*(1.-x[1])*(1.+x[1])*(1.+3.*x[1])/256.); }

FLOAT r_q3_15_01(x)
FLOAT x[DIM];
{ return(81.*(9.*x[0]*x[0]-2.*x[0]-3.)*(-9.*x[1]*x[1]-2.*x[1]+3.)/256.); }

FLOAT r_q3_15_11(x)
FLOAT x[DIM];
{ return(81.*(1.-x[0])*(1.+x[0])*(1.-3.*x[0])*(-18.*x[1]-2.)/256.); }

FLOAT rd_q3_0(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,-1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_1(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,-1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_2(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_3(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_4(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-0.3333333333333333333333,-1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_5(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,-0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_6(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.3333333333333333333333,1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_7(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_8(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.3333333333333333333333,-1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_9(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={1.,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_10(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-0.3333333333333333333333,1.};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_11(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-1.,-0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_12(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-0.3333333333333333333333,-0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_13(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.3333333333333333333333,-0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_14(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={0.3333333333333333333333,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

FLOAT rd_q3_15(v,ref_map)
FLOAT (*v)();
REF_MAPPING *ref_map;
{ FLOAT x[DIM]={-0.3333333333333333333333,0.3333333333333333333333};
  return(fcn_ref_map_value(x,ref_map,v)); }

/* q1rot basis funtions */

DOUBLE r_phi0(x)
DOUBLE *x;
{
   return( 0.25 - 0.5*x[1] - (Q1ROT_T(x[0])-Q1ROT_T(x[1]))/Q1ROT_K );
}

DOUBLE r_phi0_0(x)  /* x[0] derivative  */
DOUBLE *x;
{
   return(  -Q1ROT_TD(x[0])/Q1ROT_K );
}

DOUBLE r_phi0_1(x)  /* x[1] derivative  */
DOUBLE *x;
{
   return(  -0.5 + Q1ROT_TD(x[1])/Q1ROT_K );
}

DOUBLE r_phi1(x)
DOUBLE *x;
{
   return( 0.25 + 0.5*x[0] + (Q1ROT_T(x[0])-Q1ROT_T(x[1]))/Q1ROT_K );
}

DOUBLE r_phi1_0(x)
DOUBLE *x;
{
   return( 0.5 + Q1ROT_TD(x[0])/Q1ROT_K );
}

DOUBLE r_phi1_1(x)
DOUBLE *x;
{
   return( -Q1ROT_TD(x[1])/Q1ROT_K );
}

DOUBLE r_phi2(x)
DOUBLE *x;
{
   return( 0.25 + 0.5*x[1] - (Q1ROT_T(x[0])-Q1ROT_T(x[1]))/Q1ROT_K );
}

DOUBLE r_phi2_0(x)
DOUBLE *x;
{
   return(  -Q1ROT_TD(x[0])/Q1ROT_K );
}

DOUBLE r_phi2_1(x)
DOUBLE *x;
{
   return( 0.5 + Q1ROT_TD(x[1])/Q1ROT_K );
}

DOUBLE r_phi3(x)
DOUBLE *x;
{
   return( 0.25 - 0.5*x[0] + (Q1ROT_T(x[0])-Q1ROT_T(x[1]))/Q1ROT_K );
}

DOUBLE r_phi3_0(x)
DOUBLE *x;
{
   return( -0.5 + Q1ROT_TD(x[0])/Q1ROT_K );
}

DOUBLE r_phi3_1(x)
DOUBLE *x;
{
   return(  -Q1ROT_TD(x[1])/Q1ROT_K );
}

DOUBLE sr_phi_0(x,f0,f1,f2,f3)
DOUBLE *x, f0, f1, f2, f3;
{
   return( f0*r_phi0_0(x) + f1*r_phi1_0(x) +
           f2*r_phi2_0(x) + f3*r_phi3_0(x) );
}

DOUBLE sr_phi_1(x,f0,f1,f2,f3)
DOUBLE *x, f0, f1, f2, f3;
{
   return( f0*r_phi0_1(x) + f1*r_phi1_1(x) +
           f2*r_phi2_1(x) + f3*r_phi3_1(x) );
}

DOUBLE phi0(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return(r_phi0(y));
}

DOUBLE phi0_0(h,xc,x)  /* x[0] derivative  */
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi0_0(y));
}

DOUBLE phi0_1(h,xc,x)  /* x[1] derivative  */
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi0_1(y));
}

DOUBLE phi1(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return(r_phi1(y));
}

DOUBLE phi1_0(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi1_0(y));
}

DOUBLE phi1_1(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi1_1(y));
}

DOUBLE phi2(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return(r_phi2(y));
}

DOUBLE phi2_0(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi2_0(y));
}

DOUBLE phi2_1(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi2_1(y));
}

DOUBLE phi3(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return(r_phi3(y));
}

DOUBLE phi3_0(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi3_0(y));
}

DOUBLE phi3_1(h,xc,x)
DOUBLE h, *xc, *x;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*r_phi3_1(y));
}

DOUBLE sphi_0(h,xc,x,f0,f1,f2,f3)
DOUBLE h, *xc, *x, f0, f1, f2, f3;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*sr_phi_0(y,f0,f1,f2,f3));
}

DOUBLE sphi_1(h,xc,x,f0,f1,f2,f3)
DOUBLE h, *xc, *x, f0, f1, f2, f3;
{
   FLOAT y[2]={2.*(x[0]-xc[0])/h, 2.*(x[1]-xc[1])/h};
   return((2./h)*sr_phi_1(y,f0,f1,f2,f3));
}

/* P1C ============= */
DOUBLE_FUNC p1c_basis[3]    = { r_p1_0,   r_p1_1,   r_p1_2 };
DOUBLE_FUNC p1c_basis_0[3]  = { r_p1_0_0, r_p1_1_0, r_p1_2_0 };
DOUBLE_FUNC p1c_basis_1[3]  = { r_p1_0_1, r_p1_1_1, r_p1_2_1 };
DOUBLE_FUNC p1c_basis_00[3] = { zero_fcn, zero_fcn, zero_fcn };
DOUBLE_FUNC p1c_basis_01[3] = { zero_fcn, zero_fcn, zero_fcn };
DOUBLE_FUNC p1c_basis_11[3] = { zero_fcn, zero_fcn, zero_fcn };
DOUBLE_FUNC p1c_dofs[3]    = { rd_p1_0,  rd_p1_1,  rd_p1_2 };

FINITE_ELEMENT p1c_element = { SIMPLEX, 1, 0, 0, 
                               p1c_basis,    NULL, NULL, 
                               p1c_basis_0,  NULL, NULL,
                               p1c_basis_1,  NULL, NULL,
                               p1c_basis_00, NULL, NULL,
                               p1c_basis_01, NULL, NULL,
                               p1c_basis_11, NULL, NULL,
                               p1c_dofs,    NULL, NULL };

/* P1X3C ============= */
DOUBLE_FUNC p1x3_ebasis[1]    = {r_p1x3_b};
DOUBLE_FUNC p1x3_ebasis_0[1]  = {r_p1x3_b_0};
DOUBLE_FUNC p1x3_ebasis_1[1]  = {r_p1x3_b_1};
DOUBLE_FUNC p1x3_ebasis_00[1] = {zero_fcn};
DOUBLE_FUNC p1x3_ebasis_01[1] = {zero_fcn};
DOUBLE_FUNC p1x3_ebasis_11[1] = {zero_fcn};
DOUBLE_FUNC p1x3_edofs[1]    = {rd_p1x3_b};

FINITE_ELEMENT p1x3c_element = { SIMPLEX, 1, 0, 1, 
                                 p1c_basis,    NULL, p1x3_ebasis, 
                                 p1c_basis_0,  NULL, p1x3_ebasis_0,
                                 p1c_basis_1,  NULL, p1x3_ebasis_1,
                                 p1c_basis_00, NULL, p1x3_ebasis_00,
                                 p1c_basis_01, NULL, p1x3_ebasis_01,
                                 p1c_basis_11, NULL, p1x3_ebasis_11,
                                 p1c_dofs,     NULL, p1x3_edofs};

/* P1C_ELBUB ============= */
DOUBLE_FUNC p1cb_ebasis[1] = {r_l0l1l2};
DOUBLE_FUNC p1cb_ebasis_0[1] = {r_l0l1l2_0};
DOUBLE_FUNC p1cb_ebasis_1[1] = {r_l0l1l2_1};
DOUBLE_FUNC p1cb_ebasis_00[1] = {r_l0l1l2_00};
DOUBLE_FUNC p1cb_ebasis_01[1] = {r_l0l1l2_01};
DOUBLE_FUNC p1cb_ebasis_11[1] = {r_l0l1l2_11};
DOUBLE_FUNC p1cb_edofs[1] = {rd_l0l1l2};

FINITE_ELEMENT p1cb_element = { SIMPLEX, 1, 0, 1, 
                                p1c_basis,    NULL, p1cb_ebasis, 
                                p1c_basis_0,  NULL, p1cb_ebasis_0,
                                p1c_basis_1,  NULL, p1cb_ebasis_1,
                                p1c_basis_00, NULL, p1cb_ebasis_00,
                                p1c_basis_01, NULL, p1cb_ebasis_01,
                                p1c_basis_11, NULL, p1cb_ebasis_11,
                                p1c_dofs,     NULL, p1cb_edofs};

/* P2C ============= */
DOUBLE_FUNC p2c_fbasis[3]    = { r_l1l2,    r_l0l2,    r_l0l1 };
DOUBLE_FUNC p2c_fbasis_0[3]  = { r_l1l2_0,  r_l0l2_0,  r_l0l1_0 };
DOUBLE_FUNC p2c_fbasis_1[3]  = { r_l1l2_1,  r_l0l2_1,  r_l0l1_1 };
DOUBLE_FUNC p2c_fbasis_00[3] = { r_l1l2_00, r_l0l2_00, r_l0l1_00 };
DOUBLE_FUNC p2c_fbasis_01[3] = { r_l1l2_01, r_l0l2_01, r_l0l1_01 };
DOUBLE_FUNC p2c_fbasis_11[3] = { r_l1l2_11, r_l0l2_11, r_l0l1_11 };
DOUBLE_FUNC p2c_fdofs[3]     = { rd_p2_f12, rd_p2_f02, rd_p2_f01 };

FINITE_ELEMENT p2c_element = { SIMPLEX, 1, 1, 0, 
                               p1c_basis,    p2c_fbasis,    NULL, 
                               p1c_basis_0,  p2c_fbasis_0,  NULL,
                               p1c_basis_1,  p2c_fbasis_1,  NULL,
                               p1c_basis_00, p2c_fbasis_00, NULL,
                               p1c_basis_01, p2c_fbasis_01, NULL,
                               p1c_basis_11, p2c_fbasis_11, NULL,
                               p1c_dofs,     p2c_fdofs,     NULL };

/* P2X3C ============= */
DOUBLE_FUNC p2x3_ebasis[4]    = { r_p2x3_0,    r_p2x3_1,    r_p2x3_2,    r_p1x3_b };
DOUBLE_FUNC p2x3_ebasis_0[4]  = { r_p2x3_0_0,  r_p2x3_1_0,  r_p2x3_2_0,  r_p1x3_b_0 };
DOUBLE_FUNC p2x3_ebasis_1[4]  = { r_p2x3_0_1,  r_p2x3_1_1,  r_p2x3_2_1,  r_p1x3_b_1 };
DOUBLE_FUNC p2x3_ebasis_00[4] = { r_p2x3_0_00, r_p2x3_1_00, r_p2x3_2_00, zero_fcn };
DOUBLE_FUNC p2x3_ebasis_01[4] = { r_p2x3_0_01, r_p2x3_1_01, r_p2x3_2_01, zero_fcn };
DOUBLE_FUNC p2x3_ebasis_11[4] = { r_p2x3_0_11, r_p2x3_1_11, r_p2x3_2_11, zero_fcn };
DOUBLE_FUNC p2x3_edofs[4]    = { rd_p2x3_0,   rd_p2x3_1,   rd_p2x3_2,   rd_p2x3_b };

FINITE_ELEMENT p2x3c_element = { SIMPLEX, 1, 1, 4,
                                 p1c_basis,    p2c_fbasis,    p2x3_ebasis, 
                                 p1c_basis_0,  p2c_fbasis_0,  p2x3_ebasis_0,
                                 p1c_basis_1,  p2c_fbasis_1,  p2x3_ebasis_1,
                                 p1c_basis_00, p2c_fbasis_00, p2x3_ebasis_00,
                                 p1c_basis_01, p2c_fbasis_01, p2x3_ebasis_01,
                                 p1c_basis_11, p2c_fbasis_11, p2x3_ebasis_11,
                                 p1c_dofs,     p2c_fdofs,     p2x3_edofs};

/* P2C_3ELBUB ============= */
DOUBLE_FUNC p2cnb_ebasis[3]    = { r_l0l0l1l2,    r_l0l1l1l2,    r_l0l1l2l2 };
DOUBLE_FUNC p2cnb_ebasis_0[3]  = { r_l0l0l1l2_0,  r_l0l1l1l2_0,  r_l0l1l2l2_0 };
DOUBLE_FUNC p2cnb_ebasis_1[3]  = { r_l0l0l1l2_1,  r_l0l1l1l2_1,  r_l0l1l2l2_1 };
DOUBLE_FUNC p2cnb_ebasis_00[3] = { r_l0l0l1l2_00, r_l0l1l1l2_00, r_l0l1l2l2_00};
DOUBLE_FUNC p2cnb_ebasis_01[3] = { r_l0l0l1l2_01, r_l0l1l1l2_01, r_l0l1l2l2_01};
DOUBLE_FUNC p2cnb_ebasis_11[3] = { r_l0l0l1l2_11, r_l0l1l1l2_11, r_l0l1l2l2_11};
DOUBLE_FUNC p2cnb_edofs[3]     = {rd_l0l0l1l2,   rd_l0l1l1l2,   rd_l0l1l2l2 };

FINITE_ELEMENT p2c3b_element = { SIMPLEX, 1, 1, 3, 
                               p1c_basis,    p2c_fbasis,    p2cnb_ebasis, 
                               p1c_basis_0,  p2c_fbasis_0,  p2cnb_ebasis_0,
                               p1c_basis_1,  p2c_fbasis_1,  p2cnb_ebasis_1,
                               p1c_basis_00, p2c_fbasis_00, p2cnb_ebasis_00,
                               p1c_basis_01, p2c_fbasis_01, p2cnb_ebasis_01,
                               p1c_basis_11, p2c_fbasis_11, p2cnb_ebasis_11,
                               p1c_dofs,     p2c_fdofs,     p2cnb_edofs };

/* P2C_6ELBUB ============= */
DOUBLE_FUNC p2c6b_ebasis[6]    = { r_l0l0l1l2,    r_l0l1l1l2,    r_l0l1l2l2,
                                   r_l0l0l0l1l2,  r_l0l1l1l1l2,  r_l0l1l2l2l2 };
DOUBLE_FUNC p2c6b_ebasis_0[6]  = { r_l0l0l1l2_0,  r_l0l1l1l2_0,  r_l0l1l2l2_0,
                                   r_l0l0l0l1l2_0,r_l0l1l1l1l2_0,r_l0l1l2l2l2_0 };
DOUBLE_FUNC p2c6b_ebasis_1[6]  = { r_l0l0l1l2_1,  r_l0l1l1l2_1,  r_l0l1l2l2_1,
                                   r_l0l0l0l1l2_1,r_l0l1l1l1l2_1,r_l0l1l2l2l2_1 };
DOUBLE_FUNC p2c6b_ebasis_00[6] = { r_l0l0l1l2_00, r_l0l1l1l2_00, r_l0l1l2l2_00,
                                   r_l0l0l0l1l2_00,r_l0l1l1l1l2_00,r_l0l1l2l2l2_00 };
DOUBLE_FUNC p2c6b_ebasis_01[6] = { r_l0l0l1l2_01, r_l0l1l1l2_01, r_l0l1l2l2_01,
                                   r_l0l0l0l1l2_01,r_l0l1l1l1l2_01,r_l0l1l2l2l2_01 };
DOUBLE_FUNC p2c6b_ebasis_11[6] = { r_l0l0l1l2_11, r_l0l1l1l2_11, r_l0l1l2l2_11,
                                   r_l0l0l0l1l2_11,r_l0l1l1l1l2_11,r_l0l1l2l2l2_11 };
DOUBLE_FUNC p2c6b_edofs[6]     = {rd_l0l0l1l2,   rd_l0l1l1l2,   rd_l0l1l2l2,
                                  rd_l0l0l0l1l2, rd_l0l1l1l1l2, rd_l0l1l2l2l2 };

FINITE_ELEMENT p2c6b_element = { SIMPLEX, 1, 1, 6, 
                               p1c_basis,    p2c_fbasis,    p2c6b_ebasis, 
                               p1c_basis_0,  p2c_fbasis_0,  p2c6b_ebasis_0,
                               p1c_basis_1,  p2c_fbasis_1,  p2c6b_ebasis_1,
                               p1c_basis_00, p2c_fbasis_00, p2c6b_ebasis_00,
                               p1c_basis_01, p2c_fbasis_01, p2c6b_ebasis_01,
                               p1c_basis_11, p2c_fbasis_11, p2c6b_ebasis_11,
                               p1c_dofs,     p2c_fdofs,     p2c6b_edofs };

/* P3S ============= */

DOUBLE_FUNC p3s0_basis[3] = { r_p3_01, r_p3_12, r_p3_20 };
DOUBLE_FUNC p3s0_basis_0[3] = { r_p3_01_0, r_p3_12_0, r_p3_20_0 };
DOUBLE_FUNC p3s0_basis_1[3] = { r_p3_01_1, r_p3_12_1, r_p3_20_1 };

DOUBLE_FUNC p3s1_basis[3] = { r_p3_02, r_p3_10, r_p3_21 };
DOUBLE_FUNC p3s1_basis_0[3] = { r_p3_02_0, r_p3_10_0, r_p3_21_0 };
DOUBLE_FUNC p3s1_basis_1[3] = { r_p3_02_1, r_p3_10_1, r_p3_21_1 };

/* P4S ============= */

DOUBLE_FUNC p4s_basis[3] = { r_0012, r_1012, r_2012 };
DOUBLE_FUNC p4s_basis_0[3] = { r_0012_0, r_1012_0, r_2012_0 };
DOUBLE_FUNC p4s_basis_1[3] = { r_0012_1, r_1012_1, r_2012_1 };

/* Q1C ============= */
DOUBLE_FUNC q1c_basis[4]    = { r_q1_0,    r_q1_1,    r_q1_2,    r_q1_3 };
DOUBLE_FUNC q1c_basis_0[4]  = { r_q1_0_0,  r_q1_1_0,  r_q1_2_0,  r_q1_3_0 };
DOUBLE_FUNC q1c_basis_1[4]  = { r_q1_0_1,  r_q1_1_1,  r_q1_2_1,  r_q1_3_1 };
DOUBLE_FUNC q1c_basis_00[4] = { zero_fcn,  zero_fcn,  zero_fcn,  zero_fcn };
DOUBLE_FUNC q1c_basis_01[4] = { r_q1_0_01, r_q1_1_01, r_q1_2_01, r_q1_3_01 };
DOUBLE_FUNC q1c_basis_11[4] = { zero_fcn,  zero_fcn,  zero_fcn,  zero_fcn };
DOUBLE_FUNC q1c_dofs[4]    = { rd_q1_0,   rd_q1_1,   rd_q1_2,   rd_q1_3 };

FINITE_ELEMENT q1c_element = { CUBE, 1, 0, 0, 
                               q1c_basis,    NULL, NULL, 
                               q1c_basis_0,  NULL, NULL,
                               q1c_basis_1,  NULL, NULL,
                               q1c_basis_00, NULL, NULL,
                               q1c_basis_01, NULL, NULL,
                               q1c_basis_11, NULL, NULL,
                               q1c_dofs,     NULL, NULL };

/* Q1X4C ============= */
DOUBLE_FUNC q1x4_nbasis[4]    = { r_q1x4_0,    r_q1x4_1,    r_q1x4_2,    r_q1x4_3 };
DOUBLE_FUNC q1x4_nbasis_0[4]  = { r_q1x4_0_0,  r_q1x4_1_0,  r_q1x4_2_0,  r_q1x4_3_0 };
DOUBLE_FUNC q1x4_nbasis_1[4]  = { r_q1x4_0_1,  r_q1x4_1_1,  r_q1x4_2_1,  r_q1x4_3_1 };
DOUBLE_FUNC q1x4_nbasis_00[4] = { zero_fcn,    zero_fcn,    zero_fcn,    zero_fcn };
DOUBLE_FUNC q1x4_nbasis_01[4] = { r_q1x4_0_01, r_q1x4_1_01, r_q1x4_2_01, r_q1x4_3_01 };
DOUBLE_FUNC q1x4_nbasis_11[4] = { zero_fcn,    zero_fcn,    zero_fcn,    zero_fcn };
DOUBLE_FUNC q1x4_ndofs[4]    = { rd_q1x4_0,   rd_q1x4_1,   rd_q1x4_2,   rd_q1x4_3 };
DOUBLE_FUNC q1x4_fbasis[4]    = { r_q1x4_01,    r_q1x4_12,    r_q1x4_23,    r_q1x4_30 };
DOUBLE_FUNC q1x4_fbasis_0[4]  = { r_q1x4_01_0,  r_q1x4_12_0,  r_q1x4_23_0,  r_q1x4_30_0 };
DOUBLE_FUNC q1x4_fbasis_1[4]  = { r_q1x4_01_1,  r_q1x4_12_1,  r_q1x4_23_1,  r_q1x4_30_1 };
DOUBLE_FUNC q1x4_fbasis_00[4] = { zero_fcn,     zero_fcn,     zero_fcn,     zero_fcn };
DOUBLE_FUNC q1x4_fbasis_01[4] = { r_q1x4_01_01, r_q1x4_12_01, r_q1x4_23_01, r_q1x4_30_01 };
DOUBLE_FUNC q1x4_fbasis_11[4] = { zero_fcn,     zero_fcn,     zero_fcn,     zero_fcn };
DOUBLE_FUNC q1x4_fdofs[4]    = { rd_q1x4_01,   rd_q1x4_12,   rd_q1x4_23,   rd_q1x4_30 };
DOUBLE_FUNC q1x4_ebasis[1]    = {r_q1x4_b};
DOUBLE_FUNC q1x4_ebasis_0[1]  = {r_q1x4_b_0};
DOUBLE_FUNC q1x4_ebasis_1[1]  = {r_q1x4_b_1};
DOUBLE_FUNC q1x4_ebasis_00[1] = {zero_fcn};
DOUBLE_FUNC q1x4_ebasis_01[1] = {r_q1x4_b_01};
DOUBLE_FUNC q1x4_ebasis_11[1] = {zero_fcn};
DOUBLE_FUNC q1x4_edofs[1]    = {rd_q1x4_b};

FINITE_ELEMENT q1x4c_element = { CUBE, 1, 1, 1, 
                               q1x4_nbasis,    q1x4_fbasis,    q1x4_ebasis, 
                               q1x4_nbasis_0,  q1x4_fbasis_0,  q1x4_ebasis_0,
                               q1x4_nbasis_1,  q1x4_fbasis_1,  q1x4_ebasis_1,
                               q1x4_nbasis_00, q1x4_fbasis_00, q1x4_ebasis_00,
                               q1x4_nbasis_01, q1x4_fbasis_01, q1x4_ebasis_01,
                               q1x4_nbasis_11, q1x4_fbasis_11, q1x4_ebasis_11,
                               q1x4_ndofs,     q1x4_fdofs,     q1x4_edofs };

/* Q1C_ELBUB ============= */
DOUBLE_FUNC q2c_ebasis[1]    = {r_q2_b};
DOUBLE_FUNC q2c_ebasis_0[1]  = {r_q2_b_0};
DOUBLE_FUNC q2c_ebasis_1[1]  = {r_q2_b_1};
DOUBLE_FUNC q2c_ebasis_00[1] = {r_q2_b_00};
DOUBLE_FUNC q2c_ebasis_01[1] = {r_q2_b_01};
DOUBLE_FUNC q2c_ebasis_11[1] = {r_q2_b_11};
DOUBLE_FUNC q1cb_edofs[1]    = {rd_q1cb_b};

FINITE_ELEMENT q1cb_element = { CUBE, 1, 0, 1,
                                q1c_basis,    NULL, q2c_ebasis,
                                q1c_basis_0,  NULL, q2c_ebasis_0,
                                q1c_basis_1,  NULL, q2c_ebasis_1,
                                q1c_basis_00, NULL, q2c_ebasis_00,
                                q1c_basis_01, NULL, q2c_ebasis_01,
                                q1c_basis_11, NULL, q2c_ebasis_11,
                                q1c_dofs,     NULL, q1cb_edofs };

/* Q2C ============= */
DOUBLE_FUNC q2c_fbasis[4]    = { r_q2_0,    r_q2_1,    r_q2_2,    r_q2_3 };
DOUBLE_FUNC q2c_fbasis_0[4]  = { r_q2_0_0,  r_q2_1_0,  r_q2_2_0,  r_q2_3_0 };
DOUBLE_FUNC q2c_fbasis_1[4]  = { r_q2_0_1,  r_q2_1_1,  r_q2_2_1,  r_q2_3_1 };
DOUBLE_FUNC q2c_fbasis_00[4] = { r_q2_0_00, r_q2_1_00, r_q2_2_00, r_q2_3_00 };
DOUBLE_FUNC q2c_fbasis_01[4] = { r_q2_0_01, r_q2_1_01, r_q2_2_01, r_q2_3_01 };
DOUBLE_FUNC q2c_fbasis_11[4] = { r_q2_0_11, r_q2_1_11, r_q2_2_11, r_q2_3_11 };
DOUBLE_FUNC q2c_fdofs[4]    = { rd_q2_0,   rd_q2_1,   rd_q2_2,   rd_q2_3 };
DOUBLE_FUNC q2c_edofs[1]    = {rd_q2_b};

FINITE_ELEMENT q2c_element = { CUBE, 1, 1, 1,
                               q1c_basis,    q2c_fbasis,    q2c_ebasis,
                               q1c_basis_0,  q2c_fbasis_0,  q2c_ebasis_0,
                               q1c_basis_1,  q2c_fbasis_1,  q2c_ebasis_1,
                               q1c_basis_00, q2c_fbasis_00, q2c_ebasis_00,
                               q1c_basis_01, q2c_fbasis_01, q2c_ebasis_01,
                               q1c_basis_11, q2c_fbasis_11, q2c_ebasis_11,
                               q1c_dofs,     q2c_fdofs,     q2c_edofs };

/* Q2X4C ============= */
DOUBLE_FUNC q2x4_fbasis[12]    = { r_q2x4_00,    r_q1x4_01,    r_q2x4_01,  
                                   r_q2x4_11,    r_q1x4_12,    r_q2x4_12, 
                                   r_q2x4_22,    r_q1x4_23,    r_q2x4_23, 
                                   r_q2x4_33,    r_q1x4_30,    r_q2x4_30 };
DOUBLE_FUNC q2x4_fbasis_0[12]  = { r_q2x4_00_0,  r_q1x4_01_0,  r_q2x4_01_0,
                                   r_q2x4_11_0,  r_q1x4_12_0,  r_q2x4_12_0,
                                   r_q2x4_22_0,  r_q1x4_23_0,  r_q2x4_23_0,
                                   r_q2x4_33_0,  r_q1x4_30_0,  r_q2x4_30_0 };
DOUBLE_FUNC q2x4_fbasis_1[12]  = { r_q2x4_00_1,  r_q1x4_01_1,  r_q2x4_01_1,
                                   r_q2x4_11_1,  r_q1x4_12_1,  r_q2x4_12_1,
                                   r_q2x4_22_1,  r_q1x4_23_1,  r_q2x4_23_1,
                                   r_q2x4_33_1,  r_q1x4_30_1,  r_q2x4_30_1 };
DOUBLE_FUNC q2x4_fbasis_00[12] = { r_q2x4_00_00, zero_fcn,     r_q2x4_01_00,
                                   r_q2x4_11_00, zero_fcn,     r_q2x4_12_00,
                                   r_q2x4_22_00, zero_fcn,     r_q2x4_23_00,
                                   r_q2x4_33_00, zero_fcn,     r_q2x4_30_00 };
DOUBLE_FUNC q2x4_fbasis_01[12] = { r_q2x4_00_01, r_q1x4_01_01, r_q2x4_01_01,
                                   r_q2x4_11_01, r_q1x4_12_01, r_q2x4_12_01,
                                   r_q2x4_22_01, r_q1x4_23_01, r_q2x4_23_01,
                                   r_q2x4_33_01, r_q1x4_30_01, r_q2x4_30_01 };
DOUBLE_FUNC q2x4_fbasis_11[12] = { r_q2x4_00_11, zero_fcn,     r_q2x4_01_11,
                                   r_q2x4_11_11, zero_fcn,     r_q2x4_12_11,
                                   r_q2x4_22_11, zero_fcn,     r_q2x4_23_11,
                                   r_q2x4_33_11, zero_fcn,     r_q2x4_30_11 };
DOUBLE_FUNC q2x4_fdofs[12]    = { rd_q2x4_00,   rd_q1x4_01,   rd_q2x4_01,
                                  rd_q2x4_11,   rd_q1x4_12,   rd_q2x4_12,
                                  rd_q2x4_22,   rd_q1x4_23,   rd_q2x4_23,
                                  rd_q2x4_33,   rd_q1x4_30,   rd_q2x4_30 };
DOUBLE_FUNC q2x4_ebasis[9]    = { r_q2x4_0,     r_q2x4_1,     r_q2x4_2,     r_q2x4_3,
                                  r_q2x4_b0,    r_q2x4_b1,    r_q2x4_b2,    r_q2x4_b3,
                                  r_q1x4_b };
DOUBLE_FUNC q2x4_ebasis_0[9]  = { r_q2x4_0_0,   r_q2x4_1_0,   r_q2x4_2_0,   r_q2x4_3_0,
                                  r_q2x4_b0_0,  r_q2x4_b1_0,  r_q2x4_b2_0,  r_q2x4_b3_0,
                                  r_q1x4_b_0 };
DOUBLE_FUNC q2x4_ebasis_1[9]  = { r_q2x4_0_1,   r_q2x4_1_1,   r_q2x4_2_1,   r_q2x4_3_1,
                                  r_q2x4_b0_1,  r_q2x4_b1_1,  r_q2x4_b2_1,  r_q2x4_b3_1,
                                  r_q1x4_b_1 };
DOUBLE_FUNC q2x4_ebasis_00[9] = { r_q2x4_0_00,  r_q2x4_1_00,  r_q2x4_2_00,  r_q2x4_3_00,
                                  r_q2x4_b0_00, r_q2x4_b1_00, r_q2x4_b2_00, r_q2x4_b3_00,
                                  zero_fcn };
DOUBLE_FUNC q2x4_ebasis_01[9] = { r_q2x4_0_01,  r_q2x4_1_01,  r_q2x4_2_01,  r_q2x4_3_01,
                                  r_q2x4_b0_01, r_q2x4_b1_01, r_q2x4_b2_01, r_q2x4_b3_01,
                                  r_q1x4_b_01};
DOUBLE_FUNC q2x4_ebasis_11[9] = { r_q2x4_0_11,  r_q2x4_1_11,  r_q2x4_2_11,  r_q2x4_3_11,
                                  r_q2x4_b0_11, r_q2x4_b1_11, r_q2x4_b2_11, r_q2x4_b3_11,
                                  zero_fcn};
DOUBLE_FUNC q2x4_edofs[9]    = { rd_q2x4_0,    rd_q2x4_1,    rd_q2x4_2,    rd_q2x4_3,
                                 rd_q2x4_b0,   rd_q2x4_b1,   rd_q2x4_b2,   rd_q2x4_b3,
                                 rd_q1x4_b };

FINITE_ELEMENT q2x4c_element = { CUBE, 1, 3, 9, 
                               q1x4_nbasis,    q2x4_fbasis,    q2x4_ebasis, 
                               q1x4_nbasis_0,  q2x4_fbasis_0,  q2x4_ebasis_0,
                               q1x4_nbasis_1,  q2x4_fbasis_1,  q2x4_ebasis_1,
                               q1x4_nbasis_00, q2x4_fbasis_00, q2x4_ebasis_00,
                               q1x4_nbasis_01, q2x4_fbasis_01, q2x4_ebasis_01,
                               q1x4_nbasis_11, q2x4_fbasis_11, q2x4_ebasis_11,
                               q1x4_ndofs,     q2x4_fdofs,     q2x4_edofs };

/* Q2C_2ELBUB ============= */
DOUBLE_FUNC q2cnb_ebasis[3]    = {r_q2_b,    r_q2nb0,    r_q2nb1 };
DOUBLE_FUNC q2cnb_ebasis_0[3]  = {r_q2_b_0,  r_q2nb0_0,  r_q2nb1_0 };
DOUBLE_FUNC q2cnb_ebasis_1[3]  = {r_q2_b_1,  r_q2nb0_1,  r_q2nb1_1 };
DOUBLE_FUNC q2cnb_ebasis_00[3] = {r_q2_b_00, r_q2nb0_00, r_q2nb1_00 };
DOUBLE_FUNC q2cnb_ebasis_01[3] = {r_q2_b_01, r_q2nb0_01, r_q2nb1_01 };
DOUBLE_FUNC q2cnb_ebasis_11[3] = {r_q2_b_11, r_q2nb0_11, r_q2nb1_11 };
DOUBLE_FUNC q2cnb_edofs[3]    = {rd_q2_b,   rd_q2nb0,   rd_q2nb1 };

FINITE_ELEMENT q2c2b_element = { CUBE, 1, 1, 3,
                                 q1c_basis,    q2c_fbasis,    q2cnb_ebasis,
                                 q1c_basis_0,  q2c_fbasis_0,  q2cnb_ebasis_0,
                                 q1c_basis_1,  q2c_fbasis_1,  q2cnb_ebasis_1,
                                 q1c_basis_00, q2c_fbasis_00, q2cnb_ebasis_00,
                                 q1c_basis_01, q2c_fbasis_01, q2cnb_ebasis_01,
                                 q1c_basis_11, q2c_fbasis_11, q2cnb_ebasis_11,
                                 q1c_dofs,     q2c_fdofs,     q2cnb_edofs };

/* Q2C_3ELBUB ============= */
DOUBLE_FUNC q2c3b_ebasis[4]    = {r_q2_b,    r_q2nb0,    r_q2nb1,    r_q2nb01 };
DOUBLE_FUNC q2c3b_ebasis_0[4]  = {r_q2_b_0,  r_q2nb0_0,  r_q2nb1_0,  r_q2nb01_0 };
DOUBLE_FUNC q2c3b_ebasis_1[4]  = {r_q2_b_1,  r_q2nb0_1,  r_q2nb1_1,  r_q2nb01_1 };
DOUBLE_FUNC q2c3b_ebasis_00[4] = {r_q2_b_00, r_q2nb0_00, r_q2nb1_00, r_q2nb01_00 };
DOUBLE_FUNC q2c3b_ebasis_01[4] = {r_q2_b_01, r_q2nb0_01, r_q2nb1_01, r_q2nb01_01 };
DOUBLE_FUNC q2c3b_ebasis_11[4] = {r_q2_b_11, r_q2nb0_11, r_q2nb1_11, r_q2nb01_11 };
DOUBLE_FUNC q2c3b_edofs[4]    = {rd_q2_b,   rd_q2nb0,   rd_q2nb1,   rd_q2nb01 };

FINITE_ELEMENT q2c3b_element = { CUBE, 1, 1, 4,
                                 q1c_basis,    q2c_fbasis,    q2c3b_ebasis,
                                 q1c_basis_0,  q2c_fbasis_0,  q2c3b_ebasis_0,
                                 q1c_basis_1,  q2c_fbasis_1,  q2c3b_ebasis_1,
                                 q1c_basis_00, q2c_fbasis_00, q2c3b_ebasis_00,
                                 q1c_basis_01, q2c_fbasis_01, q2c3b_ebasis_01,
                                 q1c_basis_11, q2c_fbasis_11, q2c3b_ebasis_11,
                                 q1c_dofs,     q2c_fdofs,     q2c3b_edofs };

/* Q2B3C ============= */
DOUBLE_FUNC q2b3_fbasis[8]    = { r_q2_0,    r_q3_4,  
                                  r_q2_1,    r_q3_5, 
                                  r_q2_2,    r_q3_10, 
                                  r_q2_3,    r_q3_11 };
DOUBLE_FUNC q2b3_fbasis_0[8]  = { r_q2_0_0,  r_q3_4_0,
                                  r_q2_1_0,  r_q3_5_0,
                                  r_q2_2_0,  r_q3_10_0,
                                  r_q2_3_0,  r_q3_11_0 };
DOUBLE_FUNC q2b3_fbasis_1[8]  = { r_q2_0_1,  r_q3_4_1,
                                  r_q2_1_1,  r_q3_5_1,
                                  r_q2_2_1,  r_q3_10_1,
                                  r_q2_3_1,  r_q3_11_1 };
DOUBLE_FUNC q2b3_fbasis_00[8] = { r_q2_0_00, r_q3_4_00,
                                  r_q2_1_00, r_q3_5_00,
                                  r_q2_2_00, r_q3_10_00,
                                  r_q2_3_00, r_q3_11_00 };
DOUBLE_FUNC q2b3_fbasis_01[8] = { r_q2_0_01, r_q3_4_01,
                                  r_q2_1_01, r_q3_5_01,
                                  r_q2_2_01, r_q3_10_01,
                                  r_q2_3_01, r_q3_11_01 };
DOUBLE_FUNC q2b3_fbasis_11[8] = { r_q2_0_11, r_q3_4_11,
                                  r_q2_1_11, r_q3_5_11,
                                  r_q2_2_11, r_q3_10_11,
                                  r_q2_3_11, r_q3_11_11 };
DOUBLE_FUNC q2b3_fdofs[8]     = { no_rd, no_rd, no_rd, no_rd, 
                                  no_rd, no_rd, no_rd, no_rd };
DOUBLE_FUNC q2b3_ebasis[4]    = { r_q2_b,    r_q3_13,    r_q3_14,    r_q3_15 };
DOUBLE_FUNC q2b3_ebasis_0[4]  = { r_q2_b_0,  r_q3_13_0,  r_q3_14_0,  r_q3_15_0 };
DOUBLE_FUNC q2b3_ebasis_1[4]  = { r_q2_b_1,  r_q3_13_1,  r_q3_14_1,  r_q3_15_1 };
DOUBLE_FUNC q2b3_ebasis_00[4] = { r_q2_b_00, r_q3_13_00, r_q3_14_00, r_q3_15_00};
DOUBLE_FUNC q2b3_ebasis_01[4] = { r_q2_b_01, r_q3_13_01, r_q3_14_01, r_q3_15_01};
DOUBLE_FUNC q2b3_ebasis_11[4] = { r_q2_b_11, r_q3_13_11, r_q3_14_11, r_q3_15_11};
DOUBLE_FUNC q2b3_edofs[4]     = { no_rd, no_rd, no_rd, no_rd };

FINITE_ELEMENT q2b3c_element = { CUBE, 1, 2, 4, 
                                 q1c_basis,    q2b3_fbasis,    q2b3_ebasis, 
                                 q1c_basis_0,  q2b3_fbasis_0,  q2b3_ebasis_0,
                                 q1c_basis_1,  q2b3_fbasis_1,  q2b3_ebasis_1,
                                 q1c_basis_00, q2b3_fbasis_00, q2b3_ebasis_00,
                                 q1c_basis_01, q2b3_fbasis_01, q2b3_ebasis_01,
                                 q1c_basis_11, q2b3_fbasis_11, q2b3_ebasis_11,
                                 q1c_dofs,     q2b3_fdofs,     q2b3_edofs };

/* Q3C ============= */
DOUBLE_FUNC q3_nbasis[4]    = { r_q3_0,    r_q3_1,    r_q3_2,    r_q3_3 };
DOUBLE_FUNC q3_nbasis_0[4]  = { r_q3_0_0,  r_q3_1_0,  r_q3_2_0,  r_q3_3_0 };
DOUBLE_FUNC q3_nbasis_1[4]  = { r_q3_0_1,  r_q3_1_1,  r_q3_2_1,  r_q3_3_1 };
DOUBLE_FUNC q3_nbasis_00[4] = { r_q3_0_00, r_q3_1_00, r_q3_2_00, r_q3_3_00 };
DOUBLE_FUNC q3_nbasis_01[4] = { r_q3_0_01, r_q3_1_01, r_q3_2_01, r_q3_3_01 };
DOUBLE_FUNC q3_nbasis_11[4] = { r_q3_0_11, r_q3_1_11, r_q3_2_11, r_q3_3_11 };
DOUBLE_FUNC q3_ndofs[4]    = { rd_q3_0,   rd_q3_1,   rd_q3_2,   rd_q3_3 };
DOUBLE_FUNC q3_fbasis[8]    = { r_q3_4,    r_q3_8,  
                                r_q3_5,    r_q3_9, 
                                r_q3_6,    r_q3_10, 
                                r_q3_7,    r_q3_11 };
DOUBLE_FUNC q3_fbasis_0[8]  = { r_q3_4_0,  r_q3_8_0,  
                                r_q3_5_0,  r_q3_9_0,
                                r_q3_6_0,  r_q3_10_0,
                                r_q3_7_0,  r_q3_11_0 };
DOUBLE_FUNC q3_fbasis_1[8]  = { r_q3_4_1,  r_q3_8_1,  
                                r_q3_5_1,  r_q3_9_1,
                                r_q3_6_1,  r_q3_10_1,
                                r_q3_7_1,  r_q3_11_1 };
DOUBLE_FUNC q3_fbasis_00[8] = { r_q3_4_00, r_q3_8_00,  
                                r_q3_5_00, r_q3_9_00,
                                r_q3_6_00, r_q3_10_00,
                                r_q3_7_00, r_q3_11_00 };
DOUBLE_FUNC q3_fbasis_01[8] = { r_q3_4_01, r_q3_8_01,  
                                r_q3_5_01, r_q3_9_01,
                                r_q3_6_01, r_q3_10_01,
                                r_q3_7_01, r_q3_11_01 };
DOUBLE_FUNC q3_fbasis_11[8] = { r_q3_4_11, r_q3_8_11,  
                                r_q3_5_11, r_q3_9_11,
                                r_q3_6_11, r_q3_10_11,
                                r_q3_7_11, r_q3_11_11 };
DOUBLE_FUNC q3_fdofs[8]    = { rd_q3_4,   rd_q3_8,  
                               rd_q3_5,   rd_q3_9, 
                               rd_q3_6,   rd_q3_10, 
                               rd_q3_7,   rd_q3_11 };
DOUBLE_FUNC q3_ebasis[4]    = { r_q3_12,    r_q3_13,    r_q3_14,    r_q3_15 };
DOUBLE_FUNC q3_ebasis_0[4]  = { r_q3_12_0,  r_q3_13_0,  r_q3_14_0,  r_q3_15_0 };
DOUBLE_FUNC q3_ebasis_1[4]  = { r_q3_12_1,  r_q3_13_1,  r_q3_14_1,  r_q3_15_1 };
DOUBLE_FUNC q3_ebasis_00[4] = { r_q3_12_00, r_q3_13_00, r_q3_14_00, r_q3_15_00};
DOUBLE_FUNC q3_ebasis_01[4] = { r_q3_12_01, r_q3_13_01, r_q3_14_01, r_q3_15_01};
DOUBLE_FUNC q3_ebasis_11[4] = { r_q3_12_11, r_q3_13_11, r_q3_14_11, r_q3_15_11};
DOUBLE_FUNC q3_edofs[4]    = { rd_q3_12,   rd_q3_13,   rd_q3_14,   rd_q3_15 };

FINITE_ELEMENT q3c_element = { CUBE, 1, 2, 4, 
                               q3_nbasis,    q3_fbasis,    q3_ebasis, 
                               q3_nbasis_0,  q3_fbasis_0,  q3_ebasis_0,
                               q3_nbasis_1,  q3_fbasis_1,  q3_ebasis_1,
                               q3_nbasis_00, q3_fbasis_00, q3_ebasis_00,
                               q3_nbasis_01, q3_fbasis_01, q3_ebasis_01,
                               q3_nbasis_11, q3_fbasis_11, q3_ebasis_11,
                               q3_ndofs,     q3_fdofs,     q3_edofs };

#if N_DATA & SCALAR_NODE_DATA

FLOAT sp1c_value_ref(pel,x,u)
ELEMENT *pel;
FLOAT *x;
INT u;
{
   return(REF_LIN(x,NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u)));
}

FLOAT sp1c_value(n0,n1,n2,x,b,u)
NODE *n0, *n1, *n2;
FLOAT *x, b[DIM2][DIM2];
INT u;
{
   return(NDS(n0,u)*LINV(b[0],x)+NDS(n1,u)*LINV(b[1],x)+NDS(n2,u)*LINV(b[2],x));
}

void sp1c_grad(n0,n1,n2,b,u,grad)
NODE *n0, *n1, *n2;
FLOAT b[DIM2][DIM2], *grad;
INT u;
{
   grad[0] = NDS(n0,u)*b[0][0]+NDS(n1,u)*b[1][0]+NDS(n2,u)*b[2][0];
   grad[1] = NDS(n0,u)*b[0][1]+NDS(n1,u)*b[1][1]+NDS(n2,u)*b[2][1];
}

#else

FLOAT sp1c_value_ref(pel,x,u)
ELEMENT *pel; FLOAT *x; INT u;
{  eprintf("Error: sp1c_value_ref not available.\n"); return(0.);  }

FLOAT sp1c_value(n0,n1,n2,x,b,u)
NODE *n0, *n1, *n2; FLOAT *x, b[DIM2][DIM2]; INT u;
{  eprintf("Error: sp1c_value not available.\n"); return(0.);  }

void sp1c_grad(n0,n1,n2,b,u,grad)
NODE *n0, *n1, *n2; FLOAT b[DIM2][DIM2], *grad; INT u;
{  eprintf("Error: sp1c_grad not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (DIM == 2)

void vp1c_value(n0,n1,n2,x,b,u,val)
NODE *n0, *n1, *n2;
FLOAT *x, b[DIM2][DIM2], *val;
INT u;
{
   DOUBLE l0=LINV(b[0],x), l1=LINV(b[1],x), l2=LINV(b[2],x);

   val[0] = ND(n0,u,0)*l0 + ND(n1,u,0)*l1 + ND(n2,u,0)*l2;
   val[1] = ND(n0,u,1)*l0 + ND(n1,u,1)*l1 + ND(n2,u,1)*l2;
}

#else

void vp1c_value(n0,n1,n2,x,b,u,val)
NODE *n0, *n1, *n2; FLOAT *x, b[DIM2][DIM2], *val; INT u;
{  eprintf("Error: vp1c_value not available.\n");  }

#endif

#if F_DATA & SCALAR_FACE_DATA

FLOAT sp1nc_value(fa0,fa1,fa2,x,b,u)
FACE *fa0, *fa1, *fa2;
FLOAT *x, b[DIM2][DIM2];
INT u;
{
   return(FD(fa0,u)*(1.-2.*LINV(b[0],x))+FD(fa1,u)*(1.-2.*LINV(b[1],x))
                                        +FD(fa2,u)*(1.-2.*LINV(b[2],x)));
}

void sp1nc_grad(fa0,fa1,fa2,b,u,grad)
FACE *fa0, *fa1, *fa2;
FLOAT b[DIM2][DIM2], *grad;
INT u;
{
   grad[0] = -2.*(FD(fa0,u)*b[0][0]+FD(fa1,u)*b[1][0]+FD(fa2,u)*b[2][0]);
   grad[1] = -2.*(FD(fa0,u)*b[0][1]+FD(fa1,u)*b[1][1]+FD(fa2,u)*b[2][1]);
}

#else

FLOAT sp1nc_value(fa0,fa1,fa2,x,b,u)
FACE *fa0, *fa1, *fa2; FLOAT *x, b[DIM2][DIM2]; INT u;
{  eprintf("Error: sp1nc_value not available.\n");  }

void sp1nc_grad(fa0,fa1,fa2,b,u,grad)
FACE *fa0, *fa1, *fa2; FLOAT b[DIM2][DIM2], *grad; INT u;
{  eprintf("Error: sp1nc_grad not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT sp2c_value_ref(pel,x,u)
ELEMENT *pel;
FLOAT *x;
INT u;
{
   return(REF_QUADR(x,NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u),
                       FD(pel->f[2],u), FD(pel->f[1],u), FD(pel->f[0],u)));
}

FLOAT sp2c_value(n0,n1,n2,fa0,fa1,fa2,x,b,u)
NODE *n0, *n1, *n2;
FACE *fa0, *fa1, *fa2;
FLOAT *x, b[DIM2][DIM2];
INT u;
{
   return(NDS(n0,u)*LINV(b[0],x)+NDS(n1,u)*LINV(b[1],x)+NDS(n2,u)*LINV(b[2],x)
      +FD(fa0,u)*LINV(b[1],x)*LINV(b[2],x)+FD(fa1,u)*LINV(b[0],x)*LINV(b[2],x)
                                          +FD(fa2,u)*LINV(b[0],x)*LINV(b[1],x));
}

void sp2c_grad(n0,n1,n2,fa0,fa1,fa2,x,b,u,grad)
NODE *n0, *n1, *n2;
FACE *fa0, *fa1, *fa2;
FLOAT *x, b[DIM2][DIM2], *grad;
INT u;
{
   grad[0] = NDS(n0,u)*b[0][0]+NDS(n1,u)*b[1][0]+NDS(n2,u)*b[2][0]
            +FD(fa0,u)*(b[2][0]*LINV(b[1],x) + b[1][0]*LINV(b[2],x))
            +FD(fa1,u)*(b[2][0]*LINV(b[0],x) + b[0][0]*LINV(b[2],x))
            +FD(fa2,u)*(b[1][0]*LINV(b[0],x) + b[0][0]*LINV(b[1],x));
   grad[1] = NDS(n0,u)*b[0][1]+NDS(n1,u)*b[1][1]+NDS(n2,u)*b[2][1]
            +FD(fa0,u)*(b[2][1]*LINV(b[1],x) + b[1][1]*LINV(b[2],x))
            +FD(fa1,u)*(b[2][1]*LINV(b[0],x) + b[0][1]*LINV(b[2],x))
            +FD(fa2,u)*(b[1][1]*LINV(b[0],x) + b[0][1]*LINV(b[1],x));
}

FLOAT sp2c_Lapl(fa0,fa1,fa2,b,u)
FACE *fa0, *fa1, *fa2;
FLOAT b[DIM2][DIM2];
INT u;
{
   return(2.*(FD(fa0,u)*DOT(b[1],b[2])+FD(fa1,u)*DOT(b[0],b[2])
                                                    +FD(fa2,u)*DOT(b[0],b[1])));
}

#else

FLOAT sp2c_value_ref(pel,x,u)
ELEMENT *pel; FLOAT *x; INT u;
{  eprintf("Error: sp2c_value_ref not available.\n");  }

FLOAT sp2c_value(n0,n1,n2,fa0,fa1,fa2,x,b,u)
NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; FLOAT *x, b[DIM2][DIM2]; INT u;
{  eprintf("Error: sp2c_value not available.\n");  }

void sp2c_grad(n0,n1,n2,fa0,fa1,fa2,x,b,u,grad)
NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; FLOAT *x, b[DIM2][DIM2], *grad; INT u;
{  eprintf("Error: sp2c_grad not available.\n");  }

FLOAT sp2c_Lapl(fa0,fa1,fa2,b,u)
FACE *fa0, *fa1, *fa2; FLOAT b[DIM2][DIM2]; INT u;
{  eprintf("Error: sp2c_Lapl not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

FLOAT sq1c_value_ref(pel,x,u)
ELEMENT *pel;
FLOAT *x;
INT u;
{
   return(REF_BILIN(x,NDS(pel->n[0],u),NDS(pel->n[1],u),
                      NDS(pel->n[2],u),NDS(pel->n[3],u)));
}

#else

FLOAT sq1c_value_ref(pel,x,u)
ELEMENT *pel; FLOAT *x; INT u;
{  eprintf("Error: sq1c_value_ref not available.\n");  }

#endif

#if (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

FLOAT sq1rot_value_ref(pel,x,u)
ELEMENT *pel;
FLOAT *x;
INT u;
{
   return(REF_Q1ROT(x,FD(pel->f[0],u),FD(pel->f[1],u),
                      FD(pel->f[2],u),FD(pel->f[3],u)));
}

#else

FLOAT sq1rot_value_ref(pel,x,u)
ELEMENT *pel; FLOAT *x; INT u;
{  eprintf("Error: sq1rot_value_ref not available.\n");  }

#endif

FLOAT sfunc_value_ref(pel,x,u,space)
ELEMENT *pel;
FLOAT *x;
INT u, space;
{
   switch(space){
   case P1_NC:   eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case P1_MOD:  eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case Q1ROT:   return(sq1rot_value_ref(pel,x,u));
        break;
   case P1C:     return(sp1c_value_ref(pel,x,u));
        break;
   case P1C_FBUB: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case P1C_NEW_FBUB: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case P1C_ELBUB: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case MINI_L_DIV_FR: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case Q1C:     return(sq1c_value_ref(pel,x,u));
        break;
   case P2C:     return(sp2c_value_ref(pel,x,u));
        break;
   case IP2C:    eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case P2C_ELBUB: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   case IP2C_ELBUB: eprintf("Error: sfunc_value_ref not available.\n");
                 return(0.);
        break;
   default:
        eprintf("Error: sfunc_value_ref not available.\n");
        break;
   }
}

FLOAT sfunc_value(pel,x,b,u,space)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2];
INT u, space;
{
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
   switch(space){
   case P1_NC:   return(sp1nc_value(fa0,fa1,fa2,x,b,u));
        break;
   case P1_MOD:  if (0) 
                    return(sp1nc_value(fa0,fa1,fa2,x,b,u));
                 else{
                    eprintf("Error: sfunc_value not available.\n");
                    return(0.);
                 }
        break;
   case Q1ROT:   eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case P1C:     return(sp1c_value(n0,n1,n2,x,b,u));
        break;
   case P1C_FBUB: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case P1C_NEW_FBUB: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case P1C_ELBUB: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case MINI_L_DIV_FR: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case P2C:     return(sp2c_value(n0,n1,n2,fa0,fa1,fa2,x,b,u));
        break;
   case IP2C:    eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case P2C_ELBUB: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   case IP2C_ELBUB: eprintf("Error: sfunc_value not available.\n");
                 return(0.);
        break;
   default:
        eprintf("Error: sfunc_value not available.\n");
        break;
   }
}

void vfunc_value(pel,x,b,u,space,val)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], *val;
INT u, space;
{
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
   switch(space){
   case P1_NC:     eprintf("Error: vfunc_value not available.\n");
        break;
   case P1_MOD:    eprintf("Error: vfunc_value not available.\n");
        break;
   case Q1ROT:     eprintf("Error: vfunc_value not available.\n");
        break;
   case P1C:       vp1c_value(n0,n1,n2,x,b,u,val);
        break;
   case P1C_FBUB:  eprintf("Error: vfunc_value not available.\n");
        break;
   case P1C_NEW_FBUB: eprintf("Error: vfunc_value not available.\n");
        break;
   case P1C_ELBUB: eprintf("Error: vfunc_value not available.\n");
        break;
   case MINI_L_DIV_FR: eprintf("Error: vfunc_value not available.\n");
        break;
   case P2C:       eprintf("Error: vfunc_value not available.\n");
        break;
   case IP2C:      eprintf("Error: vfunc_value not available.\n");
        break;
   case P2C_ELBUB: eprintf("Error: vfunc_value not available.\n");
        break;
   case IP2C_ELBUB: eprintf("Error: vfunc_value not available.\n");
        break;
   default:
        eprintf("Error: vfunc_value not available.\n");
        break;
   }
}

void sgrad_value(pel,x,b,u,space,grad)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], *grad;
INT u, space;
{
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
   switch(space){
   case P1_NC:   sp1nc_grad(fa0,fa1,fa2,b,u,grad);
        break;
   case P1_MOD:  if (0) 
                    sp1nc_grad(fa0,fa1,fa2,b,u,grad);
                 else
                    eprintf("Error: sgrad_value not available.\n");
        break;
   case Q1ROT:   eprintf("Error: sgrad_value not available.\n");
        break;
   case P1C:     sp1c_grad(n0,n1,n2,b,u,grad);
        break;
   case P1C_FBUB: eprintf("Error: sgrad_value not available.\n");
        break;
   case P1C_NEW_FBUB: eprintf("Error: sgrad_value not available.\n");
        break;
   case P1C_ELBUB: eprintf("Error: sgrad_value not available.\n");
        break;
   case MINI_L_DIV_FR: eprintf("Error: sgrad_value not available.\n");
        break;
   case P2C:     sp2c_grad(n0,n1,n2,fa0,fa1,fa2,x,b,u,grad);
        break;
   case IP2C:    eprintf("Error: sgrad_value not available.\n");
        break;
   case P2C_ELBUB: eprintf("Error: sgrad_value not available.\n");
        break;
   case IP2C_ELBUB: eprintf("Error: sgrad_value not available.\n");
        break;
   default:
        eprintf("Error: sgrad_value not available.\n");
        break;
   }
}

FLOAT sLapl_value(pel,x,b,u,space)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2];
INT u, space;
{
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
   switch(space){
   case P1_NC:   return(0.);
        break;
   case P1_MOD:  if (0) 
                    return(0.);
                 else{
                    eprintf("Error: sLapl_value not available.\n");
                    return(0.);
                 }
        break;
   case Q1ROT:   eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case P1C:     return(0.);
        break;
   case P1C_FBUB: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case P1C_NEW_FBUB: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case P1C_ELBUB: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case MINI_L_DIV_FR: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case P2C:     return(sp2c_Lapl(fa0,fa1,fa2,b,u));
        break;
   case IP2C:    eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case P2C_ELBUB: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   case IP2C_ELBUB: eprintf("Error: sLapl_value not available.\n");
                 return(0.);
        break;
   default:
        eprintf("Error: sLapl_value not available.\n");
        break;
   }
}

FLOAT scalar_conv_diff_res(pel,x,b,u,space,eps,bb0,bb1,react,rhs)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u, space;
{
   FLOAT grad_u[DIM], bb[DIM];
   
   bb[0] = bb0(x);
   bb[1] = bb1(x);
   sgrad_value(pel,x,b,u,space,grad_u);
   return(-eps*sLapl_value(pel,x,b,u,space) + DOT(bb,grad_u) + 
           react(x)*sfunc_value(pel,x,b,u,space) - rhs(x));
}

double L2norm2_of_cd_res(pelem,b,u,space,eps,bb0,bb1,react,rhs)
ELEMENT *pelem;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u, space;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], 
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pelem);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pelem,x1,  b,u,space,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pelem,x2,  b,u,space,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pelem,x3,  b,u,space,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pelem,x112,b,u,space,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pelem,x113,b,u,space,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pelem,x221,b,u,space,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pelem,x223,b,u,space,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pelem,x331,b,u,space,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pelem,x332,b,u,space,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pelem,x123,b,u,space,eps,bb0,bb1,react,rhs);
   return(p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)*VOLUME(pelem));
}

#if N_DATA & SCALAR_NODE_DATA

double cd_res_times_test_cd_res(pel,b,u,v,i,eps,bb0,bb1,react,rhs)
ELEMENT *pel;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u, v, i;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], 
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123,
         q1, q2, q3, q112, q113, q221, q223, q331, q332, q123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pel);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pel,x1,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pel,x2,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pel,x3,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pel,x112,b,u,P1C,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pel,x113,b,u,P1C,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pel,x221,b,u,P1C,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pel,x223,b,u,P1C,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pel,x331,b,u,P1C,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pel,x332,b,u,P1C,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pel,x123,b,u,P1C,eps,bb0,bb1,react,rhs);
   NDS(pel->n[0],v) = NDS(pel->n[1],v) = NDS(pel->n[2],v) = 0.;
   NDS(pel->n[i],v) = 1.;
   q1   = scalar_conv_diff_res(pel,x1,  b,v,P1C,eps,bb0,bb1,react,zero_func);
   q2   = scalar_conv_diff_res(pel,x2,  b,v,P1C,eps,bb0,bb1,react,zero_func);
   q3   = scalar_conv_diff_res(pel,x3,  b,v,P1C,eps,bb0,bb1,react,zero_func);
   q112 = scalar_conv_diff_res(pel,x112,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q113 = scalar_conv_diff_res(pel,x113,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q221 = scalar_conv_diff_res(pel,x221,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q223 = scalar_conv_diff_res(pel,x223,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q331 = scalar_conv_diff_res(pel,x331,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q332 = scalar_conv_diff_res(pel,x332,b,v,P1C,eps,bb0,bb1,react,zero_func);
   q123 = scalar_conv_diff_res(pel,x123,b,v,P1C,eps,bb0,bb1,react,zero_func);
   return(pq3(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel));
}

#else

double cd_res_times_test_cd_res(pel,b,u,v,i,eps,bb0,bb1,react,rhs)
ELEMENT *pel; FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(); INT u, v, i;
{  eprintf("Error: cd_res_times_test_cd_res not available.\n");  }

#endif

#if TAU_SPACE == P0

void cd_res_times_sd(pel,tot_der,b,u,v,eps,bb0,bb1,react,rhs)
ELEMENT *pel;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT tot_der, u, v;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], grad[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123,
         q1, q2, q3, q112, q113, q221, q223, q331, q332, q123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pel);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pel,x1,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pel,x2,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pel,x3,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pel,x112,b,u,P1C,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pel,x113,b,u,P1C,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pel,x221,b,u,P1C,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pel,x223,b,u,P1C,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pel,x331,b,u,P1C,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pel,x332,b,u,P1C,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pel,x123,b,u,P1C,eps,bb0,bb1,react,rhs);
   sp1c_grad(pel->n[0],pel->n[1],pel->n[2],b,v,grad);
   q1   = bb0(x1  )*grad[0] + bb1(x1  )*grad[1];
   q2   = bb0(x2  )*grad[0] + bb1(x2  )*grad[1];
   q3   = bb0(x3  )*grad[0] + bb1(x3  )*grad[1];
   q112 = bb0(x112)*grad[0] + bb1(x112)*grad[1];
   q113 = bb0(x113)*grad[0] + bb1(x113)*grad[1];
   q221 = bb0(x221)*grad[0] + bb1(x221)*grad[1];
   q223 = bb0(x223)*grad[0] + bb1(x223)*grad[1];
   q331 = bb0(x331)*grad[0] + bb1(x331)*grad[1];
   q332 = bb0(x332)*grad[0] + bb1(x332)*grad[1];
   q123 = bb0(x123)*grad[0] + bb1(x123)*grad[1];
   ED(pel,tot_der) = -pq3(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
}

#elif TAU_SPACE == P1C

void cd_res_times_sd(pel,tot_der,b,u,v,eps,bb0,bb1,react,rhs)
ELEMENT *pel;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT tot_der, u, v;
{
   NODE *n0, *n1, *n2;
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], grad[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123,
         q1, q2, q3, q112, q113, q221, q223, q331, q332, q123;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   VERTICES_OF_ELEMENT(x1,x2,x3,pel);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pel,x1,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pel,x2,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pel,x3,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pel,x112,b,u,P1C,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pel,x113,b,u,P1C,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pel,x221,b,u,P1C,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pel,x223,b,u,P1C,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pel,x331,b,u,P1C,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pel,x332,b,u,P1C,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pel,x123,b,u,P1C,eps,bb0,bb1,react,rhs);
   sp1c_grad(pel->n[0],pel->n[1],pel->n[2],b,v,grad);
   q1   = bb0(x1  )*grad[0] + bb1(x1  )*grad[1];
   q2   = bb0(x2  )*grad[0] + bb1(x2  )*grad[1];
   q3   = bb0(x3  )*grad[0] + bb1(x3  )*grad[1];
   q112 = bb0(x112)*grad[0] + bb1(x112)*grad[1];
   q113 = bb0(x113)*grad[0] + bb1(x113)*grad[1];
   q221 = bb0(x221)*grad[0] + bb1(x221)*grad[1];
   q223 = bb0(x223)*grad[0] + bb1(x223)*grad[1];
   q331 = bb0(x331)*grad[0] + bb1(x331)*grad[1];
   q332 = bb0(x332)*grad[0] + bb1(x332)*grad[1];
   q123 = bb0(x123)*grad[0] + bb1(x123)*grad[1];
   NDS(n0,tot_der)=-pq3(p1,0.,0.,2.*p112/3.,2.*p113/3.,p221/3.,0.,p331/3.,0.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
   NDS(n1,tot_der)=-pq3(0.,p2,0.,p112/3.,0.,2.*p221/3.,2.*p223/3.,0.,p332/3.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
   NDS(n2,tot_der)=-pq3(0.,0.,p3,0.,p113/3.,0.,p223/3.,2.*p331/3.,2.*p332/3.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
}

#elif TAU_SPACE == P1_NC

void cd_res_times_sd(pel,tot_der,b,u,v,eps,bb0,bb1,react,rhs)
ELEMENT *pel;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT tot_der, u, v;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], grad[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123,
         q1, q2, q3, q112, q113, q221, q223, q331, q332, q123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pel);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pel,x1,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pel,x2,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pel,x3,  b,u,P1C,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pel,x112,b,u,P1C,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pel,x113,b,u,P1C,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pel,x221,b,u,P1C,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pel,x223,b,u,P1C,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pel,x331,b,u,P1C,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pel,x332,b,u,P1C,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pel,x123,b,u,P1C,eps,bb0,bb1,react,rhs);
   sp1c_grad(pel->n[0],pel->n[1],pel->n[2],b,v,grad);
   q1   = bb0(x1  )*grad[0] + bb1(x1  )*grad[1];
   q2   = bb0(x2  )*grad[0] + bb1(x2  )*grad[1];
   q3   = bb0(x3  )*grad[0] + bb1(x3  )*grad[1];
   q112 = bb0(x112)*grad[0] + bb1(x112)*grad[1];
   q113 = bb0(x113)*grad[0] + bb1(x113)*grad[1];
   q221 = bb0(x221)*grad[0] + bb1(x221)*grad[1];
   q223 = bb0(x223)*grad[0] + bb1(x223)*grad[1];
   q331 = bb0(x331)*grad[0] + bb1(x331)*grad[1];
   q332 = bb0(x332)*grad[0] + bb1(x332)*grad[1];
   q123 = bb0(x123)*grad[0] + bb1(x123)*grad[1];
   EDSN(pel,tot_der,0)=-pq3(p1,0.,0.,2.*p112/3.,2.*p113/3.,p221/3.,0.,p331/3.,0.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
   EDSN(pel,tot_der,1)=-pq3(0.,p2,0.,p112/3.,0.,2.*p221/3.,2.*p223/3.,0.,p332/3.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
   EDSN(pel,tot_der,2)=-pq3(0.,0.,p3,0.,p113/3.,0.,p223/3.,2.*p331/3.,2.*p332/3.,p123/3.,
              q1,q2,q3,q112,q113,q221,q223,q331,q332,q123)*VOLUME(pel);
}

#else

void cd_res_times_sd(pel,tot_der,b,u,v,eps,bb0,bb1,react,rhs)
ELEMENT *pel; FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT tot_der, u, v;
{  eprintf("Error: cd_res_times_sd not available.\n");  }

#endif

void L2norms_of_cd_res_and_grad_u(pelem,b,u,space,eps,bb0,bb1,react,rhs,res,grad)
ELEMENT *pelem;
FLOAT b[DIM2][DIM2], eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), *res, *grad;
INT u, space;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], q1[DIM], q2[DIM], q3[DIM], q112[DIM], 
         q113[DIM], q221[DIM], q223[DIM], q331[DIM], q332[DIM], q123[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pelem);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = scalar_conv_diff_res(pelem,x1,  b,u,space,eps,bb0,bb1,react,rhs);
   p2   = scalar_conv_diff_res(pelem,x2,  b,u,space,eps,bb0,bb1,react,rhs);
   p3   = scalar_conv_diff_res(pelem,x3,  b,u,space,eps,bb0,bb1,react,rhs);
   p112 = scalar_conv_diff_res(pelem,x112,b,u,space,eps,bb0,bb1,react,rhs);
   p113 = scalar_conv_diff_res(pelem,x113,b,u,space,eps,bb0,bb1,react,rhs);
   p221 = scalar_conv_diff_res(pelem,x221,b,u,space,eps,bb0,bb1,react,rhs);
   p223 = scalar_conv_diff_res(pelem,x223,b,u,space,eps,bb0,bb1,react,rhs);
   p331 = scalar_conv_diff_res(pelem,x331,b,u,space,eps,bb0,bb1,react,rhs);
   p332 = scalar_conv_diff_res(pelem,x332,b,u,space,eps,bb0,bb1,react,rhs);
   p123 = scalar_conv_diff_res(pelem,x123,b,u,space,eps,bb0,bb1,react,rhs);
   *res = sqrt(p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)*VOLUME(pelem));
   sgrad_value(pelem,x1,  b,u,space,q1);
   sgrad_value(pelem,x2,  b,u,space,q2);
   sgrad_value(pelem,x3,  b,u,space,q3);
   sgrad_value(pelem,x112,b,u,space,q112);
   sgrad_value(pelem,x113,b,u,space,q113);
   sgrad_value(pelem,x221,b,u,space,q221);
   sgrad_value(pelem,x223,b,u,space,q223);
   sgrad_value(pelem,x331,b,u,space,q331);
   sgrad_value(pelem,x332,b,u,space,q332);
   sgrad_value(pelem,x123,b,u,space,q123);
   *grad = sqrt((p32(q1[0],q2[0],q3[0],q112[0],q113[0],q221[0],q223[0],
                               q331[0],q332[0],q123[0]) +
                 p32(q1[1],q2[1],q3[1],q112[1],q113[1],q221[1],q223[1],
	                       q331[1],q332[1],q123[1]))*VOLUME(pelem));
}

FLOAT L2norm_of_u(pelem,b,u,space)
ELEMENT *pelem;
FLOAT b[DIM2][DIM2];
INT u, space;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
     x332[DIM], x123[DIM], p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;

   VERTICES_OF_ELEMENT(x1, x2, x3, pelem);
   points(x1, x2, x3, x112, x113, x221, x223, x331, x332, x123);
   p1   = sfunc_value(pelem,x1,  b,u,space);
   p2   = sfunc_value(pelem,x2,  b,u,space);
   p3   = sfunc_value(pelem,x3,  b,u,space);
   p112 = sfunc_value(pelem,x112,b,u,space);
   p113 = sfunc_value(pelem,x113,b,u,space);
   p221 = sfunc_value(pelem,x221,b,u,space);
   p223 = sfunc_value(pelem,x223,b,u,space);
   p331 = sfunc_value(pelem,x331,b,u,space);
   p332 = sfunc_value(pelem,x332,b,u,space);
   p123 = sfunc_value(pelem,x123,b,u,space);
   return(sqrt(p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)*VOLUME(pelem)));
}

void set_directions_of_edges(pel,dir)
ELEMENT *pel;
INT *dir;
{
   INT i, j;

   dir[0] = dir[1] = dir[2] = dir[3] = 1;
   if (DIRECTION_OF_EDGES == YES){
      if (ELEMENT_TYPE == SIMPLEX)
         dir[1] = 1;
      else if (ELEMENT_TYPE == CUBE)
         for (i=0; i < NVERT; i++){
            j = i + 1;
            if (j == NVERT)
               j = 0;
            if (pel->n[i]->index > pel->n[j]->index)
               dir[i] = 0;
         }
   }
}

INT f_i(ind,nn,dir)
INT ind, nn, dir;
{
   if (nn == 1 || dir)
      return(ind);
   else
      return(nn - 1 - ind);
}

NODE *ltop_node();
FACE *ltop_face();

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

FLOAT gfunc_value_ref(tGrid,pel,u,x,finite_el)
GRID *tGrid;
ELEMENT *pel;
FLOAT x[DIM];
INT u;
FINITE_ELEMENT finite_el;
{
   NODE *pnode;
   FACE *pface;
   DOUBLE_FUNC *nq, *fq, *eq;
   FLOAT val=0.;
   INT i, k, m, kk, nn, mm, dir[4];

   nq = finite_el.nbasis;
   fq = finite_el.fbasis;
   eq = finite_el.ebasis;
   kk = finite_el.k;
   nn = finite_el.n;
   mm = finite_el.m;
   if (kk)
      for (i=0; i < NVERT; i++){
         pnode = ltop_node(pel->n[i],tGrid);
         m = i*kk;
         for (k = m; k < m+kk; k++)
            val += NDMV(pnode,u,k-m)*(nq[k])(x);
      }
   if (nn){
      set_directions_of_edges(pel,dir);
      for (i=0; i < NVERT; i++){
         pface = ltop_face(pel->f[i],tGrid);
         m = i*nn;
         for (k = m; k < m+nn; k++) 
            val += FDMV(pface,u,f_i(k-m,nn,dir[i]))*(fq[k])(x);
      }
   }
   for (i = 0; i < mm; i++)
      val += EDMV(pel,u,i)*(eq[i])(x);
   return(val);
}

void ggrad_value_ref(tGrid,pel,u,x,finite_el,grad)
GRID *tGrid;
ELEMENT *pel;
FLOAT x[DIM], grad[DIM];
INT u;
FINITE_ELEMENT finite_el;
{
   NODE *pnode;
   FACE *pface;
   DOUBLE_FUNC *nq0, *nq1, *fq0, *fq1, *eq0, *eq1;
   INT i, k, m, kk, nn, mm, dir[4];

   SET7(grad,0.)
   nq0 = finite_el.nbasis_0;
   nq1 = finite_el.nbasis_1;
   fq0 = finite_el.fbasis_0;
   fq1 = finite_el.fbasis_1;
   eq0 = finite_el.ebasis_0;
   eq1 = finite_el.ebasis_1;
   kk = finite_el.k;
   nn = finite_el.n;
   mm = finite_el.m;
   if (kk)
      for (i=0; i < NVERT; i++){
         pnode = ltop_node(pel->n[i],tGrid);
         m = i*kk;
         for (k = m; k < m+kk; k++){
            grad[0] += NDMV(pnode,u,k-m)*(nq0[k])(x);
            grad[1] += NDMV(pnode,u,k-m)*(nq1[k])(x);
         }
      }
   if (nn){
      set_directions_of_edges(pel,dir);
      for (i=0; i < NVERT; i++){
         pface = ltop_face(pel->f[i],tGrid);
         m = i*nn;
         for (k = m; k < m+nn; k++){
            grad[0] += FDMV(pface,u,f_i(k-m,nn,dir[i]))*(fq0[k])(x);
            grad[1] += FDMV(pface,u,f_i(k-m,nn,dir[i]))*(fq1[k])(x);
         }
      }
   }
   for (i = 0; i < mm; i++){
      grad[0] += EDMV(pel,u,i)*(eq0[i])(x);
      grad[1] += EDMV(pel,u,i)*(eq1[i])(x);
   }
}

#else

FLOAT gfunc_value_ref(tGrid,pel,u,x,finite_el)
GRID *tGrid; ELEMENT *pel; FLOAT x[DIM]; INT u; FINITE_ELEMENT finite_el;
{  eprintf("Error: gfunc_value_ref not available.\n");  }

void ggrad_value_ref(tGrid,pel,u,x,finite_el,grad)
GRID *tGrid; ELEMENT *pel; FLOAT x[DIM], grad[DIM]; INT u; FINITE_ELEMENT finite_el;
{  eprintf("Error: ggrad_value_ref not available.\n");  }

#endif

/* value of grad(u) at the point F(x_ref) */
void ggrad_value(tGrid,pel,u,x_ref,finite_el,ref_map,ggrad)
GRID *tGrid;
ELEMENT *pel;
FLOAT x_ref[DIM], ggrad[DIM];
INT u;
FINITE_ELEMENT finite_el;
REF_MAPPING ref_map;
{
   FLOAT b[DIM][DIM], grad[DIM], z;

   z = jacobian(x_ref,&ref_map);
   inv_of_jac_matr_times_jacobian(x_ref,b,&ref_map);
   ggrad_value_ref(tGrid,pel,u,x_ref,finite_el,grad);
   ggrad[0] = (grad[0]*b[0][0] + grad[1]*b[1][0])/z;
   ggrad[1] = (grad[0]*b[0][1] + grad[1]*b[1][1])/z;
}

/* value of Res(u) at the point F(x_ref) */
FLOAT gscalar_conv_diff_res(tGrid,pel,x_ref,u,eps,bb0,bb1,react,rhs,
                                                              finite_el,ref_map)
GRID *tGrid;
ELEMENT *pel;
FLOAT *x_ref, eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)();
INT u;
FINITE_ELEMENT finite_el;
REF_MAPPING ref_map;
{
   FLOAT grad_u[DIM], bb[DIM], x[DIM];
   
   ref_map_point(x_ref,x,&ref_map);
   bb[0] = bb0(x);
   bb[1] = bb1(x);
   ggrad_value(tGrid,pel,u,x_ref,finite_el,ref_map,grad_u);
   return(
//-eps*gLapl_value + 
      DOT(bb,grad_u) + 
      react(x)*gfunc_value_ref(tGrid,pel,u,x_ref,finite_el) - rhs(x));
}

void gL2norms_of_cd_res_and_grad_u(pel,u,space,eps,bb0,bb1,react,rhs,res,grad)
ELEMENT *pel; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), *res, *grad; INT u, space;
{  eprintf("Error: gL2norms_of_cd_res_and_grad_u not available.\n");  }

FLOAT gL2norm_of_u(pelem,u,space)
ELEMENT *pelem; INT u, space;
{  eprintf("Error: gL2norm_of_u not available.\n");  }

INT approx_order(space)
INT space;
{
   INT order=1;

   if (space == P1C || space == IP1C || space == P1C_FBUB || space == GP1C ||
       space == GP1X3C || space == P1C_NEW_FBUB || space == P1C_FBUB_LG || 
       space == P1C_ELBUB || space == GP1C_ELBUB ||
       space == P1_NC || space == P1_MOD || space == Q1C || space == GQ1C || 
       space == GQ1X4C || space == GQ1C_ELBUB || space == Q1ROT)
      order = 1;
   else if (space == P2C || space == IP2C || space == P2C_ELBUB || 
            space == IP2C_ELBUB || space == GP2C || space == GP2X3C || 
            space == GP2C_3ELBUB || space == GP2C_6ELBUB || 
            space == GQ2C || space == GQ2X4C || space == GQ2B3C ||
            space == GQ2C_2ELBUB || space == GQ2C_3ELBUB)
      order = 2;
   else if (space == P3C || space == Q3C || space == GP3C || space == GQ3C)
      order = 3;
   else if (space == P4C || space == Q4C || space == GP4C || space == GQ4C)
      order = 4;
   else
      eprintf("Error: approx_order not defined for space used.\n");
   return(order);
}

#else /*  DIM == 3  */

FINITE_ELEMENT p1c_element = { SIMPLEX, 1, 0, 0, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL };

FLOAT fcn_ref_map_value(x,ref_map,f)
FLOAT *x, (*f)(); REF_MAPPING *ref_map;
{  eprintf("Error: fcn_ref_map_value not available.\n"); return(0.);  }

FLOAT gfunc_value_ref(tGrid,pel,u,x,finite_el)
GRID *tGrid; ELEMENT *pel; FLOAT x[DIM]; INT u; FINITE_ELEMENT finite_el;
{  eprintf("Error: gfunc_value_ref not available.\n"); return(0.);  }

void ggrad_value_ref(tGrid,pel,u,x,finite_el,grad)
GRID *tGrid; ELEMENT *pel; FLOAT x[DIM], grad[DIM]; INT u; FINITE_ELEMENT finite_el;
{  eprintf("Error: ggrad_value_ref not available.\n");  }

#endif

