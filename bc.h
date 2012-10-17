 double u01(x)  /* first component of the Dirichlet boundary condition */
 double x[3];
 {
 return( x[0]*x[0]*x[0]+0.1);
 }
 
 double u02(x)  /*  second component of the Dirichlet boundary condition */
 double x[3];
 {
   return( x[1]*x[1]*x[1]+0.1 );
 }
 
 double u03(x)  /*  third component of the Dirichlet boundary condition */
 double x[3];
 {
   return( x[2]*x[2]*x[2]+0.1 );
 }
 
 
 double f01(x)  /* first component of the right hand sight */
 double x[3];
 {
    return(-6.0*x[0]);
 }
 
  double f02(x)  /* second component of the right hand sight */
 double x[3];
 {
  return(-6.0*x[1]);
 }
  double f03(x)  /* third component of the right hand sight */
 double x[3];
 {
  return(-6.0*x[2]);
 }
 
 double u011(x)  
 double x[3];
 {
 return( 3.0*x[0]*x[0] );
 }
 double u012(x)  
 double x[3];
 {
 return( 0.0 );
 }
 double u013(x)  
 double x[3];
 {
 return( 0.0 );
 }

 double u021(x)  
 double x[3];
 {
   return( 0.0 );
 }
 double u022(x)  
 double x[3];
 {
   return( 3.0*x[1]*x[1] );
 }
 double u023(x)  
 double x[3];
 {
   return( 0.0 );
 }
 
 double u031(x)  
 double x[3];
 {
   return( 0.0 );
 }
 double u032(x)  
 double x[3];
 {
   return( 0.0 );
 }
 double u033(x)  
 double x[3];
 {
   return( 3.0*x[2]*x[2] );
 }
 
 double p0(x)  
 double x[3];
 {
   return( 3.0*x[2]*x[2] );
 }
 
 double t0(x)
 double x[3];
 {
    if (x[2] > 0.62)
       return(-0.5);
    else
       return( (-3. + 40.*(x[0]*x[0] + x[1]*x[1]))/74. );
 }
 
 double t01(x)
 double x[3];
 {
    return( 0. );
 }
 
 double t02(x)
 double x[3];
 {
    return( 0. );
 }
 
 double t03(x)
 double x[3];
 {
    return( 0. );
 }
 
 double ft0(x)
 double x[3];
 {
   return( 0. );
 }
 
 double nft0(x)
 double x[3];
 {
   return( TNU*2.923882 );
 }
 
 double fi1(x)
 double x;
 {
   return( TNU*(1.30472e-6*x*x*x + 2.39998e-4*x*x + 0.0165550*x + 0.507537) );
 }
