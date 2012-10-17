#if DIM == 3

#else  /* if DIM == 2 */

 double u01(x)  /* first component of the Dirichlet boundary condition */
 double x[2];
 {
    return((4.*x[1]-2.)*x[0]*x[0]*(x[0]-1.)*(x[0]-1.)*x[1]*(x[1]-1.));
 }
 
 double u02(x)  /*  second component of the Dirichlet boundary condition */
 double x[2];
 {
    return(-(4.*x[0]-2.)*x[0]*(x[0]-1.)*x[1]*x[1]*(x[1]-1.)*(x[1]-1.));
 }
 
 double f01(x)  /* first component of the right hand sight */
 double x[2];
 {
  return(-NU*((12.*x[0]*x[0]-12.*x[0]+2.)*x[1]*(x[1]-1.)*(4.*x[1]-2.)+
               x[0]*x[0]*(x[0]-1.)*(x[0]-1.)*(24.*x[1]-12.)) +
         3.*x[0]*x[0]);
 }

 
 double f02(x)  /* second component of the right hand sight */
 double x[2];
 {
 return(NU*((24.*x[0]-12.)*x[1]*x[1]*(x[1]-1.)*(x[1]-1.)+
            (4.*x[0]-2.)*x[0]*(x[0]-1.)*(12.*x[1]*x[1]-12.*x[1]+2.)) +
	3.*x[1]*x[1] ); 
 }
 
 double p0(x)    /* pressure */
 double x[2];
 {
    return( x[0]*x[0]*x[0] + x[1]*x[1]*x[1] - 0.5);
 }

 double u011(x)  
 double x[2];
 {
    return(x[0]*(x[0]-1.)*(4.*x[0]-2.)*x[1]*(x[1]-1.)*(4.*x[1]-2.));
 }

 double u012(x)  
 double x[2];
 {
    return(x[0]*x[0]*(x[0]-1.)*(x[0]-1.)*(12.*x[1]*x[1]-12.*x[1]+2.));
 }

 double u021(x)  
 double x[2];
 {
    return(-(12.*x[0]*x[0]-12.*x[0]+2.)*x[1]*x[1]*(x[1]-1.)*(x[1]-1.));
 }

 double u022(x)  
 double x[2];
 {
    return(-(4.*x[0]-2.)*x[0]*(x[0]-1.)*(4.*x[1]-2.)*x[1]*(x[1]-1.));
 }

 double t0(x)
 double x[2];
 {
    return( 7.*x[0] - 2.*x[1] + 0.1 );
 }
 
 double t01(x)
 double x[2];
 {
    return( 2.0*x[0] + x[1] );
 }
 
 double t02(x)
 double x[2];
 {
    return( x[1]*x[1] );
 }
 
 double t03(x)
 double x[2];
 {
    return( 0.0 );
 }
 
 double ft0(x)
 double x[2];
 {
    return( 0.0 );
 }
 
 double nft0(x)
 double x[2];
 {
    return( -2.0*x[0] );
 }
 
 double fi1(x)
 double x;
 {
    return( x*x*x/3 );
 }

#endif
