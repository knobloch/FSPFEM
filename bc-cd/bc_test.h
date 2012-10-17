#define K0   1.
#define K1   -1.
#define K2   3.
#define K3   1.
#define K4   5.
#define K5   -7.
#define L0   3.
#define L1   1.
#define L2   -2.
#define M0   5.
#define M1   1.
#define M2   8.
#define N0   3.
#define N1   -7.
#define N2   7.

#define CC0   1.

#if DIM == 3

#else  /* if DIM == 2 */

 double div_bb(x)
 double x[DIM];
 {
    return(0.);
 }

 double react(x)  /*  function from the reactive term  */
 double x[DIM];
 {
    return( 1. );
 }

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
  return(-((12.*x[0]*x[0]-12.*x[0]+2.)*x[1]*(x[1]-1.)*(4.*x[1]-2.)+
               x[0]*x[0]*(x[0]-1.)*(x[0]-1.)*(24.*x[1]-12.)) +
         3.*x[0]*x[0]);
 }

 
 double f02(x)  /* second component of the right hand sight */
 double x[2];
 {
 return(((24.*x[0]-12.)*x[1]*x[1]*(x[1]-1.)*(x[1]-1.)+
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
    return(K0*x[0]*x[0] + K1*x[1]*x[1] + K2*x[0]*x[1] + K3*x[0] + K4*x[1] + K5);
 }
 
 double t01(x)
 double x[2];
 {
    return(2.*K0*x[0] + K2*x[1] + K3);
 }
 
 double t02(x)
 double x[2];
 {
    return(2.*K1*x[1] + K2*x[0] + K4);
 }
 
 double t03(x)
 double x[2];
 {
    return( 0.0 );
 }
 
 double q0(x)
 double x[2];
 {
    return(L0*x[0] + L1*x[1] + L2);
 }
 
 double bb0(x)
 double x[2];
 {
    return(M0*x[0] + M1*x[1] + M2);
 }
 
 double bb1(x)
 double x[2];
 {
    return(N0*x[0] + N1*x[1] + N2);
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
