#define K0   1.
#define K1  -1.
#define K2   1.
#define K3   3.
#define K4  -2.
#define K5   8.
#define M0  -1.
#define M1  -2.
#define M2   3.
#define M3  -3.
#define M4   4.
#define M5  -1.
#define L0   2.
#define L1  -5.
#define L2   3.
#define L3 -15.
#define L4  -5.
#define L5   3.
#define R0  -0.e3
#define R1   0.e3
#define R2  -0.e3
#define R3  -3.e3
#define R4   2.e3
#define R5   1.e3

#define CC0   1.

#if DIM == 3

#else  /* if DIM == 2 */

 double bb0(x)
 double x[2];
 {
    return(M0*x[0]*x[0] + M1*x[1]*x[1] + M2*x[0]*x[1] +M3*x[0] + M4*x[1] + M5);
 }
 
 double bb1(x)
 double x[2];
 {
    return(L0*x[0]*x[0] + L1*x[1]*x[1] + L2*x[0]*x[1] +L3*x[0] + L4*x[1] + L5);
 }
 
 double div_bb(x)
 double x[DIM];
 {
    return(2.*M0*x[0] + M2*x[1] + M3 + 2.*L1*x[1] + L2*x[0] + L4);
 }

 double react(x)  /*  function from the reactive term  */
 double x[DIM];
 {
    return(R0*x[0]*x[0] + R1*x[1]*x[1] + R2*x[0]*x[1] + R3*x[0] + R4*x[1] + R5);
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
 
 double ft0(x)
 double x[2];
 {
    return( -TNU*2.*(K0+K1) + bb0(x)*t01(x) + bb1(x)*t02(x) + react(x)*t0(x) );
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

#endif
