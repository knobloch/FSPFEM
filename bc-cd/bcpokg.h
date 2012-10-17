#define CC0   1.   /*   min( react(x) - 0.5*div_bb(x) ); if min < 0, CC0=0  */

/* t0 */
#define KOEF_116  -5.
#define KOEF_117   7.
#define KOEF_118   3.
#define KOEF_120  -2.
/* recat */
#define KOEF_110  0.
#define KOEF_111  0.
#define KOEF_112  0.
#define KOEF_114  0.  
/* flow field bb0 */
#define KOEF_210 0.
#define KOEF_211 0.
#define KOEF_212 0.
#define KOEF_214 0.
/* flow field bb1 */
#define KOEF_216 0.
#define KOEF_217 0.
#define KOEF_218 0.
#define KOEF_220 0.

#if DIM == 3

#else  /* if DIM == 2 */

 double bb0(x)  /* first component of the Dirichlet boundary condition */
 double x[2];
 {
    return( KOEF_210 + KOEF_211*x[0] + KOEF_212*x[1] + KOEF_214*x[0]*x[1]);
 }
 
 double bb1(x)  /*  second component of the Dirichlet boundary condition */
 double x[2];
 {
    return( KOEF_216 + KOEF_217*x[0] + KOEF_218*x[1] + KOEF_220*x[0]*x[1]);
 }
 
 double div_bb(x)
 double x[DIM];
 {
    return( KOEF_211 + KOEF_214*x[1] + KOEF_218 + KOEF_220*x[0]);
 }

 double react(x)  /* first component of the Dirichlet boundary condition */
 double x[2];
 {
    return( KOEF_110 + KOEF_111*x[0] + KOEF_112*x[1] + KOEF_114*x[0]*x[1]);
 }
 
 double t0(x)  /*  second component of the Dirichlet boundary condition */
 double x[2];
 {
    return( KOEF_116 + KOEF_117*x[0] + KOEF_118*x[1] + KOEF_120*x[0]*x[1]);
 }
 
 double t01(x)
 double x[2];
 {
    return( KOEF_117 + KOEF_120*x[1] ); 
 }
 
 double t02(x)
 double x[2];
 {
    return( KOEF_118 + KOEF_120*x[0]);
 }
 
 double t03(x)
 double x[2];
 {
    return( 0.0 );
 }
 
 double ft0(x)
 double x[2];
 {
    return( bb0(x)*t01(x) + bb1(x)*t02(x) + react(x)*t0(x) );
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
    return(0.);
 }
 
 double u02(x)  /*  second component of the Dirichlet boundary condition */
 double x[2];
 {
    return(0.);
 }
 
 double f01(x)  /* first component of the right hand sight */
 double x[2];
 {
    return(0.);
 }

 
 double f02(x)  /* second component of the right hand sight */
 double x[2];
 {
    return(0.);
 }
 
 double p0(x)    /* pressure */
 double x[2];
 {
    return(0.);
 }

 double u011(x)  
 double x[2];
 {
    return(0.);
 }

 double u012(x)  
 double x[2];
 {
    return(0.);
 }

 double u021(x)  
 double x[2];
 {
    return(0.);
 }

 double u022(x)  
 double x[2];
 {
    return(0.);
 }

#endif
