/*  Example 3 from JMT97  */

#define CC0   1.   /*   min( react(x) - 0.5*div_bb(x) ); if min < 0, CC0=0  */

#if DIM == 3

#else  /* if DIM == 2 */

 double bb0(x)
 double x[2];
 {
    return(2.);
 }

 double bb1(x)
 double x[2];
 {
    return(3.);
 }

 double div_bb(x)
 double x[DIM];
 {
    return(0.);
 }

 double react(x)  /*  function from the reactive term  */
 double x[DIM];
 {
    return( 0. );
 }
 
 double exp_f1(x)
 double x[2];
 {
    return( exp( 2.*(x[0]-1.)/SNU ) );
 }
 
 double exp_f2(x)
 double x[2];
 {
    return( exp( 3.*(x[1]-1.)/SNU ) );
 }
 
 double exp_f3(x)
 double x[2];
 {
    return( exp( ( 2.*(x[0]-1.) + 3.*(x[1]-1.) )/SNU ) );
 }
 
 double t0(x)
 double x[2];
 {
    return( x[0]*x[1]*x[1] - x[1]*x[1]*exp_f1(x) - x[0]*exp_f2(x) + exp_f3(x) );
 }
 
 double t01(x)
 double x[2];
 {
    return( x[1]*x[1] - 2./SNU*x[1]*x[1]*exp_f1(x) - exp_f2(x) + 
                                                             2./SNU*exp_f3(x) );
 }
 
 double t02(x)
 double x[2];
 {
    return( 2.*x[0]*x[1] - 2.*x[1]*exp_f1(x) - 3./SNU*x[0]*exp_f2(x) + 
                                                             3./SNU*exp_f3(x) );
 }
 
 double t03(x)
 double x[2];
 {
    return( 0.0 );
 }
 
 double ft0(x)
 double x[2];
 {
    return( -SNU*( 2.*x[0] - (2.+4./SNU/SNU*x[1]*x[1])*exp_f1(x) - 
                   9./SNU/SNU*x[0]*exp_f2(x) + 13./SNU/SNU*exp_f3(x) ) +
             bb0(x)*t01(x) + bb1(x)*t02(x) + react(x)*t0(x) );
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
