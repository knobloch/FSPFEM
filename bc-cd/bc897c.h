#define CC0   0.   /*   min( react(x) - 0.5*div_bb(x) ); if min < 0, CC0=0  */

#if DIM == 3

#else  /* if DIM == 2 */

 double bb0(x)
 double x[2];
 {
    return(1.);
 }

 double bb1(x)
 double x[2];
 {
    return(-0.0001);
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
 
 double t0(x)
 double x[2];
 {
    if (x[0] > 1.e-6 && x[0] < 0.999999 && x[1] > 1.e-6 && x[1] < 0.999999)
       return(x[0]);
    else
       return(0.);
 }
 
 double t01(x)
 double x[2];
 {
    return(0.);
 }
 
 double t02(x)
 double x[2];
 {
    return(0.);
 }
 
 double t03(x)
 double x[2];
 {
    return(0.);
 }
 
 double ft0(x)
 double x[2];
 {
   return(1.);
 }
 
 double nft0(x)
 double x[2];
 {
    return(0.);
 }
 
 double fi1(x)
 double x;
 {
    return(0.);
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
