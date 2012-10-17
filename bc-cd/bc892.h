#define CC0   0.   /*   min( react(x) - 0.5*div_bb(x) ); if min < 0, CC0=0  */

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
 
 double t0(x)
 double x[2];
 {
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
 
 double ft0(z)
 double z[2];
 {
  double x,y,d,r0,x0,y0,rxy,cxy,hxy,g1xy,g1xy_x,g1xy_xx,g1xy_y,g1xy_yy;
  double scale,g2xy,g2xy_x,g2xy_xx,g2xy_y,g2xy_yy,gxy,gxy_x,gxy_xx,gxy_y,gxy_yy;
  double gxy_l;

  d   =  2.0 * sqrt(1.0/TNU);
  r0  =  0.25;
  x0  =  0.5;
  y0  =  0.5;
  
    x = z[0];
    y = z[1];

    rxy =  (x-x0)*(x-x0) + (y-y0)*(y-y0);
    cxy =  d*(r0*r0 - rxy);
    hxy =  2.0*d / (1 + cxy*cxy);
    g1xy    =  atan(cxy) + PI/2.0;
    g1xy_x  = -hxy*(x-x0);
    g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
    g1xy_y  = -hxy*(y-y0);
    g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
    g1xy    =  g1xy / PI; 
    g1xy_x  =  g1xy_x / PI;
    g1xy_xx =  g1xy_xx / PI;
    g1xy_y  =  g1xy_y / PI;
    g1xy_yy =  g1xy_yy / PI;
    
    scale   =  16.0;
    g2xy    =  scale*x*(x - 1)*y*(y - 1);
    g2xy_x  =  scale*(2*x - 1)*y*(y - 1);
    g2xy_xx =  scale*2*y*(y - 1);
    g2xy_y  =  scale*(2*y - 1)*x*(x - 1);
    g2xy_yy =  scale*2*x*(x - 1);

    gxy     =  g1xy*g2xy;
    gxy_x   =  g1xy_x*g2xy + g1xy*g2xy_x;
    gxy_xx  =  g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx;
    gxy_y   =  g1xy_y*g2xy + g1xy*g2xy_y;
    gxy_yy  =  g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy;

    gxy_l  =  gxy_xx + gxy_yy;
  
    return(-TNU*gxy_l+bb0(z)*gxy_x+bb1(z)*gxy_y+react(z)*gxy);
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
