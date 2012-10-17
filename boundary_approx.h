/******************************************************************************/
/*                                                                            */
/*                       approximation of the boundary                        */
/*                                                                            */
/******************************************************************************/

void eprintf();

#if (DATA_STR & BOUNDARY_APPROX) && (DIM == 3)

INT NW, NLG, NLS;
FLOAT RR[5300], ZZ[5300];

FLOAT bf1(r)
FLOAT r;
{
   return(5.-sqrt(25.-r*r));
}

FLOAT bf2(z)
FLOAT z;
{
   z = z - 0.5999825;
   return(0.57596 - sqrt(0.02-z*z) - 0.07*log(5.*z/(sqrt(0.02)+sqrt(0.02-z*z))));
}

void boundary_points()  /* the boundary has to contain the point [0,0,0] */
{
   FLOAT d, eps, r, r1, rr, z, z1, zz;
   INT i=1, nlg;
   
   d = eps = 0.0005;
   rr = 1. - eps;
   zz = bf1(1.) - eps;
   RR[0] = r = 0.;
   ZZ[0] = z = 0.;
   while(r < rr || z < zz){
      r1 = r;
      z1 = z;
      r = r1 + d;
      while ((z=bf1(r)) > z1 + eps){
         d *= 0.7;
         r = r1 + d;
      }
      RR[i] = r;
      ZZ[i] = z;
      i++;
   }
   RR[i] = 1.;
   ZZ[i] = bf1(1.);
   z = ZZ[i] + eps;
   i++;
   while(z < 0.6){
      RR[i] = 1.;
      ZZ[i] = z;
      z += eps;
      i++;
   }
   RR[i] = 1.;
   ZZ[i] = 0.6;
   i++;
   d = 0.0005;
   rr = 1. - eps;
   zz = 0.6 + eps;
   RR[i] = r = 0.45;
   ZZ[i] = z = 0.65;
   while(r < rr || z > zz){
      r1 = r;
      z1 = z;
      z = z1 - d;
      while ((r=bf2(z)) > r1 + eps){
         d *= 0.7;
         z = z1 - d;
      }
      i++;
      RR[i] = r;
      ZZ[i] = z;
   }
   nlg = i;
   NLS = nlg + 900;
   for(i=nlg+1; i <= NLS; i++){
      RR[i] = 0.45 - (i-nlg)*0.0005;
      ZZ[i] = 0.65;
   }
   printf("NLS = %i\n",NLS);
}

void change_coord(x)
FLOAT *x;
{
   FLOAT r, s, d;
   INT i, j;
   
   r = sqrt(x[0]*x[0] + x[1]*x[1]);
   d = r*r + x[2]*x[2] + 1.;
   for (i = 0; i <= NLS; i++){
      s = (r-RR[i])*(r-RR[i]) + (x[2]-ZZ[i])*(x[2]-ZZ[i]);
      if (s < d){
         d = s;
         j = i;
      }
   }
   if (r > EPSA){
      x[0] = x[0]*RR[j]/r;
      x[1] = x[1]*RR[j]/r;
   }
   else if (RR[j] < EPSA)
      x[0] = x[1] = 0.;
   else
      eprintf("Error in change_coord.\n");
   x[2] = ZZ[j];
}

void shift_curve_point(pnode,r)   
NODE *pnode;
FLOAT r;
{
   FLOAT q;
   
   q = r/sqrt(MYVERTEX(pnode)->x[0]*MYVERTEX(pnode)->x[0] + 
              MYVERTEX(pnode)->x[1]*MYVERTEX(pnode)->x[1]);
   MYVERTEX(pnode)->x[0] *= q;
   MYVERTEX(pnode)->x[1] *= q;
}

void shift_boundary_vertices(theGrid)  /*  boundary vertices are moved onto the    */
GRID *theGrid;                         /*                                boundary  */
{
   NODE *pnode;
   
   for (pnode = FDBN(theGrid); pnode != FIRSTNODE(theGrid); pnode = pnode->succ)
      if (pnode->father == NULL && IS_BN(pnode))
         if (NTYPE(pnode) & GW && NTYPE(pnode) & GLG)
            shift_curve_point(pnode,1.);
         else if (NTYPE(pnode) & GLG && NTYPE(pnode) & GLS)
            shift_curve_point(pnode,0.45);
         else if (NTYPE(pnode) & GWW)
            shift_curve_point(pnode,1.);
         else
            change_coord(MYVERTEX(pnode)->x);
}

void nLG(r,n,np)  /*  r \in (0.45,1.>  */
FLOAT r, n[3], np[3];
{
   FLOAT d, f1, f2, g1, g2, s, z, z1, z2;
   INT i, j;
   
   if (r >= 1.)
      z = 0.6;
   else{
      d = r + 1.;
      for (i = NW+1; i <= NLG; i++){
         s = fabs(r - RR[i]);
         if (s < d){
            d = s;
            j = i;
         }
      }
      if (j == NLG && r > RR[j]){
         z1 = ZZ[j];
         z2 = 0.6;
      }
      else if (j == NW+1 || r > RR[j]){
         z1 = ZZ[j];
         z2 = ZZ[j+1];
      }
      else{
         z1 = ZZ[j-1];
         z2 = ZZ[j];
      }
      if (z2 >= z1){
         z1 = 0.65;
         z2 = 0.6;
      }
      while (r < bf2(z1)) z1 += 0.0001;
      while (r > bf2(z2)) z2 -= 0.0001;
      z = (z1 + z2)/2.;
      while (fabs( (d=bf2(z)) - r ) > 1.e-5){
         if (r > d)
            z1 = z;
         else
            z2 = z;
         z = (z1 + z2)/2.;
      }
   }
   z = z - 0.5999825;
   d = sqrt(0.02 - z*z);
   s = sqrt(0.02) + d;
   r = z/s/d;
   g1 = -0.07/z + z*(1. - 0.07/s)/d;
   g2 = 0.07/z/z + (1. - 0.07/s)*0.02/d/d/d - 0.07*r*r;
   f1 = 1./g1;
   f2 = -g2/g1/g1/g1;
   d = sqrt(1. + f1*f1);
   n[1] = -f1/d;
   n[2] = 1./d;
   np[1] = -f2/d/d/d;
   np[2] = f1*np[1];
}

void values_of_m(x,m,divm,gm)
FLOAT x[3], m[3], *divm, gm[3][3];
{
   FLOAT a = 79.414609748,  b = -160.220945, 
         c = -1.7335173717, d =   17.519410642,
         r, n[3], np[3];
   
   r = sqrt(x[0]*x[0] + x[1]*x[1]);
   if (r <= 0.45){
      m[0] = (a*r*r + b*r*r*r)*x[0];
      m[1] = (a*r*r + b*r*r*r)*x[1];
      m[2] = -0.5 + c*r*r + d*r*r*r;
      gm[0][0] = (2.*a + 3.*b*r)*x[0]*x[0] + a*r*r + b*r*r*r;
      gm[0][1] = (2.*a + 3.*b*r)*x[0]*x[1];
      gm[0][2] = 0.;
      gm[1][0] = (2.*a + 3.*b*r)*x[0]*x[1];
      gm[1][1] = (2.*a + 3.*b*r)*x[1]*x[1] + a*r*r + b*r*r*r;
      gm[1][2] = 0.;
      gm[2][0] = (2.*c + 3.*d*r)*x[0];
      gm[2][1] = (2.*c + 3.*d*r)*x[1];
      gm[2][2] = 0.;
   }
   else{
      nLG(r,n,np);
      m[0] = n[1]*x[0]/r;
      m[1] = n[1]*x[1]/r;
      m[2] = n[2];
      gm[0][0] = (np[1] - n[1]/r)*x[0]*x[0]/r/r + n[1]/r;
      gm[0][1] = (np[1] - n[1]/r)*x[0]*x[1]/r/r;
      gm[0][2] = 0.;
      gm[1][0] = (np[1] - n[1]/r)*x[0]*x[1]/r/r;
      gm[1][1] = (np[1] - n[1]/r)*x[1]*x[1]/r/r + n[1]/r;
      gm[1][2] = 0.;
      gm[2][0] = np[2]*x[0]/r;
      gm[2][1] = np[2]*x[1]/r;
      gm[2][2] = 0.;
   }
   *divm = gm[0][0] + gm[1][1] + gm[2][2];
}

#else

void boundary_points()
{}

void shift_boundary_vertices(theGrid)
GRID *theGrid;
{}

void values_of_m(x,m,divm,gm)
FLOAT x[3], m[3], *divm, gm[3][3];
{
   m[0] = 0.;
   m[1] = 0.;
   m[2] = 1.;
   *divm = 0.;
   gm[0][0] = gm[0][1] = gm[0][2] = 0.;
   gm[1][0] = gm[1][1] = gm[1][2] = 0.;
   gm[2][0] = gm[2][1] = gm[2][2] = 0.;
}

#endif

