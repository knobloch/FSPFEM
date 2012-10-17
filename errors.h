/******************************************************************************/
/*                                                                            */
/*                         evaluation of the results                          */
/*                                                                            */
/******************************************************************************/

FLOAT sd_delta();

#if DIM == 3

/* (\int_K p^2 \dx)/|K|; exact for p\in\P_3 */
FLOAT p32(p1,p2,p3,p4,p112,p113,p114,
                      p221,p223,p224,
                      p331,p332,p334,
                      p441,p442,p443,
                      p123,p124,p134,p234)
                    
FLOAT p1,p2,p3,p4,p112,p113,p114,
                  p221,p223,p224,
                  p331,p332,p334,
                  p441,p442,p443,
                  p123,p124,p134,p234;
{
   return(( 2.0*( p1*(p2+p3+p4)+p2*(p3+p4)+p3*p4 ) +
            8.0*( p1*p1+p2*p2+p3*p3+p4*p4 ) +
            3.0*( p1*(p223+p224+p332+p334+p442+p443)+
                  p2*(p113+p114+p331+p334+p441+p443)+
                  p3*(p112+p114+p221+p224+p441+p442)+
                  p4*(p112+p113+p221+p223+p331+p332) ) -
           12.0*( p1*(p112+p113+p114)+
                  p2*(p221+p223+p224)+
                  p3*(p331+p332+p334)+
                  p4*(p441+p442+p443) ) +
            6.0*( p1*(p221+p331+p441)+
                  p2*(p112+p332+p442)+
                  p3*(p113+p223+p443)+
                  p4*(p114+p224+p334) ) +
           36.0*( p1*p234+p2*p134+p3*p124+p4*p123 ) +
           18.0*( p123*(p1+p2+p3)+p124*(p1+p2+p4)+p134*(p1+p3+p4)+
                  p234*(p2+p3+p4) ) +
           54.0*( p112*p112+p113*p113+p114*p114+
                  p221*p221+p223*p223+p224*p224+
                  p331*p331+p332*p332+p334*p334+
                  p441*p441+p442*p442+p443*p443 ) -
           54.0*( p112*p221+p113*p331+p114*p441+p223*p332+p224*p442+p334*p443 )+
           54.0*( p112*(p113+p114)+p113*p114+p221*(p223+p224)+p223*p224+
                  p331*(p332+p334)+p332*p334+p441*(p442+p443)+p442*p443 ) -
           27.0*( p112*(p223+p224)+p113*(p332+p334)+p114*(p442+p443)+
                  p221*(p113+p114)+p223*(p331+p334)+p224*(p441+p443)+
                  p331*(p112+p114)+p332*(p221+p224)+p334*(p441+p442)+
                  p441*(p112+p113)+p442*(p221+p223)+p443*(p331+p332) ) -
           54.0*( p123*(p441+p442+p443)+p124*(p331+p332+p334)+
                  p134*(p221+p223+p224)+p234*(p112+p113+p114) ) +
          216.0*( p123*p123+p124*p124+p134*p134+p234*p234 ) + 
          216.0*( p123*(p124+p134+p234)+p124*(p134+p234)+p134*p234 ) )/2240.0 );
}

void points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                        x334,x441,x442,x443,x123,x124,x134,x234)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, 
     *x331, *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234;
{
   POINT2(x1,x2,x112,2.0,3.0);
   POINT2(x1,x3,x113,2.0,3.0);
   POINT2(x1,x4,x114,2.0,3.0);
   POINT2(x2,x1,x221,2.0,3.0);
   POINT2(x2,x3,x223,2.0,3.0);
   POINT2(x2,x4,x224,2.0,3.0);
   POINT2(x3,x1,x331,2.0,3.0);
   POINT2(x3,x2,x332,2.0,3.0);
   POINT2(x3,x4,x334,2.0,3.0);
   POINT2(x4,x1,x441,2.0,3.0);
   POINT2(x4,x2,x442,2.0,3.0);
   POINT2(x4,x3,x443,2.0,3.0);
   POINT3(x1,x2,x3,x123);
   POINT3(x1,x2,x4,x124);
   POINT3(x1,x3,x4,x134);
   POINT3(x2,x3,x4,x234);
}
  
void lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                   x331,x332,x334,x441,x442,x443,
                       p1,p2,p3,p4,p112,p113,p114,p221,p223,p224,
                                   p331,p332,p334,p441,p442,p443,
                       v1,v2,v3,v4)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224,
                          *x331, *x332, *x334, *x441, *x442, *x443,
      *p1, *p2, *p3, *p4, *p112, *p113, *p114, *p221, *p223, *p224,
                          *p331, *p332, *p334, *p441, *p442, *p443,
      v1, v2, v3, v4;                       
{
   *p1 = v1;
   *p2 = v2;
   *p3 = v3;
   *p4 = v4;
   *p112 = (2.0*v1 + v2)/3.0;
   *p113 = (2.0*v1 + v3)/3.0;
   *p114 = (2.0*v1 + v4)/3.0;
   *p221 = (2.0*v2 + v1)/3.0;
   *p223 = (2.0*v2 + v3)/3.0;
   *p224 = (2.0*v2 + v4)/3.0;
   *p331 = (2.0*v3 + v1)/3.0;
   *p332 = (2.0*v3 + v2)/3.0;
   *p334 = (2.0*v3 + v4)/3.0;
   *p441 = (2.0*v4 + v1)/3.0;
   *p442 = (2.0*v4 + v2)/3.0;
   *p443 = (2.0*v4 + v3)/3.0;
}

void L2_lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                      x331,x332,x334,x441,x442,x443,
                          p1,p2,p3,p4,p112,p113,p114,p221,p223,p224,
                                      p331,p332,p334,p441,p442,p443,
                          v1,v2,v3,v4,u0)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224,
                          *x331, *x332, *x334, *x441, *x442, *x443,
      *p1, *p2, *p3, *p4, *p112, *p113, *p114, *p221, *p223, *p224,
                          *p331, *p332, *p334, *p441, *p442, *p443,
      v1, v2, v3, v4, (*u0)();                       
{
   *p1 = v1 - (*u0)(x1);
   *p2 = v2 - (*u0)(x2);
   *p3 = v3 - (*u0)(x3);
   *p4 = v4 - (*u0)(x4);
   *p112 = (2.0*v1 + v2)/3.0 - (*u0)(x112);
   *p113 = (2.0*v1 + v3)/3.0 - (*u0)(x113);
   *p114 = (2.0*v1 + v4)/3.0 - (*u0)(x114);
   *p221 = (2.0*v2 + v1)/3.0 - (*u0)(x221);
   *p223 = (2.0*v2 + v3)/3.0 - (*u0)(x223);
   *p224 = (2.0*v2 + v4)/3.0 - (*u0)(x224);
   *p331 = (2.0*v3 + v1)/3.0 - (*u0)(x331);
   *p332 = (2.0*v3 + v2)/3.0 - (*u0)(x332);
   *p334 = (2.0*v3 + v4)/3.0 - (*u0)(x334);
   *p441 = (2.0*v4 + v1)/3.0 - (*u0)(x441);
   *p442 = (2.0*v4 + v2)/3.0 - (*u0)(x442);
   *p443 = (2.0*v4 + v3)/3.0 - (*u0)(x443);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT oldL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + face bubbles                 */
FLOAT (*u0)();           /* j = 0, 1, 2                                       */
INT u, j;
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                         x221[3],x223[3],x224[3],
                         x331[3],x332[3],x334[3],
                         x441[3],x442[3],x443[3],
                         x123[3],x124[3],x134[3],x234[3],
                         nn1[3],nn2[3],nn3[3],nn4[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234; 
   NODE *n1, *n2, *n3, *n4;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   L2_lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                    x331,x332,x334,x441,x442,x443,
              &p1,&p2,&p3,&p4,&p112,&p113,&p114,&p221,&p223,&p224,
                              &p331,&p332,&p334,&p441,&p442,&p443,
              ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),ND(n4,u,j),u0);
   p123 = (ND(n1,u,j)+ND(n2,u,j)+ND(n3,u,j))/3.0 + 
          FD(ltop_face(pelement->f[3],tGrid),u)*FMULT*nn4[j]/27.0 -(*u0)(x123);
   p124 = (ND(n1,u,j)+ND(n2,u,j)+ND(n4,u,j))/3.0 + 
          FD(ltop_face(pelement->f[2],tGrid),u)*FMULT*nn3[j]/27.0 -(*u0)(x124);
   p134 = (ND(n1,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0 + 
          FD(ltop_face(pelement->f[1],tGrid),u)*FMULT*nn2[j]/27.0 -(*u0)(x134);
   p234 = (ND(n2,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0 + 
          FD(ltop_face(pelement->f[0],tGrid),u)*FMULT*nn1[j]/27.0 -(*u0)(x234);
     
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}
             
#else

FLOAT oldL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: oldL2s not available.\n"); return(0.);  }

#endif

#if F_DATA & SCALAR_FACE_DATA

FLOAT oldbubbleL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of */
GRID *tGrid;                           /* the function saved in u             */
ELEMENT *pelement;                     /* u = face bubbles                    */
FLOAT (*u0)();                         /* j = 0, 1, 2                         */
INT u, j;
{
   FLOAT *x1, *x2, *x3, *x4, nn1[3], nn2[3], nn3[3], nn4[3], 
         p123, p124, p134, p234; 
                     
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   p123 = FD(ltop_face(pelement->f[3],tGrid),u)*FMULT*nn4[j]/27.0;
   p124 = FD(ltop_face(pelement->f[2],tGrid),u)*FMULT*nn3[j]/27.0;
   p134 = FD(ltop_face(pelement->f[1],tGrid),u)*FMULT*nn2[j]/27.0;
   p234 = FD(ltop_face(pelement->f[0],tGrid),u)*FMULT*nn1[j]/27.0;
   return( p32(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}

#else

FLOAT oldbubbleL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: oldbubbleL2s not available.\n"); return(0.);  }

#endif

FLOAT L2s_for_subtetr(x1,x2,x3,x4,v1,v2,v3,v4,nn4j,face_value,u0)
FLOAT *x1, *x2, *x3, *x4, v1, v2, v3, v4, nn4j, face_value, (*u0)();
{                                                     /* x4 is the barycentre */
   FLOAT x112[3], x113[3], x114[3], x221[3], x223[3], x224[3], x331[3], x332[3],
         x334[3], x441[3], x442[3], x443[3], x123[3], x124[3], x134[3], x234[3],
         p1, p2, p3, p4, p112, p113, p114, p221, p223, p224, p331, p332, p334,
                                     p441, p442, p443, p123, p124, p134, p234; 
                     
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   L2_lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                    x331,x332,x334,x441,x442,x443,
              &p1,&p2,&p3,&p4,&p112,&p113,&p114,&p221,&p223,&p224,
                              &p331,&p332,&p334,&p441,&p442,&p443,
              v1,v2,v3,v4,u0);
   p123 = (v1 + v2 + v3)/3.0 + face_value*FMULT*nn4j/27.0 - (*u0)(x123);
   p124 = (v1 + v2 + v4)/3.0 - (*u0)(x124);
   p134 = (v1 + v3 + v4)/3.0 - (*u0)(x134);
   p234 = (v2 + v3 + v4)/3.0 - (*u0)(x234);
     
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT newL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + face bubbles                 */
FLOAT (*u0)();           /* j = 0, 1, 2                                       */
INT u, j;
{
   FLOAT *x1, *x2, *x3, *x4, xc[3], nn1[3], nn2[3], nn3[3], nn4[3],
         v1, v2, v3, v4, vc; 
   NODE *n1, *n2, *n3, *n4;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   POINT4(x1,x2,x3,x4,xc);
   v1 = ND(n1,u,j);
   v2 = ND(n2,u,j);
   v3 = ND(n3,u,j);
   v4 = ND(n4,u,j);
   vc = (v1 + v2 + v3 + v4)/4.;
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   return((L2s_for_subtetr(x1,x2,x3,xc,v1,v2,v3,vc,nn4[j],
                                     FD(ltop_face(pelement->f[3],tGrid),u),u0)+
           L2s_for_subtetr(x1,x2,x4,xc,v1,v2,v4,vc,nn3[j],
                                     FD(ltop_face(pelement->f[2],tGrid),u),u0)+
           L2s_for_subtetr(x1,x3,x4,xc,v1,v3,v4,vc,nn2[j],
                                     FD(ltop_face(pelement->f[1],tGrid),u),u0)+
           L2s_for_subtetr(x2,x3,x4,xc,v2,v3,v4,vc,nn1[j],
            FD(ltop_face(pelement->f[0],tGrid),u),u0))*volume(x1,x2,x3,x4)/4.);
}

#else

FLOAT newL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: newL2s not available.\n"); return(0.);  }

#endif

#if F_DATA & SCALAR_FACE_DATA

FLOAT newbubbleL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of */
GRID *tGrid;                           /* the function saved in u             */
ELEMENT *pelement;                     /* u = face bubbles                    */
FLOAT (*u0)();                         /* j = 0, 1, 2                         */
INT u, j;
{
   FLOAT *x1, *x2, *x3, *x4, nn1[3], nn2[3], nn3[3], nn4[3], fv1, fv2, fv3, fv4;
                     
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   fv1 = FD(ltop_face(pelement->f[3],tGrid),u)*nn4[j];
   fv2 = FD(ltop_face(pelement->f[2],tGrid),u)*nn3[j];
   fv3 = FD(ltop_face(pelement->f[1],tGrid),u)*nn2[j];
   fv4 = FD(ltop_face(pelement->f[0],tGrid),u)*nn1[j];
   return( (fv1*fv1 + fv2*fv2 + fv3*fv3 + fv4*fv4)*
              volume(x1,x2,x3,x4)*FMULT*FMULT/30240.);
}

#else

FLOAT newbubbleL2s(tGrid,pelement,u0,u,j) 
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: newbubbleL2s not available.\n"); return(0.);  }

#endif

FLOAT vp1cL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: vp1cL2s not available.\n"); return(0.);  }

FLOAT miniL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: miniL2s not available.\n"); return(0.);  }

FLOAT minibubbleL2s(tGrid,pelement,u0,u,j) 
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: minibubbleL2s not available.\n"); return(0.);  }

FLOAT L2sq(tGrid,pelement,u0,u,j)  /* L2 norm of the j-th component of the   */
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: L2sq not available.\n"); return(0.);  }

FLOAT L2sq_iso(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: L2sq_iso not available.\n"); return(0.);  }

FLOAT H10q_iso(tGrid,pelement,u011,u012,u021,u022,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u011)(), (*u012)(), (*u021)(), (*u022)(); INT u;
{  eprintf("Error: H10q_iso not available.\n"); return(0.);  }

#if (N_DATA & VECTOR_NODE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT linL2norms(tGrid,pelement,u,j)    /* L2 norm of the j-th component of   */
GRID *tGrid;                            /* the function saved in u            */
ELEMENT *pelement;                      /* u = cont. pw. lin.                 */
INT u, j;                               /* j = 0, 1, 2                        */
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234; 
   NODE *n1, *n2, *n3, *n4;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                 x331,x332,x334,x441,x442,x443,
           &p1,&p2,&p3,&p4,&p112,&p113,&p114,&p221,&p223,&p224,
                           &p331,&p332,&p334,&p441,&p442,&p443,
                    ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),ND(n4,u,j));
   p123 = (ND(n1,u,j)+ND(n2,u,j)+ND(n3,u,j))/3.0;
   p124 = (ND(n1,u,j)+ND(n2,u,j)+ND(n4,u,j))/3.0;
   p134 = (ND(n1,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0;
   p234 = (ND(n2,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0;
     
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}
 
FLOAT linL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin.                                */
FLOAT (*u0)();           /* j = 0, 1, 2                                       */
INT u, j;
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234; 
   NODE *n1, *n2, *n3, *n4;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   L2_lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                    x331,x332,x334,x441,x442,x443,
              &p1,&p2,&p3,&p4,&p112,&p113,&p114,&p221,&p223,&p224,
                              &p331,&p332,&p334,&p441,&p442,&p443,
              ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),ND(n4,u,j),u0);
   p123 = (ND(n1,u,j)+ND(n2,u,j)+ND(n3,u,j))/3.0 - (*u0)(x123);
   p124 = (ND(n1,u,j)+ND(n2,u,j)+ND(n4,u,j))/3.0 - (*u0)(x124);
   p134 = (ND(n1,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0 - (*u0)(x134);
   p234 = (ND(n2,u,j)+ND(n3,u,j)+ND(n4,u,j))/3.0 - (*u0)(x234);
     
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}
 
FLOAT lin_div(tGrid,pelement,u) /*  || div u_lin ||_{0,pel}^2  */
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT c, detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3, *n4;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   detB = barycentric_coordinates(n4->myvertex->x,n1->myvertex->x,
                                  n2->myvertex->x,n3->myvertex->x,b);
   c = DOT(NDD(n1,u),b[1]) + DOT(NDD(n2,u),b[2]) +
       DOT(NDD(n3,u),b[3]) + DOT(NDD(n4,u),b[0]);
   return(c*c*detB);
}

#else

FLOAT linL2norms(tGrid,pelement,u,j)
GRID *tGrid; ELEMENT *pelement; INT u, j;
{  eprintf("Error: linL2norms not available.\n"); return(0.);  }
 
FLOAT linL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: linL2s not available.\n"); return(0.);  }
 
FLOAT lin_div(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u; 
{  eprintf("Error: lin_div not available.\n"); return(0.);  }

#endif

FLOAT pL2s(pelement,p0,c) /*  ||c-p0||_{0,pelement}^2  */
ELEMENT *pelement;
FLOAT (*p0)(), c;
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                         x221[3],x223[3],x224[3],
                         x331[3],x332[3],x334[3],
                         x441[3],x442[3],x443[3],
                         x123[3],x124[3],x134[3],x234[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234; 

   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   p1 = c - (*p0)(x1);
   p2 = c - (*p0)(x2);
   p3 = c - (*p0)(x3);
   p4 = c - (*p0)(x4);
   p112 = c - (*p0)(x112);
   p113 = c - (*p0)(x113);
   p114 = c - (*p0)(x114);
   p221 = c - (*p0)(x221);
   p223 = c - (*p0)(x223);
   p224 = c - (*p0)(x224);
   p331 = c - (*p0)(x331);
   p332 = c - (*p0)(x332);
   p334 = c - (*p0)(x334);
   p441 = c - (*p0)(x441);
   p442 = c - (*p0)(x442);
   p443 = c - (*p0)(x443);
   p123 = c - (*p0)(x123);
   p124 = c - (*p0)(x124);
   p134 = c - (*p0)(x134);
   p234 = c - (*p0)(x234);
    
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}

FLOAT pL2s_lin(pelement,p0,c,q1,q2,q3,q4) /* pressure on pelement consists of */
ELEMENT *pelement;               /* a constant part c and a pw. lin. part     */
FLOAT (*p0)(), c, q1, q2, q3, q4;/* q_i=value of lin. part in pelement->n[i-1]*/
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234;

   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   p1 = c + q1 - (*p0)(x1);
   p2 = c + q2 - (*p0)(x2);
   p3 = c + q3 - (*p0)(x3);
   p4 = c + q4 - (*p0)(x4);
   p112 = c + (2.*q1+q2)/3. - (*p0)(x112);
   p113 = c + (2.*q1+q3)/3. - (*p0)(x113);
   p114 = c + (2.*q1+q4)/3. - (*p0)(x114);
   p221 = c + (2.*q2+q1)/3. - (*p0)(x221);
   p223 = c + (2.*q2+q3)/3. - (*p0)(x223);
   p224 = c + (2.*q2+q4)/3. - (*p0)(x224);
   p331 = c + (2.*q3+q1)/3. - (*p0)(x331);
   p332 = c + (2.*q3+q2)/3. - (*p0)(x332);
   p334 = c + (2.*q3+q4)/3. - (*p0)(x334);
   p441 = c + (2.*q4+q1)/3. - (*p0)(x441);
   p442 = c + (2.*q4+q2)/3. - (*p0)(x442);
   p443 = c + (2.*q4+q3)/3. - (*p0)(x443);
   p123 = c + (q1+q2+q3)/3. - (*p0)(x123);
   p124 = c + (q1+q2+q4)/3. - (*p0)(x124);
   p134 = c + (q1+q3+q4)/3. - (*p0)(x134);
   p234 = c + (q2+q3+q4)/3. - (*p0)(x234);
    
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
} 

FLOAT pL2s_q(pelement,p0,c0)
ELEMENT *pelement; FLOAT (*p0)(), c0;
{  eprintf("Error: pL2s_q not available.\n"); return(0.);  }

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/* (L2-error of j-th derivative of i-th comp.)^2 / |K| on one element */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.;
   nn1i is the i-th component of the normal vector to f1 etc. */  
FLOAT oldH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
           x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
           nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i) /* i=0,1,2 */
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1,p2,p3,p4,p112,p113,p114,p221,p223,p224,p331,p332,p334,
         p441,p442,p443,p123,p124,p134,p234,c;
     
   c =ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j + ND(n4,u,i)*gl4j;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p4 = c - (*u0ij)(x4);
   p112 = Q29*FMULT*( FD(f3,u)*nn3i*gl4j + FD(f4,u)*nn4i*gl3j ) + c;
   p113 = Q29*FMULT*( FD(f2,u)*nn2i*gl4j + FD(f4,u)*nn4i*gl2j ) + c;
   p114 = Q29*FMULT*( FD(f2,u)*nn2i*gl3j + FD(f3,u)*nn3i*gl2j ) + c;
   p223 = Q29*FMULT*( FD(f1,u)*nn1i*gl4j + FD(f4,u)*nn4i*gl1j ) + c;
   p224 = Q29*FMULT*( FD(f1,u)*nn1i*gl3j + FD(f3,u)*nn3i*gl1j ) + c;
   p334 = Q29*FMULT*( FD(f1,u)*nn1i*gl2j + FD(f2,u)*nn2i*gl1j ) + c;
   p221 = p112;
   p331 = p113;
   p332 = p223;
   p441 = p114;
   p442 = p224;
   p443 = p334;
   p112 = p112 - (*u0ij)(x112);
   p113 = p113 - (*u0ij)(x113);
   p114 = p114 - (*u0ij)(x114);
   p221 = p221 - (*u0ij)(x221);
   p223 = p223 - (*u0ij)(x223);
   p224 = p224 - (*u0ij)(x224);
   p331 = p331 - (*u0ij)(x331);
   p332 = p332 - (*u0ij)(x332);
   p334 = p334 - (*u0ij)(x334);
   p441 = p441 - (*u0ij)(x441);
   p442 = p442 - (*u0ij)(x442);
   p443 = p443 - (*u0ij)(x443);
   p123 = ( FD(f1,u)*nn1i + FD(f2,u)*nn2i + FD(f3,u)*nn3i - 
            FD(f4,u)*nn4i )*FMULT*gl4j/9.0 + c  - (*u0ij)(x123);
   p124 = ( FD(f1,u)*nn1i + FD(f2,u)*nn2i + FD(f4,u)*nn4i - 
            FD(f3,u)*nn3i )*FMULT*gl3j/9.0 + c  - (*u0ij)(x124);
   p134 = ( FD(f1,u)*nn1i + FD(f3,u)*nn3i + FD(f4,u)*nn4i - 
            FD(f2,u)*nn2i )*FMULT*gl2j/9.0 + c  - (*u0ij)(x134);
   p234 = ( FD(f2,u)*nn2i + FD(f3,u)*nn3i + FD(f4,u)*nn4i - 
            FD(f1,u)*nn1i )*FMULT*gl1j/9.0 + c  - (*u0ij)(x234);
   
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}  

#else

FLOAT oldH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
           x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
           nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i) /* i=0,1,2 */
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldH10ij_3 not available.\n"); return(0.);  }

#endif

FLOAT oldH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i) 
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldH10ij_2 not available.\n"); return(0.);  }

FLOAT oldbubbleH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldbubbleH10ij_2 not available.\n"); return(0.);  }
 
#if F_DATA & SCALAR_FACE_DATA

FLOAT oldbubbleH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                  x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                  nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p112,p113,p114,p221,p223,p224,p331,p332,p334,
         p441,p442,p443,p123,p124,p134,p234;
     
   p112 = Q29*FMULT*( FD(f3,u)*nn3i*gl4j + FD(f4,u)*nn4i*gl3j );
   p113 = Q29*FMULT*( FD(f2,u)*nn2i*gl4j + FD(f4,u)*nn4i*gl2j );
   p114 = Q29*FMULT*( FD(f2,u)*nn2i*gl3j + FD(f3,u)*nn3i*gl2j );
   p223 = Q29*FMULT*( FD(f1,u)*nn1i*gl4j + FD(f4,u)*nn4i*gl1j );
   p224 = Q29*FMULT*( FD(f1,u)*nn1i*gl3j + FD(f3,u)*nn3i*gl1j );
   p334 = Q29*FMULT*( FD(f1,u)*nn1i*gl2j + FD(f2,u)*nn2i*gl1j );
   p221 = p112;
   p331 = p113;
   p332 = p223;
   p441 = p114;
   p442 = p224;
   p443 = p334;
   p123 = ( FD(f1,u)*nn1i + FD(f2,u)*nn2i + FD(f3,u)*nn3i - 
            FD(f4,u)*nn4i )*FMULT*gl4j/9.0;
   p124 = ( FD(f1,u)*nn1i + FD(f2,u)*nn2i + FD(f4,u)*nn4i - 
            FD(f3,u)*nn3i )*FMULT*gl3j/9.0;
   p134 = ( FD(f1,u)*nn1i + FD(f3,u)*nn3i + FD(f4,u)*nn4i - 
            FD(f2,u)*nn2i )*FMULT*gl2j/9.0;
   p234 = ( FD(f2,u)*nn2i + FD(f3,u)*nn3i + FD(f4,u)*nn4i - 
            FD(f1,u)*nn1i )*FMULT*gl1j/9.0;
   
   return( p32(0.,0.,0.,0.,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}

#else

FLOAT oldbubbleH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                  x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                  nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldbubbleH10ij_3 not available.\n"); return(0.);  }

#endif

FLOAT newH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
            x441,x442,x443,x123,x124,x134,x234,v1,v2,v3,v4,face_val,
            gl1j,gl2j,gl3j,gl4j,u0ij)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      v1, v2, v3, v4, face_val, gl1j, gl2j, gl3j, gl4j;  /* x4 is barycentre */
FLOAT (*u0ij)();
{
   FLOAT p1, p2, p3, p4, p112, p113, p114, p221, p223, p224, p331, p332, p334,
         p441, p442, p443, p123, p124, p134, p234, c;
     
   c = v1*gl1j + v2*gl2j + v3*gl3j + v4*gl4j;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p4 = c - (*u0ij)(x4);
   p114 = c - (*u0ij)(x114);
   p224 = c - (*u0ij)(x224);
   p334 = c - (*u0ij)(x334);
   p441 = c - (*u0ij)(x441);
   p442 = c - (*u0ij)(x442);
   p443 = c - (*u0ij)(x443);
   p112 = Q29*face_val*gl3j + c;
   p113 = Q29*face_val*gl2j + c;
   p223 = Q29*face_val*gl1j + c;
   p221 = p112;
   p331 = p113;
   p332 = p223;
   p112 = p112 - (*u0ij)(x112);
   p113 = p113 - (*u0ij)(x113);
   p221 = p221 - (*u0ij)(x221);
   p223 = p223 - (*u0ij)(x223);
   p331 = p331 - (*u0ij)(x331);
   p332 = p332 - (*u0ij)(x332);
   p123 = -face_val*gl4j/9.0 + c  - (*u0ij)(x123);
   p124 =  face_val*gl3j/9.0 + c  - (*u0ij)(x124);
   p134 =  face_val*gl2j/9.0 + c  - (*u0ij)(x134);
   p234 =  face_val*gl1j/9.0 + c  - (*u0ij)(x234);
   
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}  

FLOAT newbubbleH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                  x441,x442,x443,x123,x124,x134,x234,v1,v2,v3,v4,face_val,
                  gl1j,gl2j,gl3j,gl4j,u0ij)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      v1, v2, v3, v4, face_val, gl1j, gl2j, gl3j, gl4j;  /* x4 is barycentre */
FLOAT (*u0ij)();
{
   FLOAT p112, p113, p221, p223, p331, p332, p123, p124, p134, p234;
     
   p112 = Q29*face_val*gl3j;
   p113 = Q29*face_val*gl2j;
   p223 = Q29*face_val*gl1j;
   p221 = p112;
   p331 = p113;
   p332 = p223;
   p123 = -face_val*gl4j/9.0;
   p124 =  face_val*gl3j/9.0;
   p134 =  face_val*gl2j/9.0;
   p234 =  face_val*gl1j/9.0;
   
   return( p32(0.,0.,0.,0.,p112,p113,0.,
                           p221,p223,0.,
                           p331,p332,0.,
                           0.,0.,0.,p123,p124,p134,p234) );
}  

/* H10 norm on a subelement                                                   */
FLOAT H10_for_sub(x1,x2,x3,x4,v10,v20,v30,v40,v11,v21,v31,v41,v12,v22,v32,v42,
                                                           face_val,nn4,H10ij_3)
FLOAT *x1, *x2, *x3, *x4, v10, v20, v30, v40, v11, v21, v31, v41, 
      v12, v22, v32, v42, face_val, *nn4, (*H10ij_3)();
{
   FLOAT x112[3], x113[3], x114[3], x221[3], x223[3], x224[3], x331[3], x332[3],
         x334[3], x441[3], x442[3], x443[3], x123[3], x124[3], x134[3], x234[3],
         face_val0, face_val1, face_val2, detB, b[DIM2][DIM2];
   
   detB = barycentric_coordinates(x4,x1,x2,x3,b);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   face_val0 = face_val*FMULT*nn4[0];
   face_val1 = face_val*FMULT*nn4[1];
   face_val2 = face_val*FMULT*nn4[2];
   
   return((H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v10,v20,v30,v40,face_val0,
          b[1][0],b[2][0],b[3][0],b[0][0],u011) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v10,v20,v30,v40,face_val0,
          b[1][1],b[2][1],b[3][1],b[0][1],u012) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v10,v20,v30,v40,face_val0,
          b[1][2],b[2][2],b[3][2],b[0][2],u013) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v11,v21,v31,v41,face_val1,
          b[1][0],b[2][0],b[3][0],b[0][0],u021) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v11,v21,v31,v41,face_val1,
          b[1][1],b[2][1],b[3][1],b[0][1],u022) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v11,v21,v31,v41,face_val1,
          b[1][2],b[2][2],b[3][2],b[0][2],u023) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v12,v22,v32,v42,face_val2,
          b[1][0],b[2][0],b[3][0],b[0][0],u031) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v12,v22,v32,v42,face_val2,
          b[1][1],b[2][1],b[3][1],b[0][1],u032) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,v12,v22,v32,v42,face_val2,
          b[1][2],b[2][2],b[3][2],b[0][2],u033))*detB);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT H10s(tGrid,pelement,u)
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1, *x2, *x3, *x4, xc[DIM], nn1[DIM], nn2[DIM], nn3[DIM], nn4[DIM],
         v10, v20, v30, v40, vc0, 
         v11, v21, v31, v41, vc1, 
         v12, v22, v32, v42, vc2;
   NODE *n1, *n2, *n3, *n4;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   POINT4(x1,x2,x3,x4,xc);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   v10 = ND(n1,u,0);
   v20 = ND(n2,u,0);
   v30 = ND(n3,u,0);
   v40 = ND(n4,u,0);
   v11 = ND(n1,u,1);
   v21 = ND(n2,u,1);
   v31 = ND(n3,u,1);
   v41 = ND(n4,u,1);
   v12 = ND(n1,u,2);
   v22 = ND(n2,u,2);
   v32 = ND(n3,u,2);
   v42 = ND(n4,u,2);
   vc0 = (v10 + v20 + v30 + v40)/4.;
   vc1 = (v11 + v21 + v31 + v41)/4.;
   vc2 = (v12 + v22 + v32 + v42)/4.;
   
   return(H10_for_sub(x1,x2,x3,xc,v10,v20,v30,vc0,v11,v21,v31,vc1,v12,v22,v32,
                     vc2,FD(ltop_face(pelement->f[3],tGrid),u),nn4,newH10ij_3) +
          H10_for_sub(x1,x2,x4,xc,v10,v20,v40,vc0,v11,v21,v41,vc1,v12,v22,v42,
                     vc2,FD(ltop_face(pelement->f[2],tGrid),u),nn3,newH10ij_3) +
          H10_for_sub(x1,x3,x4,xc,v10,v30,v40,vc0,v11,v31,v41,vc1,v12,v32,v42,
                     vc2,FD(ltop_face(pelement->f[1],tGrid),u),nn2,newH10ij_3) +
          H10_for_sub(x2,x3,x4,xc,v20,v30,v40,vc0,v21,v31,v41,vc1,v22,v32,v42,
                     vc2,FD(ltop_face(pelement->f[0],tGrid),u),nn1,newH10ij_3));
}   

FLOAT bubbleH10(tGrid,pelement,u)
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1, *x2, *x3, *x4, xc[3], nn1[3], nn2[3], nn3[3], nn4[3];
   NODE *n1, *n2, *n3, *n4;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   POINT4(x1,x2,x3,x4,xc);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   return(H10_for_sub(x1,x2,x3,xc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[3],tGrid),u),nn4,newbubbleH10ij_3) +
          H10_for_sub(x1,x2,x4,xc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[2],tGrid),u),nn3,newbubbleH10ij_3) +
          H10_for_sub(x1,x3,x4,xc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[1],tGrid),u),nn2,newbubbleH10ij_3) +
          H10_for_sub(x2,x3,x4,xc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[0],tGrid),u),nn1,newbubbleH10ij_3));
}   

#else

FLOAT H10s(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{  eprintf("Error: H10s not available.\n"); return(0.);  }


FLOAT bubbleH10(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{  eprintf("Error: bubbleH10 not available.\n"); return(0.);  }

#endif 

FLOAT miniH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                  x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                  nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: miniH10ij_3 not available.\n"); return(0.);  }

FLOAT minibubbleH10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                  x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                  nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: minibubbleH10ij_3 not available.\n"); return(0.);  }

#if N_DATA & VECTOR_NODE_DATA

FLOAT linH10ijnorm_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                     x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                     nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT c;
     
   c =ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j + ND(n4,u,i)*gl4j;
   return( c*c );
}  

FLOAT H10ijlin_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                 x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                 nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1,p2,p3,p4,p112,p113,p114,p221,p223,p224,p331,p332,p334,
         p441,p442,p443,p123,p124,p134,p234,c;
     
   c =ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j + ND(n4,u,i)*gl4j;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p4 = c - (*u0ij)(x4);
   p112 = c - (*u0ij)(x112);
   p113 = c - (*u0ij)(x113);
   p114 = c - (*u0ij)(x114);
   p221 = c - (*u0ij)(x221);
   p223 = c - (*u0ij)(x223);
   p224 = c - (*u0ij)(x224);
   p331 = c - (*u0ij)(x331);
   p332 = c - (*u0ij)(x332);
   p334 = c - (*u0ij)(x334);
   p441 = c - (*u0ij)(x441);
   p442 = c - (*u0ij)(x442);
   p443 = c - (*u0ij)(x443);
   p123 = c  - (*u0ij)(x123);
   p124 = c  - (*u0ij)(x124);
   p134 = c  - (*u0ij)(x134);
   p234 = c  - (*u0ij)(x234);
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}  

#else

FLOAT linH10ijnorm_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                     x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                     nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: linH10ijnorm_3 not available.\n"); return(0.);  }

FLOAT H10ijlin_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                 x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
                 nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      nn1i,nn2i,nn3i,nn4i,gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4; FACE *f1, *f2, *f3, *f4; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: H10ijlin_3 not available.\n"); return(0.);  }

#endif

FLOAT linH10ijnorm_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: linH10ijnorm_2 not available.\n"); return(0.);  }

FLOAT H10ijlin_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: H10ijlin_2 not available.\n"); return(0.);  }

/* (H10-error on one element)^2 */
FLOAT H10_2(tGrid,pelement,u,H10ij_2)
GRID *tGrid; ELEMENT *pelement; INT u; FLOAT (*H10ij_2)();
{  eprintf("Error: H10_2 not available.\n"); return(0.);  }

FLOAT H10_3(tGrid,pelement,u,H10ij_3)
GRID *tGrid;
ELEMENT *pelement;
INT u;
FLOAT (*H10ij_3)();
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
                     nn1[3],nn2[3],nn3[3],nn4[3],
                     detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3, *n4;
   FACE *f1, *f2, *f3, *f4;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,f4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   NORMAL_VECTORS(nn1,nn2,nn3,nn4,pelement);
   detB = barycentric_coordinates(x4,x1,x2,x3,b);
   
   return((H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[0],nn2[0],nn3[0],nn4[0],b[1][0],b[2][0],b[3][0],b[0][0],
          u011,u,0) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[0],nn2[0],nn3[0],nn4[0],b[1][1],b[2][1],b[3][1],b[0][1],
          u012,u,0) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[0],nn2[0],nn3[0],nn4[0],b[1][2],b[2][2],b[3][2],b[0][2],
          u013,u,0) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[1],nn2[1],nn3[1],nn4[1],b[1][0],b[2][0],b[3][0],b[0][0],
          u021,u,1) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[1],nn2[1],nn3[1],nn4[1],b[1][1],b[2][1],b[3][1],b[0][1],
          u022,u,1) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[1],nn2[1],nn3[1],nn4[1],b[1][2],b[2][2],b[3][2],b[0][2],
          u023,u,1) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[2],nn2[2],nn3[2],nn4[2],b[1][0],b[2][0],b[3][0],b[0][0],
          u031,u,2) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[2],nn2[2],nn3[2],nn4[2],b[1][1],b[2][1],b[3][1],b[0][1],
          u032,u,2) +
   H10ij_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
          x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,f1,f2,f3,f4,
          nn1[2],nn2[2],nn3[2],nn4[2],b[1][2],b[2][2],b[3][2],b[0][2],
          u033,u,2))*detB);
}
 
/* (\int_K g \dx)/|K|; exact for g\in\P_3 */
FLOAT integr5(pelement,g)
ELEMENT *pelement;
FLOAT (*g)();
{
   FLOAT *x1, *x2, *x3, *x4, x123[3], x124[3], x134[3], x234[3];
                    
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   POINT3(x1,x2,x3,x123);
   POINT3(x1,x2,x4,x124);
   POINT3(x1,x3,x4,x134);
   POINT3(x2,x3,x4,x234);
   return(( (*g)(x1) + (*g)(x2) + (*g)(x3) + (*g)(x4) +
            9.0*((*g)(x123) + (*g)(x124) + (*g)(x134) + (*g)(x234)) )/40.0);
}

FLOAT integr_q5(pelement,g)
ELEMENT *pelement; FLOAT (*g)();
{  eprintf("Error: integr_q5 not available.\n"); return(0.);  }
 
#if N_DATA & SCALAR_NODE_DATA

FLOAT sL2s(tGrid,pelement,u0,u) /* L2 norm of the difference between u0       */
GRID *tGrid;                    /* and a scalar pw. lin. function saved in u  */
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
         p1,p2,p3,p4,p112,p113,p114,
                     p221,p223,p224,
                     p331,p332,p334,
                     p441,p442,p443,
                     p123,p124,p134,p234; 
   NODE *n1, *n2, *n3, *n4;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   L2_lin_node_values_3(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,
                                    x331,x332,x334,x441,x442,x443,
              &p1,&p2,&p3,&p4,&p112,&p113,&p114,&p221,&p223,&p224,
                              &p331,&p332,&p334,&p441,&p442,&p443,
                      NDS(n1,u),NDS(n2,u),NDS(n3,u),NDS(n4,u),u0);
   p123 = (NDS(n1,u) + NDS(n2,u) + NDS(n3,u))/3.0 - (*u0)(x123);
   p124 = (NDS(n1,u) + NDS(n2,u) + NDS(n4,u))/3.0 - (*u0)(x124);
   p134 = (NDS(n1,u) + NDS(n3,u) + NDS(n4,u))/3.0 - (*u0)(x134);
   p234 = (NDS(n2,u) + NDS(n3,u) + NDS(n4,u))/3.0 - (*u0)(x234);
     
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234)*volume(x1,x2,x3,x4) );
}
             
/* (L2-error of j-th derivative)^2 / |K| on one element                       */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.          */
FLOAT sH10j(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
            x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,
            gl1j,gl2j,gl3j,gl4j,u0j,u)
FLOAT *x1, *x2, *x3, *x4, *x112, *x113, *x114, *x221, *x223, *x224, *x331, 
      *x332, *x334, *x441, *x442, *x443, *x123, *x124, *x134, *x234, 
      gl1j,gl2j,gl3j,gl4j;
NODE *n1, *n2, *n3, *n4;
FLOAT (*u0j)();
INT u;
{
   FLOAT p1,p2,p3,p4,p112,p113,p114,p221,p223,p224,p331,p332,p334,
         p441,p442,p443,p123,p124,p134,p234,c;
     
   c = NDS(n1,u)*gl1j + NDS(n2,u)*gl2j + NDS(n3,u)*gl3j + NDS(n4,u)*gl4j;
   p1 = c - (*u0j)(x1);
   p2 = c - (*u0j)(x2);
   p3 = c - (*u0j)(x3);
   p4 = c - (*u0j)(x4);
   p112 = c - (*u0j)(x112);
   p113 = c - (*u0j)(x113);
   p114 = c - (*u0j)(x114);
   p221 = c - (*u0j)(x221);
   p223 = c - (*u0j)(x223);
   p224 = c - (*u0j)(x224);
   p331 = c - (*u0j)(x331);
   p332 = c - (*u0j)(x332);
   p334 = c - (*u0j)(x334);
   p441 = c - (*u0j)(x441);
   p442 = c - (*u0j)(x442);
   p443 = c - (*u0j)(x443);
   p123 = c - (*u0j)(x123);
   p124 = c - (*u0j)(x124);
   p134 = c - (*u0j)(x134);
   p234 = c - (*u0j)(x234);
   return( p32(p1,p2,p3,p4,p112,p113,p114,
                           p221,p223,p224,
                           p331,p332,p334,
                           p441,p442,p443,
                           p123,p124,p134,p234) );
}  
 
/* (H10-error on one element)^2 */
FLOAT sH10(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u; 
{
   FLOAT *x1,*x2,*x3,*x4,x112[3],x113[3],x114[3],
                     x221[3],x223[3],x224[3],
                     x331[3],x332[3],x334[3],
                     x441[3],x442[3],x443[3],
                     x123[3],x124[3],x134[3],x234[3],
                     detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3, *n4;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,n4,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,x4,pelement);
   points(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,
                      x334,x441,x442,x443,x123,x124,x134,x234);
   detB = barycentric_coordinates(x4,x1,x2,x3,b);
   
   return((sH10j(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                 x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,
                 b[1][0],b[2][0],b[3][0],b[0][0],u01,u) +
           sH10j(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                 x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,
                 b[1][1],b[2][1],b[3][1],b[0][1],u02,u) +
           sH10j(x1,x2,x3,x4,x112,x113,x114,x221,x223,x224,x331,x332,x334,
                 x441,x442,x443,x123,x124,x134,x234,n1,n2,n3,n4,
                 b[1][2],b[2][2],b[3][2],b[0][2],u03,u))*detB);   
}
 
#else

FLOAT sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2s not available.\n"); return(0.);  }

FLOAT sH10(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u; 
{  eprintf("Error: sH10 not available.\n"); return(0.);  }

#endif  /* N_DATA & SCALAR_NODE_DATA */

FLOAT p1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{ eprintf("Error: p1c_L2conv not available.\n"); return(0.); }

FLOAT p1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: p1c_Linf not available.\n"); return(0.);  }

FLOAT sL2_q1(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2_q1 not available.\n"); return(0.);  }

FLOAT sH10_q1(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: sH10_q1 not available.\n"); return(0.);  }

FLOAT q1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: q1c_L2conv not available.\n"); return(0.);  }

FLOAT q1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: q1c_Linf not available.\n"); return(0.);  }

FLOAT sL2sq(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2sq not available.\n"); return(0.);  }

FLOAT sH10q(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: sH10q not available.\n"); return(0.);  }

FLOAT p2c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: p2c_L2conv not available.\n"); return(0.);  }

FLOAT p2c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: p2c_Linf not available.\n"); return(0.);  }

FLOAT nc_sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: nc_sL2s not available.\n"); return(0.);  }

FLOAT nc_L2conv(tGrid,pelement,bb0,bb1,u01,u02,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*u01)(), (*u02)(); INT u;
{  eprintf("Error: nc_L2conv not available.\n"); return(0.);  }

FLOAT nc_sH10j(x1,x2,x3,x12,x13,x23,f1,f2,f3,gl1j,gl2j,gl3j,u0j,u)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
FACE *f1, *f2, *f3; FLOAT (*u0j)(); INT u;
{  eprintf("Error: nc_sH10j not available.\n"); return(0.);  }

FLOAT nc_L2s(tGrid,pelement,u0,u,i)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, i;
{ eprintf("Error: nc_L2s not available.\n"); return(0.);  }

FLOAT nc_H10(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{ eprintf("Error: nc_H10 not available.\n"); return(0.); }

FLOAT nc_lin_div(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{ eprintf("Error: nc_lin_div not available.\n"); return(0.); }

void nc_L_inf_error(tGrid,u0,u,err)
GRID *tGrid; FLOAT (*u0)(), *err; INT u;
{  eprintf("Error: nc_L_inf_error not available.\n");  }

FLOAT P1disc_sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{ eprintf("Error: P1disc_sL2s not available.\n"); return(0.);  }

/******************************************************************************/

#else  /* if DIM == 2 */

/* (\int_K p^2 \dx)/|K|; exact for p\in\P_2 */
FLOAT p22(p1,p2,p3,p12,p13,p23)
FLOAT p1,p2,p3,p12,p13,p23; 
{
   return(( 6.*(p1*p1+p2*p2+p3*p3) +
           32.*(p12*p12+p13*p13+p23*p23+
                p12*p13+p12*p23+p13*p23) -
            8.*(p1*p23+p2*p13+p3*p12) -
               (p1*(p2+p3)+p2*(p1+p3)+p3*(p1+p2)) )/180.);
}

/* (\int_K p^2 \dx)/|K|; exact for p\in\P_3 */
FLOAT p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)
FLOAT p1,p2,p3,p112,p113,p221,p223,p331,p332,p123; 
{
   return((76.*( p1*p1+p2*p2+p3*p3 ) +
           22.*( p1*p2+p1*p3+p2*p3 ) +
          540.*( p112*p112+p113*p113+p221*p221+p223*p223+p331*p331+p332*p332+
                 p112*p113+p221*p223+p331*p332 ) -
          378.*( p112*p221+p113*p331+p223*p332 ) -
          135.*( p112*p331+p113*p221+p221*p332+p223*p112+p331*p223+p332*p113+
                 p112*p223+p113*p332+p221*p113+p223*p331+p331*p112+p332*p221 ) -
          108.*( p112*p332+p113*p223+p221*p331 ) +
         1944.*p123*p123 +
           36.*( p1*(p112+p113) + p2*(p221+p223) + p3*(p331+p332) ) +
           54.*( p1*(p223+p332) + p2*(p113+p331) + p3*(p112+p221) ) +
          324.*p123*(p112+p113+p221+p223+p331+p332) +
           72.*p123*(p1+p2+p3))/6720.);
}

/* (\int_K p^2 \dx)/|K|; exact for p\in\P_4 */
FLOAT p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332)
FLOAT p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332;
{
   return((145.*( p1*p1+p2*p2+p3*p3 ) +
          5376.*( p11*p11+p22*p22+p33*p33 ) + 
          1584.*( p12*p12+p13*p13+p23*p23 ) + 
           384.*( p11*(p12+p13)+p22*(p12+p23)+p33*(p13+p23) ) + 
           512.*( p112*p223+p221*p332+p223*p331+p332*p113+p331*p112+p113*p221 )+
           256.*( p11*(p221+p331)+p22*(p112+p332)+p33*(p113+p223)+
                  p113*p223+p112*p332+p221*p331 ) + 
            64.*( p12*(p331+p332)+p13*(p221+p223)+p23*(p112+p113) ) + 
            48.*( p12*p13+p12*p23+p13*p23 ) + 
           160.*( p1*(p112+p113+p11)+p2*(p221+p223+p22)+p3*(p331+p332+p33) ) + 
           768.*( p112*p221+p113*p331+p223*p332 ) + 
          1280.*( p11*(p112+p113)+p22*(p223+p221)+p33*(p331+p332)+
                  p112*p112+p113*p113+p221*p221+p223*p223+p331*p331+
                  p332*p332+p221*p223+p112*p113+p331*p332 ) - 
         (1536.*( p11*p22+p11*p33+p22*p33 ) + 
           768.*( p11*p23+p22*p13+p33*p12 ) + 
            27.*( p1*p2+p1*p3+p2*p3 ) + 
          1280.*( p12*(p112+p221)+p13*(p113+p331)+p23*(p223+p332) ) + 
            80.*( p1*(p12+p13)+p2*(p12+p23)+p3*(p13+p23) ) + 
           160.*( p1*(p22+p33)+p2*(p11+p33)+p3*(p11+p22) ) + 
            12.*( p1*p23+p2*p13+p3*p12 ) + 
           112.*( p1*(p223+p332)+p2*(p113+p331)+p3*(p112+p221) ) + 
           256.*( p11*(p223+p332)+p22*(p113+p331)+p33*(p112+p221) ) + 
           960.*( p12*(p113+p223)+p13*(p112+p332)+p23*(p221+p331) )) )/28350.);
}

void points(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123;
{
   POINT2(x1,x2,x112,2.0,3.0);
   POINT2(x1,x3,x113,2.0,3.0);
   POINT2(x2,x1,x221,2.0,3.0);
   POINT2(x2,x3,x223,2.0,3.0);
   POINT2(x3,x1,x331,2.0,3.0);
   POINT2(x3,x2,x332,2.0,3.0);
   POINT3(x1,x2,x3,x123);
}

/* The element (x1,x2,x3) is the image of the reference element under a 
bilinear mapping. Thus only the edge (x2,x3) is curved. */
void mapped_points(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,a,c,alpha)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, 
      a[][DIM], *c, *alpha;
{
   POINT2(x1,x2,x112,2.0,3.0);
   POINT2(x1,x3,x113,2.0,3.0);
   POINT2(x2,x1,x221,2.0,3.0);
   POINT2(x3,x1,x331,2.0,3.0);
   V_BILIN_VALUE(x223,a,c,alpha,2./3.,1./3.)
   V_BILIN_VALUE(x332,a,c,alpha,1./3.,2./3.)
   V_BILIN_VALUE(x123,a,c,alpha,1./3.,1./3.)
}

void points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332)
FLOAT *x1, *x2, *x3, *x11, *x22, *x33, *x12, *x13, *x23, 
      *x112, *x113, *x221, *x223, *x331, *x332;
{
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   MIDPOINTS(x12,x23,x13,x22,x11,x33)
   POINT2(x1,x2,x112,3.0,4.0);
   POINT2(x1,x3,x113,3.0,4.0);
   POINT2(x2,x1,x221,3.0,4.0);
   POINT2(x2,x3,x223,3.0,4.0);
   POINT2(x3,x1,x331,3.0,4.0);
   POINT2(x3,x2,x332,3.0,4.0);
}

void test_p42()
{
   FLOAT x1[DIM], x2[DIM], x3[DIM], x11[DIM], x22[DIM], x33[DIM],
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], x221[DIM], 
         x223[DIM], x331[DIM], x332[DIM];

   x1[0] = 0.;
   x1[1] = 0.;
   x2[0] = 1.;
   x2[1] = 0.;
   x3[0] = 0.;
   x3[1] = 1.;
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   printf("integral = %e\n",
     p42(t0(x1),t0(x2),t0(x3),t0(x11),t0(x22),t0(x33),t0(x12),t0(x13),t0(x23),
         t0(x112),t0(x113),t0(x221),t0(x223),t0(x331),t0(x332)));
}

FLOAT abs_max_of_6(p1,p2,p3,p12,p13,p23)
FLOAT p1, p2, p3, p12, p13, p23;
{
   FLOAT max = 0.;

   max = MAX(max,fabs(p1));
   max = MAX(max,fabs(p2));
   max = MAX(max,fabs(p3));
   max = MAX(max,fabs(p12));
   max = MAX(max,fabs(p13));
   max = MAX(max,fabs(p23));
   return(max);
}

FLOAT abs_max_of_10(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)
FLOAT p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;
{
   FLOAT max = 0.;

   max = MAX(max,fabs(p1));
   max = MAX(max,fabs(p2));
   max = MAX(max,fabs(p3));
   max = MAX(max,fabs(p112));
   max = MAX(max,fabs(p113));
   max = MAX(max,fabs(p221));
   max = MAX(max,fabs(p223));
   max = MAX(max,fabs(p331));
   max = MAX(max,fabs(p332));
   max = MAX(max,fabs(p123));
   return(max);
}

FLOAT abs_max_of_15(p1,p2,p3,p11,p22,p33,p12,p13,p23,
                                                  p112,p113,p221,p223,p331,p332)
FLOAT p1, p2, p3, p11, p22, p33, p12, p13, p23, 
                                             p112, p113, p221, p223, p331, p332;
{
   FLOAT max = 0.;

   max = MAX(max,fabs(p1));
   max = MAX(max,fabs(p2));
   max = MAX(max,fabs(p3));
   max = MAX(max,fabs(p11));
   max = MAX(max,fabs(p22));
   max = MAX(max,fabs(p33));
   max = MAX(max,fabs(p12));
   max = MAX(max,fabs(p13));
   max = MAX(max,fabs(p23));
   max = MAX(max,fabs(p112));
   max = MAX(max,fabs(p113));
   max = MAX(max,fabs(p221));
   max = MAX(max,fabs(p223));
   max = MAX(max,fabs(p331));
   max = MAX(max,fabs(p332));
   return(max);
}

void lin_node_values_2(x1,x2,x3,x12,x13,x23,
		       p1,p2,p3,p12,p13,p23,v1,v2,v3)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23,
      *p1, *p2, *p3, *p12, *p13, *p23, v1, v2, v3;                       
{
   *p1 = v1;
   *p2 = v2;
   *p3 = v3;
   *p12 = (v1 + v2)/2.0;
   *p13 = (v1 + v3)/2.0;
   *p23 = (v2 + v3)/2.0;
}

void lin_node_values_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,
                       p1,p2,p3,p112,p113,p221,p223,p331,p332,v1,v2,v3)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, 
      *p1, *p2, *p3, *p112, *p113, *p221, *p223, *p331, *p332, v1, v2, v3;
{
   *p1 = v1;
   *p2 = v2;
   *p3 = v3;
   *p112 = (2.0*v1 + v2)/3.0;
   *p113 = (2.0*v1 + v3)/3.0;
   *p221 = (2.0*v2 + v1)/3.0;
   *p223 = (2.0*v2 + v3)/3.0;
   *p331 = (2.0*v3 + v1)/3.0;
   *p332 = (2.0*v3 + v2)/3.0;
}

void L2_lin_node_values_2(x1,x2,x3,x12,x13,x23,
			  p1,p2,p3,p12,p13,p23,v1,v2,v3,u0)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23,
      *p1, *p2, *p3, *p12, *p13, *p23, v1, v2, v3, (*u0)(); 
{
   *p1 = v1 - (*u0)(x1);
   *p2 = v2 - (*u0)(x2);
   *p3 = v3 - (*u0)(x3);
   *p12 = (v1 + v2)/2.0 - (*u0)(x12);
   *p13 = (v1 + v3)/2.0 - (*u0)(x13);
   *p23 = (v2 + v3)/2.0 - (*u0)(x23);
}

void L2_lin_face_values(x1,x2,x3,x12,x13,x23,
			p1,p2,p3,p12,p13,p23,v12,v13,v23,u0)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23,
      *p1, *p2, *p3, *p12, *p13, *p23, v12, v13, v23, (*u0)(); 
{
   *p1 = v12 + v13 - v23 - (*u0)(x1);
   *p2 = v12 + v23 - v13 - (*u0)(x2);
   *p3 = v13 + v23 - v12 - (*u0)(x3);
   *p12 = v12 - (*u0)(x12);
   *p13 = v13 - (*u0)(x13);
   *p23 = v23 - (*u0)(x23);
}

void L2_lin_node_values_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,
                          p1,p2,p3,p112,p113,p221,p223,p331,p332,v1,v2,v3,u0)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, 
      *p1, *p2, *p3, *p112, *p113, *p221, *p223, *p331, *p332, v1, v2, v3, 
      (*u0)();                       
{
   *p1 = v1 - (*u0)(x1);
   *p2 = v2 - (*u0)(x2);
   *p3 = v3 - (*u0)(x3);
   *p112 = (2.0*v1 + v2)/3.0 - (*u0)(x112);
   *p113 = (2.0*v1 + v3)/3.0 - (*u0)(x113);
   *p221 = (2.0*v2 + v1)/3.0 - (*u0)(x221);
   *p223 = (2.0*v2 + v3)/3.0 - (*u0)(x223);
   *p331 = (2.0*v3 + v1)/3.0 - (*u0)(x331);
   *p332 = (2.0*v3 + v2)/3.0 - (*u0)(x332);
}

void L2_lin_node_values_4(
          x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
          p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332,
          v1,v2,v3,u0)
FLOAT *x1, *x2, *x3, *x11, *x22, *x33, *x12, *x13, *x23, 
      *x112, *x113, *x221, *x223, *x331, *x332,
      *p1, *p2, *p3, *p11, *p22, *p33, *p12, *p13, *p23, 
      *p112, *p113, *p221, *p223, *p331, *p332, v1, v2, v3, (*u0)();                       
{
   FLOAT s;

   s = 0.25*(v1 + v2 + v3);
   *p1 = v1 - (*u0)(x1);
   *p2 = v2 - (*u0)(x2);
   *p3 = v3 - (*u0)(x3);
   *p12 = (v1 + v2)/2.0 - (*u0)(x12);
   *p13 = (v1 + v3)/2.0 - (*u0)(x13);
   *p23 = (v2 + v3)/2.0 - (*u0)(x23);
   *p11 = s + 0.25*v1 - (*u0)(x11);
   *p22 = s + 0.25*v2 - (*u0)(x22);
   *p33 = s + 0.25*v3 - (*u0)(x33);
   *p112 = (3.*v1 + v2)/4. - (*u0)(x112);
   *p113 = (3.*v1 + v3)/4. - (*u0)(x113);
   *p221 = (3.*v2 + v1)/4. - (*u0)(x221);
   *p223 = (3.*v2 + v3)/4. - (*u0)(x223);
   *p331 = (3.*v3 + v1)/4. - (*u0)(x331);
   *p332 = (3.*v3 + v2)/4. - (*u0)(x332);
}

void L2_mult_lin_node_values_4(
          x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
          p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332,
          v1,v2,v3,u0,b0)
FLOAT *x1, *x2, *x3, *x11, *x22, *x33, *x12, *x13, *x23, 
      *x112, *x113, *x221, *x223, *x331, *x332,
      *p1, *p2, *p3, *p11, *p22, *p33, *p12, *p13, *p23, 
      *p112, *p113, *p221, *p223, *p331, *p332, v1, v2, v3, (*u0)(), (*b0)();                       
{
   FLOAT s;

   s = 0.25*(v1 + v2 + v3);
   *p1 = (v1 - (*u0)(x1))*b0(x1);
   *p2 = (v2 - (*u0)(x2))*b0(x2);
   *p3 = (v3 - (*u0)(x3))*b0(x3);
   *p12 = ((v1 + v2)/2.0 - (*u0)(x12))*b0(x12);
   *p13 = ((v1 + v3)/2.0 - (*u0)(x13))*b0(x13);
   *p23 = ((v2 + v3)/2.0 - (*u0)(x23))*b0(x23);
   *p11 = (s + 0.25*v1 - (*u0)(x11))*b0(x11);
   *p22 = (s + 0.25*v2 - (*u0)(x22))*b0(x22);
   *p33 = (s + 0.25*v3 - (*u0)(x33))*b0(x33);
   *p112 = ((3.*v1 + v2)/4. - (*u0)(x112))*b0(x112);
   *p113 = ((3.*v1 + v3)/4. - (*u0)(x113))*b0(x113);
   *p221 = ((3.*v2 + v1)/4. - (*u0)(x221))*b0(x221);
   *p223 = ((3.*v2 + v3)/4. - (*u0)(x223))*b0(x223);
   *p331 = ((3.*v3 + v1)/4. - (*u0)(x331))*b0(x331);
   *p332 = ((3.*v3 + v2)/4. - (*u0)(x332))*b0(x332);
}

void L2_quad_node_values_4(
          x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
          p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332,
          v1,v2,v3,v12,v13,v23,u0)
FLOAT *x1, *x2, *x3, *x11, *x22, *x33, *x12, *x13, *x23, 
      *x112, *x113, *x221, *x223, *x331, *x332,
      *p1, *p2, *p3, *p11, *p22, *p33, *p12, *p13, *p23, 
      *p112, *p113, *p221, *p223, *p331, *p332, v1, v2, v3, v12, v13, v23, 
      (*u0)();                       
{
   FLOAT s;

   s = 0.5*(v12 + v13 + v23) - 0.125*(v1 + v2 + v3);
   *p1 = v1 - (*u0)(x1);
   *p2 = v2 - (*u0)(x2);
   *p3 = v3 - (*u0)(x3);
   *p12 = v12 - (*u0)(x12);
   *p13 = v13 - (*u0)(x13);
   *p23 = v23 - (*u0)(x23);
   *p11 = s + 0.125*v1 - 0.25*v23 - (*u0)(x11);
   *p22 = s + 0.125*v2 - 0.25*v13 - (*u0)(x22);
   *p33 = s + 0.125*v3 - 0.25*v12 - (*u0)(x33);
   *p112 = (3.*v1 - v2)/8. + 0.75*v12 - (*u0)(x112);
   *p113 = (3.*v1 - v3)/8. + 0.75*v13 - (*u0)(x113);
   *p221 = (3.*v2 - v1)/8. + 0.75*v12 - (*u0)(x221);
   *p223 = (3.*v2 - v3)/8. + 0.75*v23 - (*u0)(x223);
   *p331 = (3.*v3 - v1)/8. + 0.75*v13 - (*u0)(x331);
   *p332 = (3.*v3 - v2)/8. + 0.75*v23 - (*u0)(x332);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT oldL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + face bubbles                 */
FLOAT (*u0)();           /* j = 0, 1                                          */
INT u, j;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],nn1[DIM],nn2[DIM],nn3[DIM],
         p1,p2,p3,p12,p13,p23;
   NODE *n1, *n2, *n3;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   L2_lin_node_values_2(x1,x2,x3,x12,x13,x23,
                        &p1,&p2,&p3,&p12,&p13,&p23,
                        ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),u0);
   p12 += FD(ltop_face(pelement->f[2],tGrid),u)*FMULT*nn3[j]/4.;
   p13 += FD(ltop_face(pelement->f[1],tGrid),u)*FMULT*nn2[j]/4.;
   p23 += FD(ltop_face(pelement->f[0],tGrid),u)*FMULT*nn1[j]/4.;
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}
             
#else

FLOAT oldL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: oldL2s not available.\n"); return(0.);  }

#endif

#if (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT oldbubbleL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of */
GRID *tGrid;                           /* the function saved in u             */
ELEMENT *pelement;                     /* u = face bubbles                    */
FLOAT (*u0)();                         /* j = 0, 1                            */
INT u, j;
{
   FLOAT *x1,*x2,*x3,nn1[DIM],nn2[DIM],nn3[DIM],p12,p13,p23;
                     
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   p12 = FD(ltop_face(pelement->f[2],tGrid),u)*FMULT*nn3[j]/4.;
   p13 = FD(ltop_face(pelement->f[1],tGrid),u)*FMULT*nn2[j]/4.;
   p23 = FD(ltop_face(pelement->f[0],tGrid),u)*FMULT*nn1[j]/4.;
   return( p22(0.,0.,0.,p12,p13,p23)*volume(x1,x2,x3) );
}

#else

FLOAT oldbubbleL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: oldbubbleL2s not available.\n"); return(0.);  }

#endif

FLOAT L2s_for_subtriangle(x1,x2,x3,v1,v2,v3,nn3j,face_value,u0)
FLOAT *x1, *x2, *x3, v1, v2, v3, nn3j, face_value, (*u0)();
{
   FLOAT x12[DIM],x13[DIM],x23[DIM],p1,p2,p3,p12,p13,p23;
                     
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_node_values_2(x1,x2,x3,x12,x13,x23,
                        &p1,&p2,&p3,&p12,&p13,&p23,v1,v2,v3,u0);
   p12 += face_value*FMULT*nn3j/4.;
   return( p22(p1,p2,p3,p12,p13,p23) );
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT newL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + face bubbles                 */
FLOAT (*u0)();           /* j = 0, 1                                          */
INT u, j;
{
   FLOAT *x1,*x2,*x3,xc[DIM],nn1[DIM],nn2[DIM],nn3[DIM],v1,v2,v3,vc;
   NODE *n1, *n2, *n3;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   POINT3(x1,x2,x3,xc);
   v1 = ND(n1,u,j);
   v2 = ND(n2,u,j);
   v3 = ND(n3,u,j);
   vc = (v1 + v2 + v3)/3.;
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   return((L2s_for_subtriangle(x1,x2,xc,v1,v2,vc,nn3[j],
                                     FD(ltop_face(pelement->f[2],tGrid),u),u0)+
           L2s_for_subtriangle(x1,x3,xc,v1,v3,vc,nn2[j],
                                     FD(ltop_face(pelement->f[1],tGrid),u),u0)+
           L2s_for_subtriangle(x2,x3,xc,v2,v3,vc,nn1[j],
               FD(ltop_face(pelement->f[0],tGrid),u),u0))*volume(x1,x2,x3)/3.);
}

#else

FLOAT newL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: newL2s not available.\n"); return(0.);  }

#endif

#if (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT newbubbleL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of */
GRID *tGrid;                           /* the function saved in u             */
ELEMENT *pelement;                     /* u = face bubbles                    */
FLOAT (*u0)();                         /* j = 0, 1                            */
INT u, j;
{
   FLOAT *x1, *x2, *x3, nn1[DIM], nn2[DIM], nn3[DIM], fv1, fv2, fv3;
                     
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   fv1 = FD(ltop_face(pelement->f[2],tGrid),u)*nn3[j];
   fv2 = FD(ltop_face(pelement->f[1],tGrid),u)*nn2[j];
   fv3 = FD(ltop_face(pelement->f[0],tGrid),u)*nn1[j];
   return( (fv1*fv1 + fv2*fv2 + fv3*fv3)*volume(x1,x2,x3)*FMULT*FMULT/270.);
}

#else

FLOAT newbubbleL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: newbubbleL2s not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT vp1cL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + element bubbles              */
FLOAT (*u0)();           /* j = 0, 1                                          */
INT u, j;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123);
   L2_lin_node_values_3(x1, x2, x3, x112, x113, x221, x223, x331, x332,
                        &p1,&p2,&p3,&p112,&p113,&p221,&p223,&p331,&p332,
                        ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),u0);
   p123 = (ND(n1,u,j)+ND(n2,u,j)+ND(n3,u,j))/3.0 - (*u0)(x123);
   return( p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)*volume(x1,x2,x3) );
}

#else

FLOAT vp1cL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: vp1cL2s not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT miniL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin. + element bubbles              */
FLOAT (*u0)();           /* j = 0, 1                                          */
INT u, j;
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM],
         p1, p2, p3, p112, p113, p221, p223, p331, p332, p123;
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123);
   L2_lin_node_values_3(x1, x2, x3, x112, x113, x221, x223, x331, x332,
                        &p1,&p2,&p3,&p112,&p113,&p221,&p223,&p331,&p332,
                        ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),u0);
   p123 = (ND(n1,u,j)+ND(n2,u,j)+ND(n3,u,j))/3.0 + 
          EDV(pelement,u,j)*FMULT/27.0 -(*u0)(x123);
   return( p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123)*volume(x1,x2,x3) );
}

#else

FLOAT miniL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the  */
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: miniL2s not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT L2sq(tGrid,pelement,u0,u,j)  /* L2 norm of the j-th component of the   */
GRID *tGrid;            /* difference between u0 and a function saved in u   */
ELEMENT *pelement;
FLOAT (*u0)();
INT u, j;
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], 
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], 
         x221[DIM], x223[DIM], x331[DIM], x332[DIM], p1, p2, p3, 
         p11, p22, p33, p12, p13, p23, p112, p113, p221, p223, p331, p332,
         v1, v2, v3, v12, v13, v23;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   v1 = ND(n1,u,j); 
   v2 = ND(n2,u,j); 
   v3 = ND(n3,u,j); 
   v12 = 0.5*(v1+v2) + 0.25*FDV(f3,u,j);
   v13 = 0.5*(v1+v3) + 0.25*FDV(f2,u,j);
   v23 = 0.5*(v2+v3) + 0.25*FDV(f1,u,j);
   L2_quad_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,
                  x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
                  &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,
                  v1,v2,v3,v12,v13,v23,u0);
   return( p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,
               p112,p113,p221,p223,p331,p332)*volume(x1,x2,x3) );
}

#else

FLOAT L2sq(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: L2sq not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (F_DATA & CURVED_FACE_MIDDLE)

FLOAT L2sq_iso(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;
FLOAT (*u0)();
INT u, j;
{
   INT k;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT v0, v1, v2, v01, v02, v12, a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2],
         *x0, *x1, *x2, s=0, z;

   LTOP_C_DATA_OF_ELEMENT(n0,n1,n2,fa0,fa1,fa2,x0,x1,x2,pelement,tGrid)
   P2_reference_mapping0(n0,n1,n2,fa0,a,c,alpha,jac);
   v0 = ND(n0,u,j); 
   v1 = ND(n1,u,j); 
   v2 = ND(n2,u,j); 
   v01 = FDV(fa2,u,j);
   v02 = FDV(fa1,u,j);
   v12 = FDV(fa0,u,j);
   for (k=0; k < QR8_N; k++){
      z = REF_QUADR(QR8_P[k],v0,v1,v2,v01,v02,v12) - 
          fcn_bilin_value(QR8_P[k],a,c,alpha,u0);
      s += QR8_W[k]*z*z*fabs(LINV(jac,QR8_P[k]));
   }
   return(s*QR_VOL);
}

FLOAT H10q_iso(tGrid,pelement,u011,u012,u021,u022,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u011)(), (*u012)(), (*u021)(), (*u022)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT v0, v1, v2, v01, v02, v12, w0, w1, w2, w01, w02, w12, b[2][2][DIM2], 
         a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], *x0, *x1, *x2, x[DIM], 
         z, z11, z12, z21, z22, s=0;

   LTOP_C_DATA_OF_ELEMENT(n0,n1,n2,fa0,fa1,fa2,x0,x1,x2,pelement,tGrid)
   P2_reference_mapping0_with_inverse(n0,n1,n2,fa0,b,a,c,alpha,jac);
   v0 = ND(n0,u,0); 
   v1 = ND(n1,u,0); 
   v2 = ND(n2,u,0); 
   w0 = ND(n0,u,1); 
   w1 = ND(n1,u,1); 
   w2 = ND(n2,u,1); 
   v01 = FDV(fa2,u,0);
   v02 = FDV(fa1,u,0);
   v12 = FDV(fa0,u,0);
   w01 = FDV(fa2,u,1);
   w02 = FDV(fa1,u,1);
   w12 = FDV(fa0,u,1);
   for (k=0; k < QR8_N; k++){
      SET1(x,QR8_P[k])
      z = LINV(jac,x);
      z11 = (DX_REF_QUADR(x,v0,v1,v2,v01,v02,v12)*LINV(b[0][0],x) +
             DY_REF_QUADR(x,v0,v1,v2,v01,v02,v12)*LINV(b[1][0],x))/z
            - fcn_bilin_value(x,a,c,alpha,u011);
      z12 = (DX_REF_QUADR(x,v0,v1,v2,v01,v02,v12)*LINV(b[0][1],x) +
             DY_REF_QUADR(x,v0,v1,v2,v01,v02,v12)*LINV(b[1][1],x))/z
            - fcn_bilin_value(x,a,c,alpha,u012);
      z21 = (DX_REF_QUADR(x,w0,w1,w2,w01,w02,w12)*LINV(b[0][0],x) +
             DY_REF_QUADR(x,w0,w1,w2,w01,w02,w12)*LINV(b[1][0],x))/z
            - fcn_bilin_value(x,a,c,alpha,u021);
      z22 = (DX_REF_QUADR(x,w0,w1,w2,w01,w02,w12)*LINV(b[0][1],x) +
             DY_REF_QUADR(x,w0,w1,w2,w01,w02,w12)*LINV(b[1][1],x))/z
            - fcn_bilin_value(x,a,c,alpha,u022);
      s += QR8_W[k]*(z11*z11+z12*z12+z21*z21+z22*z22)*fabs(z);
   }
   return(s*QR_VOL);
}

#else

FLOAT L2sq_iso(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: L2sq_iso not available.\n"); return(0.);  }

FLOAT H10q_iso(tGrid,pelement,u011,u012,u021,u022,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u011)(), (*u012)(), (*u021)(), (*u022)(); INT u;
{  eprintf("Error: H10q_iso not available.\n"); return(0.);  }

#endif

#if (E_DATA & VECTOR_ELEMENT_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT minibubbleL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of*/
GRID *tGrid;                           /* the function saved in u             */
ELEMENT *pelement;                     /* u = element bubbles                 */
FLOAT (*u0)();                         /* j = 0, 1                            */
INT u, j;
{
   FLOAT *x1, *x2, *x3, a;

   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   a = EDV(pelement,u,j)*FMULT;
   return( a*a/2520.*volume(x1,x2,x3) ); 
}

#else

FLOAT minibubbleL2s(tGrid,pelement,u0,u,j) 
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: minibubbleL2s not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT linL2norms(tGrid,pelement,u,j)    /* L2 norm of the j-th component of   */
GRID *tGrid;                            /* the function saved in u            */
ELEMENT *pelement;                      /* u = cont. pw. lin.                 */
INT u, j;                               /* j = 0, 1                           */
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],
         p1,p2,p3,p12,p13,p23;
   NODE *n1, *n2, *n3;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   lin_node_values_2(x1,x2,x3,x12,x13,x23,
                     &p1,&p2,&p3,&p12,&p13,&p23,
                     ND(n1,u,j),ND(n2,u,j),ND(n3,u,j));
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

FLOAT linL2s(tGrid,pelement,u0,u,j) /* L2 norm of the j-th component of the   */
GRID *tGrid;             /* difference between u0 and a function saved in u   */
ELEMENT *pelement;       /* u = cont. pw. lin.                                */
FLOAT (*u0)();           /* j = 0, 1                                          */
INT u, j;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],
         p1,p2,p3,p12,p13,p23;
   NODE *n1, *n2, *n3;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_node_values_2(x1,x2,x3,x12,x13,x23,
                        &p1,&p2,&p3,&p12,&p13,&p23,
                        ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),u0);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

FLOAT lin_div(tGrid,pelement,u) /*  || div u_lin ||_{0,pel}^2  */
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT c, detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   detB = barycentric_coordinates(n3->myvertex->x,n1->myvertex->x,
                                                  n2->myvertex->x,b);
   c = DOT(NDD(n1,u),b[1]) + DOT(NDD(n2,u),b[2]) + DOT(NDD(n3,u),b[0]);
   return(c*c*detB);
}

#else

FLOAT linL2norms(tGrid,pelement,u,j)
GRID *tGrid; ELEMENT *pelement; INT u, j;
{  eprintf("Error: linL2norms not available.\n"); return(0.);  }

FLOAT linL2s(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: linL2s not available.\n"); return(0.);  }

FLOAT lin_div(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{  eprintf("Error: lin_div not available.\n"); return(0.);  }

#endif

#if ELEMENT_TYPE == SIMPLEX

FLOAT pL2s(pelement,p0,c) /*  ||c-p0||_{0,pelement}^2  */
ELEMENT *pelement;
FLOAT (*p0)(), c;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],
         p1,p2,p3,p12,p13,p23;

   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   p1 = c - (*p0)(x1);
   p2 = c - (*p0)(x2);
   p3 = c - (*p0)(x3);
   p12 = c - (*p0)(x12);
   p13 = c - (*p0)(x13);
   p23 = c - (*p0)(x23);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

FLOAT pL2s_lin(pelement,p0,c,q1,q2,q3) /* pressure on pelement consists of a  */
ELEMENT *pelement;                     /* constant part c and a pw. lin. part */
FLOAT (*p0)(), c, q1, q2, q3; /* q_i = value of lin. part in pelement->n[i-1] */
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],
         p1,p2,p3,p12,p13,p23;

   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   p1 = c + q1 - (*p0)(x1);
   p2 = c + q2 - (*p0)(x2);
   p3 = c + q3 - (*p0)(x3);
   p12 = c + (q1+q2)/2. - (*p0)(x12);
   p13 = c + (q1+q3)/2. - (*p0)(x13);
   p23 = c + (q2+q3)/2. - (*p0)(x23);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

#else

FLOAT pL2s(pelement,p0,c)
ELEMENT *pelement; FLOAT (*p0)(), c;
{  eprintf("Error: pL2s not available.\n"); return(0.);  }

FLOAT pL2s_lin(pelement,p0,c,q1,q2,q3)
ELEMENT *pelement; FLOAT (*p0)(), c, q1, q2, q3;
{  eprintf("Error: pL2s_lin not available.\n"); return(0.);  }

#endif

FLOAT pL2s_q(pelement,p0,c0) /*  ||c0-p0||_{0,pelement}^2  */
ELEMENT *pelement;           /*   pelement=quadrilateral   */
FLOAT (*p0)(), c0;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], s=0, z;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac);
   for (k=0; k < QS5_N; k++){
      z = c0 - fcn_bilin_value(QS5_P[k],a,c,alpha,p0);
      s += QS5_W[k]*z*z*fabs(LINV(jac,QS5_P[k]));
   }
   return(s*QR_VOL);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

/* (L2-error of j-th derivative of i-th comp.)^2 / |K| on one element */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.;
   nn1i is the i-th component of the normal vector to f1 etc. */  
FLOAT oldH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,
	    gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   p1 = c + FMULT*( FD(f2,u)*nn2i*gl3j + FD(f3,u)*nn3i*gl2j ) - (*u0ij)(x1);
   p2 = c + FMULT*( FD(f1,u)*nn1i*gl3j + FD(f3,u)*nn3i*gl1j ) - (*u0ij)(x2);
   p3 = c + FMULT*( FD(f1,u)*nn1i*gl2j + FD(f2,u)*nn2i*gl1j ) - (*u0ij)(x3);
   p12 = c + 0.5*FMULT*( (FD(f1,u)*nn1i + FD(f2,u)*nn2i)*gl3j + 
                          FD(f3,u)*nn3i*(gl1j+gl2j) ) - (*u0ij)(x12);
   p13 = c + 0.5*FMULT*( (FD(f1,u)*nn1i + FD(f3,u)*nn3i)*gl2j + 
                          FD(f2,u)*nn2i*(gl1j+gl3j) ) - (*u0ij)(x13);
   p23 = c + 0.5*FMULT*( (FD(f2,u)*nn2i + FD(f3,u)*nn3i)*gl1j + 
                          FD(f1,u)*nn1i*(gl2j+gl3j) ) - (*u0ij)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

#else

FLOAT oldH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i) 
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldH10ij_2 not available.\n"); return(0.);  }

#endif

FLOAT oldH10ij_3(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldH10ij_3 not available.\n"); return(0.);  }
 
#if F_DATA & SCALAR_FACE_DATA

FLOAT oldbubbleH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,
	          gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1,p2,p3,p12,p13,p23;
     
   p1 = FMULT*( FD(f2,u)*nn2i*gl3j + FD(f3,u)*nn3i*gl2j );
   p2 = FMULT*( FD(f1,u)*nn1i*gl3j + FD(f3,u)*nn3i*gl1j );
   p3 = FMULT*( FD(f1,u)*nn1i*gl2j + FD(f2,u)*nn2i*gl1j );
   p12 = 0.5*FMULT*( (FD(f1,u)*nn1i + FD(f2,u)*nn2i)*gl3j + 
                      FD(f3,u)*nn3i*(gl1j+gl2j) );
   p13 = 0.5*FMULT*( (FD(f1,u)*nn1i + FD(f3,u)*nn3i)*gl2j + 
                      FD(f2,u)*nn2i*(gl1j+gl3j) );
   p23 = 0.5*FMULT*( (FD(f2,u)*nn2i + FD(f3,u)*nn3i)*gl1j + 
                      FD(f1,u)*nn1i*(gl2j+gl3j) );
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

#else

FLOAT oldbubbleH10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldbubbleH10ij_2 not available.\n"); return(0.);  }

#endif

FLOAT oldbubbleH10ij_3(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: oldbubbleH10ij_3 not available.\n"); return(0.);  }

FLOAT newH10ij_2(x1,x2,x3,x12,x13,x23,v1,v2,v3,face_val,gl1j,gl2j,gl3j,u0ij)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, v1, v2, v3, face_val, gl1j, gl2j, gl3j;
FLOAT (*u0ij)();  /* x3 is barycentre,  face_val = FMULT*FD(f3,u)*nn3i */
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = v1*gl1j + v2*gl2j + v3*gl3j;
   p1 = c + face_val*gl2j - (*u0ij)(x1);
   p2 = c + face_val*gl1j - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p12 = c + 0.5*face_val*(gl1j+gl2j) - (*u0ij)(x12);
   p13 = c + 0.5*face_val*gl2j - (*u0ij)(x13);
   p23 = c + 0.5*face_val*gl1j - (*u0ij)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

FLOAT newbubbleH10ij_2(x1,x2,x3,x12,x13,x23,v1,v2,v3,face_val,gl1j,gl2j,gl3j,u0ij)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, v1, v2, v3, face_val, gl1j, gl2j, gl3j;
FLOAT (*u0ij)();  /* x3 is barycentre,  face_val = FMULT*FD(f3,u)*nn3i */
{
   FLOAT p1,p2,p3,p12,p13,p23;
     
   p1 = face_val*gl2j;
   p2 = face_val*gl1j;
   p3 = 0.;
   p12 = 0.5*face_val*(gl1j+gl2j);
   p13 = 0.5*face_val*gl2j;
   p23 = 0.5*face_val*gl1j;
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

/* H10 norm on a subelement                                                   */
FLOAT H10_for_sub(x1,x2,x3,v10,v20,v30,v11,v21,v31,face_val,nn3,H10ij_2)
FLOAT *x1, *x2, *x3, v10, v20, v30, v11, v21, v31, face_val, *nn3, (*H10ij_2)();
{
   FLOAT x12[DIM], x13[DIM], x23[DIM], 
         face_val0, face_val1, detB, b[DIM2][DIM2];

   detB = barycentric_coordinates(x3,x1,x2,b);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   face_val0 = face_val*FMULT*nn3[0];
   face_val1 = face_val*FMULT*nn3[1];

   return((H10ij_2(x1,x2,x3,x12,x13,x23,
           v10,v20,v30,face_val0,b[1][0],b[2][0],b[0][0],u011) +
   H10ij_2(x1,x2,x3,x12,x13,x23,
           v10,v20,v30,face_val0,b[1][1],b[2][1],b[0][1],u012) +
   H10ij_2(x1,x2,x3,x12,x13,x23,
           v11,v21,v31,face_val1,b[1][0],b[2][0],b[0][0],u021) +
   H10ij_2(x1,x2,x3,x12,x13,x23,
           v11,v21,v31,face_val1,b[1][1],b[2][1],b[0][1],u022))*detB);
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT H10s(tGrid,pelement,u)
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1, *x2, *x3, xc[DIM], nn1[DIM], nn2[DIM], nn3[DIM], 
         v10, v20, v30, vc0, v11, v21, v31, vc1;
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   POINT3(x1,x2,x3,xc);
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   v10 = ND(n1,u,0);
   v20 = ND(n2,u,0);
   v30 = ND(n3,u,0);
   v11 = ND(n1,u,1);
   v21 = ND(n2,u,1);
   v31 = ND(n3,u,1);
   vc0 = (v10 + v20 + v30)/3.;
   vc1 = (v11 + v21 + v31)/3.;
   
   return(H10_for_sub(x1,x2,xc,v10,v20,vc0,v11,v21,vc1,
                         FD(ltop_face(pelement->f[2],tGrid),u),nn3,newH10ij_2) +
          H10_for_sub(x1,x3,xc,v10,v30,vc0,v11,v31,vc1,
                         FD(ltop_face(pelement->f[1],tGrid),u),nn2,newH10ij_2) +
          H10_for_sub(x2,x3,xc,v20,v30,vc0,v21,v31,vc1,
                         FD(ltop_face(pelement->f[0],tGrid),u),nn1,newH10ij_2));
}

FLOAT bubbleH10(tGrid,pelement,u)
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1, *x2, *x3, xc[DIM], nn1[DIM], nn2[DIM], nn3[DIM];
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   POINT3(x1,x2,x3,xc);
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   return(H10_for_sub(x1,x2,xc,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[2],tGrid),u),nn3,newbubbleH10ij_2) +
          H10_for_sub(x1,x3,xc,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[1],tGrid),u),nn2,newbubbleH10ij_2) +
          H10_for_sub(x2,x3,xc,0.,0.,0.,0.,0.,0.,
                   FD(ltop_face(pelement->f[0],tGrid),u),nn1,newbubbleH10ij_2));
}

#else

FLOAT H10s(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{  eprintf("Error: H10s not available.\n"); return(0.);  }

FLOAT bubbleH10(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{  eprintf("Error: bubbleH10 not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

FLOAT miniH10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
	    gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, 
      gl1j, gl2j, gl3j;
NODE *n1, *n2, *n3;
ELEMENT *pelement;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1, p2, p3, p112, p113, p221, p223, p331, p332, p123, c, d;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   d = EDV(pelement,u,i)*2./9.;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p123 = c - (*u0ij)(x123);
   p112 = c + d*gl3j - (*u0ij)(x112);
   p221 = c + d*gl3j - (*u0ij)(x221);
   p223 = c + d*gl1j - (*u0ij)(x223);
   p332 = c + d*gl1j - (*u0ij)(x332);
   p113 = c + d*gl2j - (*u0ij)(x113);
   p331 = c + d*gl2j - (*u0ij)(x331);
   return( p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123) );
}  
 
#else

FLOAT miniH10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, gl1j, gl2j, gl3j; NODE *n1, *n2, *n3; ELEMENT *pelement; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: miniH10ij_3 not available.\n"); return(0.);  }

#endif

#if E_DATA & VECTOR_ELEMENT_DATA

FLOAT minibubbleH10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,
                      pelement,gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, 
      gl1j, gl2j, gl3j;
NODE *n1, *n2, *n3;
ELEMENT *pelement;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p112, p113, p221, p223, p331, p332, d;
     
   d = EDV(pelement,u,i)*2./9.;
   p112 = p221 = d*gl3j;
   p223 = p332 = d*gl1j;
   p113 = p331 = d*gl2j;
   return( p32(0.,0.,0.,p112,p113,p221,p223,p331,p332,0.) );
}  

#else

FLOAT minibubbleH10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, gl1j, gl2j, gl3j; NODE *n1, *n2, *n3; ELEMENT *pelement; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: minibubbleH10ij_3 not available.\n"); return(0.);  }

#endif

#if N_DATA & VECTOR_NODE_DATA

FLOAT linH10ijnorm_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,
	             gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT c;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   return( c*c );
}  
 
FLOAT linH10ijnorm_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
	           gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, 
      gl1j, gl2j, gl3j;
NODE *n1, *n2, *n3;
ELEMENT *pelement;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT c;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   return( c*c );
}  
 
FLOAT H10ijlin_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,
	    gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p12 = c - (*u0ij)(x12);
   p13 = c - (*u0ij)(x13);
   p23 = c - (*u0ij)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

FLOAT H10ijlin_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
	       gl1j,gl2j,gl3j,u0ij,u,i) /* i=0,1 */
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, 
      gl1j, gl2j, gl3j;
NODE *n1, *n2, *n3;
ELEMENT *pelement;
FLOAT (*u0ij)();
INT u,i;
{
   FLOAT p1, p2, p3, p112, p113, p221, p223, p331, p332, p123, c;
     
   c = ND(n1,u,i)*gl1j + ND(n2,u,i)*gl2j + ND(n3,u,i)*gl3j;
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p123 = c - (*u0ij)(x123);
   p112 = c - (*u0ij)(x112);
   p221 = c - (*u0ij)(x221);
   p223 = c - (*u0ij)(x223);
   p332 = c - (*u0ij)(x332);
   p113 = c - (*u0ij)(x113);
   p331 = c - (*u0ij)(x331);
   return( p32(p1,p2,p3,p112,p113,p221,p223,p331,p332,p123) );
}  
 
#else

FLOAT linH10ijnorm_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j; NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: linH10ijnorm_2 not available.\n"); return(0.);  }

FLOAT linH10ijnorm_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, gl1j, gl2j, gl3j; NODE *n1, *n2, *n3; ELEMENT *pelement; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: linH10ijnorm_3 not available.\n"); return(0.);  }

FLOAT H10ijlin_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,nn1i,nn2i,nn3i,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, nn1i,nn2i,nn3i,gl1j,gl2j,gl3j;
NODE *n1, *n2, *n3; FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: H10ijlin_2 not available.\n"); return(0.);  }

FLOAT H10ijlin_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x112, *x113, *x221, *x223, *x331, *x332, *x123, gl1j, gl2j, gl3j; NODE *n1, *n2, *n3; ELEMENT *pelement; FLOAT (*u0ij)(); INT u,i;
{  eprintf("Error: H10ijlin_3 not available.\n"); return(0.);  }

#endif

FLOAT H10_2(tGrid,pelement,u,H10ij_2)
GRID *tGrid;
ELEMENT *pelement;
INT u;
FLOAT (*H10ij_2)();
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],
                     nn1[DIM],nn2[DIM],nn3[DIM],
                     detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   NORMAL_VECTORS(nn1,nn2,nn3,pelement);
   detB = barycentric_coordinates(x3,x1,x2,b);
   
   return((H10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,
          nn1[0],nn2[0],nn3[0],b[1][0],b[2][0],b[0][0],u011,u,0) +
   H10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,
          nn1[0],nn2[0],nn3[0],b[1][1],b[2][1],b[0][1],u012,u,0) +
   H10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,
          nn1[1],nn2[1],nn3[1],b[1][0],b[2][0],b[0][0],u021,u,1) +
   H10ij_2(x1,x2,x3,x12,x13,x23,n1,n2,n3,f1,f2,f3,
          nn1[1],nn2[1],nn3[1],b[1][1],b[2][1],b[0][1],u022,u,1))*detB);
}
 
FLOAT H10_3(tGrid,pelement,u,H10ij_3)
GRID *tGrid;
ELEMENT *pelement;
INT u;
FLOAT (*H10ij_3)();
{
   FLOAT *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], x123[DIM], detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123);
   detB = barycentric_coordinates(x3,x1,x2,b);
   
   return((H10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
          b[1][0],b[2][0],b[0][0],u011,u,0) +
   H10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
          b[1][1],b[2][1],b[0][1],u012,u,0) +
   H10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
          b[1][0],b[2][0],b[0][0],u021,u,1) +
   H10ij_3(x1,x2,x3,x112,x113,x221,x223,x331,x332,x123,n1,n2,n3,pelement,
          b[1][1],b[2][1],b[0][1],u022,u,1))*detB);
}
 
/* (\int_K g \dx)/|K|; exact for g\in\P_3 */
FLOAT integr5(pelement,g)
ELEMENT *pelement;
FLOAT (*g)();
{
   FLOAT *x1, *x2, *x3, x12[DIM], x13[DIM], x23[DIM], x123[DIM];

   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   POINT3(x1,x2,x3,x123);
   return(( 3.*((*g)(x1) + (*g)(x2) + (*g)(x3)) + 
            8.0*((*g)(x12) + (*g)(x13) + (*g)(x23)) + 
            27.0*(*g)(x123) )/60.0);
}

FLOAT integr_q5(pelement,g)
ELEMENT *pelement;
FLOAT (*g)();
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], s=0;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac);
   for (k=0; k < QS5_N; k++)
      s += QS5_W[k]*fcn_bilin_value(QS5_P[k],a,c,alpha,g)*
                                                       fabs(LINV(jac,QS5_P[k]));
   return(s*QR_VOL/VOLUME(pelement));
}

#if (N_DATA & SCALAR_NODE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT sL2s(tGrid,pelement,u0,u) /* L2 norm of the difference between u0       */
GRID *tGrid;                    /* and a scalar pw. lin. function saved in u  */
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],p1,p2,p3,p12,p13,p23;
   NODE *n1, *n2, *n3;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_node_values_2(x1,x2,x3,x12,x13,x23,
                        &p1,&p2,&p3,&p12,&p13,&p23,
                        NDS(n1,u),NDS(n2,u),NDS(n3,u),u0);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

/* (L2-error of j-th derivative)^2 / |K| on one element                       */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.          */
FLOAT sH10j(x1,x2,x3,x12,x13,x23,n1,n2,n3,gl1j,gl2j,gl3j,u0j,u)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
NODE *n1, *n2, *n3;
FLOAT (*u0j)();
INT u;
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = NDS(n1,u)*gl1j + NDS(n2,u)*gl2j + NDS(n3,u)*gl3j;
   p1 = c - (*u0j)(x1);
   p2 = c - (*u0j)(x2);
   p3 = c - (*u0j)(x3);
   p12 = c - (*u0j)(x12);
   p13 = c - (*u0j)(x13);
   p23 = c - (*u0j)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

/* (H10-error on one element)^2 */
FLOAT sH10(tGrid,pelement,u01,u02,u03,u)   /*  u03 not used  */
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u; 
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM], detB, b[DIM2][DIM2];
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   detB = barycentric_coordinates(x3,x1,x2,b);
   
   return(
     (sH10j(x1,x2,x3,x12,x13,x23,n1,n2,n3,b[1][0],b[2][0],b[0][0],u01,u) +
      sH10j(x1,x2,x3,x12,x13,x23,n1,n2,n3,b[1][1],b[2][1],b[0][1],u02,u))*detB);
}

/* (error of  bb*grad u  on one element)^2 */
FLOAT p1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid;                                      /*  bb02 and u03 not used  */
ELEMENT *pelement;
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
INT u;
{
   FLOAT *x1, *x2, *x3, x12[DIM], x13[DIM], x23[DIM], detB, b[DIM2][DIM2],
         p1, p2, p3, p12, p13, p23, g0, g1;
   NODE *n1, *n2, *n3;
   
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   detB = barycentric_coordinates(x3,x1,x2,b);

   g0 = NDS(n1,u)*b[1][0] + NDS(n2,u)*b[2][0] + NDS(n3,u)*b[0][0];
   g1 = NDS(n1,u)*b[1][1] + NDS(n2,u)*b[2][1] + NDS(n3,u)*b[0][1];

   p1  = bb0(x1 )*(g0 - (*u01)(x1 )) + bb1(x1 )*(g1 - (*u02)(x1));
   p2  = bb0(x2 )*(g0 - (*u01)(x2 )) + bb1(x2 )*(g1 - (*u02)(x2));
   p3  = bb0(x3 )*(g0 - (*u01)(x3 )) + bb1(x3 )*(g1 - (*u02)(x3));
   p12 = bb0(x12)*(g0 - (*u01)(x12)) + bb1(x12)*(g1 - (*u02)(x12));
   p13 = bb0(x13)*(g0 - (*u01)(x13)) + bb1(x13)*(g1 - (*u02)(x13));
   p23 = bb0(x23)*(g0 - (*u01)(x23)) + bb1(x23)*(g1 - (*u02)(x23));

   return( p22(p1,p2,p3,p12,p13,p23)*detB );
}

FLOAT p1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, d, max;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   max = fabs(NDS(n1,u) - (*u0)(x1));
   d   = fabs(NDS(n2,u) - (*u0)(x2));
   if (d > max) max = d;
   d   = fabs(NDS(n3,u) - (*u0)(x3));
   if (d > max) max = d;
   return(max);
}

#else

FLOAT sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2s not available.\n"); return(0.);  }

FLOAT sH10(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u; 
{  eprintf("Error: sH10 not available.\n"); return(0.);  }

FLOAT p1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{ eprintf("Error: p1c_L2conv not available.\n"); return(0.); }

FLOAT p1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{ eprintf("Error: p1c_Linf not available.\n"); return(0.); }

#endif  /* N_DATA & SCALAR_NODE_DATA */
 
/* (L2-error of j-th derivative)^2 / |K| on one element                       */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.          */
FLOAT sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
             u1,u2,u3,u12,u13,u23,gl1j,gl2j,gl3j,u0j,u)
FLOAT *x1, *x2, *x3, *x11, *x22, *x33, *x12, *x13, *x23, 
      *x112, *x113, *x221, *x223, *x331, *x332, u1, u2, u3, u12, u13, u23, 
      gl1j, gl2j, gl3j, (*u0j)();
INT u;
{
   FLOAT p1, p2, p3, p11, p22, p33, p12, p13, p23, 
         p112, p113, p221, p223, p331, p332, v1, v2, v3, c;
     
   c = u1*gl1j + u2*gl2j + u3*gl3j;
   v1 = c + u12*gl2j + u13*gl3j;
   v2 = c + u12*gl1j + u23*gl3j;
   v3 = c + u13*gl1j + u23*gl2j;
   L2_lin_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,
               x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
               &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,v1,v2,v3,u0j);
   return( p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332) );
}  

#if (N_DATA & SCALAR_NODE_DATA) && (ELEMENT_TYPE == CUBE)

FLOAT sL2_q1(tGrid,pelement,u0,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], s=0, z;

   LTOPNODES_OF_4ELEMENT(n0,n1,n2,n3,pelement,tGrid);
   Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac);
   for (k=0; k < QS5_N; k++){
      z = REF_BILIN(QS5_P[k],NDS(n0,u),NDS(n1,u),NDS(n2,u),NDS(n3,u)) - 
          fcn_bilin_value(QS5_P[k],a,c,alpha,u0);
      s += QS5_W[k]*z*z*fabs(LINV(jac,QS5_P[k]));
   }
   return(s*QR_VOL);
}

FLOAT sH10_q1(tGrid,pelement,u01,u02,u03,u)   /*  u03 not used  */
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FLOAT v0, v1, v2, v3, b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], 
         jac[DIM2], x[DIM], z, z0, z1, s=0;

   LTOPNODES_OF_4ELEMENT(n0,n1,n2,n3,pelement,tGrid);
   Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac);
   v0 = NDS(n0,u); 
   v1 = NDS(n1,u); 
   v2 = NDS(n2,u); 
   v3 = NDS(n3,u); 
   for (k=0; k < QS5_N; k++){
      SET1(x,QS5_P[k])
      z = LINV(jac,x);
      z0 = (DX_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[0][0],x) +
            DY_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[1][0],x))/z
            - fcn_bilin_value(x,a,c,alpha,u01);
      z1 = (DX_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[0][1],x) +
            DY_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[1][1],x))/z
            - fcn_bilin_value(x,a,c,alpha,u02);
      s += QS5_W[k]*(z0*z0+z1*z1)*fabs(z);
   }
   return(s*QR_VOL);
}

/* (error of  bb*grad u  on one element)^2 */
FLOAT q1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid;                                        /*  bb2 and u03 not used  */
ELEMENT *pelement;
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FLOAT v0, v1, v2, v3, b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], 
         jac[DIM2], x[DIM], z, z0, s=0;

   LTOPNODES_OF_4ELEMENT(n0,n1,n2,n3,pelement,tGrid);
   Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac);
   v0 = NDS(n0,u); 
   v1 = NDS(n1,u); 
   v2 = NDS(n2,u); 
   v3 = NDS(n3,u); 
   for (k=0; k < QS5_N; k++){
      SET1(x,QS5_P[k])
      z = LINV(jac,x);
      z0 = ((DX_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[0][0],x) +
             DY_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[1][0],x))/z
           - fcn_bilin_value(x,a,c,alpha,u01))*fcn_bilin_value(x,a,c,alpha,bb0)
         + ((DX_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[0][1],x) +
             DY_REF_BILIN(x,v0,v1,v2,v3)*LINV(b[1][1],x))/z
           - fcn_bilin_value(x,a,c,alpha,u02))*fcn_bilin_value(x,a,c,alpha,bb1);
      s += QS5_W[k]*z0*z0*fabs(z);
   }
   return(s*QR_VOL);
}

FLOAT q1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   NODE *n0, *n1, *n2, *n3;
   FLOAT *x0, *x1, *x2, *x3, d, max;
                     
   LTOPNODES_OF_4ELEMENT(n0,n1,n2,n3,pelement,tGrid);
   VERTICES_OF_4ELEMENT(x0,x1,x2,x3,pelement);
   max = fabs(NDS(n0,u) - (*u0)(x0));
   d   = fabs(NDS(n1,u) - (*u0)(x1));
   if (d > max) max = d;
   d   = fabs(NDS(n2,u) - (*u0)(x2));
   if (d > max) max = d;
   d   = fabs(NDS(n3,u) - (*u0)(x3));
   if (d > max) max = d;
   return(max);
}


#else

FLOAT sL2_q1(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2_q1 not available.\n"); return(0.);  }

FLOAT sH10_q1(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: sH10_q1 not available.\n"); return(0.);  }

FLOAT q1c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: q1c_L2conv not available.\n"); return(0.);  }

FLOAT q1c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: q1c_Linf not available.\n"); return(0.);  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT sL2sq(tGrid,pelement,u0,u) /* L2 norm of the difference between u0      */
GRID *tGrid;                     /* and a scalar pw. lin. function saved in u */
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], 
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], 
         x221[DIM], x223[DIM], x331[DIM], x332[DIM], p1, p2, p3, 
         p11, p22, p33, p12, p13, p23, p112, p113, p221, p223, p331, p332,
         v1, v2, v3, v12, v13, v23;
                     
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   v1 = NDS(n1,u); 
   v2 = NDS(n2,u); 
   v3 = NDS(n3,u); 
   v12 = 0.5*(v1+v2) + 0.25*FD(f3,u);
   v13 = 0.5*(v1+v3) + 0.25*FD(f2,u);
   v23 = 0.5*(v2+v3) + 0.25*FD(f1,u);
   L2_quad_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,
                  x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
                  &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,
                  v1,v2,v3,v12,v13,v23,u0);
   return( p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,
               p112,p113,p221,p223,p331,p332)*volume(x1,x2,x3) );
}

/* (H10-error on one element)^2 */
FLOAT sH10q(tGrid,pelement,u01,u02,u03,u)    /*  u03 not used  */
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u;
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], 
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], 
         x221[DIM], x223[DIM], x331[DIM], x332[DIM], detB, b[DIM2][DIM2];
 
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   detB = barycentric_coordinates(x3,x1,x2,b);
   return((sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
             NDS(n1,u),NDS(n2,u),NDS(n3,u),FD(f3,u),FD(f2,u),FD(f1,u),
             b[1][0],b[2][0],b[0][0],u01,u) +
           sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
             NDS(n1,u),NDS(n2,u),NDS(n3,u),FD(f3,u),FD(f2,u),FD(f1,u),
             b[1][1],b[2][1],b[0][1],u02,u))*detB);
}

/* (error of  bb*grad u  on one element)^2 */
FLOAT p2c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid;                                        /*  bb2 and u03 not used  */
ELEMENT *pelement;
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
INT u;
{
   NODE *n0, *n1, *n2;
   FACE *f0, *f1, *f2;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], x12[DIM], x13[DIM], 
         x23[DIM], x112[DIM], x113[DIM], x221[DIM], x223[DIM], x331[DIM], 
         x332[DIM], p1, p2, p3, p11, p22, p33, p12, p13, p23, p112, p113, p221, 
         p223, p331, p332, v0[DIM], v1[DIM], v2[DIM], c0, c1, detB, 
         b[DIM2][DIM2];
 
   LTOPNODES_OF_ELEMENT(n0,n1,n2,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f0,f1,f2,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   detB = barycentric_coordinates(x1,x2,x3,b);
   c0 = NDS(n0,u)*b[0][0] + NDS(n1,u)*b[1][0] + NDS(n2,u)*b[2][0];
   c1 = NDS(n0,u)*b[0][1] + NDS(n1,u)*b[1][1] + NDS(n2,u)*b[2][1];
   v0[0] = c0 + FD(f2,u)*b[1][0] + FD(f1,u)*b[2][0];
   v0[1] = c1 + FD(f2,u)*b[1][1] + FD(f1,u)*b[2][1];
   v1[0] = c0 + FD(f2,u)*b[0][0] + FD(f0,u)*b[2][0];
   v1[1] = c1 + FD(f2,u)*b[0][1] + FD(f0,u)*b[2][1];
   v2[0] = c0 + FD(f1,u)*b[0][0] + FD(f0,u)*b[1][0];
   v2[1] = c1 + FD(f1,u)*b[0][1] + FD(f0,u)*b[1][1];
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   L2_mult_lin_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,
               x112,x113,x221,x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
               &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,
               v0[0],v1[0],v2[0],u01,bb0);
   c0 = p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332);
   L2_mult_lin_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,
               x112,x113,x221,x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
               &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,
               v0[1],v1[1],v2[1],u02,bb1);
   c1 = p42(p1,p2,p3,p11,p22,p33,p12,p13,p23,p112,p113,p221,p223,p331,p332);
   return((c0+c1)*detB);
}  

FLOAT p2c_Linf(tGrid,pelement,u0,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], u01, u02, u12;
                     
   LTOPNODES_OF_ELEMENT(n0,n1,n2,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(fa0,fa1,fa2,pelement,tGrid);
   VERTICES_OF_ELEMENT(x0,x1,x2,pelement);
   MIDPOINTS(x0,x1,x2,x01,x02,x12)
   u01 = 0.5*(NDS(n0,u)+NDS(n1,u)) + 0.25*FD(fa2,u);
   u02 = 0.5*(NDS(n0,u)+NDS(n2,u)) + 0.25*FD(fa1,u);
   u12 = 0.5*(NDS(n1,u)+NDS(n2,u)) + 0.25*FD(fa0,u);
   return(abs_max_of_6(NDS(n0,u)-(*u0)(x0),
                       NDS(n1,u)-(*u0)(x1),
                       NDS(n2,u)-(*u0)(x2),
                       u01 - (*u0)(x01),
                       u02 - (*u0)(x02),
                       u12 - (*u0)(x12)));
}

#else

FLOAT sL2sq(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2sq not available.\n"); return(0.);  }

FLOAT sH10q(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: sH10q not available.\n"); return(0.);  }

FLOAT p2c_L2conv(tGrid,pelement,bb0,bb1,bb2,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: p2c_L2conv not available.\n"); return(0.);  }

FLOAT p2c_Linf(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: p2c_Linf not available.\n"); return(0.);  }

#endif  /*  (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)  */
 
#if (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT nc_sL2s(tGrid,pelement,u0,u)   /* L2 norm of the difference between u0  */
GRID *tGrid;      /* and a nonconforming scalar pw. lin. function saved in u  */
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],p1,p2,p3,p12,p13,p23;
   FACE *f1, *f2, *f3;
                     
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_face_values(x1,x2,x3,x12,x13,x23,
                      &p1,&p2,&p3,&p12,&p13,&p23,
                      FD(f3,u),FD(f2,u),FD(f1,u),u0);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

/* (error of  bb*grad u  on one element)^2 */
FLOAT nc_L2conv(tGrid,pelement,bb0,bb1,u01,u02,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*bb0)(), (*bb1)(), (*u01)(), (*u02)();
INT u;
{
   FLOAT *x1, *x2, *x3, x12[DIM], x13[DIM], x23[DIM], detB, b[DIM2][DIM2],
         p1, p2, p3, p12, p13, p23, g0, g1;
   FACE *f1, *f2, *f3;
   
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   detB = barycentric_coordinates(x3,x1,x2,b);

   g0 = -2.*(FD(f1,u)*b[1][0] + FD(f2,u)*b[2][0] + FD(f3,u)*b[0][0]);
   g1 = -2.*(FD(f1,u)*b[1][1] + FD(f2,u)*b[2][1] + FD(f3,u)*b[0][1]);

   p1 = bb0(x1)*(g0 - (*u01)(x1)) + bb1(x1)*(g1 - (*u02)(x1));
   p2 = bb0(x2)*(g0 - (*u01)(x2)) + bb1(x2)*(g1 - (*u02)(x2));
   p3 = bb0(x3)*(g0 - (*u01)(x3)) + bb1(x3)*(g1 - (*u02)(x3));
   p12 = bb0(x12)*(g0 - (*u01)(x12)) + bb1(x12)*(g1 - (*u02)(x12));
   p13 = bb0(x13)*(g0 - (*u01)(x13)) + bb1(x13)*(g1 - (*u02)(x13));
   p23 = bb0(x23)*(g0 - (*u01)(x23)) + bb1(x23)*(g1 - (*u02)(x23));

   return( p22(p1,p2,p3,p12,p13,p23)*detB );
}

/* (L2-error of j-th derivative)^2 / |K| on one element                       */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.          */
FLOAT nc_sH10j(x1,x2,x3,x12,x13,x23,f1,f2,f3,gl1j,gl2j,gl3j,u0j,u)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
FACE *f1, *f2, *f3;
FLOAT (*u0j)();
INT u;
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = -2.*(FD(f1,u)*gl1j + FD(f2,u)*gl2j + FD(f3,u)*gl3j);
   p1 = c - (*u0j)(x1);
   p2 = c - (*u0j)(x2);
   p3 = c - (*u0j)(x3);
   p12 = c - (*u0j)(x12);
   p13 = c - (*u0j)(x13);
   p23 = c - (*u0j)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}

void nc_L_inf_error(tGrid,u0,u,err)
GRID *tGrid;
FLOAT (*u0)(), *err;
INT u;
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM], d, max=0.;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pel);
         VERTICES_OF_ELEMENT(x1,x2,x3,pel);
         MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
         d = fabs(FD(f1,u) - (*u0)(xc1));
         if (d > max) max = d;
         d = fabs(FD(f2,u) - (*u0)(xc2));
         if (d > max) max = d;
         d = fabs(FD(f3,u) - (*u0)(xc3));
         if (d > max) max = d;
      }
   printf("L-infinity error: %12.6e\n",(*err=max));
}

#else

FLOAT nc_sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: nc_sL2s not available.\n"); return(0.);  }

FLOAT nc_L2conv(tGrid,pelement,bb0,bb1,u01,u02,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*bb0)(), (*bb1)(), (*u01)(), (*u02)(); INT u;
{  eprintf("Error: nc_L2conv not available.\n"); return(0.);  }

FLOAT nc_sH10j(x1,x2,x3,x12,x13,x23,f1,f2,f3,gl1j,gl2j,gl3j,u0j,u)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
FACE *f1, *f2, *f3; FLOAT (*u0j)(); INT u;
{  eprintf("Error: nc_sH10j not available.\n"); return(0.);  }

void nc_L_inf_error(tGrid,u0,u,err)
GRID *tGrid; FLOAT (*u0)(), *err; INT u;
{  eprintf("Error: nc_L_inf_error not available.\n");  }

#endif

/* (H10-error on one element)^2 */
FLOAT nc_sH10(tGrid,pelement,u01,u02,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)();
INT u;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM], detB, b[DIM2][DIM2];
   FACE *f1, *f2, *f3;
   
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   detB = barycentric_coordinates(x3,x1,x2,b);
   
   return(
  (nc_sH10j(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][0],b[2][0],b[0][0],u01,u) +
   nc_sH10j(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][1],b[2][1],b[0][1],u02,u))*detB);
}

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT nc_L2s(tGrid,pelement,u0,u,i) /* L2 norm of the i-th component of the   */
GRID *tGrid;            /* difference between u0 and a nonconforming scalar   */
ELEMENT *pelement;      /* pw. lin. function saved in u                       */
FLOAT (*u0)();          /* i = 0, 1                                           */
INT u, i;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],p1,p2,p3,p12,p13,p23;
   FACE *f1, *f2, *f3;
                     
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_face_values(x1,x2,x3,x12,x13,x23,
                      &p1,&p2,&p3,&p12,&p13,&p23,
                      FDV(f3,u,i),FDV(f2,u,i),FDV(f1,u,i),u0);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

/* (L2-error of j-th derivative of the i-th component)^2 / |K| on one element */
/* gl1j is the j-th derivative of the barycentric coordinate l1 etc.          */
FLOAT nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,gl1j,gl2j,gl3j,u0ij,u,i)/*i=0,1*/
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
FACE *f1, *f2, *f3;
FLOAT (*u0ij)();
INT u, i;
{
   FLOAT p1,p2,p3,p12,p13,p23,c;
     
   c = -2.*(FDV(f1,u,i)*gl1j + FDV(f2,u,i)*gl2j + FDV(f3,u,i)*gl3j);
   p1 = c - (*u0ij)(x1);
   p2 = c - (*u0ij)(x2);
   p3 = c - (*u0ij)(x3);
   p12 = c - (*u0ij)(x12);
   p13 = c - (*u0ij)(x13);
   p23 = c - (*u0ij)(x23);
   return( p22(p1,p2,p3,p12,p13,p23) );
}  

FLOAT nc_lin_div(tGrid,pelement,u) /*  || div u_lin ||_{0,pel}^2  */
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1,*x2,*x3, detB, b[DIM2][DIM2], c;
   FACE *f1, *f2, *f3;
   
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   detB = barycentric_coordinates(x3,x1,x2,b);

   c = DOT(FDVP(f1,u),b[1]) + DOT(FDVP(f2,u),b[2]) + DOT(FDVP(f3,u),b[0]);
   return(4.*c*c*detB);
}

#else  /*  !(F_DATA & VECTOR_FACE_DATA)  */

FLOAT nc_L2s(tGrid,pelement,u0,u,i)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, i;
{ eprintf("Error: nc_L2s not available.\n"); return(0.);  }

FLOAT nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,gl1j,gl2j,gl3j,u0ij,u,i)
FLOAT *x1, *x2, *x3, *x12, *x13, *x23, gl1j, gl2j, gl3j;
FACE *f1, *f2, *f3; FLOAT (*u0ij)(); INT u, i;
{ eprintf("Error: nc_H10ij_2 not available.\n"); return(0.);  }

FLOAT nc_lin_div(tGrid,pelement,u)
GRID *tGrid; ELEMENT *pelement; INT u;
{ eprintf("Error: nc_lin_div not available.\n"); return(0.);  }

#endif

/* (H10-error on one element)^2 */
FLOAT nc_H10(tGrid,pelement,u)
GRID *tGrid;
ELEMENT *pelement;
INT u;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM], detB, b[DIM2][DIM2];
   FACE *f1, *f2, *f3;
   
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   detB = barycentric_coordinates(x3,x1,x2,b);
   
   return(
  (nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][0],b[2][0],b[0][0],u011,u,0) +
   nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][1],b[2][1],b[0][1],u012,u,0) +
   nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][0],b[2][0],b[0][0],u021,u,1) +
   nc_H10ij_2(x1,x2,x3,x12,x13,x23,f1,f2,f3,b[1][1],b[2][1],b[0][1],u022,u,1))*detB);
}

#if (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (ELEMENT_TYPE == SIMPLEX)

FLOAT P1disc_sL2s(tGrid,pel,u0,u)  /* L2 norm of the difference between u0    */
GRID *tGrid;      /* and a discontinuous scalar pw. lin. function saved in u  */
ELEMENT *pel;     /* VALUES CORRESPOND TO MIDPOINTS OF EDGES                  */
FLOAT (*u0)();
INT u;
{
   FLOAT *x1,*x2,*x3,x12[DIM],x13[DIM],x23[DIM],p1,p2,p3,p12,p13,p23;
   FACE *f1, *f2, *f3;
                     
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pel,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pel);
   MIDPOINTS(x1,x2,x3,x12,x13,x23)
   L2_lin_face_values(x1,x2,x3,x12,x13,x23,
                      &p1,&p2,&p3,&p12,&p13,&p23,
                      EDSN(pel,u,2),EDSN(pel,u,1),EDSN(pel,u,0),u0);
   return( p22(p1,p2,p3,p12,p13,p23)*volume(x1,x2,x3) );
}

#else

FLOAT P1disc_sL2s(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{ eprintf("Error: P1disc_sL2s not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

/* (H10-error on one element)^2 */
FLOAT H10q(tGrid,pelement,u011,u012,u021,u022,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u011)(), (*u012)(), (*u021)(), (*u022)();
INT u;
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], 
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], 
         x221[DIM], x223[DIM], x331[DIM], x332[DIM], detB, b[DIM2][DIM2];
 
   LTOPNODES_OF_ELEMENT(n1,n2,n3,pelement,tGrid);
   LTOPFACES_OF_ELEMENT(f1,f2,f3,pelement,tGrid);
   VERTICES_OF_ELEMENT(x1,x2,x3,pelement);
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   detB = barycentric_coordinates(x3,x1,x2,b);
   return((sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
           ND(n1,u,0),ND(n2,u,0),ND(n3,u,0),FDV(f3,u,0),FDV(f2,u,0),FDV(f1,u,0),
           b[1][0],b[2][0],b[0][0],u011,u) +
           sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
           ND(n1,u,0),ND(n2,u,0),ND(n3,u,0),FDV(f3,u,0),FDV(f2,u,0),FDV(f1,u,0),
           b[1][1],b[2][1],b[0][1],u012,u) +
           sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
           ND(n1,u,1),ND(n2,u,1),ND(n3,u,1),FDV(f3,u,1),FDV(f2,u,1),FDV(f1,u,1),
           b[1][0],b[2][0],b[0][0],u021,u) +
           sH10qj(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332,
           ND(n1,u,1),ND(n2,u,1),ND(n3,u,1),FDV(f3,u,1),FDV(f2,u,1),FDV(f1,u,1),
           b[1][1],b[2][1],b[0][1],u022,u))*detB);
}

#endif  /*  (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)  */
     
#endif  /* DIM == 2 */

#if (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

FLOAT sL2_q1rot(tGrid,pelement,u0,u)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u0)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], s=0, z;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   LTOPFACES_OF_4ELEMENT(fa0,fa1,fa2,fa3,pelement,tGrid);
   Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac);
   for (k=0; k < QS5_N; k++){
      z = REF_Q1ROT(QS5_P[k],FD(fa0,u),FD(fa1,u),FD(fa2,u),FD(fa3,u)) - 
          fcn_bilin_value(QS5_P[k],a,c,alpha,u0);
      s += QS5_W[k]*z*z*fabs(LINV(jac,QS5_P[k]));
   }
   return(s*QR_VOL);
}

FLOAT sH10_q1rot(tGrid,pelement,u01,u02,u03,u)   /*  u03 not used  */
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   FLOAT v0, v1, v2, v3, b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], 
         jac[DIM2], x[DIM], z, z0, z1, s=0;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   LTOPFACES_OF_4ELEMENT(fa0,fa1,fa2,fa3,pelement,tGrid);
   Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac);
   v0 = FD(fa0,u); 
   v1 = FD(fa1,u); 
   v2 = FD(fa2,u); 
   v3 = FD(fa3,u); 
   for (k=0; k < QS5_N; k++){
      SET1(x,QS5_P[k])
      z = LINV(jac,x);
      z0 = (DX_REF_Q1ROT(x,v0,v1,v2,v3)*LINV(b[0][0],x) +
            DY_REF_Q1ROT(x,v0,v1,v2,v3)*LINV(b[1][0],x))/z
            - fcn_bilin_value(x,a,c,alpha,u01);
      z1 = (DX_REF_Q1ROT(x,v0,v1,v2,v3)*LINV(b[0][1],x) +
            DY_REF_Q1ROT(x,v0,v1,v2,v3)*LINV(b[1][1],x))/z
            - fcn_bilin_value(x,a,c,alpha,u02);
      s += QS5_W[k]*(z0*z0+z1*z1)*fabs(z);
   }
   return(s*QR_VOL);
}

#else

FLOAT sL2_q1rot(tGrid,pelement,u0,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u;
{  eprintf("Error: sL2_q1rot not available.\n");  }

FLOAT sH10_q1rot(tGrid,pelement,u01,u02,u03,u)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u01)(), (*u02)(), (*u03)(); INT u;
{  eprintf("Error: sH10_q1rot not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

FLOAT L2s_q1rot(tGrid,pelement,u0,u,j)/* L2 norm of the j-th component of the */
GRID *tGrid;               /* difference between u0 and a function saved in u */
ELEMENT *pelement;
FLOAT (*u0)();
INT u, j;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], s=0, z;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   LTOPFACES_OF_4ELEMENT(fa0,fa1,fa2,fa3,pelement,tGrid);
   Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac);
   for (k=0; k < QS5_N; k++){
      z = REF_Q1ROT(QS5_P[k],FDV(fa0,u,j),FDV(fa1,u,j),
                             FDV(fa2,u,j),FDV(fa3,u,j)) - 
          fcn_bilin_value(QS5_P[k],a,c,alpha,u0);
      s += QS5_W[k]*z*z*fabs(LINV(jac,QS5_P[k]));
   }
   return(s*QR_VOL);
}

void q5_values(h,xc,sphi,u0,f0,f1,f2,f3,x11,x12,x13,x21,x22,x23,x31,x32,x33,
                                        v11,v12,v13,v21,v22,v23,v31,v32,v33)
DOUBLE h, *xc, (*sphi)(), (*u0)(), f0, f1, f2, f3,
       *x11, *x12, *x13, *x21, *x22, *x23, *x31, *x32, *x33,
       *v11, *v12, *v13, *v21, *v22, *v23, *v31, *v32, *v33;
{
   *v11 = sphi(h,xc,x11,f0,f1,f2,f3) - u0(x11);
   *v12 = sphi(h,xc,x12,f0,f1,f2,f3) - u0(x12);
   *v13 = sphi(h,xc,x13,f0,f1,f2,f3) - u0(x13);
   *v21 = sphi(h,xc,x21,f0,f1,f2,f3) - u0(x21);
   *v22 = sphi(h,xc,x22,f0,f1,f2,f3) - u0(x22);
   *v23 = sphi(h,xc,x23,f0,f1,f2,f3) - u0(x23);
   *v31 = sphi(h,xc,x31,f0,f1,f2,f3) - u0(x31);
   *v32 = sphi(h,xc,x32,f0,f1,f2,f3) - u0(x32);
   *v33 = sphi(h,xc,x33,f0,f1,f2,f3) - u0(x33);
}

/*  integration formula exact for Q5 functions  */
FLOAT q1rot_H10_ar(tGrid,pelem,u,u00,u01,u10,u11)
GRID *tGrid;
ELEMENT *pelem;
DOUBLE (*u00)(), (*u01)(), (*u10)(), (*u11)();
INT u;
{
   FACE *fa0, *fa1, *fa2, *fa3;
   DOUBLE xc[DIM], x11[DIM], x12[DIM], x13[DIM], x21[DIM], x22[DIM], x23[DIM], 
         x31[DIM], x32[DIM], x33[DIM], 
         v11, v12, v13, v21, v22, v23, v31, v32, v33, h, e00, e01, e10, e11;

   FACES_OF_4ELEMENT(fa0,fa1,fa2,fa3,pelem);
   xc[0] = 0.5*(pelem->n[0]->myvertex->x[0] + pelem->n[1]->myvertex->x[0]);
   xc[1] = 0.5*(pelem->n[1]->myvertex->x[1] + pelem->n[2]->myvertex->x[1]);
   h     =      pelem->n[1]->myvertex->x[0] - pelem->n[0]->myvertex->x[0];
   q5_points(h,xc,x11,x12,x13,x21,x22,x23,x31,x32,x33);
   q5_values(h,xc,sphi_0,u00,
             FDV(fa0,u,0),FDV(fa1,u,0),FDV(fa2,u,0),FDV(fa3,u,0),
             x11, x12, x13, x21, x22, x23, x31, x32, x33,
            &v11,&v12,&v13,&v21,&v22,&v23,&v31,&v32,&v33);
   e00 = ( 25.*(v11*v11 + v13*v13 + v31*v31 + v33*v33)+
           40.*(v12*v12 + v21*v21 + v23*v23 + v32*v32)+ 64.*v22*v22 );
   q5_values(h,xc,sphi_1,u01,
             FDV(fa0,u,0),FDV(fa1,u,0),FDV(fa2,u,0),FDV(fa3,u,0),
             x11, x12, x13, x21, x22, x23, x31, x32, x33,
            &v11,&v12,&v13,&v21,&v22,&v23,&v31,&v32,&v33);
   e01 = ( 25.*(v11*v11 + v13*v13 + v31*v31 + v33*v33)+
           40.*(v12*v12 + v21*v21 + v23*v23 + v32*v32)+ 64.*v22*v22 );
   q5_values(h,xc,sphi_0,u10,
             FDV(fa0,u,1),FDV(fa1,u,1),FDV(fa2,u,1),FDV(fa3,u,1),
             x11, x12, x13, x21, x22, x23, x31, x32, x33,
            &v11,&v12,&v13,&v21,&v22,&v23,&v31,&v32,&v33);
   e10 = ( 25.*(v11*v11 + v13*v13 + v31*v31 + v33*v33)+
           40.*(v12*v12 + v21*v21 + v23*v23 + v32*v32)+ 64.*v22*v22 );
   q5_values(h,xc,sphi_1,u11,
             FDV(fa0,u,1),FDV(fa1,u,1),FDV(fa2,u,1),FDV(fa3,u,1),
             x11, x12, x13, x21, x22, x23, x31, x32, x33,
            &v11,&v12,&v13,&v21,&v22,&v23,&v31,&v32,&v33);
   e11 = ( 25.*(v11*v11 + v13*v13 + v31*v31 + v33*v33)+
           40.*(v12*v12 + v21*v21 + v23*v23 + v32*v32)+ 64.*v22*v22 );
   return((e00 + e01 + e10 + e11)/324.*h*h);
}

FLOAT q1rot_H10(tGrid,pelement,u,u00,u01,u10,u11)
GRID *tGrid;
ELEMENT *pelement;
FLOAT (*u00)(), (*u01)(), (*u10)(), (*u11)();
INT u;
{
   INT k;
   NODE *n0, *n1, *n2, *n3;
   FACE *fa0, *fa1, *fa2, *fa3;
   FLOAT v0_0, v0_1, v0_2, v0_3, v1_0, v1_1, v1_2, v1_3, 
         b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2], x[DIM], 
         z, z0_0, z0_1, z1_0, z1_1, s=0;

   NODES_OF_4ELEMENT(n0,n1,n2,n3,pelement);
   LTOPFACES_OF_4ELEMENT(fa0,fa1,fa2,fa3,pelement,tGrid);
   Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac);
   v0_0 = FDV(fa0,u,0); 
   v0_1 = FDV(fa1,u,0); 
   v0_2 = FDV(fa2,u,0); 
   v0_3 = FDV(fa3,u,0); 
   v1_0 = FDV(fa0,u,1); 
   v1_1 = FDV(fa1,u,1); 
   v1_2 = FDV(fa2,u,1); 
   v1_3 = FDV(fa3,u,1); 
   for (k=0; k < QS5_N; k++){
      SET1(x,QS5_P[k])
      z = LINV(jac,x);
      z0_0 = (DX_REF_Q1ROT(x,v0_0,v0_1,v0_2,v0_3)*LINV(b[0][0],x) +
              DY_REF_Q1ROT(x,v0_0,v0_1,v0_2,v0_3)*LINV(b[1][0],x))/z
              - fcn_bilin_value(x,a,c,alpha,u00);
      z0_1 = (DX_REF_Q1ROT(x,v0_0,v0_1,v0_2,v0_3)*LINV(b[0][1],x) +
              DY_REF_Q1ROT(x,v0_0,v0_1,v0_2,v0_3)*LINV(b[1][1],x))/z
              - fcn_bilin_value(x,a,c,alpha,u01);
      z1_0 = (DX_REF_Q1ROT(x,v1_0,v1_1,v1_2,v1_3)*LINV(b[0][0],x) +
              DY_REF_Q1ROT(x,v1_0,v1_1,v1_2,v1_3)*LINV(b[1][0],x))/z
              - fcn_bilin_value(x,a,c,alpha,u10);
      z1_1 = (DX_REF_Q1ROT(x,v1_0,v1_1,v1_2,v1_3)*LINV(b[0][1],x) +
              DY_REF_Q1ROT(x,v1_0,v1_1,v1_2,v1_3)*LINV(b[1][1],x))/z
              - fcn_bilin_value(x,a,c,alpha,u11);
      s += QS5_W[k]*(z0_0*z0_0 + z0_1*z0_1 + z1_0*z1_0 + z1_1*z1_1)*fabs(z);
   }
   return(s*QR_VOL);
}

void q1rot_H10test(tGrid,u,err)
GRID *tGrid;
FLOAT *err;
INT u;
{
   ELEMENT *pel;
   FLOAT sum=0.;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      sum += q1rot_H10(tGrid,pel,u,u011,u012,u021,u022);
   printf("H10 error: %12.6e\n",(*err=sqrt(sum)));
}

void L_infinity_error_q1rot(tGrid,u0,u1,u,err,vset_edge_value)
GRID *tGrid;
INT u;
FLOAT (*u0)(), (*u1)(), *err;
void (*vset_edge_value)();
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3, *f4;
   FLOAT *x1, *x2, *x3, *x4, a[DIM], d[DIM], max=0.;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
         FACES_OF_4ELEMENT(f1,f2,f3,f4,pel);
         VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);

         vset_edge_value(a,x1,x2,u01,u02);
         SET25(d,FDVP(f1,u),a)
         max = MAX(max,MAX(d[0],d[1]));

         vset_edge_value(a,x2,x3,u01,u02);
         SET25(d,FDVP(f2,u),a)
         max = MAX(max,MAX(d[0],d[1]));

         vset_edge_value(a,x3,x4,u01,u02);
         SET25(d,FDVP(f3,u),a)
         max = MAX(max,MAX(d[0],d[1]));

         vset_edge_value(a,x4,x1,u01,u02);
         SET25(d,FDVP(f4,u),a)
         max = MAX(max,MAX(d[0],d[1]));
      }
   printf("L-infinity error: %12.6e\n",(*err=max));
}

#else

FLOAT L2s_q1rot(tGrid,pelement,u0,u,j)
GRID *tGrid; ELEMENT *pelement; FLOAT (*u0)(); INT u, j;
{  eprintf("Error: L2s_q1rot not available.\n");  }

void L_infinity_error_q1rot(tGrid,u0,u1,u,err,vset_edge_value)
GRID *tGrid; INT u; FLOAT (*u0)(), (*u1)(), *err; void (*vset_edge_value)();
{  eprintf("Error: L_infinity_error_q1rot not available.\n");  }

void q1rot_H10test(tGrid,u,err)
GRID *tGrid; FLOAT *err; INT u;
{  eprintf("Error: q1rot_H10test not available.\n");  }

#endif

#if N_DATA & SCALAR_NODE_DATA

void p1c_L_inf_error(tGrid,u0,u,err)
GRID *tGrid;
FLOAT (*u0)(), *err;
INT u;
{
   GRID *theGrid;
   NODE *theNode;
   FLOAT d, max=0.;

   for (theNode=FIRSTN(tGrid); theNode != NULL; theNode=SUCC(theNode)){
      d = fabs(NDS(theNode,u) - (*u0)(MYVERTEX(theNode)->x));
      if (d > max) max = d;
   }

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRSTN(theGrid); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid)){
            d = fabs(NDS(theNode,u) - (*u0)(MYVERTEX(theNode)->x));
            if (d > max) max = d;
         }
   printf("L-infinity error: %12.6e\n",(*err=max));
}

#else

void p1c_L_inf_error(tGrid,u0,u,err)
GRID *tGrid; FLOAT (*u0)(), *err; INT u;
{  eprintf("Error: p1c_L_inf_error not available.\n");  }

#endif

FLOAT H10(tGrid,pelement,u,H10ij,degree)
GRID *tGrid;
ELEMENT *pelement;
INT u, degree;
FLOAT (*H10ij)();
{
   if (degree == 2) return(H10_2(tGrid,pelement,u,H10ij));
   else if (degree == 3) return(H10_3(tGrid,pelement,u,H10ij));
   else{
      eprintf("Error: H10 not available.\n");
      return(0.);
   }
}

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void general_L_inf_error(tGrid,u0,u,err,rm_type,finite_el)
GRID *tGrid;
FLOAT *err, (*u0)();
INT u, rm_type;
FINITE_ELEMENT finite_el;
{
   ELEMENT *pel;
   FLOAT d, max=0.;
   REF_MAPPING ref_map;
   DOUBLE_FUNC *nd, *fd, *ed;
   INT i, k, m, kk, nn, mm;

   nd = finite_el.ndofs;
   fd = finite_el.fdofs;
   ed = finite_el.edofs;
   kk = finite_el.k;
   nn = finite_el.n;
   mm = finite_el.m;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      reference_mapping(pel,rm_type,&ref_map);
      if (kk)
         for (i=0; i < NVERT; i++){
            m = i*kk;
            for (k = m; k < m+kk; k++){
               d = fabs(NDMV(pel->n[i],u,k-m) - (nd[k])(u0,&ref_map));
               if (d > max) max = d;
            }
         }
      if (nn)
         for (i=0; i < NVERT; i++){
            m = i*nn;
            for (k = m; k < m+nn; k++){
               d = fabs(FDMV(pel->f[i],u,k-m) - (fd[k])(u0,&ref_map));
               if (d > max) max = d;
            }
         }
      for (i = 0; i < mm; i++){
         d = fabs(EDMV(pel,u,i) - (ed[i])(u0,&ref_map));
         if (d > max) max = d;
      }
   }
   printf("L-infinity error: %12.6e\n",(*err=max));
}

#else

void general_L_inf_error(tGrid,u0,u,err,rm_type,finite_el)
GRID *tGrid; FLOAT *err, (*u0)(); INT u, rm_type; FINITE_ELEMENT finite_el;
{  eprintf("Error: general_L_inf_error not available.\n");  }

#endif

DOUBLE gs_max_error_on_elem(tGrid,pel,u,u0,points,n,finite_el,rmtype)
GRID *tGrid;
ELEMENT *pel;
INT u, n, rmtype;
DOUBLE (*u0)(), points[][DIM];
FINITE_ELEMENT finite_el;
{
   INT k;
   DOUBLE err=0., s;
   REF_MAPPING ref_map;

   reference_mapping(pel,rmtype,&ref_map);
   for (k = 0; k < n; k++){
      s = fabs(gfunc_value_ref(tGrid,pel,u,points[k],finite_el) -
               fcn_ref_map_value(points[k],&ref_map,u0));
      err = MAX(err,s);
   }
   return(err);
}

void gs_max_error(tGrid,u,u0,err,points,n,finite_el,rmtype)
GRID *tGrid;
INT u, n, rmtype;
DOUBLE (*u0)(), points[][DIM], *err;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   DOUBLE s, max=0.;

   *err = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      s = gs_max_error_on_elem(tGrid,pel,u,u0,points,n,finite_el,rmtype);
      *err = MAX(*err,s);
   }
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            s = gs_max_error_on_elem(tGrid,pel,u,u0,points,n,finite_el,rmtype);
            *err = MAX(*err,s);
         }
   printf("L-infinity error: %12.6e\n",*err);
}

void gs_max_error_wb(tGrid,u,u0,err,points,n,finite_el,rmtype)
GRID *tGrid;
INT u, n, rmtype;
DOUBLE (*u0)(), points[][DIM], *err;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   DOUBLE s, max=0.;
   INT i, j;

   *err = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = 1;
      for (i = 0; i < NVERT; i++)
         if (pel->n[i]->myvertex->x[1] > 0.99999)
            j = 0;
      if (j){
         s = gs_max_error_on_elem(tGrid,pel,u,u0,points,n,finite_el,rmtype);
         *err = MAX(*err,s);
      }
   }
   printf("L-infinity error: %12.6e\n",*err);
}

/******************************************************************************/

#if DIM == 3
#define SUM_L2S       sum  += L2s(tGrid,pel,u01,u,0) + L2s(tGrid,pel,u02,u,1) +\
                                                       L2s(tGrid,pel,u03,u,2);
#define SUM_LINL2S    sum1 += linL2s(tGrid,pel,u01,u,0) +                      \
                              linL2s(tGrid,pel,u02,u,1) +                      \
                              linL2s(tGrid,pel,u03,u,2);
#define SUM_BUBBLEL2S sum2 += bubbleL2s(tGrid,pel,u01,ub,0) +                  \
                              bubbleL2s(tGrid,pel,u02,ub,1) +                  \
                              bubbleL2s(tGrid,pel,u03,ub,2);
#define SUM_LINL2N    sum  += linL2norms(tGrid,pel,u,0) +                      \
                              linL2norms(tGrid,pel,u,1) +                      \
                              linL2norms(tGrid,pel,u,2);
#define SUM_PL2SLIN   sum += pL2s_lin(pel,p0,ED(pel,eu)+diff,ND(pel->n[0],p,0),\
                         ND(pel->n[1],p,0),ND(pel->n[2],p,0),ND(pel->n[3],p,0));
#define SUM_PL2SLS    sum += pL2s_lin(pel,p0,diff,NDS(pel->n[0],p),            \
                         NDS(pel->n[1],p),NDS(pel->n[2],p),NDS(pel->n[3],p));
#define SUM_PL2SLN    sum += ((NDS(pel->n[0],p) + NDS(pel->n[1],p) +           \
                               NDS(pel->n[2],p) + NDS(pel->n[3],p) - 4.*diff)* \
                              (NDS(pel->n[0],p) + NDS(pel->n[1],p) +           \
                               NDS(pel->n[2],p) + NDS(pel->n[3],p) - 4.*diff)+ \
                              (NDS(pel->n[0],p)-diff)*(NDS(pel->n[0],p)-diff)+ \
                              (NDS(pel->n[1],p)-diff)*(NDS(pel->n[1],p)-diff)+ \
                              (NDS(pel->n[2],p)-diff)*(NDS(pel->n[2],p)-diff)+ \
                              (NDS(pel->n[3],p)-diff)*(NDS(pel->n[3],p)-diff)) \
                              *VOLUME(pel)/20.;
/*  #define SUM_H10Q      sum += H10q(tGrid,pel,u011,u012,u013,u021,u022,u023,u031,u032,u033,u);  */
#define SUM_H10Q       
#define SUM_NC_SL2S   sum1 += e0 = -1.;
#define SUM_NC_SH10   sum2 += e1 = -1.;
#define SUM_NC_L2CONV sum3 += ec = -1.;
#define SUM_NC_L2S    sum1 += -1.;
#define SUM_P1DISC_SL2S sum1 += -1.;
#else
#define SUM_L2S       sum  += L2s(tGrid,pel,u01,u,0) + L2s(tGrid,pel,u02,u,1);
#define SUM_LINL2S    sum1 += linL2s(tGrid,pel,u01,u,0) +                      \
                              linL2s(tGrid,pel,u02,u,1);
#define SUM_BUBBLEL2S sum2 += bubbleL2s(tGrid,pel,u01,ub,0) +                  \
                              bubbleL2s(tGrid,pel,u02,ub,1);
#define SUM_LINL2N    sum  += linL2norms(tGrid,pel,u,0) +                      \
                              linL2norms(tGrid,pel,u,1);
#define SUM_PL2SLIN   sum += pL2s_lin(pel,p0,ED(pel,eu)+diff,ND(pel->n[0],p,0),\
                                           ND(pel->n[1],p,0),ND(pel->n[2],p,0));
#define SUM_PL2SLS    sum += pL2s_lin(pel,p0,diff,NDS(pel->n[0],p),            \
                                           NDS(pel->n[1],p),NDS(pel->n[2],p));
#define SUM_PL2SLN    sum += ((NDS(pel->n[0],p) + NDS(pel->n[1],p) - 2.*diff)* \
                              (NDS(pel->n[0],p) + NDS(pel->n[1],p) - 2.*diff)+ \
                              (NDS(pel->n[0],p) + NDS(pel->n[2],p) - 2.*diff)* \
                              (NDS(pel->n[0],p) + NDS(pel->n[2],p) - 2.*diff)+ \
                              (NDS(pel->n[1],p) + NDS(pel->n[2],p) - 2.*diff)* \
                              (NDS(pel->n[1],p) + NDS(pel->n[2],p) - 2.*diff)) \
                              *VOLUME(pel)/12.;
#define SUM_H10Q      sum += H10q(tGrid,pel,u011,u012,u021,u022,u);
#define SUM_NC_SL2S   sum1 += e0 = nc_sL2s(tGrid,pel,u0,u);
#define SUM_NC_SH10   sum2 += e1 = nc_sH10(tGrid,pel,u01,u02,u);
#define SUM_NC_L2CONV sum3 += ec = nc_L2conv(tGrid,pel,bb0,bb1,u01,u02,u);
#define SUM_NC_L2S    sum1 += nc_L2s(tGrid,pel,u01,u,0) +                      \
                              nc_L2s(tGrid,pel,u02,u,1);
#define SUM_P1DISC_SL2S sum1 += P1disc_sL2s(tGrid,pel,u0,u);
#endif

void L2test(tGrid,i,u,err,L2s)
GRID *tGrid;
INT i, u;
FLOAT *err, (*L2s)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_L2S

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_L2S
   *err=sqrt(sum);
   if (i) printf("L2-error: %12.6e\n",*err);
}

void L2errors(tGrid,i,u,ub,err,errlin,errbubble,L2s,linL2s,bubbleL2s)
GRID *tGrid;
INT i, u, ub;
FLOAT *err, *errlin, *errbubble, (*L2s)(), (*linL2s)(), (*bubbleL2s)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum1=0.0, sum2=0.0;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
     SUM_L2S
     SUM_LINL2S
     SUM_BUBBLEL2S
   }
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            SUM_L2S
            SUM_LINL2S
            SUM_BUBBLEL2S
         }
   *err=sqrt(sum);
   *errlin=sqrt(sum1);
   *errbubble=sqrt(sum2);
   if (i) printf("L2-error: %12.6e (lin: %12.6e, bubble: %12.6e)\n",
                 (*err=sqrt(sum)),(*errlin=sqrt(sum1)),(*errbubble=sqrt(sum2)));
}

void linL2norm(tGrid,i,u,err)
GRID *tGrid;
INT i, u;
FLOAT *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
     SUM_LINL2N
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_LINL2N
   *err=sqrt(sum);
   if (i) printf("L2 norm: %12.6e\n",*err);
}

/* L2 norm of the divergence of the linear part of velocity */
void lin_div_norm(tGrid,u,err)
GRID *tGrid;
INT u;
FLOAT *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
     sum += lin_div(tGrid,pel,u);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += lin_div(tGrid,pel,u);
   *err=sqrt(sum);
   printf("|| div u_lin ||: %12.6e\n",*err);
}
 
void H10err(tGrid,u,err,H10ij,degree)
GRID *tGrid;
INT u, degree;
FLOAT *err, (*H10ij)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
                    
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum  += H10(tGrid,pel,u,H10ij,degree);

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum  += H10(tGrid,pel,u,H10ij,degree);

   printf("H10-error: %12.6e\n",(*err=sqrt(sum)));
}

void H10err_p1c_new_fbub(tGrid,u,err)
GRID *tGrid;
INT u;
FLOAT *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
                    
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum  += H10s(tGrid,pel,u);

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum  += H10s(tGrid,pel,u);

   printf("H10-error: %12.6e\n",(*err=sqrt(sum)));
}

void H10errors(tGrid,u,err,errlin,errbubble,H10ij,H10ijlin,bubbleH10ij,degree)
GRID *tGrid;
INT u, degree;
FLOAT *err, *errlin, *errbubble, (*H10ij)(), (*H10ijlin)(), (*bubbleH10ij)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum1=0.0, sum2=0.0;
                    
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      sum  += H10(tGrid,pel,u,H10ij,degree);
      sum1 += H10(tGrid,pel,u,H10ijlin,degree);
      sum2 += H10(tGrid,pel,u,bubbleH10ij,degree);
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            sum  += H10(tGrid,pel,u,H10ij,degree);
            sum1 += H10(tGrid,pel,u,H10ijlin,degree);
            sum2 += H10(tGrid,pel,u,bubbleH10ij,degree);
         } 
   printf("H10-error: %12.6e (lin: %12.6e, bubble: %12.6e)\n",
                 (*err=sqrt(sum)),(*errlin=sqrt(sum1)),(*errbubble=sqrt(sum2)));
}

void H10errors_p1c_new_fbub(tGrid,u,err,errlin,errbubble,H10ijlin,degree)
GRID *tGrid;
INT u, degree;
FLOAT *err, *errlin, *errbubble, (*H10ijlin)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum1=0.0, sum2=0.0;
                    
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      sum  += H10s(tGrid,pel,u);
      sum1 += H10(tGrid,pel,u,H10ijlin,degree);
      sum2 += bubbleH10(tGrid,pel,u);
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            sum  += H10s(tGrid,pel,u);
            sum1 += H10(tGrid,pel,u,H10ijlin,degree);
            sum2 += bubbleH10(tGrid,pel,u);
         } 
   printf("H10-error: %12.6e (lin: %12.6e, bubble: %12.6e)\n",
                 (*err=sqrt(sum)),(*errlin=sqrt(sum1)),(*errbubble=sqrt(sum2)));
}

void linH10norm(tGrid,u,err,space)
GRID *tGrid;
INT u, space;
FLOAT *err;
{
   printf("The following is a H10 norm and not an error\n");
   if ((space == P1C_ELBUB) || (DIM == 3))
      H10err(tGrid,u,err,linH10ijnorm_3,3);
   else
      H10err(tGrid,u,err,linH10ijnorm_2,2);
}

void H10qtest_iso(tGrid,u,err)
GRID *tGrid;
FLOAT  *err;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum += H10q_iso(tGrid,pel,u011,u012,u021,u022,u);
 
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += H10q_iso(tGrid,pel,u011,u012,u021,u022,u);
   printf("H10 error: %12.6e\n",(*err=sqrt(sum)));
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void H10qtest(tGrid,u,err)
GRID *tGrid;
FLOAT  *err;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_H10Q
  
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_H10Q
   printf("H10 error: %12.6e\n",(*err=sqrt(sum)));
}

#else

void H10qtest(tGrid,u,err)
GRID *tGrid; FLOAT  *err; INT u;
{  eprintf("Error: H10qtest not available.\n");  }

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void pL2norm(tGrid,u,err,t)
GRID *tGrid;
INT u, t;
FLOAT *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum2=0.0, vol;
  
   if (t & WITHOUT_FIRST)
      ED(FIRSTELEMENT(tGrid),u) = 0.0;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){ 
      vol   = VOLUME(pel);
      sum2 += ED(pel,u)*vol;
      sum  += vol;
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            vol   = VOLUME(pel);
            sum2 += ED(pel,u)*vol;
            sum  += vol;
         }
  
   sum2 /= sum;
   sum = 0.0;  
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      sum += (ED(pel,u)-sum2)*(ED(pel,u)-sum2)*VOLUME(pel);
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += (ED(pel,u)-sum2)*(ED(pel,u)-sum2)*VOLUME(pel);
   printf("L2 norm of pressure: %12.6e\n",(*err=sqrt(sum)));
}

FLOAT P0_average(pel,p,u)
ELEMENT *pel;
INT p, u;
{
    return(ED(pel,p));
}

#else  /*  !(E_DATA & SCALAR_ELEMENT_DATA)  */

void pL2norm(tGrid,u,err,t)
GRID *tGrid; INT u, t; FLOAT *err;
{  eprintf("Error: pL2norm not available.\n");  }

FLOAT P0_average(pel,p,u)
ELEMENT *pel; INT p, u;
{  eprintf("Error: P0_average not available.\n"); return(0.);  }

#endif

#if (E_DATA & SCALAR_ELEMENT_DATA) && (N_DATA & VECTOR_NODE_DATA)

FLOAT P1P0_average(pel,p,u)
ELEMENT *pel;
INT p, u;
{
    return(ED(pel,u)+AVER_N(pel,p,0));
}

#else

FLOAT P1P0_average(pel,p,u)
ELEMENT *pel; INT p, u; 
{  eprintf("Error: P1P0_average not available.\n"); return(0.);  }

#endif

#if F_DATA & SCALAR_FACE_DATA

FLOAT P1nc_average(pel,p)
ELEMENT *pel;
INT p;
{
    return(AVER_FS(pel,p));
}

#else

FLOAT P1nc_average(pel,p)
ELEMENT *pel; INT p;
{  eprintf("Error: P1nc_average not available.\n"); return(0.);  }

#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES

FLOAT P1disc_average(pel,p)
ELEMENT *pel;
INT p;
{
    return(AVER_ESN(pel,p));
}

#else

FLOAT P1disc_average(pel,p)
ELEMENT *pel; INT p;
{ eprintf("Error: P1disc_averag not available.\n"); return(0.);  }

#endif

FLOAT mean_values(tGrid,p,u,p0,average,integr)
GRID *tGrid;
FLOAT (*p0)(), (*average)(), (*integr)();
INT p, u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT mean_p=0., mean_p0=0., sum=0., vol;
  
   mean_p = mean_p0 = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){ 
      vol   = VOLUME(pel);
      mean_p  += average(pel,p,u,0)*vol;
      mean_p0 += integr(pel,p0)*vol;
      sum  += vol;
   }
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            vol   = VOLUME(pel);
            mean_p  += average(pel,p,u,0)*vol;
            mean_p0 += integr(pel,p0)*vol;
            sum  += vol;
         }
   return((mean_p0 - mean_p)/sum);
}

#if E_DATA & SCALAR_ELEMENT_DATA

void pL2test(tGrid,p0,u,err,t,average,integr,p_L2s)
GRID *tGrid;
FLOAT *err, (*p0)(), (*p_L2s)(), (*average)(), (*integr)();
INT u, t;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, diff;
  
   if (t & WITHOUT_FIRST)
      ED(FIRSTELEMENT(tGrid),u) = 0.0;
  
   diff = mean_values(tGrid,u,0,p0,average,integr);

   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      sum += p_L2s(pel,p0,ED(pel,u)+diff);
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
             sum += p_L2s(pel,p0,ED(pel,u)+diff);
    
   printf("L2-error in pressure: %12.6e\n",(*err=sqrt(sum)));
}

#else  /*  !(E_DATA & SCALAR_ELEMENT_DATA)  */

void pL2test(tGrid,p0,u,err,t,average,integr,p_L2s)
GRID *tGrid; FLOAT *err, (*p0)(), (*p_L2s)(), (*average)(), (*integr)(); INT u, t;
{  eprintf("Error: pL2test not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void pL2test_lin(tGrid,err,p,eu,t)/* pressure consists of a piecewise constant*/
GRID *tGrid;                      /* part stored in eu and a piecewise linear */
FLOAT *err;                       /* part stored in p                         */
INT p, t;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum1, sum2, diff;
  
   if (t & WITHOUT_FIRST)
      ED(FIRSTELEMENT(tGrid),eu) = 0.0;
   diff = mean_values(tGrid,p,eu,p0,P1P0_average,integr5);
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      SUM_PL2SLIN
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_PL2SLIN
    
   printf("L2-error in pressure: %12.6e\n",(*err=sqrt(sum)));
}

#else

void pL2test_lin(tGrid,err,p,eu,t)
GRID *tGrid; FLOAT *err; INT p, t;
{  eprintf("Error: pL2test_lin not available.\n");  }

#endif

void sL2test(tGrid,u0,u,err,sL2)
GRID *tGrid;
INT u;
FLOAT (*sL2)(), (*u0)(), *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum += sL2(tGrid,pel,u0,u);
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += sL2(tGrid,pel,u0,u);
   printf("L2 error: %12.6e\n",(*err=sqrt(sum)));
}

void sH10test(tGrid,u01,u02,u03,u,err,sH10) /*  if DIM == 2, u03 is not used  */
GRID *tGrid;
FLOAT (*sH10)(), (*u01)(), (*u02)(), (*u03)(), *err;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum2=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum2 += sH10(tGrid,pel,u01,u02,u03,u);
  
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum2 += sH10(tGrid,pel,u01,u02,u03,u);
   printf("H10 error: %12.6e\n",(*err=sqrt(sum2)));
}

void s_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err,sL2conv) 
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(), *err, 
      (*sL2conv)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum3=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum3 += sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum3 += sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);
   printf("conv. term error: %12.6e\n",(*err=sqrt(sum3)));
}

void s_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err4,H10,sL2conv,sL2) 
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), (*u03)(), 
      *err4, (*H10)(), (*sL2conv)(), (*sL2)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT e1, ec, sum2=0., sum3=0., sum4=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      sum2 += e1 = H10(tGrid,pel,u01,u02,u03,u);
      sum3 += ec = sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);
      sum4 += eps*e1 + CC0*sL2(tGrid,pel,u0,u) + sd_delta(pel,eps,bb0,bb1)*ec;
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            sum2 += e1 = H10(tGrid,pel,u01,u02,u03,u);
            sum3 += ec = sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);
            sum4 += eps*e1 + CC0*sL2(tGrid,pel,u0,u) + 
                                                   sd_delta(pel,eps,bb0,bb1)*ec;
         }
   printf("SD error: %12.6e\n",(*err4=sqrt(sum4)));
}

void s_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,
                            err1,err2,err3,err4,err5,area,sL2,H10,sL2conv,sLinf)
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT xm, ym, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), 
      (*u03)(), *err1, *err2, *err3, *err4, *err5, *area,
      (*sL2)(), (*H10)(), (*sL2conv)(), (*sLinf)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   NODE *theNode;
   FLOAT e0, e1, ec, sum1=0., sum2=0., sum3=0., sum4=0., max=0., ar=0., d;
                     
   printf("square for errors bounded by xm = %5.3f, ym = %5.3f\n",xm,ym);
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      if (IS_IN_SQUARE(pel,xm,ym)){
         sum1 += e0 = sL2(tGrid,pel,u0,u);
         sum2 += e1 = H10(tGrid,pel,u01,u02,u03,u);
         sum3 += ec = sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);
         sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;
         ar += VOLUME(pel);
         if ((d=sLinf(tGrid,pel,u0,u)) > max) max = d;
      }
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid) && IS_IN_SQUARE(pel,xm,ym)){
            sum1 += e0 = sL2(tGrid,pel,u0,u);
            sum2 += e1 = H10(tGrid,pel,u01,u02,u03,u);
            sum3 += ec = sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u);
            sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;
            ar += VOLUME(pel);
            if ((d=sLinf(tGrid,pel,u0,u)) > max) max = d;
         }
   }
   printf("L2 error: %12.6e\n",(*err1=sqrt(sum1)));
   printf("H10 error: %12.6e\n",(*err2=sqrt(sum2)));
   printf("conv. term error: %12.6e\n",(*err3=sqrt(sum3)));
   printf("SD error: %12.6e\n",(*err4=sqrt(sum4)));
   printf("L-infinity error: %12.6e\n",(*err5=max));
   printf("Computed area: %12.6e\n",(*area=ar));
}

FLOAT general_sL2(tGrid,pel,u0,u,qr,rm_type,finite_el)
GRID *tGrid;
ELEMENT *pel;
FLOAT (*u0)();
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   INT k;
   FLOAT s=0., z;
   REF_MAPPING ref_map;

   reference_mapping(pel,rm_type,&ref_map);
   for (k=0; k < qr.n; k++){
      z = gfunc_value_ref(tGrid,pel,u,qr.points[k],finite_el) -
          fcn_ref_map_value(qr.points[k],&ref_map,u0);
      s += qr.weights[k]*z*z*fabs(jacobian(qr.points[k],&ref_map));
   }
   return(s*QR_VOL);
}

void general_sL2test(tGrid,u0,u,qr,rm_type,finite_el,err,sL2)
GRID *tGrid;
FLOAT (*sL2)(), (*u0)(), *err;
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      sum += sL2(tGrid,pel,u0,u,qr,rm_type,finite_el);
    
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += sL2(tGrid,pel,u0,u,qr,rm_type,finite_el);
   printf("L2 error: %12.6e\n",(*err=sqrt(sum)));
}

void general_sL2test_without_belems(tGrid,u0,u,qr,rm_type,finite_el,err,sL2)
GRID *tGrid;
FLOAT (*sL2)(), (*u0)(), *err;
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.;
   INT i, j;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = 1;
      for (i = 0; i < NVERT; i++)
         if (pel->n[i]->myvertex->x[1] > 0.99999)
            j = 0;
      if (j)
         sum += sL2(tGrid,pel,u0,u,qr,rm_type,finite_el);
   }
   printf("L2 error: %12.6e\n",(*err=sqrt(sum)));
}

FLOAT general_sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el)
GRID *tGrid;
ELEMENT *pel;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   INT k;
   FLOAT b[DIM][DIM], grad[DIM], z, z0, z1, s=0.;
   REF_MAPPING ref_map;

   reference_mapping_with_inverse(pel,rm_type,&ref_map);
   for (k=0; k < qr.n; k++){
      z = jacobian(qr.points[k],&ref_map);
      inv_of_jac_matr_times_jacobian(qr.points[k],b,&ref_map);
      ggrad_value_ref(tGrid,pel,u,qr.points[k],finite_el,grad);
      z0 = (grad[0]*b[0][0] + grad[1]*b[1][0])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u01);
      z1 = (grad[0]*b[0][1] + grad[1]*b[1][1])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u02);
      s += qr.weights[k]*(z0*z0+z1*z1)*fabs(z);
   }
   return(s*QR_VOL);
}

FLOAT general_sH10_test(tGrid,u01,u02,u03,u,qr,rm_type,finite_el,err,sH10)
GRID *tGrid;
FLOAT (*sH10)(), (*u01)(), (*u02)(), (*u03)(), *err;
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum2=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      sum2 += sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el);
  
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum2 += sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el);
   printf("H10 error: %12.6e\n",(*err=sqrt(sum2)));
}

FLOAT general_sH10_test_without_belems(tGrid,u01,u02,u03,u,qr,rm_type,finite_el,err,sH10)
GRID *tGrid;
FLOAT (*sH10)(), (*u01)(), (*u02)(), (*u03)(), *err;
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum2=0.;
   INT i, j;
                     
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = 1;
      for (i = 0; i < NVERT; i++)
         if (pel->n[i]->myvertex->x[1] > 0.99999)
            j = 0;
      if (j)
         sum2 += sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el);
   }
   printf("H10 error: %12.6e\n",(*err=sqrt(sum2)));
}

/* (error of  bb*grad u  on one element)^2 */
FLOAT general_sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u,qr,rm_type,finite_el)
GRID *tGrid;                                        /*  bb2 and u03 not used  */
ELEMENT *pel;
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   INT k;
   FLOAT b[DIM][DIM], grad[DIM], r, z, z0, z1, s=0.;
   REF_MAPPING ref_map;

   reference_mapping_with_inverse(pel,rm_type,&ref_map);
   for (k=0; k < qr.n; k++){
      z = jacobian(qr.points[k],&ref_map);
      inv_of_jac_matr_times_jacobian(qr.points[k],b,&ref_map);
      ggrad_value_ref(tGrid,pel,u,qr.points[k],finite_el,grad);
      z0 = (grad[0]*b[0][0] + grad[1]*b[1][0])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u01);
      z1 = (grad[0]*b[0][1] + grad[1]*b[1][1])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u02);
      r = z0*fcn_ref_map_value(qr.points[k],&ref_map,bb0) +
          z1*fcn_ref_map_value(qr.points[k],&ref_map,bb1);
      s += qr.weights[k]*r*r*fabs(z);
   }
   return(s*QR_VOL);
}

FLOAT general_sL2conv_test(tGrid,bb0,bb1,bb2,u01,u02,u03,u,
                                               qr,rm_type,finite_el,err,sL2conv)
GRID *tGrid;
FLOAT (*sL2conv)(), (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(), 
      *err;
INT u, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum2=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      sum2 += sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u,qr,rm_type,finite_el);
  
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum2 += sL2conv(tGrid,pel,bb0,bb1,bb2,u01,u02,u03,u,qr,rm_type,finite_el);
   printf("L2conv error: %12.6e\n",(*err=sqrt(sum2)));
}

FLOAT general_s_sd_test(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,space,
                                      qr,rm_type,finite_el,err,sH10,sL2conv,sL2)
GRID *tGrid;
FLOAT (*sH10)(), (*sL2conv)(), (*sL2)(), eps, (*bb0)(), (*bb1)(), (*bb2)(), 
      (*u0)(), (*u01)(), (*u02)(), (*u03)(), *err;
INT u, space, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.;
                     
   if(PW_CONST_PAR == YES){
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
         sum += eps*sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el) +
                gsd_tau(qr.points[0],pel,eps,bb0,bb1,space)*sL2conv(tGrid,
                          pel,bb0,bb1,bb2,u01,u02,u03,u,qr,rm_type,finite_el) +
                CC0*sL2(tGrid,pel,u0,u,qr,rm_type,finite_el);
  
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
         for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
            if (IS_LTOP_ELEMENT(pel,tGrid))
               sum += eps*sH10(tGrid,pel,u01,u02,u03,u,qr,rm_type,finite_el) +
                      gsd_tau(qr.points[0],pel,eps,bb0,bb1,space)*sL2conv(tGrid,
                          pel,bb0,bb1,bb2,u01,u02,u03,u,qr,rm_type,finite_el) +
                      CC0*sL2(tGrid,pel,u0,u,qr,rm_type,finite_el);
      printf("SD error: %12.6e\n",(*err=sqrt(sum)));
   }
   else{
      eprintf("Error: general_s_sd_test not available.\n");
      *err = -1.;
   }
}

FLOAT general_lp_err(tGrid,pel,u01,u02,u03,u,
                     eps,n,bb0,bb1,f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el)
GRID *tGrid;
ELEMENT *pel;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   INT i, j, k;
   DOUBLE_FUNC f[8] = { f2, f3, f4, f5, f6, f7, f8, f9 };
   DOUBLE a[8][N_SGGEM], cq[2][N_SGGEM], cr[2][8], x[2][N_SGGEM], b[DIM][DIM],
          grad[DIM], fi, z, z0, z1, s=0.;
   REF_MAPPING ref_map;

   if (n > 8)
      eprintf("Error: too large n in general_lp_err.\n");
   for (i=0; i < n; i++){
      cq[0][i] = cq[1][i] = 0.;
      for (j=0; j < n; j++)
         a[i][j] = 0.;
   }
   reference_mapping_with_inverse(pel,rm_type,&ref_map);
   for (k=0; k < qr.n; k++){
      z = jacobian(qr.points[k],&ref_map);
      inv_of_jac_matr_times_jacobian(qr.points[k],b,&ref_map);
      ggrad_value_ref(tGrid,pel,u,qr.points[k],finite_el,grad);
      z0 = (grad[0]*b[0][0] + grad[1]*b[1][0])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u01);
      z1 = (grad[0]*b[0][1] + grad[1]*b[1][1])/z
            - fcn_ref_map_value(qr.points[k],&ref_map,u02);
      s += qr.weights[k]*(z0*z0+z1*z1)*fabs(z);
      for (i=0; i < n; i++){
         fi = qr.weights[k]*f[i](qr.points[k])*fabs(z);
         cq[0][i] += fi*z0;
         cq[1][i] += fi*z1;
         for (j=0; j < n; j++)
            a[i][j] += fi*f[j](qr.points[k]);
      }
   }
   for (i=0; i < n; i++){
      cr[0][i] = cq[0][i];
      cr[1][i] = cq[1][i];
   }
   sggem(a,cq,x,n,2);
   for (i=0; i < n; i++)
      s -= x[0][i]*cr[0][i] + x[1][i]*cr[1][i];
   return(lp_delta(pel,eps,bb0,bb1)*s*QR_VOL);
}

FLOAT general_conv_lp_err(tGrid,pel,u01,u02,u03,u,
                     eps,n,bb0,bb1,f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el)
GRID *tGrid;
ELEMENT *pel;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   INT i, j, k;
   DOUBLE_FUNC f[8] = { f2, f3, f4, f5, f6, f7, f8, f9 };
   DOUBLE a[8][N_GGEM], cq[8], cr[8], x[8], b[DIM][DIM], xc[DIM],
          grad[DIM], b0x, b1x, fi, z, qk, s=0.;
   REF_MAPPING ref_map;

   if (n > 8)
      eprintf("Error: too large n in general_conv_lp_err.\n");
   for (i=0; i < n; i++){
      cq[i] = 0.;
      for (j=0; j < n; j++)
         a[i][j] = 0.;
   }
   if (USE_BM == YES){
      coord_of_barycentre(pel,xc);
      b0x = bb0(xc);
      b1x = bb1(xc);
   }
   reference_mapping_with_inverse(pel,rm_type,&ref_map);
   for (k=0; k < qr.n; k++){
      z = jacobian(qr.points[k],&ref_map);
      inv_of_jac_matr_times_jacobian(qr.points[k],b,&ref_map);
      ggrad_value_ref(tGrid,pel,u,qr.points[k],finite_el,grad);
      if (USE_BM == NO){
         b0x = fcn_ref_map_value(qr.points[k],&ref_map,bb0);
         b1x = fcn_ref_map_value(qr.points[k],&ref_map,bb1);
      }
      qk = b0x*((grad[0]*b[0][0] + grad[1]*b[1][0])/z
                 - fcn_ref_map_value(qr.points[k],&ref_map,u01)) +
           b1x*((grad[0]*b[0][1] + grad[1]*b[1][1])/z
                 - fcn_ref_map_value(qr.points[k],&ref_map,u02));
      s += qr.weights[k]*qk*qk*fabs(z);
      for (i=0; i < n; i++){
         fi = qr.weights[k]*f[i](qr.points[k])*fabs(z);
         cq[i] += fi*qk;
         for (j=0; j < n; j++)
            a[i][j] += fi*f[j](qr.points[k]);
      }
   }
   for (i=0; i < n; i++)
      cr[i] = cq[i];
   ggem(a,cq,x,n);
   for (i=0; i < n; i++)
      s -= x[i]*cr[i];
   return(lp_delta(pel,eps,bb0,bb1)*s*QR_VOL);
}

#if DATA_S & N_LINK_TO_ELEMENTS

FLOAT general_macro_conv_lp_err(tGrid,pn,u01,u02,u03,u,
                     eps,n,bb0,bb1,f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el)
GRID *tGrid;
NODE *pn;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   NELINK *pnel;
   ELEMENT *pel;
   INT i, j, k;
   DOUBLE_FUNC f[8] = { f2, f3, f4, f5, f6, f7, f8, f9 };
   DOUBLE a[8][N_GGEM], cq[8], cr[8], x[8], b[DIM][DIM],
          grad[DIM], b0x, b1x, fi, z, qk, s=0.;
   REF_MAPPING ref_map;

   if (n > 8)
      eprintf("Error: too large n in general_macro_conv_lp_err.\n");
   for (i=0; i < n; i++){
      cq[i] = 0.;
      for (j=0; j < n; j++)
         a[i][j] = 0.;
   }
   if (USE_BM == YES){
      b0x = bb0(pn->myvertex->x);
      b1x = bb1(pn->myvertex->x);
   }
   for (pnel = NESTART(pn); pnel; pnel = pnel->next){
      pel = pnel->nbel;
      reference_mapping_with_inverse(pel,rm_type,&ref_map);
      for (k=0; k < qr.n; k++){
         z = jacobian(qr.points[k],&ref_map);
         inv_of_jac_matr_times_jacobian(qr.points[k],b,&ref_map);
         ggrad_value_ref(tGrid,pel,u,qr.points[k],finite_el,grad);
         if (USE_BM == NO){
            b0x = fcn_ref_map_value(qr.points[k],&ref_map,bb0);
            b1x = fcn_ref_map_value(qr.points[k],&ref_map,bb1);
         }
         qk = b0x*((grad[0]*b[0][0] + grad[1]*b[1][0])/z
                    - fcn_ref_map_value(qr.points[k],&ref_map,u01)) +
              b1x*((grad[0]*b[0][1] + grad[1]*b[1][1])/z
                    - fcn_ref_map_value(qr.points[k],&ref_map,u02));
         s += qr.weights[k]*qk*qk*fabs(z);
         for (i=0; i < n; i++){
            fi = qr.weights[k]*f[i](qr.points[k])*fabs(z);
            cq[i] += fi*qk;
            for (j=0; j < n; j++)
               a[i][j] += fi*f[j](qr.points[k]);
         }
      }
   }
   for (i=0; i < n; i++)
      cr[i] = cq[i];
   ggem(a,cq,x,n);
   for (i=0; i < n; i++)
      s -= x[i]*cr[i];
   return(lpm_delta(pn,eps,bb0,bb1)*s*QR_VOL);
}

#else

FLOAT general_macro_conv_lp_err(tGrid,pn,u01,u02,u03,u,eps,n,bb0,bb1,f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el)
GRID *tGrid; NODE *pn; FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)(); INT u, n, rm_type; QUADRATURE_RULE qr; FINITE_ELEMENT finite_el;
{  eprintf("Error: general_macro_conv_lp_err not available.\n");  }

#endif

FLOAT general_lp_test(tGrid,u01,u02,u03,u,eps,n,bb0,bb1,
                      f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el,err,lp_err)
GRID *tGrid;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)(), *err, (*lp_err)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      sum += lp_err(tGrid,pel,u01,u02,u03,u,eps,n,bb0,bb1,
                    f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el);
  
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += lp_err(tGrid,pel,u01,u02,u03,u,eps,n,bb0,bb1,
                          f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el);
   if (sum < 0.) *err=-sqrt(-sum);
   else *err=sqrt(sum);
   printf("LP error: %12.6e\n",*err);
}

FLOAT general_macro_lp_test(tGrid,u01,u02,u03,u,eps,n,bb0,bb1,
                      f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el,err,lp_err)
GRID *tGrid;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)(), *err, (*lp_err)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   NODE *pn;
   FLOAT sum=0.;
                     
   for (pn = FIRSTNODE(tGrid); pn; pn = pn->succ)
      if (is_macro_node(pn))
         sum += lp_err(tGrid,pn,u01,u02,u03,u,eps,n,bb0,bb1,
                       f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el);
   if (sum < 0.) *err=-sqrt(-sum);
   else *err=sqrt(sum);
   printf("LP error: %12.6e\n",*err);
}

FLOAT general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,n,bb0,bb1,
                      f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el,err,lp_err)
GRID *tGrid;
FLOAT (*u01)(), (*u02)(), (*u03)(), eps, (*bb0)(), (*bb1)(), (*f2)(), (*f3)(), 
      (*f4)(), (*f5)(), (*f6)(), (*f7)(), (*f8)(), (*f9)(), *err, (*lp_err)();
INT u, n, rm_type;
QUADRATURE_RULE qr;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.;
   INT i, j;
                     
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      j = 1;
      for (i = 0; i < NVERT; i++)
         if (pel->n[i]->myvertex->x[1] > 0.99999)
            j = 0;
    if (j)
         sum += lp_err(tGrid,pel,u01,u02,u03,u,eps,n,bb0,bb1,
                       f2,f3,f4,f5,f6,f7,f8,f9,qr,rm_type,finite_el);
   }
   if (sum < 0.) *err=-sqrt(-sum);
   else *err=sqrt(sum);
   printf("LP error: %12.6e\n",*err);
}

#if (N_DATA & SCALAR_NODE_DATA) && (ELEMENT_TYPE == SIMPLEX)

void pL2test_cont(tGrid,err,p,p0,t) 
GRID *tGrid;
FLOAT *err, (*p0)();
INT p, t;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum1=0.0, sum2=0.0, diff, vol;
  
   if (t & WITHOUT_FIRST)
      NDS(FIRSTN(tGrid),p) = 0.;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){ 
      vol   = VOLUME(pel);
      sum1 += integr5(pel,p0)*vol;
      sum2 += AVER_NS(pel,p)*vol;
      sum  += vol;
   }
  
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            vol   = VOLUME(pel);
            sum1 += integr5(pel,p0)*vol;
            sum2 += AVER_NS(pel,p)*vol;
            sum  += vol;
         }
    
   diff = (sum1 - sum2)/sum;
   sum = 0.0;  
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      SUM_PL2SLS
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
             SUM_PL2SLS
    
   printf("L2-error in pressure: %12.6e\n",(*err=sqrt(sum)));
}

void pL2norm_cont(tGrid,err,p,t) 
GRID *tGrid;
FLOAT *err;
INT p, t;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0, sum2=0.0, diff, vol;
  
   if (t & WITHOUT_FIRST)
      NDS(FIRSTN(tGrid),p) = 0.;
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){ 
      vol   = VOLUME(pel);
      sum2 += AVER_NS(pel,p)*vol;
      sum  += vol;
   }
  
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            vol   = VOLUME(pel);
            sum2 += AVER_NS(pel,p)*vol;
            sum  += vol;
         }
    
   diff = sum2/sum;
   sum = 0.0;  
  
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      SUM_PL2SLN
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
             SUM_PL2SLN
    
   printf("L2-norm of pressure: %12.6e\n",(*err=sqrt(sum)));
}

#else

void pL2test_cont(tGrid,err,p,p0,t) 
GRID *tGrid; FLOAT *err, (*p0)(); INT p, t;
{  eprintf("Error: pL2test_cont not available.\n");  }

void pL2norm_cont(tGrid,err,p,t) 
GRID *tGrid; FLOAT *err; INT p, t;
{  eprintf("Error: pL2norm_cont not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DIM == 2)

FLOAT max_p2c_error_at_15_points_of_element(pel,u0,u)
ELEMENT *pel;
FLOAT (*u0)();
INT u;
{
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, x11[DIM], x22[DIM], x33[DIM], 
         x12[DIM], x13[DIM], x23[DIM], x112[DIM], x113[DIM], 
         x221[DIM], x223[DIM], x331[DIM], x332[DIM], p1, p2, p3, 
         p11, p22, p33, p12, p13, p23, p112, p113, p221, p223, p331, p332,
         v1, v2, v3, v12, v13, v23;
                     
   NODES_OF_ELEMENT(n1,n2,n3,pel);
   FACES_OF_ELEMENT(f1,f2,f3,pel);
   VERTICES_OF_ELEMENT(x1,x2,x3,pel);
   points4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,x223,x331,x332);
   v1 = NDS(n1,u); 
   v2 = NDS(n2,u); 
   v3 = NDS(n3,u); 
   v12 = 0.5*(v1+v2) + 0.25*FD(f3,u);
   v13 = 0.5*(v1+v3) + 0.25*FD(f2,u);
   v23 = 0.5*(v2+v3) + 0.25*FD(f1,u);
   L2_quad_node_values_4(x1,x2,x3,x11,x22,x33,x12,x13,x23,x112,x113,x221,
                  x223,x331,x332,&p1,&p2,&p3,&p11,&p22,&p33,
                  &p12,&p13,&p23,&p112,&p113,&p221,&p223,&p331,&p332,
                  v1,v2,v3,v12,v13,v23,u0);
   return( abs_max_of_15(p1,p2,p3,p11,p22,p33,p12,p13,p23,
                            p112,p113,p221,p223,p331,p332) );
}

void p2c_L_inf_error(tGrid,u0,u,err)
GRID *tGrid;
FLOAT (*u0)(), *err;
INT u;
{
   ELEMENT *pel;
   NODE *theNode;
   FLOAT d, max=0.;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      d = p2c_Linf(tGrid,pel,u0,u);
      if (d > max) max = d;
   }

   printf("L-infinity error: %12.6e\n",(*err=max));
}

#else

void p2c_L_inf_error(tGrid,u0,u,err)
GRID *tGrid; FLOAT (*u0)(), *err; INT u;
{  eprintf("Error: p2c_L_inf_error not available.\n");  }

#endif

void nc_L2err(tGrid,u0,u,err1) 
GRID *tGrid;
FLOAT (*u0)(), *err1;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT e0, sum1=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_NC_SL2S

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_NC_SL2S
   printf("L2 error: %12.6e\n",(*err1=sqrt(sum1)));
}

void nc_H10err(tGrid,u01,u02,u03,u,err2) 
GRID *tGrid;                     /*  if DIM == 2, u03 is not used    */
FLOAT (*u01)(), (*u02)(), (*u03)(), *err2;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT e1, sum2=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_NC_SH10

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_NC_SH10
   printf("H10 error: %12.6e\n",(*err2=sqrt(sum2)));
}

void nc_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err3) 
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)(), *err3;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT ec, sum3=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_NC_L2CONV

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_NC_L2CONV
   printf("conv. term error: %12.6e\n",(*err3=sqrt(sum3)));
}

void nc_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err4) 
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), (*u03)(), 
      *err4;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT e0, e1, ec, sum1=0., sum2=0., sum3=0., sum4=0.;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      SUM_NC_SL2S
      SUM_NC_SH10
      SUM_NC_L2CONV
      sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;
   }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid)){
            SUM_NC_SL2S
            SUM_NC_SH10
            SUM_NC_L2CONV
            sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;
         }
   printf("SD error: %12.6e\n",(*err4=sqrt(sum4)));
}

#if F_DATA & SCALAR_FACE_DATA

void nc_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,
                                                  err1,err2,err3,err4,err5,area)
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT xm, ym, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), 
      (*u03)(), *err1, *err2, *err3, *err4, *err5, *area;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT e0, e1, ec, sum1=0., sum2=0., sum3=0., sum4=0., max=0., ar=0.,
         *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM], d;
                     
   printf("square for errors bounded by xm = %5.3f, ym = %5.3f\n",xm,ym);
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      if (IS_IN_SQUARE(pel,xm,ym)){
         SUM_NC_SL2S
         SUM_NC_SH10
         SUM_NC_L2CONV
         sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;

         FACES_OF_ELEMENT(f1,f2,f3,pel);
         VERTICES_OF_ELEMENT(x1,x2,x3,pel);
         MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
         d = fabs(FD(f1,u) - (*u0)(xc1));
         if (d > max) max = d;
         d = fabs(FD(f2,u) - (*u0)(xc2));
         if (d > max) max = d;
         d = fabs(FD(f3,u) - (*u0)(xc3));
         if (d > max) max = d;
         ar += VOLUME(pel);
      }

   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid) && IS_IN_SQUARE(pel,xm,ym)){
            SUM_NC_SL2S
            SUM_NC_SH10
            SUM_NC_L2CONV
            sum4 += eps*e1 + CC0*e0 + sd_delta(pel,eps,bb0,bb1)*ec;

            FACES_OF_ELEMENT(f1,f2,f3,pel);
            VERTICES_OF_ELEMENT(x1,x2,x3,pel);
            MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
            d = fabs(FD(f1,u) - (*u0)(xc1));
            if (d > max) max = d;
            d = fabs(FD(f2,u) - (*u0)(xc2));
            if (d > max) max = d;
            d = fabs(FD(f3,u) - (*u0)(xc3));
            if (d > max) max = d;
            ar += VOLUME(pel);
         }
   printf("L2 error: %12.6e\n",(*err1=sqrt(sum1)));
   printf("H10 error: %12.6e\n",(*err2=sqrt(sum2)));
   printf("conv. term error: %12.6e\n",(*err3=sqrt(sum3)));
   printf("SD error: %12.6e\n",(*err4=sqrt(sum4)));
   printf("L-infinity error: %12.6e\n",(*err5=max));
   printf("Computed area: %12.6e\n",(*area=ar));
}

#else

void nc_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err1,err2,err3,err4,err5,area)
GRID *tGrid; FLOAT xm, ym, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), (*u03)(), *err1, *err2, *err3, *err4, *err5, *area; INT u;
{  eprintf("Error: nc_errors_on_square not available.\n");  }

#endif

void nc_L2test(tGrid,u,err)
GRID *tGrid;
INT u;
FLOAT *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum1=0.0;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_NC_L2S
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_NC_L2S
   printf("L2 error: %12.6e\n",(*err=sqrt(sum1)));
}

void nc_H10test(tGrid,u,err)
GRID *tGrid;
FLOAT *err;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT e1, sum2=0.0;
                     
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      sum2 += e1 = nc_H10(tGrid,pel,u);
  
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum2 += e1 = nc_H10(tGrid,pel,u);
   printf("H10 error: %12.6e\n",(*err=sqrt(sum2)));
}

/* L2 norm of the divergence of a nonconforming pw. linear velocity */
void nc_lin_div_norm(tGrid,u,err)
GRID *tGrid;
FLOAT *err;
INT u;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum=0.0;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
     sum += nc_lin_div(tGrid,pel,u);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            sum += nc_lin_div(tGrid,pel,u);
   *err=sqrt(sum);
   printf("|| div u_lin ||: %12.6e\n",*err);
}
 
void P1disc_sL2test(tGrid,u0,u,err)
GRID *tGrid;
INT u;
FLOAT (*u0)(), *err;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT sum1=0.0;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      SUM_P1DISC_SL2S
    
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            SUM_P1DISC_SL2S
   printf("L2 error: %12.6e\n",(*err=sqrt(sum1)));
}

void sL2_error(tGrid,u,v,err,u0,space,structure,t,type)
GRID *tGrid;
INT u, v, space, structure, t, type;
FLOAT *err, (*u0)();
{
   switch(space){
   case P0:      if (structure == P_SCALAR){
                    if (ELEMENT_TYPE==SIMPLEX)
                       pL2test(tGrid,u0,u,err,t,P0_average,integr5,pL2s);
                    else if (ELEMENT_TYPE==CUBE)
                       pL2test(tGrid,u0,u,err,t,P0_average,integr_q5,pL2s_q);
                    else
                       eprintf("Error: sL2_error not available.\n");
                 }
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1_NC:   if (structure == P_SCALAR){
                    if (t & WITHOUT_FIRST)
                       set_first(tGrid,0.,u,t,type);
                    set_value(tGrid,
                       mean_values(tGrid,u,0,u0,P1nc_average,integr5),v,0,type);
                    add(tGrid,u,v,v,0,type);
                    printf("Pressure:   ");
                    nc_L2err(tGrid,u0,v,err);
                 }
                 else if (structure == SCALAR)
                    nc_L2err(tGrid,u0,u,err);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1_MOD:  if (structure == SCALAR){
                    vector_to_scalar_f(tGrid,u,v,0);
                    nc_L2err(tGrid,u0,v,err);
                 }
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1_DISC: if (structure == P_SCALAR){
                    if (t & WITHOUT_FIRST)
                       set_first(tGrid,0.,u,t,type);
                    set_value(tGrid,
                     mean_values(tGrid,u,0,u0,P1disc_average,integr5),v,0,type);
                    add(tGrid,u,v,v,0,type);
                    printf("Pressure:   ");
                    P1disc_sL2test(tGrid,u0,v,err);
                 }
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1C:     if (structure == P_SCALAR)
                    pL2test_cont(tGrid,err,u,u0,t);
                 else if (structure == SCALAR)
                    sL2test(tGrid,u0,u,err,sL2s);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1C_P2L: if (structure == P_SCALAR)
                    pL2test_cont(tGrid,err,u,u0,t);
                 else if (structure == SCALAR)
                    sL2test(tGrid,u0,u,err,sL2s);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P1C_ELBUB: if (structure == SCALAR)
                    eprintf("Error: sL2_error not available.\n");
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case Q1C:     if (structure == P_SCALAR)
                    eprintf("Error: sL2_error not available.\n");
                 else if (structure == SCALAR)
                    sL2test(tGrid,u0,u,err,sL2_q1);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case Q1ROT:   if (structure == SCALAR)
                    sL2test(tGrid,u0,u,err,sL2_q1rot);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P2C:     if (structure == SCALAR)
                    sL2test(tGrid,u0,u,err,sL2sq);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case IP2C:    if (structure == SCALAR)
                    eprintf("Error: sL2_error not available.\n");
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case P2C_ELBUB: if (structure == SCALAR)
                    eprintf("Error: sL2_error not available.\n");
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   case GP1C:
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1C:
   case GQ1X4C:
   case GQ1C_ELBUB:
   case GP2C:
   case GP2X3C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:
   case GQ2C:
   case GQ2X4C:
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB:
                 if (structure == P_SCALAR)
                    eprintf("Error: sL2_error not available.\n");
                 else if (structure == SCALAR)
                    general_sL2test(tGrid,u0,u,L2_Q_RULE,REF_MAP,ELEM,
                                    err,general_sL2);
                 else
                    eprintf("Error: sL2_error not available.\n");
        break;
   default:
        eprintf("Error: sL2_error not available.\n");
        break;
   }
}

void vL2_errors(tGrid,u,v,err,errlin,errbubble,print,choice,space)
GRID *tGrid;
INT u, v, print, choice, space;
FLOAT *err, *errlin, *errbubble;
{
   switch(space){
   case P1_NC:   nc_L2test(tGrid,u,err);
        break;
   case P1_MOD:  dvector_to_vector_f(tGrid,u,v,0);
                 nc_L2test(tGrid,v,err);
        break;
   case P1C:     L2test(tGrid,print,u,err,vp1cL2s);
        break;
   case P1C_FBUB: switch(choice){
                    case 1: L2test(tGrid,print,u,err,oldL2s);
                         break;
                    case ALL: L2errors(tGrid,print,u,u,err,errlin,errbubble,
                                       oldL2s,linL2s,oldbubbleL2s);
                         break;
                    default:
                         eprintf("Error: vL2_errors not available.\n");
                         break;
                    }
        break;
   case P1C_NEW_FBUB: switch(choice){
                    case 1: L2test(tGrid,print,u,err,newL2s);
                         break;
                    case ALL: L2errors(tGrid,print,u,u,err,errlin,errbubble,
                                       newL2s,linL2s,newbubbleL2s);
                         break;
                    default:
                         eprintf("Error: vL2_errors not available.\n");
                         break;
                    }
        break;
   case P1C_ELBUB: switch(choice){
                    case 1: L2test(tGrid,print,u,err,miniL2s);
                         break;
                    case ALL: L2errors(tGrid,print,u,u,err,errlin,errbubble,
                                       miniL2s,linL2s,minibubbleL2s);
                         break;
                    default:
                         eprintf("Error: vL2_errors not available.\n");
                         break;
                    }
        break;
   case Q1ROT:   L2test(tGrid,print,u,err,L2s_q1rot);
        break;
   case P2C:     L2test(tGrid,print,u,err,L2sq);
        break;
   case IP2C:    L2test(tGrid,print,u,err,L2sq_iso);
        break;
   case P2C_ELBUB: L2test(tGrid,print,u,err,L2sq);
        break;
   case IP2C_ELBUB: L2test(tGrid,print,u,err,L2sq_iso);
        break;
   default:
        eprintf("Error: vL2_errors not available.\n");
        break;
   }
}

void sH10_error(tGrid,u,v,err,u01,u02,u03,space)
GRID *tGrid;
INT u, v, space;
FLOAT *err, (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case P1_NC:   nc_H10err(tGrid,u01,u02,u03,u,err);
        break;
   case P1_MOD:  vector_to_scalar_f(tGrid,u,v,0);
                 nc_H10err(tGrid,u01,u02,u03,v,err);
        break;
   case P1C:     sH10test(tGrid,u01,u02,u03,u,err,sH10);
        break;
   case Q1C:     sH10test(tGrid,u01,u02,u03,u,err,sH10_q1);
        break;
   case Q1ROT:   sH10test(tGrid,u01,u02,u03,u,err,sH10_q1rot);
        break;
   case P2C:     sH10test(tGrid,u01,u02,u03,u,err,sH10q);
        break;
   case GP1C:
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1C:
   case GQ1X4C:
   case GQ1C_ELBUB:
   case GP2C:
   case GP2X3C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:
   case GQ2C:
   case GQ2X4C:
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB:
                 general_sH10_test(tGrid,u01,u02,u03,u,
                                   H10_Q_RULE,REF_MAP,ELEM,err,general_sH10);
        break;
   default:
        eprintf("Error: sH10_error not available.\n");
        break;
   }
}

void vH10_errors(tGrid,u,v,err,errlin,errbubble,choice,space)
GRID *tGrid;
INT u, v, choice, space;
FLOAT *err, *errlin, *errbubble;
{
   switch(space){
   case P1_NC:   nc_H10test(tGrid,u,err);
        break;
   case P1_MOD:  dvector_to_vector_f(tGrid,u,v,0);
                 nc_H10test(tGrid,v,err);
        break;
   case Q1ROT:   q1rot_H10test(tGrid,u,err);
        break;
   case P1C:     H10err(tGrid,u,err,H10ijlin_3,3);
        break;
   case P1C_FBUB: switch(choice){
                    case 1:
                         if(DIM == 2)
                            H10err(tGrid,u,err,oldH10ij_2,2);
                         else if(DIM == 3)
                            H10err(tGrid,u,err,oldH10ij_3,3);
                         break;
                    case ALL: 
                         if(DIM == 2)
                            H10errors(tGrid,u,err,errlin,errbubble,
                                      oldH10ij_2,H10ijlin_2,oldbubbleH10ij_2,2);
                         else if(DIM == 3)
                            H10errors(tGrid,u,err,errlin,errbubble,
                                      oldH10ij_3,H10ijlin_3,oldbubbleH10ij_3,3);
                         break;
                    default:
                         eprintf("Error: vH10_errors not available.\n");
                         break;
                    }
        break;
   case P1C_NEW_FBUB: switch(choice){
                    case 1: H10err_p1c_new_fbub(tGrid,u,err);
                         break;
                    case ALL: 
                         if(DIM == 2)
                            H10errors_p1c_new_fbub(tGrid,u,err,errlin,errbubble,
                                                                  H10ijlin_2,2);
                         else if(DIM == 3)
                            H10errors_p1c_new_fbub(tGrid,u,err,errlin,errbubble,
                                                                  H10ijlin_3,3);
                         break;
                    default:
                         eprintf("Error: vH10_errors not available.\n");
                         break;
                    }
        break;
   case P1C_ELBUB: switch(choice){
                    case 1:   H10err(tGrid,u,err,miniH10ij_3,3);
                         break;
                    case 2:   H10err(tGrid,u,errlin,H10ijlin_3,3);
                         break;
                    case ALL: H10errors(tGrid,u,err,errlin,errbubble,
                                    miniH10ij_3,H10ijlin_3,minibubbleH10ij_3,3);
                         break;
                    default:
                         eprintf("Error: vH10_errors not available.\n");
                         break;
                    }
        break;
   case P2C:     H10qtest(tGrid,u,err);
        break;
   case IP2C:    H10qtest_iso(tGrid,u,err);
        break;
   case P2C_ELBUB: H10qtest(tGrid,u,err);
        break;
   case IP2C_ELBUB: H10qtest_iso(tGrid,u,err);
        break;
   default:
        eprintf("Error: vH10_errors not available.\n");
        break;
   }
}

void s_conv_error(tGrid,u,v,err,bb0,bb1,bb2,u01,u02,u03,space) 
GRID *tGrid;
INT u, v, space;
FLOAT *err, (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case P1_NC:   nc_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err);
        break;
   case P1_MOD:  vector_to_scalar_f(tGrid,u,v,0);
                 nc_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,v,err);
        break;
   case P1C:     s_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err,p1c_L2conv);
        break;
   case Q1C:     s_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err,q1c_L2conv);
        break;
   case P2C:     s_conv_err(tGrid,bb0,bb1,bb2,u01,u02,u03,u,err,p2c_L2conv);
        break;
   case GP1C:
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1C:
   case GQ1X4C:
   case GQ1C_ELBUB:
   case GP2C:
   case GP2X3C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:
   case GQ2C:
   case GQ2X4C:
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB:
                 general_sL2conv_test(tGrid,bb0,bb1,bb2,u01,u02,u03,u,
                                  CONV_Q_RULE,REF_MAP,ELEM,err,general_sL2conv);
        break;
   default:
        eprintf("Error: s_conv_error not available.\n");
        break;
   }
}

void s_sd_error(tGrid,u,v,err,eps,bb0,bb1,bb2,u0,u01,u02,u03,space) 
GRID *tGrid;
INT u, v, space;
FLOAT *err, eps, (*bb0)(), (*bb1)(), (*bb2)(), 
      (*u0)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case P1_NC:   nc_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err);
        break;
   case P1_MOD:  vector_to_scalar_f(tGrid,u,v,0);
                 nc_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,v,err);
        break;
   case P1C:     s_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err,
                                                          sH10,p1c_L2conv,sL2s);
        break;
   case Q1C:     s_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err,
                                                     sH10_q1,q1c_L2conv,sL2_q1);
        break;
   case P2C:     s_sd_err(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,err,
                                                        sH10q,p2c_L2conv,sL2sq);
        break;
   case GP1C:
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1C:
   case GQ1X4C:
   case GQ1C_ELBUB:
   case GP2C:
   case GP2X3C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:
   case GQ2C:
   case GQ2X4C:
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB:
                 general_s_sd_test(tGrid,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,space,
                                   SD_Q_RULE,REF_MAP,ELEM,err,
                                   general_sH10,general_sL2conv,general_sL2);
        break;
   default:
        eprintf("Error: s_sd_error not available.\n");
        break;
   }
}

void sL_infinity_error(tGrid,u,v,err,u0,space)
GRID *tGrid;
INT u, v, space;
FLOAT *err, (*u0)();
{
   switch(space){
   case P1_NC:   nc_L_inf_error(tGrid,u0,u,err);
        break;
   case P1_MOD:  vector_to_scalar_f(tGrid,u,v,0);
                 nc_L_inf_error(tGrid,u0,v,err);
        break;
   case P1C:     p1c_L_inf_error(tGrid,u0,u,err);
        break;
   case Q1C:     p1c_L_inf_error(tGrid,u0,u,err);
        break;
   case P2C:     p2c_L_inf_error(tGrid,u0,u,err);
        break;
   case GP1C:    gs_max_error(tGrid,u,u0,err,LP1,3,p1c_element,REF_MAP);
        break;
   case GP1X3C:  gs_max_error(tGrid,u,u0,err,LP3,10,p1x3c_element,REF_MAP);
        break;
   case GP1C_ELBUB: gs_max_error(tGrid,u,u0,err,LP3,10,p1cb_element,REF_MAP);
        break;
   case GQ1C:    gs_max_error(tGrid,u,u0,err,LQ1,4,q1c_element,REF_MAP);
        break;
   case GQ1X4C:  gs_max_error(tGrid,u,u0,err,LQ2,9,q1x4c_element,REF_MAP);
        break;
   case GQ1C_ELBUB: gs_max_error(tGrid,u,u0,err,LQ4,25,q1cb_element,REF_MAP);
        break;
   case GP2C:    gs_max_error(tGrid,u,u0,err,LP2,6,p2c_element,REF_MAP);
        break;
   case GP2X3C:  gs_max_error(tGrid,u,u0,err,LP6,28,p2x3c_element,REF_MAP);
        break;
   case GP2C_3ELBUB: gs_max_error(tGrid,u,u0,err,LP6,28,p2c3b_element,REF_MAP);
        break;
   case GP2C_6ELBUB: gs_max_error(tGrid,u,u0,err,LP6,28,p2c6b_element,REF_MAP);
        break;
   case GQ2C:    gs_max_error(tGrid,u,u0,err,LQ2,9,q2c_element,REF_MAP);
        break;
   case GQ2X4C:  gs_max_error(tGrid,u,u0,err,LQ4,25,q2x4c_element,REF_MAP);
        break;
   case GQ2C_2ELBUB: gs_max_error(tGrid,u,u0,err,LQ4,25,q2c2b_element,REF_MAP);
        break;
   case GQ2C_3ELBUB: gs_max_error(tGrid,u,u0,err,LQ4,25,q2c3b_element,REF_MAP);
        break;
//               general_L_inf_error(tGrid,u0,u,err,REF_MAP,ELEM);
   default:
        eprintf("Error: sL_infinity_error not available.\n");
        break;
   }
}

void sL_infinity_error_without_belems(tGrid,u,v,err,u0,space)
GRID *tGrid;
INT u, v, space;
FLOAT *err, (*u0)();
{
   switch(space){
   case GP1C:    gs_max_error_wb(tGrid,u,u0,err,LP1,3,p1c_element,REF_MAP);
        break;
   case GP1X3C:  gs_max_error_wb(tGrid,u,u0,err,LP3,10,p1x3c_element,REF_MAP);
        break;
   case GP1C_ELBUB: gs_max_error_wb(tGrid,u,u0,err,LP3,10,p1cb_element,REF_MAP);
        break;
   case GQ1C:    gs_max_error_wb(tGrid,u,u0,err,LQ1,4,q1c_element,REF_MAP);
        break;
   case GQ1X4C:  gs_max_error_wb(tGrid,u,u0,err,LQ2,9,q1x4c_element,REF_MAP);
        break;
   case GQ1C_ELBUB: gs_max_error_wb(tGrid,u,u0,err,LQ4,25,q1cb_element,REF_MAP);
        break;
   case GP2C:    gs_max_error_wb(tGrid,u,u0,err,LP2,6,p2c_element,REF_MAP);
        break;
   case GP2X3C:  gs_max_error_wb(tGrid,u,u0,err,LP6,28,p2x3c_element,REF_MAP);
        break;
   case GP2C_3ELBUB: gs_max_error_wb(tGrid,u,u0,err,LP6,28,p2c3b_element,REF_MAP);
        break;
   case GP2C_6ELBUB: gs_max_error_wb(tGrid,u,u0,err,LP6,28,p2c6b_element,REF_MAP);
        break;
   case GQ2C:    gs_max_error_wb(tGrid,u,u0,err,LQ2,9,q2c_element,REF_MAP);
        break;
   case GQ2X4C:  gs_max_error_wb(tGrid,u,u0,err,LQ4,25,q2x4c_element,REF_MAP);
        break;
   case GQ2C_2ELBUB: gs_max_error_wb(tGrid,u,u0,err,LQ4,25,q2c2b_element,REF_MAP);
        break;
   case GQ2C_3ELBUB: gs_max_error_wb(tGrid,u,u0,err,LQ4,25,q2c3b_element,REF_MAP);
        break;
   default:
        eprintf("Error: sL_infinity_error_without_belems not available.\n");
        break;
   }
}

void vL_infinity_error(tGrid,u,v,err,space)
GRID *tGrid;
INT u, v, space;
FLOAT *err;
{
   switch(space){
   case Q1ROT:   if (DOF == MIDPOINT)
                    L_infinity_error_q1rot(tGrid,u01,u02,u,err,vset_edge_value);
                 else if (DOF == MEAN_VALUE)
                    L_infinity_error_q1rot(tGrid,u01,u02,u,err,vset_edge_mean_value);
                 else
                    eprintf("Error: vL_infinity_error not available.\n");
        break;
   default:
        eprintf("Error: vL_infinity_error not available.\n");
        break;
   }
}

void errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,u0,u01,u02,u03,u,v,
                                            err1,err2,err3,err4,err5,area,space)
GRID *tGrid;                     /*  if DIM == 2, bb2 and u03 are not used    */
FLOAT xm, ym, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u0)(), (*u01)(), (*u02)(), 
      (*u03)(), *err1, *err2, *err3, *err4, *err5, *area;
INT u, v, space;
{
   switch(space){
   case P1_NC:   nc_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,
                                u0,u01,u02,u03,u,err1,err2,err3,err4,err5,area);
        break;
   case P1_MOD:  vector_to_scalar_f(tGrid,u,v,0);
                 nc_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,
                                u0,u01,u02,u03,v,err1,err2,err3,err4,err5,area);
        break;
   case P1C:     s_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,
                                 u0,u01,u02,u03,u,err1,err2,err3,err4,err5,area,
                                                 sL2s,sH10,p1c_L2conv,p1c_Linf);
        break;
   case Q1C:     s_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,
                                 u0,u01,u02,u03,u,err1,err2,err3,err4,err5,area,
                                            sL2_q1,sH10_q1,q1c_L2conv,q1c_Linf);
        break;
   case P2C:     s_errors_on_square(tGrid,xm,ym,eps,bb0,bb1,bb2,
                                 u0,u01,u02,u03,u,err1,err2,err3,err4,err5,area,
                                               sL2sq,sH10q,p2c_L2conv,p2c_Linf);
        break;
   default:
        *err1 = *err2 = *err3 = *err4 = *err5 = *area = -1.;
        eprintf("Error: errors_on_square not available.\n");
        break;
   }
}

void s_lp_error(tGrid,u,err,eps,bb0,bb1,bb2,u01,u02,u03,space)
GRID *tGrid;
INT u, space;
FLOAT *err, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1X4C:
   case GQ1C_ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,1,
                        bb0,bb1,one_fcn,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GP2X3C:
   case GP2C_3ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GP2C_6ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GQ2X4C:      general_lp_test(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
/*
   case GQ2X4C:      general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
*/
        break;
   case GQ2C_2ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GQ2C_3ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   default:
        eprintf("Error: s_lp_error not available.\n");
        break;
   }
}

void s_lp_error_without_belems(tGrid,u,err,eps,bb0,bb1,bb2,u01,u02,u03,space)
GRID *tGrid;
INT u, space;
FLOAT *err, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1X4C:
   case GQ1C_ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,1,
                        bb0,bb1,one_fcn,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GP2X3C:
   case GP2C_3ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GP2C_6ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GQ2X4C:      general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
/*
   case GQ2X4C:      general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
*/
        break;
   case GQ2C_2ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   case GQ2C_3ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_lp_err);
        break;
   default:
        eprintf("Error: s_lp_error not available.\n");
        break;
   }
}

void s_conv_lp_error(tGrid,u,err,eps,bb0,bb1,bb2,u01,u02,u03,space)
GRID *tGrid;
INT u, space;
FLOAT *err, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case GP1C:
   case GQ1C:       general_macro_lp_test(tGrid,u01,u02,u03,u,eps,1,
                        bb0,bb1,one_fcn,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_macro_conv_lp_err);
        break;
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1X4C:
   case GQ1C_ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,1,
                        bb0,bb1,one_fcn,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GP2C:        general_macro_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_macro_conv_lp_err);
        break;
   case GP2X3C:
   case GP2C_3ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GP2C_6ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GQ2X4C:      general_lp_test(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
/*
   case GQ2X4C:      general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
*/
        break;
   case GQ2C_2ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GQ2C_3ELBUB: general_lp_test(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   default:
        eprintf("Error: s_conv_lp_error not available.\n");
        break;
   }
}

void s_conv_lp_error_without_belems(tGrid,u,err,eps,bb0,bb1,bb2,u01,u02,u03,space)
GRID *tGrid;
INT u, space;
FLOAT *err, eps, (*bb0)(), (*bb1)(), (*bb2)(), (*u01)(), (*u02)(), (*u03)();
{
   switch(space){
   case GP1X3C:
   case GP1C_ELBUB:
   case GQ1X4C:
   case GQ1C_ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,1,
                        bb0,bb1,one_fcn,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GP2X3C:
   case GP2C_3ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GP2C_6ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,r_l0,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GQ2X4C:      general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
/*
   case GQ2X4C:      general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
*/
        break;
   case GQ2C_2ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,3,
                        bb0,bb1,one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   case GQ2C_3ELBUB: general_lp_test_without_belems(tGrid,u01,u02,u03,u,eps,4,
                        bb0,bb1,one_fcn,r_l1,r_l2,r_l1l2,NULL,NULL,NULL,NULL,
                        L2_Q_RULE,REF_MAP,ELEM,err,general_conv_lp_err);
        break;
   default:
        eprintf("Error: s_conv_lp_error not available.\n");
        break;
   }
}

/*

pL2test_lin(tGrid,err,p,eu,t)

L2 div-norm:
============

   if (U_SPACE == P1_NC || U_SPACE == P1_MOD){
      if (U_SPACE == P1_MOD)
         dvector_to_vector_f(topGrid,U,U,0);
      nc_lin_div_norm(topGrid,U,&err); 
   }
   else if (U_SPACE == P1C_FBUB){
      lin_div_norm(TOP_GRID(mg),U,&err);
   }

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

L2 norm:
========

   if (U_SPACE == P1C_FBUB){
      linL2norm(TOP_GRID(mg),1,U,&err);
   }

   if (P_SPACE == P0){
      pL2norm(TOP_GRID(mg),U,&err,T_FOR_P);
   }
   else if (P_SPACE == P1C){
      pL2norm_cont(topGrid,&err,U,T_FOR_P);
   }

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

H10 norm:
=========

   if (U_SPACE == P1C_FBUB){
      linH10norm(TOP_GRID(mg),U,&err,U_SPACE);  
   }

*/
