/******************************************************************************/
/*                                                                            */
/*                             boundary condition                             */
/*                                                                            */
/******************************************************************************/

void set_value();
void mv_set_value_n();
void mv_set_value_f();

FLOAT integrate_using_multiple_Simpson(x1,x2,f,n)   /*  integral of f over a  */
FLOAT *x1, *x2, (*f)();                /*  line connecting x1 and x2 divided  */
INT n;                                 /*  by the distance between x1 and x2  */
{
   FLOAT d[DIM], d1[DIM], x[DIM], s=0., z=0.;
   INT i;

   SET12(d,x2,x1,n)
   SET2(d1,d,0.5)
   for (i=1; i < n; i++){
      SET23(x,x1,d,i)
      s += f(x);
      SET6(x,d1)
      z += f(x);
   }
   SET11(x,x2,d1)
   z += f(x);
   return((f(x1)+f(x2)+2.*s+4.*z)/6./n);
}

#if N_DATA & SCALAR_NODE_DATA

void sbound_val1(mg,u,u0)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   NODE *pnode;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            NDS(pnode,u) = u0(pnode->myvertex->x);
}

#else

void sbound_val1(mg,u,u0)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: sbound_val1 not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void vbound_val1(mg,u,u0,u1)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)();
INT u;
{
   GRID *theGrid;
   NODE *pnode;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
         }
}

#else

void vbound_val1(mg,u,u0,u1)
MULTIGRID *mg; FLOAT (*u0)(), (*u1)(); INT u;
{  eprintf("Error: vbound_val1 not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void sbound_val2(mg,u,u0)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            NDS(pnode,u) = u0(pnode->myvertex->x);

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            if (NOT_FF(fa0)){
               AVERAGE(x1,x2,x12);
               FD(fa0,u) = 4.*u0(x12) - 2.*(NDS(n1,u)+NDS(n2,u));
            }
            if (NOT_FF(fa1)){
               AVERAGE(x0,x2,x02);
               FD(fa1,u) = 4.*u0(x02) - 2.*(NDS(n0,u)+NDS(n2,u));
            }
            if (NOT_FF(fa2)){
               AVERAGE(x0,x1,x01);
               FD(fa2,u) = 4.*u0(x01) - 2.*(NDS(n0,u)+NDS(n1,u));
            }
         }
}

#else

void sbound_val2(mg,u,u0)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: sbound_val2 not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (F_DATA & CURVED_FACE_MIDDLE) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void sbound_val2_iso(mg,u,u0)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            NDS(pnode,u) = u0(pnode->myvertex->x);

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            if (NOT_FF(fa0)){
               if (fa0->c_midpoint)
                  SET1(x12,fa0->c_midpoint->x)
               else
                  AVERAGE(x1,x2,x12);
               FD(fa0,u) = 4.*u0(x12) - 2.*(NDS(n1,u)+NDS(n2,u));
            }
            if (NOT_FF(fa1)){
               if (fa1->c_midpoint)
                  SET1(x02,fa1->c_midpoint->x)
               else
                  AVERAGE(x0,x2,x02);
               FD(fa1,u) = 4.*u0(x02) - 2.*(NDS(n0,u)+NDS(n2,u));
            }
            if (NOT_FF(fa2)){
               if (fa2->c_midpoint)
                  SET1(x01,fa2->c_midpoint->x)
               else
                  AVERAGE(x0,x1,x01);
               FD(fa2,u) = 4.*u0(x01) - 2.*(NDS(n0,u)+NDS(n1,u));
            }
         }
}

#else

void sbound_val2_iso(mg,u,u0)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: sbound_val2_iso not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void bound_val_vn_vf_2(mg,u,u0,u1)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
         }

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            if (NOT_FF(fa0)){
               AVERAGE(x1,x2,x12);
               FDV(fa0,u,0) = 4.*u0(x12) - 2.*(ND(n1,u,0)+ND(n2,u,0));
               FDV(fa0,u,1) = 4.*u1(x12) - 2.*(ND(n1,u,1)+ND(n2,u,1));
            }
            if (NOT_FF(fa1)){
               AVERAGE(x0,x2,x02);
               FDV(fa1,u,0) = 4.*u0(x02) - 2.*(ND(n0,u,0)+ND(n2,u,0));
               FDV(fa1,u,1) = 4.*u1(x02) - 2.*(ND(n0,u,1)+ND(n2,u,1));
            }
            if (NOT_FF(fa2)){
               AVERAGE(x0,x1,x01);
               FDV(fa2,u,0) = 4.*u0(x01) - 2.*(ND(n0,u,0)+ND(n1,u,0));
               FDV(fa2,u,1) = 4.*u1(x01) - 2.*(ND(n0,u,1)+ND(n1,u,1));
            }
         }
}

void correct_edge_2_5(f0,n1,n2,u)
FACE *f0;
NODE *n1, *n2;  /*  nodes of f0  */
INT u;
{
   FLOAT *x1, *x2, x11112[DIM], x11122[DIM], x11222[DIM], x12222[DIM], nn[DIM], s;
   
   if (NOT_FF(f0)){
      x1 = n1->myvertex->x;
      x2 = n2->myvertex->x;
      normal_vector(x1,x2,f0,nn); 
      ADDMULT(.8,x1,.2,x2,x11112);
      ADDMULT(.6,x1,.4,x2,x11122);
      ADDMULT(.4,x1,.6,x2,x11222);
      ADDMULT(.2,x1,.8,x2,x12222);
      s = (( nn[0]*( 19.*(u01(x1)+u01(x2)) + 75.*(u01(x11112)+u01(x12222)) +
                     50.*(u01(x11122)+u01(x11222)) ) +
             nn[1]*( 19.*(u02(x1)+u02(x2)) + 75.*(u02(x11112)+u02(x12222)) +
                     50.*(u02(x11122)+u02(x11222)) ))/288.
             - ( DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn) )*.5  
             - DOT(FDVP(f0,u),nn)/6. )*6.;
      SET4(FDVP(f0,u),nn,s)
   }
}

void correct_boundary_values_2_5(mg,u)
MULTIGRID *mg;
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2;
   FACE *f0, *f1, *f2;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer){
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(f0,f1,f2,pelem);
            correct_edge_2_5(f0,n1,n2,u);
            correct_edge_2_5(f1,n0,n2,u);
            correct_edge_2_5(f2,n0,n1,u);
         }
   }
}

#else

void bound_val_vn_vf_2(mg,u,u0,u1)
MULTIGRID *mg; FLOAT (*u0)(), (*u1)(); INT u;
{  eprintf("Error: bound_val_vn_vf_2 not available.\n");  }

void correct_boundary_values_2_5(mg,u)
MULTIGRID *mg; INT u;
{  eprintf("Error: correct_boundary_values_2_5 not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (F_DATA & CURVED_FACE_MIDDLE) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void bound_val_vn_vf_2_iso(mg,u,u0,u1)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
         }

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            if (NOT_FF(fa0)){
               if (fa0->c_midpoint)
                  SET1(x12,fa0->c_midpoint->x)
               else
                  AVERAGE(x1,x2,x12);
               FDV(fa0,u,0) = 4.*u0(x12) - 2.*(ND(n1,u,0)+ND(n2,u,0));
               FDV(fa0,u,1) = 4.*u1(x12) - 2.*(ND(n1,u,1)+ND(n2,u,1));
            }
            if (NOT_FF(fa1)){
               if (fa1->c_midpoint)
                  SET1(x02,fa1->c_midpoint->x)
               else
                  AVERAGE(x0,x2,x02);
               FDV(fa1,u,0) = 4.*u0(x02) - 2.*(ND(n0,u,0)+ND(n2,u,0));
               FDV(fa1,u,1) = 4.*u1(x02) - 2.*(ND(n0,u,1)+ND(n2,u,1));
            }
            if (NOT_FF(fa2)){
               if (fa2->c_midpoint)
                  SET1(x01,fa2->c_midpoint->x)
               else
                  AVERAGE(x0,x1,x01);
               FDV(fa2,u,0) = 4.*u0(x01) - 2.*(ND(n0,u,0)+ND(n1,u,0));
               FDV(fa2,u,1) = 4.*u1(x01) - 2.*(ND(n0,u,1)+ND(n1,u,1));
            }
         }
}

#else

void bound_val_vn_vf_2_iso(mg,u,u0,u1)
MULTIGRID *mg; FLOAT (*u0)(), (*u1)(); INT u;
{  eprintf("Error: bound_val_vn_vf_2_iso not available.\n");  }

#endif

#if DIM == 3

#if N_DATA & VECTOR_NODE_DATA

FLOAT l_face_flux(n1,n2,n3,f4,u)
NODE *n1, *n2, *n3;
FACE *f4;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f4)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,
                           n3->myvertex->x,f4,nn);
      return((DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn)+DOT(NDD(n3,u),nn))*area/3.);
   }
   else
      return(0.);
}

#else

FLOAT l_face_flux(n1,n2,n3,f4,u)
NODE *n1, *n2, *n3; FACE *f4; INT u;
{  eprintf("Error: l_face_flux not available.\n"); return(0.);  } 

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT lb_face_flux(n1,n2,n3,f4,u)
NODE *n1, *n2, *n3;
FACE *f4;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f4)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,
                           n3->myvertex->x,f4,nn);
      return((DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn)+DOT(NDD(n3,u),nn))*area/3.
              + FMULT*FD(f4,u)*area/60. );
   }
   else
      return(0.);
}

#else

FLOAT lb_face_flux(n1,n2,n3,f4,u)
NODE *n1, *n2, *n3; FACE *f4; INT u;
{  eprintf("Error: lb_face_flux not available.\n"); return(0.);  }

#endif

#if ELEMENT_TYPE == SIMPLEX

FLOAT boundary_flux(mg,u,face_flux)
MULTIGRID *mg;
INT u;
FLOAT (*face_flux)();
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3, *n4;
   FACE *f1, *f2, *f3, *f4;
   FLOAT sum=0.;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n1,n2,n3,n4,pelem);
            FACES_OF_ELEMENT(f1,f2,f3,f4,pelem);
            sum += face_flux(n2,n3,n4,f1,u) + face_flux(n1,n3,n4,f2,u) +
                   face_flux(n1,n2,n4,f3,u) + face_flux(n1,n2,n3,f4,u);
         }
   return(sum);
}

#else

FLOAT boundary_flux(mg,u,face_flux)
MULTIGRID *mg; INT u; FLOAT (*face_flux)();
{  eprintf("Error: boundary_flux not available.\n"); return(-1.);  }

#endif

FLOAT vlb_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: vlb_face_flux not available.\n"); return(0.);  }

#else  /*  DIM != 3  */

#if N_DATA & VECTOR_NODE_DATA

FLOAT l_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return((DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn))*area/2.);
   }
   else
      return(0.);
}

#else

FLOAT l_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: l_face_flux not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

FLOAT lb_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return((DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn))*area/2.
              + FMULT*FD(f3,u)*area/6. );
   }
   else
      return(0.);
}

#else

FLOAT lb_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: lb_face_flux not available.\n"); return(0.);  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

FLOAT vlb_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return((DOT(NDD(n1,u),nn)+DOT(NDD(n2,u),nn))*area/2.
              + DOT(FDVP(f3,u),nn)*area/6. );
   }
   else
      return(0.);
}

#else

FLOAT vlb_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: vlb_face_flux not available.\n"); return(0.);  }

#endif

#if ELEMENT_TYPE == SIMPLEX

FLOAT boundary_flux(mg,u,face_flux)
MULTIGRID *mg;
INT u;
FLOAT (*face_flux)();
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;
   FLOAT sum=0.;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n1,n2,n3,pelem);
            FACES_OF_ELEMENT(f1,f2,f3,pelem);
            sum += face_flux(n2,n3,f1,u) + face_flux(n1,n3,f2,u) +
                                           face_flux(n1,n2,f3,u);
         }
   return(sum);
}

#else

FLOAT boundary_flux(mg,u,face_flux)
MULTIGRID *mg; INT u; FLOAT (*face_flux)();
{  eprintf("Error: boundary_flux not available.\n"); return(-1.);  }

#endif

#endif

#if F_DATA & SCALAR_FACE_DATA

void bubble_correction(mg,tGrid,sum1,sum2,u)
MULTIGRID *mg;
GRID *tGrid;
FLOAT sum1, sum2;
INT u;
{
   GRID *theGrid;
   FACE *pface;
   FLOAT flux, e;
   
   printf("Boundary condition with bubbles.\n");
   printf("Flux of lin. int. of u0 through the boundary: %e\n",
                                               boundary_flux(mg,u,l_face_flux));
   printf("Flux of lin. int. + bubbles: %e\n",boundary_flux(mg,u,lb_face_flux));
   printf("Approximative flux of u0 through the boundary: %e\n",sum2);
   /* correction of coefficients at face basis functions in order to obtain the 
      correct total flux through the boundary */
   flux = 0.0;
   if (fabs(sum1)<EPSC)
      printf("Correction cannot be done.\n");
   else{   
      e = (flux - sum2)/sum1;
      for (theGrid = tGrid; theGrid != NULL; theGrid = theGrid->finer)
        for (pface=FDBF(theGrid); pface !=FIRSTFACE(theGrid); pface=pface->succ)
           if (IS_TOP_FACE(pface))
              FD(pface,u) += FD(pface,u)*e;
      printf("Correction done ( e = %e, flux = %e ).\n",e,
                                              boundary_flux(mg,u,lb_face_flux));
   }
}

#else

void bubble_correction(mg,tGrid,sum1,sum2,u)
MULTIGRID *mg; GRID *tGrid; FLOAT sum1, sum2; INT u;
{  eprintf("Error: bubble_correction not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void node_correction(mg,u,v,w)
MULTIGRID *mg;
INT u, v, w;
{
   GRID *theGrid;
   NODE *pnode;
   FLOAT e, q, sp = 0., sm = 0.;
   
   printf("Boundary condition without bubbles.\n");
   printf("Flux of lin. int. of u0 through the boundary: %e\n",
                                               boundary_flux(mg,u,l_face_flux));
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            q = ND(pnode,w,0)/DOT(NDD(pnode,v),NDD(pnode,v));
            SET4(NDD(pnode,u),NDD(pnode,v),q)
            q = DOT(NDD(pnode,u),NDD(pnode,v));
            if (q < 0.){
               ND(pnode,w,0) = -1.;
               sm += q;
            }
            else{
               ND(pnode,w,0) = 1.;
               sp += q;
            }
         }
   printf("Flux of modified lin. int. of u0 through the boundary: %e\n",
                                               boundary_flux(mg,u,l_face_flux));
   printf("Approximative flux of u0 through the boundary: %e\n",sp+sm);
   if (fabs(sp+sm)<EPSC)
      printf("Correction cannot be done.\n");
   else{   
      e = (sp+sm)/(sp-sm);
      for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
         for (pnode = FDBN(theGrid); pnode != FIRSTNODE(theGrid);
                                                            pnode = pnode->succ)
            if (IS_TOP_NODE(pnode)){
               if (ND(pnode,w,0) < 0.)
                   SET4(NDD(pnode,u),NDD(pnode,u),e)
               else
                   SET9(NDD(pnode,u),NDD(pnode,u),e)
            }
      printf("Correction done ( e = %e, flux = %e ).\n",e,
                                               boundary_flux(mg,u,l_face_flux));
   }
}

#else

void node_correction(mg,u,v,w)
MULTIGRID *mg; INT u, v, w;
{  eprintf("Error: node_correction not available.\n");  }

#endif

#if DIM == 3

 /* s1 is a linear approximation of (\int_K' u0\cdot nn \ds)/K'; s2 is a better 
    (quadratic) approximation of this integral. */
 void integr2(x1,x2,x3,nn,s1,s2)
 FLOAT *x1, *x2, *x3, *nn, *s1, *s2;
 {
   FLOAT xs1[DIM], xs2[DIM], xs3[DIM];
   
   AVERAGE(x1,x2,xs1);
   AVERAGE(x2,x3,xs2);
   AVERAGE(x3,x1,xs3);
   
   *s1 = ( nn[0]*(u01(x1)+u01(x2)+u01(x3))+
           nn[1]*(u02(x1)+u02(x2)+u02(x3))+
           nn[2]*(u03(x1)+u03(x2)+u03(x3)) )/3.0;
          
   *s2 = ( nn[0]*(u01(xs1)+u01(xs2)+u01(xs3))+
           nn[1]*(u02(xs1)+u02(xs2)+u02(xs3))+
           nn[2]*(u03(xs1)+u03(xs2)+u03(xs3)) )/3.0;       
 }
 
#if PROBLEM & BUBBLE_BC

#if !(DATA_STR & LG_DATA)

#if F_DATA & SCALAR_FACE_DATA

 /* x1, x2, x3 are vertexes of the i-th face, x0 is the remaining vertex */
 void face_val(pelement,i,x1,x2,x3,x0,u,sum1,sum2)
 ELEMENT *pelement;
 INT i, u;
 FLOAT  *x1, *x2, *x3, *x0, *sum1, *sum2;  
 {                                /* sum1=\int_{\partial\om} p\cdot\nn \ds    */
    FLOAT nn[DIM], k0, k1, area;  /* sum2 ~ \int_{\partial\om} u0\cdot\nn \ds */
    
    if (NOT_FF(pelement->f[i])){
       area = normal_vector(x1,x2,x3,pelement->f[i],nn); 
       integr2(x1,x2,x3,nn,&k0,&k1);  
       FD(pelement->f[i],u) = 60.0*(k1-k0)/FMULT;
       test_normal_vector_direction(nn,x0,x1,&area);
       *sum1 += area*FD(pelement->f[i],u);
       *sum2 += area*k1;
    }
 }

#else

 void face_val(pelement,i,x1,x2,x3,x0,u,sum1,sum2)
 ELEMENT *pelement; INT i, u; FLOAT  *x1, *x2, *x3, *x0, *sum1, *sum2;  
 {  eprintf("Error: face_val not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

 void bound_val(mg,u,v,w) /* v, w not used */
 MULTIGRID *mg;
 INT u, v, w;
 {
   GRID *theGrid;
   NODE *pnode;
   ELEMENT *pelem;
   FLOAT *x0, *x1, *x2, *x3, sum1=0.0, sum2=0.0;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer){
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u01(pnode->myvertex->x);
            ND(pnode,u,1) = u02(pnode->myvertex->x);
            ND(pnode,u,2) = u03(pnode->myvertex->x);
         }
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            VERTICES_OF_ELEMENT(x0,x1,x2,x3,pelem);
            face_val(pelem,0,x1,x2,x3,x0,u,&sum1,&sum2);
            face_val(pelem,1,x0,x2,x3,x1,u,&sum1,&sum2);
            face_val(pelem,2,x0,x1,x3,x2,u,&sum1,&sum2);
            face_val(pelem,3,x0,x1,x2,x3,u,&sum1,&sum2);
         }
   }
   sum1 = sum1/60.0*FMULT; /* flux of the bubble part through the boundary */
   bubble_correction(mg,FIRSTGRID(mg),sum1,sum2,u);
 }

#else

 void bound_val(mg,u,v,w)
 MULTIGRID *mg; INT u, v, w;
 {  eprintf("Error: bound_val not available.\n");  }

#endif

#endif

#if DATA_STR & LG_DATA

 /* x1, x2, x3 are vertexes of the i-th face, x0 is the remaining vertex */
 void face_val(pel,i,x1,x2,x3,x0,u,sum1,sum2)
 ELEMENT *pel;
 INT i, u;
 FLOAT  *x1, *x2, *x3, *x0, *sum1, *sum2;
 {                                /* sum1 = \int_{\partial\om} p \cdot\nn \ds */
    FLOAT nn[DIM], k0, k1, area;  /* sum2 ~ \int_{\partial\om} u0\cdot\nn \ds */

    if (NOT_FF(pel->f[i]))
       if (pel->n[0]->lgd && i != 0 || pel->n[1]->lgd && i != 1 || 
           pel->n[2]->lgd && i != 2 || pel->n[3]->lgd && i != 3)
          FD(pel->f[i],u) = 0.;
       else{
          area = normal_vector(x1,x2,x3,pel->f[i],nn); 
          integr2(x1,x2,x3,nn,&k0,&k1);  
          FD(pel->f[i],u) = 60.0*(k1-k0)/FMULT;
          test_normal_vector_direction(nn,x0,x1,&area);
          *sum1 += area*FD(pel->f[i],u);
          *sum2 += area*k1;
       }
 }

 void bound_val(mg,u,v,w) /* v, w not used */
 MULTIGRID *mg;
 INT u, v, w;
 {
   GRID *theGrid;
   NODE *pnode;
   ELEMENT *pelem;
   FLOAT *x0, *x1, *x2, *x3, sum1=0.0, sum2=0.0;
   
   theGrid = TOP_GRID(mg);
   for (pnode = FDBN(theGrid); pnode != FIRSTNODE(theGrid); pnode = pnode->succ)
      if (pnode->lgd)
         ND(pnode,u,0) = ND(pnode,u,1) = ND(pnode,u,2) = 0.;
      else{
         ND(pnode,u,0) = u01(pnode->myvertex->x);
         ND(pnode,u,1) = u02(pnode->myvertex->x);
         ND(pnode,u,2) = u03(pnode->myvertex->x);
      }
   for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,x3,pelem);
      face_val(pelem,0,x1,x2,x3,x0,u,&sum1,&sum2);
      face_val(pelem,1,x0,x2,x3,x1,u,&sum1,&sum2);
      face_val(pelem,2,x0,x1,x3,x2,u,&sum1,&sum2);
      face_val(pelem,3,x0,x1,x2,x3,u,&sum1,&sum2);
   }
   sum1 = sum1/60.0*FMULT;
   bubble_correction(mg,theGrid,sum1,sum2,u);
 }

#endif
#endif

#if !(PROBLEM & BUBBLE_BC)

#if N_DATA & VECTOR_NODE_DATA

void face_val(n1,n2,n3,f0,u,v,w)
NODE *n1, *n2, *n3;
FACE *f0;
INT u, v, w;
{
   FLOAT nn[DIM], area, s, q, *x1, *x2, *x3, x112[DIM], x113[DIM], x221[DIM], 
                                  x223[DIM], x331[DIM], x332[DIM], x123[DIM];
    
   if (NOT_FF(f0)){
      x1 = n1->myvertex->x;
      x2 = n2->myvertex->x;
      x3 = n3->myvertex->x;
      area = normal_vector(x1,x2,x3,f0,nn);
      POINT2(x1,x2,x112,2.0,3.0);
      POINT2(x1,x3,x113,2.0,3.0);
      POINT2(x2,x1,x221,2.0,3.0);
      POINT2(x2,x3,x223,2.0,3.0);
      POINT2(x3,x1,x331,2.0,3.0);
      POINT2(x3,x2,x332,2.0,3.0);
      POINT3(x1,x2,x3,x123);
      s = nn[0]*( 2.*u01(x123) - ND(n1,u,0) - ND(n2,u,0) - ND(n3,u,0) ) + 
          nn[1]*( 2.*u02(x123) - ND(n1,u,1) - ND(n2,u,1) - ND(n3,u,1) ) +
          nn[2]*( 2.*u03(x123) - ND(n1,u,2) - ND(n2,u,2) - ND(n3,u,2) );
      q = area*3./40.;
      ND(n1,w,0) += q*( nn[0]*( u01(x112) + u01(x113) - ND(n1,u,0) ) + 
                        nn[1]*( u02(x112) + u02(x113) - ND(n1,u,1) ) + 
                        nn[2]*( u03(x112) + u03(x113) - ND(n1,u,2) ) + s );
      ND(n2,w,0) += q*( nn[0]*( u01(x221) + u01(x223) - ND(n2,u,0) ) + 
                        nn[1]*( u02(x221) + u02(x223) - ND(n2,u,1) ) + 
                        nn[2]*( u03(x221) + u03(x223) - ND(n2,u,2) ) + s );
      ND(n3,w,0) += q*( nn[0]*( u01(x331) + u01(x332) - ND(n3,u,0) ) + 
                        nn[1]*( u02(x331) + u02(x332) - ND(n3,u,1) ) + 
                        nn[2]*( u03(x331) + u03(x332) - ND(n3,u,2) ) + s );
      q = area/3.;
      SET4(NDD(n1,v),nn,q)
      SET4(NDD(n2,v),nn,q)
      SET4(NDD(n3,v),nn,q)
   }
}

void bound_val(mg,u,v,w) /* Dirichlet b.c. has to be prescribed on the whole  */
MULTIGRID *mg;           /*                                          boundary */
INT u, v, w;
{
   GRID *theGrid;
   NODE *pnode;
   FACE *pface;
   ELEMENT *pelem;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u01(pnode->myvertex->x);
            ND(pnode,u,1) = u02(pnode->myvertex->x);
            ND(pnode,u,2) = u03(pnode->myvertex->x);
            ND(pnode,v,0) = ND(pnode,v,1) = ND(pnode,v,2) = ND(pnode,w,0) = 0.;
         }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            face_val(pelem->n[1],pelem->n[2],pelem->n[3],pelem->f[0],u,v,w);
            face_val(pelem->n[0],pelem->n[2],pelem->n[3],pelem->f[1],u,v,w);
            face_val(pelem->n[0],pelem->n[1],pelem->n[3],pelem->f[2],u,v,w);
            face_val(pelem->n[0],pelem->n[1],pelem->n[2],pelem->f[3],u,v,w);
         }
   node_correction(mg,u,v,w);
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pface=FDBF(theGrid); pface != FIRSTFACE(theGrid); pface=pface->succ)
         if (IS_TOP_FACE(pface))
            FD(pface,u) = 0.;
}

#else

void face_val(n1,n2,n3,f0,u,v,w)
NODE *n1, *n2, *n3; FACE *f0; INT u, v, w;
{  eprintf("Error: face_val not available.\n");  }

void bound_val(mg,u,v,w)
MULTIGRID *mg; INT u, v, w;
{  eprintf("Error: bound_val not available.\n");  }

#endif


#endif

#endif

/* -------------------------------------------------------------------------- */

#if DIM == 2

 /* s1 is a linear approximation of (\int_K' u0\cdot nn \ds)/K'; s2 is a better 
    (quadratic) approximation of this integral. */
 void integr2(x1,x2,nn,s1,s2)
 FLOAT *x1, *x2, *nn, *s1, *s2;
 {
   FLOAT xs[DIM];
   
   AVERAGE(x1,x2,xs);
   
   *s1 = ( nn[0]*(u01(x1) + u01(x2))+
           nn[1]*(u02(x1) + u02(x2)) )/2.0;
          
   *s2 = ( nn[0]*(u01(x1) + 4.*u01(xs) + u01(x2))+
           nn[1]*(u02(x1) + 4.*u02(xs) + u02(x2)) )/6.0;
 }
 
#if PROBLEM & BUBBLE_BC

#if !(DATA_STR & LG_DATA)

#if F_DATA & SCALAR_FACE_DATA

 /* x1, x2 are vertexes of the i-th face, x0 is the remaining vertex */
 void face_val(pelement,i,x1,x2,x0,u,sum1,sum2)
 ELEMENT *pelement;
 INT i, u;
 FLOAT  *x1, *x2, *x0, *sum1, *sum2;  
 {                              /* sum1 = \int_{\partial\om} p\cdot\nn \ds  */
    FLOAT nn[DIM], k0, k1, area;/* sum2 ~ \int_{\partial\om} u0\cdot\nn \ds */
    
    if (NOT_FF(pelement->f[i])){
       area = normal_vector(x1,x2,pelement->f[i],nn); 
       integr2(x1,x2,nn,&k0,&k1);  
       FD(pelement->f[i],u) = 6.0*(k1-k0)/FMULT;
       test_normal_vector_direction(nn,x0,x1,&area);
       *sum1 += area*FD(pelement->f[i],u);
       *sum2 += area*k1;
    }
 }

#else

 void face_val(pelement,i,x1,x2,x0,u,sum1,sum2)
 ELEMENT *pelement; INT i, u; FLOAT  *x1, *x2, *x0, *sum1, *sum2;  
 {  eprintf("Error: face_val not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

 void bound_val(mg,u,v,w) /* v, w not used */
 MULTIGRID *mg;
 INT u, v, w;
 {
   GRID *theGrid;
   NODE *pnode;
   ELEMENT *pelem;
   FLOAT *x0, *x1, *x2, sum1=0.0, sum2=0.0;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer){
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u01(pnode->myvertex->x);
            ND(pnode,u,1) = u02(pnode->myvertex->x);
         }
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            face_val(pelem,0,x1,x2,x0,u,&sum1,&sum2);
            face_val(pelem,1,x0,x2,x1,u,&sum1,&sum2);
            face_val(pelem,2,x0,x1,x2,u,&sum1,&sum2);
         }
   }
   sum1 = sum1/6.0*FMULT;
   bubble_correction(mg,FIRSTGRID(mg),sum1,sum2,u);
 }

#else

 void bound_val(mg,u,v,w)
 MULTIGRID *mg; INT u, v, w;
 {  eprintf("Error: bound_val not available.\n");  }

#endif

#endif

#if DATA_STR & LG_DATA

void bound_val(mg,u,v,w) 
MULTIGRID *mg;
INT u, v, w;
{  eprintf("Error: bound_val not available.\n");  }

#endif
#endif

#if !(PROBLEM & BUBBLE_BC)

#if N_DATA & VECTOR_NODE_DATA

void face_val(n1,n2,f0,u,v,w)
NODE *n1, *n2;
FACE *f0;
INT u, v, w;
{
   FLOAT nn[DIM], area, s, q, *x1, *x2, x112[DIM], x221[DIM];
    
   if (NOT_FF(f0)){
      x1 = n1->myvertex->x;
      x2 = n2->myvertex->x;
      area = normal_vector(x1,x2,f0,nn);
      POINT2(x1,x2,x112,2.0,3.0);
      POINT2(x2,x1,x221,2.0,3.0);
      s = 2.*(nn[0]*(ND(n1,u,0)+ND(n2,u,0)) + nn[1]*(ND(n1,u,1)+ND(n2,u,1)));
      q = area*3./40.;
      ND(n1,w,0) += q*( nn[0]*( 4.*u01(x112) + u01(x221) - ND(n1,u,0) ) + 
                        nn[1]*( 4.*u02(x112) + u02(x221) - ND(n1,u,1) ) - s );
      ND(n2,w,0) += q*( nn[0]*( 4.*u01(x221) + u01(x112) - ND(n2,u,0) ) + 
                        nn[1]*( 4.*u02(x221) + u02(x112) - ND(n2,u,1) ) - s );
      q = area/2.;
      SET4(NDD(n1,v),nn,q)
      SET4(NDD(n2,v),nn,q)
   }
}

void P1_bound_val(mg,u,v,w,u01,u02) /*  Dirichlet b.c. has to be prescribed  */
MULTIGRID *mg;                      /*                on the whole boundary  */
FLOAT (*u01)(), (*u02)();
INT u, v, w;
{
   GRID *theGrid;
   NODE *pnode;
   ELEMENT *pelem;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FDBN(theGrid); pnode !=FIRSTNODE(theGrid); pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u01(pnode->myvertex->x);
            ND(pnode,u,1) = u02(pnode->myvertex->x);
            ND(pnode,v,0) = ND(pnode,v,1) = ND(pnode,w,0) = 0.;
         }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            face_val(pelem->n[1],pelem->n[2],pelem->f[0],u,v,w);
            face_val(pelem->n[0],pelem->n[2],pelem->f[1],u,v,w);
            face_val(pelem->n[0],pelem->n[1],pelem->f[2],u,v,w);
         }
   node_correction(mg,u,v,w);
}

#else

void face_val(n1,n2,f0,u,v,w)
NODE *n1, *n2; FACE *f0; INT u, v, w;
{  eprintf("Error: face_val not available.\n");  }

void P1_bound_val(mg,u,v,w,u01,u02)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u, v, w;
{  eprintf("Error: P1_bound_val not available.\n");  }

#endif

void bound_val(mg,u,v,w) /* Dirichlet b.c. has to be prescribed on the whole  */
MULTIGRID *mg;           /*                                          boundary */
INT u, v, w;
{
   P1_bound_val(mg,u,v,w,u01,u02);
   set_value(TOP_GRID(mg),0.,u,STOP_IS_FIRST_INNER,Q_SF);
}

#endif

#else  /*  DIM != 2  */

void P1_bound_val(mg,u,v,w,u01,u02)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u, v, w;
{  eprintf("Error: P1_bound_val not available.\n");  }

#endif

/******************************************************************************/
/*                                                                            */
/*                     nonconforming boundary conditions                      */
/*                                                                            */
/******************************************************************************/

#if DIM == 2

#if F_DATA & SCALAR_FACE_DATA

void set_bound_edge_value(f12,x1,x2,u0,u)
FACE *f12;
FLOAT *x1, *x2, (*u0)();
INT u;
{
   FLOAT x12[DIM]; 

   if (NOT_FF(f12)){
      AVERAGE(x1,x2,x12);
      FD(f12,u) = u0(x12);
   }
}

#else

void set_bound_edge_value(f12,x1,x2,u0,u)
FACE *f12; FLOAT *x1, *x2, (*u0)(); INT u;
{  eprintf("Error: set_bound_edge_value not available.\n");  }

#endif

#if ELEMENT_TYPE == SIMPLEX

void face_sbound_val(mg,u0,u)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         VERTICES_OF_ELEMENT(x1,x2,x3,pelem);
         set_bound_edge_value(f1,x2,x3,u0,u);
         set_bound_edge_value(f2,x3,x1,u0,u);
         set_bound_edge_value(f3,x1,x2,u0,u);
      }
}     

#else

void face_sbound_val(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_sbound_val not available.\n");  }

#endif

#if F_DATA & VECTOR_FACE_DATA

void vset_bound_edge_value(f12,x1,x2,u01,u02,u)
FACE *f12;
FLOAT *x1, *x2, (*u01)(), (*u02)();
INT u;
{
   FLOAT x12[DIM]; 

   if (NOT_FF(f12)){
/*
      AVERAGE(x1,x2,x12);
      FDV(f12,u,0) = (u01(x1) + 4.*u01(x12) + u01(x2))/6.;
      FDV(f12,u,1) = (u02(x1) + 4.*u02(x12) + u02(x2))/6.;
*/
      FDV(f12,u,0) = integrate_using_multiple_Simpson(x1,x2,u01,500);
      FDV(f12,u,1) = integrate_using_multiple_Simpson(x1,x2,u02,500);
   }
}

#else

void vset_bound_edge_value(f12,x1,x2,u01,u02,u)
FACE *f12; FLOAT *x1, *x2, (*u01)(), (*u02)(); INT u;
{  eprintf("Error: vset_bound_edge_value not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX)

FLOAT nc_l_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return(DOT(FDVP(f3,u),nn)*area);
   }
   else
      return(0.);
}

FLOAT tnc_l_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return(FDV(f3,u,0)*area);
   }
   else
      return(0.);
}

void tvset_edge_value(f12,x1,x2,u01,u02,u)
FACE *f12;
FLOAT *x1, *x2, (*u01)(), (*u02)();
INT u;
{
   FLOAT x12[DIM], nn[DIM], tt[DIM], uu[DIM]; 

   if (IS_FF(f12)){
      normal_vector(x1,x2,f12,nn);
      tt[0] = -nn[1];
      tt[1] =  nn[0];
      AVERAGE(x1,x2,x12);
      uu[0] = u01(x12);
      uu[1] = u02(x12);
      FDV(f12,u,0) = DOT(uu,nn);
      FDV(f12,u,1) = DOT(uu,tt);
   }
}

void tvset_bound_edge_value(f12,x1,x2,u01,u02,u)
FACE *f12;
FLOAT *x1, *x2, (*u01)(), (*u02)();
INT u;
{
   FLOAT x12[DIM], nn[DIM], tt[DIM], uu[DIM];

   if (NOT_FF(f12)){
      normal_vector(x1,x2,f12,nn);
      tt[0] = -nn[1];
      tt[1] =  nn[0];
      AVERAGE(x1,x2,x12);
      uu[0] = (u01(x1) + 4.*u01(x12) + u01(x2))/6.;
      uu[1] = (u02(x1) + 4.*u02(x12) + u02(x2))/6.;
      FDV(f12,u,0) = DOT(uu,nn);
      FDV(f12,u,1) = DOT(uu,tt);
   }
}

void tvface_init(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         VERTICES_OF_ELEMENT(x1,x2,x3,pelem);
         tvset_edge_value(f1,x2,x3,u01,u02,u);
         tvset_edge_value(f2,x1,x3,u01,u02,u);
         tvset_edge_value(f3,x1,x2,u01,u02,u);
      }
}

void vface_bound_val(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         VERTICES_OF_ELEMENT(x1,x2,x3,pelem);
         vset_bound_edge_value(f1,x2,x3,u01,u02,u);
         vset_bound_edge_value(f2,x3,x1,u01,u02,u);
         vset_bound_edge_value(f3,x1,x2,u01,u02,u);
      }
}     

void tvface_bound_val(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         VERTICES_OF_ELEMENT(x1,x2,x3,pelem);
         tvset_bound_edge_value(f1,x2,x3,u01,u02,u);
         tvset_bound_edge_value(f2,x1,x3,u01,u02,u);
         tvset_bound_edge_value(f3,x1,x2,u01,u02,u);
      }
      printf("Flux of lin. int. of u0 through the boundary: %e\n",
                                           boundary_flux(mg,u,tnc_l_face_flux));
}     

void set_bound_edge_values_mbub(f12,n1,n2,u0,u)
FACE *f12;
NODE *n1, *n2;
FLOAT (*u0)();
INT u;
{
   FLOAT v1, v2; 

   if (NOT_FF(f12)){
      v1 = u0(n1->myvertex->x);
      v2 = u0(n2->myvertex->x);
      FDV(f12,u,0) = 0.5*(v1 + v2);
      FDV(f12,u,1) =  5.*(v1 - v2)*NINDI(n2,n1);
   }
}

void face_sbound_val_mbub(mg,u0,u)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         set_bound_edge_values_mbub(f1,n2,n3,u0,u);
         set_bound_edge_values_mbub(f2,n3,n1,u0,u);
         set_bound_edge_values_mbub(f3,n1,n2,u0,u);
      }
}     

#else  /*  !(F_DATA & VECTOR_FACE_DATA)  */

FLOAT nc_l_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: nc_l_face_flux not available.\n"); return(0.);  }

void vface_bound_val(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: vface_bound_val not available.\n");  }

void face_sbound_val_mbub(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_sbound_val_mbub not available.\n");  }

#endif

#else  /*  DIM != 2  */

FLOAT nc_l_face_flux(n1,n2,f3,u)
NODE *n1, *n2; FACE *f3; INT u;
{  eprintf("Error: nc_l_face_flux not available.\n"); return(0.);  }

void vface_bound_val(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: vface_bound_val not available.\n");  }

void face_sbound_val(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_sbound_val not available.\n");  }

void face_sbound_val_mbub(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_sbound_val_mbub not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

void vface_bound_val_q1rot(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f0, *f1, *f2, *f3;
   FLOAT *x0, *x1, *x2, *x3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_4ELEMENT(f0,f1,f2,f3,pelem);
         VERTICES_OF_4ELEMENT(x0,x1,x2,x3,pelem);
         vset_bound_edge_value(f0,x0,x1,u01,u02,u);
         vset_bound_edge_value(f1,x1,x2,u01,u02,u);
         vset_bound_edge_value(f2,x2,x3,u01,u02,u);
         vset_bound_edge_value(f3,x3,x0,u01,u02,u);
      }
}     

#else

void vface_bound_val_q1rot(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: vface_bound_val_q1rot not available.\n");  }

#endif

#if (F_DATA & DVECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

FLOAT dv_nc_l_face_flux(n1,n2,f3,u)
NODE *n1, *n2;
FACE *f3;
INT u;
{
   FLOAT nn[DIM], area;

   if (NOT_FF(f3)){
      area = normal_vector(n1->myvertex->x,n2->myvertex->x,f3,nn);
      return(DOT(FDDVP(f3,u,0),nn)*area);
   }
   else
      return(0.);
}

void dvset_bound_edge_values_mbub(f12,n1,n2,u01,u02,u)
FACE *f12;
NODE *n1, *n2;
FLOAT (*u01)(), (*u02)();
INT u;
{
   FLOAT v1, v2; 

   if (NOT_FF(f12)){
      v1 = u01(n1->myvertex->x);
      v2 = u01(n2->myvertex->x);
      FDDV(f12,u,0,0) = 0.5*(v1 + v2);
      FDDV(f12,u,1,0) =  5.*(v1 - v2)*NINDI(n2,n1);
      v1 = u02(n1->myvertex->x);
      v2 = u02(n2->myvertex->x);
      FDDV(f12,u,0,1) = 0.5*(v1 + v2);
      FDDV(f12,u,1,1) =  5.*(v1 - v2)*NINDI(n2,n1);
   }
}

void sdvset_bound_edge_values_mbub(f12,n1,n2,u1,u2,u)
FACE *f12;
NODE *n1, *n2;
FLOAT u1[DIM], u2[DIM];
INT u;
{
   if (NOT_FF(f12)){
      FDDV(f12,u,0,0) = 0.5*(u1[0] + u2[0]);
      FDDV(f12,u,1,0) =  5.*(u1[0] - u2[0])*NINDI(n2,n1);
      FDDV(f12,u,0,1) = 0.5*(u1[1] + u2[1]);
      FDDV(f12,u,1,1) =  5.*(u1[1] - u2[1])*NINDI(n2,n1);
   }
}

void dvface_bound_val_mbub(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         dvset_bound_edge_values_mbub(f1,n2,n3,u01,u02,u);
         dvset_bound_edge_values_mbub(f2,n3,n1,u01,u02,u);
         dvset_bound_edge_values_mbub(f3,n1,n2,u01,u02,u);
      }
}     

void correct_dvface_bound_val_mbub(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   if (boundary_flux(mg,u,dv_nc_l_face_flux) > 1.e-6){
      eprintf("*****************************************************\n");
      eprintf("*****************************************************\n");
      eprintf("***                                               ***\n");
      eprintf("*** boundary condition not approximated correctly ***\n");
      eprintf("***                                               ***\n");
      eprintf("*****************************************************\n");
      eprintf("*****************************************************\n");
/*
      bound_val(mg,U,F,D);
      for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
         for (pelem = FDBE(theGrid); pelem != NULL; pelem = pelem->succ){
            dvset_bound_edge_values_mbub(f1,n2,n3,NDD(n2,U),NDD(n3,U),u);
            dvset_bound_edge_values_mbub(f2,n3,n1,NDD(n3,U),NDD(n1,U),u);
            dvset_bound_edge_values_mbub(f3,n1,n2,NDD(n1,U),NDD(n2,U),u);
         }
*/
   }
   printf("Flux of lin. int. of u0 through the boundary: %e\n",
                                         boundary_flux(mg,u,dv_nc_l_face_flux));
}     

#else  /*  !(F_DATA & DVECTOR_FACE_DATA)  */

void dvface_bound_val_mbub(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: dvface_bound_val_mbub not available.\n");  }

void correct_dvface_bound_val_mbub(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: correct_dvface_bound_val_mbub not available.\n");  }

#endif

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (DIM == 2)

void general_bound_val(mg,u,u0,rm_type,finite_el)
MULTIGRID *mg;
FLOAT (*u0)();
INT u, rm_type;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   DOUBLE_FUNC *nd, *fd;
   REF_MAPPING ref_map;
   INT i, k, m, kk, nn, nv=NVERT-1, dir[4];

   nd = finite_el.ndofs;
   fd = finite_el.fdofs;
   kk = finite_el.k;
   nn = finite_el.n;
   for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel) && (NOT_FF(pel->f[0]) || NOT_FF(pel->f[1]) 
                                  || NOT_FF(pel->f[2]) || NOT_FF(pel->f[nv]))){
            reference_mapping(pel,rm_type,&ref_map);
            if (kk)
               for (i=0; i < NVERT; i++){
                  m = i*kk;
                  if (NOT_FN(pel->n[i]))
                     for (k = m; k < m+kk; k++) 
                        NDMV(pel->n[i],u,k-m) = (nd[k])(u0,&ref_map);
               }
            if (nn){
               set_directions_of_edges(pel,dir);
               for (i=0; i < NVERT; i++){
                  m = i*nn;
                  if (NOT_FF(pel->f[i]))
                     for (k = m; k < m+nn; k++) 
                        FDMV(pel->f[i],u,f_i(k-m,nn,dir[i])) = 
                                                           (fd[k])(u0,&ref_map);
               }
            }
         }
}

#else

void general_bound_val(mg,u,u0,rm_type,finite_el)
MULTIGRID *mg; FLOAT (*u0)(); INT u, rm_type; FINITE_ELEMENT finite_el;
{  eprintf("Error: general_bound_val not available.\n");  }

#endif


/******************************************************************************/

void boundary_values(mg,u01,u02,u03,u,v,w,space,structure)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u, v, w, space;                  /*  v, w ... auxiliary variables  */
{
   switch(space){
   case P1C:       if (structure == SCALAR)
                      sbound_val1(mg,u,u01);
                   else if (structure == VECTOR)
                      vbound_val1(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case P1C_ELBUB: if (structure == SCALAR)
                      sbound_val1(mg,u,u01);
                   else if (structure == VECTOR)
                      vbound_val1(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case P1_NC:     if (structure == SCALAR)
                      face_sbound_val(mg,u01,u);
                   else if (structure == VECTOR)
                      vface_bound_val(mg,u01,u02,u);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case P1_MOD:    if (structure == SCALAR)
                      face_sbound_val_mbub(mg,u01,u);
                   else if (structure == VECTOR)
                      dvface_bound_val_mbub(mg,u01,u02,u);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case MINI_L_DIV_FR: if (structure == SCALAR)
                      sbound_val1(mg,u,u01);
                   else if (structure == VECTOR)
                      vbound_val1(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case Q1C:       if (structure == SCALAR)
                      sbound_val1(mg,u,u01);
                   else if (structure == VECTOR)
                      vbound_val1(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case P2C:       if (structure == SCALAR)
                      sbound_val2(mg,u,u01);
                   else if (structure == VECTOR && DIM == 2)
                      bound_val_vn_vf_2(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case IP2C:      if (structure == SCALAR)
                      sbound_val2_iso(mg,u,u01);
                   else if (structure == VECTOR && DIM == 2)
                      bound_val_vn_vf_2_iso(mg,u,u01,u02);
                   else
                      eprintf("Error: boundary_values not available.\n");
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
                   if (structure == SCALAR)
                      general_bound_val(mg,u,u01,REF_MAP,ELEM);
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   case GQ2B3C:    if (structure == SCALAR){
                      mv_set_value_n(TOP_GRID(mg),0.,u,STOP_IS_FIRST_INNER);
                      mv_set_value_f(TOP_GRID(mg),0.,u,STOP_IS_FIRST_INNER);
                      general_bound_val(mg,u,u01,REF_MAP,q2c_element);
                   }
                   else
                      eprintf("Error: boundary_values not available.\n");
        break;
   default:
        eprintf("Error: boundary_values not available.\n");
        break;
   }
}

void boundary_values_of_velocity(mg,u01,u02,u03,u,v,w,space)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)(), (*u03)();
INT u, v, w, space;                  /*  v, w ... auxiliary variables  */
{
   switch(space){
   case P1_NC:     vface_bound_val(mg,u01,u02,u);
                   printf("Flux of lin. int. of u0 through the boundary: %e\n",
                                            boundary_flux(mg,u,nc_l_face_flux));
        break;
   case P1_MOD:    dvface_bound_val_mbub(mg,u01,u02,u);
                   correct_dvface_bound_val_mbub(mg,u01,u02,u);
        break;
   case P1C:       P1_bound_val(mg,u,v,w,u01,u02);
        break;
   case P1C_FBUB:  if (!(PROBLEM & BUBBLE_BC)){
                      P1_bound_val(mg,u,v,w,u01,u02);
                      set_value(TOP_GRID(mg),0.,u,STOP_IS_FIRST_INNER,Q_SF);
                   }
                   else
                      eprintf("Error: boundary_values_of_velocity not available.\n");
        break;
   case P1C_NEW_FBUB:  if (!(PROBLEM & BUBBLE_BC)){
                      P1_bound_val(mg,u,v,w,u01,u02);
                      set_value(TOP_GRID(mg),0.,u,STOP_IS_FIRST_INNER,Q_SF);
                   }
                   else
                      eprintf("Error: boundary_values_of_velocity not available.\n");
        break;
   case P1C_ELBUB: P1_bound_val(mg,u,v,w,u01,u02);
        break;
   case Q1ROT:     vface_bound_val_q1rot(mg,u01,u02,u);
        break;
   case P2C:       bound_val_vn_vf_2(mg,u,u01,u02);
                   printf("Flux of the discrete b.c.: %e\n",
                                             boundary_flux(mg,u,vlb_face_flux));
                   correct_boundary_values_2_5(mg,u);
                   printf("Flux of the discrete b.c.: %e\n",
                                             boundary_flux(mg,u,vlb_face_flux));
        break;
   case IP2C:      bound_val_vn_vf_2_iso(mg,u,u01,u02);
        break;
   case P2C_ELBUB: bound_val_vn_vf_2(mg,u,u01,u02);
                   printf("Flux of the discrete b.c.: %e\n",
                                             boundary_flux(mg,u,vlb_face_flux));
                   correct_boundary_values_2_5(mg,u);
                   printf("Flux of the discrete b.c.: %e\n",
                                             boundary_flux(mg,u,vlb_face_flux));
        break;
   case IP2C_ELBUB: bound_val_vn_vf_2(mg,u,u01,u02);
        break;
   default:
        eprintf("Error: boundary_values_of_velocity not available.\n");
        break;
   }
}

