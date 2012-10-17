/******************************************************************************/
/*                                                                            */
/*                               postprocessing                               */
/*                                                                            */
/******************************************************************************/

INT cmp_vertices(pel,a0,a1,b0,b1,c0,c1,j0,j1,j2,i0,i1,i2)
ELEMENT *pel;
FLOAT a0, a1, b0, b1, c0, c1;
INT j0, j1, j2, *i0, *i1, *i2;
{
   FLOAT *x0, *x1, *x2, eps=1.e-12;

   x0 = pel->n[j0]->myvertex->x;
   x1 = pel->n[j1]->myvertex->x;
   x2 = pel->n[j2]->myvertex->x;
   if (fabs(x0[0]-a0) < eps && fabs(x0[1]-a1) < eps &&
       fabs(x1[0]-b0) < eps && fabs(x1[1]-b1) < eps &&
       fabs(x2[0]-c0) < eps && fabs(x2[1]-c1) < eps){
         *i0 = j0;
         *i1 = j1;
         *i2 = j2;
         return(1);
   }
   return(0);
}
       
INT has_vertices(pel,a0,a1,b0,b1,c0,c1,i0,i1,i2)
ELEMENT *pel;
FLOAT a0, a1, b0, b1, c0, c1;
INT *i0, *i1, *i2;
{
   if (cmp_vertices(pel,a0,a1,b0,b1,c0,c1,0,1,2,i0,i1,i2) ||
       cmp_vertices(pel,a0,a1,b0,b1,c0,c1,1,2,0,i0,i1,i2) ||
       cmp_vertices(pel,a0,a1,b0,b1,c0,c1,2,0,1,i0,i1,i2) ||
       cmp_vertices(pel,a0,a1,b0,b1,c0,c1,0,2,1,i0,i1,i2) ||
       cmp_vertices(pel,a0,a1,b0,b1,c0,c1,1,0,2,i0,i1,i2) ||
       cmp_vertices(pel,a0,a1,b0,b1,c0,c1,2,1,0,i0,i1,i2))
      return(1);
   else
      return(0);
}

INT find_element_with_vertices(tGrid,a0,a1,b0,b1,c0,c1,i0,i1,i2,pelem)
GRID *tGrid;
ELEMENT **pelem;
FLOAT a0, a1, b0, b1, c0, c1;
INT *i0, *i1, *i2;
{
   INT i=0;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel && !i; pel=pel->succ)
      if(i=has_vertices(pel,a0,a1,b0,b1,c0,c1,i0,i1,i2))
         *pelem = pel;
   return(i);
}

double azimuthal_coordinate(x)
double *x;
{
   double phi;

   if (fabs(x[0]) < 1.e-10 && fabs(x[1]) < 1.e-10)
      return(0.);
   else if (fabs(x[0]) < 1.e-10){
      if (x[1] > 0.)
         return(0.5*PI);
      else
         return(1.5*PI);
   }
   else{
      phi = atan(x[1]/x[0]);
      if (x[0] > 0. && x[1] < 0.)
         return(2.*PI + phi);
      else if (x[0] < 0.)
         return(PI + phi);
      else
         return(phi);
   }
}

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void edge_value_to_vertices(f1,n2,n3,u,v,w)
FACE *f1;
NODE *n2, *n3;
INT u, v, w;
{
   if (IS_FF(f1)){
      NDS(n2,v) += FD(f1,u);
      NDS(n3,v) += FD(f1,u);
      NDS(n2,w) += 1.;
      NDS(n3,w) += 1.;
   }
   else{
      NDS(n2,v) += 2.*FD(f1,u);
      NDS(n3,v) += 2.*FD(f1,u);
      NDS(n2,w) += 2.;
      NDS(n3,w) += 2.;
   }
}

#else  /*  !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void edge_value_to_vertices(f1,n2,n3,u,v,w)
FACE *f1; NODE *n2, *n3; INT u, v, w;
{  eprintf("Error: edge_value_to_vertices not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (DIM == 2)

void nc_conf_averaging(tGrid,u,v,w,t0)
GRID *tGrid;
INT u, v, w;
FLOAT (*t0)();
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      NDS(pnode,v) = NDS(pnode,w) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      edge_value_to_vertices(f1,n2,n3,u,v,w);
      edge_value_to_vertices(f2,n3,n1,u,v,w);
      edge_value_to_vertices(f3,n1,n2,u,v,w);
   }
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      NDS(pnode,v) /= NDS(pnode,w);
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      if (IS_DN(pnode))
         NDS(pnode,v) = t0(pnode->myvertex->x);
}

void make_averages(tGrid,u_old,u_new)
GRID *tGrid;
INT u_old, u_new;
{
   NODE *pnode;
   LINK *pli;
   FLOAT sum;
   INT n;
  
   for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ){
      sum = 0.;
      n = 0;
      for (pli = START(pnode); pli; pli = pli->next){
         n++;
         sum += NDS(NBNODE(pli),u_old);
      }
      NDS(pnode,u_new) = 0.5*(NDS(pnode,u_old) + sum/n); 
   }
}

#else  /*  !(N_DATA & VECTOR_NODE_DATA)  */

void nc_conf_averaging(tGrid,u,v,w,t0)
GRID *tGrid; INT u, v, w; FLOAT (*t0)();
{  eprintf("Error: nc_conf_averaging not available.\n");  } 

void make_averages(tGrid,u_old,u_new)
GRID *tGrid; INT u_old, u_new;
{  eprintf("Error: make_averages not available.\n");  } 

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void v_edge_value_to_vertices(f1,n2,n3,u,v)
FACE *f1;
NODE *n2, *n3;
INT u, v;
{
   if (IS_FF(f1)){
      ND(n2,u,0) += FDV(f1,u,0);
      ND(n3,u,0) += FDV(f1,u,0);
      ND(n2,u,1) += FDV(f1,u,1);
      ND(n3,u,1) += FDV(f1,u,1);
      ND(n2,v,1) += 1.;
      ND(n3,v,1) += 1.;
   }
   else{
      ND(n2,u,0) += 2.*FDV(f1,u,0);
      ND(n3,u,0) += 2.*FDV(f1,u,0);
      ND(n2,u,1) += 2.*FDV(f1,u,1);
      ND(n3,u,1) += 2.*FDV(f1,u,1);
      ND(n2,v,1) += 2.;
      ND(n3,v,1) += 2.;
   }
}

void v_nc_conf_averaging(tGrid,u,v)
GRID *tGrid;
INT u, v;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      ND(pnode,u,0) = ND(pnode,u,1) = ND(pnode,v,1) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      v_edge_value_to_vertices(f1,n2,n3,u,v);
      v_edge_value_to_vertices(f2,n3,n1,u,v);
      v_edge_value_to_vertices(f3,n1,n2,u,v);
   }
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ){
      ND(pnode,u,0) /= ND(pnode,v,1);
      ND(pnode,u,1) /= ND(pnode,v,1);
   }
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      if (IS_DN(pnode)){
         ND(pnode,u,0) = u01(pnode->myvertex->x);
         ND(pnode,u,1) = u02(pnode->myvertex->x);
      }
}

#else

void v_edge_value_to_vertices(f1,n2,n3,u,v)
FACE *f1; NODE *n2, *n3; INT u, v;
{  eprintf("Error: v_edge_value_to_vertices not available.\n");  }

void v_nc_conf_averaging(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: v_nc_conf_averaging not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void dv_edge_value_to_vertices(f1,n2,n3,u,v)
FACE *f1;
NODE *n2, *n3;
INT u, v;
{
   if (IS_FF(f1)){
      ND(n2,u,0) += FDDV(f1,u,0,0);
      ND(n3,u,0) += FDDV(f1,u,0,0);
      ND(n2,u,1) += FDDV(f1,u,0,1);
      ND(n3,u,1) += FDDV(f1,u,0,1);
      ND(n2,v,1) += 1.;
      ND(n3,v,1) += 1.;
   }
   else{
      ND(n2,u,0) += 2.*FDDV(f1,u,0,0);
      ND(n3,u,0) += 2.*FDDV(f1,u,0,0);
      ND(n2,u,1) += 2.*FDDV(f1,u,0,1);
      ND(n3,u,1) += 2.*FDDV(f1,u,0,1);
      ND(n2,v,1) += 2.;
      ND(n3,v,1) += 2.;
   }
}

void dv_nc_conf_averaging(tGrid,u,v)
GRID *tGrid;
INT u, v;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *pnode;
   FACE *f1, *f2, *f3;

   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      ND(pnode,u,0) = ND(pnode,u,1) = ND(pnode,v,1) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      dv_edge_value_to_vertices(f1,n2,n3,u,v);
      dv_edge_value_to_vertices(f2,n3,n1,u,v);
      dv_edge_value_to_vertices(f3,n1,n2,u,v);
   }
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ){
      ND(pnode,u,0) /= ND(pnode,v,1);
      ND(pnode,u,1) /= ND(pnode,v,1);
   }
   for (pnode = FDBN(tGrid); pnode != NULL; pnode = pnode->succ)
      if (IS_DN(pnode)){
         ND(pnode,u,0) = u01(pnode->myvertex->x);
         ND(pnode,u,1) = u02(pnode->myvertex->x);
      }
}

#else

void dv_edge_value_to_vertices(f1,n2,n3,u,v)
FACE *f1; NODE *n2, *n3; INT u, v;
{  eprintf("Error: dv_edge_value_to_vertices not available.\n");  }

void dv_nc_conf_averaging(tGrid,u,v)
GRID *tGrid; INT u, v;
{  eprintf("Error: dv_nc_conf_averaging not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void p1nc_update_drag_and_lift(pel,n0,n1,n2,fa0,fa1,fa2,u,lift,drag)
ELEMENT *pel;
NODE *n0, *n1, *n2;   
FACE *fa0, *fa1, *fa2;  /* fa0 lies on the circle  */
INT u;
FLOAT *lift, *drag;
{
   FLOAT nn[DIM], d, s, gu0_0, gu0_1, gu1_0, gu1_1, p, b[DIM2][DIM2];

   barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,b);
   d = normal_vector(n1->myvertex->x,n2->myvertex->x,fa0,nn);
   gu0_0 = -2.*(FDV(fa0,u,0)*b[0][0]+FDV(fa1,u,0)*b[1][0]+FDV(fa2,u,0)*b[2][0]);
   gu0_1 = -2.*(FDV(fa0,u,0)*b[0][1]+FDV(fa1,u,0)*b[1][1]+FDV(fa2,u,0)*b[2][1]);
   gu1_0 = -2.*(FDV(fa0,u,1)*b[0][0]+FDV(fa1,u,1)*b[1][0]+FDV(fa2,u,1)*b[2][0]);
   gu1_1 = -2.*(FDV(fa0,u,1)*b[0][1]+FDV(fa1,u,1)*b[1][1]+FDV(fa2,u,1)*b[2][1]);
   p = ED(pel,U);
/*   s = gu0_1-gu1_0; */
   s = (gu0_0*nn[1] - gu1_0*nn[0])*nn[0] + (gu0_1*nn[1] - gu1_1*nn[0])*nn[1];
   *lift += d*(NU*s*nn[0]+p*nn[1]);
   *drag += d*(NU*s*nn[1]-p*nn[0]);
}

#else

void p1nc_update_drag_and_lift(pel,n0,n1,n2,fa0,fa1,fa2,u,lift,drag)
ELEMENT *pel; NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; INT u; FLOAT *lift, *drag;
{  eprintf("Error: p1nc_update_drag_and_lift not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (N_DATA & SCALAR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void p2c_update_drag_and_lift(n0,n1,n2,fa0,fa1,fa2,u,lift,drag)
NODE *n0, *n1, *n2;   
FACE *fa0, *fa1, *fa2;  /* fa0 lies on the circle  */
INT u;
FLOAT *lift, *drag;
{
   FLOAT nn[DIM], d, s, gu0_0, gu0_1, gu1_0, gu1_1, p, b[DIM2][DIM2];

   barycentric_coordinates(n0->myvertex->x,n1->myvertex->x,n2->myvertex->x,b);
   d = normal_vector(n1->myvertex->x,n2->myvertex->x,fa0,nn);
   gu0_0 = ND(n0,u,0)*b[0][0]+ND(n1,u,0)*b[1][0]+ND(n2,u,0)*b[2][0]+
           0.5*(FDV(fa1,u,0)+FDV(fa2,u,0)-FDV(fa0,u,0))*b[0][0];
   gu0_1 = ND(n0,u,0)*b[0][1]+ND(n1,u,0)*b[1][1]+ND(n2,u,0)*b[2][1]+
           0.5*(FDV(fa1,u,0)+FDV(fa2,u,0)-FDV(fa0,u,0))*b[0][1];
   gu1_0 = ND(n0,u,1)*b[0][0]+ND(n1,u,1)*b[1][0]+ND(n2,u,1)*b[2][0]+
           0.5*(FDV(fa1,u,1)+FDV(fa2,u,1)-FDV(fa0,u,1))*b[0][0];
   gu1_1 = ND(n0,u,1)*b[0][1]+ND(n1,u,1)*b[1][1]+ND(n2,u,1)*b[2][1]+
           0.5*(FDV(fa1,u,1)+FDV(fa2,u,1)-FDV(fa0,u,1))*b[0][1];
   p = 0.5*(NDS(n1,u)+NDS(n2,u));
/*   s = gu0_1-gu1_0;  */
   s = (gu0_0*nn[1] - gu1_0*nn[0])*nn[0] + (gu0_1*nn[1] - gu1_1*nn[0])*nn[1];
   *lift += d*(NU*s*nn[0]+p*nn[1]);
   *drag += d*(NU*s*nn[1]-p*nn[0]);
}

void define_v(tGrid,x,i)
GRID *tGrid;
INT x,i;
{
   NODE *theNode;

   set_value(tGrid,0.,x,0,Q_VNVF);
   for (theNode=FIRSTN(tGrid); theNode!=FIRSTNODE(tGrid); theNode=SUCC(theNode))
   if (theNode->myvertex->x[0] > 0.05 && theNode->myvertex->x[0] < 0.5 &&
       theNode->myvertex->x[1] > 0.1  && theNode->myvertex->x[1] < 0.3){
/*
      printf("%e\n",(theNode->myvertex->x[0]-0.2)*(theNode->myvertex->x[0]-0.2)+
        (theNode->myvertex->x[1]-0.2)*(theNode->myvertex->x[1]-0.2));
*/
      ND(theNode,x,i) = 1.;
   }
}

void lift_and_drag_for_p2c_p1c_2(tGrid,ZA,ZB,u,d,r,q,t_for_p,u_type,p_type,
                                 a_struct,b_struct)
GRID *tGrid;
INT ZA,ZB,u,d,r,q,t_for_p,u_type,p_type,a_struct,b_struct;
{
   set_value(tGrid,0.,r,0,u_type);
   Stokes_defect(tGrid,ZA,ZB,r,r,u,u,d,d,q,
                                   0,t_for_p,u_type,p_type,a_struct,b_struct,0);
   define_v(tGrid,q,0);
   printf("drag = %e\n",-500*dot(tGrid,d,q,0,u_type));
   define_v(tGrid,q,1);
   printf("lift = %e\n",-500*dot(tGrid,d,q,0,u_type));
}

#else

void p2c_update_drag_and_lift(n0,n1,n2,fa0,fa1,fa2,u,lift,drag)
NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; INT u; FLOAT *lift, *drag;
{  eprintf("Error: p2c_update_drag_and_lift not available.\n");  }

void lift_and_drag_for_p2c_p1c_2(tGrid,ZA,ZB,u,d,r,q,t_for_p,u_type,p_type,
                                 a_struct,b_struct)
GRID *tGrid; INT ZA,ZB,u,d,r,q,t_for_p,u_type,p_type,a_struct,b_struct;
{  eprintf("Error: lift_and_drag_for_p2c_p1c_2 not available.\n");  }

#endif

void update_drag_and_lift(pel,n0,n1,n2,fa0,fa1,fa2,u,space,lift,drag) 
ELEMENT *pel;
NODE *n0, *n1, *n2;   
FACE *fa0, *fa1, *fa2;  /* fa0 lies on the circle  */
INT u, space;
FLOAT *lift, *drag;
{
   switch(space){
   case P1_NC: p1nc_update_drag_and_lift(pel,n0,n1,n2,fa0,fa1,fa2,u,lift,drag);
        break;
   case P2C: p2c_update_drag_and_lift(n0,n1,n2,fa0,fa1,fa2,u,lift,drag);
        break;
   default:
        eprintf("Error: update_drag_and_lift not available.\n");
        break;
   }
}

INT is_on_circle(theNode,c,r)
NODE *theNode;
FLOAT *c, r;
{
   if (fabs((theNode->myvertex->x[0]-c[0])*(theNode->myvertex->x[0]-c[0])+
            (theNode->myvertex->x[1]-c[1])*(theNode->myvertex->x[1]-c[1])
            -r*r) < 1.e-10) return(1);
   else
      return(0);
}

#if DIM == 2

void compute_drag_and_lift_for_benchmark_channel(tGrid,u,space)
GRID *tGrid;
INT u, space;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT drag=0, lift=0, c[DIM], r=0.05;

   c[0] = c[1] = 0.2;
   for (pel = FIRSTELEMENT(tGrid);pel != NULL;pel = pel->succ){
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      if (is_on_circle(n0,c,r)+is_on_circle(n1,c,r)+is_on_circle(n2,c,r)>1){
         FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
         if (!is_on_circle(n0,c,r))
            update_drag_and_lift(pel,n0,n1,n2,fa0,fa1,fa2,u,space,&lift,&drag);
         else if (!is_on_circle(n1,c,r))
            update_drag_and_lift(pel,n1,n0,n2,fa1,fa0,fa2,u,space,&lift,&drag);
         else if (!is_on_circle(n2,c,r))
            update_drag_and_lift(pel,n2,n0,n1,fa2,fa0,fa1,u,space,&lift,&drag);
         else
            eprintf("Error in compute_drag_and_lift.\n");
      }
   }
   printf("drag = %e\n",-500.*drag);
   printf("lift = %e\n",500.*lift);
}

#else

void compute_drag_and_lift_for_benchmark_channel(tGrid,u,space)
GRID *tGrid; INT u, space;
{  eprintf("Error: compute_drag_and_lift_for_benchmark_channel not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (N_DATA & IxD_NODE_MATR) && (N_DATA & Dx1_NODE_FACE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (N_DATA & SCALAR_NODE_DATA) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void lift_and_drag_for_p2c_p1c(tGrid,ZA,ZB,x,y)
GRID *tGrid;
INT ZA,ZB,x,y;
{
   NODE *theNode, *pnode, *pn;
   FACE *theFace, *pface;
   LINK *pl;
   NFLINK *pnfl;
   DOUBLE sum[DIM];
	
   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ)
      SET_VALUE(pnode,0.,y)
   for (pnode = FDBN(tGrid); pnode; pnode = pnode->succ){
      SET4(NDD(pnode,y),COEFFBP(pnode,ZB),NDS(pnode,x))
      for (pl = pnode->tstart; pl!=pnode->start; pl = pl->next){
         pn=NBNODE(pl);
         SET4(NDD(pn,y),COEFFBLP(pl,ZB),NDS(pnode,x))
      }
   }
   sum[0] = sum[1] = 0.;
   for (theNode=FIRSTN(tGrid); theNode!=FIRSTNODE(tGrid); theNode=SUCC(theNode))
   if (theNode->myvertex->x[0] > 0.05 && theNode->myvertex->x[0] < 0.5 &&
       theNode->myvertex->x[1] > 0.1  && theNode->myvertex->x[1] < 0.3){
      SET5(sum,NDD(theNode,y))
      SET4(sum,NDD(theNode,x),COEFFN(theNode,ZA))
      for (pl=TSTART(theNode); pl; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,ZA))
      for (pnfl=TNFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
         SET4(sum,FDVP(NBFACE(pnfl),x),COEFFL(pnfl,ZA))
   }
   printf("drag = %e, lift = %e\n",-20/0.04*sum[0],-20/0.04*sum[1]);
}

#else

void lift_and_drag_for_p2c_p1c(tGrid,ZA,ZB,x,y)
GRID *tGrid; INT ZA,ZB,x,y;
{  eprintf("Error: lift_and_drag_for_p2c_p1c not available.\n");  }

#endif

#if (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExF_MATR) && (E_DATA & FxE_MATR) && (E_DATA & ExE_MATR) && (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA) && (E_DATA & ExFxDN_MATR) && (E_DATA & ExFxDF_MATR) && (E_DATA & ExDF_MATR) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)

void lift_and_drag_for_p2cb_p1disc(tGrid,ZA,ZB,x,y)
GRID *tGrid;
INT ZA,ZB,x,y;
{
   ELEMENT *pel;
   NODE *theNode, *pnode;
   LINK *pl;
   NFLINK *pnfl;
   DOUBLE sum[DIM];
   INT i, j;

   for (pnode = FIRSTN(tGrid); pnode != NULL; pnode = pnode->succ)
      SET_VALUE(pnode,0.,y)
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      for (i = 0; i < SIDES; i++)
         for (j = 0; j < NVERT; j++)
            SET4(NDD(pel->n[j],y),COEFF_BFDNP(pel,ZB,i,j),EDSN(pel,x,i))
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      ADD_VMULT(pel->n,y,COEFF_NEP(pel,ZA),EDVP(pel,x))
   sum[0] = sum[1] = 0.;
   for (theNode=FIRSTN(tGrid); theNode!=FIRSTNODE(tGrid); theNode=SUCC(theNode))
   if (theNode->myvertex->x[0] > 0.05 && theNode->myvertex->x[0] < 0.5 &&
       theNode->myvertex->x[1] > 0.1  && theNode->myvertex->x[1] < 0.3){
      SET5(sum,NDD(theNode,y))
      SET4(sum,NDD(theNode,x),COEFFN(theNode,ZA))
      for (pl=TSTART(theNode); pl; pl=NEXT(pl))
         SET4(sum,NDD(NBNODE(pl),x),COEFFL(pl,ZA))
      for (pnfl=TNFSTART(theNode); pnfl; pnfl=NEXT(pnfl))
         SET4(sum,FDVP(NBFACE(pnfl),x),COEFFL(pnfl,ZA))
   }
   printf("drag = %e, lift = %e\n",-20/0.04*sum[0],-20/0.04*sum[1]);
}

#else

void lift_and_drag_for_p2cb_p1disc(tGrid,ZA,ZB,x,y)
GRID *tGrid; INT ZA,ZB,x,y;
{  eprintf("Error: lift_and_drag_for_p2cb_p1disc not available.\n");  }

#endif

void lift_and_drag_for_benchmark_channel(tGrid,ZA,ZB,x,y,u_space,p_space)
GRID *tGrid;
INT ZA,ZB,x,y,u_space,p_space;
{
   if ((u_space == P2C && p_space == P1C) || 
       (u_space == IP2C && p_space == IP1C))
      lift_and_drag_for_p2c_p1c(tGrid,ZA,ZB,x,y);
   else if ((u_space == P2C_ELBUB && p_space == P1_DISC) || 
            (u_space == IP2C_ELBUB && p_space == IP1_DISC))
      lift_and_drag_for_p2cb_p1disc(tGrid,ZA,ZB,x,y);
   else
      eprintf("Error: lift_and_drag_for_benchmark_channel not available.\n");
}

#if N_DATA & SCALAR_NODE_DATA

FLOAT pressure_difference_for_p1c(vertexes,u)
VERTEX *vertexes;
INT u;
{
   return(NDS(vertexes[3].topnode,u)-NDS(vertexes[0].topnode,u));
}

#else

FLOAT pressure_difference_for_p1c(vertexes,u)
VERTEX *vertexes; INT u;
{  eprintf("Error: pressure_difference_for_p1c not available.\n");  }

#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES

FLOAT pressure_difference_for_p1disc(mg,vertexes,u)
MULTIGRID *mg;
VERTEX *vertexes;
INT u;
{
   ELEMENT *pel;
   NODE *n1=vertexes[3].topnode, *n2=vertexes[0].topnode;
   FLOAT s1=0., s2=0.;
   INT k1=0, k2=0;

   for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel != NULL; pel = pel->succ){
      if (pel->n[0] == n1){
         s1 += -EDSN(pel,u,0) + EDSN(pel,u,1) + EDSN(pel,u,2);
         k1++;
      }
      else if (pel->n[1] == n1){
         s1 += EDSN(pel,u,0) - EDSN(pel,u,1) + EDSN(pel,u,2);
         k1++;
      }
      else if (pel->n[2] == n1){
         s1 += EDSN(pel,u,0) + EDSN(pel,u,1) - EDSN(pel,u,2);
         k1++;
      }
      if (pel->n[0] == n2){
         s2 += -EDSN(pel,u,0) + EDSN(pel,u,1) + EDSN(pel,u,2);
         k2++;
      }
      else if (pel->n[1] == n2){
         s2 += EDSN(pel,u,0) - EDSN(pel,u,1) + EDSN(pel,u,2);
         k2++;
      }
      else if (pel->n[2] == n2){
         s2 += EDSN(pel,u,0) + EDSN(pel,u,1) - EDSN(pel,u,2);
         k2++;
      }
   }
   return(s1/k1 - s2/k2);
}

#else

FLOAT pressure_difference_for_p1disc(mg,vertexes,u)
MULTIGRID *mg; VERTEX *vertexes; INT u;
{  eprintf("Error: pressure_difference_for_p1disc not available.\n");  }

#endif

FLOAT pressure_difference_for_benchmark_channel(mg,vertexes,u,space)
MULTIGRID *mg;
VERTEX *vertexes;
INT u, space;
{
   switch(space){
   case P1C:      return(pressure_difference_for_p1c(vertexes,u));
        break;
   case IP1C:     return(pressure_difference_for_p1c(vertexes,u));
        break;
   case P1_DISC:  return(pressure_difference_for_p1disc(mg,vertexes,u));
        break;
   case IP1_DISC: return(pressure_difference_for_p1disc(mg,vertexes,u));
        break;
   default:
        eprintf("Error: pressure_difference_for_benchmark_channel not available.\n");
        break;
   }
}

