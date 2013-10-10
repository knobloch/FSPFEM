/******************************************************************************/
/*                                                                            */
/*                               initializations                              */
/*                                                                            */
/******************************************************************************/

#if E_DATA & SCALAR_ELEMENT_DATA

void set_ed_value(pelem,u,x)
ELEMENT *pelem;
FLOAT x;
INT u;
{
   ED(pelem,u) = x;
}

FLOAT get_ed_value(pelem,u)
ELEMENT *pelem;
INT u;
{
   return(ED(pelem,u));
}

#else

void set_ed_value(pelem,u,x)
ELEMENT *pelem; FLOAT x; INT u;
{  eprintf("Error: set_ed_value not available.\n");  }

FLOAT get_ed_value(pelem,u)
ELEMENT *pelem; INT u;
{  eprintf("Error: get_ed_value not available.\n");  }

#endif

#if E_DATA & VECTOR_ELEMENT_DATA

void set_edv_value(pelem,u,x)
ELEMENT *pelem;
FLOAT x[DIM];
INT u;
{
   SET1(EDVP(pelem,u),x)
}

void get_edv_value(pelem,u,x)
ELEMENT *pelem;
FLOAT x[DIM];
INT u;
{
   SET1(x,EDVP(pelem,u))
}

#else

void set_edv_value(pelem,u,x)
ELEMENT *pelem; FLOAT x[DIM]; INT u;
{  eprintf("Error: set_edv_value not available.\n");  }

void get_edv_value(pelem,u,x)
ELEMENT *pelem; FLOAT x[DIM]; INT u;
{  eprintf("Error: get_edv_value not available.\n");  }

#endif

#if N_DATA & SCALAR_NODE_DATA

void init_sn(mg,u,u0)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   NODE *pnode;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode !=NULL; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            NDS(pnode,u) = u0(pnode->myvertex->x);
}

#else

void init_sn(mg,u,u0)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: init_sn not available.\n");  }

#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES

void P1disc_init(mg,p,p0)
MULTIGRID *mg;
FLOAT (*p0)();
INT p;
{
   GRID *tGrid;
   ELEMENT *pel;
   FLOAT x01[DIM], x12[DIM], x02[DIM];

   for (tGrid = FIRSTGRID(mg); tGrid != NULL; tGrid = tGrid->finer)
      for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ){
         AVERAGE(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,x01);
         AVERAGE(pel->n[1]->myvertex->x,pel->n[2]->myvertex->x,x12);
         AVERAGE(pel->n[0]->myvertex->x,pel->n[2]->myvertex->x,x02);
         EDSN(pel,p,0) = p0(x12);
         EDSN(pel,p,1) = p0(x02);
         EDSN(pel,p,2) = p0(x01);
/*
         FLOAT v0, v1, v2, v01, v02, v12, c;

         Interpolation at vertices:

         v0 = p0(pel->n[0]->myvertex->x);
         v1 = p0(pel->n[1]->myvertex->x);
         v2 = p0(pel->n[2]->myvertex->x);
         EDSN(pel,p,0) = (v1+v2)*0.5;
         EDSN(pel,p,1) = (v0+v2)*0.5;
         EDSN(pel,p,2) = (v0+v1)*0.5;

         L2 projection:

         v12 = p0(x12);
         v02 = p0(x02);
         v01 = p0(x01);
         c = v0 + v1 + v2 + 2.*(v01 + v02 + v12);
         EDSN(pel,p,0) = 0.1*(c - 3.*v0 + 4.*v12);
         EDSN(pel,p,1) = 0.1*(c - 3.*v1 + 4.*v02);
         EDSN(pel,p,2) = 0.1*(c - 3.*v2 + 4.*v01);
*/
      }
}

#else

void P1disc_init(mg,p,p0)
MULTIGRID *mg; FLOAT (*p0)(); INT p;
{  eprintf("Error: P1disc_init not available.\n");  }

#endif

#if F_DATA & SCALAR_FACE_DATA

void set_edge_value(f12,x1,x2,u0,u)
FACE *f12;
FLOAT *x1, *x2, (*u0)();
INT u;
{
   FLOAT x12[DIM]; 

   if (IS_FF(f12)){
      AVERAGE(x1,x2,x12);
      FD(f12,u) = u0(x12);
   }
}

#else

void set_edge_value(f12,x1,x2,u0,u)
FACE *f12; FLOAT *x1, *x2, (*u0)(); INT u;
{  eprintf("Error: set_edge_value not available.\n");  }

#endif

#if (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void face_init(mg,u0,u)
MULTIGRID *mg;
FLOAT (*u0)();
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
         set_edge_value(f1,x2,x3,u0,u);
         set_edge_value(f2,x3,x1,u0,u);
         set_edge_value(f3,x1,x2,u0,u);
      }
}

#else

void face_init(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_init not available.\n");  }

#endif

#if DIM == 2

void vset_edge_value(a,x1,x2,u01,u02)
FLOAT *a, *x1, *x2, (*u01)(), (*u02)();
{
   FLOAT x12[DIM]; 

   AVERAGE(x1,x2,x12);
   a[0] = u01(x12);
   a[1] = u02(x12);
}

void vset_edge_mean_value(a,x1,x2,u01,u02)
FLOAT *a, *x1, *x2, (*u01)(), (*u02)();
{
   FLOAT x12[DIM], x11[DIM], x22[DIM], d[DIM], q=sqrt(15.)/10.;

   AVERAGE(x1,x2,x12);
   SET11(d,x2,x1)
   SET23(x11,x12,d,-q)
   SET23(x22,x12,d, q)
   a[0] = (5.*(u01(x11)+u01(x22)) + 8.*u01(x12))/18.;
   a[1] = (5.*(u02(x11)+u02(x22)) + 8.*u02(x12))/18.;
}

#else

void vset_edge_value(a,x1,x2,u01,u02)
FLOAT *a, *x1, *x2, (*u01)(), (*u02)();
{  eprintf("Error: vset_edge_value not available.\n");  }

void vset_edge_mean_value(a,x1,x2,u01,u02)
FLOAT *a, *x1, *x2, (*u01)(), (*u02)();
{  eprintf("Error: vset_edge_mean_value not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void vface_init(mg,u01,u02,u)
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
         vset_edge_value(FDVP(f1,u),x2,x3,u01,u02);
         vset_edge_value(FDVP(f2,u),x3,x1,u01,u02);
         vset_edge_value(FDVP(f3,u),x1,x2,u01,u02);
      }
}

void set_edge_values_mbub(f12,n1,n2,u0,u)
FACE *f12;
NODE *n1, *n2;
FLOAT (*u0)();
INT u;
{
   FLOAT v1, v2; 

   if (IS_FF(f12)){
      v1 = u0(n1->myvertex->x);
      v2 = u0(n2->myvertex->x);
      FDV(f12,u,0) = 0.5*(v1 + v2);
      FDV(f12,u,1) = 0.;
   }
}

void face_init_mbub(mg,u0,u)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         set_edge_values_mbub(f1,n2,n3,u0,u);
         set_edge_values_mbub(f2,n3,n1,u0,u);
         set_edge_values_mbub(f3,n1,n2,u0,u);
      }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         if (IS_FF(f1)) 
            FDV(f1,u,1) += 5.*(FDV(f2,u,0)-FDV(f3,u,0))*NINDI(n2,n3);
         if (IS_FF(f2))
            FDV(f2,u,1) += 5.*(FDV(f1,u,0)-FDV(f3,u,0))*NINDI(n1,n3);
         if (IS_FF(f3))
            FDV(f3,u,1) += 5.*(FDV(f1,u,0)-FDV(f2,u,0))*NINDI(n1,n2);
      }
}     

#else  /*  !( (F_DATA & VECTOR_FACE_DATA) && (DIM == 2) )  */

void vface_init(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: vface_init not available.\n");  }

void face_init_mbub(mg,u0,u)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: face_init_mbub not available.\n");  }

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == CUBE) && (DIM == 2)

void vq1rot_init(mg,u01,u02,u,vset_edge_value)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
void (*vset_edge_value)();
{
   GRID *theGrid;
   ELEMENT *pelem;
   FACE *f1, *f2, *f3, *f4;
   FLOAT *x1, *x2, *x3, *x4;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_4ELEMENT(f1,f2,f3,f4,pelem);
         VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pelem);
         vset_edge_value(FDVP(f1,u),x1,x2,u01,u02);
         vset_edge_value(FDVP(f2,u),x2,x3,u01,u02);
         vset_edge_value(FDVP(f3,u),x3,x4,u01,u02);
         vset_edge_value(FDVP(f4,u),x4,x1,u01,u02);
      }
}

#else

void vq1rot_init(mg,u01,u02,u,vset_edge_value)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u; void (*vset_edge_value)();
{  eprintf("Error: vq1rot_init not available.\n");  }

#endif

#if (F_DATA & DVECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void dvset_edge_values_mbub(f12,n1,n2,u01,u02,u)
FACE *f12;
NODE *n1, *n2;
FLOAT (*u01)(), (*u02)();
INT u;
{
   FLOAT v1, v2; 

   if (IS_FF(f12)){
      v1 = u01(n1->myvertex->x);
      v2 = u01(n2->myvertex->x);
      FDDV(f12,u,0,0) = 0.5*(v1 + v2);
      FDDV(f12,u,1,0) = 0.;
      v1 = u02(n1->myvertex->x);
      v2 = u02(n2->myvertex->x);
      FDDV(f12,u,0,1) = 0.5*(v1 + v2);
      FDDV(f12,u,1,1) = 0.;
   }
}

void dvface_init_mbub(mg,u01,u02,u)
MULTIGRID *mg;
FLOAT (*u01)(), (*u02)();
INT u;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n1, *n2, *n3;
   FACE *f1, *f2, *f3;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         dvset_edge_values_mbub(f1,n2,n3,u01,u02,u);
         dvset_edge_values_mbub(f2,n3,n1,u01,u02,u);
         dvset_edge_values_mbub(f3,n1,n2,u01,u02,u);
      }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
         FACES_OF_ELEMENT(f1,f2,f3,pelem);
         NODES_OF_ELEMENT(n1,n2,n3,pelem);
         if (IS_FF(f1)){
            FDDV(f1,u,1,0) += 5.*(FDDV(f2,u,0,0)-FDDV(f3,u,0,0))*NINDI(n2,n3);
            FDDV(f1,u,1,1) += 5.*(FDDV(f2,u,0,1)-FDDV(f3,u,0,1))*NINDI(n2,n3);
         }
         if (IS_FF(f2)){
            FDDV(f2,u,1,0) += 5.*(FDDV(f1,u,0,0)-FDDV(f3,u,0,0))*NINDI(n1,n3);
            FDDV(f2,u,1,1) += 5.*(FDDV(f1,u,0,1)-FDDV(f3,u,0,1))*NINDI(n1,n3);
         }
         if (IS_FF(f3)){
            FDDV(f3,u,1,0) += 5.*(FDDV(f1,u,0,0)-FDDV(f2,u,0,0))*NINDI(n1,n2);
            FDDV(f3,u,1,1) += 5.*(FDDV(f1,u,0,1)-FDDV(f2,u,0,1))*NINDI(n1,n2);
         }
      }
}     

#else  /*  if !((F_DATA & DVECTOR_FACE_DATA) && (DIM == 2))  */

void dvface_init_mbub(mg,u01,u02,u)
MULTIGRID *mg; FLOAT (*u01)(), (*u02)(); INT u;
{  eprintf("Error: dvface_init_mbub not available.\n");  }

#endif

#if DIM == 2

void p2cs_dofs(pel,v0,v1,v2,v01,v02,v12,f)
ELEMENT *pel;
FLOAT *v0, *v1, *v2, *v01, *v02, *v12, (*f)();
{
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM];

   VERTICES_OF_ELEMENT(x0,x1,x2,pel);
   AVERAGE(x0,x1,x01);
   AVERAGE(x0,x2,x02);
   AVERAGE(x1,x2,x12);
   *v0 = f(x0);
   *v1 = f(x1);
   *v2 = f(x2);
   *v01 = 4.*f(x01) - 2.*(*v0 + *v1);
   *v02 = 4.*f(x02) - 2.*(*v0 + *v2);
   *v12 = 4.*f(x12) - 2.*(*v1 + *v2);
}

#else

void p2cs_dofs(pel,v0,v1,v2,v01,v02,v12,f)
ELEMENT *pel; FLOAT *v0, *v1, *v2, *v01, *v02, *v12, (*f)();
{  eprintf("Error: p2cs_dofs not available.\n");  }

#endif

void p2cv_dofs(pel,v0,v1,v2,v01,v02,v12,f0,f1)
ELEMENT *pel;
FLOAT v0[DIM], v1[DIM], v2[DIM], v01[DIM], v02[DIM], v12[DIM], (*f0)(), (*f1)();
{
   p2cs_dofs(pel,&v0[0],&v1[0],&v2[0],&v01[0],&v02[0],&v12[0],f0);
   p2cs_dofs(pel,&v0[1],&v1[1],&v2[1],&v01[1],&v02[1],&v12[1],f1);
}

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void init_sn_sf(mg,u,u0,bubble)
MULTIGRID *mg;
FLOAT (*u0)();
INT u, bubble;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], xc[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode !=NULL; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            NDS(pnode,u) = u0(pnode->myvertex->x);

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            TOPNODES_OF_ELEMENT(n0,n1,n2,pelem);
            TOPFACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            AVERAGE(x1,x2,x12);
            FD(fa0,u) = 4.*u0(x12) - 2.*(NDS(n1,u)+NDS(n2,u));
            AVERAGE(x0,x2,x02);
            FD(fa1,u) = 4.*u0(x02) - 2.*(NDS(n0,u)+NDS(n2,u));
            AVERAGE(x0,x1,x01);
            FD(fa2,u) = 4.*u0(x01) - 2.*(NDS(n0,u)+NDS(n1,u));
            if (bubble){
               POINT3(x0,x1,x2,xc);
               set_ed_value(pelem,u,27.*u0(xc)
                                    -9.*(NDS(n0,u)+NDS(n1,u)+NDS(n2,u))
                                    -3.*(FD(fa0,u)+FD(fa1,u)+FD(fa2,u)));
            }
         }
}

#else

void init_sn_sf(mg,u,u0,bubble)
MULTIGRID *mg; FLOAT (*u0)(); INT u, bubble;
{  eprintf("Error: init_sn_sf not available.\n");  }

#endif

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & VECTOR_FACE_DATA) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void init_vn_vf_2(mg,u,u0,u1,u2,bubble)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)(), (*u2)();
INT u, bubble;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *pnode;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], xc[DIM], vxc[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode !=NULL; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
         }

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem)){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            FACES_OF_ELEMENT(fa0,fa1,fa2,pelem);
            VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
            AVERAGE(x1,x2,x12);
            FDV(fa0,u,0) = 4.*u0(x12) - 2.*(ND(n1,u,0)+ND(n2,u,0));
            FDV(fa0,u,1) = 4.*u1(x12) - 2.*(ND(n1,u,1)+ND(n2,u,1));
            AVERAGE(x0,x2,x02);
            FDV(fa1,u,0) = 4.*u0(x02) - 2.*(ND(n0,u,0)+ND(n2,u,0));
            FDV(fa1,u,1) = 4.*u1(x02) - 2.*(ND(n0,u,1)+ND(n2,u,1));
            AVERAGE(x0,x1,x01);
            FDV(fa2,u,0) = 4.*u0(x01) - 2.*(ND(n0,u,0)+ND(n1,u,0));
            FDV(fa2,u,1) = 4.*u1(x01) - 2.*(ND(n0,u,1)+ND(n1,u,1));
            if (bubble){
               POINT3(x0,x1,x2,xc);
               vxc[0] = 27.*u0(xc)-9.*( ND( n0,u,0)+ ND( n1,u,0)+ ND( n2,u,0))
                                  -3.*(FDV(fa0,u,0)+FDV(fa1,u,0)+FDV(fa2,u,0));
               vxc[1] = 27.*u1(xc)-9.*( ND( n0,u,1)+ ND( n1,u,1)+ ND( n2,u,1))
                                  -3.*(FDV(fa0,u,1)+FDV(fa1,u,1)+FDV(fa2,u,1));
               set_edv_value(pelem,u,vxc);
            }
         }
}

#else

void init_vn_vf_2(mg,u,u0,u1,u2,bubble)
MULTIGRID *mg; FLOAT (*u0)(), (*u1)(), (*u2)(); INT u, bubble;
{  eprintf("Error: init_vn_vf_2 not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

#if DIM == 2

void init_vn(mg,u,u0,u1,u2)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)(), (*u2)();
INT u;
{
   GRID *theGrid;
   NODE *pnode;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode !=NULL; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
         }
}

#else  /*  if (DIM == 3)  */

void init_vn(mg,u,u0,u1,u2)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)(), (*u2)();
INT u;
{
   GRID *theGrid;
   NODE *n0, *n1, *n2, *pnode;

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode !=NULL; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode)){
            ND(pnode,u,0) = u0(pnode->myvertex->x);
            ND(pnode,u,1) = u1(pnode->myvertex->x);
            ND(pnode,u,2) = u2(pnode->myvertex->x);
         }
}

#endif  /*  DIM == 2  */

#else

void init_vn(mg,u,u0,u1,u2)
MULTIGRID *mg; FLOAT (*u0)(), (*u1)(), (*u2)(); INT u;
{  eprintf("Error: init_vn not available.\n");  }

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void coord_of_barycentre();

#if DIM == 3

void P0_init(mg,ue,p0)
MULTIGRID *mg;
FLOAT (*p0)();
INT ue;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT xc[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            POINT4(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                   pel->n[2]->myvertex->x,pel->n[3]->myvertex->x,xc);
            ED(pel,ue) = p0(xc);
         }
}

#else  /* if DIM == 2 */

void P0_init(mg,ue,p0)
MULTIGRID *mg;
FLOAT (*p0)();
INT ue;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT xc[DIM];

   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            coord_of_barycentre(pel,xc);
            ED(pel,ue) = p0(xc);
         }
}

#endif  /*  DIM == 2  */

#else

void P0_init(mg,ue,p0)
MULTIGRID *mg; FLOAT (*p0)(); INT ue;
{  eprintf("Error: P0_init not available.\n");  }

#endif

#if DATA_S & SPECIAL_NODES_AND_FACES

void p2l_init(mg,u,u0)
MULTIGRID *mg;
FLOAT (*u0)();
INT u;
{
   GRID *theGrid;
   SFACE *psf;
   DOUBLE x12[DIM], u1, u2;

   for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
      for (psf = FIRSTSF(theGrid); psf; psf = psf->succ)
         if (IS_LF(psf->f)){
            if (IS_LN(psf->n1))
               u1 = SNDS(psf->n1->s_node,u) = u0(psf->n1->myvertex->x);
            else
               u1 = 0.;
            if (IS_LN(psf->n2))
               u2 = SNDS(psf->n2->s_node,u) = u0(psf->n2->myvertex->x);
            else
               u2 = 0.;
            AVERAGE(psf->n1->myvertex->x,psf->n2->myvertex->x,x12);
            SFDS(psf,u) = 4.*u0(x12) - 2.*(u1 + u2);
         }
} 

#else

void p2l_init(mg,u,u0)
MULTIGRID *mg; FLOAT (*u0)(); INT u;
{  eprintf("Error: p2l_init not available.\n");  }

#endif

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void general_init(mg,u,u0,rm_type,finite_el)
MULTIGRID *mg;
FLOAT (*u0)();
INT u, rm_type;
FINITE_ELEMENT finite_el;
{
   GRID *theGrid;
   ELEMENT *pel;
   REF_MAPPING ref_map;
   DOUBLE_FUNC *nd, *fd, *ed;
   INT i, k, m, kk, nn, mm, dir[4];

   nd = finite_el.ndofs;
   fd = finite_el.fdofs;
   ed = finite_el.edofs;
   kk = finite_el.k;
   nn = finite_el.n;
   mm = finite_el.m;
   for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            reference_mapping(pel,rm_type,&ref_map);
            if (kk)
               for (i=0; i < NVERT; i++){
                  m = i*kk;
                  for (k = m; k < m+kk; k++) 
                     NDMV(pel->n[i],u,k-m) = (nd[k])(u0,&ref_map);
               }
            if (nn){
               set_directions_of_edges(pel,dir);
               for (i=0; i < NVERT; i++){
                  m = i*nn;
                  for (k = m; k < m+nn; k++) 
                     FDMV(pel->f[i],u,f_i(k-m,nn,dir[i])) = (fd[k])(u0,&ref_map);
               }
            }
            for (i = 0; i < mm; i++)
               EDMV(pel,u,i) = (ed[i])(u0,&ref_map);
         }
}

#else

void general_init(mg,u,u0,rm_type,finite_el)
MULTIGRID *mg; FLOAT (*u0)(); INT u, rm_type; FINITE_ELEMENT finite_el;
{  eprintf("Error: general_init not available.\n");  }

#endif

void initialize(mg,u,u0,u1,u2,space,structure)
MULTIGRID *mg;
FLOAT (*u0)(), (*u1)(), (*u2)();
INT u, space, structure;
{
   switch(space){
   case P0:   if (structure == SCALAR || structure == P_SCALAR)
                 P0_init(mg,u,u0);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P1_NC:if (structure == SCALAR)
                 face_init(mg,u0,u);
              else if (structure == VECTOR)
                 vface_init(mg,u0,u1,u);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P1_MOD:if (structure == SCALAR)
                 face_init_mbub(mg,u0,u);
               else if (structure == VECTOR)
                 dvface_init_mbub(mg,u0,u1,u);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case Q1ROT:if (structure == SCALAR)
                 eprintf("Error: initialize not available.\n");
              else if (structure == VECTOR || DOF == MIDPOINT)
                 vq1rot_init(mg,u0,u1,u,vset_edge_value);
              else if (structure == VECTOR || DOF == MEAN_VALUE)
                 vq1rot_init(mg,u0,u1,u,vset_edge_mean_value);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P1C:  if (structure == SCALAR || structure == P_SCALAR)
                 init_sn(mg,u,u0);
              else if (structure == VECTOR)
                 init_vn(mg,u,u0,u1,u2);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P1_DISC:  if (structure == SCALAR || structure == P_SCALAR)
                 P1disc_init(mg,u,u0);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case MINI_L_DIV_FR: if (structure == SCALAR)
                 init_sn(mg,u,u0);
              else if (structure == VECTOR)
                 init_vn(mg,u,u0,u1,u2);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case Q1C:  if (structure == SCALAR || structure == P_SCALAR)
                 init_sn(mg,u,u0);
              else if (structure == VECTOR)
                 init_vn(mg,u,u0,u1,u2);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P2C:  if (structure == SCALAR  || structure == P_SCALAR)
                 init_sn_sf(mg,u,u0,0);
              else if (structure == VECTOR)
                 init_vn_vf_2(mg,u,u0,u1,u2,0);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P2C_ELBUB:  if (structure == SCALAR  || structure == P_SCALAR)
                 init_sn_sf(mg,u,u0,1);
              else if (structure == VECTOR)
                 init_vn_vf_2(mg,u,u0,u1,u2,1);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case P2L:  p2l_init(mg,u,u0);
        break;
   case GP1C: 
              if (structure == SCALAR || structure == P_SCALAR)
                 general_init(mg,u,u0,REF_MAP,p1c_element);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case GQ1C: 
              if (structure == SCALAR || structure == P_SCALAR)
                 general_init(mg,u,u0,REF_MAP,q1c_element);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   case GQ2C: 
              if (structure == SCALAR || structure == P_SCALAR)
                 general_init(mg,u,u0,REF_MAP,q2c_element);
              else
                 eprintf("Error: initialize not available.\n");
        break;
   default:
        eprintf("Error: initialize not available.\n");
        break;
   }
}

