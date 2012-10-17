/******************************************************************************/
/*                                                                            */
/*                             vector operations                              */
/*                                                                            */
/******************************************************************************/

#if N_DATA & SCALAR_NODE_DATA

void sb_copy(tGrid,x,z,t)  /* z := x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid;
INT x, z, t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   if (t & USE_IS_DIR){
      for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); 
                                                          theNode=SUCC(theNode))
         if(IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
         for (theNode=FDBN(theGrid); theNode != FIRSTNODE(theGrid); 
                                                          theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x);
   }
   else{
      for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); 
                                                          theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
         for (theNode=FDBN(theGrid); theNode != FIRSTNODE(theGrid); 
                                                          theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x);
   }
}

#else  /*  if !(N_DATA & SCALAR_NODE_DATA)  */

void sb_copy(tGrid,x,z,t)  /* z := x;  t=1 ... temp., t=2 ... conc. */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: sb_copy not available.\n");  }

#endif  /*  !(N_DATA & SCALAR_NODE_DATA)  */

#if F_DATA & SCALAR_FACE_DATA

void b_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FDBF(tGrid); theFace != FIRSTFACE(tGrid); theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FDBF(theGrid); theFace != FIRSTFACE(theGrid); 
                                                          theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x);
}

#else  /*  if !(F_DATA & SCALAR_FACE_DATA)  */

void b_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: b_copy_f not available.\n");  }

#endif  /*  !(F_DATA & SCALAR_FACE_DATA)  */

#if N_DATA & VECTOR_NODE_DATA

void v_b_copy_n(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      COPY(theNode,x,z)
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theNode=FDBN(theGrid); theNode != FIRSTNODE(theGrid); 
                                                        theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            COPY(theNode,x,z)
}

#else  /* if !(N_DATA & VECTOR_NODE_DATA)  */

void v_b_copy_n(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: v_b_copy_n not available.\n");  }

#endif  /* if !(N_DATA & VECTOR_NODE_DATA)  */

#if F_DATA & VECTOR_FACE_DATA

void vb_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FDBF(tGrid); theFace != FIRSTFACE(tGrid); theFace=SUCC(theFace))
      COPY(theFace,x,z)
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FDBF(theGrid); theFace != FIRSTFACE(theGrid); 
                                                          theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            COPY(theFace,x,z)
}

#else  /* if !(F_DATA & VECTOR_FACE_DATA)  */

void vb_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: vb_copy_f not available.\n");  }

#endif  /*  !(F_DATA & VECTOR_FACE_DATA)  */

#if F_DATA & DVECTOR_FACE_DATA

void dvb_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FDBF(tGrid); theFace != FIRSTFACE(tGrid); theFace=SUCC(theFace))
      D_COPY(theFace,x,z)
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser)
      for (theFace=FDBF(theGrid); theFace != FIRSTFACE(theGrid); 
                                                          theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_COPY(theFace,x,z)
}

#else  /*  if !(F_DATA & DVECTOR_FACE_DATA)  */

void dvb_copy_f(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: dvb_copy_f not available.\n");  }

#endif  /*  !(F_DATA & DVECTOR_FACE_DATA)  */

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void vs_b_copy_nf(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      COPY(theNode,x,z)
   for (theFace=FDBF(tGrid); theFace != FIRSTFACE(tGrid); theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid != NULL; theGrid = theGrid->coarser){
      for (theNode=FDBN(theGrid); theNode != FIRSTNODE(theGrid); 
                                                        theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            COPY(theNode,x,z)
      for (theFace=FDBF(theGrid); theFace != FIRSTFACE(theGrid); 
                                                        theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x);
   }	
}

#else  /* if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void vs_b_copy_nf(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: vs_b_copy_nf not available.\n");  }

#endif  /* if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DATA_STR & LG_DATA)

void vs_b_copy_nflg(tGrid,x,z)  /* z := x */
GRID *tGrid;
INT x,z;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FDBN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      COPY(theNode,x,z)
   for (theFace=FDBF(tGrid); theFace != FIRSTFACE(tGrid); theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
}

#else  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                (DATA_STR & LG_DATA))  */

void vs_b_copy_nflg(tGrid,x,z)  /* z := x */
GRID *tGrid; INT x,z;
{  eprintf("Error: vs_b_copy_nflg not available.\n");  }

#endif  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                 (DATA_STR & LG_DATA))  */

void b_copy(tGrid,x,z,t,type)  /* z := x */
GRID *tGrid;
INT x, z, t, type;
{ 
   if(operation_available(t,&type,"b_copy"))
   switch(type){
   case Q_SN: sb_copy(tGrid,x,z,t);
        break;
   case Q_VN: v_b_copy_n(tGrid,x,z);
        break;
   case Q_SF: b_copy_f(tGrid,x,z);
        break;
   case Q_VF: vb_copy_f(tGrid,x,z);
        break;
   case Q_DVF: dvb_copy_f(tGrid,x,z);
        break;
   case Q_SNSF: sb_copy(tGrid,x,z,t);
                b_copy_f(tGrid,x,z);
        break;
   case Q_VNSF: vs_b_copy_nf(tGrid,x,z);
        break;
   case Q_VNVF: v_b_copy_n(tGrid,x,z);
                vb_copy_f(tGrid,x,z);
        break;
   case Q_VNSFLG: vs_b_copy_nflg(tGrid,x,z);
        break;
   default:
        eprintf("Error: b_copy not available.\n");
        break;
   }
}

