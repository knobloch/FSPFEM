/******************************************************************************/
/*                                                                            */
/*                             vector operations                              */
/*                                                                            */
/******************************************************************************/

NODE *ltop_node(n,tGrid)
NODE *n;
GRID *tGrid;
{
   while (n->son && INDEX(n->son) <= tGrid->maxNodeIndex) n = n->son;
   return(n);
}

FACE *ltop_face(f,tGrid)
FACE *f;
GRID *tGrid;
{
   if (f->sons[0] == f)
      return(f);
   while (f->sons[0] && INDEX(f->sons[0]) <= tGrid->maxFaceIndex
                     && INDEX(f->sons[0]) > 0) f = f->sons[0];
   return(f);
}

#if N_DATA & SCALAR_NODE_DATA

void sset_value(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = r;
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = r;
      }
   }
}

DOUBLE smax_abs_value(tGrid,z,t)  /* max := max(z_i) */
GRID *tGrid;
INT z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   DOUBLE max=0.;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            max = MAX(max,fabs(NDS(theNode,z)));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               max = MAX(max,fabs(NDS(theNode,z)));
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         max = MAX(max,fabs(NDS(theNode,z)));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               max = MAX(max,fabs(NDS(theNode,z)));
      }
   }
   return(max);
}

void sadd_value(tGrid,r,z,t)  /* z_i := z_i + r */
GRID *tGrid;
FLOAT r;
INT z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) += r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) += r;
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) += r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) += r;
      }
   }
}

DOUBLE ssum(tGrid,z,t)
GRID *tGrid;
INT z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   DOUBLE sum=0.;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            sum += NDS(theNode,z);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               sum += NDS(theNode,z);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         sum += NDS(theNode,z);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               sum += NDS(theNode,z);
      }
   }
   return(sum);
}

DOUBLE ssum_n(tGrid,z,n,t)
GRID *tGrid;
INT z, *n, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   DOUBLE sum=0.;
   	
   *n = 0;
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t)){
            sum += NDS(theNode,z);
            (*n)++;
         }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t)){
               sum += NDS(theNode,z);
               (*n)++;
            }
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode)){
         sum += NDS(theNode,z);
         (*n)++;
      }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid)){
               sum += NDS(theNode,z);
               (*n)++;
            }
      }
   }
   return(sum);
}

void sset_first(tGrid,r,z)  /* z_1 := r */
GRID *tGrid;
FLOAT r;
INT z;
{
   NDS(FIRSTN(tGrid),z) = r;
}

void ssubtr_first(tGrid,z)  /* z_i := z_i - z_1 */
GRID *tGrid;
INT z;
{
   GRID *theGrid;
   NODE *theNode;
   FLOAT f;
   	
   f = NDS(FIRSTN(tGrid),z);
   for (theNode=FIRSTN(tGrid); theNode != NULL; theNode=SUCC(theNode))
      NDS(theNode,z) -= f;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRSTN(theGrid); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NDS(theNode,z) -= f;
}

DOUBLE smean_value(tGrid,z)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid;
INT z;
{
   GRID *theGrid;
   ELEMENT *pelem;
   FLOAT sum=0., vol=0., evol;
  
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      evol = VOLUME(pelem);
      sum += AVER_NS(pelem,z)*evol;
      vol += evol;
   }

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            evol = VOLUME(pelem);
            sum += AVER_NS(pelem,z)*evol;
            vol += evol;
         }
   return(sum/vol);
}

void sset_zero_mean(tGrid,z)
GRID *tGrid;
INT z;
{
   GRID *theGrid;
   NODE *theNode;
   DOUBLE mean;
  
   mean = smean_value(tGrid,z);

   for (theNode=FIRSTN(tGrid); theNode != NULL; theNode=SUCC(theNode))
      NDS(theNode,z) -= mean;
  
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRSTN(theGrid); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NDS(theNode,z) -= mean;
   }
}

void scopy(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x);
      }
   }
}

void mso_scopy(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
		 if (t & STOP_IS_FIRST_INNER) {
			if(IS_DIR(theNode,t))
				NDS(theNode,z) = NDS(theNode,x);
		 } else {
			 if(!IS_DIR(theNode,t))
				NDS(theNode,z) = NDS(theNode,x);
		 }
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
			 if (t & STOP_IS_FIRST_INNER) {
				if (IS_LTOP_NODE(theNode,tGrid) && IS_DIR(theNode,t))
					NDS(theNode,z) = NDS(theNode,x);
			 } else {
				 if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
					NDS(theNode,z) = NDS(theNode,x);
			 }
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x);
      }
   }
}

void sinv(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = -NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = -NDS(theNode,x);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = -NDS(theNode,x);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = -NDS(theNode,x);
      }
   }
}

void sadd(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x) + NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x) + NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x) + NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x) + NDS(theNode,y);
      }
   }
}

void ssubtr(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x) - NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x) - NDS(theNode,y);
      }
   }
}

void sadd_and_subtr(tGrid,w,x,y,z,t)  /*  z := w + x - y  */
GRID *tGrid;
INT w, x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,w) + NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,w) + NDS(theNode,x) - NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,w) + NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,w) + NDS(theNode,x) - NDS(theNode,y);
      }
   }
}

void smult_and_add(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = r*NDS(theNode,x) + NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = r*NDS(theNode,x) + NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = r*NDS(theNode,x) + NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = r*NDS(theNode,x) + NDS(theNode,y);
      }
   }
}

void smult_and_subtr(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = r*NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = r*NDS(theNode,x) - NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = r*NDS(theNode,x) - NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = r*NDS(theNode,x) - NDS(theNode,y);
      }
   }
}

void sdamp(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x) + r*(NDS(theNode,y) - NDS(theNode,x));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x) + r*(NDS(theNode,y) - NDS(theNode,x));
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x) + r*(NDS(theNode,y) - NDS(theNode,x));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x) + r*(NDS(theNode,y) - NDS(theNode,x));
      }
   }
}

void smult(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x)*r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x)*r;
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x)*r;
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x)*r;
      }
   }
}

DOUBLE sdot(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x, y, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   DOUBLE sum=0.0;
	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            sum +=  NDS(theNode,x)*NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               sum +=  NDS(theNode,x)*NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         sum +=  NDS(theNode,x)*NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               sum +=  NDS(theNode,x)*NDS(theNode,y);
      }
   }
   return(sum);
}

DOUBLE sdot_without_first(tGrid,x,y)
GRID *tGrid;
INT x, y;
{
   GRID *theGrid;
   NODE *theNode;
   DOUBLE sum=0.0;

   for (theNode=SUCC(FIRSTN(tGrid)); theNode != NULL; theNode=SUCC(theNode))
      sum +=  NDS(theNode,x)*NDS(theNode,y);

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRSTN(theGrid); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            sum +=  NDS(theNode,x)*NDS(theNode,y);
   return(sum);
}

void sdivide(tGrid,x,y,z,t)  /* z_i := x_i/y_i */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x)/NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x)/NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x)/NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x)/NDS(theNode,y);
      }
   }
}

void smultiply(tGrid,x,y,z,t)  /* z_i := x_i*y_i */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = NDS(theNode,x)*NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = NDS(theNode,x)*NDS(theNode,y);
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = NDS(theNode,x)*NDS(theNode,y);
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = NDS(theNode,x)*NDS(theNode,y);
      }
   }
}

void smake_sqrt(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         if(!IS_DIR(theNode,t))
            NDS(theNode,z) = sqrt(NDS(theNode,x));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid) && !IS_DIR(theNode,t))
               NDS(theNode,z) = sqrt(NDS(theNode,x));
      }
   }
   else{
      for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
         NDS(theNode,z) = sqrt(NDS(theNode,x));
   
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid))
               NDS(theNode,z) = sqrt(NDS(theNode,x));
      }
   }
}

void smult_and_add1_for_GMRES(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{
   GRID *theGrid;
   NODE *pnode, *stop;
   INT i;
   
   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (pnode = FIRST_NODE(tGrid,t); pnode != stop; pnode = pnode->succ)
         if (!IS_DIR(pnode,t))
            for (i = 0; i <= j; i++)
               NDS(pnode,l) -= h[i][j]*NDS(pnode,k+i);
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (pnode = FIRST_NODE(theGrid,t); pnode != stop; pnode = pnode->succ)
            if (IS_LTOP_NODE(pnode,tGrid) && !IS_DIR(pnode,t))
               for (i = 0; i <= j; i++)
                  NDS(pnode,l) -= h[i][j]*NDS(pnode,k+i);
      }
   }
   else{
      for (pnode = FIRST_NODE(tGrid,t); pnode != stop; pnode = pnode->succ)
         for (i = 0; i <= j; i++)
            NDS(pnode,l) -= h[i][j]*NDS(pnode,k+i);
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (pnode = FIRST_NODE(theGrid,t); pnode != stop; pnode = pnode->succ)
            if (IS_LTOP_NODE(pnode,tGrid))
               for (i = 0; i <= j; i++)
                  NDS(pnode,l) -= h[i][j]*NDS(pnode,k+i);
      }
   }
}

void smult_and_add2_for_GMRES(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{
   GRID *theGrid;
   NODE *pnode, *stop;
   INT i;

   stop = STOP_NODE(tGrid,t);
   if (t & USE_IS_DIR){
      for (pnode = FIRST_NODE(tGrid,t); pnode != stop; pnode = pnode->succ)
         if (!IS_DIR(pnode,t))
            for (i = 0; i <= j; i++)
               NDS(pnode,x) += g[i]*NDS(pnode,k+i);
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (pnode = FIRST_NODE(theGrid,t); pnode != stop; pnode = pnode->succ)
            if (IS_LTOP_NODE(pnode,tGrid) && !IS_DIR(pnode,t))
               for (i = 0; i <= j; i++)
                  NDS(pnode,x) += g[i]*NDS(pnode,k+i);
      }
   }
   else{
      for (pnode = FIRST_NODE(tGrid,t); pnode != stop; pnode = pnode->succ)
         for (i = 0; i <= j; i++)
            NDS(pnode,x) += g[i]*NDS(pnode,k+i);
      for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
         stop = STOP_NODE(theGrid,t);
         for (pnode = FIRST_NODE(theGrid,t); pnode != stop; pnode = pnode->succ)
            if (IS_LTOP_NODE(pnode,tGrid))
               for (i = 0; i <= j; i++)
                  NDS(pnode,x) += g[i]*NDS(pnode,k+i);
      }
   }
}

void save_sn_vector(tGrid,z,m,fp,t)
GRID *tGrid;
INT z, m, t;
FILE *fp;
{
   NODE *theNode, *stop=FIRSTNODE(tGrid);

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      fprintf(fp,"%i %e\n",m+theNode->index,NDS(theNode,z));
   if (!(t & NSTART_FROM_INNER))
      for (theNode=FIRSTN(tGrid); theNode != stop; theNode=SUCC(theNode))
         fprintf(fp,"%i %e\n",m+theNode->index,NDS(theNode,z));
}

void sn_set_zero_on_2nd_periodic_boundary(tGrid,z)
GRID *tGrid;
INT z;
{
   P_NODE *pp_node;

   for (pp_node = FIRSTPN(tGrid); pp_node; pp_node = pp_node->succ)
      NDS(pp_node->n2,z) = 0.;
}

void sn_copy_1st_to_2nd_periodic_boundary(tGrid,z)
GRID *tGrid;
INT z;
{
   P_NODE *pp_node;

   for (pp_node = FIRSTPN(tGrid); pp_node; pp_node = pp_node->succ)
      NDS(pp_node->n2,z) = NDS(pp_node->n1,z);
}

void sn_add_2nd_to_1st_periodic_boundary(tGrid,z)
GRID *tGrid;
INT z;
{
   P_NODE *pp_node;

   for (pp_node = FIRSTPN(tGrid); pp_node; pp_node = pp_node->succ)
      NDS(pp_node->n1,z) += NDS(pp_node->n2,z);
}

#else  /*  if !(N_DATA & SCALAR_NODE_DATA)  */

void sset_value(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z, t;
{  eprintf("Error: sset_value not available.\n");  }

DOUBLE smax_abs_value(tGrid,z,t)  /* max := max(z_i) */
GRID *tGrid; INT z, t;
{  eprintf("Error: smax_abs_value not available.\n"); return(0.);  }

void sadd_value(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z, t;
{  eprintf("Error: sadd_value not available.\n");  }

DOUBLE ssum(tGrid,z,t)
GRID *tGrid; INT z, t;
{  eprintf("Error: ssum not available.\n"); return(0.);  }

DOUBLE ssum_n(tGrid,z,n,t)
GRID *tGrid; INT z, *n, t;
{  eprintf("Error: ssum_n not available.\n"); return(0.);  }

void sset_first(tGrid,r,z)  /* z_1 := r */
GRID *tGrid; FLOAT r; INT z;
{  eprintf("Error: sset_first not available.\n");  }

void ssubtr_first(tGrid,z)  /* z_i := z_i - z_1 */
GRID *tGrid; INT z;
{  eprintf("Error: ssubtr_first not available.\n");  }

DOUBLE smean_value(tGrid,z)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid; INT z;
{  eprintf("Error: smean_value not available.\n"); return(0.);  }

void sset_zero_mean(tGrid,z)
GRID *tGrid; INT z;
{  eprintf("Error: sset_zero_mean not available.\n");  }

void scopy(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: scopy not available.\n");  }

void mso_scopy(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: mso_scopy not available.\n");  }

void sinv(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: sinv not available.\n");  }

void sadd(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: sadd not available.\n");  }

void ssubtr(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: ssubtr not available.\n");  }

void sadd_and_subtr(tGrid,w,x,y,z,t)  /*  z := w + x - y  */
GRID *tGrid; INT w, x, y, z, t;
{  eprintf("Error: sadd_and_subtr not available.\n");  }

void smult_and_add(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x, y, z, t;
{  eprintf("Error: smult_and_add not available.\n");  }

void smult_and_subtr(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x, y, z, t;
{  eprintf("Error: smult_and_subtr not available.\n");  }

void sdamp(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x, y, z, t;
{  eprintf("Error: sdamp not available.\n");  }

void smult(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x, z, t;
{  eprintf("Error: smult not available.\n");  }

DOUBLE sdot(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x, y, t;
{  eprintf("Error: sdot not available.\n"); return(0.);  }

DOUBLE sdot_without_first(tGrid,x,y)
GRID *tGrid; INT x, y;
{  eprintf("Error: sdot_without_first not available.\n"); return(0.);  }

void sdivide(tGrid,x,y,z,t)  /* z_i := x_i/y_i */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: sdivide not available.\n");  }

void smultiply(tGrid,x,y,z,t)  /* z_i := x_i*y_i */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: smultiply not available.\n");  }

void smake_sqrt(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: smake_sqrt not available.\n");  }

void smult_and_add1_for_GMRES(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: smult_and_add1_for_GMRES not available.\n");  }

void smult_and_add2_for_GMRES(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: smult_and_add2_for_GMRES not available.\n");  }

void save_sn_vector(tGrid,z,m,fp,t)
GRID *tGrid; INT z, m, t; FILE *fp;
{  eprintf("Error: save_sn_vector not available.\n");  }

void sn_set_zero_on_2nd_periodic_boundary(tGrid,z)
GRID *tGrid; INT z;
{  eprintf("Error: sn_set_zero_on_2nd_periodic_boundary not available.\n");  }

void sn_copy_1st_to_2nd_periodic_boundary(tGrid,z)
GRID *tGrid; INT z;
{  eprintf("Error: sn_copy_1st_to_2nd_periodic_boundary not available.\n");  }

void sn_add_2nd_to_1st_periodic_boundary(tGrid,z)
GRID *tGrid; INT z;
{  eprintf("Error: sn_add_2nd_to_1st_periodic_boundary not available.\n");  }

#endif  /*  !(N_DATA & SCALAR_NODE_DATA)  */

#if F_DATA & SCALAR_FACE_DATA

void set_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      FD(theFace,z) = r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r;
   }
}

void set_first_f(tGrid,r,z)  /* z_1 := r */
GRID *tGrid;
FLOAT r;
INT z;
{ 
   FD(FIRSTF(tGrid),z) = r;
}

void subtr_first_f(tGrid,z)  /* z_i := z_i - z_1 */
GRID *tGrid;
INT z;
{ 
   GRID *theGrid;
   FACE *theFace;
   FLOAT f;
   	
   f = FD(FIRSTF(tGrid),z);
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) -= f;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) -= f;
}

void copy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x);
   }
}

void inv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = -FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = -FD(theFace,x);
}

void add_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) + FD(theFace,y);
}

void subtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) - FD(theFace,y);
   }
}

void add_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid;
INT w,x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,w) + FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,w) + FD(theFace,x) - FD(theFace,y);
}

void mult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = r*FD(theFace,x) + FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r*FD(theFace,x) + FD(theFace,y);
}

void mult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = r*FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r*FD(theFace,x) - FD(theFace,y);
}

void damp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + r*(FD(theFace,y) - FD(theFace,x));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) + r*(FD(theFace,y) - FD(theFace,x));
}

void mult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)*r;
}

DOUBLE dot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      sum += FD(theFace,x)*FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            sum += FD(theFace,x)*FD(theFace,y);
   return(sum);
}

DOUBLE dot_f_without_first(tGrid,x,y)  /* := x.y */
GRID *tGrid;
INT x,y;
{ 
   GRID *theGrid;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theFace = SUCC(FIRSTF(tGrid)); theFace != NULL; theFace=SUCC(theFace))
      sum += FD(theFace,x)*FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            sum += FD(theFace,x)*FD(theFace,y);
   return(sum);
}

void divide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)/FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)/FD(theFace,y);
}

void multiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)*FD(theFace,y);
}

void make_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = sqrt(FD(theFace,x));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = sqrt(FD(theFace,x));
}

void mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;
   
   for (pface = FIRST_FACE(tGrid,t); pface!=NULL; pface=pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,l) -=  h[i][j]*FD(pface,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               FD(pface,l) -=  h[i][j]*FD(pface,k+i);
}

void mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;

   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,x) += g[i]*FD(pface,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)   
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               FD(pface,x) += g[i]*FD(pface,k+i);
}

void save_sf_vector(tGrid,z,m,fp,t)
GRID *tGrid;
INT z, m, t;
FILE *fp;
{
   FACE *theFace, *stop=FIRSTFACE(tGrid);

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FD(theFace,z));
   if (!(t & FSTART_FROM_INNER))
      for (theFace=FIRSTF(tGrid); theFace != stop; theFace=SUCC(theFace))
         fprintf(fp,"%i %e\n",m+theFace->index,FD(theFace,z));
}

#else  /*  if !(F_DATA & SCALAR_FACE_DATA)  */

void set_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid; FLOAT r; INT z, t;
{  eprintf("Error: set_value_f not available.\n");  }

void set_first_f(tGrid,r,z)  /* z_1 := r */
GRID *tGrid; FLOAT r; INT z;
{  eprintf("Error: set_first_f not available.\n");  }
 
void subtr_first_f(tGrid,z)  /* z_i := z_i - z_1 */
GRID *tGrid; INT z;
{  eprintf("Error: subtr_first_f not available.\n");  }

void copy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: copy_f not available.\n");  }

void inv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: inv_f not available.\n");  }

void add_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: add_f not available.\n");  }

void subtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: subtr_f not available.\n");  }

void add_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: add_and_subtr_f not available.\n");  }

void mult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mult_and_add_f not available.\n");  }

void mult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mult_and_subtr_f not available.\n");  }

void damp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: damp_f not available.\n");  }

void mult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: mult_f not available.\n");  }

DOUBLE dot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: dot_f not available.\n"); return(0.);  }

DOUBLE dot_f_without_first(tGrid,x,y)  /* := x.y */
GRID *tGrid; INT x,y;
{  eprintf("Error: dot_f_without_first not available.\n"); return(0.);  }

void divide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: divide_f not available.\n");  }

void multiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: multiply_f not available.\n");  }

void make_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: make_sqrt_f not available.\n");  }

void mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: mult_and_add1_for_GMRES_f not available.\n");  }

void mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: mult_and_add2_for_GMRES_f not available.\n");  }

void save_sf_vector(tGrid,z,m,fp,t)
GRID *tGrid; INT z, m, t; FILE *fp;
{  eprintf("Error: save_sf_vector not available.\n");  }

#endif  /*  !(F_DATA & SCALAR_FACE_DATA)  */

#if E_DATA & SCALAR_ELEMENT_DATA

void set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid;
FLOAT r;
INT ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = r;
}

DOUBLE max_abs_value_e(tGrid,ze)  /* max := max(z_e) */
GRID *tGrid;
INT ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE max=0.;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      max = MAX(max,fabs(ED(pelem,ze)));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            max = MAX(max,fabs(ED(pelem,ze)));
   return(max);
}

void add_value_e(tGrid,r,ze)  /* ze_i := ze_i + r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) += r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) += r;
}

DOUBLE sum_e(tGrid,ze)
GRID *tGrid;
INT ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.; 
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      sum += ED(pelem,ze);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += ED(pelem,ze);
   return(sum);
}

DOUBLE sum_e_n(tGrid,ze,n)
GRID *tGrid;
INT ze, *n;
{
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.; 
   
   *n = 0;
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      sum += ED(pelem,ze);
      (*n)++;
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            sum += ED(pelem,ze);
            (*n)++;
         }
   return(sum);
}

DOUBLE mean_value_e(tGrid,ze)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   FLOAT sum=0., vol=0., evol;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      evol = VOLUME(pelem);
      sum += ED(pelem,ze)*evol;
      vol += evol;
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            evol = VOLUME(pelem);
            sum += ED(pelem,ze)*evol;
            vol += evol;
         }
   return(sum/vol);
}

void set_zero_mean_e(tGrid,ze)
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE mean;
   
   mean = mean_value_e(tGrid,ze);

   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) -= mean;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) -= mean;
}

void set_first_element(tGrid,r,ze)  /* z_1 := r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   ED(FIRSTELEMENT(tGrid),ze) = r;
}

void subtr_first_element(tGrid,ze)  /* ze_i := ze_i - ze_1 */
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   FLOAT fel;
   
   fel = ED(FIRSTELEMENT(tGrid),ze);
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) -= fel;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) -= fel;
}

void copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe);
}

void inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = -ED(pelem,xe);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = -ED(pelem,xe);
}

void add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe) + ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe) + ED(pelem,ye);
}

void subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe) - ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe) - ED(pelem,ye);
}

void add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid;
INT we,xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,we) + ED(pelem,xe) - ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,we) + ED(pelem,xe) - ED(pelem,ye);
}

void mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = r*ED(pelem,xe) + ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = r*ED(pelem,xe) + ED(pelem,ye);
}

void mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = r*ED(pelem,xe) - ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = r*ED(pelem,xe) - ED(pelem,ye);
}

void damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe) + r*(ED(pelem,ye) - ED(pelem,xe));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe) + r*(ED(pelem,ye) - ED(pelem,xe));
}

void mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid;
FLOAT r;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe)*r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe)*r;
}

void mult2_e(tGrid,r,xe,ze,t)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze,t;
{  mult_e(tGrid,r,xe,ze);  }

void add_sum_of_multiples_e(tGrid,r,xe,n,ze)
GRID *tGrid;  /* ze += sum_{i=0}^{n-1} r[i]*(xe+i) */
DOUBLE *r;
INT xe, n, ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum;
   INT i;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      sum = 0.;
      for (i = 0; i < n; i++)
         sum += ED(pelem,xe+i)*r[i];
      ED(pelem,ze) += sum;
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            sum = 0.;
            for (i = 0; i < n; i++)
               sum += ED(pelem,xe+i)*r[i];
            ED(pelem,ze) += sum;
         }
}

DOUBLE dot_e_with_first(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.0;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      sum += ED(pelem,xe)*ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += ED(pelem,xe)*ED(pelem,ye);
   return(sum);
}

DOUBLE dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.0;
   
   for (pelem = SUCC(FIRSTELEMENT(tGrid)); pelem!=NULL; pelem=pelem->succ)
      sum += ED(pelem,xe)*ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += ED(pelem,xe)*ED(pelem,ye);
   return(sum);
}

DOUBLE dot2_e(tGrid,xe,ye,t)  /* := xe*ye */
GRID *tGrid; INT xe,ye,t;
{  return(dot_e(tGrid,xe,ye));  }

void divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe)/ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe)/ED(pelem,ye);
}

void multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = ED(pelem,xe)*ED(pelem,ye);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = ED(pelem,xe)*ED(pelem,ye);
}

void make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ED(pelem,ze) = sqrt(ED(pelem,xe));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ED(pelem,ze) = sqrt(ED(pelem,xe));
}

void mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         ED(pel,l) -=  h[i][j]*ED(pel,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               ED(pel,l) -= h[i][j]*ED(pel,k+i);
}

void mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid;
FLOAT g[GMN];
INT xe, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         ED(pel,xe) += g[i]*ED(pel,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)   
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               ED(pel,xe) += g[i]*ED(pel,k+i);
}

void save_e_vector(tGrid,z,m,fp)
GRID *tGrid;
INT z, m;
FILE *fp;
{
   ELEMENT *pel;

   for (pel=FIRSTELEMENT(tGrid); pel; pel=SUCC(pel))
      fprintf(fp,"%i %e\n",m+pel->eflag,ED(pel,z));
}

#else  /*  if !(E_DATA & SCALAR_ELEMENT_DATA)  */

void set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: set_value_e not available.\n");  }

DOUBLE max_abs_value_e(tGrid,ze)  /* max := max(z_e) */
GRID *tGrid; INT ze;
{  eprintf("Error: max_abs_value_e not available.\n"); return(0.);  }

void add_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: add_value_e not available.\n");  }

DOUBLE sum_e(tGrid,ze)
GRID *tGrid; INT ze;
{  eprintf("Error: sum_e not available.\n"); return(0.);  }

DOUBLE sum_e_n(tGrid,ze,n)
GRID *tGrid; INT ze, *n;
{  eprintf("Error: sum_e_n not available.\n"); return(0.);  }

DOUBLE mean_value_e(tGrid,ze)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid; INT ze;
{  eprintf("Error: mean_value_e not available.\n"); return(0.);  }
 
void set_zero_mean_e(tGrid,ze)
GRID *tGrid; INT ze;
{  eprintf("Error: set_zero_mean_e not available.\n");  }

void set_first_element(tGrid,r,ze)  /* z_1 := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: set_first_element not available.\n");  }

void subtr_first_element(tGrid,ze)  /* ze_i := ze_i - ze_1 */
GRID *tGrid; INT ze;
{  eprintf("Error: subtr_first_element not available.\n");  }

void copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: copy_e not available.\n");  }

void inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: inv_e not available.\n");  }

void add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: add_e not available.\n");  }

void subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: subtr_e not available.\n");  }

void add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid; INT we,xe,ye,ze;
{  eprintf("Error: add_and_subtr_e not available.\n");  }

void mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: mult_and_add_e not available.\n");  }

void mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: mult_and_subtr_e not available.\n");  }

void damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: damp_e not available.\n");  }

void mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze;
{  eprintf("Error: mult_e not available.\n");  }

void mult2_e(tGrid,r,xe,ze,t)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze,t;
{  eprintf("Error: mult2_e not available.\n");  }

void add_sum_of_multiples_e(tGrid,r,xe,n,ze)
GRID *tGrid; DOUBLE *r; INT xe, n, ze;
{  eprintf("Error: add_sum_of_multiples_e not available.\n");  }

DOUBLE dot_e_with_first(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: dot_e_with_first not available.\n"); return(0.);  }

DOUBLE dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: dot_e not available.\n"); return(0.);  }

DOUBLE dot2_e(tGrid,xe,ye,t)  /* := xe*ye */
GRID *tGrid; INT xe,ye,t;
{  eprintf("Error: dot2_e not available.\n"); return(0.);  }

void divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: divide_e not available.\n");  }

void multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: multiply_e not available.\n");  }

void make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: make_sqrt_e not available.\n");  }

void mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k;
{  eprintf("Error: mult_and_add1_e_for_GMRES not available.\n");  }

void mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid; FLOAT g[GMN]; INT xe, j, k;
{  eprintf("Error: mult_and_add2_e_for_GMRES not available.\n");  }

void save_e_vector(tGrid,z,m,fp)
GRID *tGrid; INT z, m; FILE *fp;
{  eprintf("Error: save_e_vector not available.\n");  }

#endif  /*  if !(E_DATA & SCALAR_ELEMENT_DATA)  */

#if N_DATA & VECTOR_NODE_DATA

void v_set_value_n(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      SET_VALUE(theNode,r,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            SET_VALUE(theNode,r,z)
   }
}

void v_copy_n(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      COPY(theNode,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            COPY(theNode,x,z)
   }
}

void v_inv_n(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      INV(theNode,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            INV(theNode,x,z)
}

void v_add_n(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      ADD(theNode,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            ADD(theNode,x,y,z)
}

void v_subtr_n(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      NSUBTR(theNode,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NSUBTR(theNode,x,y,z)
   }
}

void v_add_and_subtr_n(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      NASUBTR(theNode,w,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NASUBTR(theNode,w,x,y,z)
}

void v_mult_and_add_n(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MADD(theNode,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MADD(theNode,r,x,y,z)
}

void v_mult_and_subtr_n(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MSUBTR(theNode,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MSUBTR(theNode,r,x,y,z)
}

void v_damp_n(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DAMP(theNode,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            DAMP(theNode,r,x,y,z)
}

void v_mult_n(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULT(theNode,r,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MULT(theNode,r,x,z)
}

DOUBLE v_dot_n(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{
   GRID *theGrid;
   NODE *theNode;
   DOUBLE sum=0.0;
	
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      sum += DOT(NDD(theNode,x),NDD(theNode,y));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            sum += DOT(NDD(theNode,x),NDD(theNode,y));
   return(sum);
}

void v_divide_n(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DIVIDE(theNode,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            DIVIDE(theNode,x,y,z)
}

void v_multiply_n(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULTIPLY(theNode,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MULTIPLY(theNode,x,y,z)
}

void v_make_sqrt_n(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MSQRT(theNode,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MSQRT(theNode,x,z)
}

void v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   NODE *theNode;
   INT i;
   
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,(-h[i][j]),k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(theNode,(-h[i][j]),k+i,l)
}

void v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i;

   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,g[i],k+i,x)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(theNode,g[i],k+i,x)
}

void save_vn_vector(tGrid,z,m,ni,fp)
GRID *tGrid;
INT z, m, ni;
FILE *fp;
{
   NODE *theNode;

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      fprintf(fp,"%i %e\n",m+theNode->index,ND(theNode,z,0));
   m += ni;
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      fprintf(fp,"%i %e\n",m+theNode->index,ND(theNode,z,1));
   m += ni;
   if (DIM == 3)
      for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
         fprintf(fp,"%i %e\n",m+theNode->index,ND(theNode,z,2));
}

#else  /* if !(N_DATA & VECTOR_NODE_DATA)  */

void v_set_value_n(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: v_set_value_n not available.\n");  }

void v_copy_n(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: v_copy_n not available.\n");  }

void v_inv_n(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: v_inv_n not available.\n");  }

void v_add_n(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: v_add_n not available.\n");  }

void v_subtr_n(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: v_subtr_n not available.\n");  }

void v_add_and_subtr_n(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: v_add_and_subtr_n not available.\n");  }

void v_mult_and_add_n(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: v_mult_and_add_n not available.\n");  }

void v_mult_and_subtr_n(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: v_mult_and_subtr_n not available.\n");  }

void v_damp_n(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: v_damp_n not available.\n");  }

void v_mult_n(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: v_mult_n not available.\n");  }

DOUBLE v_dot_n(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: v_dot_n not available.\n"); return(0.);  }

void v_divide_n(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: v_divide_n not available.\n");  }

void v_multiply_n(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: v_multiply_n not available.\n");  }

void v_make_sqrt_n(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: v_make_sqrt_n not available.\n");  }

void v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: v_mult_and_add1_for_GMRES_n not available.\n");  }

void v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: v_mult_and_add2_for_GMRES_n not available.\n");  }

void save_vn_vector(tGrid,z,m,ni,fp)
GRID *tGrid; INT z, m, ni; FILE *fp;
{  eprintf("Error: save_vn_vector not available.\n");  }

#endif  /* if !(N_DATA & VECTOR_NODE_DATA)  */

#if N_DATA & MVECTOR_NODE_DATA

void mv_set_value_n(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   NODE *theNode, *stop;
   INT i, kk=N_OF_NODE_FUNC;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      GSET7(NDMVP(theNode,z),r,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET7(NDMVP(theNode,z),r,i,kk)
   }
}

void mv_copy_n(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode, *stop;
   INT i, kk=N_OF_NODE_FUNC;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      GSET1(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET1(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
   }
}

void mv_inv_n(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET8(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET8(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
}

void mv_add_n(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET10(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET10(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
}

void mv_subtr_n(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode, *stop;
   INT i, kk=N_OF_NODE_FUNC;
   	
   stop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != stop; theNode=SUCC(theNode))
      GSET11(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != stop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET11(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   }
}

void mv_add_and_subtr_n(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET30(NDMVP(theNode,z),
             NDMVP(theNode,w),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET30(NDMVP(theNode,z),
                   NDMVP(theNode,w),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
}

void mv_mult_and_add_n(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET23(NDMVP(theNode,z),NDMVP(theNode,y),NDMVP(theNode,x),r,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET23(NDMVP(theNode,z),NDMVP(theNode,y),NDMVP(theNode,x),r,i,kk)
}

void mv_mult_and_subtr_n(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET24(NDMVP(theNode,z),NDMVP(theNode,y),NDMVP(theNode,x),r,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET24(NDMVP(theNode,z),NDMVP(theNode,y),NDMVP(theNode,x),r,i,kk)
}

void mv_damp_n(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET31(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),r,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET31(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),r,i,kk)
}

void mv_mult_n(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET2(NDMVP(theNode,z),NDMVP(theNode,x),r,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET2(NDMVP(theNode,z),NDMVP(theNode,x),r,i,kk)
}

DOUBLE mv_dot_n(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{
   GRID *theGrid;
   NODE *theNode;
   DOUBLE sum=0.0;
   INT i, kk=N_OF_NODE_FUNC;
	
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GDOTS(NDMVP(theNode,x),NDMVP(theNode,y),sum,i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GDOTS(NDMVP(theNode,x),NDMVP(theNode,y),sum,i,kk)
   return(sum);
}

void mv_divide_n(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET32(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET32(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
}

void mv_multiply_n(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET33(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET33(NDMVP(theNode,z),NDMVP(theNode,x),NDMVP(theNode,y),i,kk)
}

void mv_make_sqrt_n(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i, kk=N_OF_NODE_FUNC;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      GSET34(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            GSET34(NDMVP(theNode,z),NDMVP(theNode,x),i,kk)
}

void mv_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   NODE *theNode;
   INT i, m, kk=N_OF_NODE_FUNC;
   
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         GSET4(NDMVP(theNode,l),NDMVP(theNode,k+i),(-h[i][j]),m,kk)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               GSET4(NDMVP(theNode,l),NDMVP(theNode,k+i),(-h[i][j]),m,kk)
}

void mv_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{
   GRID *theGrid;
   NODE *theNode;
   INT i, m, kk=N_OF_NODE_FUNC;

   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         GSET4(NDMVP(theNode,x),NDMVP(theNode,k+i),g[i],m,kk)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               GSET4(NDMVP(theNode,x),NDMVP(theNode,k+i),g[i],m,kk)
}

#else  /* if !(N_DATA & MVECTOR_NODE_DATA)  */

void mv_set_value_n(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: mv_set_value_n not available.\n");  }

void mv_copy_n(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: mv_copy_n not available.\n");  }

void mv_inv_n(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: mv_inv_n not available.\n");  }

void mv_add_n(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: mv_add_n not available.\n");  }

void mv_subtr_n(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: mv_subtr_n not available.\n");  }

void mv_add_and_subtr_n(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: mv_add_and_subtr_n not available.\n");  }

void mv_mult_and_add_n(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_mult_and_add_n not available.\n");  }

void mv_mult_and_subtr_n(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_mult_and_subtr_n not available.\n");  }

void mv_damp_n(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_damp_n not available.\n");  }

void mv_mult_n(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: mv_mult_n not available.\n");  }

DOUBLE mv_dot_n(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: mv_dot_n not available.\n"); return(0.);  }

void mv_divide_n(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: mv_divide_n not available.\n");  }

void mv_multiply_n(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: mv_multiply_n not available.\n");  }

void mv_make_sqrt_n(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: mv_make_sqrt_n not available.\n");  }

void mv_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: mv_mult_and_add1_for_GMRES_n not available.\n");  }

void mv_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: mv_mult_and_add2_for_GMRES_n not available.\n");  }

#endif  /* if !(N_DATA & MVECTOR_NODE_DATA)  */

#if F_DATA & VECTOR_FACE_DATA

void vset_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      SET_VALUE(theFace,r,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            SET_VALUE(theFace,r,z)
   }
}

void vcopy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      COPY(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            COPY(theFace,x,z)
   }
}

void vinv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      INV(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            INV(theFace,x,z)
}

void vadd_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      ADD(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            ADD(theFace,x,y,z)
}

void vsubtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      NSUBTR(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            NSUBTR(theFace,x,y,z)
   }
}

void vadd_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      NASUBTR(theFace,w,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            NASUBTR(theFace,w,x,y,z)
}

void vmult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      MADD(theFace,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MADD(theFace,r,x,y,z)
}

void vmult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      MSUBTR(theFace,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MSUBTR(theFace,r,x,y,z)
}

void vdamp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      DAMP(theFace,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            DAMP(theFace,r,x,y,z)
}

void vmult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      MULT(theFace,r,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MULT(theFace,r,x,z)
}

DOUBLE vdot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      sum += DOT(FDVP(theFace,x),FDVP(theFace,y));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            sum += DOT(FDVP(theFace,x),FDVP(theFace,y));
   return(sum);
}

void vdivide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      DIVIDE(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            DIVIDE(theFace,x,y,z)
}

void vmultiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      MULTIPLY(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MULTIPLY(theFace,x,y,z)
}

void vmake_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      MSQRT(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            MSQRT(theFace,x,z)
}

void vmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;
   
   for (pface = FIRST_FACE(tGrid,t); pface!=NULL; pface=pface->succ)
      for (i = 0; i <= j; i++)
         SUBTR_MULT(pface,h[i][j],k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               SUBTR_MULT(pface,h[i][j],k+i,l)
}

void vmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;

   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(pface,g[i],k+i,x)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(pface,g[i],k+i,x)
}

void save_vf_vector(tGrid,z,m,ni,fp)
GRID *tGrid;
INT z, m, ni;
FILE *fp;
{
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDV(theFace,z,0));
   m += ni;
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDV(theFace,z,1));
   m += ni;
   if (DIM == 3)
      for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
         fprintf(fp,"%i %e\n",m+theFace->index,FDV(theFace,z,2));
}

#else  /* if !(F_DATA & VECTOR_FACE_DATA)  */

void vset_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: vset_value_f not available.\n");  }

void vcopy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vcopy_f not available.\n");  }

void vinv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vinv_f not available.\n");  }

void vadd_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vadd_f not available.\n");  }

void vsubtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vsubtr_f not available.\n");  }

void vadd_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: vadd_and_subtr_f not available.\n");  }

void vmult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vmult_and_add_f not available.\n");  }

void vmult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vmult_and_subtr_f not available.\n");  }

void vdamp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vdamp_f not available.\n");  }

void vmult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: vmult_f not available.\n");  }

DOUBLE vdot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: vdot_f not available.\n"); return(0.);  }

void vdivide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: vdivide_f not available.\n");  }

void vmultiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: vmultiply_f not available.\n");  }

void vmake_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: vmake_sqrt_f not available.\n");  }

void vmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: vmult_and_add1_for_GMRES_f not available.\n");  }

void vmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: vmult_and_add2_for_GMRES_f not available.\n");  }

void save_vf_vector(tGrid,z,m,ni,fp)
GRID *tGrid; INT z, m, ni; FILE *fp;
{  eprintf("Error: save_vf_vector not available.\n");  }

#endif  /*  !(F_DATA & VECTOR_FACE_DATA)  */

#if F_DATA & DVECTOR_FACE_DATA

void dvset_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      D_SET_VALUE(theFace,r,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_SET_VALUE(theFace,r,z)
   }
}

void dvcopy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      D_COPY(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_COPY(theFace,x,z)
   }
}

void dvinv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_INV(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_INV(theFace,x,z)
}

void dvadd_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_ADD(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_ADD(theFace,x,y,z)
}

void dvsubtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace, *stop;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      D_NSUBTR(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_NSUBTR(theFace,x,y,z)
   }
}

void dvadd_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_NASUBTR(theFace,w,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_NASUBTR(theFace,w,x,y,z)
}

void dvmult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_MADD(theFace,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_MADD(theFace,r,x,y,z)
}

void dvmult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: dvmult_and_subtr_f not available.\n");  }

void dvdamp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_DAMP(theFace,r,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_DAMP(theFace,r,x,y,z)
}

void dvmult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_MULT(theFace,r,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_MULT(theFace,r,x,z)
}

DOUBLE dvdot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      sum += DOT(FDDVP(theFace,x,0),FDDVP(theFace,y,0)) +
             DOT(FDDVP(theFace,x,1),FDDVP(theFace,y,1));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            sum += DOT(FDDVP(theFace,x,0),FDDVP(theFace,y,0)) +
                   DOT(FDDVP(theFace,x,1),FDDVP(theFace,y,1));
   return(sum);
}

void dvdivide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_DIVIDE(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_DIVIDE(theFace,x,y,z)
}

void dvmultiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_MULTIPLY(theFace,x,y,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_MULTIPLY(theFace,x,y,z)
}

void dvmake_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      D_MSQRT(theFace,x,z)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            D_MSQRT(theFace,x,z)
}

void dvmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;
   
   for (pface = FIRST_FACE(tGrid,t); pface!=NULL; pface=pface->succ)
      for (i = 0; i <= j; i++)
         D_SUBTR_MULT(pface,h[i][j],k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               D_SUBTR_MULT(pface,h[i][j],k+i,l)
}

void dvmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i;

   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         D_ADD_MULT(pface,g[i],k+i,x)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               D_ADD_MULT(pface,g[i],k+i,x)
}

void save_dvf_vector(tGrid,z,m,ni,fp)
GRID *tGrid;
INT z, m, ni;
FILE *fp;
{
   FACE *theFace;

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,0,0));
   m += ni;
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,1,0));
   m += ni;
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,0,1));
   m += ni;
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,1,1));
   if (DIM == 3){
      m += ni;
      for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
         fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,0,2));
      m += ni;
      for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
         fprintf(fp,"%i %e\n",m+theFace->index,FDDV(theFace,z,1,2));
   }
}

#else  /*  if !(F_DATA & DVECTOR_FACE_DATA)  */

void dvset_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: dvset_value_f not available.\n");  }

void dvcopy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: dvcopy_f not available.\n");  }

void dvinv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: dvinv_f not available.\n");  }

void dvadd_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: dvadd_f not available.\n");  }

void dvsubtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: dvsubtr_f not available.\n");  }

void dvadd_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: dvadd_and_subtr_f not available.\n");  }

void dvmult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: dvmult_and_add_f not available.\n");  }

void dvmult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: dvmult_and_subtr_f not available.\n");  }

void dvdamp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: dvdamp_f not available.\n");  }

void dvmult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: dvmult_f not available.\n");  }

DOUBLE dvdot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: dvdot_f not available.\n"); return(0.);  }

void dvdivide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: dvdivide_f not available.\n");  }

void dvmultiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: dvmultiply_f not available.\n");  }

void dvmake_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: dvmake_sqrt_f not available.\n");  }

void dvmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: dvmult_and_add1_for_GMRES_f not available.\n");  }

void dvmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: dvmult_and_add2_for_GMRES_f not available.\n");  }

void save_dvf_vector(tGrid,z,m,ni,fp)
GRID *tGrid; INT z, m, ni; FILE *fp;
{  eprintf("Error: save_dvf_vector not available.\n");  }

#endif  /*  !(F_DATA & DVECTOR_FACE_DATA)  */

#if F_DATA & MVECTOR_FACE_DATA

void mv_set_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   INT i, nn=N_OF_FACE_FUNC;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      GSET7(FDMVP(theFace,z),r,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET7(FDMVP(theFace,z),r,i,nn)
   }
}

void mv_copy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace, *stop;
   INT i, nn=N_OF_FACE_FUNC;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      GSET1(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET1(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
   }
}

void mv_inv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET8(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET8(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
}

void mv_add_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET10(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET10(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
}

void mv_subtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace, *stop;
   INT i, nn=N_OF_FACE_FUNC;
   	
   stop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != stop; theFace=SUCC(theFace))
      GSET11(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      stop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != stop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET11(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   }
}

void mv_add_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET30(FDMVP(theFace,z),
             FDMVP(theFace,w),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET30(FDMVP(theFace,z),
                   FDMVP(theFace,w),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
}

void mv_mult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET23(FDMVP(theFace,z),FDMVP(theFace,y),FDMVP(theFace,x),r,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET23(FDMVP(theFace,z),FDMVP(theFace,y),FDMVP(theFace,x),r,i,nn)
}

void mv_mult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET24(FDMVP(theFace,z),FDMVP(theFace,y),FDMVP(theFace,x),r,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET24(FDMVP(theFace,z),FDMVP(theFace,y),FDMVP(theFace,x),r,i,nn)
}

void mv_damp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET31(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),r,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET31(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),r,i,nn)
}

void mv_mult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET2(FDMVP(theFace,z),FDMVP(theFace,x),r,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET2(FDMVP(theFace,z),FDMVP(theFace,x),r,i,nn)
}

DOUBLE mv_dot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   GRID *theGrid;
   FACE *theFace;
   DOUBLE sum=0.0;
   INT i, nn=N_OF_FACE_FUNC;
	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GDOTS(FDMVP(theFace,x),FDMVP(theFace,y),sum,i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GDOTS(FDMVP(theFace,x),FDMVP(theFace,y),sum,i,nn)
   return(sum);
}

void mv_divide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET32(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET32(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
}

void mv_multiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET33(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET33(FDMVP(theFace,z),FDMVP(theFace,x),FDMVP(theFace,y),i,nn)
}

void mv_make_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   FACE *theFace;
   INT i, nn=N_OF_FACE_FUNC;
   	
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      GSET34(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            GSET34(FDMVP(theFace,z),FDMVP(theFace,x),i,nn)
}

void mv_mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i, m, nn=N_OF_FACE_FUNC;
   
   for (pface = FIRST_FACE(tGrid,t); pface!=NULL; pface=pface->succ)
      for (i = 0; i <= j; i++)
         GSET9(FDMVP(pface,l),FDMVP(pface,k+i),h[i][j],m,nn)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               GSET9(FDMVP(pface,l),FDMVP(pface,k+i),h[i][j],m,nn)
}

void mv_mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{   
   GRID *theGrid;
   FACE *pface;
   INT i, m, nn=N_OF_FACE_FUNC;

   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         GSET4(FDMVP(pface,x),FDMVP(pface,k+i),g[i],m,nn)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               GSET4(FDMVP(pface,x),FDMVP(pface,k+i),g[i],m,nn)
}

#else  /* if !(F_DATA & MVECTOR_FACE_DATA)  */

void mv_set_value_f(tGrid,r,z,t)  /* zf := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: mv_set_value_f not available.\n");  }

void mv_copy_f(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: mv_copy_f not available.\n");  }

void mv_inv_f(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: mv_inv_f not available.\n");  }

void mv_add_f(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: mv_add_f not available.\n");  }

void mv_subtr_f(tGrid,x,y,z,t)  /* zf := xf - yf */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: mv_subtr_f not available.\n");  }

void mv_add_and_subtr_f(tGrid,w,x,y,z,t)  /* zf := wf + xf - yf */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: mv_add_and_subtr_f not available.\n");  }

void mv_mult_and_add_f(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_mult_and_add_f not available.\n");  }

void mv_mult_and_subtr_f(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_mult_and_subtr_f not available.\n");  }

void mv_damp_f(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: mv_damp_f not available.\n");  }

void mv_mult_f(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: mv_mult_f not available.\n");  }

DOUBLE mv_dot_f(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: mv_dot_f not available.\n"); return(0.);  }

void mv_divide_f(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: mv_divide_f not available.\n");  }

void mv_multiply_f(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: mv_multiply_f not available.\n");  }

void mv_make_sqrt_f(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: mv_make_sqrt_f not available.\n");  }

void mv_mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: mv_mult_and_add1_for_GMRES_f not available.\n");  }

void mv_mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: mv_mult_and_add2_for_GMRES_f not available.\n");  }

#endif  /*  !(F_DATA & MVECTOR_FACE_DATA)  */

#if E_DATA & VECTOR_ELEMENT_DATA

void vset_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      SET_VALUE(pelem,r,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            SET_VALUE(pelem,r,ze)
}

void vcopy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      COPY(pelem,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            COPY(pelem,xe,ze)
}

void vinv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      INV(pelem,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            INV(pelem,xe,ze)
}

void vadd_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      ADD(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            ADD(pelem,xe,ye,ze)
}

void vsubtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid;
INT xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      NSUBTR(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            NSUBTR(pelem,xe,ye,ze)
}

void vadd_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid;
INT we,xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      NASUBTR(pelem,we,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            NASUBTR(pelem,we,xe,ye,ze)
}

void vmult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      MADD(pelem,r,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            MADD(pelem,r,xe,ye,ze)
}

void vmult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      MSUBTR(pelem,r,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            MSUBTR(pelem,r,xe,ye,ze)
}

void vdamp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      DAMP(pel,r,xe,ye,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            DAMP(pel,r,xe,ye,ze)
}

void vmult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid;
FLOAT r;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      MULT(pelem,r,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            MULT(pelem,r,xe,ze)
}

void vmult_e_s(tGrid,r,xe,ze)  /* ze_i := r_i*xe_i, i=1,...,DIM */
GRID *tGrid;
FLOAT r[DIM];
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      MULT_S(pelem,r,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            MULT_S(pelem,r,xe,ze)
}

void vadd_sum_of_multiples_e(tGrid,r,xe,n,ze)
GRID *tGrid;  /* ze += sum_{i=0}^{n-1} r[i]*(xe+i) */
DOUBLE *r;
INT xe, n, ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   INT i;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      for (i = 0; i < n; i++)
         ADD_MULT(pelem,r[i],xe+i,ze)
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            for (i = 0; i < n; i++)
               ADD_MULT(pelem,r[i],xe+i,ze)
         }
}


DOUBLE vdot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.0;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      sum += DOT(NDD(pelem,xe),NDD(pelem,ye));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += DOT(NDD(pelem,xe),NDD(pelem,ye)); 
   return(sum);
}

void vdivide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vdivide_e not available.\n");  }

void vmultiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vmultiply_e not available.\n");  }
   
void vmake_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: vmake_sqrt_e not available.\n");  }

void vmult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(pel,(-h[i][j]),k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(pel,(-h[i][j]),k+i,l)
}

void vmult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid;
FLOAT g[GMN];
INT xe, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(pel,g[i],k+i,xe)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser) 
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(pel,g[i],k+i,xe)
}

void save_ve_vector(tGrid,z,m,ne,fp)
GRID *tGrid;
INT z, m, ne;
FILE *fp;
{
   ELEMENT *pelem;

   for (pelem=FIRSTELEMENT(tGrid); pelem; pelem=SUCC(pelem))
      fprintf(fp,"%i %e\n",m+pelem->eflag,EDV(pelem,z,0));
   m += ne;
   for (pelem=FIRSTELEMENT(tGrid); pelem; pelem=SUCC(pelem))
      fprintf(fp,"%i %e\n",m+pelem->eflag,EDV(pelem,z,1));
   m += ne;
   if (DIM == 3)
      for (pelem=FIRSTELEMENT(tGrid); pelem; pelem=SUCC(pelem))
         fprintf(fp,"%i %e\n",m+pelem->eflag,EDV(pelem,z,2));
}

#else  /*  if !(E_DATA & VECTOR_ELEMENT_DATA)  */

void vset_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: vset_value_e not available.\n");  }

void vcopy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: vcopy_e not available.\n");  }

void vinv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: vinv_e not available.\n");  }

void vadd_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vadd_e not available.\n");  }

void vsubtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vsubtr_e not available.\n");  }

void vadd_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid; INT we,xe,ye,ze;
{  eprintf("Error: vadd_and_subtr_e not available.\n");  }

void vmult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: vmult_and_add_e not available.\n");  }

void vmult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: vmult_and_subtr_e not available.\n");  }

void vdamp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: vdamp_e not available.\n");  }

void vmult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze;
{  eprintf("Error: vmult_e not available.\n");  }

void vmult_e_s(tGrid,r,xe,ze)  /* ze_i := r_i*xe_i, i=1,...,DIM */
GRID *tGrid; FLOAT r[DIM]; INT xe,ze;
{  eprintf("Error: vmult_e_s not available.\n");  }

void vadd_sum_of_multiples_e(tGrid,r,xe,n,ze)
GRID *tGrid; DOUBLE *r; INT xe, n, ze;
{  eprintf("Error: vadd_sum_of_multiples_e not available.\n");  }

DOUBLE vdot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: vdot_e not available.\n"); return(0.);  }

void vdivide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vdivide_e not available.\n");  }

void vmultiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vmultiply_e not available.\n");  }
   
void vmake_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: vmake_sqrt_e not available.\n");  }

void vmult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k;
{  eprintf("Error: vmult_and_add1_e_for_GMRES not available.\n");  }

void vmult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid; FLOAT g[GMN]; INT xe, j, k;
{  eprintf("Error: vmult_and_add2_e_for_GMRES not available.\n");  }

void save_ve_vector(tGrid,z,m,ne,fp)
GRID *tGrid; INT z, m, ne; FILE *fp;
{  eprintf("Error: save_ve_vector not available.\n");  }

#endif  /*  !(E_DATA & VECTOR_ELEMENT_DATA)  */

#if E_DATA & MVECTOR_ELEMENT_DATA

void mv_set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET7(EDMVP(pelem,ze),r,i,mm)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET7(EDMVP(pelem,ze),r,i,mm)
}

void mv_set_value_es(tGrid,r,ze)
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem=pelem->succ)
      for(i = 1; i < mm; i++)
         EDMV(pelem,ze,i) = r;
}

void mv_copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET1(EDMVP(pelem,ze),EDMVP(pelem,xe),i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET1(EDMVP(pelem,ze),EDMVP(pelem,xe),i,mm)
}

void mv_inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET8(EDMVP(pelem,ze),EDMVP(pelem,xe),i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET8(EDMVP(pelem,ze),EDMVP(pelem,xe),i,mm)
}

void mv_add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET10(EDMVP(pelem,ze),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET10(EDMVP(pelem,ze),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
}

void mv_subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid;
INT xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET11(EDMVP(pelem,ze),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET11(EDMVP(pelem,ze),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
}

void mv_add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid;
INT we,xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET30(EDMVP(pelem,ze),
             EDMVP(pelem,we),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET30(EDMVP(pelem,ze),
                   EDMVP(pelem,we),EDMVP(pelem,xe),EDMVP(pelem,ye),i,mm)
}

void mv_mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET23(EDMVP(pelem,ze),EDMVP(pelem,ye),EDMVP(pelem,xe),r,i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET23(EDMVP(pelem,ze),EDMVP(pelem,ye),EDMVP(pelem,xe),r,i,mm)
}

void mv_mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET24(EDMVP(pelem,ze),EDMVP(pelem,ye),EDMVP(pelem,xe),r,i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET24(EDMVP(pelem,ze),EDMVP(pelem,ye),EDMVP(pelem,xe),r,i,mm)
}

void mv_damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      GSET31(EDMVP(pel,ze),EDMVP(pel,xe),EDMVP(pel,ye),r,i,mm)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            GSET31(EDMVP(pel,ze),EDMVP(pel,xe),EDMVP(pel,ye),r,i,mm)
}

void mv_mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid;
FLOAT r;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GSET2(EDMVP(pelem,ze),EDMVP(pelem,xe),r,i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GSET2(EDMVP(pelem,ze),EDMVP(pelem,xe),r,i,mm)
}

DOUBLE mv_dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.0;
   INT i, mm=N_OF_ELEM_FUNC;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      GDOTS(EDMVP(pelem,xe),EDMVP(pelem,ye),sum,i,mm)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            GDOTS(EDMVP(pelem,xe),EDMVP(pelem,ye),sum,i,mm)
   return(sum);
}

void mv_divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_divide_e not available.\n");  }

void mv_multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_multiply_e not available.\n");  }
   
void mv_make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: mv_make_sqrt_e not available.\n");  }

void mv_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i, m, mm=N_OF_ELEM_FUNC;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         GSET4(EDMVP(pel,l),EDMVP(pel,k+i),(-h[i][j]),m,mm)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               GSET4(EDMVP(pel,l),EDMVP(pel,k+i),(-h[i][j]),m,mm)
}

void mv_mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid;
FLOAT g[GMN];
INT xe, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i, m, mm=N_OF_ELEM_FUNC;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         GSET4(EDMVP(pel,xe),EDMVP(pel,k+i),g[i],m,mm)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser) 
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               GSET4(EDMVP(pel,xe),EDMVP(pel,k+i),g[i],m,mm)
}

#else  /*  if !(E_DATA & MVECTOR_ELEMENT_DATA)  */

void mv_set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: mv_set_value_e not available.\n");  }

void mv_set_value_es(tGrid,r,ze)
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: mv_set_value_es not available.\n");  }

void mv_copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: mv_copy_e not available.\n");  }

void mv_inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: mv_inv_e not available.\n");  }

void mv_add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_add_e not available.\n");  }

void mv_subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_subtr_e not available.\n");  }

void mv_add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid; INT we,xe,ye,ze;
{  eprintf("Error: mv_add_and_subtr_e not available.\n");  }

void mv_mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: mv_mult_and_add_e not available.\n");  }

void mv_mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: mv_mult_and_subtr_e not available.\n");  }

void mv_damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: mv_damp_e not available.\n");  }

void mv_mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze;
{  eprintf("Error: mv_mult_e not available.\n");  }

DOUBLE mv_dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: mv_dot_e not available.\n"); return(0.);  }

void mv_divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_divide_e not available.\n");  }

void mv_multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: mv_multiply_e not available.\n");  }
   
void mv_make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: mv_make_sqrt_e not available.\n");  }

void mv_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k;
{  eprintf("Error: mv_mult_and_add1_e_for_GMRES not available.\n");  }

void mv_mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid; FLOAT g[GMN]; INT xe, j, k;
{  eprintf("Error: mv_mult_and_add2_e_for_GMRES not available.\n");  }

#endif  /*  !(E_DATA & MVECTOR_ELEMENT_DATA)  */

#if (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)

void set_vect_to_zero_for_q2b3(tGrid,u)
GRID *tGrid;
INT u;
{
   FACE *theFace;
   ELEMENT *pel;
   INT k;

   for (theFace=FIRSTF(tGrid); theFace; theFace=SUCC(theFace))
      FDMV(theFace,u,1) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      EDMV(pel,u,1) = EDMV(pel,u,2) = EDMV(pel,u,3) = 0.;
}

void partial_set_vect_to_zero_for_q2b3(tGrid,u)
GRID *tGrid;
INT u;
{
   NODE *n1, *n2;
   ELEMENT *pel;
   INT k;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
   if (IS_B_EL(pel)){
      for (k=0; k < NVERT; k++){
         n1 = pel->n[k];
         n2 = pel->n[(k+1)%NVERT];
         if (IS_BF(pel->f[k]) || (IS_FN(n1) && IS_FN(n2)))
            FDMV(pel->f[k],u,1) = 0.;
      }
   }
   else
   {
      EDMV(pel,u,1) = EDMV(pel,u,2) = EDMV(pel,u,3) = 0.;
      for (k=0; k < NVERT; k++)
         FDMV(pel->f[k],u,1) = 0.;
   }
}

#else

void set_vect_to_zero_for_q2b3(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: set_vect_to_zero_for_q2b3 not available.\n");  }

void partial_set_vect_to_zero_for_q2b3(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: partial_set_vect_to_zero_for_q2b3 not available.\n");  }

#endif  /* !((F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)) */

#if (E_DATA & SCALAR_ELEMENT_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

void vdivide_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      SDIVIDE(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            SDIVIDE(pelem,xe,ye,ze)
}

void vmultiply_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      SMULTIPLY(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            SMULTIPLY(pelem,xe,ye,ze)
}

#else

void vdivide_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vdivide_1e not available.\n");  }

void vmultiply_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: vmultiply_1e not available.\n");  }

#endif

#if E_DATA & SCALAR_DATA_IN_ELEMENT_NODES

void sn_set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_SET_VALUE(pelem,r,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_SET_VALUE(pelem,r,ze)
}

void sn_add_value_e(tGrid,r,ze)  /* ze_i := ze_i + r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_ADD_VALUE(pelem,r,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_ADD_VALUE(pelem,r,ze)
}

DOUBLE sn_sum_e(tGrid,ze)
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      sum += E_SUM(pelem,ze);

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += E_SUM(pelem,ze);
   return(sum);
}

DOUBLE sn_sum_e_n(tGrid,ze,n)
GRID *tGrid;
INT ze, *n;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.;
   
   *n = 0;
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      sum += E_SUM(pelem,ze);
      *n += DIM2;
   }

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            sum += E_SUM(pelem,ze);
            *n += DIM2;
         }
   return(sum);
}

void sn_set_first(tGrid,r,ze)  /* z_1 := r */
GRID *tGrid;
FLOAT r;
INT ze;
{ 
   EDSN(FIRSTELEMENT(tGrid),ze,DIM) = r;
}

void sn_subtr_one(tGrid,ze)
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   FLOAT f;
   
   f = EDSN(FIRSTELEMENT(tGrid),ze,DIM);
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_SUBTR_C(pelem,f,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_SUBTR_C(pelem,f,ze)
}

DOUBLE sn_mean_value_e(tGrid,ze)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   FLOAT sum=0., vol=0., evol;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      evol = VOLUME(pelem);
      sum += AVER_ESN(pelem,ze)*evol;
      vol += evol;
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            evol = VOLUME(pelem);
            sum += AVER_ESN(pelem,ze)*evol;
            vol += evol;
         }
   return(sum/vol);
}

void sn_set_zero_mean_e(tGrid,ze)
GRID *tGrid;
INT ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE mean;
   
   mean = sn_mean_value_e(tGrid,ze);

   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ){
      EDSN(pelem,ze,0) -= mean;
      EDSN(pelem,ze,1) -= mean;
      EDSN(pelem,ze,2) -= mean;
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid)){
            EDSN(pelem,ze,0) -= mean;
            EDSN(pelem,ze,1) -= mean;
            EDSN(pelem,ze,2) -= mean;
         }
}

void sn_copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_COPY(pelem,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_COPY(pelem,xe,ze)
}

void sn_inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_INV(pelem,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_INV(pelem,xe,ze)
}

void sn_add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_ADD(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_ADD(pelem,xe,ye,ze)
}

void sn_subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid;
INT xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_NSUBTR(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_NSUBTR(pelem,xe,ye,ze)
}

void sn_add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid;
INT we,xe,ye,ze;
{
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_NASUBTR(pelem,we,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_NASUBTR(pelem,we,xe,ye,ze)
}

void sn_mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_MADD(pelem,r,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_MADD(pelem,r,xe,ye,ze)
}

void sn_mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: sn_mult_and_subtr_e not available.\n");  }

void sn_damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid;
FLOAT r;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pel;
   
   for (pel = FIRSTELEMENT(tGrid); pel!=NULL; pel=pel->succ)
      E_DAMP(pel,r,xe,ye,ze)

   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel!=NULL; pel=pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            E_DAMP(pel,r,xe,ye,ze)
}

void sn_mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid;
FLOAT r;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_MULT(pelem,r,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_MULT(pelem,r,xe,ze)
}

DOUBLE sn_dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum=0.0;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      sum += DOTS(EDSNP(pelem,xe),EDSNP(pelem,ye));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += DOTS(EDSNP(pelem,xe),EDSNP(pelem,ye)); 
   return(sum);
}

DOUBLE sn_dot_e_without_one(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid;
INT xe,ye;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   DOUBLE sum;
   
   pelem = FIRSTELEMENT(tGrid);
   sum = DOT(EDSNP(pelem,xe),EDSNP(pelem,ye));
   for (pelem = SUCC(FIRSTELEMENT(tGrid)); pelem!=NULL; pelem=pelem->succ)
      sum += DOTS(EDSNP(pelem,xe),EDSNP(pelem,ye));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            sum += DOTS(EDSNP(pelem,xe),EDSNP(pelem,ye)); 
   return(sum);
}

void sn_divide_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_DIVIDE(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_DIVIDE(pelem,xe,ye,ze)
}

void sn_multiply_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid;
INT xe,ye,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_MULTIPLY(pelem,xe,ye,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_MULTIPLY(pelem,xe,ye,ze)
}
   
void sn_divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_divide_e not available.\n");  }

void sn_multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_multiply_e not available.\n");  }
   
void sn_make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid;
INT xe,ze;
{ 
   GRID *theGrid;
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem!=NULL; pelem=pelem->succ)
      E_MSQRT(pelem,xe,ze)
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pelem = FIRSTELEMENT(theGrid); pelem!=NULL; pelem=pelem->succ)
         if (IS_LTOP_ELEMENT(pelem,tGrid))
            E_MSQRT(pelem,xe,ze)
}

void sn_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         E_ADD_MULT(pel,(-h[i][j]),k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               E_ADD_MULT(pel,(-h[i][j]),k+i,l)
}

void sn_mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid;
FLOAT g[GMN];
INT xe, j, k;
{   
   GRID *theGrid;
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i <= j; i++)
         E_ADD_MULT(pel,g[i],k+i,xe)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser) 
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_LTOP_ELEMENT(pel,tGrid))
            for (i = 0; i <= j; i++)
               E_ADD_MULT(pel,g[i],k+i,xe)
}

void save_sne_vector(tGrid,z,m,fp)
GRID *tGrid;
INT z, m;
FILE *fp;
{
   ELEMENT *pelem;
   INT k;

   for (pelem=FIRSTELEMENT(tGrid); pelem; pelem=SUCC(pelem)){
      k = m + DIM2*pelem->eflag;
      fprintf(fp,"%i %e\n",k,EDSN(pelem,z,0));
      fprintf(fp,"%i %e\n",k+1,EDSN(pelem,z,1));
      fprintf(fp,"%i %e\n",k+2,EDSN(pelem,z,2));
      if (DIM == 3)
         fprintf(fp,"%i %e\n",k+3,EDSN(pelem,z,3));
   }
}

#else  /*  if !(E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)  */

void sn_set_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: sn_set_value_e not available.\n");  }

void sn_add_value_e(tGrid,r,ze)  /* ze := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: sn_add_value_e not available.\n");  }

DOUBLE sn_sum_e(tGrid,ze)
GRID *tGrid; INT ze;
{  eprintf("Error: sn_sum_e not available.\n"); return(0.);  }
 
DOUBLE sn_sum_e_n(tGrid,ze,n)
GRID *tGrid; INT ze, *n;
{  eprintf("Error: sn_sum_e_n not available.\n"); return(0.);  }
 
void sn_set_first(tGrid,r,ze)  /* z_1 := r */
GRID *tGrid; FLOAT r; INT ze;
{  eprintf("Error: sn_set_first not available.\n");  }

void sn_subtr_one(tGrid,ze)
GRID *tGrid; INT ze;
{  eprintf("Error: sn_subtr_one not available.\n");  }
 
DOUBLE sn_mean_value_e(tGrid,ze)  /*  (\int_omega ze dx)/|omega|  */
GRID *tGrid; INT ze;
{  eprintf("Error: sn_mean_value_e not available.\n"); return(0.);  }
 
void sn_set_zero_mean_e(tGrid,ze)
GRID *tGrid; INT ze;
{  eprintf("Error: sn_set_zero_mean_e not available.\n");  }
 
void sn_copy_e(tGrid,xe,ze)  /* ze := xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: sn_copy_e not available.\n");  }

void sn_inv_e(tGrid,xe,ze)  /* ze := -xe */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: sn_inv_e not available.\n");  }

void sn_add_e(tGrid,xe,ye,ze)  /* ze := xe + ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_add_e not available.\n");  }

void sn_subtr_e(tGrid,xe,ye,ze)  /* ze := xe - ye */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_subtr_e not available.\n");  }

void sn_add_and_subtr_e(tGrid,we,xe,ye,ze)  /* ze := we + xe - ye */
GRID *tGrid; INT we,xe,ye,ze;
{  eprintf("Error: sn_add_and_subtr_e not available.\n");  }

void sn_mult_and_add_e(tGrid,r,xe,ye,ze)  /* ze := r*xe + ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: sn_mult_and_add_e not available.\n");  }

void sn_mult_and_subtr_e(tGrid,r,xe,ye,ze)  /* ze := r*xe - ye */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: sn_mult_and_subtr_e not available.\n");  }

void sn_damp_e(tGrid,r,xe,ye,ze)  /* ze := xe + r*(ye - xe) */
GRID *tGrid; FLOAT r; INT xe,ye,ze;
{  eprintf("Error: sn_damp_e not available.\n");  }

void sn_mult_e(tGrid,r,xe,ze)  /* ze := r*xe */
GRID *tGrid; FLOAT r; INT xe,ze;
{  eprintf("Error: sn_mult_e not available.\n");  }

DOUBLE sn_dot_e(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: sn_dot_e not available.\n"); return(0.);  }

DOUBLE sn_dot_e_without_one(tGrid,xe,ye)  /* := xe*ye */
GRID *tGrid; INT xe,ye;
{  eprintf("Error: sn_dot_e_without_one not available.\n"); return(0.);  }

void sn_divide_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_divide_1e not available.\n");  }

void sn_multiply_1e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_multiply_1e not available.\n");  }

void sn_divide_e(tGrid,xe,ye,ze)  /* ze_i := xe_i/ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_divide_e not available.\n");  }

void sn_multiply_e(tGrid,xe,ye,ze)  /* ze_i := xe_i * ye_i */
GRID *tGrid; INT xe,ye,ze;
{  eprintf("Error: sn_multiply_e not available.\n");  }
   
void sn_make_sqrt_e(tGrid,xe,ze)  /* ze_i := sqrt(xe_i) */
GRID *tGrid; INT xe,ze;
{  eprintf("Error: sn_make_sqrt_e not available.\n");  }

void sn_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k;
{  eprintf("Error: sn_mult_and_add1_e_for_GMRES not available.\n");  }

void sn_mult_and_add2_e_for_GMRES(tGrid,g,xe,j,k)
GRID *tGrid; FLOAT g[GMN]; INT xe, j, k;
{  eprintf("Error: sn_mult_and_add2_e_for_GMRES not available.\n");  }

void save_sne_vector(tGrid,z,m,fp)
GRID *tGrid; INT z, m; FILE *fp;
{  eprintf("Error: save_sne_vector not available.\n");  }

#endif  /*  !(E_DATA & SCALAR_DATA_IN_ELEMENT_NODES)  */

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void vs_set_value_nf(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   GRID *theGrid;
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      SET_VALUE(theNode,r,z)
   fstop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      nstop = STOP_NODE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != nstop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            SET_VALUE(theNode,r,z)
      fstop = STOP_FACE(theGrid,t);
      for (theFace=FIRST_FACE(theGrid,t); theFace != fstop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r;
   }	
}

void vs_copy_nf(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   fstop = STOP_FACE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      COPY(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      nstop = STOP_NODE(theGrid,t);
      fstop = STOP_FACE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != nstop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            COPY(theNode,x,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != fstop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x);
   }	
}

void vs_inv_nf(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      INV(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = -FD(theFace,x);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            INV(theNode,x,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = -FD(theFace,x);
   }	
}

void vs_add_nf(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      ADD(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            ADD(theNode,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) + FD(theFace,y);
   }	
}

void vs_subtr_nf(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   fstop = STOP_FACE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      NSUBTR(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      nstop = STOP_NODE(theGrid,t);
      fstop = STOP_FACE(theGrid,t);
      for (theNode=FIRST_NODE(theGrid,t); theNode != nstop; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NSUBTR(theNode,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != fstop; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) - FD(theFace,y);
   }	
}

void vs_add_and_subtr_nf(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid;
INT w,x,y,z,t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      NASUBTR(theNode,w,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,w) + FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            NASUBTR(theNode,w,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,w) + FD(theFace,x) - FD(theFace,y);
   }	
}

void vs_mult_and_add_nf(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MADD(theNode,r,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = r*FD(theFace,x) + FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MADD(theNode,r,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r*FD(theFace,x) + FD(theFace,y);
   }	
}

void vs_mult_and_subtr_nf(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MSUBTR(theNode,r,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = r*FD(theFace,x) - FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MSUBTR(theNode,r,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = r*FD(theFace,x) - FD(theFace,y);
   }	
}

void vs_damp_nf(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DAMP(theNode,r,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + r*(FD(theFace,y) - FD(theFace,x));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            DAMP(theNode,r,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x) + r*(FD(theFace,y) - FD(theFace,x));
   }	
}

void vs_mult_nf(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULT(theNode,r,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*r;
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MULT(theNode,r,x,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)*r;
   }	
}

DOUBLE vs_dot_nf(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      sum += DOT(NDD(theNode,x),NDD(theNode,y));
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      sum += FD(theFace,x)*FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            sum += DOT(NDD(theNode,x),NDD(theNode,y));
      for (theFace = FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            sum += FD(theFace,x)*FD(theFace,y);
   }	
   return(sum);
}

void vs_divide_nf(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DIVIDE(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)/FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            DIVIDE(theNode,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)/FD(theFace,y);
   }
}

void vs_multiply_nf(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULTIPLY(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*FD(theFace,y);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MULTIPLY(theNode,x,y,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FD(theFace,x)*FD(theFace,y);
   }
}

void vs_make_sqrt_nf(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MSQRT(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = sqrt(FD(theFace,x));
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser){
      for (theNode=FIRST_NODE(theGrid,t); theNode != NULL; theNode=SUCC(theNode))
         if (IS_LTOP_NODE(theNode,tGrid))
            MSQRT(theNode,x,z)
      for (theFace=FIRST_FACE(theGrid,t); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = sqrt(FD(theFace,x));
   }
}

void vs_mult_and_add1_for_GMRES_nf(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   GRID *theGrid;
   NODE *theNode;
   FACE *pface;
   INT i;
   
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,(-h[i][j]),k+i,l)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(theNode,(-h[i][j]),k+i,l)
   for (pface = FIRST_FACE(tGrid,t); pface!=NULL; pface=pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,l) -=  h[i][j]*FD(pface,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               FD(pface,l) -=  h[i][j]*FD(pface,k+i);
}

void vs_mult_and_add2_for_GMRES_nf(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *pface;
   INT i;

   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,g[i],k+i,x)
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theNode = FIRST_NODE(theGrid,t); theNode != NULL; 
                                                      theNode = theNode->succ)
         if (IS_LTOP_NODE(theNode,tGrid))
            for (i = 0; i <= j; i++)
               ADD_MULT(theNode,g[i],k+i,x)
   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,x) += g[i]*FD(pface,k+i);
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (pface = FIRST_FACE(theGrid,t); pface != NULL; pface = pface->succ)  
         if (IS_LTOP_FACE(pface,tGrid))
            for (i = 0; i <= j; i++)
               FD(pface,x) += g[i]*FD(pface,k+i);
}

void print_vector(mg,tGrid,v)
MULTIGRID *mg;
GRID *tGrid;
INT v;
{
   GRID *theGrid;
   NODE *theNode;
   FACE *theFace;
   INT i,k;
  
   printf("\n");
   for (k=0;k<DIM;k++){
      i=1;
      printf("%i. component:\n",k+1);
      for (theGrid = FIRSTGRID(mg); theGrid->level <= tGrid->level; 
                                                       theGrid = theGrid->finer)
         for (theNode = FIRSTNODE(theGrid); theNode != NULL; 
                                                          theNode=SUCC(theNode))
            if (IS_LTOP_NODE(theNode,tGrid)) 
               printf("%3i:   %11.4e\n",i++,ND(theNode,v,k));
   }
   i=1;
   for (theGrid = FIRSTGRID(mg); theGrid->level <= tGrid->level; 
                                                       theGrid = theGrid->finer)
      for (theFace = FIRSTFACE(tGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)) 
            printf("%3i:   %11.4e\n",i++,FD(theFace,v));
}

#else  /* if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

void vs_set_value_nf(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z, t;
{  eprintf("Error: vs_set_value_nf not available.\n");  }

void vs_copy_nf(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vs_copy_nf not available.\n");  }

void vs_inv_nf(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vs_inv_nf not available.\n");  }

void vs_add_nf(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_add_nf not available.\n");  }

void vs_subtr_nf(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_subtr_nf not available.\n");  }

void vs_add_and_subtr_nf(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: vs_add_and_subtr_nf not available.\n");  }

void vs_mult_and_add_nf(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_mult_and_add_nf not available.\n");  }

void vs_mult_and_subtr_nf(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_mult_and_subtr_nf not available.\n");  }

void vs_damp_nf(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_damp_nf not available.\n");  }

void vs_mult_nf(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: vs_mult_nf not available.\n");  }

DOUBLE vs_dot_nf(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: vs_dot_nf not available.\n"); return(0.);  }

void vs_divide_nf(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: vs_divide_nf not available.\n");  }

void vs_multiply_nf(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x, y, z, t;
{  eprintf("Error: vs_multiply_nf not available.\n");  }

void vs_make_sqrt_nf(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x, z, t;
{  eprintf("Error: vs_make_sqrt_nf not available.\n");  }

void vs_mult_and_add1_for_GMRES_nf(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: vs_mult_and_add1_for_GMRES_nf not available.\n");  }

void vs_mult_and_add2_for_GMRES_nf(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: vs_mult_and_add2_for_GMRES_nf not available.\n");  }

void print_vector(mg,tGrid,v)
MULTIGRID *mg; GRID *tGrid; INT v;
{  eprintf("Error: print_vector not available.\n");  }

#endif  /* if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA))  */

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DATA_STR & LG_DATA)

void vs_set_value_nflg(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      SET_VALUE(theNode,r,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_SET_VALUE(theNode,r,z)
   fstop = STOP_FACE(tGrid,t);
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = r;
}

void vs_copy_nflg(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   fstop = STOP_FACE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      COPY(theNode,x,z)
   if (nstop == NULL)
      for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
         if (theNode->lgd)
            LG_COPY(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x);
}

void vs_inv_nflg(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      INV(theNode,x,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_INV(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = -FD(theFace,x);
}

void vs_add_nflg(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      ADD(theNode,x,y,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_ADD(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + FD(theFace,y);
}

void vs_subtr_nflg(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{ 
   NODE *theNode, *nstop;
   FACE *theFace, *fstop;
   	
   nstop = STOP_NODE(tGrid,t);
   fstop = STOP_FACE(tGrid,t);
   for (theNode=FIRST_NODE(tGrid,t); theNode != nstop; theNode=SUCC(theNode))
      NSUBTR(theNode,x,y,z)
   if (nstop == NULL)
      for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
         if (theNode->lgd)
            LG_SUBTR(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != fstop; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) - FD(theFace,y);
}

void vs_add_and_subtr_nflg(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: vs_add_and_subtr_nflg not available.\n");  }

void vs_mult_and_add_nflg(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MADD(theNode,r,x,y,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_MADD(theNode,r,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = r*FD(theFace,x) + FD(theFace,y);
}

void vs_mult_and_subtr_nflg(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_mult_and_subtr_nflg not available.\n");  }

void vs_damp_nflg(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DAMP(theNode,r,x,y,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_DAMP(theNode,r,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x) + r*(FD(theFace,y) - FD(theFace,x));
}

void vs_mult_nflg(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULT(theNode,r,x,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_MULT(theNode,r,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*r;
}

DOUBLE vs_dot_nflg(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   NODE *theNode;
   FACE *theFace;
   DOUBLE sum=0.0;
	
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      sum += DOT(NDD(theNode,x),NDD(theNode,y));
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         sum += LG_DOT(theNode,x,y)
   for (theFace = FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      sum += FD(theFace,x)*FD(theFace,y);
   return(sum);
}

void vs_divide_nflg(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x,y,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      DIVIDE(theNode,x,y,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_DIVIDE(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)/FD(theFace,y);
}

void vs_multiply_nflg(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x,y,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MULTIPLY(theNode,x,y,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_MULTIPLY(theNode,x,y,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FD(theFace,x)*FD(theFace,y);
}

void vs_make_sqrt_nflg(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x,z,t;
{ 
   NODE *theNode;
   FACE *theFace;
   	
   for (theNode=FIRST_NODE(tGrid,t); theNode != NULL; theNode=SUCC(theNode))
      MSQRT(theNode,x,z)
   for (theNode=FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); theNode=SUCC(theNode))
      if (theNode->lgd)
         LG_MSQRT(theNode,x,z)
   for (theFace=FIRST_FACE(tGrid,t); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = sqrt(FD(theFace,x));
}

void vs_mult_and_add1_for_GMRES_nflg(tGrid,h,l,j,k,t)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t;
{   
   NODE *theNode;
   FACE *pface;
   INT i;
   
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,(-h[i][j]),k+i,l)
   for (theNode = FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); 
                                                    theNode = theNode->succ)
      if (theNode->lgd)
         for (i = 0; i <= j; i++)
            LG_MADD(theNode,(-h[i][j]),k+i,l,l)
   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,l) -=  h[i][j]*FD(pface,k+i);
}

void vs_mult_and_add2_for_GMRES_nflg(tGrid,g,x,j,k,t)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t;
{   
   NODE *theNode;
   FACE *pface;
   INT i, t;
   
   for (theNode = FIRST_NODE(tGrid,t); theNode != NULL; theNode = theNode->succ)
      for (i = 0; i <= j; i++)
         ADD_MULT(theNode,g[i],k+i,x)
   for (theNode = FIRSTN(tGrid); theNode != FIRSTNODE(tGrid); 
                                                     theNode = theNode->succ)
      if (theNode->lgd)
         for (i = 0; i <= j; i++)
            LG_MADD(theNode,g[i],k+i,x,x)
   for (pface = FIRST_FACE(tGrid,t); pface != NULL; pface = pface->succ)
      for (i = 0; i <= j; i++)
         FD(pface,x) += g[i]*FD(pface,k+i);
}

#else  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                (DATA_STR & LG_DATA))  */

void vs_set_value_nflg(tGrid,r,z,t)  /* z := r */
GRID *tGrid; FLOAT r; INT z,t;
{  eprintf("Error: vs_set_value_nflg not available.\n");  }

void vs_copy_nflg(tGrid,x,z,t)  /* z := x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vs_copy_nflg not available.\n");  }

void vs_inv_nflg(tGrid,x,z,t)  /* z := -x */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vs_inv_nflg not available.\n");  }

void vs_add_nflg(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_add_nflg not available.\n");  }

void vs_subtr_nflg(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_subtr_nflg not available.\n");  }

void vs_add_and_subtr_nflg(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid; INT w,x,y,z,t;
{  eprintf("Error: vs_add_and_subtr_nflg not available.\n");  }

void vs_mult_and_add_nflg(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_mult_and_add_nflg not available.\n");  }

void vs_mult_and_subtr_nflg(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_mult_and_subtr_nflg not available.\n");  }

void vs_damp_nflg(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid; FLOAT r; INT x,y,z,t;
{  eprintf("Error: vs_damp_nflg not available.\n");  }

void vs_mult_nflg(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid; FLOAT r; INT x,z,t;
{  eprintf("Error: vs_mult_nflg not available.\n");  }

DOUBLE vs_dot_nflg(tGrid,x,y,t)  /* := x.y */
GRID *tGrid; INT x,y,t;
{  eprintf("Error: vs_dot_nflg not available.\n"); return(0.);  }

void vs_divide_nflg(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_divide_nflg not available.\n");  }

void vs_multiply_nflg(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid; INT x,y,z,t;
{  eprintf("Error: vs_multiply_nflg not available.\n");  }

void vs_make_sqrt_nflg(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid; INT x,z,t;
{  eprintf("Error: vs_make_sqrt_nflg not available.\n");  }

void vs_mult_and_add1_for_GMRES_nflg(tGrid,h,l,j,k,t)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k, t;
{  eprintf("Error: vs_mult_and_add1_for_GMRES_nflg not available.\n");  }

void vs_mult_and_add2_for_GMRES_nflg(tGrid,g,x,j,k,t)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k, t;
{  eprintf("Error: vs_mult_and_add2_for_GMRES_nflg not available.\n");  }

#endif  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                 (DATA_STR & LG_DATA))  */

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (DATA_STR & LG_DATA) && (DIM == 3)

void print_vector_lg(mg,tGrid,u)
MULTIGRID *mg;
GRID *tGrid;
INT u;
{    
   NODE *pnode;
   FACE *pface;
   INT i=0;
   
   printf("Nodes on GLG:\n");  
   for (pnode = FIRSTN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (pnode->lgd)
         printf("%4i:  %9.2e   %9.2e\n",++i,NDLG(pnode,u,0),NDLG(pnode,u,1));
   i = 0;
   printf("\nInner nodes:\n");
   for (pnode = FIRSTNODE(tGrid); pnode != NULL ; pnode = pnode->succ)
      printf("%4i:  %9.2e   %9.2e   %9.2e\n",
                              ++i,ND(pnode,u,0),ND(pnode,u,1),ND(pnode,u,2));
   i = 0;
   printf("\nInner faces:\n");
   for (pface = FIRSTFACE(tGrid); pface != NULL ; pface = pface->succ)
      printf("%4i:  %9.2e\n",++i,FD(pface,u));
}

void set_LG_values(tGrid,u)
GRID *tGrid;
INT u;
{
   NODE *pnode;
   
   for (pnode = FIRSTN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (pnode->lgd){
         ND(pnode,u,0) = pnode->lgd->t1[0]*NDLG(pnode,u,0) + 
                         pnode->lgd->t2[0]*NDLG(pnode,u,1);
         ND(pnode,u,1) = pnode->lgd->t1[1]*NDLG(pnode,u,0) + 
                         pnode->lgd->t2[1]*NDLG(pnode,u,1);
         ND(pnode,u,2) = pnode->lgd->t1[2]*NDLG(pnode,u,0) + 
                         pnode->lgd->t2[2]*NDLG(pnode,u,1);
      }
}

#else  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                (DATA_STR & LG_DATA) && (DIM == 3))  */

void print_vector_lg(mg,tGrid,u)
MULTIGRID *mg; GRID *tGrid; INT u;
{  eprintf("Error: print_vector_lg not available.\n");  }

void set_LG_values(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: set_LG_values not available.\n");  }

#endif  /*  if !((N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && 
                 (DATA_STR & LG_DATA) && (DIM == 3))  */

#if (N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA)

void sset_value_BD(tGrid,r,z,t)	//z := r
GRID *tGrid;
FLOAT r;
INT z, t;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = r;
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = r;

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = r;
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = r;
	}
}

void scopy_BD(tGrid,x,z)	//z = x
GRID *tGrid;
INT x, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = NDSBD(pnode,x);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = FDBD(pface,x);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = NDSBD(pnode,x);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = FDBD(pface,x);
	}
}

void sinv_BD(tGrid,x,z)		//z = -x
GRID *tGrid;
INT x, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = -NDSBD(pnode,x);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = -FDBD(pface,x);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = -NDSBD(pnode,x);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = -FDBD(pface,x);
	}
}

void ssubtr_BD(tGrid,x,y,z)	//z = x - y
GRID *tGrid;
INT x, y, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = NDSBD(pnode,x) - NDSBD(pnode,y);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = FDBD(pface,x) - FDBD(pface,y);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = NDSBD(pnode,x) - NDSBD(pnode,y);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = FDBD(pface,x) - FDBD(pface,y);
	}
}

void smult_BD(tGrid,r,x,z)	//z = r*x
GRID *tGrid;
FLOAT r;
INT x, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = NDSBD(pnode,x)*r;
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = FDBD(pface,x)*r;

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = NDSBD(pnode,x)*r;
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = FDBD(pface,x)*r;
	}
}

void smult_and_add_BD(tGrid,r,x,y,z)	//z = r*x + y
GRID *tGrid;
FLOAT r;
INT x, y, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = r*NDSBD(pnode,x) + NDSBD(pnode,y);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = r*FDBD(pface,x) + FDBD(pface,y);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = r*NDSBD(pnode,x) + NDSBD(pnode,y);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = r*FDBD(pface,x) + FDBD(pface,y);
	}
}

void smult_and_subtr_BD(tGrid,r,x,y,z)	//z = r*x - y
GRID *tGrid;
FLOAT r;
INT x, y, z;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			NDSBD(pnode,z) = r*NDSBD(pnode,x) - NDSBD(pnode,y);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			FDBD(pface,z) = r*FDBD(pface,x) - FDBD(pface,y);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				NDSBD(pnode,z) = r*NDSBD(pnode,x) - NDSBD(pnode,y);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				FDBD(pface,z) = r*FDBD(pface,x) - FDBD(pface,y);
	}
}

DOUBLE sdot_BD(tGrid,x,y)	// = x.y
GRID *tGrid;
INT x, y;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;
	DOUBLE sum = 0.0;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			sum += NDSBD(pnode,x)*NDSBD(pnode,y);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			sum += FDBD(pface,x)*FDBD(pface,y);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				sum += NDSBD(pnode,x)*NDSBD(pnode,y);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				sum += FDBD(pface,x)*FDBD(pface,y);
	}
	return (sum);
}

void smult_and_add1_for_GMRES_BD(tGrid,h,l,j,k)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;
	INT i;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			for (i = 0; i <= j; i++)
				NDSBD(pnode,l) -= h[i][j]*NDSBD(pnode,k+i);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			for (i = 0; i <= j; i++)
				FDBD(pface,l) -= h[i][j]*FDBD(pface,k+i);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				for (i = 0; i <= j; i++)
					NDSBD(pnode,l) -= h[i][j]*NDSBD(pnode,k+i);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				for (i = 0; i <= j; i++)
					FDBD(pface,l) -= h[i][j]*FDBD(pface,k+i);
	}
}

void smult_and_add2_for_GMRES_BD(tGrid,g,x,j,k)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k;
{
	GRID *pg;
	NODE *pnode;
	FACE *pface;
	INT i;

	for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
		if (IS_ZNN(pnode))
			for (i = 0; i <= j; i++)
				NDSBD(pnode,x) += g[i]*NDSBD(pnode,k+i);
	for (pface = FIRSTFACE(tGrid); pface; pface = pface->succ)
		if (IS_ZNF(pface))
			for (i = 0; i <= j; i++)
				FDBD(pface,x) += g[i]*FDBD(pface,k+i);

	for (pg = tGrid->coarser; pg; pg = pg->coarser) {
		for (pnode = FIRSTNODE(pg); pnode; pnode = pnode->succ)
			if (IS_ZNN(pnode) && IS_LTOP_NODE(pnode,tGrid))
				for (i = 0; i <= j; i++)
					NDSBD(pnode,x) += g[i]*NDSBD(pnode,k+i);
		for (pface = FIRSTFACE(pg); pface; pface = pface->succ)
			if (IS_ZNF(pface) && IS_LTOP_FACE(pface,tGrid))
				for (i = 0; i <= j; i++)
					FDBD(pface,x) += g[i]*FDBD(pface,k+i);
	}
}

#else	/* if !((N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA)) */

void sset_value_BD(tGrid,r,z,t)	//z := r
GRID *tGrid;
FLOAT r;
INT z, t;
{ eprintf("Error: sset_value_BD not available.\n"); }

void scopy_BD(tGrid,x,z)	//z = x
GRID *tGrid; INT x, z;
{ eprintf("Error: scopy_BD not available.\n"); }

void sinv_BD(tGrid,x,z)		//z = -x
GRID *tGrid; INT x, z;
{ eprintf("Error: sinv_BD not available.\n"); }

void ssubtr_BD(tGrid,x,y,z)	//z = x - y
GRID *tGrid; INT x, y, z;
{ eprintf("Error: ssubtr_BD not available.\n"); }

void smult_BD(tGrid,r,x,z)	//z = r*x
GRID *tGrid; FLOAT r; INT x, z;
{ eprintf("Error: smult_BD not available.\n"); }

void smult_and_add_BD(tGrid,r,x,y,z)	//z = r*x + y
GRID *tGrid; FLOAT r; INT x, y, z;
{ eprintf("Error: smult_and_add_BD not available.\n"); }

void smult_and_subtr_BD(tGrid,r,x,y,z)	//z = r*x - y
GRID *tGrid; FLOAT r; INT x, y, z;
{ eprintf("Error: smult_and_subtr_BD not available.\n"); }

DOUBLE sdot_BD(tGrid,x,y)	// = x.y
GRID *tGrid; INT x, y;
{ eprintf("Error: sdot_BD not available.\n"); return(0.); }

void smult_and_add1_for_GMRES_BD(tGrid,h,l,j,k)
GRID *tGrid; FLOAT h[GMN][GMN]; INT l, j, k;
{ eprintf("Error: smult_and_add1_for_GMRES_BD not available.\n"); }

void smult_and_add2_for_GMRES_BD(tGrid,g,x,j,k)
GRID *tGrid; FLOAT g[GMN]; INT x, j, k;
{ eprintf("Error: smult_and_add2_for_GMRES_BD not available.\n"); }

#endif	/* !((N_DATA & NODE_BD_DATA) && (F_DATA & FACE_BD_DATA)) */

void set_value_ssn(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = r;
}

void copy_ssn(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x);
}

void inv_ssn(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = -SNDS(psnode,x);
}

void add_ssn(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x) + SNDS(psnode,y);
}

void subtr_ssn(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x) - SNDS(psnode,y);
}

void add_and_subtr_ssn(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid;
INT w,x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,w) + SNDS(psnode,x) - SNDS(psnode,y);
}

void mult_and_add_ssn(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = r*SNDS(psnode,x) + SNDS(psnode,y);
}

void mult_and_subtr_ssn(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = r*SNDS(psnode,x) - SNDS(psnode,y);
}

void damp_ssn(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x) + r*(SNDS(psnode,y) - SNDS(psnode,x));
}

void mult_ssn(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x)*r;
}

DOUBLE dot_ssn(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   SNODE *psnode;
   DOUBLE sum=0.0;
	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      sum += SNDS(psnode,x)*SNDS(psnode,y);
   return(sum);
}

void divide_ssn(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x)/SNDS(psnode,y);
}

void multiply_ssn(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = SNDS(psnode,x)*SNDS(psnode,y);
}

void make_sqrt_ssn(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   SNODE *psnode;
   	
   for (psnode = FIRSTSN(tGrid); psnode; psnode = SUCC(psnode))
      SNDS(psnode,z) = sqrt(SNDS(psnode,x));
}

void set_value_ssf(tGrid,r,z,t)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = r;
}

void copy_ssf(tGrid,x,z,t)  /* z := x */
GRID *tGrid;
INT x,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x);
}

void inv_ssf(tGrid,x,z,t)  /* z := -x */
GRID *tGrid;
INT x,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = -SFDS(psface,x);
}

void add_ssf(tGrid,x,y,z,t)  /* z := x + y */
GRID *tGrid;
INT x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x) + SFDS(psface,y);
}

void subtr_ssf(tGrid,x,y,z,t)  /* z := x - y */
GRID *tGrid;
INT x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x) - SFDS(psface,y);
}

void add_and_subtr_ssf(tGrid,w,x,y,z,t)  /* z := w + x - y */
GRID *tGrid;
INT w,x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,w) + SFDS(psface,x) - SFDS(psface,y);
}

void mult_and_add_ssf(tGrid,r,x,y,z,t)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = r*SFDS(psface,x) + SFDS(psface,y);
}

void mult_and_subtr_ssf(tGrid,r,x,y,z,t)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = r*SFDS(psface,x) - SFDS(psface,y);
}

void damp_ssf(tGrid,r,x,y,z,t)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x,y,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x) + r*(SFDS(psface,y) - SFDS(psface,x));
}

void mult_ssf(tGrid,r,x,z,t)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x,z,t;
{ 
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x)*r;
}

DOUBLE dot_ssf(tGrid,x,y,t)  /* := x.y */
GRID *tGrid;
INT x,y,t;
{ 
   SFACE *psface;
   DOUBLE sum=0.0;
	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      sum += SFDS(psface,x)*SFDS(psface,y);
   return(sum);
}

void divide_ssf(tGrid,x,y,z,t)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x)/SFDS(psface,y);
}

void multiply_ssf(tGrid,x,y,z,t)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t;
{
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = SFDS(psface,x)*SFDS(psface,y);
}

void make_sqrt_ssf(tGrid,x,z,t)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t;
{
   SFACE *psface;
   	
   for (psface = FIRSTSF(tGrid); psface; psface = SUCC(psface))
      SFDS(psface,z) = sqrt(SFDS(psface,x));
}

#if (N_DATA & SCALAR_NODE_DATA) && (DIM == 2)

DOUBLE p1c_normal_derivative(n0,n1,n2,u)  /*  n0, n1, n2 are nodes of a       */
NODE *n0, *n1, *n2;        /*  triangle, u denotes a scalar linear function.  */
INT u;                     /*  The result is the derivative of u w.r.t. the   */
{                          /*  normal vector to [n0,n1] pointing to n2.       */
   DOUBLE n[DIM], t[DIM], s[DIM], d, u3;

   SUBTR(n1->myvertex->x,n0->myvertex->x,t);
   SUBTR(n2->myvertex->x,n0->myvertex->x,s);
   d = DOT(t,t);
   u3 = NDS(n0,u) + (NDS(n1,u)-NDS(n0,u))*DOT(t,s)/d;
   ORT_VECT(n,t)
   d = fabs(DOT(n,s))/sqrt(d);
   return((NDS(n2,u)-u3)/d);
}

#else

DOUBLE p1c_normal_derivative(n0,n1,n2,u)
NODE *n0, *n1, *n2; INT u;
{  eprintf("Error: p1c_normal_derivative not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (DATA_S & N_LINK_TO_NODES)

DOUBLE p1c_other_normal_derivative(n0,n1,n2,u)
NODE *n0, *n1, *n2;
INT u;
{
   LINK *pl, *pli;
   INT i=0;
   DOUBLE res;

   for (pl = n0->tstart; pl; pl = pl->next)
      for (pli = n1->tstart; pli; pli = pli->next)
         if ((pl->nbnode == pli->nbnode) && (pl->nbnode != n2)){
            res = p1c_normal_derivative(n0,n1,pl->nbnode,u);
            i++;
         }
   if (i == 0)
      eprintf("Error in p1c_other_normal_derivative (i==0).\n");
   else if (i > 1)
      eprintf("Error in p1c_other_normal_derivative (i>1).\n");
   return(res);
}

#else

DOUBLE p1c_other_normal_derivative(n0,n1,n2,u)
NODE *n0, *n1, *n2; INT u;
{  eprintf("Error: p1c_other_normal_derivative not available.\n");  }

#endif

DOUBLE p1c_jump_of_normal_derivative(n0,n1,n2,f2,u)
NODE *n0, *n1, *n2;
FACE *f2;
INT u;
{
   if (IS_BF(f2)){
      if (IS_FN(n0) || IS_FN(n1))
         return(0.);
      else
         return(p1c_normal_derivative(n0,n1,n2,u));
   }
   else
      return(p1c_normal_derivative(n0,n1,n2,u) +
             p1c_other_normal_derivative(n0,n1,n2,u));
}

#if (F_DATA & SCALAR_FACE_DATA) && (F_DATA & VECTOR_FACE_DATA)

void vector_to_scalar_f(tGrid,x,z,i)  /* z := x_i */
GRID *tGrid;
INT x,z,i;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FDV(theFace,x,i);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FDV(theFace,x,i);
}

#else  /*  if !((F_DATA & SCALAR_FACE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

void vector_to_scalar_f(tGrid,x,z,i)  /* z := x_i */
GRID *tGrid; INT x,z,i;
{  eprintf("Error: vector_to_scalar_f not available.\n");  }

#endif  /*  !((F_DATA & SCALAR_FACE_DATA) && (F_DATA & VECTOR_FACE_DATA))  */

#if (F_DATA & SCALAR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void dvector_to_scalar_f(tGrid,x,z,i,j)  /* z := x_(i,j) */
GRID *tGrid;
INT x,z,i,j;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace))
      FD(theFace,z) = FDDV(theFace,x,i,j);
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            FD(theFace,z) = FDDV(theFace,x,i,j);
}

#else  /*  if !((F_DATA & SCALAR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void dvector_to_scalar_f(tGrid,x,z,i,j)  /* z := x_(i,j) */
GRID *tGrid; INT x,z,i,j;
{  eprintf("Error: dvector_to_scalar_f not available.\n");  }

#endif  /*  !((F_DATA & SCALAR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

#if (F_DATA & VECTOR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA)

void dvector_to_vector_f(tGrid,x,z,i)
GRID *tGrid;
INT x,z,i;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace))
      SET1(FDVP(theFace,z),FDDVP(theFace,x,i))
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid))
            SET1(FDVP(theFace,z),FDDVP(theFace,x,i))
}

void dvector_to_dscalar_f(tGrid,x,z,i)
GRID *tGrid;
INT x,z,i;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      FDV(theFace,z,0) = FDDV(theFace,x,0,i);
      FDV(theFace,z,1) = FDDV(theFace,x,1,i);
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            FDV(theFace,z,0) = FDDV(theFace,x,0,i);
            FDV(theFace,z,1) = FDDV(theFace,x,1,i);
         }
}

void dscalar_to_dvector_f(tGrid,x,z,i)
GRID *tGrid;
INT x,z,i;
{ 
   GRID *theGrid;
   FACE *theFace;
   	
   for (theFace=FIRSTF(tGrid); theFace != NULL; theFace=SUCC(theFace)){
      FDDV(theFace,z,0,i) = FDV(theFace,x,0);
      FDDV(theFace,z,1,i) = FDV(theFace,x,1);
   }
   
   for (theGrid = tGrid->coarser; theGrid; theGrid = theGrid->coarser)
      for (theFace=FIRSTF(theGrid); theFace != NULL; theFace=SUCC(theFace))
         if (IS_LTOP_FACE(theFace,tGrid)){
            FDDV(theFace,z,0,i) = FDV(theFace,x,0);
            FDDV(theFace,z,1,i) = FDV(theFace,x,1);
         }
}

#else  /* if !((F_DATA & VECTOR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

void dvector_to_vector_f(tGrid,x,z,i)
GRID *tGrid; INT x,z,i;
{  eprintf("Error: dvector_to_vector_f not available.\n");  }

void dvector_to_dscalar_f(tGrid,x,z,i)
GRID *tGrid; INT x,z,i;
{  eprintf("Error: dvector_to_dscalar_f not available.\n");  }

void dscalar_to_dvector_f(tGrid,x,z,i)
GRID *tGrid; INT x,z,i;
{  eprintf("Error: dscalar_to_dvector_f not available.\n");  }

#endif  /* !((F_DATA & VECTOR_FACE_DATA) && (F_DATA & DVECTOR_FACE_DATA))  */

#if N_DATA & SCALAR_NODE_DATA

DOUBLE ssn_node_value(pnode,u)
NODE *pnode;
INT u;
{
   return(NDS(pnode,u));
}

void add_to_ssn_node_value(pnode,u,r)
NODE *pnode;
DOUBLE r;
INT u;
{
   NDS(pnode,u) += r;
}

#else

DOUBLE ssn_node_value(pnode,u)
NODE *pnode; INT u;
{  eprintf("Error: ssn_node_value not available.\n");  }

void add_to_ssn_node_value(pnode,u,r)
NODE *pnode; DOUBLE r; INT u;
{  eprintf("Error: add_to_ssn_node_value not available.\n");  }

#endif

#if N_DATA & MVECTOR_NODE_DATA

DOUBLE gsn_node_value(pnode,u)
NODE *pnode;
INT u;
{
   return(NDMV(pnode,u,0));
}

void add_to_gsn_node_value(pnode,u,r)
NODE *pnode;
DOUBLE r;
INT u;
{
   NDMV(pnode,u,0) += r;
}

#else

DOUBLE gsn_node_value(pnode,u)
NODE *pnode; INT u;
{  eprintf("Error: gsn_node_value not available.\n");  }

void add_to_gsn_node_value(pnode,u,r)
NODE *pnode; DOUBLE r; INT u; 
{  eprintf("Error: add_to_gsn_node_value not available.\n");  }

#endif

DOUBLE sn_node_value(pnode,u,space)
NODE *pnode;
INT u, space;
{
      if (space < 500)
         return(ssn_node_value(pnode,u));
      else
         return(gsn_node_value(pnode,u));
}

void add_to_sn_node_value(pnode,u,r,space)
NODE *pnode;
DOUBLE r;
INT u, space;
{
      if (space < 500)
         add_to_ssn_node_value(pnode,u,r);
      else
         add_to_gsn_node_value(pnode,u,r);
}

#if F_DATA & SCALAR_FACE_DATA

DOUBLE ssf_face_value(pface,u)
FACE *pface;
INT u;
{
   return(FD(pface,u));
}

void add_to_ssf_face_value(pface,u,r)
FACE *pface;
DOUBLE r;
INT u;
{
   FD(pface,u) += r;
}

#else

DOUBLE ssf_face_value(pface,u)
FACE *pface; INT u;
{  eprintf("Error: ssf_face_value not available.\n");  }

void add_to_ssf_face_value(pface,u,r)
FACE *pface; DOUBLE r; INT u;
{  eprintf("Error: add_to_ssf_face_value not available.\n");  }

#endif

#if F_DATA & MVECTOR_FACE_DATA

DOUBLE gsf_face_value(pface,u)
FACE *pface;
INT u;
{
   return(FDMV(pface,u,0));
}

void add_to_gsf_face_value(pface,u,r)
FACE *pface;
DOUBLE r;
INT u;
{
   FDMV(pface,u,0) += r;
}

#else

DOUBLE gsf_face_value(pface,u)
FACE *pface; INT u;
{  eprintf("Error: gsf_face_value not available.\n");  }

void add_to_gsf_face_value(pface,u,r)
FACE *pface; DOUBLE r; INT u;
{  eprintf("Error: add_to_gsf_face_value not available.\n");  }

#endif

DOUBLE sf_face_value(pface,u,space)
FACE *pface;
INT u, space;
{
      if (space < 500)
         return(ssf_face_value(pface,u));
      else
         return(gsf_face_value(pface,u));
}

void add_to_sf_face_value(pface,u,r,space)
FACE *pface;
DOUBLE r;
INT u, space;
{
      if (space < 500)
         add_to_ssf_face_value(pface,u,r);
      else
         add_to_gsf_face_value(pface,u,r);
}

#if E_DATA & SCALAR_ELEMENT_DATA

DOUBLE sse_element_value(pel,u)
ELEMENT *pel;
INT u;
{
   return(ED(pel,u));
}

void add_to_sse_element_value(pel,u,r)
ELEMENT *pel;
DOUBLE r;
INT u;
{
   ED(pel,u) += r;
}

#else

DOUBLE sse_element_value(pel,u)
ELEMENT *pel; INT u;
{  eprintf("Error: sse_element_value not available.\n");  }

void add_to_sse_element_value(pel,u,r)
ELEMENT *pel; DOUBLE r; INT u;
{  eprintf("Error: add_to_sse_element_value not available.\n");  }

#endif

#if E_DATA & MVECTOR_ELEMENT_DATA

DOUBLE gse_element_value(pel,u)
ELEMENT *pel;
INT u;
{
   return(EDMV(pel,u,0));
}

void add_to_gse_element_value(pel,u,r)
ELEMENT *pel;
DOUBLE r;
INT u;
{
   EDMV(pel,u,0) += r;
}

#else

DOUBLE gse_element_value(pel,u)
ELEMENT *pel; INT u;
{  eprintf("Error: gse_element_value not available.\n");  }

void add_to_gse_element_value(pel,u,r)
ELEMENT *pel; DOUBLE r; INT u;
{  eprintf("Error: add_to_gse_element_value not available.\n");  }

#endif

DOUBLE se_element_value(pel,u,space)
ELEMENT *pel;
INT u, space;
{
      if (space < 500)
         return(sse_element_value(pel,u));
      else
         return(gse_element_value(pel,u));
}

void add_to_se_element_value(pel,u,r,space)
ELEMENT *pel;
DOUBLE r;
INT u, space;
{
      if (space < 500)
         add_to_sse_element_value(pel,u,r);
      else
         add_to_gse_element_value(pel,u,r);
}

INT eqstr(s,t)
char s[], t[];
{
   INT i=0;

   while (s[i] == t[i])
      if (s[i++] == '\0')
         return(1);
   return(0);
}

INT operation_available(t,type,name)
INT t, *type;
char name[];
{
   if ( ( (t & USE_IS_DIR) && !(*type == Q_SN || *type == Q_SNSE) ) ||
        ( (t & USE_IS_DIR) &&
                (*type == Q_SE || *type == Q_SNE || *type == Q_VE) ) ||
        ( (t & (STOP_IS_FIRST_INNER | ONLY_INNER)) && (*type == Q_SNE) ) ||
        ( (t & STOP_IS_FIRST_INNER) && (t & ONLY_INNER) ) ||
        ( (t & STOP_IS_FIRST_INNER) && !(eqstr(name,"set_value") || 
                   eqstr(name,"copy") || eqstr(name,"subtr")) ) ){
      eeprintf("Error: %s not available for t = %i and type = %i.\n",
              name,t,*type);
      return(0);
   }
   else if ( (t & STOP_IS_FIRST_INNER) && (*type == Q_SE || *type == Q_VE) ) 
      return(0);
   else{
      if (t & STOP_IS_FIRST_INNER)
         switch(*type){
            case Q_SNSE: *type = Q_SN;
                 break;
            case Q_SFSE: *type = Q_SF;
                 break;
            case Q_VNVE: *type = Q_VN;
                 break;
            case Q_VNVFVE: *type = Q_VNVF;
                 break;
            case Q_GENERAL: *type = Q_MVNMVF;
                 break;
            default:
                 break;
         }
      return(1);
   }
}

void set_value(tGrid,r,z,t,type)  /* z := r */
GRID *tGrid;
FLOAT r;
INT z, t, type;
{
   if(operation_available(t,&type,"set_value"))
   switch(type){
   case Q_SN: sset_value(tGrid,r,z,t);
              if (t & ADD_ZN_CONDITION)
                 sset_value_BD(tGrid,r,z,t);
        break;
   case Q_VN: v_set_value_n(tGrid,r,z,t);
        break;
   case Q_SF: set_value_f(tGrid,r,z,t);
        break;
   case Q_VF: vset_value_f(tGrid,r,z,t);
        break;
   case Q_DVF: dvset_value_f(tGrid,r,z,t);
        break;
   case Q_SE: set_value_e(tGrid,r,z);
        break;
   case Q_SNE: sn_set_value_e(tGrid,r,z);
        break;
   case Q_VE: vset_value_e(tGrid,r,z);
        break;
   case Q_SNSF: sset_value(tGrid,r,z,t);
                set_value_f(tGrid,r,z,t);
        break;
   case Q_SNSE: sset_value(tGrid,r,z,t);
                set_value_e(tGrid,r,z);
        break;
   case Q_SFSE: set_value_f(tGrid,r,z,t);
                set_value_e(tGrid,r,z);
        break;
   case Q_VNSF: vs_set_value_nf(tGrid,r,z,t);
        break;
   case Q_VNVF: v_set_value_n(tGrid,r,z,t);
                vset_value_f(tGrid,r,z,t);
        break;
   case Q_VNVE: v_set_value_n(tGrid,r,z,t);
                vset_value_e(tGrid,r,z);
        break;
   case Q_VNVFVE: v_set_value_n(tGrid,r,z,t);
                  vset_value_f(tGrid,r,z,t);
                  vset_value_e(tGrid,r,z);
        break;
   case Q_MVNMVF: mv_set_value_n(tGrid,r,z,t);
                  mv_set_value_f(tGrid,r,z,t);
        break;
   case Q_GENERAL: mv_set_value_n(tGrid,r,z,t);
                   mv_set_value_f(tGrid,r,z,t);
                   mv_set_value_e(tGrid,r,z);
        break;
   case Q_SSN: set_value_ssn(tGrid,r,z,t);
        break;
   case Q_SSF: set_value_ssf(tGrid,r,z,t);
        break;
   case Q_SNSSNSSF: sset_value(tGrid,r,z,t);
                    set_value_ssn(tGrid,r,z,t);
                    set_value_ssf(tGrid,r,z,t);
        break;
   case Q_VNSFLG: vs_set_value_nflg(tGrid,r,z,t);
        break;
   default:
        eprintf("Error: set_value not available.\n");
        break;
   }
}

DOUBLE max_abs_value(tGrid,z,t,type)  /* max := max(z_i) */
GRID *tGrid;
INT z, t, type;
{
   if(operation_available(t,&type,"max_abs_value"))
   switch(type){
   case Q_SN: if (t & ADD_ZN_CONDITION)
                 eprintf("Error: max_abs_value not available for ADD_ZN_CONDITION.\n");
              return(smax_abs_value(tGrid,z,t));
        break;
   case Q_SE: return(max_abs_value_e(tGrid,z));
        break;
   default:
        eprintf("Error: max_abs_value not available.\n");
        return(0.);
        break;
   }
}

void add_value(tGrid,r,z,t,type)  /* z_i := z_i + r */
GRID *tGrid;
FLOAT r;
INT z, t, type;
{
   if(operation_available(t,&type,"add_value"))
   switch(type){
   case Q_SN: sadd_value(tGrid,r,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: add_value not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SE: add_value_e(tGrid,r,z);
        break;
   case Q_SNE: sn_add_value_e(tGrid,r,z);
        break;
   case Q_SNSSNSSF: sadd_value(tGrid,r,z,t);
        break;
   default:
        eprintf("Error: add_value not available.\n");
        break;
   }
}

DOUBLE sum(tGrid,z,t,type)  /* sum of the components (inner product with 1) */
GRID *tGrid;
INT z, t, type;
{
   if(operation_available(t,&type,"sum"))
   switch(type){
   case Q_SN: return(ssum(tGrid,z,t));
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: sum not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SE: return(sum_e(tGrid,z));
        break;
   case Q_SNE: return(sn_sum_e(tGrid,z));
        break;
   case Q_SNSSNSSF: return(ssum(tGrid,z,t));
        break;
   default:
        eprintf("Error: sum not available.\n");
        return(0.);
        break;
   }
}

DOUBLE sum_n(tGrid,z,n,t,type)/* sum of the components (inner product with 1) */
GRID *tGrid;                  /* n is the number of values which are summed   */
INT z, *n, t, type;
{
   if(operation_available(t,&type,"sum_n"))
   switch(type){
   case Q_SN: return(ssum_n(tGrid,z,n,t));
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: sum_n not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SE: return(sum_e_n(tGrid,z,n));
        break;
   case Q_SNE: return(sn_sum_e_n(tGrid,z,n));
        break;
   case Q_SNSSNSSF: return(ssum_n(tGrid,z,n,t));
        break;
   default:
        eprintf("Error: sum_n not available.\n");
        return(0.);
        break;
   }
}

DOUBLE mean_value(tGrid,z,t,type)  /* integral mean of the resp. FE function */
GRID *tGrid;
INT z, t, type;
{
   switch(type){
   case Q_SN: return(smean_value(tGrid,z));
        break;
   case Q_SE: return(mean_value_e(tGrid,z));
        break;
   case Q_SNE: return(sn_mean_value_e(tGrid,z));
        break;
   default:
        eprintf("Error: mean_value not available.\n");
        return(0.);
        break;
   }
}

void set_zero_mean(tGrid,z,t,type)  /* integral mean := 0 */
GRID *tGrid;
INT z, t, type;
{
   switch(type){
   case Q_SN: sset_zero_mean(tGrid,z);
        break;
   case Q_SE: set_zero_mean_e(tGrid,z);
        break;
   case Q_SNE: sn_set_zero_mean_e(tGrid,z);
        break;
   case Q_SNSSNSSF: sset_zero_mean(tGrid,z);
        break;
   default:
        eprintf("Error: set_zero_mean not available.\n");
        break;
   }
}

void set_first(tGrid,r,z,t,type)  /* z_1 := r */
GRID *tGrid;
FLOAT r;
INT z, t, type;
{
   switch(type){
   case Q_SN: sset_first(tGrid,r,z);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: set_first not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SF: set_first_f(tGrid,r,z);
        break;
   case Q_SE: set_first_element(tGrid,r,z);
        break;
   case Q_SNE: sn_set_first(tGrid,r,z);
               sn_subtr_one(tGrid,z);
        break;
   case Q_SNSE: sset_first(tGrid,r,z);
                set_first_element(tGrid,r,z);
        break;
   case Q_SFSE: set_first_f(tGrid,r,z);
                set_first_element(tGrid,r,z);
        break;
   default:
        eprintf("Error: set_first not available.\n");
        break;
   }
}

void subtr_first(tGrid,z,t,type)  /* z_i := z_i - z_1 */
GRID *tGrid;
INT z, t, type;
{
   switch(type){
   case Q_SN: ssubtr_first(tGrid,z);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: subtr_first not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SF: subtr_first_f(tGrid,z);
        break;
   case Q_SE: subtr_first_element(tGrid,z);
        break;
   case Q_SNE: sn_subtr_one(tGrid,z);
        break;
   case Q_SNSE: ssubtr_first(tGrid,z);
                subtr_first_element(tGrid,z);
        break;
   case Q_SFSE: subtr_first_f(tGrid,z);
                subtr_first_element(tGrid,z);
        break;
   default:
        eprintf("Error: subtr_first not available.\n");
        break;
   }
}

void copy(tGrid,x,z,t,type)  /* z := x */
GRID *tGrid;
INT x, z, t, type;
{
   if(operation_available(t,&type,"copy"))
   switch(type){
   case Q_SN: scopy(tGrid,x,z,t);
              if (t & ADD_ZN_CONDITION)
                 scopy_BD(tGrid,x,z);
        break;
   case Q_VN: v_copy_n(tGrid,x,z,t);
        break;
   case Q_SF: copy_f(tGrid,x,z,t);
        break;
   case Q_VF: vcopy_f(tGrid,x,z,t);
        break;
   case Q_DVF: dvcopy_f(tGrid,x,z,t);
        break;
   case Q_SE: copy_e(tGrid,x,z);
        break;
   case Q_SNE: sn_copy_e(tGrid,x,z);
        break;
   case Q_VE: vcopy_e(tGrid,x,z);
        break;
   case Q_SNSF: scopy(tGrid,x,z,t);
                copy_f(tGrid,x,z,t);
        break;
   case Q_SNSE: scopy(tGrid,x,z,t);
                copy_e(tGrid,x,z);
        break;
   case Q_SFSE: copy_f(tGrid,x,z,t);
                copy_e(tGrid,x,z);
        break;
   case Q_VNSF: vs_copy_nf(tGrid,x,z,t);
        break;
   case Q_VNVF: v_copy_n(tGrid,x,z,t);
                vcopy_f(tGrid,x,z,t);
        break;
   case Q_VNVE: v_copy_n(tGrid,x,z,t);
                vcopy_e(tGrid,x,z);
        break;
   case Q_VNVFVE: v_copy_n(tGrid,x,z,t);
                  vcopy_f(tGrid,x,z,t);
                  vcopy_e(tGrid,x,z);
        break;
   case Q_MVNMVF: mv_copy_n(tGrid,x,z,t);
                  mv_copy_f(tGrid,x,z,t);
        break;
   case Q_GENERAL: mv_copy_n(tGrid,x,z,t);
                   mv_copy_f(tGrid,x,z,t);
                   mv_copy_e(tGrid,x,z);
        break;
   case Q_SSN: copy_ssn(tGrid,x,z,t);
        break;
   case Q_SSF: copy_ssf(tGrid,x,z,t);
        break;
   case Q_SNSSNSSF: scopy(tGrid,x,z,t);
                    copy_ssn(tGrid,x,z,t);
                    copy_ssf(tGrid,x,z,t);
        break;
   case Q_VNSFLG: vs_copy_nflg(tGrid,x,z,t);
        break;
   default:
        eprintf("Error: copy not available.\n");
        break;
   }
}

void inv(tGrid,x,z,t,type)  /* z := -x */
GRID *tGrid;
INT x, z, t, type;
{
   if(operation_available(t,&type,"inv"))
   switch(type){
   case Q_SN: sinv(tGrid,x,z,t);
              if (t & ADD_ZN_CONDITION)
                 sinv_BD(tGrid,x,z);
        break;
   case Q_VN: v_inv_n(tGrid,x,z,t);
        break;
   case Q_SF: inv_f(tGrid,x,z,t);
        break;
   case Q_VF: vinv_f(tGrid,x,z,t);
        break;
   case Q_DVF: dvinv_f(tGrid,x,z,t);
        break;
   case Q_SE: inv_e(tGrid,x,z);
        break;
   case Q_SNE: sn_inv_e(tGrid,x,z);
        break;
   case Q_VE: vinv_e(tGrid,x,z);
        break;
   case Q_SNSF: sinv(tGrid,x,z,t);
                inv_f(tGrid,x,z,t);
        break;
   case Q_SNSE: sinv(tGrid,x,z,t);
                inv_e(tGrid,x,z);
        break;
   case Q_SFSE: inv_f(tGrid,x,z,t);
                inv_e(tGrid,x,z);
        break;
   case Q_VNSF: vs_inv_nf(tGrid,x,z,t);
        break;
   case Q_VNVF: v_inv_n(tGrid,x,z,t);
                vinv_f(tGrid,x,z,t);
        break;
   case Q_VNVE: v_inv_n(tGrid,x,z,t);
                vinv_e(tGrid,x,z);
        break;
   case Q_VNVFVE: v_inv_n(tGrid,x,z,t);
                  vinv_f(tGrid,x,z,t);
                  vinv_e(tGrid,x,z);
        break;
   case Q_MVNMVF: mv_inv_n(tGrid,x,z,t);
                  mv_inv_f(tGrid,x,z,t);
        break;
   case Q_GENERAL: mv_inv_n(tGrid,x,z,t);
                   mv_inv_f(tGrid,x,z,t);
                   mv_inv_e(tGrid,x,z);
        break;
   case Q_SSN: inv_ssn(tGrid,x,z,t);
        break;
   case Q_SSF: inv_ssf(tGrid,x,z,t);
        break;
   case Q_SNSSNSSF: sinv(tGrid,x,z,t);
                    inv_ssn(tGrid,x,z,t);
                    inv_ssf(tGrid,x,z,t);
        break;
   case Q_VNSFLG: vs_inv_nflg(tGrid,x,z,t);
        break;
   default:
        eprintf("Error: inv not available.\n");
        break;
   }
}

void add(tGrid,x,y,z,t,type)  /* z := x + y */
GRID *tGrid;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"add"))
   switch(type){
   case Q_SN: sadd(tGrid,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: add not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_add_n(tGrid,x,y,z,t);
        break;
   case Q_SF: add_f(tGrid,x,y,z,t);
        break;
   case Q_VF: vadd_f(tGrid,x,y,z,t);
        break;
   case Q_DVF: dvadd_f(tGrid,x,y,z,t);
        break;
   case Q_SE: add_e(tGrid,x,y,z);
        break;
   case Q_SNE: sn_add_e(tGrid,x,y,z);
        break;
   case Q_VE: vadd_e(tGrid,x,y,z);
        break;
   case Q_SNSF: sadd(tGrid,x,y,z,t);
                add_f(tGrid,x,y,z,t);
        break;
   case Q_SNSE: sadd(tGrid,x,y,z,t);
                add_e(tGrid,x,y,z);
        break;
   case Q_SFSE: add_f(tGrid,x,y,z,t);
                add_e(tGrid,x,y,z);
        break;
   case Q_VNSF: vs_add_nf(tGrid,x,y,z,t);
        break;
   case Q_VNVF: v_add_n(tGrid,x,y,z,t);
                vadd_f(tGrid,x,y,z,t);
        break;
   case Q_VNVE: v_add_n(tGrid,x,y,z,t);
                vadd_e(tGrid,x,y,z);
        break;
   case Q_VNVFVE: v_add_n(tGrid,x,y,z,t);
                  vadd_f(tGrid,x,y,z,t);
                  vadd_e(tGrid,x,y,z);
        break;
   case Q_MVNMVF: mv_add_n(tGrid,x,y,z,t);
                  mv_add_f(tGrid,x,y,z,t);
        break;
   case Q_GENERAL: mv_add_n(tGrid,x,y,z,t);
                   mv_add_f(tGrid,x,y,z,t);
                   mv_add_e(tGrid,x,y,z);
        break;
   case Q_SSN: add_ssn(tGrid,x,y,z,t);
        break;
   case Q_SSF: add_ssf(tGrid,x,y,z,t);
        break;
   case Q_SNSSNSSF: sadd(tGrid,x,y,z,t);
                    add_ssn(tGrid,x,y,z,t);
                    add_ssf(tGrid,x,y,z,t);
        break;
   case Q_VNSFLG: vs_add_nflg(tGrid,x,y,z,t);
        break;
   default:
        eprintf("Error: add not available.\n");
        break;
   }
}

void subtr(tGrid,x,y,z,t,type)  /* z := x - y */
GRID *tGrid;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"subtr"))
   switch(type){
   case Q_SN: ssubtr(tGrid,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 ssubtr_BD(tGrid,x,y,z);
        break;
   case Q_VN: v_subtr_n(tGrid,x,y,z,t);
        break;
   case Q_SF: subtr_f(tGrid,x,y,z,t);
        break;
   case Q_VF: vsubtr_f(tGrid,x,y,z,t);
        break;
   case Q_DVF: dvsubtr_f(tGrid,x,y,z,t);
        break;
   case Q_SE: subtr_e(tGrid,x,y,z);
        break;
   case Q_SNE: sn_subtr_e(tGrid,x,y,z);
        break;
   case Q_VE: vsubtr_e(tGrid,x,y,z);
        break;
   case Q_SNSF: ssubtr(tGrid,x,y,z,t);
                subtr_f(tGrid,x,y,z,t);
        break;
   case Q_SNSE: ssubtr(tGrid,x,y,z,t);
                subtr_e(tGrid,x,y,z);
        break;
   case Q_SFSE: subtr_f(tGrid,x,y,z,t);
                subtr_e(tGrid,x,y,z);
        break;
   case Q_VNSF: vs_subtr_nf(tGrid,x,y,z,t);
        break;
   case Q_VNVF: v_subtr_n(tGrid,x,y,z,t);
                vsubtr_f(tGrid,x,y,z,t);
        break;
   case Q_VNVE: v_subtr_n(tGrid,x,y,z,t);
                vsubtr_e(tGrid,x,y,z);
        break;
   case Q_VNVFVE: v_subtr_n(tGrid,x,y,z,t);
                  vsubtr_f(tGrid,x,y,z,t);
                  vsubtr_e(tGrid,x,y,z);
        break;
   case Q_MVNMVF: mv_subtr_n(tGrid,x,y,z,t);
                  mv_subtr_f(tGrid,x,y,z,t);
        break;
   case Q_GENERAL: mv_subtr_n(tGrid,x,y,z,t);
                   mv_subtr_f(tGrid,x,y,z,t);
                   mv_subtr_e(tGrid,x,y,z);
        break;
   case Q_SSN: subtr_ssn(tGrid,x,y,z,t);
        break;
   case Q_SSF: subtr_ssf(tGrid,x,y,z,t);
        break;
   case Q_SNSSNSSF: ssubtr(tGrid,x,y,z,t);
                    subtr_ssn(tGrid,x,y,z,t);
                    subtr_ssf(tGrid,x,y,z,t);
        break;
   case Q_VNSFLG: vs_subtr_nflg(tGrid,x,y,z,t);
        break;
   default:
        eprintf("Error: subtr not available.\n");
        break;
   }
}

void add_and_subtr(tGrid,w,x,y,z,t,type)  /* z := w + x - y */
GRID *tGrid;
INT w, x, y, z, t, type;
{
   if(operation_available(t,&type,"add_and_subtr"))
   switch(type){
   case Q_SN: sadd_and_subtr(tGrid,w,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: add_and_subtr not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_add_and_subtr_n(tGrid,w,x,y,z,t);
        break;
   case Q_SF: add_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_VF: vadd_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_DVF: dvadd_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_SE: add_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_SNE: sn_add_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_VE: vadd_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_SNSF: sadd_and_subtr(tGrid,w,x,y,z,t);
                add_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_SNSE: sadd_and_subtr(tGrid,w,x,y,z,t);
                add_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_SFSE: add_and_subtr_f(tGrid,w,x,y,z,t);
                add_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_VNSF: vs_add_and_subtr_nf(tGrid,w,x,y,z,t);
        break;
   case Q_VNVF: v_add_and_subtr_n(tGrid,w,x,y,z,t);
                vadd_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_VNVE: v_add_and_subtr_n(tGrid,w,x,y,z,t);
                vadd_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_VNVFVE: v_add_and_subtr_n(tGrid,w,x,y,z,t);
                  vadd_and_subtr_f(tGrid,w,x,y,z,t);
                  vadd_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_MVNMVF: mv_add_and_subtr_n(tGrid,w,x,y,z,t);
                  mv_add_and_subtr_f(tGrid,w,x,y,z,t);
        break;
   case Q_GENERAL: mv_add_and_subtr_n(tGrid,w,x,y,z,t);
                   mv_add_and_subtr_f(tGrid,w,x,y,z,t);
                   mv_add_and_subtr_e(tGrid,w,x,y,z);
        break;
   case Q_SSN: add_and_subtr_ssn(tGrid,w,x,y,z,t);
        break;
   case Q_SSF: add_and_subtr_ssf(tGrid,w,x,y,z,t);
        break;
   case Q_SNSSNSSF: sadd_and_subtr(tGrid,w,x,y,z,t);
                    add_and_subtr_ssn(tGrid,w,x,y,z,t);
                    add_and_subtr_ssf(tGrid,w,x,y,z,t);
        break;
   case Q_VNSFLG: vs_add_and_subtr_nflg(tGrid,w,x,y,z,t);
        break;
   default:
        eprintf("Error: add_and_subtr not available.\n");
        break;
   }
}

void mult_and_add(tGrid,r,x,y,z,t,type)  /* z := r*x + y */
GRID *tGrid;
FLOAT r;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"mult_and_add"))
   switch(type){
   case Q_SN: smult_and_add(tGrid,r,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 smult_and_add_BD(tGrid,r,x,y,z);
        break;
   case Q_VN: v_mult_and_add_n(tGrid,r,x,y,z,t);
        break;
   case Q_SF: mult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_VF: vmult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_DVF: dvmult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_SE: mult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_SNE: sn_mult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_VE: vmult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_SNSF: smult_and_add(tGrid,r,x,y,z,t);
                mult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_SNSE: smult_and_add(tGrid,r,x,y,z,t);
                mult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_SFSE: mult_and_add_f(tGrid,r,x,y,z,t);
                mult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_VNSF: vs_mult_and_add_nf(tGrid,r,x,y,z,t);
        break;
   case Q_VNVF: v_mult_and_add_n(tGrid,r,x,y,z,t);
                vmult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_VNVE: v_mult_and_add_n(tGrid,r,x,y,z,t);
                vmult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_VNVFVE: v_mult_and_add_n(tGrid,r,x,y,z,t);
                  vmult_and_add_f(tGrid,r,x,y,z,t);
                  vmult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_MVNMVF: mv_mult_and_add_n(tGrid,r,x,y,z,t);
                  mv_mult_and_add_f(tGrid,r,x,y,z,t);
        break;
   case Q_GENERAL: mv_mult_and_add_n(tGrid,r,x,y,z,t);
                   mv_mult_and_add_f(tGrid,r,x,y,z,t);
                   mv_mult_and_add_e(tGrid,r,x,y,z);
        break;
   case Q_SSN: mult_and_add_ssn(tGrid,r,x,y,z,t);
        break;
   case Q_SSF: mult_and_add_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_SNSSNSSF: smult_and_add(tGrid,r,x,y,z,t);
                    mult_and_add_ssn(tGrid,r,x,y,z,t);
                    mult_and_add_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_VNSFLG: vs_mult_and_add_nflg(tGrid,r,x,y,z,t);
        break;
   default:
        eprintf("Error: mult_and_add not available.\n");
        break;
   }
}

void mult_and_subtr(tGrid,r,x,y,z,t,type)  /* z := r*x - y */
GRID *tGrid;
FLOAT r;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"mult_and_subtr"))
   switch(type){
   case Q_SN: smult_and_subtr(tGrid,r,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 smult_and_subtr_BD(tGrid,r,x,y,z);
        break;
   case Q_VN: v_mult_and_subtr_n(tGrid,r,x,y,z,t);
        break;
   case Q_SF: mult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_VF: vmult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_DVF: dvmult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_SE: mult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_SNE: sn_mult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_VE: vmult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_SNSF: smult_and_subtr(tGrid,r,x,y,z,t);
                mult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_SNSE: smult_and_subtr(tGrid,r,x,y,z,t);
                mult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_SFSE: mult_and_subtr_f(tGrid,r,x,y,z,t);
                mult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_VNSF: vs_mult_and_subtr_nf(tGrid,r,x,y,z,t);
        break;
   case Q_VNVF: v_mult_and_subtr_n(tGrid,r,x,y,z,t);
                vmult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_VNVE: v_mult_and_subtr_n(tGrid,r,x,y,z,t);
                vmult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_VNVFVE: v_mult_and_subtr_n(tGrid,r,x,y,z,t);
                  vmult_and_subtr_f(tGrid,r,x,y,z,t);
                  vmult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_MVNMVF: mv_mult_and_subtr_n(tGrid,r,x,y,z,t);
                  mv_mult_and_subtr_f(tGrid,r,x,y,z,t);
        break;
   case Q_GENERAL: mv_mult_and_subtr_n(tGrid,r,x,y,z,t);
                   mv_mult_and_subtr_f(tGrid,r,x,y,z,t);
                   mv_mult_and_subtr_e(tGrid,r,x,y,z);
        break;
   case Q_SSN: mult_and_subtr_ssn(tGrid,r,x,y,z,t);
        break;
   case Q_SSF: mult_and_subtr_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_SNSSNSSF: smult_and_subtr(tGrid,r,x,y,z,t);
                    mult_and_subtr_ssn(tGrid,r,x,y,z,t);
                    mult_and_subtr_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_VNSFLG: vs_mult_and_subtr_nflg(tGrid,r,x,y,z,t);
        break;
   default:
        eprintf("Error: mult_and_subtr not available.\n");
        break;
   }
}

void damp(tGrid,r,x,y,z,t,type)  /* z := x + r*(y - x) */
GRID *tGrid;
FLOAT r;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"damp"))
   switch(type){
   case Q_SN: sdamp(tGrid,r,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: damp not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_damp_n(tGrid,r,x,y,z,t);
        break;
   case Q_SF: damp_f(tGrid,r,x,y,z,t);
        break;
   case Q_VF: vdamp_f(tGrid,r,x,y,z,t);
        break;
   case Q_DVF: dvdamp_f(tGrid,r,x,y,z,t);
        break;
   case Q_SE: damp_e(tGrid,r,x,y,z);
        break;
   case Q_SNE: sn_damp_e(tGrid,r,x,y,z);
        break;
   case Q_VE: vdamp_e(tGrid,r,x,y,z);
        break;
   case Q_SNSF: sdamp(tGrid,r,x,y,z,t);
                damp_f(tGrid,r,x,y,z,t);
        break;
   case Q_SNSE: sdamp(tGrid,r,x,y,z,t);
                damp_e(tGrid,r,x,y,z);
        break;
   case Q_SFSE: damp_f(tGrid,r,x,y,z,t);
                damp_e(tGrid,r,x,y,z);
        break;
   case Q_VNSF: vs_damp_nf(tGrid,r,x,y,z,t);
        break;
   case Q_VNVF: v_damp_n(tGrid,r,x,y,z,t);
                vdamp_f(tGrid,r,x,y,z,t);
        break;
   case Q_VNVE: v_damp_n(tGrid,r,x,y,z,t);
                vdamp_e(tGrid,r,x,y,z);
        break;
   case Q_VNVFVE: v_damp_n(tGrid,r,x,y,z,t);
                  vdamp_f(tGrid,r,x,y,z,t);
                  vdamp_e(tGrid,r,x,y,z);
        break;
   case Q_MVNMVF: mv_damp_n(tGrid,r,x,y,z,t);
                  mv_damp_f(tGrid,r,x,y,z,t);
        break;
   case Q_GENERAL: mv_damp_n(tGrid,r,x,y,z,t);
                   mv_damp_f(tGrid,r,x,y,z,t);
                   mv_damp_e(tGrid,r,x,y,z);
        break;
   case Q_SSN: damp_ssn(tGrid,r,x,y,z,t);
        break;
   case Q_SSF: damp_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_SNSSNSSF: sdamp(tGrid,r,x,y,z,t);
                    damp_ssn(tGrid,r,x,y,z,t);
                    damp_ssf(tGrid,r,x,y,z,t);
        break;
   case Q_VNSFLG: vs_damp_nflg(tGrid,r,x,y,z,t);
        break;
   default:
        eprintf("Error: dmap not available.\n");
        break;
   }
}

void mult(tGrid,r,x,z,t,type)  /* z := r*x */
GRID *tGrid;
FLOAT r;
INT x, z, t, type;
{
   if(operation_available(t,&type,"mult"))
   switch(type){
   case Q_SN: smult(tGrid,r,x,z,t);
              if (t & ADD_ZN_CONDITION)
                 smult_BD(tGrid,r,x,z);
        break;
   case Q_VN: v_mult_n(tGrid,r,x,z,t);
        break;
   case Q_SF: mult_f(tGrid,r,x,z,t);
        break;
   case Q_VF: vmult_f(tGrid,r,x,z,t);
        break;
   case Q_DVF: dvmult_f(tGrid,r,x,z,t);
        break;
   case Q_SE: mult_e(tGrid,r,x,z);
        break;
   case Q_SNE: sn_mult_e(tGrid,r,x,z);
        break;
   case Q_VE: vmult_e(tGrid,r,x,z);
        break;
   case Q_SNSF: smult(tGrid,r,x,z,t);
                mult_f(tGrid,r,x,z,t);
        break;
   case Q_SNSE: smult(tGrid,r,x,z,t);
                mult_e(tGrid,r,x,z);
        break;
   case Q_SFSE: mult_f(tGrid,r,x,z,t);
                mult_e(tGrid,r,x,z);
        break;
   case Q_VNSF: vs_mult_nf(tGrid,r,x,z,t);
        break;
   case Q_VNVF: v_mult_n(tGrid,r,x,z,t);
                vmult_f(tGrid,r,x,z,t);
        break;
   case Q_VNVE: v_mult_n(tGrid,r,x,z,t);
                vmult_e(tGrid,r,x,z);
        break;
   case Q_VNVFVE: v_mult_n(tGrid,r,x,z,t);
                  vmult_f(tGrid,r,x,z,t);
                  vmult_e(tGrid,r,x,z);
        break;
   case Q_MVNMVF: mv_mult_n(tGrid,r,x,z,t);
                  mv_mult_f(tGrid,r,x,z,t);
        break;
   case Q_GENERAL: mv_mult_n(tGrid,r,x,z,t);
                   mv_mult_f(tGrid,r,x,z,t);
                   mv_mult_e(tGrid,r,x,z);
        break;
   case Q_SSN: mult_ssn(tGrid,r,x,z,t);
        break;
   case Q_SSF: mult_ssf(tGrid,r,x,z,t);
        break;
   case Q_SNSSNSSF: smult(tGrid,r,x,z,t);
                    mult_ssn(tGrid,r,x,z,t);
                    mult_ssf(tGrid,r,x,z,t);
        break;
   case Q_VNSFLG: vs_mult_nflg(tGrid,r,x,z,t);
        break;
   default:
        eprintf("Error: mult not available.\n");
        break;
   }
}

void add_sum_of_multiples(tGrid,r,x,n,z,t,type)  
GRID *tGrid;   /* z += sum_{i=0}^{n-1} r[i]*(x+i) */
DOUBLE *r;
INT x, n, z, t, type;
{
   if(operation_available(t,&type,"add_sum_of_multiples"))
   switch(type){
   case Q_SE: add_sum_of_multiples_e(tGrid,r,x,n,z);
        break;
   case Q_VE: vadd_sum_of_multiples_e(tGrid,r,x,n,z);
        break;
   default:
        eprintf("Error: add_sum_of_multiples not available.\n");
        break;
   }
}

void mult_s(tGrid,r,x,z,t,type)  /* z_i := r_i*x_i, i=1,...,DIM */
GRID *tGrid;
FLOAT r[DIM];
INT x, z, t, type;
{
   if(operation_available(t,&type,"mult_s"))
   switch(type){
   case Q_VE: vmult_e_s(tGrid,r,x,z);
        break;
   default:
        eprintf("Error: mult_s not available.\n");
        break;
   }
}

DOUBLE dot_without_first(tGrid,x,y,t,type)  /* := x.y */
GRID *tGrid;
INT x, y, t, type;
{
   switch(type){
   case Q_SN: return(sdot_without_first(tGrid,x,y));
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: dot_without_first not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_SF: return(dot_f_without_first(tGrid,x,y));
        break;
   case Q_SE: return(dot_e(tGrid,x,y));
        break;
   case Q_SNE: return(sn_dot_e_without_one(tGrid,x,y));
        break;
   case Q_SNSF: return(sdot_without_first(tGrid,x,y)+dot_f(tGrid,x,y,0));
        break;
   case Q_SNSE: return(sdot_without_first(tGrid,x,y)+dot_e(tGrid,x,y));
        break;
   case Q_SFSE: return(dot_f_without_first(tGrid,x,y)+dot_e(tGrid,x,y));
        break;
   default:
        eprintf("Error: dot_without_first not available.\n");
        return(0.);
        break;
   }
}

DOUBLE dot(tGrid,x,y,t,type)  /* := x.y */
GRID *tGrid;
INT x, y, t, type;
{
   if (t & WITHOUT_FIRST)
      return(dot_without_first(tGrid,x,y,t,type));
   if(operation_available(t,&type,"dot"))
   switch(type){
   case Q_SN: if (t & ADD_ZN_CONDITION)
                 return(sdot(tGrid,x,y,t) + sdot_BD(tGrid,x,y));
              else
                 return(sdot(tGrid,x,y,t));
        break;
   case Q_VN: return(v_dot_n(tGrid,x,y,t));
        break;
   case Q_SF: return(dot_f(tGrid,x,y,t));
        break;
   case Q_VF: return(vdot_f(tGrid,x,y,t));
        break;
   case Q_DVF: return(dvdot_f(tGrid,x,y,t));
        break;
   case Q_SE: return(dot_e_with_first(tGrid,x,y));
        break;
   case Q_SNE: return(sn_dot_e(tGrid,x,y));
        break;
   case Q_VE: return(vdot_e(tGrid,x,y));
        break;
   case Q_SNSF: return(sdot(tGrid,x,y,t)+dot_f(tGrid,x,y,t));
        break;
   case Q_SNSE: return(sdot(tGrid,x,y,t)+dot_e_with_first(tGrid,x,y));
        break;
   case Q_SFSE: return(dot_f(tGrid,x,y,t)+dot_e_with_first(tGrid,x,y));
        break;
   case Q_VNSF: return(vs_dot_nf(tGrid,x,y,t));
        break;
   case Q_VNVF: return(v_dot_n(tGrid,x,y,t)+vdot_f(tGrid,x,y,t));
        break;
   case Q_VNVE: return(v_dot_n(tGrid,x,y,t)+vdot_e(tGrid,x,y));
        break;
   case Q_VNVFVE: return(v_dot_n(tGrid,x,y,t)+vdot_f(tGrid,x,y,t)+
                                              vdot_e(tGrid,x,y));
        break;
   case Q_MVNMVF: return(mv_dot_n(tGrid,x,y,t) + mv_dot_f(tGrid,x,y,t));
        break;
   case Q_GENERAL: return(mv_dot_n(tGrid,x,y,t) + mv_dot_f(tGrid,x,y,t) +
                                                  mv_dot_e(tGrid,x,y));
        break;
   case Q_SSN: return(dot_ssn(tGrid,x,y,t));
        break;
   case Q_SSF: return(dot_ssf(tGrid,x,y,t));
        break;
   case Q_SNSSNSSF: return(sdot(tGrid,x,y,t)
                      + dot_ssn(tGrid,x,y,t)
                      + dot_ssf(tGrid,x,y,t));
        break;
   case Q_VNSFLG: return(vs_dot_nflg(tGrid,x,y,t));
        break;
   default:
        eprintf("Error: dot not available.\n");
        return(0.);
        break;
   }
}

void divide(tGrid,x,y,z,t,type)  /*  z_i := x_i/y_i  */
GRID *tGrid;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"divide"))
   switch(type){
   case Q_SN: sdivide(tGrid,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: divide not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_divide_n(tGrid,x,y,z,t);
        break;
   case Q_SF: divide_f(tGrid,x,y,z,t);
        break;
   case Q_VF: vdivide_f(tGrid,x,y,z,t);
        break;
   case Q_DVF: dvdivide_f(tGrid,x,y,z,t);
        break;
   case Q_SE: divide_e(tGrid,x,y,z);
        break;
   case Q_SNE: sn_divide_e(tGrid,x,y,z);
        break;
   case Q_VE: vdivide_e(tGrid,x,y,z);
        break;
   case Q_SNSF: sdivide(tGrid,x,y,z,t);
                divide_f(tGrid,x,y,z,t);
        break;
   case Q_SNSE: sdivide(tGrid,x,y,z,t);
                divide_e(tGrid,x,y,z);
        break;
   case Q_SFSE: divide_f(tGrid,x,y,z,t);
                divide_e(tGrid,x,y,z);
        break;
   case Q_VNSF: vs_divide_nf(tGrid,x,y,z,t);
        break;
   case Q_VNVF: v_divide_n(tGrid,x,y,z,t);
                vdivide_f(tGrid,x,y,z,t);
        break;
   case Q_VNVE: v_divide_n(tGrid,x,y,z,t);
                vdivide_e(tGrid,x,y,z);
        break;
   case Q_VNVFVE: v_divide_n(tGrid,x,y,z,t);
                  vdivide_f(tGrid,x,y,z,t);
                  vdivide_e(tGrid,x,y,z);
        break;
   case Q_MVNMVF: mv_divide_n(tGrid,x,y,z,t);
                  mv_divide_f(tGrid,x,y,z,t);
        break;
   case Q_GENERAL: mv_divide_n(tGrid,x,y,z,t);
                   mv_divide_f(tGrid,x,y,z,t);
                   mv_divide_e(tGrid,x,y,z);
        break;
   case Q_SSN: divide_ssn(tGrid,x,y,z,t);
        break;
   case Q_SSF: divide_ssf(tGrid,x,y,z,t);
        break;
   case Q_SNSSNSSF: sdivide(tGrid,x,y,z,t);
                    divide_ssn(tGrid,x,y,z,t);
                    divide_ssf(tGrid,x,y,z,t);
        break;
   case Q_VNSFLG: vs_divide_nflg(tGrid,x,y,z,t);
        break;
   default:
        eprintf("Error: divide not available.\n");
        break;
   }
}

void multiply(tGrid,x,y,z,t,type)  /*  z_i := x_i*y_i  */
GRID *tGrid;
INT x, y, z, t, type;
{
   if(operation_available(t,&type,"multiply"))
   switch(type){
   case Q_SN: smultiply(tGrid,x,y,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: multiply not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_multiply_n(tGrid,x,y,z,t);
        break;
   case Q_SF: multiply_f(tGrid,x,y,z,t);
        break;
   case Q_VF: vmultiply_f(tGrid,x,y,z,t);
        break;
   case Q_DVF: dvmultiply_f(tGrid,x,y,z,t);
        break;
   case Q_SE: multiply_e(tGrid,x,y,z);
        break;
   case Q_SNE: sn_multiply_e(tGrid,x,y,z);
        break;
   case Q_VE: vmultiply_e(tGrid,x,y,z);
        break;
   case Q_SNSF: smultiply(tGrid,x,y,z,t);
                multiply_f(tGrid,x,y,z,t);
        break;
   case Q_SNSE: smultiply(tGrid,x,y,z,t);
                multiply_e(tGrid,x,y,z);
        break;
   case Q_SFSE: multiply_f(tGrid,x,y,z,t);
                multiply_e(tGrid,x,y,z);
        break;
   case Q_VNSF: vs_multiply_nf(tGrid,x,y,z,t);
        break;
   case Q_VNVF: v_multiply_n(tGrid,x,y,z,t);
                vmultiply_f(tGrid,x,y,z,t);
        break;
   case Q_VNVE: v_multiply_n(tGrid,x,y,z,t);
                vmultiply_e(tGrid,x,y,z);
        break;
   case Q_VNVFVE: v_multiply_n(tGrid,x,y,z,t);
                  vmultiply_f(tGrid,x,y,z,t);
                  vmultiply_e(tGrid,x,y,z);
        break;
   case Q_MVNMVF: mv_multiply_n(tGrid,x,y,z,t);
                  mv_multiply_f(tGrid,x,y,z,t);
        break;
   case Q_GENERAL: mv_multiply_n(tGrid,x,y,z,t);
                   mv_multiply_f(tGrid,x,y,z,t);
                   mv_multiply_e(tGrid,x,y,z);
        break;
   case Q_SSN: multiply_ssn(tGrid,x,y,z,t);
        break;
   case Q_SSF: multiply_ssf(tGrid,x,y,z,t);
        break;
   case Q_SNSSNSSF: smultiply(tGrid,x,y,z,t);
                    multiply_ssn(tGrid,x,y,z,t);
                    multiply_ssf(tGrid,x,y,z,t);
        break;
   case Q_VNSFLG: vs_multiply_nflg(tGrid,x,y,z,t);
        break;
   default:
        eprintf("Error: multiply not available.\n");
        break;
   }
}

void make_sqrt(tGrid,x,z,t,type)  /*  z_i := sqrt(x_i)  */
GRID *tGrid;
INT x, z, t, type;
{
   if(operation_available(t,&type,"make_sqrt"))
   switch(type){
   case Q_SN: smake_sqrt(tGrid,x,z,t);
              if (t & ADD_ZN_CONDITION)
                 eprintf("Error: make_sqrt not available for ADD_ZN_CONDITION.\n");
        break;
   case Q_VN: v_make_sqrt_n(tGrid,x,z,t);
        break;
   case Q_SF: make_sqrt_f(tGrid,x,z,t);
        break;
   case Q_VF: vmake_sqrt_f(tGrid,x,z,t);
        break;
   case Q_DVF: dvmake_sqrt_f(tGrid,x,z,t);
        break;
   case Q_SE: make_sqrt_e(tGrid,x,z);
        break;
   case Q_SNE: sn_make_sqrt_e(tGrid,x,z);
        break;
   case Q_VE: vmake_sqrt_e(tGrid,x,z);
        break;
   case Q_SNSF: smake_sqrt(tGrid,x,z,t);
                make_sqrt_f(tGrid,x,z,t);
        break;
   case Q_SNSE: smake_sqrt(tGrid,x,z,t);
                make_sqrt_e(tGrid,x,z);
        break;
   case Q_SFSE: make_sqrt_f(tGrid,x,z,t);
                make_sqrt_e(tGrid,x,z);
        break;
   case Q_VNSF: vs_make_sqrt_nf(tGrid,x,z,t);
        break;
   case Q_VNVF: v_make_sqrt_n(tGrid,x,z,t);
                vmake_sqrt_f(tGrid,x,z,t);
        break;
   case Q_VNVE: v_make_sqrt_n(tGrid,x,z,t);
                vmake_sqrt_e(tGrid,x,z);
        break;
   case Q_VNVFVE: v_make_sqrt_n(tGrid,x,z,t);
                  vmake_sqrt_f(tGrid,x,z,t);
                  vmake_sqrt_e(tGrid,x,z);
        break;
   case Q_MVNMVF: mv_make_sqrt_n(tGrid,x,z,t);
                  mv_make_sqrt_f(tGrid,x,z,t);
        break;
   case Q_GENERAL: mv_make_sqrt_n(tGrid,x,z,t);
                   mv_make_sqrt_f(tGrid,x,z,t);
                   mv_make_sqrt_e(tGrid,x,z);
        break;
   case Q_SSN: make_sqrt_ssn(tGrid,x,z,t);
        break;
   case Q_SSF: make_sqrt_ssf(tGrid,x,z,t);
        break;
   case Q_SNSSNSSF: smake_sqrt(tGrid,x,z,t);
                    make_sqrt_ssn(tGrid,x,z,t);
                    make_sqrt_ssf(tGrid,x,z,t);
        break;
   case Q_VNSFLG: vs_make_sqrt_nflg(tGrid,x,z,t);
        break;
   default:
        eprintf("Error: make_sqrt not available.\n");
        break;
   }
}

void mult_and_add1_for_GMRES(tGrid,h,l,j,k,t,type)
GRID *tGrid;
FLOAT h[GMN][GMN];
INT l, j, k, t, type;
{
   if(operation_available(t,&type,"mult_and_add1_for_GMRES"))
   switch(type){
   case Q_SN: smult_and_add1_for_GMRES(tGrid,h,l,j,k,t);
              if (t & ADD_ZN_CONDITION)
                 smult_and_add1_for_GMRES_BD(tGrid,h,l,j,k);
        break;
   case Q_VN: v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
        break;
   case Q_SF: mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_VF: vmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_DVF: dvmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_SE: mult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_SNE: sn_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_VE: vmult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_SNSF: smult_and_add1_for_GMRES(tGrid,h,l,j,k,t);
                mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_SNSE: smult_and_add1_for_GMRES(tGrid,h,l,j,k,t);
                mult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_SFSE: mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
                mult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_VNSF: vs_mult_and_add1_for_GMRES_nf(tGrid,h,l,j,k,t);
        break;
   case Q_VNVF: v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
                vmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_VNVE: v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
                vmult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_VNVFVE: v_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
                  vmult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
                  vmult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_MVNMVF: mv_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
                  mv_mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
        break;
   case Q_GENERAL: mv_mult_and_add1_for_GMRES_n(tGrid,h,l,j,k,t);
                   mv_mult_and_add1_for_GMRES_f(tGrid,h,l,j,k,t);
                   mv_mult_and_add1_e_for_GMRES(tGrid,h,l,j,k);
        break;
   case Q_VNSFLG: vs_mult_and_add1_for_GMRES_nflg(tGrid,h,l,j,k,t);
        break;
   default:
        eprintf("Error: mult_and_add1_for_GMRES not available.\n");
        break;
   }
}

void mult_and_add2_for_GMRES(tGrid,g,x,j,k,t,type)
GRID *tGrid;
FLOAT g[GMN];
INT x, j, k, t, type;
{
   if(operation_available(t,&type,"mult_and_add2_for_GMRES"))
   switch(type){
   case Q_SN: smult_and_add2_for_GMRES(tGrid,g,x,j,k,t);
              if (t & ADD_ZN_CONDITION)
                 smult_and_add2_for_GMRES_BD(tGrid,g,x,j,k);
        break;
   case Q_VN: v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
        break;
   case Q_SF: mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_VF: vmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_DVF: dvmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_SE: mult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_SNE: sn_mult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_VE: vmult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_SNSF: smult_and_add2_for_GMRES(tGrid,g,x,j,k,t);
                mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_SNSE: smult_and_add2_for_GMRES(tGrid,g,x,j,k,t);
                mult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_SFSE: mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
                mult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_VNSF: vs_mult_and_add2_for_GMRES_nf(tGrid,g,x,j,k,t);
        break;
   case Q_VNVF: v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
                vmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_VNVE: v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
                vmult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_VNVFVE: v_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
                  vmult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
                  vmult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_MVNMVF: mv_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
                  mv_mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
        break;
   case Q_GENERAL: mv_mult_and_add2_for_GMRES_n(tGrid,g,x,j,k,t);
                   mv_mult_and_add2_for_GMRES_f(tGrid,g,x,j,k,t);
                   mv_mult_and_add2_e_for_GMRES(tGrid,g,x,j,k);
        break;
   case Q_VNSFLG: vs_mult_and_add2_for_GMRES_nflg(tGrid,g,x,j,k,t);
        break;
   default:
        eprintf("Error: mult_and_add2_for_GMRES not available.\n");
        break;
   }
}

void save_vector(tGrid,z,m,nn,nni,nf,nfi,ne,t,type,name)
GRID *tGrid;
INT z, m, nn, nni, nf, nfi, ne, t, type;
char name[];
{
   FILE *fp;

   fp = fopen(name,"w");
   switch(type){
   case Q_SN: save_sn_vector(tGrid,z,m,fp,t);
              save_sn_vector(tGrid,z,m,fp,t);
        break;
   case Q_VN: save_vn_vector(tGrid,z,m,nni,fp);
        break;
   case Q_SF: save_sf_vector(tGrid,z,m,fp,t);
        break;
   case Q_VF: save_vf_vector(tGrid,z,m,nfi,fp);
        break;
   case Q_DVF: save_dvf_vector(tGrid,z,m,nfi,fp);
        break;
   case Q_SE: save_e_vector(tGrid,z,m,fp);
        break;
   case Q_SNE: save_sne_vector(tGrid,z,m,fp);
        break;
   case Q_VE: save_ve_vector(tGrid,z,m,ne,fp);
        break;
   case Q_SNSF: save_sn_vector(tGrid,z,m,fp,t);
                if (t & NSTART_FROM_INNER)
                   m += nni;
                else
                   m += nn;
                save_sf_vector(tGrid,z,m,fp,t);
        break;
   case Q_SNSE: save_sn_vector(tGrid,z,m,fp,t);
                if (t & NSTART_FROM_INNER)
                   m += nni;
                else
                   m += nn;
                save_e_vector(tGrid,z,m,fp);
        break;
   case Q_SFSE: save_sf_vector(tGrid,z,m,fp,t);
                if (t & FSTART_FROM_INNER)
                   m += nfi;
                else
                   m += nf;
                save_e_vector(tGrid,z,m,fp);
        break;
   case Q_VNSF: save_vn_vector(tGrid,z,m,nni,fp);
                m += DIM*nni;
                save_sf_vector(tGrid,z,m,fp,t);
        break;
   case Q_VNVF: save_vn_vector(tGrid,z,m,nni,fp);
                m += DIM*nni;
                save_vf_vector(tGrid,z,m,nfi,fp);
        break;
   case Q_VNVE: save_vn_vector(tGrid,z,m,nni,fp);
                m += DIM*nni;
                save_ve_vector(tGrid,z,m,ne,fp);
        break;
   case Q_VNVFVE: save_vn_vector(tGrid,z,m,nni,fp);
                m += DIM*nni;
                save_vf_vector(tGrid,z,m,nfi,fp);
                m += DIM*nfi;
                save_ve_vector(tGrid,z,m,ne,fp);
        break;
   default:
        eprintf("Error: save_vector not available.\n");
        break;
   }
   fclose(fp);
}

void set_zero_on_2nd_periodic_boundary(tGrid,z,type)
GRID *tGrid;
INT z, type;
{
   P_NODE *pp_node;

   if (type == Q_SN)
      sn_set_zero_on_2nd_periodic_boundary(tGrid,z);
   else
      eprintf("Error: set_zero_on_2nd_periodic_boundary not available.\n");
}

void copy_1st_to_2nd_periodic_boundary(tGrid,z,type)
GRID *tGrid;
INT z, type;
{
   P_NODE *pp_node;

   if (type == Q_SN)
      sn_copy_1st_to_2nd_periodic_boundary(tGrid,z);
   else
      eprintf("Error: copy_to_2nd_periodic_boundary not available.\n");
}

void add_2nd_to_1st_periodic_boundary(tGrid,z,type)
GRID *tGrid;
INT z, type;
{
   P_NODE *pp_node;

   if (type == Q_SN)
      sn_add_2nd_to_1st_periodic_boundary(tGrid,z);
   else
      eprintf("Error: add_2nd_to_1st_periodic_boundary not available.\n");
}
