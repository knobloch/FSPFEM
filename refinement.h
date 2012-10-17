/******************************************************************************/
/*                                                                            */
/*                                 refinement                                 */
/*                                                                            */
/******************************************************************************/

/* refines the finest triangulation which has to consist of elements
   with status = 0 (i.e., all elements of the finest triangulation have to 
   belong to the same level).                                                 */

void Order2(n1,n2,f1,f2)
NODE **n1, **n2;
FACE **f1, **f2;
{
   NODE *n;
   FACE *f;
   
   if ( (*n1)->index > (*n2)->index ){
      EXCHANGE(*n1,*n2,n)
      EXCHANGE(*f1,*f2,f)
   }
}
 
void Order3(n1,n2,n3,f1,f2,f3)
NODE **n1, **n2, **n3;
FACE **f1, **f2, **f3;
{
   Order2(n1,n2,f1,f2);
   Order2(n1,n3,f1,f3);
   Order2(n2,n3,f2,f3);
}

void Order4(n1,n2,n3,n4,f1,f2,f3,f4)
NODE **n1, **n2, **n3, **n4;
FACE **f1, **f2, **f3, **f4;
{
   Order2(n1,n2,f1,f2);
   Order2(n1,n3,f1,f3);
   Order2(n1,n4,f1,f4);
   Order2(n2,n3,f2,f3);
   Order2(n2,n4,f2,f4);
   Order2(n3,n4,f3,f4);
}

LINK *find_link(n1,n2)  /*  n1->index < n2->index       */
NODE *n1, *n2;          /*  (n1,n2) has to be an edge   */
{            
   LINK *pli;

   for (pli=n1->tstart; pli->nbnode!=n2; pli=pli->next);
   return(pli);
}

void prepare_links(mg)
MULTIGRID *mg;
{
   NODE *pnode;
   LINK *pli;
   GRID *theGrid;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pnode = FIRSTN(theGrid); pnode != NULL; pnode = pnode->succ)
         for (pli = pnode->tstart; pli != NULL; pli = pli->next) 
            pli->flag &= 7;
}

void make_nson(n,pnode,pnodeil,pnodebl,nodeNumber)
NODE *n, **pnode, **pnodeil, **pnodebl;
INT  *nodeNumber;
{
   if (NOT_FN(n))
      *pnodebl = (*pnodebl)->succ = *pnode;
   else
      *pnodeil = (*pnodeil)->succ = *pnode;
   n->son = *pnode;
   nmake((*pnode)++,(NODE*)NULL,++(*nodeNumber),n,n->myvertex);
}
 
void new_neighbour(n,neighbour,flag,type,plink)/* type & 1 ... boundary edge  */
NODE *n, *neighbour;                           /* !(type & 1) ... other edge  */
LINK **plink;
INT flag, type;
{
   LINK *pli;
   
   if (IS_FN(neighbour))
      ADD_NEW_LINK(n->start,*plink,pli)
   else
      ADD_NEW_LINK(n->tstart,*plink,pli)
   (*plink)->nbnode = neighbour;
   (*plink)->flag = 0;
   if (flag == 1)
      MARK_EDGE_TO_MIDDLE(*plink);
   else if (flag == 2)
      MARK_EDGE_TO_OLD(*plink);
   if (type & 1)
      MARK_BND_EDGE(*plink);
   ((*plink)++)->next = NULL;
}
 
void new_2neighbours(n,neighbour1,neighbour2,type,plink)
NODE *n, *neighbour1, *neighbour2;
LINK **plink;
INT type;
{
   LINK *pli;
   
   if (IS_FN(neighbour1))
      ADD_NEW_LINK(n->start,*plink,pli)
   else
      ADD_NEW_LINK(n->tstart,*plink,pli)
   (*plink)->next = NULL;
   (*plink)->flag = 0;
   if (type)
      MARK_BND_EDGE(*plink);
   ((*plink)++)->nbnode = neighbour1;

   if (IS_FN(neighbour2))
      ADD_NEW_LINK(n->start,*plink,pli)
   else
      ADD_NEW_LINK(n->tstart,*plink,pli)
   (*plink)->next = NULL;
   (*plink)->flag = 0;
   if (type)
      MARK_BND_EDGE(*plink);
   ((*plink)++)->nbnode = neighbour2;
}
 
void neighbours_of_middle(mnode,n1,n2,plink)  /*  no neighbour of mnode has   */
NODE *mnode, *n1, *n2;                        /*  been stored up to now       */
LINK **plink;
{
   if (IS_FN(n1))
      mnode->start = *plink;
   else
      mnode->tstart = *plink;
   (*plink)->next = NULL;
   (*plink)->flag = 0;
   ((*plink)++)->nbnode = n1;
   if (IS_FN(n2))
      if (mnode->start)
         mnode->start->next = *plink;
      else
         mnode->start = *plink;
   else
      if (mnode->tstart)
         mnode->tstart->next = *plink;
      else
         mnode->tstart = *plink;
   (*plink)->next = NULL;
   (*plink)->flag = 0;
   ((*plink)++)->nbnode = n2;
}

#if DATA_S & N_LINK_TO_FACES

void new_Fneighbour(n,nf,pnflink)
NODE *n;
FACE *nf;
NFLINK **pnflink;
{
   NFLINK *pnf;
                            
   ADD_NEW_LINK(n->nfstart,*pnflink,pnf)
   (*pnflink)->nbface = nf;
   ((*pnflink)++)->next = NULL;
}

void new_F_neighbour(n,nf,pnflink)
NODE *n;
FACE *nf;
NFLINK **pnflink;
{
   NFLINK *pnf;
                            
   if (IS_FF(nf))
      ADD_NEW_LINK(n->nfstart,*pnflink,pnf)
   else
      ADD_NEW_LINK(n->tnfstart,*pnflink,pnf)
   (*pnflink)->nbface = nf;
   ((*pnflink)++)->next = NULL;
}

void new_3Fneighbours(n,nf1,nf2,nf3,pnflink)
NODE *n;
FACE *nf1, *nf2, *nf3;
NFLINK **pnflink;
{
   NFLINK *pnf;
                            
   ADD_NEW_LINK(n->nfstart,*pnflink,pnf)
   (*pnflink)->nbface = nf1;
   ADD_NEXT_FACE(nf2)
   if (IS_FF(nf3)) ADD_NEXT_FACE(nf3)
   ((*pnflink)++)->next = NULL;
   if (NOT_FF(nf3)) ADD_FACE(nf3)
}

void new_3F_neighbours(n,nf1,nf2,nf3,pnflink)
NODE *n;
FACE *nf1, *nf2, *nf3;
NFLINK **pnflink;
{
   NFLINK *pnf;
                            
   ADD_NEW_LINK(n->tnfstart,*pnflink,pnf)
   (*pnflink)->nbface = nf1;
   ADD_NEXT_FACE(nf2)
   ADD_NEXT_FACE(nf3)
   ((*pnflink)++)->next = NULL;
}

void new_4Fneighbours(n,nf1,nf2,nf3,nf4,pnflink)
NODE *n;
FACE *nf1, *nf2, *nf3, *nf4;
NFLINK **pnflink;
{
   ADD_FIRST_FACE2(nf1)
   ADD_NEXT_FACE(nf2)
   if (IS_FF(nf3)) ADD_NEXT_FACE(nf3)
   if (IS_FF(nf4)) ADD_NEXT_FACE(nf4)
   ((*pnflink)++)->next = NULL;
   if (NOT_FF(nf3)) ADD_FACE(nf3)
   if (NOT_FF(nf4)) ADD_FACE(nf4)
}

void new_5Fneighbours(n,nf1,nf2,nf3,nf4,nf5,pnflink)
NODE *n;
FACE *nf1, *nf2, *nf3, *nf4, *nf5;
NFLINK **pnflink;
{
   NFLINK *pnf;
                            
   ADD_NEW_LINK(n->nfstart,*pnflink,pnf)
   (*pnflink)->nbface = nf1;
   ADD_NEXT_FACE(nf2)
   ((*pnflink)++)->next = NULL;
   ADD_FACE(nf3)
   ADD_FACE(nf4)
   ADD_FACE(nf5)
}

void new_8Fneighbours(n,nf1,nf2,nf3,nf4,nf5,nf6,nf7,nf8,pnflink)
NODE *n;
FACE *nf1, *nf2, *nf3, *nf4, *nf5, *nf6, *nf7, *nf8;
NFLINK **pnflink;
{
   ADD_FIRST_FACE2(nf1)
   ADD_NEXT_FACE(nf2)
   ADD_NEXT_FACE(nf3)
   ADD_NEXT_FACE(nf4)
   ADD_NEXT_FACE(nf5)
   ADD_NEXT_FACE(nf6)
   if (IS_FF(nf7)) ADD_NEXT_FACE(nf7)
   if (IS_FF(nf8)) ADD_NEXT_FACE(nf8)
   ((*pnflink)++)->next = NULL;
   if (NOT_FF(nf7)) ADD_FACE(nf7)
   if (NOT_FF(nf8)) ADD_FACE(nf8)
}

#else

void new_Fneighbour(n,nf,pnflink)
NODE *n; FACE *nf; NFLINK **pnflink;
{}

void new_F_neighbour(n,nf,pnflink)
NODE *n; FACE *nf; NFLINK **pnflink;
{}

void new_3Fneighbours(n,nf1,nf2,nf3,pnflink)
NODE *n; FACE *nf1, *nf2, *nf3; NFLINK **pnflink;
{}

void new_3F_neighbours(n,nf1,nf2,nf3,pnflink)
NODE *n; FACE *nf1, *nf2, *nf3; NFLINK **pnflink;
{}

void new_4Fneighbours(n,nf1,nf2,nf3,nf4,pnflink)
NODE *n; FACE *nf1, *nf2, *nf3, *nf4; NFLINK **pnflink;
{}

void new_5Fneighbours(n,nf1,nf2,nf3,nf4,nf5,pnflink)
NODE *n; FACE *nf1, *nf2, *nf3, *nf4, *nf5; NFLINK **pnflink;
{}

void new_8Fneighbours(n,nf1,nf2,nf3,nf4,nf5,nf6,nf7,nf8,pnflink)
NODE *n; FACE *nf1, *nf2, *nf3, *nf4, *nf5, *nf6, *nf7, *nf8; NFLINK **pnflink;
{}

#endif 
 
/* pvert ... new vertex in the middle of the edge (n1,n2);
   pl12  ... link between end nodes on the father level;
   pnode ... the corresponding node.     */
     
NODE *make_vertex(n1,n2,pl12,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber)
NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK *pl12, **plink;
INT  *nodeNumber;
{
   if (IS_BND_EDGE(pl12))
      (*pvert)->type = (NTYPE(n1) & NTYPE(n2)) | 1 | 
                       ((NTYPE(n1)|NTYPE(n2)) & PNMASK);
   else
      (*pvert)->type = 0;
   AVERAGE(n1->myvertex->x,n2->myvertex->x,(*pvert)->x);
   nmake(*pnode,(NODE*)NULL,++(*nodeNumber),(NODE*)NULL,*pvert);
   if(IS_FN(*pnode))
      *pnodeil = (*pnodeil)->succ = *pnode;
   else
      *pnodebl = (*pnodebl)->succ = *pnode;
   neighbours_of_middle(*pnode,n1,n2,plink);
   new_neighbour(n1,*pnode,1,(*pvert)->type,plink);
   new_neighbour(n2,*pnode,1,(*pvert)->type,plink); 
   MARK_DIVIDED_EDGE(pl12);
   (*pnode)++;
   (*pvert)++;
   return(*pnode - 1);
}

NODE *mid_node(n1,n2) /*  looks for the middle of (n1,n2).                    */
NODE *n1, *n2;        /*  The middle has to exist!                            */
{
   LINK *pli, *plin;
   INT t=NSTART_FROM_INNER;
      
   if (NOT_FN(n1) && NOT_FN(n2)){
      if (n1->father->index < n2->father->index)
         pli = find_link(n1->father,n2->father);
      else
         pli = find_link(n2->father,n1->father);
      if (IS_BND_EDGE(pli))
         t = 0;
   }
   for (pli = T_START(n1,t); ; pli = pli->next){
      for(plin = T_START(n2,t); plin && plin->nbnode!=pli->nbnode; 
                                                               plin=plin->next);
      if (plin && plin->nbnode == pli->nbnode && IS_EDGE_TO_MIDDLE(plin) 
                                              && IS_EDGE_TO_MIDDLE(pli)) 
         return(plin->nbnode);
   }
} /*  For uniform refinement, the test IS_EDGE_TO_MIDDLE is not necessary     */
  /*  because sons always have exactly one common node.                       */

NODE *mid_node_or_null(n1,n2)
NODE *n1, *n2;
{
   LINK *pli;
      
   if (n1->father->index < n2->father->index)
      pli = find_link(n1->father,n2->father);
   else
      pli = find_link(n2->father,n1->father);
   if (IS_DIVIDED_EDGE(pli))
      return(mid_node(n1,n2));
   else
      return(NULL);
}
 
NODE *treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;  /*  n1->index < n2->index      */
VERTEX **pvert;                                /*  (n1,n2) has to be an edge  */
LINK **plink;                                  /*  function returns the middle*/
INT  *nodeNumber;                              /*  of this edge               */
{
   LINK *pl12;
   
   for (pl12 = n1->tstart; pl12->nbnode != n2; pl12 = pl12->next);
   if (IS_DIVIDED_EDGE(pl12)) 
      return(mid_node(n1->son,n2->son));
   else{
      if (n1->son == NULL) make_nson(n1,pnode,pnodeil,pnodebl,nodeNumber);
      if (n2->son == NULL) make_nson(n2,pnode,pnodeil,pnodebl,nodeNumber);
      return(make_vertex(n1->son,n2->son,pl12,pvert,pnode,pnodeil,pnodebl,plink,
                                                                   nodeNumber));
   }
}

NODE *Treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;                                /*  (n1,n2) has to be an edge  */
LINK **plink;                                  /*  function returns the middle*/
INT  *nodeNumber;                              /*  of this edge               */
{
   if (n1->index < n2->index)
      return(treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber));
   else
      return(treat_edge(n2,n1,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber));
}

#if DIM == 3

void elem2(theElement,i,status,n1,n2,n3,n4,f1,f2,f3,f4,pelemil,pelembl,pelement)
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
INT i,status;
{
   if (NOT_FF(f1) || NOT_FF(f2) || NOT_FF(f3) || NOT_FF(f4))
      *pelembl = (*pelembl)->succ = *pelement;
   else
      *pelemil = (*pelemil)->succ = *pelement;
   theElement->sons[i] = *pelement;
   Order4(&n1,&n2,&n3,&n4,&f1,&f2,&f3,&f4);
   (*pelement)->status = status;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = n1;
   (*pelement)->n[1] = n2;
   (*pelement)->n[2] = n3;
   (*pelement)->n[3] = n4;
   (*pelement)->f[0] = f1;
   (*pelement)->f[1] = f2;
   (*pelement)->f[2] = f3;
   (*pelement)->f[3] = f4;
   (*pelement)->succ = NULL;
   (*pelement)->father = theElement;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}
      
void treat_4face(f,f1,f2,f3,f0,n12,n13,n23,n11,n22,n33,faceNumber,
                                            pface,pfaceil,pfacebl,plink,pnflink)
FACE *f, **f1, **f2, **f3, **f0, **pface, **pfaceil, **pfacebl; 
NODE *n12, *n13, *n23, *n11, *n22, *n33;
LINK **plink;
NFLINK **pnflink;
INT *faceNumber;
{
   INT i;

   if (f->sons[0] == NULL){
      if (NOT_FF(f))
         *pfacebl = ( (*pfacebl)->succ = *pface ) + 3;
      else
         *pfaceil = ( (*pfaceil)->succ = *pface ) + 3;
      for (i=0; i<4; i++){
         f->sons[i] = *pface;
         fmake(*pface,*pface + 1,FTYPE(f),0,++(*faceNumber),f);
         (*pface)++;
      }
      ((*pface)-1)->succ = NULL;
      new_2neighbours(n12,n13,n23,IS_BF(f),plink);
      new_2neighbours(n13,n23,n12,IS_BF(f),plink);
      new_2neighbours(n23,n12,n13,IS_BF(f),plink);
      if (IS_FF(f)){
         new_3Fneighbours(n12,f->sons[1],f->sons[2],f->sons[0],pnflink);
         new_3Fneighbours(n13,f->sons[1],f->sons[3],f->sons[0],pnflink);
         new_3Fneighbours(n23,f->sons[2],f->sons[3],f->sons[0],pnflink);
      }
      else{
         new_3F_neighbours(n12,f->sons[1],f->sons[2],f->sons[0],pnflink);
         new_3F_neighbours(n13,f->sons[1],f->sons[3],f->sons[0],pnflink);
         new_3F_neighbours(n23,f->sons[2],f->sons[3],f->sons[0],pnflink);
      }
      new_F_neighbour(n11,f->sons[1],pnflink);
      new_F_neighbour(n22,f->sons[2],pnflink);
      new_F_neighbour(n33,f->sons[3],pnflink);
   }   
   *f1 = f->sons[1];
   *f2 = f->sons[2];
   *f3 = f->sons[3];
   *f0 = f->sons[0];
}

/* For uniform refinement, treat_4face instead of treat_4Face can be used. 
   (With treat_4Face, the position of the son faces in the father face is 
   given by nodes order in the new element).                                      */
void treat_4Face(f,f1,f2,f3,f0,n12,n13,n23,n11,n22,n33,faceNumber,
                                            pface,pfaceil,pfacebl,plink,pnflink)
FACE *f, **f1, **f2, **f3, **f0, **pface, **pfaceil, **pfacebl;
NODE *n12, *n13, *n23, *n11, *n22, *n33;
LINK **plink;
NFLINK **pnflink;
INT *faceNumber;
{
   if (n11->index < n22->index){
      if (n22->index < n33->index)        /*   n11 < n22 < n33   */
         treat_4face(f,f1,f2,f3,f0,n12,n13,n23,n11,n22,n33,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
      else if (n11->index < n33->index)   /*   n11 < n33 < n22   */
         treat_4face(f,f1,f3,f2,f0,n13,n12,n23,n11,n33,n22,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
      else                                /*   n33 < n11 < n22   */
         treat_4face(f,f3,f1,f2,f0,n13,n23,n12,n33,n11,n22,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
   }
   else if (n33->index < n22->index)      /*   n33 < n22 < n11   */
      treat_4face(f,f3,f2,f1,f0,n23,n13,n12,n33,n22,n11,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
   else if (n33->index < n11->index)      /*   n22 < n33 < n11   */
      treat_4face(f,f2,f3,f1,f0,n23,n12,n13,n22,n33,n11,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
   else                                   /*   n22 < n11 < n33   */
      treat_4face(f,f2,f1,f3,f0,n12,n23,n13,n22,n11,n33,faceNumber,
                                           pface,pfaceil,pfacebl,plink,pnflink);
}

void iifaces(f11,f22,f33,f44,n11,n22,n33,n44,faceNumber,pface,pfaceil,pnflink)
FACE **f11, **f22, **f33, **f44, **pface, **pfaceil;
NODE *n11, *n22, *n33, *n44;
NFLINK **pnflink;
INT *faceNumber;
{
   INT i;
   
   *f11 = *pface;
   *f22 = *pface + 1;
   *f33 = *pface + 2;
   *f44 = *pface + 3;
   (*pfaceil)->succ = *f11;
   *pfaceil = *f44; 
   for (i=0; i<4; i++){
      fmake(*pface,*pface + 1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
      (*pface)++;
   }
   ((*pface)-1)->succ = NULL;
   new_Fneighbour(n11,*f11,pnflink);
   new_Fneighbour(n22,*f22,pnflink);
   new_Fneighbour(n33,*f33,pnflink);
   new_Fneighbour(n44,*f44,pnflink);
}
  
 void inner_elements(theElement,f11,f22,f33,f44,f10,f20,f30,f40,m1,m2,m3,m4,
          m5,m6,faceNumber,pface,pfaceil,plink,pnflink,pelemil,pelembl,pelement)
 FACE *f11,*f22,*f33,*f44,*f10,*f20,*f30,*f40,**pface,**pfaceil;
 NODE *m1,*m2,*m3,*m4,*m5,*m6;
 LINK **plink;
 NFLINK **pnflink;
 ELEMENT *theElement, **pelemil, **pelembl,**pelement; 
 INT *faceNumber;
 {
    INT i;
    FACE *fi[5]; 
    
    (*pfaceil)->succ = *pface;
    *pfaceil = *pface + 3;
    for (i=1; i<5; i++){
       fi[i] = *pface;
       fmake(fi[i],++(*pface),0,OR_SET,++(*faceNumber),(FACE*)NULL);
    }
    fi[4]->succ = NULL;
    elem2(theElement,4,0,m5,m6,m1,m3,f10,f11,fi[3],fi[1],pelemil,pelembl,pelement);
    elem2(theElement,5,0,m5,m6,m1,m4,f22,f20,fi[4],fi[1],pelemil,pelembl,pelement);
    elem2(theElement,6,0,m5,m6,m2,m3,f33,f30,fi[3],fi[2],pelemil,pelembl,pelement);
    elem2(theElement,7,0,m5,m6,m2,m4,f40,f44,fi[4],fi[2],pelemil,pelembl,pelement);
    new_neighbour(m5,m6,0,0,plink);
    new_neighbour(m6,m5,0,0,plink);
    new_3Fneighbours(m1,fi[1],fi[3],fi[4],pnflink);
    new_3Fneighbours(m2,fi[2],fi[3],fi[4],pnflink);
    new_3Fneighbours(m3,fi[1],fi[2],fi[3],pnflink);
    new_3Fneighbours(m4,fi[1],fi[2],fi[4],pnflink);
    new_8Fneighbours(m5,fi[1],fi[2],fi[3],fi[4],f22,f33,f10,f40,pnflink);
    new_8Fneighbours(m6,fi[1],fi[2],fi[3],fi[4],f11,f44,f20,f30,pnflink);
 }

FLOAT cosA(x1,x2, x3)
FLOAT *x1, *x2, *x3;
{
   FLOAT a[3], b[3];
   
   SUBTR(x1,x2,a);
   SUBTR(x1,x3,b);
   return( DOT(a,b)/sqrt(DOT(a,a)*DOT(b,b)) );
}
 
FLOAT crit1(n1,n2,n3)
NODE *n1, *n2, *n3;
{
   return( fabs(cosA(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x)-0.5)+
           fabs(cosA(n2->myvertex->x,n3->myvertex->x,n1->myvertex->x)-0.5)+
           fabs(cosA(n3->myvertex->x,n1->myvertex->x,n2->myvertex->x)-0.5) );
}
 
 void find_order(n1,n2,n3,i1,i2,i3)  /*  nk is the ik-th node  */
 NODE *n1, *n2, *n3;
 INT *i1, *i2, *i3;
 {
    if (n1->index < n2->index){
       if (n2->index < n3->index){       /*  n1 < n2 < n3  */
          *i1 = 1;
          *i2 = 2;
          *i3 = 3;
       }
       else if (n1->index < n3->index){  /*  n1 < n3 < n2  */
          *i1 = 1;
          *i2 = 3;
          *i3 = 2;
       }
       else{                             /*  n3 < n1 < n2  */
          *i1 = 2;
          *i2 = 3;
          *i3 = 1;
       }
    }
    else if (n3->index < n2->index){     /*  n3 < n2 < n1  */
       *i1 = 3;
       *i2 = 2;
       *i3 = 1;
    }
    else if (n3->index < n1->index){     /*  n2 < n3 < n1  */   
       *i1 = 3;
       *i2 = 1;
       *i3 = 2;
    }
    else{                                /*  n2 < n1 < n3  */
       *i1 = 2;
       *i2 = 1;
       *i3 = 3;
    }
 }
 
INT norm_vect_or(m1,m2,m3,ind)  /* ind == 0  or  ind == ORIENT */
NODE *m1, *m2, *m3;
INT ind;
{
   if (m1->index < m2->index)
      if (m2->index < m3->index || m3->index < m1->index)
         return(ind);
      else if (ind == 0) 
         return(ORIENT);
      else 
         return(0);
   else if (m2->index < m3->index && m3->index < m1->index)
      return(ind);
   else if (ind == 0) 
      return(ORIENT);
   else 
      return(0);
}
 
void set1_norm_vect_or(m1,m2,m3,m4,f)
NODE *m1, *m2, *m3, *m4;
FACE *f;
{
   if (m1->index > m2->index){
      SET_ORIENT(f->sons[1],norm_vect_or(m1,m2,m3,FTYPE(f)&ORIENT));
      SET_ORIENT(f->sons[2],norm_vect_or(m1,m3,m4,FTYPE(f)&ORIENT));
   }
   else{
      SET_ORIENT(f->sons[1],norm_vect_or(m2,m4,m1,FTYPE(f)&ORIENT));
      SET_ORIENT(f->sons[2],norm_vect_or(m2,m3,m4,FTYPE(f)&ORIENT));
   }
}
 
void set2_norm_vect_or(m1,m2,m3,m4,f)
NODE *m1, *m2, *m3, *m4;
FACE *f;
{
   if (m1->index < m2->index){
      SET_ORIENT(f->sons[0],norm_vect_or(m1,m4,m3,FTYPE(f)&ORIENT));
      SET_ORIENT(f->sons[1],norm_vect_or(m4,m2,m3,FTYPE(f)&ORIENT));
   }
   else{
      SET_ORIENT(f->sons[0],norm_vect_or(m4,m2,m3,FTYPE(f)&ORIENT));
      SET_ORIENT(f->sons[1],norm_vect_or(m1,m4,m3,FTYPE(f)&ORIENT));
   }
}
  
void set_norm_vect_or(n11,n22,n33,n12,n13,n23,f)
NODE *n11, *n22, *n33, *n12, *n13, *n23;
FACE *f;
{
   INT i1, i2, i3;
   
   if ((FTYPE(f->sons[0]) & OR_SET) == 0)   
      if (f->sons[3]){
         find_order(n11,n22,n33,&i1,&i2,&i3);
         SET_ORIENT(f->sons[0], norm_vect_or(n12,n23,n13,FTYPE(f)&ORIENT));
         SET_ORIENT(f->sons[i1],norm_vect_or(n11,n12,n13,FTYPE(f)&ORIENT));
         SET_ORIENT(f->sons[i2],norm_vect_or(n12,n22,n23,FTYPE(f)&ORIENT));
         SET_ORIENT(f->sons[i3],norm_vect_or(n13,n23,n33,FTYPE(f)&ORIENT));
      }
      else if (f->sons[2]) /* for nonuniform refinement necessary */
         if (!n12){
            SET_ORIENT(f->sons[0],norm_vect_or(n13,n23,n33,FTYPE(f)&ORIENT));
            set1_norm_vect_or(n11,n22,n23,n13,f);
         }
         else if (!n13){
            SET_ORIENT(f->sons[0],norm_vect_or(n12,n22,n23,FTYPE(f)&ORIENT));
            set1_norm_vect_or(n33,n11,n12,n23,f);
         }   
         else{
            SET_ORIENT(f->sons[0],norm_vect_or(n11,n12,n13,FTYPE(f)&ORIENT));
            set1_norm_vect_or(n22,n33,n13,n12,f);
         }
      else if (f->sons[1])
         if (n12)
            set2_norm_vect_or(n11,n22,n33,n12,f);
         else if (n13)
            set2_norm_vect_or(n33,n11,n22,n13,f);
         else 
            set2_norm_vect_or(n22,n33,n11,n23,f);
      else
         SET_ORIENT(f->sons[0],norm_vect_or(n11,n22,n33,FTYPE(f)&ORIENT));
} 

/* theElement ... the element we want to refine;
   *pelement  ... free place for a new element;
   fij, i=1,...,4, j=1,...,4, j!=i  ...  subface of f[i-1] containing node njj;
   fi0, i=1,...,4  ...  subface of f[i-1] containing no vertex of f[i-1]; 
   fii, i=1,...,4  ...  no-father-face contained in subelement with node nii. */
 
 void finer6(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
          faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,
                                                               pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f11, *f12, *f13, *f14, *f10, *f21, *f22, *f23, *f24, *f20, 
         *f31, *f32, *f33, *f34, *f30, *f41, *f42, *f43, *f44, *f40;
    FLOAT c1, c2, c3;
    INT  i; 
 
    treat_4Face(f1,&f12,&f13,&f14,&f10,n23,n24,n34,n22,n33,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    treat_4Face(f2,&f21,&f23,&f24,&f20,n13,n14,n34,n11,n33,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    treat_4Face(f3,&f31,&f32,&f34,&f30,n12,n14,n24,n11,n22,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    treat_4Face(f4,&f41,&f42,&f43,&f40,n12,n13,n23,n11,n22,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    iifaces(&f11,&f22,&f33,&f44,n11,n22,n33,n44,faceNumber,pface,pfaceil,pnflink);
    elem2(theElement,0,0,n11,n12,n13,n14,f11,f21,f31,f41,pelemil,pelembl,pelement);
    elem2(theElement,1,0,n22,n12,n23,n24,f22,f12,f32,f42,pelemil,pelembl,pelement);
    elem2(theElement,2,0,n33,n13,n23,n34,f33,f13,f23,f43,pelemil,pelembl,pelement);
    elem2(theElement,3,0,n44,n14,n24,n34,f44,f14,f24,f34,pelemil,pelembl,pelement);
    new_4Fneighbours(n12,f11,f22,f12,f21,pnflink);
    new_4Fneighbours(n13,f11,f33,f13,f31,pnflink);
    new_4Fneighbours(n14,f11,f44,f14,f41,pnflink);
    new_4Fneighbours(n23,f22,f33,f23,f32,pnflink);
    new_4Fneighbours(n24,f22,f44,f24,f42,pnflink);
    new_4Fneighbours(n34,f33,f44,f34,f43,pnflink);
    c1 = crit1(n13,n24,n12)+crit1(n13,n24,n34)+crit1(n13,n24,n23)+crit1(n13,n24,n14);
    c2 = crit1(n12,n34,n13)+crit1(n12,n34,n24)+crit1(n12,n34,n14)+crit1(n12,n34,n23);
    c3 = crit1(n14,n23,n13)+crit1(n14,n23,n24)+crit1(n14,n23,n34)+crit1(n14,n23,n12);
    if ( c1<=c2 && c1<=c3 )  
      inner_elements(theElement,f11,f44,f22,f33,f30,f20,f40,f10,n14,n23,n12,n34,
       n13,n24,faceNumber,pface,pfaceil,plink,pnflink,pelemil,pelembl,pelement);
    else if ( c2<=c1 && c2<=c3 )
      inner_elements(theElement,f11,f33,f44,f22,f20,f40,f30,f10,n13,n24,n14,n23,
       n12,n34,faceNumber,pface,pfaceil,plink,pnflink,pelemil,pelembl,pelement);
    else
      inner_elements(theElement,f44,f22,f33,f11,f10,f30,f20,f40,n24,n13,n34,n12,
       n14,n23,faceNumber,pface,pfaceil,plink,pnflink,pelemil,pelembl,pelement);
    for (i=0; i<8; i++){
       nodesToFaces(pfnlink,theElement->sons[i],NMASK,FMASK);
       facesToFaces(pflink,theElement->sons[i],FMASK); 
    }    
 }

 void unif_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,
    plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 NODE **pnode, **pnodeil, **pnodebl;
 FACE **pface, **pfaceil, **pfacebl;
 VERTEX **pvert;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT  *nodeNumber, *faceNumber; 
 {
    NODE *n1, *n2, *n3, *n4, *n12, *n13, *n14, *n23, *n24, *n34,
         *n11, *n22, *n33, *n44;
    
    NODES_OF_ELEMENT(n1,n2,n3,n4,theElement);
    n12 = treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    n13 = treat_edge(n1,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    n14 = treat_edge(n1,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    n23 = treat_edge(n2,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    n24 = treat_edge(n2,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    n34 = treat_edge(n3,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    SONS_OF_NODES4(n1,n2,n3,n4,n11,n22,n33,n44)
    finer6(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,
        theElement->f[0],theElement->f[1],theElement->f[2],theElement->f[3],
        faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,
                                                              pelembl,pelement);
    set_norm_vect_or(n11,n22,n33,n12,n13,n23,theElement->f[3]); 
    set_norm_vect_or(n11,n22,n44,n12,n14,n24,theElement->f[2]); 
    set_norm_vect_or(n11,n33,n44,n13,n14,n34,theElement->f[1]); 
    set_norm_vect_or(n22,n33,n44,n23,n24,n34,theElement->f[0]); 
 }

#else /*  DIM == 2  */

#if ELEMENT_TYPE == SIMPLEX

void elem2(theElement,i,status,n1,n2,n3,f1,f2,f3,pelemil,pelembl,pelement)
NODE *n1, *n2, *n3;
FACE *f1, *f2, *f3;
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
INT i,status;
{
   if (NOT_FF(f1) || NOT_FF(f2) || NOT_FF(f3))
      *pelembl = (*pelembl)->succ = *pelement;
   else
      *pelemil = (*pelemil)->succ = *pelement;
   theElement->sons[i] = *pelement;
   Order3(&n1,&n2,&n3,&f1,&f2,&f3);
   (*pelement)->status = status;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = n1;
   (*pelement)->n[1] = n2;
   (*pelement)->n[2] = n3;
   (*pelement)->f[0] = f1;
   (*pelement)->f[1] = f2;
   (*pelement)->f[2] = f3;
   (*pelement)->succ = NULL;
   (*pelement)->father = theElement;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}

#elif ELEMENT_TYPE == CUBE

void elem2(theElement,i,status,n1,n2,n3,n4,f1,f2,f3,f4,pelemil,pelembl,pelement)
NODE *n1, *n2, *n3, *n4;
FACE *f1, *f2, *f3, *f4;
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
INT i,status;
{
   if (NOT_FF(f1) || NOT_FF(f2) || NOT_FF(f3) || NOT_FF(f4))
      *pelembl = (*pelembl)->succ = *pelement;
   else
      *pelemil = (*pelemil)->succ = *pelement;
   theElement->sons[i] = *pelement;
   (*pelement)->status = status;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = n1;
   (*pelement)->n[1] = n2;
   (*pelement)->n[2] = n3;
   (*pelement)->n[3] = n4;
   (*pelement)->f[0] = f1;
   (*pelement)->f[1] = f2;
   (*pelement)->f[2] = f3;
   (*pelement)->f[3] = f4;
   (*pelement)->succ = NULL;
   (*pelement)->father = theElement;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}

#endif
      
void treat_1face(f3,f33,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink,plink)
FACE *f3, **f33, **pface, **pfaceil, **pfacebl; 
NODE *n11, *n22;
NFLINK **pnflink;
LINK **plink;
INT *faceNumber;
{
   LINK *pl12;
   INT type=0;

   if (f3->sons[0] == NULL){
      if (NOT_FF(f3))
         *pfacebl = (*pfacebl)->succ = *pface;
      else
         *pfaceil = (*pfaceil)->succ = *pface;
      f3->sons[0] = *pface;
      fmake(*pface,NULL,FTYPE(f3),0,++(*faceNumber),f3);
      (*pface)++;
      f3->sons[1] = NULL;
      new_F_neighbour(n11,f3->sons[0],pnflink);
      new_F_neighbour(n22,f3->sons[0],pnflink);
      if (n11->father->index < n22->father->index)
         pl12 = find_link(n11->father,n22->father);
      else
         pl12 = find_link(n22->father,n11->father);
      if (IS_BND_EDGE(pl12)) type = 1;
      new_neighbour(n11,n22,0,type,plink);
      new_neighbour(n22,n11,0,type,plink);
   }   
   *f33 = f3->sons[0];
}

void treat_2face(f3,f31,f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink)
FACE *f3, **f31, **f32, **pface, **pfaceil, **pfacebl; 
NODE *n12, *n11, *n22;
NFLINK **pnflink;
INT *faceNumber;
{
   if (f3->sons[0] == NULL){
      if (NOT_FF(f3))
         *pfacebl = ( (*pfacebl)->succ = *pface ) + 1;
      else
         *pfaceil = ( (*pfaceil)->succ = *pface ) + 1;
      f3->sons[0] = *pface;
      fmake(*pface,*pface + 1,FTYPE(f3),0,++(*faceNumber),f3);
      (*pface)++;
      f3->sons[1] = *pface;
      fmake(*pface,NULL,FTYPE(f3),0,++(*faceNumber),f3);
      (*pface)++;
      new_F_neighbour(n11,f3->sons[0],pnflink);
      new_F_neighbour(n22,f3->sons[1],pnflink);
      new_F_neighbour(n12,f3->sons[0],pnflink);
      new_F_neighbour(n12,f3->sons[1],pnflink);
   }   
   *f31 = f3->sons[0];
   *f32 = f3->sons[1];
}

void treat_2Face(f3,f31,f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink)
FACE *f3, **f31, **f32, **pface, **pfaceil, **pfacebl; 
NODE *n12, *n11, *n22;
NFLINK **pnflink;
INT *faceNumber;
{
   if (n11->index < n22->index)
      treat_2face(f3,f31,f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink);
   else
      treat_2face(f3,f32,f31,n12,n22,n11,faceNumber,pface,pfaceil,pfacebl,pnflink);
}

void iifaces0(f11,f22,f33,faceNumber,pface,pfaceil,pnflink)
FACE **f11, **f22, **f33, **pface, **pfaceil;
NFLINK **pnflink;
INT *faceNumber;
{
   INT i;
   
   *f11 = *pface;
   *f22 = *pface + 1;
   *f33 = *pface + 2;
   (*pfaceil)->succ = *f11;
   *pfaceil = *f33; 
   for (i=0; i<3; i++){
      fmake(*pface,*pface + 1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
      (*pface)++;
   }
   ((*pface)-1)->succ = NULL;
}

void iifaces(f11,f22,f33,n11,n22,n33,faceNumber,pface,pfaceil,pnflink)
FACE **f11, **f22, **f33, **pface, **pfaceil;
NODE *n11, *n22, *n33;
NFLINK **pnflink;
INT *faceNumber;
{
   iifaces0(f11,f22,f33,faceNumber,pface,pfaceil,pnflink);
   new_Fneighbour(n11,*f11,pnflink);
   new_Fneighbour(n22,*f22,pnflink);
   new_Fneighbour(n33,*f33,pnflink);
}

void find_order(n1,n2,i1,i2)
NODE *n1, *n2;
INT *i1, *i2;
{
   if (n1->index < n2->index){
      *i1 = 0;
      *i2 = 1;
   }
   else{
      *i1 = 1;
      *i2 = 0; 
   }
}
 
INT norm_vect_or(m1,m2,ind)  /* ind == 0  or  ind == ORIENT */
NODE *m1, *m2;
INT ind;
{
   if (m1->index < m2->index)
      return(ind);
   else if (ind == 0) 
      return(ORIENT);
   else 
      return(0);
}
 
void set_norm_vect_or(n11,n22,n12,f)
NODE *n11, *n22, *n12;
FACE *f;
{
   NODE *n;
   INT i1, i2;
   
   if (n11->father->index > n22->father->index)
      EXCHANGE(n11,n22,n)
   if ((FTYPE(f->sons[0]) & OR_SET) == 0){
      if (f->sons[1]){
         find_order(n11,n22,&i1,&i2);
         SET_ORIENT(f->sons[i1],norm_vect_or(n11,n12,FTYPE(f)&ORIENT));
         SET_ORIENT(f->sons[i2],norm_vect_or(n12,n22,FTYPE(f)&ORIENT));
      }
      else
         SET_ORIENT(f->sons[0],norm_vect_or(n11,n22,FTYPE(f)&ORIENT));
   }
} 

#if ELEMENT_TYPE == SIMPLEX

/* theElement ... the element we want to refine;
   *pelement  ... free place for a new element;
   fij, i=1,2,3, j=1,2,3, j!=i  ...  subface of f[i-1] containing node njj;
   fii, i=1,2,3  ...  no-father-face contained in subelement with node nii.   */
 
void finer4(theElement,n11,n22,n33,n12,n13,n23,f1,f2,f3,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n12, *n13, *n23;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *faceNumber;
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   INT  i; 
 
   treat_2Face(f1,&f12,&f13,n23,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f2,&f21,&f23,n13,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f3,&f31,&f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink);
   iifaces(&f11,&f22,&f33,n11,n22,n33,faceNumber,pface,pfaceil,pnflink);
   elem2(theElement,0,0,n11,n12,n13,f11,f21,f31,pelemil,pelembl,pelement);
   elem2(theElement,1,0,n22,n12,n23,f22,f12,f32,pelemil,pelembl,pelement);
   elem2(theElement,2,0,n33,n13,n23,f33,f13,f23,pelemil,pelembl,pelement);
   elem2(theElement,3,0,n12,n23,n13,f33,f11,f22,pelemil,pelembl,pelement);
   new_2neighbours(n12,n13,n23,0,plink);
   new_2neighbours(n23,n12,n13,0,plink);
   new_2neighbours(n13,n12,n23,0,plink);
   new_5Fneighbours(n12,f11,f22,f33,f12,f21,pnflink);
   new_5Fneighbours(n13,f11,f22,f33,f13,f31,pnflink);
   new_5Fneighbours(n23,f11,f22,f33,f23,f32,pnflink);
   for (i=0; i<4; i++){
      nodesToFaces(pfnlink,theElement->sons[i],NMASK,FMASK);
      facesToFaces(pflink,theElement->sons[i],FMASK); 
   }    
}

#elif ELEMENT_TYPE == CUBE

void finer4(theElement,n11,n22,n33,n44,n12,n23,n34,n41,ni,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n44, *n12, *n23, *n34, *n41, *ni;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *faceNumber;
{
   FACE *f112, *f122, *f223, *f233, *f334, *f344, *f441, *f411, 
        *f12i, *f23i, *f34i, *f41i;
   INT  i; 
 
   treat_2Face(f1,&f112,&f122,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f2,&f223,&f233,n23,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f3,&f334,&f344,n34,n33,n44,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f4,&f441,&f411,n41,n44,n11,faceNumber,pface,pfaceil,pfacebl,pnflink);
   f12i = *pface;
   f23i = *pface + 1;
   f34i = *pface + 2;
   f41i = *pface + 3;
   (*pfaceil)->succ = f12i;
   *pfaceil = f41i; 
   for (i = 0; i < 4; i++){
      fmake(*pface,*pface + 1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
      (*pface)++;
   }
   ((*pface)-1)->succ = NULL;
   elem2(theElement,0,0,n41,n11,n12,ni,f411,f112,f12i,f41i,pelemil,pelembl,pelement);
   elem2(theElement,1,0,n12,n22,n23,ni,f122,f223,f23i,f12i,pelemil,pelembl,pelement);
   elem2(theElement,2,0,n23,n33,n34,ni,f233,f334,f34i,f23i,pelemil,pelembl,pelement);
   elem2(theElement,3,0,n34,n44,n41,ni,f344,f441,f41i,f34i,pelemil,pelembl,pelement);
   new_neighbour(n11,ni,0,0,plink);
   new_neighbour(n22,ni,0,0,plink);
   new_neighbour(n33,ni,0,0,plink);
   new_neighbour(n44,ni,0,0,plink);
   new_neighbour(n12,ni,0,0,plink);
   new_neighbour(n23,ni,0,0,plink);
   new_neighbour(n34,ni,0,0,plink);
   new_neighbour(n41,ni,0,0,plink);
   new_2neighbours(n12,n41,n23,0,plink);
   new_2neighbours(n23,n12,n34,0,plink);
   new_2neighbours(n34,n23,n41,0,plink);
   new_2neighbours(n41,n34,n12,0,plink);
   new_2neighbours(ni,n11,n12,0,plink);
   new_2neighbours(ni,n22,n23,0,plink);
   new_2neighbours(ni,n33,n34,0,plink);
   new_2neighbours(ni,n44,n41,0,plink);
   new_Fneighbour(n11,f41i,pnflink);
   new_Fneighbour(n11,f12i,pnflink);
   new_Fneighbour(n22,f12i,pnflink);
   new_Fneighbour(n22,f23i,pnflink);
   new_Fneighbour(n33,f23i,pnflink);
   new_Fneighbour(n33,f34i,pnflink);
   new_Fneighbour(n44,f34i,pnflink);
   new_Fneighbour(n44,f41i,pnflink);
   new_5Fneighbours(n12,f41i,f12i,f23i,f411,f223,pnflink);
   new_5Fneighbours(n23,f12i,f23i,f34i,f122,f334,pnflink);
   new_5Fneighbours(n34,f23i,f34i,f41i,f233,f441,pnflink);
   new_5Fneighbours(n41,f34i,f41i,f12i,f344,f112,pnflink);
   new_5Fneighbours(ni,f12i,f23i,f112,f122,f223,pnflink);
   new_5Fneighbours(ni,f34i,f41i,f233,f334,f344,pnflink);
   new_F_neighbour(ni,f441,pnflink);
   new_F_neighbour(ni,f411,pnflink);
   for (i = 0; i < 4; i++){
      nodesToFaces(pfnlink,theElement->sons[i],NMASK,FMASK);
      facesToFaces(pflink,theElement->sons[i],FMASK); 
   }    
}

#endif

#if ELEMENT_TYPE == SIMPLEX

void sub_elem1(theElem,i,status,n11,n12,n13,f11,f21,f31,
                                               pnflink,pelemil,pelembl,pelement)
ELEMENT *theElem, **pelemil, **pelembl, **pelement;
FACE *f11, *f21, *f31;
NODE *n11, *n12, *n13;
NFLINK **pnflink;
INT i, status;
{
   new_Fneighbour(n11,f11,pnflink);
   new_F_neighbour(n12,f21,pnflink);
   new_F_neighbour(n13,f31,pnflink);
   elem2(theElem,i,status,n11,n12,n13,f11,f21,f31,pelemil,pelembl,pelement);
}

void sub_elem2(theElem,i1,i2,n11,n12,n13,n112,f11,f21,f31o,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElem, **pelemil, **pelembl, **pelement;
FACE *f11, *f21, *f31o, **pface, **pfaceil, **pfacebl;
NODE *n11, *n12, *n13, *n112;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT i1, i2, *faceNumber;
{
   FACE *f311, *f312, *f1i;

   treat_2Face(f31o,&f311,&f312,n112,n11,n12,
                                      faceNumber,pface,pfaceil,pfacebl,pnflink);
   *pfaceil = (*pfaceil)->succ = f1i = (*pface)++;
   fmake(f1i,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   new_neighbour(n112,n13,0,0,plink);
   new_neighbour(n13,n112,0,0,plink);
   new_Fneighbour(n11,f1i,pnflink);
   new_Fneighbour(n12,f1i,pnflink);
   new_Fneighbour(n13,f1i,pnflink);
   new_Fneighbour(n112,f1i,pnflink);
   new_F_neighbour(n112,f11,pnflink);
   new_F_neighbour(n112,f21,pnflink);
   new_F_neighbour(n13,f311,pnflink);
   new_F_neighbour(n13,f312,pnflink);
   elem2(theElem,i1,1,n112,n13,n11,f21,f311,f1i,pelemil,pelembl,pelement);
   elem2(theElem,i2,1,n112,n13,n12,f11,f312,f1i,pelemil,pelembl,pelement);
}

void sub_elem3(theElem,n12,n22,n23,n122,n223,f22,f12o,f32o,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElem, **pelemil, **pelembl, **pelement;
FACE *f22, *f12o, *f32o, **pface, **pfaceil, **pfacebl;
NODE *n12, *n22, *n23, *n122, *n223;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *faceNumber;
{
   FACE *f321, *f322, *f122, *f123, *f222, *f2i;

   treat_2Face(f12o,&f122,&f123,n223,n22,n23,
                                      faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f32o,&f321,&f322,n122,n12,n22,
                                      faceNumber,pface,pfaceil,pfacebl,pnflink);
   (*pfaceil)->succ = f2i = (*pface)++;
   *pfaceil = f222 = (*pface)++;
   fmake(f2i,f222,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   fmake(f222,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   new_neighbour(n23,n122,0,0,plink);
   new_neighbour(n223,n122,0,0,plink);
   new_2neighbours(n122,n23,n223,0,plink);
   new_5Fneighbours(n122,f22,f2i,f222,f122,f123,pnflink);
   new_3Fneighbours(n223,f2i,f222,f322,pnflink);
   new_3Fneighbours(n23,f2i,f222,f321,pnflink);
   new_Fneighbour(n22,f222,pnflink);
   new_Fneighbour(n12,f2i,pnflink);
   elem2(theElem,0,1,n12,n122,n23,f2i,f22,f321,pelemil,pelembl,pelement);
   elem2(theElem,1,1,n122,n22,n223,f122,f222,f322,pelemil,pelembl,pelement);
   elem2(theElem,2,1,n122,n223,n23,f123,f2i,f222,pelemil,pelembl,pelement);
}

void finer4_green2(theElem1,theElem2,n11,n22,n33,n12,n13,n23,n112,n122,
                   f1o,f2o,f31o,f32o,faceNumber,pface,pfaceil,pfacebl,plink,
                   pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElem1, *theElem2, **pelemil, **pelembl, **pelement;
FACE *f1o, *f2o, *f31o, *f32o, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n12, *n13, *n23, *n112, *n122;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *faceNumber;
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f33, *f31, *f32;
   INT  i; 
 
   treat_2Face(f1o,&f12,&f13,n23,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f2o,&f21,&f23,n13,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   iifaces0(&f11,&f22,&f33,faceNumber,pface,pfaceil,pnflink);
   theElem1->sons[2] = theElem1->sons[3] = NULL;
   theElem2->sons[2] = theElem2->sons[3] = NULL;
   if (n112)
      sub_elem2(theElem1,1,2,n11,n12,n13,n112,f11,f21,f31o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   else{
      treat_1face(f31o,&f31,n11,n12,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem1(theElem1,1,0,n11,n12,n13,f11,f21,f31,
                                              pnflink,pelemil,pelembl,pelement);
   }
   if (n122)
      sub_elem2(theElem2,1,2,n12,n22,n23,n122,f12,f22,f32o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   else{
      treat_1face(f32o,&f32,n12,n22,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem1(theElem2,1,0,n22,n12,n23,f22,f12,f32,
                                              pnflink,pelemil,pelembl,pelement);
   }
   sub_elem1(theElem1,0,0,n33,n13,n23,f33,f13,f23,
                                              pnflink,pelemil,pelembl,pelement);
   elem2(theElem2,0,0,n12,n23,n13,f33,f11,f22,pelemil,pelembl,pelement);
   new_2neighbours(n12,n13,n23,0,plink);
   new_2neighbours(n23,n12,n13,0,plink);
   new_2neighbours(n13,n12,n23,0,plink);
   new_3Fneighbours(n12,f11,f22,f33,pnflink);
   new_3Fneighbours(n13,f11,f22,f33,pnflink);
   new_3Fneighbours(n23,f11,f22,f33,pnflink);
   for (i=0; i<3; i++){
      if (theElem1->sons[i]){
         nodesToFaces(pfnlink,theElem1->sons[i],NMASK,FMASK);
         facesToFaces(pflink,theElem1->sons[i],FMASK); 
      }
      if (theElem2->sons[i]){
         nodesToFaces(pfnlink,theElem2->sons[i],NMASK,FMASK);
         facesToFaces(pflink,theElem2->sons[i],FMASK); 
      }
   }
}

void finer4_green3(theElem1,theElem2,theElem3,n11,n22,n33,n12,n13,n23,n112,n122,
          n223,n233,f12o,f13o,f2o,f31o,f32o,f22o,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElem1, *theElem2, *theElem3, **pelemil, **pelembl, **pelement;
FACE *f12o, *f13o, *f2o, *f31o, *f32o, *f22o, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n12, *n13, *n23, *n112, *n122, *n223, *n233;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *faceNumber;
{
   FACE *f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33;
   INT  i; 
 
   treat_2Face(f2o,&f21,&f23,n13,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   iifaces0(&f11,&f22,&f33,faceNumber,pface,pfaceil,pnflink);
   theElem1->sons[2] = theElem1->sons[3] = NULL;
   theElem2->sons[1] = theElem2->sons[2] = theElem2->sons[3] = NULL;
   theElem3->sons[1] = theElem3->sons[2] = theElem3->sons[3] = NULL;
   if (n112)
      sub_elem2(theElem1,1,2,n11,n12,n13,n112,f11,f21,f31o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   else{
      treat_1face(f31o,&f31,n11,n12,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem1(theElem1,1,0,n11,n12,n13,f11,f21,f31,
                                              pnflink,pelemil,pelembl,pelement);
   }
   if (n233)
      sub_elem2(theElem3,0,1,n23,n33,n13,n233,f23,f33,f13o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   else{
      treat_1face(f13o,&f13,n23,n33,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem1(theElem3,0,0,n33,n13,n23,f33,f13,f23,
                                              pnflink,pelemil,pelembl,pelement);
   }
   if (n122 && n223)
      sub_elem3(theElem2,n12,n22,n23,n122,n223,f22,f12o,f32o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   else if (n122 && !n223){
      treat_1face(f12o,&f12,n22,n23,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem2(theElem2,0,1,n12,n22,n23,n122,f12,f22,f32o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   }
   else if (!n122 && n223){
      treat_1face(f32o,&f32,n12,n22,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem2(theElem2,0,1,n22,n23,n12,n223,f22,f32,f12o,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
   }
   else{
      treat_1face(f32o,&f32,n12,n22,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      treat_1face(f12o,&f12,n22,n23,
                                faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
      sub_elem1(theElem2,0,0,n22,n12,n23,f22,f12,f32,
                                              pnflink,pelemil,pelembl,pelement);
   }
   elem2(theElem1,0,0,n12,n23,n13,f33,f11,f22,pelemil,pelembl,pelement);
   new_2neighbours(n12,n13,n23,0,plink);
   new_2neighbours(n23,n12,n13,0,plink);
   new_2neighbours(n13,n12,n23,0,plink);
   new_3Fneighbours(n12,f11,f22,f33,pnflink);
   new_3Fneighbours(n13,f11,f22,f33,pnflink);
   new_3Fneighbours(n23,f11,f22,f33,pnflink);
   f22o->sons[0] = f22;
   f22->father = f22o;
   for (i=0; i<3; i++){
      if (theElem1->sons[i]){
         nodesToFaces(pfnlink,theElem1->sons[i],NMASK,FMASK);
         facesToFaces(pflink,theElem1->sons[i],FMASK); 
      }
      if (theElem2->sons[i]){
         nodesToFaces(pfnlink,theElem2->sons[i],NMASK,FMASK);
         facesToFaces(pflink,theElem2->sons[i],FMASK); 
      }
      if (theElem3->sons[i]){
         nodesToFaces(pfnlink,theElem3->sons[i],NMASK,FMASK);
         facesToFaces(pflink,theElem3->sons[i],FMASK); 
      }
   }
}

void finer3(theElement,n11,n22,n33,n12,n23,f1,f2,f3,status,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n12, *n23;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT status, *faceNumber;
{
   FACE *f12, *f13, *f22, *f2i, *f31, *f32, *f33;
 
   treat_2Face(f1,&f12,&f13,n23,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_2Face(f3,&f31,&f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink);
   treat_1face(f2,&f22,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   fmake(*pface,*pface + 1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   (*pfaceil)->succ = f33 = (*pface)++;
   fmake(*pface,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   *pfaceil = f2i = (*pface)++;
   new_Fneighbour(n11,f33,pnflink);
   new_Fneighbour(n22,f2i,pnflink);
   new_Fneighbour(n33,f33,pnflink);
   new_Fneighbour(n33,f2i,pnflink);
   elem2(theElement,0,status,n11,n12,n33,f33,f22,f31,pelemil,pelembl,pelement);
   elem2(theElement,1,status,n22,n12,n23,f2i,f12,f32,pelemil,pelembl,pelement);
   elem2(theElement,2,status,n33,n12,n23,f2i,f13,f33,pelemil,pelembl,pelement);
   theElement->sons[3] = NULL;
   new_2neighbours(n12,n33,n23,0,plink);
   new_neighbour(n23,n12,0,0,plink);
   new_neighbour(n33,n12,0,0,plink);
   new_5Fneighbours(n12,f33,f2i,f22,f12,f13,pnflink);
   new_3Fneighbours(n23,f33,f2i,f32,pnflink);
   new_F_neighbour(n33,f31,pnflink);
   nodesToFaces(pfnlink,theElement->sons[0],NMASK,FMASK);
   nodesToFaces(pfnlink,theElement->sons[1],NMASK,FMASK);
   nodesToFaces(pfnlink,theElement->sons[2],NMASK,FMASK);
   facesToFaces(pflink,theElement->sons[0],FMASK); 
   facesToFaces(pflink,theElement->sons[1],FMASK); 
   facesToFaces(pflink,theElement->sons[2],FMASK); 
}

void finer2(theElement,n11,n22,n3,n12,f1,f2,f3,status,nodeNumber,faceNumber,
            pnode,pnodeil,pnodebl,pvert,
            pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
            pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n3, *n12, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT status, *nodeNumber, *faceNumber;
{
   NODE *n33;
   FACE *f11, *f22, *f33, *f31, *f32;

   if (n3->son == NULL) make_nson(n3,pnode,pnodeil,pnodebl,nodeNumber);
   n33 = n3->son;
   treat_1face(f1,&f11,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_1face(f2,&f22,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_2Face(f3,&f31,&f32,n12,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink);
   fmake(*pface,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   *pfaceil = (*pfaceil)->succ = f33 = (*pface)++;
   elem2(theElement,0,status,n11,n12,n33,f33,f22,f31,pelemil,pelembl,pelement);
   elem2(theElement,1,status,n22,n12,n33,f33,f11,f32,pelemil,pelembl,pelement);
   theElement->sons[2] = theElement->sons[3] = NULL;
   new_Fneighbour(n11,f33,pnflink);
   new_Fneighbour(n22,f33,pnflink);
   new_Fneighbour(n33,f33,pnflink);
   new_Fneighbour(n12,f33,pnflink);
   new_neighbour(n12,n33,0,0,plink);
   new_neighbour(n33,n12,0,0,plink);
   new_F_neighbour(n33,f31,pnflink);
   new_F_neighbour(n33,f32,pnflink);
   new_F_neighbour(n12,f11,pnflink);
   new_F_neighbour(n12,f22,pnflink);
   nodesToFaces(pfnlink,theElement->sons[0],NMASK,FMASK);
   nodesToFaces(pfnlink,theElement->sons[1],NMASK,FMASK);
   facesToFaces(pflink,theElement->sons[0],FMASK); 
   facesToFaces(pflink,theElement->sons[1],FMASK); 
}

void finer1(theElement,n1,n2,n3,f1,f2,f3,nodeNumber,faceNumber,
            pnode,pnodeil,pnodebl,pvert,
            pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
            pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, **pface, **pfaceil, **pfacebl;
NODE *n1, *n2, *n3, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT *nodeNumber, *faceNumber;
{
   NODE *n11, *n22, *n33;
   FACE *f11, *f22, *f33;

   if (n1->son == NULL) make_nson(n1,pnode,pnodeil,pnodebl,nodeNumber);
   if (n2->son == NULL) make_nson(n2,pnode,pnodeil,pnodebl,nodeNumber);
   if (n3->son == NULL) make_nson(n3,pnode,pnodeil,pnodebl,nodeNumber);
   SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
   treat_1face(f1,&f11,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_1face(f2,&f22,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_1face(f3,&f33,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   elem2(theElement,0,theElement->status,n11,n22,n33,f11,f22,f33,
                                                     pelemil,pelembl,pelement);
   theElement->sons[1] = theElement->sons[2] = theElement->sons[3] = NULL;
   new_F_neighbour(n11,f11,pnflink);
   new_F_neighbour(n22,f22,pnflink);
   new_F_neighbour(n33,f33,pnflink);
   nodesToFaces(pfnlink,theElement->sons[0],NMASK,FMASK);
   facesToFaces(pflink,theElement->sons[0],FMASK); 
}

void unif_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,
    plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
{
   NODE *n1, *n2, *n3, *n12, *n13, *n23, *n11, *n22, *n33;
   
   NODES_OF_ELEMENT(n1,n2,n3,theElement);
   n12 = treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   n13 = treat_edge(n1,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   n23 = treat_edge(n2,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
   finer4(theElement,n11,n22,n33,n12,n13,n23,
          theElement->f[0],theElement->f[1],theElement->f[2],faceNumber,
          pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,
                                                                     pelement);
   set_norm_vect_or(n11,n22,n12,theElement->f[2]); 
   set_norm_vect_or(n11,n33,n13,theElement->f[1]); 
   set_norm_vect_or(n22,n33,n23,theElement->f[0]); 
}

#elif ELEMENT_TYPE == CUBE

void unif_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,
    plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
{
   NODE *n1, *n2, *n3, *n4, *n12, *n23, *n34, *n41, *n11, *n22, *n33, *n44, *ni;
   INT i;
   
   NODES_OF_4ELEMENT(n1,n2,n3,n4,theElement);
   n12 = Treat_edge(n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   n23 = Treat_edge(n2,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   n34 = Treat_edge(n3,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   n41 = Treat_edge(n4,n1,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   SONS_OF_NODES4(n1,n2,n3,n4,n11,n22,n33,n44)
   (*pvert)->type = 0;
   POINT4(n11->myvertex->x,n22->myvertex->x,
          n33->myvertex->x,n44->myvertex->x,(*pvert)->x);
   *pnodeil = (*pnodeil)->succ = ni = (*pnode)++;
   nmake(ni,(NODE*)NULL,++(*nodeNumber),(NODE*)NULL,(*pvert)++);
   finer4(theElement,n11,n22,n33,n44,n12,n23,n34,n41,ni,
          theElement->f[0],theElement->f[1],theElement->f[2],theElement->f[3],
          faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
          pelemil,pelembl,pelement);
   set_norm_vect_or(n11,n22,n12,theElement->f[0]); 
   set_norm_vect_or(n22,n33,n23,theElement->f[1]); 
   set_norm_vect_or(n33,n44,n34,theElement->f[2]); 
   set_norm_vect_or(n44,n11,n41,theElement->f[3]); 
}

#endif

#if ELEMENT_TYPE == SIMPLEX

void refine_to_middle(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
                      pvert,plink,pnflink,pfnlink,pflink,
                      pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
{
   NODE *n1, *n2, *n3, *n11, *n22, *n33, *n123;
   FACE *f1, *f2, *f3, *f11, *f22, *f33, *f1i, *f2i, *f3i;
   
   NODES_OF_ELEMENT(n1,n2,n3,theElement);
   FACES_OF_ELEMENT(f1,f2,f3,theElement);
   if (n1->son == NULL) make_nson(n1,pnode,pnodeil,pnodebl,nodeNumber);
   if (n2->son == NULL) make_nson(n2,pnode,pnodeil,pnodebl,nodeNumber);
   if (n3->son == NULL) make_nson(n3,pnode,pnodeil,pnodebl,nodeNumber);
   SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
   (*pvert)->type = 0;
   POINT3(n1->myvertex->x,n2->myvertex->x,n3->myvertex->x,(*pvert)->x);
   *pnodeil = (*pnodeil)->succ = n123 = (*pnode)++;
   nmake(n123,(NODE*)NULL,++(*nodeNumber),(NODE*)NULL,(*pvert)++);
   new_neighbour(n11,n123,0,0,plink);
   new_neighbour(n22,n123,0,0,plink);
   new_neighbour(n33,n123,0,0,plink);
   new_neighbour(n123,n11,0,0,plink);
   new_neighbour(n123,n22,0,0,plink);
   new_neighbour(n123,n33,0,0,plink);
   treat_1face(f1,&f11,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_1face(f2,&f22,n11,n33,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   treat_1face(f3,&f33,n11,n22,faceNumber,pface,pfaceil,pfacebl,pnflink,plink);
   (*pfaceil)->succ = f1i = (*pface)++;
   f2i = (*pface)++;
   *pfaceil = f3i = (*pface)++;
   fmake(f1i,f2i,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   fmake(f2i,f3i,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   fmake(f3i,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   new_3Fneighbours(n123,f1i,f2i,f3i,pnflink);
   new_F_neighbour(n123,f11,pnflink);
   new_F_neighbour(n123,f22,pnflink);
   new_F_neighbour(n123,f33,pnflink);
   new_3Fneighbours(n11,f1i,f2i,f3i,pnflink);
   new_3Fneighbours(n22,f1i,f2i,f3i,pnflink);
   new_3Fneighbours(n33,f1i,f2i,f3i,pnflink);
   elem2(theElement,0,0,n22,n33,n123,f3i,f2i,f11,pelemil,pelembl,pelement);
   elem2(theElement,1,0,n11,n33,n123,f3i,f1i,f22,pelemil,pelembl,pelement);
   elem2(theElement,2,0,n11,n22,n123,f2i,f1i,f33,pelemil,pelembl,pelement);
   theElement->sons[3] = NULL;
   nodesToFaces(pfnlink,theElement->sons[0],NMASK,FMASK);
   nodesToFaces(pfnlink,theElement->sons[1],NMASK,FMASK);
   nodesToFaces(pfnlink,theElement->sons[2],NMASK,FMASK);
   facesToFaces(pflink,theElement->sons[0],FMASK);
   facesToFaces(pflink,theElement->sons[1],FMASK);
   facesToFaces(pflink,theElement->sons[2],FMASK);
   set_norm_vect_or(n11,n22,(NODE*)NULL,f3);
   set_norm_vect_or(n11,n33,(NODE*)NULL,f2);
   set_norm_vect_or(n22,n33,(NODE*)NULL,f1);
}

LINK *longest_edge(n0,n1,n2,ni,nf)  /* n0, n1, n2 are nodes of an element  */
NODE *n0, *n1, *n2, **ni, **nf;     /* n0->index < n1->index < n2->index   */
{
   FLOAT m, x, h01, h02, h12;

   h01 = edge_length_square(n0,n1);
   h02 = edge_length_square(n0,n2);
   h12 = edge_length_square(n1,n2);
   *ni = n1;
   *nf = n2;
   if (fabs(h01-h02)/(h01+h02) < 0.05 && fabs(h01-h12)/(h01+h12) < 0.05 &&
                                   (n0->father || n1->father || n2->father)){
      if (n0->father && n1->father && n2->father == NULL){
         *ni = n0;
         *nf = n1; 
      }
      else if (n0->father && n1->father == NULL && n2->father){
         *ni = n0;
         *nf = n2; 
      }
   }
   else{
      m = h12;
      if (h12 < h02){
         m = h02;
         *ni = n0;
      }
      if (m < h01){
         *ni = n0;
         *nf = n1;
      }
   }
   return(find_link(*ni,*nf));
}

INT check_length(n0,n1,n2,pl02,i,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
NODE *n0, *n1, *n2, **pnode, **pnodeil, **pnodebl;
LINK *pl02;                                /* (n0,n1) is the longest edge  */
VERTEX **pvert;                            /* (n1,n2) is divided as well   */
LINK **plink;                              /* (n0,n2) is not divided       */
INT  *i, *nodeNumber;
{
   NODE *n01;
   FLOAT h12=1.01*edge_length_square(n1,n2);

   n01 = mid_node(n0->son,n1->son);
   if (h12 < edge_length_square(n01,n1) || h12 < edge_length_square(n01,n2)){
      make_vertex(n0->son,n2->son,pl02,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
      *i = 1;
      return(0);
   }
   return(1);
}

void check_ratio2(n0,n1,n2,pl02,pl12,i,j,
                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber,maxindex,sigma2)
NODE *n0, *n1, *n2, **pnode, **pnodeil, **pnodebl;
LINK *pl02, *pl12;                         /* (n0,n1) is the longest edge  */
VERTEX **pvert;                            /* (n0,n2) is not divided       */
LINK **plink;                              /* (n1,n2) is not divided       */
INT  *i, j, *nodeNumber, maxindex;
FLOAT sigma2;
{
   NODE *n01;

   n01 = mid_node(n0->son,n1->son);
   if (j && n01->index > maxindex)
      if (sigma_squared(n0,n01,n2) > sigma2 || 
          sigma_squared(n01,n1,n2) > sigma2){
         if (n2->son == NULL) 
            make_nson(n2,pnode,pnodeil,pnodebl,nodeNumber);
         make_vertex(n0->son,n2->son,pl02,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
         make_vertex(n1->son,n2->son,pl12,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
         *i = 1;
      }
}

void check_ratio3(n0,n1,n2,pl02,i,j,
                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber,maxindex,sigma2)
NODE *n0, *n1, *n2, **pnode, **pnodeil, **pnodebl;
LINK *pl02;                                /* (n0,n1) is the longest edge  */
VERTEX **pvert;                            /* (n1,n2) is divided as well   */
LINK **plink;                              /* (n0,n2) is not divided       */
INT  *i, j, *nodeNumber, maxindex;
FLOAT sigma2;
{
   NODE *n01, *n12;

   n01 = mid_node(n0->son,n1->son);
   n12 = mid_node(n1->son,n2->son);
   if (j && (n01->index > maxindex || n12->index > maxindex))
      if (sigma_squared(n0,n01,n2) > sigma2 || 
          sigma_squared(n01,n12,n2) > sigma2){
         make_vertex(n0->son,n2->son,pl02,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
         *i = 1;
      }
}

void find_all_brothers(pel,el,n0,n1,n2)  /*  pel->status = 1  */
ELEMENT *pel, *el[4];
NODE **n0, **n1, **n2;
{
   ELEMENT *father_el;
   NODE *m[2];
   FLOAT xc0, xc1;
   INT i=0, j=0, k=0;

   father_el = pel->father;
   while (father_el->sons[1] == NULL)
      father_el = father_el->father;
   while (father_el->sons[i]){
      if (father_el->sons[i]->status == 1){
         el[j] = father_el->sons[i];
         while (el[j]->sons[0])
            el[j] = el[j]->sons[0];
         j++;
      }
      i++;
   }
   el[j] = NULL;
   if (father_el->status == 0){
      TOPNODES_OF_ELEMENT(*n0,*n1,*n2,father_el);
      *n0 = (*n0)->father;
      *n1 = (*n1)->father;
      *n2 = (*n2)->father;
   }
   else if (j == 2){
      for (i = 0; i < 3; i++)
         if (IS_IN(el[0]->n[i],el[1]->n[0],el[1]->n[1],el[1]->n[2]))
            m[k++] = el[0]->n[i];
         else
            *n0 = el[0]->n[i];
      if (!IS_IN(el[1]->n[0],el[0]->n[0],el[0]->n[1],el[0]->n[2]))
         *n1 = el[1]->n[0];
      else if (!IS_IN(el[1]->n[1],el[0]->n[0],el[0]->n[1],el[0]->n[2]))
         *n1 = el[1]->n[1];
      else if (!IS_IN(el[1]->n[2],el[0]->n[0],el[0]->n[1],el[0]->n[2]))
         *n1 = el[1]->n[2];
      else
         eprintf("Error in find_all_brothers.\n");
      xc0 = 0.5*((*n0)->myvertex->x[0] + (*n1)->myvertex->x[0]);
      xc1 = 0.5*((*n0)->myvertex->x[1] + (*n1)->myvertex->x[1]);
      if (fabs(m[0]->myvertex->x[0] - xc0) + fabs(m[0]->myvertex->x[1] - xc1) 
        > fabs(m[1]->myvertex->x[0] - xc0) + fabs(m[1]->myvertex->x[1] - xc1))
         *n2 = m[0];
      else
         *n2 = m[1];
   }
   else{
      *n0 = el[0]->n[0];
      *n1 = el[1]->n[0];
      *n2 = el[2]->n[0];
      if (*n0 == *n2)
         *n0 = el[0]->n[1];
   }
}

FACE *opposite_face(pel,n)
ELEMENT *pel;
NODE *n;
{
   if (n == pel->n[0])
      return(pel->f[0]);
   else if (n == pel->n[1])
      return(pel->f[1]);
   else if (n == pel->n[2])
      return(pel->f[2]);
   else{
      eprintf("Error in opposite_face.\n");
      return(NULL);
   }
}

void nodes_and_faces_of_two_green_elems(el,n1,n2,n3,n11,n22,n33,n12,
                                                              f1o,f2o,f31o,f32o)
ELEMENT *el[4];
NODE *n1, *n2, *n3, **n11, **n22, **n33, **n12;
FACE **f1o, **f2o, **f31o, **f32o;
{
   if (IS_IN(n1,el[0]->n[0],el[0]->n[1],el[0]->n[2]) &&
       IS_IN(n1,el[1]->n[0],el[1]->n[1],el[1]->n[2])){
      *n33 = n1;
      n1 = n3;
   }
   else if (IS_IN(n2,el[0]->n[0],el[0]->n[1],el[0]->n[2]) &&
            IS_IN(n2,el[1]->n[0],el[1]->n[1],el[1]->n[2])){
      *n33 = n2;
      n2 = n3;
   }
   else
      *n33 = n3;
   if (IS_IN(n1,el[0]->n[0],el[0]->n[1],el[0]->n[2])){
      *n11 = n1;
      *n22 = n2;
   }
   else{
      *n11 = n2;
      *n22 = n1;
   }
   if (*n11 != el[0]->n[0] && *n33 != el[0]->n[0])
      *n12 = el[0]->n[0];
   else if (*n11 != el[0]->n[1] && *n33 != el[0]->n[1])
      *n12 = el[0]->n[1];
   else
      *n12 = el[0]->n[2];
   *f1o  = opposite_face(el[1],*n12);
   *f2o  = opposite_face(el[0],*n12);
   *f31o = opposite_face(el[0],*n33);
   *f32o = opposite_face(el[1],*n33);
   *n11 = (*n11)->son;
   *n22 = (*n22)->son;
   *n33 = (*n33)->son;
   *n12 = (*n12)->son;
}

void nodes_and_faces_of_three_green_elems(el,n1,n2,n3,n11,n22,n33,n12,n23,
                                                   f12o,f13o,f2o,f31o,f32o,f22o)
ELEMENT *el[4];
NODE *n1, *n2, *n3, **n11, **n22, **n33, **n12, **n23;
FACE **f12o, **f13o, **f2o, **f31o, **f32o, **f22o;
{
   if (IS_IN(n1,el[1]->n[0],el[1]->n[1],el[1]->n[2]))
      *n22 = n1;
   else if (IS_IN(n2,el[1]->n[0],el[1]->n[1],el[1]->n[2]))
      *n22 = n2;
   else
      *n22 = n3;
   if (IS_IN(n1,el[2]->n[0],el[2]->n[1],el[2]->n[2]))
      *n33 = n1;
   else if (IS_IN(n2,el[2]->n[0],el[2]->n[1],el[2]->n[2]))
      *n33 = n2;
   else
      *n33 = n3;
   if (n1 != *n22 && n1 != *n33)
      *n11 = n1;
   else if (n2 != *n22 && n2 != *n33)
      *n11 = n2;
   else
      *n11 = n3;
   if (*n22 == el[1]->n[0]){
      n1 = el[1]->n[1];
      n2 = el[1]->n[2];
      *f22o = el[1]->f[0];
   }
   else if (*n22 == el[1]->n[1]){
      n1 = el[1]->n[0];
      n2 = el[1]->n[2];
      *f22o = el[1]->f[1];
   }
   else{
      n1 = el[1]->n[0];
      n2 = el[1]->n[1];
      *f22o = el[1]->f[2];
   }
   if (IS_IN(n1,el[0]->n[0],el[0]->n[1],el[0]->n[2])){
      *n12 = n1;
      *n23 = n2;
   }
   else{
      *n12 = n2;
      *n23 = n1;
   }
   *f2o  = opposite_face(el[0],*n12);
   *f31o = opposite_face(el[0],*n33);
   *f32o = opposite_face(el[1],*n23);
   *f12o = opposite_face(el[1],*n12);
   *f13o = opposite_face(el[2],*n12);
   *n11 = (*n11)->son;
   *n22 = (*n22)->son;
   *n33 = (*n33)->son;
   *n12 = (*n12)->son;
   *n23 = (*n23)->son;
}

void divide_edge0(n1,n2,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,i)
NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;           /*  n1->index < n2->index                          */
LINK **plink;             /*  (n1,n2) has to be an edge on the father level  */
INT  *nodeNumber, *i;     /*  n1, n2 have to possess sons                    */
{
   LINK *pl12=find_link(n1,n2);

   if (!IS_DIVIDED_EDGE(pl12)){
      make_vertex(n1->son,n2->son,pl12,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
      *i = 1;
   }
}

INT divide_element_edges(pel,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
ELEMENT *pel;
NODE **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK **plink;
INT  *nodeNumber;
{
   ELEMENT *el[4];
   NODE *m1, *m2, *m3, *n0, *n1, *n2, *no0, *no1, *no2, 
        *n11, *n22, *n33, *n12, *n23;
   FACE *f;
   INT i=0, i0=0, i1=0, i2=0;

   NODES_OF_ELEMENT(n0,n1,n2,pel);
   if (pel->status == 1){
      find_all_brothers(pel,el,&no0,&no1,&no2);
      if (IS_IN(n0,no0,no1,no2))
         i0 = 1;
      if (IS_IN(n1,no0,no1,no2))
         i1 = 1;
      if (IS_IN(n2,no0,no1,no2))
         i2 = 1;
   }
   if (pel->status == 0 || (i0 && i1))
      divide_edge0(n0,n1,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
   if (pel->status == 0 || (i0 && i2))
      divide_edge0(n0,n2,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
   if (pel->status == 0 || (i1 && i2))
      divide_edge0(n1,n2,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
   if (pel->status == 1 && pel->eflag == 2){
      find_all_brothers(pel,el,&m1,&m2,&m3);
      if (el[2] == NULL){
         nodes_and_faces_of_two_green_elems(el,m1,m2,m3,&n11,&n22,&n33,&n12,
                                            &f,&f,&f,&f);
         if (IS_IN(n11->father,n0,n1,n2))
            divide_edge0(n11->father,n12->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
         else
            divide_edge0(n22->father,n12->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
      }
      else{
         nodes_and_faces_of_three_green_elems(el,m1,m2,m3,&n11,&n22,&n33,
                                              &n12,&n23,&f,&f,&f,&f,&f,&f);
         if (IS_IN(n11->father,n0,n1,n2))
            divide_edge0(n11->father,n12->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
         else if (IS_IN(n22->father,n0,n1,n2)){
            divide_edge0(n22->father,n12->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
            divide_edge0(n22->father,n23->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
         }
         else
            divide_edge0(n33->father,n23->father,
                               pvert,pnode,pnodeil,pnodebl,plink,nodeNumber,&i);
      }
   }
   return(i);
}

void divide_edges_for_red_green(mg,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
MULTIGRID *mg;
NODE **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK **plink;
INT  *nodeNumber;
{
   GRID *theGrid=TOP_GRID(mg);
   ELEMENT *pel, *el[4];
   NODE *n0, *n1, *n2, *pn;
   INT i, j, k;

   for (pn = FIRSTN(theGrid); pn; pn = pn->succ)
      make_nson(pn,pnode,pnodeil,pnodebl,nodeNumber);
   for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
      if (pel->eflag == 2)
         divide_element_edges(pel,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
   do{
      k = 0;
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if (pel->status == 1 && pel->eflag != 3){
            NODES_OF_ELEMENT(n0,n1,n2,pel);
            if (IS_DIVIDED_EDGE(find_link(n0,n1)) ||
                IS_DIVIDED_EDGE(find_link(n0,n2)) ||
                IS_DIVIDED_EDGE(find_link(n1,n2)) || pel->eflag == 2){
               find_all_brothers(pel,el,&n0,&n1,&n2);
               i = 0;
               while (el[i]){
                  el[i]->eflag = 3;
                  j = divide_element_edges(el[i++],
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
                  if (j) k = 1;
               }
            }
         }
   }
   while (k);
}

void red_green_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
                      pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,
                      pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
{
   ELEMENT *el[4];
   NODE *n1, *n2, *n3, *n12, *n13, *n23, 
        *n11, *n22, *n33, *n112, *n122, *n223, *n233;
   FACE *f1, *f2, *f3, *f1o, *f12o, *f13o, *f2o, *f31o, *f32o, *f22o;
   LINK *pl12, *pl13, *pl23;
   
   if (!theElement->sons[0]){
      if (theElement->status == 0){
         NODES_OF_ELEMENT(n1,n2,n3,theElement);
         FACES_OF_ELEMENT(f1,f2,f3,theElement);
         SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
         n12 = mid_node_or_null(n11,n22);
         n13 = mid_node_or_null(n11,n33);
         n23 = mid_node_or_null(n22,n33);
         if (n12 && n13 && n23)
            finer4(theElement,n11,n22,n33,n12,n13,n23,f1,f2,f3,faceNumber,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,
                                                                      pelement);
         else if (n12 && n13)
            finer3(theElement,n33,n11,n22,n13,n12,f3,f1,f2,1,faceNumber,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else if (n12 && n23)
            finer3(theElement,n11,n22,n33,n12,n23,f1,f2,f3,1,faceNumber,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else if (n13 && n23)
            finer3(theElement,n22,n33,n11,n23,n13,f2,f3,f1,1,faceNumber,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else if (n12)
            finer2(theElement,n11,n22,n3,n12,f1,f2,f3,1,nodeNumber,faceNumber,
                   pnode,pnodeil,pnodebl,pvert,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else if (n13)
            finer2(theElement,n33,n11,n2,n13,f3,f1,f2,1,nodeNumber,faceNumber,
                   pnode,pnodeil,pnodebl,pvert,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else if (n23)
            finer2(theElement,n22,n33,n1,n23,f2,f3,f1,1,nodeNumber,faceNumber,
                   pnode,pnodeil,pnodebl,pvert,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         else
            finer1(theElement,n1,n2,n3,f1,f2,f3,nodeNumber,faceNumber,
                   pnode,pnodeil,pnodebl,pvert,
                   pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                   pelemil,pelembl,pelement);
         set_norm_vect_or(n11,n22,n12,f3); 
         set_norm_vect_or(n11,n33,n13,f2); 
         set_norm_vect_or(n22,n33,n23,f1); 
      }
      else if (theElement->eflag == 0){
         NODES_OF_ELEMENT(n1,n2,n3,theElement);
         FACES_OF_ELEMENT(f1,f2,f3,theElement);
         SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
         finer1(theElement,n1,n2,n3,f1,f2,f3,nodeNumber,faceNumber,
                pnode,pnodeil,pnodebl,pvert,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
         set_norm_vect_or(n11,n22,NULL,f3); 
         set_norm_vect_or(n11,n33,NULL,f2); 
         set_norm_vect_or(n22,n33,NULL,f1); 
      }
      else{
         find_all_brothers(theElement,el,&n1,&n2,&n3);
         if (el[2] == NULL){
            nodes_and_faces_of_two_green_elems(el,n1,n2,n3,&n11,&n22,&n33,&n12,
                                               &f1o,&f2o,&f31o,&f32o);
            n13 = mid_node(n11,n33);
            n23 = mid_node(n22,n33);
            n112 = mid_node_or_null(n11,n12);
            n122 = mid_node_or_null(n12,n22);
            finer4_green2(el[0],el[1],n11,n22,n33,n12,n13,n23,n112,n122,
                   f1o,f2o,f31o,f32o,faceNumber,pface,pfaceil,pfacebl,plink,
                   pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
            set_norm_vect_or(n22,n33,n23,f1o); 
         }
         else{
            nodes_and_faces_of_three_green_elems(el,n1,n2,n3,
                   &n11,&n22,&n33,&n12,&n23,&f12o,&f13o,&f2o,&f31o,&f32o,&f22o);
            n13 = mid_node(n11,n33);
            n112 = mid_node_or_null(n11,n12);
            n122 = mid_node_or_null(n12,n22);
            n223 = mid_node_or_null(n22,n23);
            n233 = mid_node_or_null(n23,n33);
            finer4_green3(el[0],el[1],el[2],n11,n22,n33,n12,n13,n23,n112,n122,
               n223,n233,f12o,f13o,f2o,f31o,f32o,f22o,faceNumber,pface,pfaceil,
               pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
            set_norm_vect_or(n22,n23,n223,f12o); 
            set_norm_vect_or(n33,n23,n233,f13o); 
            set_norm_vect_or(n12,n23,NULL,f22o); 
         }
         set_norm_vect_or(n11,n12,n112,f31o); 
         set_norm_vect_or(n22,n12,n122,f32o); 
         set_norm_vect_or(n11,n33,n13,f2o); 
      }
   }
}

void divide_edges_for_bisection(mg,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
MULTIGRID *mg;
NODE **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;
LINK **plink;
INT  *nodeNumber;
{
   GRID *theGrid;
   ELEMENT *pelem;
   NODE *n0, *n1, *n2, *ni, *nf;
   LINK *pl01, *pl12, *pl02, *pl_l;
   INT i=1, j=0, k=0, maxi, maxindex=0;
   FLOAT sigma2 = MAX_SIGMA*MAX_SIGMA;

   for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->eflag == 2){
            NODES_OF_ELEMENT(n0,n1,n2,pelem);
            pl_l = longest_edge(n0,n1,n2,&ni,&nf);
            if (!IS_DIVIDED_EDGE(pl_l)){
               if (ni->son == NULL) 
                  make_nson(ni,pnode,pnodeil,pnodebl,nodeNumber);
               if (nf->son == NULL) 
                  make_nson(nf,pnode,pnodeil,pnodebl,nodeNumber);
               make_vertex(ni->son,nf->son,pl_l,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
            }
         }
   while(i){
      i = 0;
      for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
         for (pelem = FIRSTELEMENT(theGrid); pelem; pelem = pelem->succ)
            if (IS_TOP_ELEMENT(pelem)){
               NODES_OF_ELEMENT(n0,n1,n2,pelem);
               pl01 = find_link(n0,n1);
               pl02 = find_link(n0,n2);
               pl12 = find_link(n1,n2);
               pl_l = longest_edge(n0,n1,n2,&ni,&nf);
               if (!(   ( IS_DIVIDED_EDGE(pl01) &&  IS_DIVIDED_EDGE(pl02) 
                                                &&  IS_DIVIDED_EDGE(pl12))
                     || (!IS_DIVIDED_EDGE(pl01) && !IS_DIVIDED_EDGE(pl02) &&
                                                   !IS_DIVIDED_EDGE(pl12)) )){
                  if ((!IS_DIVIDED_EDGE(pl01) || pl01 == pl_l) &&
                      (!IS_DIVIDED_EDGE(pl02) || pl02 == pl_l) &&
                      (!IS_DIVIDED_EDGE(pl12) || pl12 == pl_l)){
                     if (IS_DIVIDED_EDGE(pl01))
                        check_ratio2(n0,n1,n2,pl02,pl12,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                     else if (IS_DIVIDED_EDGE(pl02))
                        check_ratio2(n2,n0,n1,pl12,pl01,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                     else
                        check_ratio2(n1,n2,n0,pl01,pl02,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                  }
                  else{
                     if (!IS_DIVIDED_EDGE(pl_l)){
                        if (ni->son == NULL) 
                           make_nson(ni,pnode,pnodeil,pnodebl,nodeNumber);
                        if (nf->son == NULL) 
                           make_nson(nf,pnode,pnodeil,pnodebl,nodeNumber);
                        make_vertex(ni->son,nf->son,pl_l,
                                  pvert,pnode,pnodeil,pnodebl,plink,nodeNumber);
                        i = 1;
                     }
                     if (!(IS_DIVIDED_EDGE(pl01) && IS_DIVIDED_EDGE(pl02) 
                                                 && IS_DIVIDED_EDGE(pl12))){
                        if (pl_l == pl01){
                           if (IS_DIVIDED_EDGE(pl02)){
                              if (check_length(n1,n0,n2,pl12,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n1,n0,n2,pl12,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                           else{
                              if (check_length(n0,n1,n2,pl02,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n0,n1,n2,pl02,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                        }
                        else if (pl_l == pl02){
                           if (IS_DIVIDED_EDGE(pl01)){
                              if (check_length(n2,n0,n1,pl12,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n2,n0,n1,pl12,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                           else{
                              if (check_length(n0,n2,n1,pl01,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n0,n2,n1,pl01,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                        }
                        else{
                           if (IS_DIVIDED_EDGE(pl01)){
                              if (check_length(n2,n1,n0,pl02,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n2,n1,n0,pl02,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                           else{
                              if (check_length(n1,n2,n0,pl01,&i,
                                  pnode,pnodeil,pnodebl,pvert,plink,nodeNumber))
                                check_ratio3(n1,n2,n0,pl01,&i,j,pnode,pnodeil,
                                pnodebl,pvert,plink,nodeNumber,maxindex,sigma2);
                           }
                        }
                     }
                  }
               }
            }
      k++;
      if (j == 1 && maxindex == 0)
         maxindex = maxi;
      if (i == 0 && j == 0){
         i = j = 1;
         k--;
         maxi = *nodeNumber;
      }
   }
   if (k > 1)
      printf("%i cycles in refinement\n",k);
   else
      printf("%i cycle in refinement\n",k);
}

void bisection_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
                     pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,
                     pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
{
   NODE *n1, *n2, *n3, *ni, *nf, *n12, *n13, *n23, *n11, *n22, *n33;
   FACE *f1, *f2, *f3;
   LINK *pl12, *pl13, *pl23, *pl_l;
   
   NODES_OF_ELEMENT(n1,n2,n3,theElement);
   FACES_OF_ELEMENT(f1,f2,f3,theElement);
   SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
   pl12 = find_link(n1,n2);
   pl13 = find_link(n1,n3);
   pl23 = find_link(n2,n3);
   pl_l = longest_edge(n1,n2,n3,&ni,&nf);
   if (IS_DIVIDED_EDGE(pl12)) n12 = mid_node(n11,n22);
   if (IS_DIVIDED_EDGE(pl13)) n13 = mid_node(n11,n33);
   if (IS_DIVIDED_EDGE(pl23)) n23 = mid_node(n22,n33);
   if (IS_DIVIDED_EDGE(pl12) && IS_DIVIDED_EDGE(pl13) && IS_DIVIDED_EDGE(pl23))
      finer4(theElement,n11,n22,n33,n12,n13,n23,f1,f2,f3,faceNumber,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,
                                                                      pelement);
   else if (IS_DIVIDED_EDGE(pl12) && IS_DIVIDED_EDGE(pl13)){
      if (pl12 == pl_l)
         finer3(theElement,n22,n11,n33,n12,n13,f2,f1,f3,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
      else
         finer3(theElement,n33,n11,n22,n13,n12,f3,f1,f2,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
   }
   else if (IS_DIVIDED_EDGE(pl12) && IS_DIVIDED_EDGE(pl23)){
      if (pl12 == pl_l)
         finer3(theElement,n11,n22,n33,n12,n23,f1,f2,f3,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
      else
         finer3(theElement,n33,n22,n11,n23,n12,f3,f2,f1,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
   }
   else if (IS_DIVIDED_EDGE(pl13) && IS_DIVIDED_EDGE(pl23)){
      if (pl13 == pl_l)
         finer3(theElement,n11,n33,n22,n13,n23,f1,f3,f2,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
      else
         finer3(theElement,n22,n33,n11,n23,n13,f2,f3,f1,0,faceNumber,
                pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                pelemil,pelembl,pelement);
   }
   else if (IS_DIVIDED_EDGE(pl12))
      finer2(theElement,n11,n22,n3,n12,f1,f2,f3,0,nodeNumber,faceNumber,
             pnode,pnodeil,pnodebl,pvert,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
             pelemil,pelembl,pelement);
   else if (IS_DIVIDED_EDGE(pl13))
      finer2(theElement,n33,n11,n2,n13,f3,f1,f2,0,nodeNumber,faceNumber,
             pnode,pnodeil,pnodebl,pvert,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
             pelemil,pelembl,pelement);
   else if (IS_DIVIDED_EDGE(pl23))
      finer2(theElement,n22,n33,n1,n23,f2,f3,f1,0,nodeNumber,faceNumber,
             pnode,pnodeil,pnodebl,pvert,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
             pelemil,pelembl,pelement);
   else
      finer1(theElement,n1,n2,n3,f1,f2,f3,nodeNumber,faceNumber,
             pnode,pnodeil,pnodebl,pvert,
             pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
             pelemil,pelembl,pelement);
   SONS_OF_NODES3(n1,n2,n3,n11,n22,n33)
   set_norm_vect_or(n11,n22,n12,f3); 
   set_norm_vect_or(n11,n33,n13,f2); 
   set_norm_vect_or(n22,n33,n23,f1); 
}

#else  /*  ELEMENT_TYPE != SIMPLEX  */

void refine_to_middle(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber; 
{  eprintf("Error: refine_to_middle not available.\n");  }

void divide_edges_for_red_green(mg,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
MULTIGRID *mg; NODE **pnode, **pnodeil, **pnodebl; VERTEX **pvert; LINK **plink; INT  *nodeNumber;
{  eprintf("Error: divide_edges_for_red_green not available.\n");  }

void red_green_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber; 
{  eprintf("Error: red_green_refine not available.\n");  }
 
void divide_edges_for_bisection(mg,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
MULTIGRID *mg; NODE **pnode, **pnodeil, **pnodebl; VERTEX **pvert; LINK **plink; INT  *nodeNumber;
{  eprintf("Error: divide_edges_for_bisection not available.\n");  }

void bisection_refine(theElement,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
ELEMENT *theElement, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber; 
{  eprintf("Error: bisection_refine not available.\n");  }

#endif

#endif /*  DIM == 2  */

#if DATA_S & N_LINK_TO_NODES

void add_start_to_tstart(pgrid)
GRID *pgrid;
{
   NODE *pn;
   LINK *pli;

   for (pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
      if (pn->tstart){
         for (pli = pn->tstart; pli->next != NULL; pli=pli->next);
         pli->next = pn->start;
      } 
      else
         pn->tstart = pn->start;
}

#else

void add_start_to_tstart(pgrid)
GRID *pgrid;
{}

#endif

#if DATA_S & N_LINK_TO_FACES

void add_nfstart_to_tnfstart(pgrid)
GRID *pgrid;
{
   NODE *pn;
   NFLINK *pnfl;

   for (pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
      if (pn->tnfstart){
         for (pnfl = pn->tnfstart; pnfl->next != NULL; pnfl=pnfl->next);
         pnfl->next = pn->nfstart;
      } 
      else
         pn->tnfstart = pn->nfstart;
}

#else

void add_nfstart_to_tnfstart(pgrid)
GRID *pgrid;
{}

#endif

#if DIM == 3

INT seek_face(pel,f1,nn)
ELEMENT *pel;
FACE *f1;
FLOAT nn[DIM];
{
   if(f1 == pel->f[0]){
      normal_vector(pel->n[1]->myvertex->x,pel->n[2]->myvertex->x,
                    pel->n[3]->myvertex->x,pel->f[0],nn);
      return(0);
   }
   else if(f1 == pel->f[1]){
      normal_vector(pel->n[0]->myvertex->x,pel->n[2]->myvertex->x,
                    pel->n[3]->myvertex->x,pel->f[1],nn);
      return(0);
   }
   else if(f1 == pel->f[2]){
      normal_vector(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                    pel->n[3]->myvertex->x,pel->f[2],nn);
      return(0);
   }
   else if(f1 == pel->f[3]){
      normal_vector(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                    pel->n[2]->myvertex->x,pel->f[3],nn);
      return(0);
   }
   else 
      return(1);
}
 
#else

INT seek_face(pel,f1,nn)
ELEMENT *pel;
FACE *f1;
FLOAT nn[DIM];
{
   if(f1 == pel->f[0]){
      normal_vector(pel->n[1]->myvertex->x,pel->n[2]->myvertex->x,pel->f[0],nn);
      return(0);
   }
   else if(f1 == pel->f[1]){
      normal_vector(pel->n[0]->myvertex->x,pel->n[2]->myvertex->x,pel->f[1],nn);
      return(0);
   }
   else if(f1 == pel->f[2]){
      normal_vector(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,pel->f[2],nn);
      return(0);
   }
   else 
      return(1);
}
 
#endif

void check_normal_vector(mg,pel,f1,nn1)
MULTIGRID *mg;
ELEMENT *pel;
FACE *f1;
FLOAT nn1[DIM];
{
   GRID *pg;
   ELEMENT *el[4], *father_el;
   FLOAT nn2[DIM];
   INT i=0, j, k;
   
   if (f1->father){ 
      pg = TOP_GRID(mg);
      while(pel->n[0]->index < pg->minNodeIndex) pg = pg->coarser;
      while(f1->index > pg->maxFaceIndex) f1 = f1->father;
      if (k=seek_face(pel,f1,nn2) && pel->status == 1){
         while(pel->father->sons[i] && 
                                (k=seek_face(pel->father->sons[i],f1,nn2))) i++;
         if (k){
            father_el = pel->father;
            while (father_el->sons[1] == NULL)
               father_el = father_el->father;
            i = j = 0;
            while (father_el->sons[i]){
               if (father_el->sons[i]->status == 1){
                  el[j] = father_el->sons[i];
                  while (el[j]->sons[0]->sons[0])
                     el[j] = el[j]->sons[0];
                  j++;
               }
               i++;
            }
            el[j] = NULL;
            i = 0;
            while (el[i] && (k=seek_face(el[i],f1,nn2))) i++;
         }
      }
      if (k) eprintf("Can't find the correct face.\n");
      else if (DOT(nn1,nn2) < 0.9)
        eprintf("Wrong normal vector!\n");
   }
}
 
#if DIM == 3

void test_or(mg,pel,f1,n2,n3,n4)
MULTIGRID *mg;
ELEMENT *pel;
FACE *f1;
NODE *n2, *n3, *n4;
{
   FLOAT nn1[DIM];
          
   normal_vector(n2->myvertex->x,n3->myvertex->x,n4->myvertex->x,f1,nn1);
   pel = pel->father;
   if (pel->status < 0)
      if(f1 == pel->f[0] || f1 == pel->f[1] || f1 == pel->f[2] || f1 == pel->f[3])
         check_normal_vector(mg,pel->father,f1,nn1);  
      else if (f1->father->index < 0){
         check_normal_vector(mg,pel,f1->father,nn1);      
         check_normal_vector(mg,pel->father,f1->father->father,nn1);
      }
      else
         check_normal_vector(mg,pel->father,f1,nn1);
   else 
      check_normal_vector(mg,pel,f1,nn1);
}
       
void test_orientation(mg,pel)
MULTIGRID *mg;
ELEMENT *pel;
{   
   test_or(mg,pel,pel->f[0],pel->n[1],pel->n[2],pel->n[3]);
   test_or(mg,pel,pel->f[1],pel->n[0],pel->n[2],pel->n[3]);
   test_or(mg,pel,pel->f[2],pel->n[0],pel->n[1],pel->n[3]);
   test_or(mg,pel,pel->f[3],pel->n[0],pel->n[1],pel->n[2]);       
}

#elif ELEMENT_TYPE == SIMPLEX

void test_or(mg,pel,f1,n2,n3)
MULTIGRID *mg;
ELEMENT *pel;
FACE *f1;
NODE *n2, *n3;
{
   FLOAT nn1[DIM];
          
   normal_vector(n2->myvertex->x,n3->myvertex->x,f1,nn1);
   pel = pel->father;
   if (pel->status < 0)
      if(f1 == pel->f[0] || f1 == pel->f[1] || f1 == pel->f[2])
         check_normal_vector(mg,pel->father,f1,nn1);  
      else if (f1->father->index < 0){
         check_normal_vector(mg,pel,f1->father,nn1);      
         check_normal_vector(mg,pel->father,f1->father->father,nn1);
      }
      else
         check_normal_vector(mg,pel->father,f1,nn1);
   else 
      check_normal_vector(mg,pel,f1,nn1);
}
       
void test_orientation(mg,pel)
MULTIGRID *mg;
ELEMENT *pel;
{   
   test_or(mg,pel,pel->f[0],pel->n[1],pel->n[2]);
   test_or(mg,pel,pel->f[1],pel->n[0],pel->n[2]);
   test_or(mg,pel,pel->f[2],pel->n[0],pel->n[1]);
}

#else

void test_orientation(mg,pel)
MULTIGRID *mg; ELEMENT *pel;
{}

#endif

 void check_orientation(mg)
 MULTIGRID *mg;
 {
    ELEMENT *pel;
    NODE *pn;
    LINK *pli, *pli2;
    INT j=0;
    
    for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel != NULL; pel = pel->succ)
       test_orientation(mg,pel);
    printf("Orientation checked.\n\n"); 
 }
 
#if E_DATA & E_E_NEIGHBOURS

void element_neighbours(tGrid,pelink)
GRID *tGrid;
ELINK **pelink;
{
   ELEMENT *pelem, *pel;
   ELINK *peli;
   INT i;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      for (i = 0; i < SONS; i++){
         pel = pelem->father->sons[i];
         if (pel && pelem != pel && ARE_NEIGHBOURS(pelem,pel))
            new_elem_neighbour(pelem,pel,pelink);
      }
      for (peli = pelem->father->estart; peli != NULL; peli = peli->next)
         for (i = 0; i < SONS; i++){
            pel = peli->nbel->sons[i];
            if (pel && ARE_NEIGHBOURS(pelem,pel))
               new_elem_neighbour(pelem,pel,pelink);
         }
   }
}

#else

void element_neighbours(tGrid,pelink)
GRID *tGrid;
ELINK **pelink;
{}

#endif

#if E_DATA & E_E_FNEIGHBOURS

void element_fneighbours(tGrid)
GRID *tGrid;
{
   ELEMENT *pelem, *pel;
   INT i, j;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      for (i = 0; i < SONS; i++){
         pel = pelem->father->sons[i];
         if (pel && pelem != pel) 
            add_elem_fneighbour(pelem,pel);
      }
      for (j = 0; j < SIDES; j++)
         if(pel=NB_EL(pelem->father,j))
            for (i = 0; i < SONS; i++)
               if (pel->sons[i])
                  add_elem_fneighbour(pelem,pel->sons[i]);
      for (i = 0; i < SIDES; i++)
         if (IS_BF(pelem->f[i]))
            NB_EL(pelem,i) = NULL;
   }
}

#else

void element_fneighbours(tGrid)
GRID *tGrid;
{}

#endif

#if DATA_S & SPECIAL_NODES_AND_FACES

void smidpoint_type(f,n1,n2,mask)  /*  n1, n2 are end points of f  */
FACE *f;                           /*  n1->index < n2->index       */
NODE *n1, *n2;
INT mask;
{
   NODE *n12;
   LINK *pl12;

   if (f->s_face){
      for (pl12 = n1->tstart; pl12->nbnode != n2; pl12 = pl12->next);
      if (IS_DIVIDED_EDGE(pl12)){
         n12 = mid_node(n1->son,n2->son);
         NTYPE(n12) |= mask;
      }
   }
}

void change_type_of_smidpoints(theGrid,mask)
GRID *theGrid;
INT mask;
{
   ELEMENT *pelem;
   NODE *n0, *n1, *n2;
   FACE *f0, *f1, *f2;

   for (pelem = FIRSTELEMENT(theGrid); pelem; pelem = pelem->succ){
      NODES_OF_ELEMENT(n0,n1,n2,pelem);
      FACES_OF_ELEMENT(f0,f1,f2,pelem);
      smidpoint_type(f0,n1,n2,mask);
      smidpoint_type(f1,n0,n2,mask);
      smidpoint_type(f2,n0,n1,mask);
   }
}

#else

void change_type_of_smidpoints(theGrid,mask)
GRID *theGrid; INT mask;
{}

#endif

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void set_itypes_on_refined_level(tGrid)
GRID *tGrid;
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   INT i;
   
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      ITYPE(pel) = ITYPE(pel->father);
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      ITYPE(pnode) = 0;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         ITYPE(pel->n[i]) |= ITYPE(pel);
   for (pface = FIRSTF(tGrid); pface; pface = pface->succ){
      FTYPE(pface) &= ~IBNDRY_FACE_BIT;
      if (pface->father)
         FTYPE(pface) |= FTYPE(pface->father) & IBNDRY_FACE_BIT;
   }
}

#else

void set_itypes_on_refined_level(tGrid)
GRID *tGrid;
{}

#endif

void refinement(mg,pnode,pface,pvert,plink,pelink,pnelink,
                pnflink,pfnlink,pflink,psnode,psface,pp_node,pp_face,pelement,
                pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata,ref_type)
MULTIGRID *mg;
ELEMENT **pelement;
NODE **pnode;
FACE **pface;
VERTEX **pvert;
LINK **plink;
ELINK **pelink;
NELINK **pnelink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
SNODE **psnode;
SFACE **psface;
P_NODE **pp_node;
P_FACE **pp_face;
NLGLINK **pnlglink;
LGNLINK **plgnlink;
LGFLINK **plgflink;
LGLGLINK **plglglink;
FLGLINK **pflglink;
LGDATA **plgdata;
INT ref_type;
{
   NODE lin, ldbn;
   FACE lif, ldbf, *pf;
   ELEMENT lie, ldbe, *theElement, *pelemil, *pelembl;
   GRID *pgrid, *theGrid;
   INT nodeNumber, faceNumber, elemNumber;
   
   lin.succ = NULL;
   lif.succ = NULL;
   lie.succ = NULL;
   ldbn.succ = NULL;
   ldbf.succ = NULL;
   ldbe.succ = NULL;    
   pgrid = TOP_GRID(mg) + 1;
   pgrid->first_grid = TOP_GRID(mg)->first_grid;
   pgrid->coarser = TOP_GRID(mg);
   pgrid->finer = NULL;
   pgrid->level = mg->toplevel + 1;
   pgrid->lastNode = &lin;
   pgrid->ldbn = &ldbn;
   pgrid->lastFace = &lif;
   pgrid->ldbf = &ldbf;
   pelemil = &lie;
   pelembl = &ldbe;
   nodeNumber = mg->nodeNumber;
   faceNumber = mg->faceNumber;
   prepare_links(mg);
   if (ref_type == BISECTION && DIM ==2)
      divide_edges_for_bisection(mg,pnode,&pgrid->lastNode,&pgrid->ldbn,
                                                       pvert,plink,&nodeNumber);
   else if (ref_type == RED_GREEN && DIM ==2)
      divide_edges_for_red_green(mg,pnode,&pgrid->lastNode,&pgrid->ldbn,
                                                       pvert,plink,&nodeNumber);
   for (theElement = FIRSTELEMENT(TOP_GRID(mg)); theElement != NULL; 
                                                 theElement = theElement->succ)
      if (ref_type == UNIFORM)
         unif_refine(theElement,pnode,&pgrid->lastNode,&pgrid->ldbn,pface,
              &pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,pflink,
              &pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
      else if (ref_type == REFINE_TO_MID)
         refine_to_middle(theElement,pnode,&pgrid->lastNode,&pgrid->ldbn,pface,
              &pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,pflink,
              &pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
      else if (ref_type == BISECTION && DIM ==2)
         bisection_refine(theElement,pnode,&pgrid->lastNode,&pgrid->ldbn,pface,
              &pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,pflink,
              &pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
      else if (ref_type == RED_GREEN && DIM ==2)
         red_green_refine(theElement,pnode,&pgrid->lastNode,&pgrid->ldbn,pface,
              &pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,pflink,
              &pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
      else
         eprintf("Error: element refinement not available.\n");
   pgrid->firstNode = lin.succ;
   pgrid->fdbn = ldbn.succ;
   pgrid->firstFace = lif.succ;
   pgrid->fdbf = ldbf.succ;
   if (lie.succ){
      FIRSTELEMENT(pgrid) = lie.succ;
      pelemil->succ = FDBE(pgrid) = ldbe.succ;
   }
   else
      FIRSTELEMENT(pgrid) = FDBE(pgrid) = ldbe.succ;
   if (ldbn.succ){
      pgrid->firstN = ldbn.succ;
      pgrid->ldbn->succ = pgrid->firstNode;
   }
   else{
      pgrid->firstN = pgrid->firstNode;
      pgrid->ldbn = NULL;
   }    
   if (ldbf.succ){
      pgrid->firstF = ldbf.succ;
      pgrid->ldbf->succ = pgrid->firstFace;
   }
   else{
      pgrid->firstF = pgrid->firstFace;
      pgrid->ldbf = NULL;
   }
   pgrid->minNodeIndex = (pgrid-1)->maxNodeIndex + 1;
   pgrid->maxNodeIndex = nodeNumber;
   pgrid->minFaceIndex = (pgrid-1)->maxFaceIndex + 1;
   pgrid->maxFaceIndex = faceNumber;
   TOP_GRID(mg)->finer = pgrid;
   for (theGrid = pgrid; theGrid; theGrid = theGrid->coarser)
      theGrid->top_grid = pgrid;
   elemNumber = ((long)(*pelement) - 
      (long)(TMIN(FIRSTELEMENT(pgrid),FDBE(pgrid))))/sizeof(ELEMENT);
   mg->nodeNumber = nodeNumber;
   mg->faceNumber = faceNumber;
   mg->elemNumber += elemNumber;
   mg->toplevel++;
   TOP_GRID(mg) = pgrid;
   add_start_to_tstart(pgrid);
   add_nfstart_to_tnfstart(pgrid);
   compute_number_of_neighbouring_elements_for_nodes(pgrid);
   lg_neighbours(pgrid,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata,
                                                                         U,F,D);
   element_neighbours(pgrid,pelink);
   element_fneighbours(pgrid);
   make_previous_nodes(pgrid);
   make_previous_faces(pgrid);
   for (pf = FIRSTF(pgrid); pf != NULL; pf = pf->succ)
      SET_FTOPLEVEL(pf,mg->toplevel);
   shift_boundary_vertices(pgrid); 
   check_new_level(mg);
   check_orientation(mg);
   change_type_of_smidpoints(pgrid->coarser,SNMASK);
   make_special_nodes_and_faces(pgrid,psnode,psface,SNMASK,SFMASK);
   set_itypes_on_refined_level(pgrid);
   FIRSTPN(pgrid) = NULL;
   FIRSTPF(pgrid) = NULL;
   printf("Level %1i: node number = %i, face number = %i,\n",mg->toplevel,
                                   pgrid->maxNodeIndex-pgrid->minNodeIndex + 1,
                                   pgrid->maxFaceIndex-pgrid->minFaceIndex + 1);
   printf("         element number = %i\n",elemNumber);  
   printf("Multigrid: node number = %i, face number = %i,\n",mg->nodeNumber,
                                                             mg->faceNumber);
   printf("           element number = %i\n",mg->elemNumber);
   max_and_min_edge_in_grid(TOP_GRID(mg));
   check_ratio(TOP_GRID(mg));
   make_periodic_nodes(TOP_GRID(mg),pp_node,PNMASK,0);
} 
