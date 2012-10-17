/******************************************************************************/
/*                                                                            */
/*                           nonuniform refinement                            */
/*                                                                            */
/******************************************************************************/

INT find_level(mg,n)
MULTIGRID *mg;
NODE *n;
{
   GRID *pg;
   
   pg = FIRSTGRID(mg);
   while( n->index > pg->maxNodeIndex) pg = pg->finer;
   return(pg->level);
}
 
INT can_be_edge(mg,n1,n2)/* after applying this function, level(n1)<=level(n2)*/
MULTIGRID *mg;
NODE **n1, **n2;
{
   NODE *n;
   INT l1, l2;
   
   l1 = find_level(mg,*n1);
   l2 = find_level(mg,*n2);
   if (l1 == l2) return(1);
   else if (l1 > l2){
      n = *n1;
      *n1 = *n2;
      *n2 = n;
   }
   return(0);
}

LINK *get_link(n1,n2)
NODE *n1, *n2;
{
   LINK *pli;
   NODE *n;

   if (n1->index > n2->index){
      n = n1;
      n1 = n2;
      n2 = n;
   }
   for (pli=n1->tstart; pli && pli->nbnode!=n2; pli=pli->next);
   return(pli);
}

LINK *Get_link(n1,n2)      
NODE *n1, *n2;              
{                           
   LINK *pli;                 
  
   for (pli=n1->start; pli && pli->nbnode!=n2; pli=pli->next);
   return(pli);
}
 
NFLINK *get_nflink(n,f)
NODE *n; 
FACE *f;              
{                             
   NFLINK *pnf;                
   
   for (pnf=n->nfstart; pnf!=NULL && pnf->nbface!=f; pnf=pnf->next);
   return(pnf);
}
 
LINK *lower_edge(mg,n1,n2)  /* (n1,n2) is a non-green edge on some level */
MULTIGRID *mg;
NODE **n1, **n2;
{
   LINK *pli;
   
   if ( (*n1 = (*n1)->father) == NULL || (*n2 = (*n2)->father) == NULL ) 
      return(NULL);
   while(1)
      if (can_be_edge(mg,n1,n2)){
         pli = get_link(*n1,*n2);
         if (pli)
            return(pli);
         else if ( (*n1 = (*n1)->father) == NULL || 
                   (*n2 = (*n2)->father) == NULL )
            return(NULL);
      }
      else if ( (*n2 = (*n2)->father) == NULL ) return(NULL);
}
  
LINK *upper_edge(mg,n1,n2)  /* (n1,n2) is a non-green edge on some level */
MULTIGRID *mg;
NODE **n1, **n2;
{
   LINK *pli;
   
   if ( (*n1 = (*n1)->son) == NULL || (*n2 = (*n2)->son) == NULL || 
        INDEX(*n1) > mg->nodeNumber || INDEX(*n2) > mg->nodeNumber) 
      return(NULL);
   while(1)
      if (can_be_edge(mg,n1,n2)){
         pli = get_link(*n1,*n2);
         if (pli)
            return(pli);
         else if ( (*n1 = (*n1)->son) == NULL || (*n2 = (*n2)->son) == NULL || 
                   INDEX(*n1) > mg->nodeNumber || INDEX(*n2) > mg->nodeNumber) 
            return(NULL);  
      }
      else if ( (*n1 = (*n1)->son) == NULL || (*n1)->index > mg->nodeNumber)
         return(NULL);
}

void mark_other_edges(mg,n1,n2,r)  /* (n1,n2) is an edge on some level */
MULTIGRID *mg;
NODE *n1, *n2;
INT r;
{
   LINK *pli;
   NODE *m1, *m2;
    
   m1 = n1; 
   m2 = n2;
   while( pli=lower_edge(mg,&m1,&m2) ) MARK_TYPE_OF_EDGE(pli,r);
   while( pli=upper_edge(mg,&n1,&n2) ) MARK_TYPE_OF_EDGE(pli,r);
}
 
void mark_edge_to_divide(mg,n1,n2) 
MULTIGRID *mg;
NODE *n1, *n2;
{
   LINK *pli;
   
   pli = get_link(n1,n2);
   if (pli && !IS_EDGE_TO_DIVIDE(pli)){
      MARK_EDGE_TO_DIVIDE(pli);
      mark_other_edges(mg,n1,n2,EDGE_TO_DIVIDE);
   }
}

ELEMENT *forefather(theElement)
ELEMENT *theElement;
{
   do
      theElement = theElement->father;
   while( theElement->status == 1 );
   return(theElement);
}
 
void link_or_middle(n1,n2,pl12,n12)
NODE *n1, *n2, **n12;
LINK **pl12;
{
   if ( *pl12=get_link(n1,n2) ) *n12 = NULL;
   else *n12 = mid_node(n1,n2);
}
 
void check_edge(mg,n1,n2,pl12,i,pnode,pnodeil,pnodebl,nodeNumber) 
MULTIGRID *mg;
NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;
LINK *pl12;
INT  *nodeNumber, *i;
{
   if (!IS_EDGE_TO_DIVIDE(pl12))
      if (NO_NODE_ON_NEW_LEVEL(n1,mg)){
         if (NO_NODE_ON_NEW_LEVEL(n2,mg))
            if (n1->index < n2->index)
               make_nson(TOP_NODE(n1),pnode,pnodeil,pnodebl,nodeNumber);
            else
               make_nson(TOP_NODE(n2),pnode,pnodeil,pnodebl,nodeNumber); 
         else if (n1->index < n2->index){
            mark_edge_to_divide(mg,n1,n2);
            (*i)++;
         }
      }
      else if (NO_NODE_ON_NEW_LEVEL(n2,mg)){
         if (n1->index > n2->index){
            mark_edge_to_divide(mg,n1,n2);
            (*i)++;
         }
      }
      else if( (n1->index < n2->index && 
                TOP_NODE(n1)->index > TOP_NODE(n2)->index) ||
               (n1->index > n2->index && 
                TOP_NODE(n1)->index < TOP_NODE(n2)->index) ){
         mark_edge_to_divide(mg,n1,n2);
         (*i)++;
      } 
}       
          
/******************************************************************************/

 void edge(n1,n2,flag,type,plink) /*  flag=1 <-> n2 is a middle of some edge  */
 NODE *n1, *n2;
 LINK **plink;
 INT flag, type;
 {
    if (IS_FN(n2)){
       new_neighbour(n1,n2,flag,type,plink);
       if (IS_FN(n1))
          new_neighbour(n2,n1,0,type,plink);
    }
    else{
       new_neighbour(n2,n1,0,type,plink);
       if (NOT_FN(n1))
          new_neighbour(n1,n2,flag,type,plink);
    }
 }    
 
 void new_edge(mg,n1,n2,pl12,flag,plink)
 MULTIGRID *mg;
 NODE *n1, *n2;
 LINK *pl12, **plink;
 INT flag;
 {
    edge(TOP_NODE(n1),TOP_NODE(n2),flag,IS_BND_EDGE(pl12),plink);
    MARK_TREATED_EDGE(pl12);
    mark_other_edges(mg,n1,n2,TREATED_EDGE);
 } 
 
 NODE *divide_edge(mg,n1,n2,pl12,pvert,pnode,pnodeil,pnodebl,plink,nodeNumber)
 MULTIGRID *mg;
 NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;
 VERTEX **pvert;
 LINK *pl12, **plink;
 INT  *nodeNumber;
 {
    mark_other_edges(mg,n1,n2,DIVIDED_EDGE);
    return(make_vertex(TOP_NODE(n1),TOP_NODE(n2),pl12,pvert,pnode,pnodeil,
                                                     pnodebl,plink,nodeNumber));
 }                                         
 
 NODE *nu_treat_edge(mg,n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber) 
 MULTIGRID *mg;                           /*  n1, n2 are nodes on some top-   */
 NODE *n1, *n2, **pnode, **pnodeil, **pnodebl;/* element level and (n1,n2) is */
 VERTEX **pvert;                          /*  either an edge of a topelement  */
 LINK **plink;                            /*  or an edge of a forefather of   */
 INT  *nodeNumber;                        /*  green elements.                 */
 {
    NODE *mnode, *n12;
    LINK *pl12, *pl112, *pl122, *pli;
    
    pl12 = get_link(n1,n2); /* pl12=NULL -> (n1,n2) is not an edge on father  */
    mnode = NULL;           /*                                   father level */
    if (pl12 && IS_DIVIDED_EDGE(pl12))
       mnode = mid_node(TOP_NODE(n1),TOP_NODE(n2));
    
    if ( mnode == NULL && (pl12 == NULL || !IS_TREATED_EDGE(pl12))){
       if (NO_NODE_ON_NEW_LEVEL(n1,mg))/* if pl12!=NULL, then !IS_DIVIDED_EDGE*/
          make_nson(TOP_NODE(n1),pnode,pnodeil,pnodebl,nodeNumber);
       if (NO_NODE_ON_NEW_LEVEL(n2,mg)) 
          make_nson(TOP_NODE(n2),pnode,pnodeil,pnodebl,nodeNumber);
       if (pl12 == NULL) { /* (n1,n2) is divided already on the coarser level */
          if (NO_NODE_ON_NEW_LEVEL( (n12=mid_node(n1,n2)), mg) ){
             make_nson(TOP_NODE(n12),pnode,pnodeil,pnodebl,nodeNumber);
             mnode = TOP_NODE(n12);
             pl112 = get_link(n1,n12);
             pl122 = get_link(n12,n2);
             if (IS_EDGE_TO_DIVIDE(pl112)){
                divide_edge(mg,n1,n12,pl112,pvert,pnode,pnodeil,pnodebl,plink,
                                                                    nodeNumber);
                if (IS_EDGE_TO_DIVIDE(pl122))
                   divide_edge(mg,n12,n2,pl122,pvert,pnode,pnodeil,pnodebl,
                                                              plink,nodeNumber);
                else{
                  new_neighbour(TOP_NODE(n2),mnode,0,IS_BND_EDGE(pl112),plink);
                  if (IS_FN(n2))
                   new_neighbour(mnode,TOP_NODE(n2),0,IS_BND_EDGE(pl112),plink);
                  MARK_TREATED_EDGE(pl122);                  
                  mark_other_edges(mg,n12,n2,TREATED_EDGE);
                }
             }
             else if (IS_EDGE_TO_DIVIDE(pl122)){
                new_neighbour(TOP_NODE(n1),mnode,0,IS_BND_EDGE(pl112),plink);
                if (IS_FN(n1))
                   new_neighbour(mnode,TOP_NODE(n1),0,IS_BND_EDGE(pl112),plink);
                MARK_TREATED_EDGE(pl112);                  
                mark_other_edges(mg,n1,n12,TREATED_EDGE);
                divide_edge(mg,n12,n2,pl122,pvert,pnode,pnodeil,pnodebl,plink,
                                                                    nodeNumber);
             }
             else{      
                new_neighbour(TOP_NODE(n1),mnode,1,IS_BND_EDGE(pl112),plink);
                new_neighbour(TOP_NODE(n2),mnode,1,IS_BND_EDGE(pl112),plink);
                if (IS_FN(n12))
                   neighbours_of_middle(mnode,TOP_NODE(n1),TOP_NODE(n2),plink);
                MARK_TREATED_EDGE(pl112); /* For Dirichlet nodes, neighbour   */
                MARK_TREATED_EDGE(pl122); /* nodes with lower indices are not */
                mark_other_edges(mg,n1,n12,TREATED_EDGE);  /*         stored. */
                mark_other_edges(mg,n12,n2,TREATED_EDGE);
             }
          }
          else{
             mnode = TOP_NODE(n12);
             pl112 = get_link(n1,n12);
             pl122 = get_link(n12,n2);
             if (!IS_EDGE_TO_DIVIDE(pl112) && !IS_DIVIDED_EDGE(pl112) &&
                 !IS_EDGE_TO_DIVIDE(pl122) && !IS_DIVIDED_EDGE(pl122)){
                if (!IS_TREATED_EDGE(pl112)) 
                   new_edge(mg,n1,n12,pl112,1,plink);
                else{
                   for (pli = TOP_NODE(n1)->start; pli->nbnode != mnode; 
                                                               pli = pli->next);
                   MARK_EDGE_TO_MIDDLE(pli);
                }
                if (!IS_TREATED_EDGE(pl122)) 
                   new_edge(mg,n2,n12,pl122,1,plink);
                else{
                   for (pli = TOP_NODE(n2)->start; pli->nbnode != mnode; 
                                                               pli = pli->next);
                   MARK_EDGE_TO_MIDDLE(pli);
                }
             }
             else{
                if (!IS_EDGE_TO_DIVIDE(pl112) && !IS_DIVIDED_EDGE(pl112) &&
                                                 !IS_TREATED_EDGE(pl112))
                   new_edge(mg,n1,n12,pl112,0,plink);
                else if (IS_EDGE_TO_DIVIDE(pl112)) 
                   divide_edge(mg,n1,n12,pl112,pvert,pnode,pnodeil,pnodebl,
                                                              plink,nodeNumber);
                if (!IS_EDGE_TO_DIVIDE(pl122) && !IS_DIVIDED_EDGE(pl122) &&
                                                 !IS_TREATED_EDGE(pl122))
                   new_edge(mg,n12,n2,pl122,0,plink);
                else  if (IS_EDGE_TO_DIVIDE(pl122)) 
                   divide_edge(mg,n12,n2,pl122,pvert,pnode,pnodeil,pnodebl,
                                                              plink,nodeNumber);
             }
          }
       }
       else if (!IS_EDGE_TO_DIVIDE(pl12))
          new_edge(mg,n1,n2,pl12,0,plink);
       else 
          mnode = divide_edge(mg,n1,n2,pl12,pvert,pnode,pnodeil,pnodebl,plink,
                                                                    nodeNumber);
    }
    return(mnode);
 }

 INT max2(a,b,i,j,lm)  /*    if a>=b  -> i, lm:=a;  */
 INT a,b,i,j,*lm;      /*  else       -> j, lm:=b;  */
 {
    if (a >= b) {
       *lm = a; 
       return(i);
    }
    else{
       *lm = b;
       return(j);
    }
 }  
  
 INT max3(a,b,c,i,j,k,lm)  /*       if a>=b && a>=c -> i, lm:=a;  */
 INT a,b,c,i,j,k,*lm;      /*  else if b>=a && b>=c -> j, lm:=b;  */
 {                         /*  else                 -> k, lm:=c.  */
    if (a >= c){
       if (a >= b) {
          *lm = a; 
          return(i);
       }
       else{
          *lm = b;
          return(j);
       }
    }
    else if (b >= c) {
       *lm = b;
       return(j);
    }
    *lm = c;
    return(k);
 }  
  
 INT new_2Fneighbours(n,nf1,nf2,pnflink)
 NODE *n;
 FACE *nf1, *nf2;
 NFLINK **pnflink;
 {
    ADD_FIRST_FACE1(nf1)
    ADD_NEXT_FACE(nf2)
    ((*pnflink)++)->next = NULL;
    return(0);
 }
 
 INT new_2SFneighbours(n,nf1,nf2,pnflink)
 NODE *n;
 FACE *nf1, *nf2;
 NFLINK **pnflink;
 {
    ADD_FIRST_FACE1(nf1)
    if (IS_FF(nf2)) ADD_NEXT_FACE(nf2)
    ((*pnflink)++)->next = NULL;
    return(0);
 }
 
 INT new_3SFneighbours(n,nf1,nf2,nf3,pnflink)
 NODE *n;
 FACE *nf1, *nf2, *nf3;
 NFLINK **pnflink;
 {
    ADD_FIRST_FACE1(nf1)
    if (IS_FF(nf2)) ADD_NEXT_FACE(nf2)
    if (IS_FF(nf3)) ADD_NEXT_FACE(nf3)
    ((*pnflink)++)->next = NULL;
    return(0);
 }

 INT new_7Fneighbours(n,nf1,nf2,nf3,nf4,nf5,nf6,nf7,pnflink)
 NODE *n;
 FACE *nf1, *nf2, *nf3, *nf4, *nf5, *nf6, *nf7;
 NFLINK **pnflink;
 {
    ADD_FIRST_FACE1(nf1)
    ADD_NEXT_FACE(nf2)
    ADD_NEXT_FACE(nf3)
    if (IS_FF(nf4)) ADD_NEXT_FACE(nf4)
    if (IS_FF(nf5)) ADD_NEXT_FACE(nf5)
    if (IS_FF(nf6)) {
       ADD_NEXT_FACE(nf6)
       ADD_NEXT_FACE(nf7)
    }
    ((*pnflink)++)->next = NULL;
    return(0);
 }
 
 void k_inner_faces(k,fi,faceNumber,pface,pfaceil)
 INT k, *faceNumber;
 FACE **fi, **pface, **pfaceil;
 {
    INT i;
    
    *pfaceil = ((*pfaceil)->succ = *pface) + (k-1);
    for (i=1; i<=k; i++){
       fmake(*pface,*pface+1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
       fi[i] = (*pface)++;
    }
    ((*pface)-1)->succ = NULL;
 }
  
 void nodes_and_facesToFaces(theElement,k,pfnlink,pflink)
 ELEMENT *theElement;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT k;
 {
    INT i; 

    for (i=0; i<k; i++){
       nodesToFaces(pfnlink,theElement->sons[i],NMASK,FMASK);
       facesToFaces(pflink,theElement->sons[i],FMASK);
    }
    theElement->sons[k] = NULL;
 }
 
 void treat_fictive_face(f,pfil,pfbl,fictNumber)
 FACE *f, *pfil, *pfbl;
 INT *fictNumber;
 {
    FACE *pf;
    
    if (NOT_FF(f)) pf = pfbl;
    else pf = pfil;
    while (pf->succ != f) pf = pf->succ;
    pf->succ = f->succ;
    f->index = -f->index;
    (*fictNumber)++;
 }        
             
 void treat_fictive_element(pel,peil,pebl,fictElNumber)
 ELEMENT *pel, *peil, *pebl;
 INT *fictElNumber;
 {
    while (pebl->succ != NULL && pebl->succ != pel) pebl = pebl->succ;
    if (pebl->succ != pel){
       while (peil->succ != pel) peil = peil->succ;
       peil->succ = pel->succ;
    }
    else pebl->succ = pel->succ;
    pel->status = -1;
    (*fictElNumber)++;
 }        
 
 NODE *middle(n1,n2)  /* looks for the middle of (n1,n2) */
 NODE *n1, *n2;
 {
    LINK *pli, *plin;
       
    if (get_link(n1,n2) == NULL)
       for (pli = n1->start; ; pli = pli->next){
          for(plin = n2->start; plin != NULL && plin->nbnode != pli->nbnode; 
                                                               plin=plin->next);
          if (plin->nbnode == pli->nbnode && IS_EDGE_TO_MIDDLE(plin) 
                                && IS_EDGE_TO_MIDDLE(pli)) return(plin->nbnode);
       }
    else
       return(NULL);
 }   
 
#if (ELEMENT_TYPE == SIMPLEX) && (DIM == 2)

void mark_all_element_edges(mg,n1,n2,n3) 
MULTIGRID *mg;
NODE *n1, *n2, *n3;
{
   mark_edge_to_divide(mg,n1,n2);
   mark_edge_to_divide(mg,n1,n3);      
   mark_edge_to_divide(mg,n2,n3);                                        
} 
 
void fsons(f,fac)
FACE *f, *fac[4];
{
   fac[0] = f->sons[0];
   fac[1] = f->sons[1];
   if (fac[2] = f->sons[2])
      fac[3] = f->sons[3];
   else
      fac[3] = NULL;
}
 
void fsons_of_next_level(fac)  
FACE *fac[4];
{
   INT i = -1;
   
   if (fac[0]->sons[1])  /*  geometrical form changed  */
      fsons(fac[0],fac);
   else if (fac[1])
      while (++i < 4 && fac[i]) fac[i] = fac[i]->sons[3]; 
   else 
      fac[0] = fac[0]->sons[0];
}
    
INT correct_ind(pg,f,n1,n2)
GRID *pg;
FACE *f;
NODE *n1, *n2;
{ 
   if (pg && f->index <= pg->maxFaceIndex && f->index > 0){
      while(f->index < pg->minFaceIndex) pg = pg->coarser;
      while(n1->index > pg->maxNodeIndex) n1 = n1->father;
      while(n2->index > pg->maxNodeIndex) n2 = n2->father;
      return(norm_vect_or(n1,n2,FTYPE(f)&ORIENT));
   }
   else
      return(FTYPE(f)&ORIENT);
}
       
void data1_from_forefather(theGrid,pelem,no1,no2,no3,el)   
GRID *theGrid;
ELEMENT *pelem, *el[4];     /* pelem is a green element */
NODE **no1, **no2, **no3;
{   
   ELEMENT *pel;
   INT i=0;
   
   while (el[i]=(pelem->father)->sons[i]) i++;
   pel = forefather(pelem);
   NODES_OF_ELEMENT(*no1,*no2,*no3,pel);
   while( (*no1)->index < theGrid->minNodeIndex ) *no1 = (*no1)->son;
   while( (*no2)->index < theGrid->minNodeIndex ) *no2 = (*no2)->son;
   while( (*no3)->index < theGrid->minNodeIndex ) *no3 = (*no3)->son;
} 
 
void data2_from_forefather(theGrid,pelem,no1,no2,no3,el,fac,ind) 
GRID *theGrid;
ELEMENT *pelem, *el[8];
NODE **no1, **no2, **no3;
FACE *fac[4][4];
INT ind[4];
{   
   GRID *pg;
   ELEMENT *pel;
   FACE *fi;
   INT i=0;
   
   while (el[i]=(pelem->father)->sons[i]) i++;
   pel = forefather(pelem);
   NODES_OF_ELEMENT(*no1,*no2,*no3,pel);
   pg = theGrid;
   while(pg && pel->n[0]->index <= pg->maxNodeIndex) pg = pg->coarser;
   ind[0] = correct_ind(pg,pel->f[0],*no2,*no3);
   ind[1] = correct_ind(pg,pel->f[1],*no1,*no3);
   ind[2] = correct_ind(pg,pel->f[2],*no1,*no2);
   while( (*no1)->index < theGrid->minNodeIndex ) *no1 = (*no1)->son;
   while( (*no2)->index < theGrid->minNodeIndex ) *no2 = (*no2)->son;
   while( (*no3)->index < theGrid->minNodeIndex ) *no3 = (*no3)->son;
   for (i = 0; i < 4; i++){
      if (IS_TOP_FACE((fi=pel->f[i]))){  /*  father element of pelem is       */
         fac[i][0] = fi;                 /*                    a fictive one  */
         fac[i][1] = fac[i][2] = fac[i][3] = NULL;
      }
      else{
         fsons(fi,fac[i]);
         while(!fac[i][3] && fac[i][0]->sons[0])
            fsons_of_next_level(fac[i]);
         while(fac[i][0]->index < theGrid->minFaceIndex)/* -> fac[i][3]!=NULL */
            fsons_of_next_level(fac[i]);
         while(!fac[i][0]->sons[1] && fac[i][0]->sons[3])
            fsons_of_next_level(fac[i]); 
      }                                  
   }      
}  
    
INT small_edge_to_divide(n1,n2,no1,no2,no3,ind)
NODE *n1, *n2, *no1, *no2, *no3;
INT *ind;
{
   LINK *pli;
   
   pli = get_link(n1,n2);
   
   if ( !IS_EDGE_TO_DIVIDE(pli) ) return(0);
   else if ( IS_IN(n1,no1,no2,no3) && IS_IN(n2,no1,no2,no3) ){
      *ind = 1;         
      return(0);
   }
   else return(1); 
}
 
INT unif_refinement_necessary(pel,no1,no2,no3,ind)     /* pel is some element */
ELEMENT *pel;               /*  with status=1, no1, no2, no3 correspond       */
NODE *no1, *no2, *no3;      /*  to the first forefather element with status=0.*/
INT *ind;
{
   NODE *n1, *n2, *n3;
   
   NODES_OF_ELEMENT(n1,n2,n3,pel);
   if(      small_edge_to_divide(n1,n2,no1,no2,no3,ind) ) return(1);
   else if( small_edge_to_divide(n1,n3,no1,no2,no3,ind) ) return(1); 
   else if( small_edge_to_divide(n2,n3,no1,no2,no3,ind) ) return(1); 
   else return(0); 
}
 
INT check_green_elements(mg,el,no1,no2,no3,i1,i2)    /*  i1...green refinement*/
MULTIGRID *mg;                                       /*  i2...treated         */
ELEMENT *el[4];
NODE *no1, *no2, *no3;
INT i1, i2;
{
   INT i=0, ind1, ind2=0, flag;
   
   do 
      ind1 = unif_refinement_necessary(el[i++],no1,no2,no3,&ind2);
   while( el[i]!=NULL && !ind1 );
   if (ind1) {
      flag = 4;
      mark_all_element_edges(mg,no1,no2,no3);
   }
   else if (ind2) flag = i1;
   else flag = i2;
   for (i = 0; el[i] != NULL; i++) el[i]->eflag = flag;
   return(ind1);
}  
 
void mark_primary_unif_refinement(mg)   
MULTIGRID *mg;
{   
   NODE *n1, *n2, *n3;
   ELEMENT *pelem, *el[4];
   GRID *theGrid;
   INT i;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->eflag == 2){
            if (pelem->status == 0) {
               NODES_OF_ELEMENT(n1,n2,n3,pelem);
            }
            else{
               data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,el);
               for (i = 0; el[i] != NULL; i++) el[i]->eflag = 3;
            }
            mark_all_element_edges(mg,n1,n2,n3);
         }
}
          
INT mark_secondary_unif_refinement(mg,i1,i2) 
MULTIGRID *mg;
INT i1,i2;
{
   NODE *no1, *no2, *no3;
   ELEMENT *pelem, *el[4];
   GRID *theGrid;
   INT ind=0;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->status == 1 && pelem->eflag != 3
                                   && pelem->eflag != 4 && pelem->eflag < i1){
            data1_from_forefather(theGrid,pelem,&no1,&no2,&no3,el);
            if ( check_green_elements(mg,el,no1,no2,no3,i1,i2) ) ind = 1;
         }
   return(ind);
}
 
INT test_green_faces(mg,pnode,pnodeil,pnodebl,nodeNumber,i1)
MULTIGRID *mg;
NODE **pnode, **pnodeil, **pnodebl;
INT  *nodeNumber, i1;
{   
   NODE *n1, *n2, *n3, *n12, *n13, *n23;
   LINK *pl12, *pl13, *pl23;
   ELEMENT *pelem, *el[4];
   GRID *theGrid;
   INT i=0,j;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->eflag == i1){
            data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,el);    
            link_or_middle(n1,n2,&pl12,&n12);
            link_or_middle(n1,n3,&pl13,&n13);
            link_or_middle(n2,n3,&pl23,&n23);
            for (j = 0; el[j] != NULL; j++) el[j]->eflag = 1;
         }
   if(i) printf("%i edges divided for node numbering.\n",i);
   return(i);
}
 
void mark_new_green_elements_and_double_refinement(mg,i1) 
MULTIGRID *mg;
INT i1;
{
   NODE *n1, *n2, *n3;
   ELEMENT *pelem, *el[4];
   GRID *theGrid;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem))
            if (pelem->status == 0 && !pelem->eflag){
               NODES_OF_ELEMENT(n1,n2,n3,pelem);
               if ( IS_EDGE_TO_DIVIDE(get_link(n1,n2)) || 
                    IS_EDGE_TO_DIVIDE(get_link(n1,n3)) ||
                    IS_EDGE_TO_DIVIDE(get_link(n2,n3)) )
                  pelem->eflag = 1; 
            }
            else if (pelem->status == 1)
               if (pelem->eflag == 3){
                  data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,el);
                  check_green_elements(mg,el,n1,n2,n3,3,3);
               }
               else if (pelem->eflag == i1) pelem->eflag = 1;
               else if (pelem->eflag >  i1) pelem->eflag = 0;
}
 
/******************************************************************************/
 
 void treat_3edges(mg,n1,n2,n3,n11,n22,n33,n12,n13,n23,
                                   pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
 MULTIGRID *mg;
 NODE *n1, *n2, *n3, **n11, **n22, **n33, **n12, **n13, **n23, 
      **pnode, **pnodeil, **pnodebl;
 VERTEX **pvert;
 LINK **plink;
 INT  *nodeNumber;
 {
    *n12 = nu_treat_edge(mg,n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n13 = nu_treat_edge(mg,n1,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n23 = nu_treat_edge(mg,n2,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n11 = TOP_NODE(n1);
    *n22 = TOP_NODE(n2);
    *n33 = TOP_NODE(n3);
 }
 
 void reorder(n1,n2,n3,f1,f2,f3,n12,n13,n23,
              m1,m2,m3,e1,e2,e3,m12,m13,m23)
 NODE  *n1, *n2, *n3, *n12, *n13, *n23, 
      **m1,**m2,**m3,**m12,**m13,**m23;
 FACE  *f1, *f2, *f3,
      **e1,**e2,**e3;
 {
     *m1 = n1;    *m2 = n2;    *m3 = n3;
     *e1 = f1;    *e2 = f2;    *e3 = f3;
    *m12 = n12;  *m13 = n13;  *m23 = n23;  
 }
 
/* permutes the nodes n1, n2, n3 in such a way that
   n1 is a node such that the number of divided edges containig this node is
   maximal and (n1,n2) is a divided edge. The couple n1, n2 is the first 
   possible one with the mentioned properties.                             */

INT permute(n1,n2,n3,f1,f2,f3,n12,n13,n23)
NODE **n1, **n2, **n3, **n12, **n13, **n23;
FACE **f1, **f2, **f3;
{
   INT km1, km2, km3, l1, l2, l3, 
       lm1, lm2, lm3, k, l, km, lm;
   
   if ( !(*n12) && !(*n13) && !(*n23)) 
      return(0);
   else if (*n12 && *n13 && *n23) 
      return(0);
   else if (*n12 && *n13) 
      return(0);
   else if (*n12 && *n23) 
      reorder(*n2,*n1,*n3,*f2,*f1,*f3,*n12,*n23,*n13,
               n1, n2, n3, f1, f2, f3, n12, n13, n23);
   else if (*n13 && *n23) 
      reorder(*n3,*n1,*n2,*f3,*f1,*f2,*n13,*n23,*n12,
               n1, n2, n3, f1, f2, f3, n12, n13, n23);
   else if (*n12) 
      return(0);
   else if (*n13) 
      reorder(*n1,*n3,*n2,*f1,*f3,*f2,*n13,*n12,*n23,
               n1, n2, n3, f1, f2, f3, n12, n13, n23);
   else 
      reorder(*n2,*n3,*n1,*f2,*f3,*f1,*n23,*n12,*n13,
               n1, n2, n3, f1, f2, f3, n12, n13, n23);
   return(0);
}

 void nu_refine(k,theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,
      f1,f2,f3,f4,faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                       pelemil,pelembl,pelement)
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
 FACE *f1, *f2, *f3, *f4;
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE **pface, **pfaceil, **pfacebl;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT k, *faceNumber; 
 {
/*
    permute(&n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n23,&n24,&n34);
    if (n12 == NULL)
       finer0(k,theElement,n11,n22,n33,n44,f1,f2,f3,f4,pnflink,pfnlink,pflink);
    else if (n13 == NULL){
       if (n34 == NULL)
          finer1(k,theElement,n11,n22,n33,n44,n12,f1,f2,f3,f4,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
       else 
          finer2A(theElement,n11,n22,n33,n44,n12,n34,f1,f2,f3,f4,faceNumber,
                  pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,
                                                              pelembl,pelement);
    }
    else if (n14 == NULL){
       if (n23 == NULL){
          if (n24 == NULL){
             if (n34 == NULL)
                finer2B(k,theElement,n11,n22,n33,n44,n12,n13,f1,f2,f3,f4,
                        faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
             else 
                finer3A(theElement,n11,n33,n22,n44,n13,n12,n34,f1,f3,f2,f4,
                        faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
          }
          else if (n34==NULL)
             finer3A(theElement,n11,n22,n33,n44,n12,n13,n24,f1,f2,f3,f4,
                     faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
          else
             finer4A(theElement,n11,n22,n33,n44,n12,n13,n24,n34,f1,f2,f3,f4,
                     faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
       }
       else
          finer3B(theElement,n11,n22,n33,n44,n12,n13,n23,f1,f2,f3,f4,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    }
    else if (n23 == NULL){
       if (n24 == NULL)
          finer3C(k,theElement,n11,n22,n33,n44,n12,n13,n14,f1,f2,f3,f4,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
       else 
          finer4B(theElement,n11,n22,n44,n33,n12,n14,n13,n24,f1,f2,f4,f3,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    }
    else if (n24 == NULL)
       finer4B(theElement,n11,n22,n33,n44,n12,n13,n14,n23,f1,f2,f3,f4,
               faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else if (n34 == NULL)
       finer5(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else
       finer6(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
*/
 }  
 
 void treat_subface(ref,f,i1,j1,facs,fath,n11,n12,ind,pface,pfaceil,pfacebl,
                                                             pnflink,faceNumber)
 FACE *f, *facs[2], *fath[2], **pface, **pfaceil, **pfacebl;
 NODE *n11, *n12;
 NFLINK **pnflink;
 INT ref, i1, j1, *faceNumber;
 INT ind;
 {
    if (facs[i1]){
       f->sons[j1] = facs[i1];
       if (ref == 3 && !facs[i1]->sons[0]) facs[i1]->sons[0] = facs[i1];
    }
    else if (ref == 3 && (!get_link(n11,n12))
       f->sons[j1] = fath[i1]; 
    else{
       if (NOT_FF(fath[i1]))
          *pfacebl = (*pfacebl)->succ = *pface;
       else
          *pfaceil = (*pfaceil)->succ = *pface;
       f->sons[j1] = *pface;
       fmake(f->sons[j1],++(*pface),FTYPE(fath[i1]),0,++(*faceNumber),fath[i1]);
       SET_ORIENT(f->sons[j1],norm_vect_or(n11,n12,ind));
       fath[i1]->sons[0] = f->sons[j1];
       if (ref == 2 && IS_FF(f->sons[j1])){
          new_Fneighbour(n11,f->sons[j1],pnflink);
          new_Fneighbour(n12,f->sons[j1],pnflink);
       }
    }
 }

void copy_faces(mg,pg,ref,fac,ind,ff,f,n1,n2,n11,n22,n12,pnode,
pnodeil,pnodebl,pvert,plink,pface,pfaceil,pfacebl,pnflink,nodeNumber,faceNumber)
MULTIGRID *mg;
GRID *pg;
NODE *n1, *n2, *n11, *n22, *n12, **pnode, **pnodeil, **pnodebl;
FACE *fac[2], *ff, **f, **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
NFLINK **pnflink;
INT ref, *faceNumber, *nodeNumber;
INT ind;
{
   FACE *facs[2], *fath[2];
   INT i, i1, i2, i3, j1, j2, j3;    
   
   *f = ff;
   (*f)->type = FTYPE(fac[0]);
   SET_ORIENT(*f,ind);
   (*f)->sons[1] = NULL;
   if (fac[1] == NULL){
      if(FACE_ON_NEW_LEVEL(fac[0],mg)){ 
         (*f)->sons[0] = fac[0];
         (*f)->father = *f; 
      }
      else{
         (*f)->sons[0] = NULL;
         (*f)->father = NULL;
      }
   }
   else{  /* fac[1] != NULL */
      if(FACE_ON_NEW_LEVEL(fac[0],mg) || fac[0]->index < 0){
         (*f)->sons[0] = fac[0];
         (*f)->sons[1] = fac[1];
         if (ref == 3){
            if (!fac[0]->sons[0]) fac[0]->sons[0] = fac[0];
            if (!fac[1]->sons[0]) fac[1]->sons[0] = fac[1];
         }
      }
      else{  /* the corresponding face of the refined "macroelement"  */   
         for (i = 0; i < 2; i++){         /*  consists of 2 subfaces  */
            facs[i] = fac[i];
            while (facs[i]->sons[0] && !facs[i]->sons[1]) 
               facs[i] = facs[i]->sons[0];
            if (!facs[i]->sons[0] && !FACE_ON_NEW_LEVEL(facs[i],mg)){
               fath[i] = facs[i];
               facs[i] = NULL;
            }
         }
         while (pg->maxFaceIndex < fac[0]->index) pg = pg->finer;
         while (n1->index < pg->minNodeIndex) n1 = n1->son;
         while (n2->index < pg->minNodeIndex) n2 = n2->son;
         find_order(n1,n2,n3,&i1,&i2,&i3);    /* fac[ik] lies at nk    */
         find_order(n11,n22,n33,&j1,&j2,&j3); /* nkk is the jk-th node */
         treat_subface(ref,*f,i1,j1,facs,fath,n11,n12,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         treat_subface(ref,*f,i2,j2,facs,fath,n22,n12,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         if (ref == 2){
            fac[0]->sons[3] = (*f)->sons[0]; 
            fac[1]->sons[3] = (*f)->sons[1]; 
         }
      }
      (*f)->father = *f; 
   }
}
 
 void fcopyS2(i1,i2,j1,j2,f,fac)
 INT i1, i2, j1, j2;
 FACE *f, *fac[4];
 {
    fac[i1]->sons[0] = f->sons[j1];
    fac[i2]->sons[0] = f->sons[j2];
    f->sons[j1]->father = fac[i1];
    f->sons[j2]->father = fac[i2];
 }
 
 void fcopy2(n1,n2,n11,n22,f,fac)
 NODE *n1, *n2, *n11, *n22;
 FACE *f, *fac[4];
 {
    if (n1->index < n2->index && n11->index < n22->index ||
        n1->index > n2->index && n11->index > n22->index )
       fcopyS2(0,1,0,1,f,fac);
    else
       fcopyS2(0,1,1,0,f,fac);
    fac[0]->sons[3] = f->sons[0];
    fac[1]->sons[3] = f->sons[1];   
 }
 
 void fcopy3(n2,n3,n22,n33,f,fac)
 NODE *n2, *n3, *n22, *n33;
 FACE *f, *fac[4];
 {
    fac[0]->sons[0] = f->sons[0];
    f->sons[0]->father = fac[0];
    if (n2->index < n3->index && n22->index < n33->index ||
        n2->index > n3->index && n22->index > n33->index )
       fcopyS2(1,2,1,2,f,fac);
    else
       fcopyS2(1,2,2,1,f,fac);
    fac[0]->sons[3] = f->sons[0];
    fac[1]->sons[3] = f->sons[1];   
    fac[2]->sons[3] = f->sons[2];   
 }
 
 void recopy_faces(pg,n1,n2,n11,n22,n12,f,fac)
 GRID *pg;
 NODE *n1, *n2, *n11, *n22, *n12;
 FACE *f, *fac[4];
 {
    INT i, j=0, k=0;
    
    if (f->father == NULL){
       set_norm_vect_or(n11,n22,n12,f);
       while(fac[++j]);
       while(++k < 4 && f->sons[k]);
       if (j < k){
          for (i = 0; i < k; i++){
             fac[0]->sons[i] = f->sons[i];
             f->sons[i]->father = fac[0];
          }
          for (i = 1; i < j; i++)
             fac[i]->sons[0] = f->sons[0];
       }
       else   /*  j == k  */
          if (j < 2){ 
             fac[0]->sons[0] = f->sons[0];
             f->sons[0]->father = fac[0];
          }
          else{  
             while(pg->maxFaceIndex < fac[0]->index) pg = pg->finer;
             while(n1->index < pg->minNodeIndex) n1 = n1->son;
             while(n2->index < pg->minNodeIndex) n2 = n2->son;
             if (j < 3)
                if (n12) 
                   fcopy2(n1,n2,n11,n22,f,fac);
             else if (j < 4)
                if (!n12)
                   fcopy3(n1,n2,n11,n22,f,fac);
          }
    }
 }
 
 void set_ind1(n11,n22,n12,ff,f)
 NODE *n11, *n22, *n12;
 FACE *ff, *f;
 {
    UNS type;
    
    type = FTYPE(f);
    FTYPE(f) = FTYPE(ff);
    set_norm_vect_or(n11,n22,n12,f); 
    FTYPE(f) = type;
 }
 
/* refinement of an element with status=0 */
void refine1(mg,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,
          pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
MULTIGRID *mg;
ELEMENT *theElem, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
    FACE *f1, *f2, *f3, *f4;
    
    if(FACE_ON_NEW_LEVEL(f1=top_face(theElem->f[0]),mg)) f1 = f1->father;
    if(FACE_ON_NEW_LEVEL(f2=top_face(theElem->f[1]),mg)) f2 = f2->father;
    if(FACE_ON_NEW_LEVEL(f3=top_face(theElem->f[2]),mg)) f3 = f3->father;
    if(FACE_ON_NEW_LEVEL(f4=top_face(theElem->f[3]),mg)) f4 = f4->father;
    treat_6edges(mg,theElem->n[0],theElem->n[1],theElem->n[2],theElem->n[3],
                 &n11,&n22,&n33,&n44,&n12,&n13,&n14,&n23,&n24,&n34,
                 pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    if (theElem->eflag>1)
       finer6(theElem,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else
       nu_refine(0,theElem,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
                 faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    set_ind1(n11,n22,n12,theElem->f[2],f3); 
    set_ind1(n11,n33,n13,theElem->f[1],f2);  
    set_ind1(n22,n33,n23,theElem->f[0],f1);    
    theElem->eflag = 0; 
 }
 
 /* single refinement of elements with status=1 */
 void refine2(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
              pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,
                                                          nodeNumber,faceNumber)
 MULTIGRID *mg;
 GRID *theGrid;
 ELEMENT *theElem, **pelemil, **pelembl, **pelement;
 NODE **pnode, **pnodeil, **pnodebl;
 FACE **pface, **pfaceil, **pfacebl;
 VERTEX **pvert;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT  *nodeNumber, *faceNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34,
         *no1, *no2, *no3, *no4;
    FACE ff1, ff2, ff3, ff4, *f1, *f2, *f3, *f4, *fac[4][4];
    ELEMENT *el[8];
    INT ind[4];
    INT i;
    
    data2_from_forefather(theGrid,theElem,&no1,&no2,&no3,el,fac,ind);
    treat_6edges(mg,no1,no2,no3,no4,&n11,&n22,&n33,&n44,&n12,&n13,&n14,
                   &n23,&n24,&n34,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
void copy_faces(mg,pg,ref,fac,ind,ff,f,n1,n2,n11,n22,n12,pnode,
    copy_faces(mg,theGrid,2,fac[0],ind[0],&ff1,&f1,no2,no3,
                n22,n33,n23,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,2,fac[1],ind[1],&ff2,&f2,no1,no3,
                n11,n33,n13,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,2,fac[2],ind[2],&ff3,&f3,no1,no2,
                n11,n22,n12,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    nu_refine(0,el[0],n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    for (i=1; el[i]; i++){
       el[i]->sons[0] = el[0]->sons[0];
       el[i]->eflag = 0;
    }
    el[0]->eflag = 0;
    recopy_faces(theGrid,no2,no3,n22,n33,n23,f1,fac[0]);
    recopy_faces(theGrid,no1,no3,n11,n33,n13,f2,fac[1]);
    recopy_faces(theGrid,no1,no2,n11,n22,n12,f3,fac[2]);
    for (i=0; el[i]; i++){
       if (el[i]->f[0]->sons[0] == NULL) el[i]->f[0]->sons[0] = el[i]->f[0];
       if (el[i]->f[1]->sons[0] == NULL) el[i]->f[1]->sons[0] = el[i]->f[1];
       if (el[i]->f[2]->sons[0] == NULL) el[i]->f[2]->sons[0] = el[i]->f[2];
       if (el[i]->f[3]->sons[0] == NULL) el[i]->f[3]->sons[0] = el[i]->f[3];
    }
 }
 
 void set_ind2(mg,n11,n22,n12,f)
 MULTIGRID *mg;
 NODE *n11, *n22, *n12;
 FACE *f;
 {
    UNS type;
    
    type = FTYPE(f);
    if(f->sons[1] && f->index > 0)
       SET_ORIENT(f,correct_ind(TOP_GRID(mg),f,n11,n22));
    set_norm_vect_or(n11,n22,n12,f); 
    FTYPE(f) = type;
 }
 
/* refinement of products of finer4 */
void second_refinement(mg,theElem,peil,pebl,pface,pfaceil,pfacebl,plink,pnflink,
                pfnlink,pflink,pelemil,pelembl,pelement,faceNumber,fictElNumber)
MULTIGRID *mg;
ELEMENT *theElem, *peil, *pebl, **pelemil, **pelembl, **pelement;
FACE **pface, **pfaceil, **pfacebl;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *faceNumber, *fictElNumber; 
 {
/*
  PREDELAT !!!
*/
    NODE *n11, *n22, *n33, *n12, *n13, *n23;
    
    NODES_OF_ELEMENT(n11,n22,n33,theElem);
    n12 = middle(n11,n22); 
    n13 = middle(n11,n33);  
    n23 = middle(n22,n33);  
    nu_refine(1,theElem,n11,n22,n33,n12,n13,n23,
         theElem->f[0],theElem->f[1],theElem->f[2],faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
    if (theElem->sons[0]) 
       treat_fictive_element(theElem,peil,pebl,fictElNumber);
    set_ind2(mg,n11,n22,n12,theElem->f[2]); 
    set_ind2(mg,n11,n33,n13,theElem->f[1]); 
    set_ind2(mg,n22,n33,n23,theElem->f[0]); 
 }

 /* double refinement of elements with status=1 */
 void refine3(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
         pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,
                                             faceNumber,fictNumber,fictElNumber)
 MULTIGRID *mg;
 GRID *theGrid;
 ELEMENT *theElem, **pelemil, **pelembl, **pelement;
 NODE **pnode, **pnodeil, **pnodebl;
 FACE **pface, **pfaceil, **pfacebl;
 VERTEX **pvert;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT  *nodeNumber, *faceNumber, *fictNumber, *fictElNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34,
         *no1, *no2, *no3, *no4;
    FACE *pfil, *pfbl, ff1, ff2, ff3, ff4, *f1, *f2, *f3, *f4, *fac[4][4];
    ELEMENT *peil, *pebl, *el[4], *eson;
    INT ind[2];
    INT i,j;
    
    pfil = *pfaceil;
    pfbl = *pfacebl;
    peil = *pelemil;
    pebl = *pelembl;
    data2_from_forefather(theGrid,theElem,&no1,&no2,&no3,el,fac,ind);
    treat_6edges(mg,no1,no2,no3,no4,&n11,&n22,&n33,&n44,&n12,&n13,&n14,
                   &n23,&n24,&n34,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    copy_faces(mg,theGrid,3,fac[0],ind[0],&ff1,&f1,no2,no3,no4,
               n22,n33,n44,n23,n24,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[1],ind[1],&ff2,&f2,no1,no3,no4,
               n11,n33,n44,n13,n14,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[2],ind[2],&ff3,&f3,no1,no2,no4,
               n11,n22,n44,n12,n14,n24,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[3],ind[3],&ff4,&f4,no1,no2,no3,
               n11,n22,n33,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    Finer6(mg,el[0],n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,pnode,
           pnodeil,pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,
                                              pfacebl,pelemil,pelembl,pelement);
    for (i=1; el[i]; i++){
       el[i]->sons[0] = el[0]->sons[0];
       el[i]->eflag = 0;
    }
    el[0]->eflag = 0;
    recopy_faces(theGrid,no2,no3,n22,n33,n23,f1,fac[0]);
    recopy_faces(theGrid,no1,no3,n11,n33,n13,f2,fac[1]);
    recopy_faces(theGrid,no1,no2,n11,n22,n12,f3,fac[2]);
    for (i=0; i<4; i++)
      second_refinement(mg,el[0]->sons[i],peil,pebl,pface,pfaceil,pfacebl,plink,
       pnflink,pfnlink,pflink,pelemil,pelembl,pelement,faceNumber,fictElNumber);
    for (i=0; i<4; i++){
       eson = el[0]->sons[i];
       if (eson->f[0]->sons[0] == eson->f[0])
           eson->f[0]->sons[0] = NULL;
       else if (eson->f[0]->sons[0] && 
                eson->f[0]->index > TOP_GRID(mg)->maxFaceIndex)
          treat_fictive_face(eson->f[0],pfil,pfbl,fictNumber);
       if (eson->f[1]->sons[0] == eson->f[1])
           eson->f[1]->sons[0] = NULL;
       else if (eson->f[1]->sons[0] && 
                eson->f[1]->index > TOP_GRID(mg)->maxFaceIndex)
          treat_fictive_face(eson->f[1],pfil,pfbl,fictNumber);
    }
    for (i=0; el[i]; i++){
       if (el[i]->f[0]->sons[0] == NULL) el[i]->f[0]->sons[0] = el[i]->f[0];
       if (el[i]->f[1]->sons[0] == NULL) el[i]->f[1]->sons[0] = el[i]->f[1];
       if (el[i]->f[2]->sons[0] == NULL) el[i]->f[2]->sons[0] = el[i]->f[2];
    }
 }
 
 void additional_neighbours(mg,n1,n2,n3,f1,f2,f3,plink,pnflink,pfnlink,pflink)
 MULTIGRID *mg;
 NODE *n1, *n2, *n3;
 FACE *f1, *f2, *f3;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 {
    FNLINK *pfnl;
    FLINK *pfl;
    
    if (NODE_ON_NEW_LEVEL(n1,mg) && IS_FN(n1)){
       if (IS_FN(n2) && !Get_link(n1,n2))
          new_neighbour(n1,n2,2,0,plink);
       if (IS_FN(n3) && !Get_link(n1,n3))
          new_neighbour(n1,n3,2,0,plink);
       if (IS_FF(f1))
          new_Fneighbour(n1,f1,pnflink);
       if (IS_FF(f2) && !FACE_ON_NEW_LEVEL(f2,mg) && !get_nflink(n1,f2))
          new_Fneighbour(n1,f2,pnflink);
       if (IS_FF(f3) && !FACE_ON_NEW_LEVEL(f3,mg) && !get_nflink(n1,f3))
          new_Fneighbour(n1,f3,pnflink);
    }
    if (IS_FF(f1) && FACE_ON_NEW_LEVEL(f1,mg)){
       if (IS_FN(n1)){
          if (f1->fnstart == NULL)
             f1->fnstart = *pfnlink;
          else{
             for (pfnl = f1->fnstart; pfnl->next != NULL; pfnl = pfnl->next);
             pfnl->next = *pfnlink;
          }
          (*pfnlink)->nbnode = n1;
          ((*pfnlink)++)->next = NULL;
       }
       if (IS_FF(f2) || IS_FF(f3)){
          if (f1->fstart == NULL)
             f1->fstart = *pflink;
          else{
             for (pfl = f1->fstart; pfl->next != NULL; pfl = pfl->next);
             pfl->next = *pflink;
          }
          if (IS_FF(f2)){
             (*pflink)->next = *pflink + 1;
             ((*pflink)++)->nbface = f2;
          } 
          if (IS_FF(f3)){
             (*pflink)->next = *pflink + 1;
             ((*pflink)++)->nbface = f3;
          } 
          (*pflink - 1)->next = NULL;
       }
    }
 }
          
 void add_neighbours_from_old_levels(mg,plink,pnflink,pfnlink,pflink)
 MULTIGRID *mg;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 {         
    NODE *n1, *n2, *n3;
    FACE *f1, *f2, *f3;
    ELEMENT *pel;
    GRID *theGrid;
    
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
          if (IS_TOP_ELEMENT(pel)){
             TOPNODES_OF_ELEMENT(n1,n2,n3,pel); 
             TOPFACES_OF_ELEMENT(f1,f2,f3,pel);  
             additional_neighbours(mg,n1,n2,n3,f1,f2,f3,plink,pnflink,
                                                                pfnlink,pflink);
             additional_neighbours(mg,n2,n3,n1,f2,f3,f1,plink,pnflink,
                                                                pfnlink,pflink);
             additional_neighbours(mg,n3,n1,n2,f3,f1,f2,plink,pnflink,
                                                                pfnlink,pflink);
          }
 }

#elif !(DATA_STR & REDUCED) && (ELEMENT_TYPE == SIMPLEX) && (DIM == 3)

void mark_all_element_edges(mg,n1,n2,n3,n4) 
MULTIGRID *mg;
NODE *n1, *n2, *n3, *n4;
{
   mark_edge_to_divide(mg,n1,n2);
   mark_edge_to_divide(mg,n1,n3);      
   mark_edge_to_divide(mg,n1,n4);      
   mark_edge_to_divide(mg,n2,n3);                                        
   mark_edge_to_divide(mg,n2,n4);
   mark_edge_to_divide(mg,n3,n4);
} 
 
void fsons(f,fac)
FACE *f, *fac[4];
{
   fac[0] = f->sons[0];
   fac[1] = f->sons[1];
   if (fac[2] = f->sons[2])
      fac[3] = f->sons[3];
   else
      fac[3] = NULL;
}
 
void fsons_of_next_level(fac)  
FACE *fac[4];
{
   INT i = -1;
   
   if (fac[0]->sons[1])  /*  geometrical form changed  */
      fsons(fac[0],fac);
   else if (fac[1])
      while (++i < 4 && fac[i]) fac[i] = fac[i]->sons[3]; 
   else 
      fac[0] = fac[0]->sons[0];
}
    
INT correct_ind(pg,f,n1,n2,n3)
GRID *pg;
FACE *f;
NODE *n1, *n2, *n3;
{ 
   if (pg && f->index <= pg->maxFaceIndex && f->index > 0){
      while(f->index < pg->minFaceIndex) pg = pg->coarser;
      while(n1->index > pg->maxNodeIndex) n1 = n1->father;
      while(n2->index > pg->maxNodeIndex) n2 = n2->father;
      while(n3->index > pg->maxNodeIndex) n3 = n3->father;
      return(norm_vect_or(n1,n2,n3,FTYPE(f)&ORIENT));
   }
   else
      return(FTYPE(f)&ORIENT);
}
       
void data1_from_forefather(theGrid,pelem,no1,no2,no3,no4,el)   
GRID *theGrid;
ELEMENT *pelem, *el[8];     /* pelem is a green element */
NODE **no1, **no2, **no3, **no4;
{   
   ELEMENT *pel;
   INT i=0;
   
   while (el[i]=(pelem->father)->sons[i]) i++;
   pel = forefather(pelem);
   NODES_OF_ELEMENT(*no1,*no2,*no3,*no4,pel);
   while( (*no1)->index < theGrid->minNodeIndex ) *no1 = (*no1)->son;
   while( (*no2)->index < theGrid->minNodeIndex ) *no2 = (*no2)->son;
   while( (*no3)->index < theGrid->minNodeIndex ) *no3 = (*no3)->son;
   while( (*no4)->index < theGrid->minNodeIndex ) *no4 = (*no4)->son;
} 
 
void data2_from_forefather(theGrid,pelem,no1,no2,no3,no4,el,fac,ind) 
GRID *theGrid;
ELEMENT *pelem, *el[8];
NODE **no1, **no2, **no3, **no4;
FACE *fac[4][4];
INT ind[4];
{   
   GRID *pg;
   ELEMENT *pel;
   FACE *fi;
   INT i=0;
   
   while (el[i]=(pelem->father)->sons[i]) i++;
   pel = forefather(pelem);
   NODES_OF_ELEMENT(*no1,*no2,*no3,*no4,pel);
   pg = theGrid;
   while(pg && pel->n[0]->index <= pg->maxNodeIndex) pg = pg->coarser;
   ind[0] = correct_ind(pg,pel->f[0],*no2,*no3,*no4);
   ind[1] = correct_ind(pg,pel->f[1],*no1,*no3,*no4);
   ind[2] = correct_ind(pg,pel->f[2],*no1,*no2,*no4);
   ind[3] = correct_ind(pg,pel->f[3],*no1,*no2,*no3);
   while( (*no1)->index < theGrid->minNodeIndex ) *no1 = (*no1)->son;
   while( (*no2)->index < theGrid->minNodeIndex ) *no2 = (*no2)->son;
   while( (*no3)->index < theGrid->minNodeIndex ) *no3 = (*no3)->son;
   while( (*no4)->index < theGrid->minNodeIndex ) *no4 = (*no4)->son;
   for (i = 0; i < 4; i++){
      if (IS_TOP_FACE((fi=pel->f[i]))){  /*  father element of pelem is       */
         fac[i][0] = fi;                 /*                    a fictive one  */
         fac[i][1] = fac[i][2] = fac[i][3] = NULL;
      }
      else{
         fsons(fi,fac[i]);
         while(!fac[i][3] && fac[i][0]->sons[0])
            fsons_of_next_level(fac[i]);
         while(fac[i][0]->index < theGrid->minFaceIndex)/* -> fac[i][3]!=NULL */
            fsons_of_next_level(fac[i]);
         while(!fac[i][0]->sons[1] && fac[i][0]->sons[3])
            fsons_of_next_level(fac[i]); 
      }                                  
   }      
}  
    
INT small_edge_to_divide(n1,n2,no1,no2,no3,no4,ind)
NODE *n1, *n2, *no1, *no2, *no3, *no4;
INT *ind;
{
   LINK *pli;
   
   pli = get_link(n1,n2);
   
   if ( !IS_EDGE_TO_DIVIDE(pli) ) return(0);
   else if ( IS_IN(n1,no1,no2,no3,no4) && IS_IN(n2,no1,no2,no3,no4) ){
      *ind = 1;         
      return(0);
   }
   else return(1); 
}
 
INT unif_refinement_necessary(pel,no1,no2,no3,no4,ind) /* pel is some element */
ELEMENT *pel;               /*  with status=1, no1, no2, no3, no4 correspond  */
NODE *no1, *no2, *no3, *no4;/*  to the first forefather element with status=0.*/
INT *ind;
{
   NODE *n1, *n2, *n3, *n4;
   
   NODES_OF_ELEMENT(n1,n2,n3,n4,pel);
   if(      small_edge_to_divide(n1,n2,no1,no2,no3,no4,ind) ) return(1);
   else if( small_edge_to_divide(n1,n3,no1,no2,no3,no4,ind) ) return(1); 
   else if( small_edge_to_divide(n1,n4,no1,no2,no3,no4,ind) ) return(1); 
   else if( small_edge_to_divide(n2,n3,no1,no2,no3,no4,ind) ) return(1); 
   else if( small_edge_to_divide(n2,n4,no1,no2,no3,no4,ind) ) return(1); 
   else if( small_edge_to_divide(n3,n4,no1,no2,no3,no4,ind) ) return(1); 
   else return(0); 
}
 
INT check_green_elements(mg,el,no1,no2,no3,no4,i1,i2)/*  i1...green refinement*/
MULTIGRID *mg;                                       /*  i2...treated         */
ELEMENT *el[8];
NODE *no1, *no2, *no3, *no4;
INT i1,i2;
{
   INT i=0, ind1, ind2=0, flag;
   
   do 
      ind1 = unif_refinement_necessary(el[i++],no1,no2,no3,no4,&ind2);
   while( el[i]!=NULL && !ind1 );
   if (ind1) {
      flag = 4;
      mark_all_element_edges(mg,no1,no2,no3,no4);
   }
   else if (ind2) flag = i1;
   else flag = i2;
   for (i = 0; el[i] != NULL; i++) el[i]->eflag = flag;
   return(ind1);
}  
 
void mark_primary_unif_refinement(mg)   
MULTIGRID *mg;
{   
   NODE *n1, *n2, *n3, *n4;
   ELEMENT *pelem, *el[8];
   GRID *theGrid;
   INT i;
   
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->eflag == 2){
            if (pelem->status == 0) {
               NODES_OF_ELEMENT(n1,n2,n3,n4,pelem);
            }
            else{
               data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,&n4,el);
               for (i = 0; el[i] != NULL; i++) el[i]->eflag = 3;
            }
            mark_all_element_edges(mg,n1,n2,n3,n4);
         }
}
          
INT mark_secondary_unif_refinement(mg,i1,i2) 
MULTIGRID *mg;
INT i1,i2;
{
   NODE *no1, *no2, *no3, *no4;
   ELEMENT *pelem, *el[8];
   GRID *theGrid;
   INT ind=0;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->status == 1 && pelem->eflag != 3
                                   && pelem->eflag != 4 && pelem->eflag < i1){
            data1_from_forefather(theGrid,pelem,&no1,&no2,&no3,&no4,el);
            if ( check_green_elements(mg,el,no1,no2,no3,no4,i1,i2) ) ind = 1;
         }
   return(ind);
}
 
void check_face(mg,n1,n2,n3,n12,n13,n23,pl12,pl13,pl23,i,pnode,pnodeil,pnodebl,
                                                                     nodeNumber)
MULTIGRID *mg;
NODE *n1, *n2, *n3, *n12, *n13, *n23, **pnode, **pnodeil, **pnodebl;
LINK *pl12, *pl13, *pl23;
INT  *nodeNumber, *i;
{
   if (n12){
     if (n13){
        if (!n23)  check_edge(mg,n2,n3,pl23,i,pnode,pnodeil,pnodebl,nodeNumber);
     }
     else if (n23) check_edge(mg,n1,n3,pl13,i,pnode,pnodeil,pnodebl,nodeNumber);
   }
   else if (n13 && n23) 
                   check_edge(mg,n1,n2,pl12,i,pnode,pnodeil,pnodebl,nodeNumber);
}       

INT test_green_faces(mg,pnode,pnodeil,pnodebl,nodeNumber,i1)
MULTIGRID *mg;
NODE **pnode, **pnodeil, **pnodebl;
INT  *nodeNumber, i1;
{   
   NODE *n1, *n2, *n3, *n4, *n12, *n13, *n14, *n23, *n24, *n34;
   LINK *pl12, *pl13, *pl14, *pl23, *pl24, *pl34;
   ELEMENT *pelem, *el[8];
   GRID *theGrid;
   INT i=0,j;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem) && pelem->eflag == i1){
            data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,&n4,el);    
            link_or_middle(n1,n2,&pl12,&n12);
            link_or_middle(n1,n3,&pl13,&n13);
            link_or_middle(n1,n4,&pl14,&n14);
            link_or_middle(n2,n3,&pl23,&n23);
            link_or_middle(n2,n4,&pl24,&n24);
            link_or_middle(n3,n4,&pl34,&n34);
            check_face(mg,n1,n2,n3,n12,n13,n23,pl12,pl13,pl23,&i,pnode,pnodeil,
                                                            pnodebl,nodeNumber);
            check_face(mg,n1,n2,n4,n12,n14,n24,pl12,pl14,pl24,&i,pnode,pnodeil,
                                                            pnodebl,nodeNumber);
            check_face(mg,n1,n3,n4,n13,n14,n34,pl13,pl14,pl34,&i,pnode,pnodeil,
                                                            pnodebl,nodeNumber);
            check_face(mg,n2,n3,n4,n23,n24,n34,pl23,pl24,pl34,&i,pnode,pnodeil,
                                                            pnodebl,nodeNumber);
            for (j = 0; el[j] != NULL; j++) el[j]->eflag = 1;
         }
   if(i) printf("%i edges divided for node numbering.\n",i);
   return(i);
}
 
void mark_new_green_elements_and_double_refinement(mg,i1) 
MULTIGRID *mg;
INT i1;
{
   NODE *n1, *n2, *n3, *n4;
   ELEMENT *pelem, *el[8];
   GRID *theGrid;
         
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ)
         if (IS_TOP_ELEMENT(pelem))
            if (pelem->status == 0 && !pelem->eflag){
               NODES_OF_ELEMENT(n1,n2,n3,n4,pelem);
               if ( IS_EDGE_TO_DIVIDE(get_link(n1,n2)) || 
                    IS_EDGE_TO_DIVIDE(get_link(n1,n3)) ||
                    IS_EDGE_TO_DIVIDE(get_link(n1,n4)) || 
                    IS_EDGE_TO_DIVIDE(get_link(n2,n3)) ||
                    IS_EDGE_TO_DIVIDE(get_link(n2,n4)) || 
                    IS_EDGE_TO_DIVIDE(get_link(n3,n4)) )
                  pelem->eflag = 1; 
            }
            else if (pelem->status == 1)
               if (pelem->eflag == 3){
                  data1_from_forefather(theGrid,pelem,&n1,&n2,&n3,&n4,el);
                  check_green_elements(mg,el,n1,n2,n3,n4,3,3);
               }
               else if (pelem->eflag == i1) pelem->eflag = 1;
               else if (pelem->eflag >  i1) pelem->eflag = 0;
}
 
/******************************************************************************/
 
 void treat_6edges(mg,n1,n2,n3,n4,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,
                                   pnode,pnodeil,pnodebl,pvert,plink,nodeNumber)
 MULTIGRID *mg;
 NODE *n1, *n2, *n3, *n4, **n11, **n22, **n33, **n44, 
      **n12, **n13, **n14, **n23, **n24, **n34, **pnode, **pnodeil, **pnodebl;
 VERTEX **pvert;
 LINK **plink;
 INT  *nodeNumber;
 {
    *n12 = nu_treat_edge(mg,n1,n2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n13 = nu_treat_edge(mg,n1,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n14 = nu_treat_edge(mg,n1,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n23 = nu_treat_edge(mg,n2,n3,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n24 = nu_treat_edge(mg,n2,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n34 = nu_treat_edge(mg,n3,n4,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    *n11 = TOP_NODE(n1);
    *n22 = TOP_NODE(n2);
    *n33 = TOP_NODE(n3);
    *n44 = TOP_NODE(n4);
 }
 
 void reorder(n1,n2,n3,n4,f1,f2,f3,f4,n12,n13,n14,n23,n24,n34,
             m1,m2,m3,m4,e1,e2,e3,e4,m12,m13,m14,m23,m24,m34)
 NODE  *n1, *n2, *n3, *n4, *n12, *n13, *n14, *n23, *n24, *n34,
      **m1,**m2,**m3,**m4,**m12,**m13,**m14,**m23,**m24,**m34;
 FACE  *f1, *f2, *f3, *f4,
      **e1,**e2,**e3,**e4;
 {
     *m1 = n1;    *m2 = n2;    *m3 = n3;   *m4 = n4;
     *e1 = f1;    *e2 = f2;    *e3 = f3;   *e4 = f4;
    *m12 = n12;  *m13 = n13;  *m14 = n14; 
    *m23 = n23;  *m24 = n24;  *m34 = n34;
 }
 
/* permutes the nodes n1, n2, n3, n4 in such a way that
   n1 is a node such that the number of divided edges containig this node is
   maximal;
   (n1,n2) is a divided edge and there do not exist any nodes n, m such that 
   (n,m) is a divided edge, the number of divided edges containing n is the same
   as for n1 and the number of divided edges containing m is grater than for n2.
   The couple n1, n2 is the first possible one with the mentioned properties.
   If n13 = 0 than also n14 = 0.
   If n13 = 0 than clearly n1->index < n2->index.                             */

INT permute(n1,n2,n3,n4,f1,f2,f3,f4,n12,n13,n14,n23,n24,n34)
NODE **n1, **n2, **n3, **n4, **n12, **n13, **n14, **n23, **n24, **n34;
FACE **f1, **f2, **f3, **f4;
{
   INT km1, km2, km3, km4, l1, l2, l3, l4, 
       lm1, lm2, lm3, lm4, k, l, km, lm;
   
   if ( !(*n12) && !(*n13) && !(*n14) && !(*n23) && !(*n24) && !(*n34)) 
      return(0);
   if ( *n12 && *n13 && *n14 && *n23 && *n24 && *n34) return(0);
   km1 = ((long)(*n12) > 0) + ((long)(*n13) > 0) + ((long)(*n14) > 0);/*number of  */
   km2 = ((long)(*n12) > 0) + ((long)(*n23) > 0) + ((long)(*n24) > 0);/*divided ed-*/
   km3 = ((long)(*n13) > 0) + ((long)(*n23) > 0) + ((long)(*n34) > 0);/*ges contai-*/
   km4 = ((long)(*n14) > 0) + ((long)(*n24) > 0) + ((long)(*n34) > 0);/*ning *n1   */
   l1 = max3( km2*((long)(*n12)>0), km3*((long)(*n13)>0), km4*((long)(*n14)>0),
                                                                 2, 3, 4, &lm1);
    /* *n_l1 is the first node such that (*n1,*n_l1) is a divided edge and the 
       number lm1 of divided edges containing *n_l1 is maximal. */
   l2 = max3( km1*((long)(*n12)>0), km3*((long)(*n23)>0), km4*((long)(*n24)>0),
                                                                 1, 3, 4, &lm2);
   l3 = max3( km1*((long)(*n13)>0), km2*((long)(*n23)>0), km4*((long)(*n34)>0),
                                                                 1, 2, 4, &lm3);
   l4 = max3( km1*((long)(*n14)>0), km2*((long)(*n24)>0), km3*((long)(*n34)>0), 
                                                                 1, 2, 3, &lm4);
   k = 1;
   l = l1;
   km = km1;
   lm = lm1;
   if ( km2 > km1 || (km2 == km1 && lm2 > lm1) ){
      k = 2;
      l = l2;
      km = km2;
      lm = lm2;
   }
   if ( km3 > km || (km3 == km && lm3 > lm) ){   
      k = 3;
      l = l3;
      km = km3;
      lm = lm3;
   }
   if ( km4 > km || (km4 == km && lm4 > lm) ){
      k = 4;   
      l = l4;
      km = km4;
      lm = lm4;
   }
   /* Now, k is index\in(1,2,3,4) of the node which is common to the most 
      divided edges; l is index of the neighbour node on a divided edge which  
      is common to the most (among the neighbours) divided edges.             */
   if (k == 1) {
      if (l==3)
         reorder(*n1,*n3,*n2,*n4,*f1,*f3,*f2,*f4,*n13,*n12,*n14,*n23,*n34,*n24,
                 n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==4)
         reorder(*n1,*n4,*n2,*n3,*f1,*f4,*f2,*f3,*n14,*n12,*n13,*n24,*n34,*n23,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
   }
   else if (k == 2) {
      if (l==1)
         reorder(*n2,*n1,*n3,*n4,*f2,*f1,*f3,*f4,*n12,*n23,*n24,*n13,*n14,*n34,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==3)
         reorder(*n2,*n3,*n1,*n4,*f2,*f3,*f1,*f4,*n23,*n12,*n24,*n13,*n34,*n14,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==4)
         reorder(*n2,*n4,*n1,*n3,*f2,*f4,*f1,*f3,*n24,*n12,*n23,*n14,*n34,*n13,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
   }
   else if (k == 3) {
      if (l==1)
         reorder(*n3,*n1,*n2,*n4,*f3,*f1,*f2,*f4,*n13,*n23,*n34,*n12,*n14,*n24,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==2)
         reorder(*n3,*n2,*n1,*n4,*f3,*f2,*f1,*f4,*n23,*n13,*n34,*n12,*n24,*n14,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==4)
         reorder(*n3,*n4,*n1,*n2,*f3,*f4,*f1,*f2,*n34,*n13,*n23,*n14,*n24,*n12,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
   }
   else if (k == 4) {   
      if (l==1)
         reorder(*n4,*n1,*n2,*n3,*f4,*f1,*f2,*f3,*n14,*n24,*n34,*n12,*n13,*n23,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==2)
         reorder(*n4,*n2,*n1,*n3,*f4,*f2,*f1,*f3,*n24,*n14,*n34,*n12,*n23,*n13,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
      if (l==3)
         reorder(*n4,*n3,*n1,*n2,*f4,*f3,*f1,*f2,*n34,*n14,*n24,*n13,*n23,*n12,
                  n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
   }
   if (*n13 == NULL && *n14 != NULL) 
      reorder(*n1,*n2,*n4,*n3,*f1,*f2,*f4,*f3,*n12,*n14,*n13,*n24,*n23,*n34,
               n1, n2, n3, n4, f1, f2, f3, f4, n12, n13, n14, n23, n24, n34);
   return(0);
}

 void treat_1face(f,f0,n11,n22,n33,faceNumber,pface,pfaceil,pfacebl,pnflink)
 FACE *f, **f0, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33;
 NFLINK **pnflink;
 INT *faceNumber;
 {
    if (f->sons[0] == NULL){
    
       if (NOT_FF(f))
          *pfacebl = (*pfacebl)->succ = *pface;
       else
          *pfaceil = (*pfaceil)->succ = *pface;
       f->sons[0] = *pface;
       fmake((*pface)++,(FACE*)NULL,FTYPE(f),0,++(*faceNumber),f);
       if (IS_FF(f)){
          new_Fneighbour(n11,f->sons[0],pnflink);
          new_Fneighbour(n22,f->sons[0],pnflink);
          new_Fneighbour(n33,f->sons[0],pnflink);
       }
                
    }   
    *f0 = f->sons[0];
 }
  
 void Treat_1face(f,f0,n11,n22,n33,pnflink)   /*  for double refinement only  */
 FACE *f, **f0;
 NODE *n11, *n22, *n33;
 NFLINK **pnflink;
 {
    if (f->sons[0] == NULL){
       f->sons[0] = f;       /*  single refined face (temporary indication)  */
       if (IS_FF(f)){
          new_Fneighbour(n11,f->sons[0],pnflink);
          new_Fneighbour(n22,f->sons[0],pnflink);
          new_Fneighbour(n33,f->sons[0],pnflink);
       }
    }   
    *f0 = f->sons[0];
 }
 
 void treat_2face(f,f1,f2,n12,n11,n22,n33,faceNumber,pface,pfaceil,pfacebl,
                                                                  plink,pnflink)
 FACE *f, **f1, **f2, **pface, **pfaceil, **pfacebl;     
 NODE *n12, *n11, *n22, *n33;
 LINK **plink;
 NFLINK **pnflink;
 INT *faceNumber;
 {
    if (f->sons[0] == NULL){
       if (NOT_FF(f))
          *pfacebl = ( (*pfacebl)->succ = *pface ) + 1;
       else
          *pfaceil = ( (*pfaceil)->succ = *pface ) + 1;
       f->sons[0] = *pface;
       fmake(*pface,*pface + 1,FTYPE(f),0,++(*faceNumber),f);
       f->sons[1] = ++(*pface);
       fmake((*pface)++,(FACE*)NULL,FTYPE(f),0,++(*faceNumber),f);
       new_neighbour(n12,n33,0,IS_BF(f),plink);
       new_neighbour(n33,n12,0,IS_BF(f),plink);
       if (IS_FF(f)){
          new_2Fneighbours(n12,f->sons[0],f->sons[1],pnflink);
          new_2Fneighbours(n33,f->sons[0],f->sons[1],pnflink);
          if (n11->index < n22->index){
             new_Fneighbour(n11,f->sons[0],pnflink);
             new_Fneighbour(n22,f->sons[1],pnflink);
          }
          else{
             new_Fneighbour(n11,f->sons[1],pnflink);
             new_Fneighbour(n22,f->sons[0],pnflink);
          }
       }
    }
    if (n11->index < n22->index){   
       *f1 = f->sons[0];
       *f2 = f->sons[1];
    }
    else{
       *f1 = f->sons[1];
       *f2 = f->sons[0];
    }    
 }
 
void treat_3face(f,f1,f2,f3,n12,n13,n11,n22,n33,faceNumber,pface,pfaceil,
                                                          pfacebl,plink,pnflink)
FACE *f, **f1, **f2, **f3, **pface, **pfaceil, **pfacebl; 
NODE *n12, *n13, *n11, *n22, *n33;                   /* n2->index < n3->index */
LINK **plink;
NFLINK **pnflink;
INT *faceNumber;
{
   INT i;

   if (f->sons[0] == NULL){
      if (NOT_FF(f))
         *pfacebl = ( (*pfacebl)->succ = *pface ) + 2;
      else
         *pfaceil = ( (*pfaceil)->succ = *pface ) + 2;
      for (i=0; i<3; i++){
         f->sons[i] = *pface;
         fmake(*pface,*pface + 1,FTYPE(f),0,++(*faceNumber),f);
         (*pface)++;
      }
      ((*pface)-1)->succ = NULL;
      new_2neighbours(n12,n13,n33,IS_BF(f),plink);
      new_neighbour(n13,n12,0,IS_BF(f),plink);
      new_neighbour(n33,n12,0,IS_BF(f),plink);
      if (IS_FF(f)){
         new_3Fneighbours(n12,f->sons[0],f->sons[1],f->sons[2],pnflink);
         new_2Fneighbours(n13,f->sons[0],f->sons[2],pnflink);
         new_Fneighbour(n11,f->sons[0],pnflink);
         new_Fneighbour(n22,f->sons[1],pnflink);
         new_2Fneighbours(n33,f->sons[1],f->sons[2],pnflink);
      }
   }   
   *f1 = f->sons[0];
   *f2 = f->sons[1];
   *f3 = f->sons[2];
}
 
 void pyramid(theElement,k1,k2,n1,n2,n3,n4,n5,f1,f2,f3,f4,f5,f52,f54,
                                               pnflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, *f5, *f52, *f54;
 NODE *n1, *n2, *n3, *n4, *n5;
 NFLINK **pnflink;
 INT k1, k2;
 {
    new_3SFneighbours(n1,f5,f2,f3,pnflink); /* Faces lying on the surface of  */
    new_3SFneighbours(n3,f5,f1,f4,pnflink); /* the pyramid are not stored as  */
    new_3SFneighbours(n5,f5,f52,f54,pnflink);/*neighbours of their nodes here.*/
    new_Fneighbour(n2,f5,pnflink);              
    new_Fneighbour(n4,f5,pnflink);              
    elem2(theElement,k1,1,n1,n2,n3,n5,f2,f5,f1,f52,pelemil,pelembl,pelement);
    elem2(theElement,k2,1,n3,n4,n1,n5,f4,f5,f3,f54,pelemil,pelembl,pelement);
 }
 
 void finer0(k,theElement,n11,n22,n33,n44,f1,f2,f3,f4,pnflink,pfnlink,pflink)
 ELEMENT *theElement;
 FACE *f1, *f2, *f3, *f4;
 NODE *n11, *n22, *n33, *n44;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT k;
 {
    FACE *f11, *f22, *f33, *f44; 
 
    if (k) {
       Treat_1face(f1,&f11,n22,n33,n44,pnflink);
       Treat_1face(f2,&f22,n11,n33,n44,pnflink);
       Treat_1face(f3,&f33,n11,n22,n44,pnflink);
       Treat_1face(f4,&f44,n11,n22,n33,pnflink);
       if (IS_FF(f11)) new_Fneighbour(n11,f11,pnflink);
       if (IS_FF(f22)) new_Fneighbour(n22,f22,pnflink);
       if (IS_FF(f33)) new_Fneighbour(n33,f33,pnflink);
       if (IS_FF(f44)) new_Fneighbour(n44,f44,pnflink);       
       nodesToFaces(pfnlink,theElement,NMASK,FMASK);
       facesToFaces(pflink,theElement,FMASK);
    }
 } 
   
void finer1(k,theElement,n11,n22,n33,n44,n12,f1,f2,f3,f4,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n44, *n12;
LINK **plink;
NFLINK **pnflink;
FNLINK **pfnlink;
FLINK **pflink;
INT k, *faceNumber;
{
   FACE *f11, *f22, *f31, *f32, *f41, *f42, *fi; 

   if (k) {
      Treat_1face(f1,&f11,n22,n33,n44,pnflink);
      Treat_1face(f2,&f22,n11,n33,n44,pnflink);
   }
   else{
      treat_1face(f1,&f11,n22,n33,n44,faceNumber,pface,pfaceil,pfacebl,pnflink);
      treat_1face(f2,&f22,n11,n33,n44,faceNumber,pface,pfaceil,pfacebl,pnflink);
   }
   treat_2face(f3,&f31,&f32,n12,n11,n22,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
   treat_2face(f4,&f41,&f42,n12,n11,n22,n33,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
   *pfaceil = (*pfaceil)->succ = fi = *pface;
   fmake((*pface)++,(FACE*)NULL,0,OR_SET,++(*faceNumber),(FACE*)NULL);
   new_Fneighbour(n11,fi,pnflink);
   new_Fneighbour(n22,fi,pnflink);
   new_3SFneighbours(n12,fi,f11,f22,pnflink);
   new_3SFneighbours(n33,fi,f31,f32,pnflink);
   new_3SFneighbours(n44,fi,f41,f42,pnflink);
   elem2(theElement,0,1,n11,n12,n33,n44,fi,f22,f31,f41,pelemil,pelembl,pelement);
   elem2(theElement,1,1,n12,n22,n33,n44,f11,fi,f32,f42,pelemil,pelembl,pelement);
   nodes_and_facesToFaces(theElement,2,pfnlink,pflink);
}
 
 void finer2A(theElement,n11,n22,n33,n44,n12,n34,f1,f2,f3,f4,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n34;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f13, *f14, *f23, *f24, *f31, *f32, *f41, *f42, *fi[5]; 
 
    treat_2face(f1,&f13,&f14,n34,n33,n44,n22,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f2,&f23,&f24,n34,n33,n44,n11,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f3,&f31,&f32,n12,n11,n22,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f4,&f41,&f42,n12,n11,n22,n33,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    k_inner_faces(4,fi,faceNumber,pface,pfaceil);
    new_3Fneighbours(n11,fi[1],fi[3],fi[4],pnflink);
    new_3Fneighbours(n22,fi[2],fi[3],fi[4],pnflink);
    new_3Fneighbours(n33,fi[3],fi[1],fi[2],pnflink);
    new_3Fneighbours(n44,fi[4],fi[1],fi[2],pnflink);
    new_4Fneighbours(n12,fi[1],fi[2],f13,f14,pnflink);
    new_4Fneighbours(n12,fi[3],fi[4],f23,f24,pnflink);
    new_4Fneighbours(n34,fi[1],fi[2],f31,f32,pnflink);
    new_4Fneighbours(n34,fi[3],fi[4],f41,f42,pnflink);
    edge(n12,n34,0,0,plink);
    elem2(theElement,0,1,n11,n12,n33,n34,fi[3],f23,fi[1],f41,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,1,1,n11,n12,n34,n44,fi[4],f24,f31,fi[1],pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,2,1,n12,n22,n33,n34,f13,fi[3],fi[2],f42,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,3,1,n12,n22,n34,n44,f14,fi[4],f32,fi[2],pelemil,pelembl,
                                                                      pelement);
    nodes_and_facesToFaces(theElement,4,pfnlink,pflink);
 }
  
 void finer2B(k,theElement,n11,n22,n33,n44,n12,n13,f1,f2,f3,f4,faceNumber,pface,
          pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT k, *faceNumber;
 {
    FACE *f11, *f21, *f23, *f31, *f32, *f41, *f42, *f43, *fi[3]; 
    NODE *n;
 
    if (n33->index < n22->index)
       reorder( n11, n33, n22, n44, f1, f3, f2, f4, n13, n12, n, n, n, n,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n,&n,&n,&n);
    if (k) 
      Treat_1face(f1,&f11,n22,n33,n44,pnflink);
    else
      treat_1face(f1,&f11,n22,n33,n44,faceNumber,pface,pfaceil,pfacebl,pnflink);
    treat_2face(f2,&f21,&f23,n13,n11,n33,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f3,&f31,&f32,n12,n11,n22,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_3face(f4,&f41,&f42,&f43,n12,n13,n11,n22,n33,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    k_inner_faces(2,fi,faceNumber,pface,pfaceil);
    new_Fneighbour(n11,fi[1],pnflink);
    new_2SFneighbours(n12,fi[1],f21,pnflink);
    new_2SFneighbours(n13,fi[1],f31,pnflink);
    new_2SFneighbours(n44,fi[1],f41,pnflink);
    elem2(theElement,0,1,n11,n12,n13,n44,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
    pyramid(theElement,1,2,n12,n22,n33,n13,n44,f32,f11,f23,fi[1],fi[2],f42,f43,
                                              pnflink,pelemil,pelembl,pelement);
    nodes_and_facesToFaces(theElement,3,pfnlink,pflink);
 }
  
 void finer3A(theElement,n11,n22,n33,n44,n12,n13,n24,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n24;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f12, *f14, *f21, *f23, *f31, *f32, *f34, *f41, *f42, *f43, *fi[6]; 
 
    treat_2face(f1,&f12,&f14,n24,n22,n44,n33,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f2,&f21,&f23,n13,n11,n33,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    edge(n13,n24,0,0,plink);
    k_inner_faces(5,fi,faceNumber,pface,pfaceil);
    if (n44->index > n11->index)
       treat_3face(f3,&f32,&f31,&f34,n12,n24,n22,n11,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    else
       treat_3face(f3,&f32,&f34,&f31,n24,n12,n22,n44,n11,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    if (n33->index > n22->index)
       treat_3face(f4,&f41,&f42,&f43,n12,n13,n11,n22,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    else
       treat_3face(f4,&f41,&f43,&f42,n13,n12,n11,n33,n22,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    new_Fneighbour(n12,fi[2],pnflink);
    new_2Fneighbours(n33,fi[3],fi[4],pnflink);
    new_2Fneighbours(n44,fi[3],fi[4],pnflink);
    new_4Fneighbours(n13,fi[2],fi[3],fi[4],f14,pnflink);
    new_4Fneighbours(n24,fi[2],fi[3],fi[4],f23,pnflink);
    if (n44->index > n11->index)
       pyramid(theElement,0,1,n44,n11,n12,n24,n13,f21,f41,fi[2],fi[4],fi[1],f31,
                                          f34,pnflink,pelemil,pelembl,pelement);
    else
       pyramid(theElement,0,1,n11,n12,n24,n44,n13,f41,fi[2],fi[4],f21,fi[1],f31,
                                          f34,pnflink,pelemil,pelembl,pelement);
    if (n33->index > n22->index)
       pyramid(theElement,2,3,n12,n22,n33,n13,n24,f32,f12,fi[3],fi[2],fi[5],f42,
                                          f43,pnflink,pelemil,pelembl,pelement);
    else
       pyramid(theElement,2,3,n22,n33,n13,n12,n24,f12,fi[3],fi[2],f32,fi[5],f43,
                                          f42,pnflink,pelemil,pelembl,pelement);
    elem2(theElement,4,1,n24,n33,n13,n44,f23,fi[4],f14,fi[3],pelemil,pelembl,
                                                                      pelement);
    nodes_and_facesToFaces(theElement,5,pfnlink,pflink);
 }
 
 void finer3B(theElement,n11,n22,n33,n44,n12,n13,n23,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n23;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f12, *f13, *f21, *f23, *f31, *f32, *f41, *f42, *f43, *f40, *fi[4]; 
 
    treat_2face(f1,&f12,&f13,n23,n22,n33,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f2,&f21,&f23,n13,n11,n33,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_2face(f3,&f31,&f32,n12,n11,n22,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_4Face(f4,&f41,&f42,&f43,&f40,n12,n13,n23,n11,n22,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    k_inner_faces(3,fi,faceNumber,pface,pfaceil);
    new_Fneighbour(n11,fi[1],pnflink);
    new_Fneighbour(n22,fi[2],pnflink);
    new_Fneighbour(n33,fi[3],pnflink);
    new_7Fneighbours(n44,fi[1],fi[2],fi[3],f41,f42,f43,f40,pnflink);
    new_5Fneighbours(n12,fi[1],fi[2],fi[3],f12,f21,pnflink);
    new_5Fneighbours(n23,fi[1],fi[2],fi[3],f23,f32,pnflink);
    new_5Fneighbours(n13,fi[1],fi[2],fi[3],f13,f31,pnflink);
    elem2(theElement,0,1,n11,n12,n13,n44,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,1,1,n22,n23,n12,n44,fi[2],f32,f12,f42,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,2,1,n33,n13,n23,n44,fi[3],f13,f23,f43,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,3,1,n12,n23,n13,n44,fi[3],fi[1],fi[2],f40,pelemil,pelembl,
                                                                      pelement);
    nodes_and_facesToFaces(theElement,4,pfnlink,pflink);
 }
             
 void finer3C(k,theElement,n11,n22,n33,n44,n12,n13,n14,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT k, *faceNumber;
 {
    FACE *f11, *f21, *f23, *f24, *f31, *f32, *f34, *f41, *f42, *f43, *fi[4]; 
    NODE *n;
    
    if (n33->index > n22->index && n33->index > n44->index)
       reorder( n11, n22, n44, n33, f1, f2, f4, f3, n12, n14, n13, n, n, n,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n,&n,&n);
    else if (n22->index > n33->index && n22->index > n44->index)
       reorder( n11, n33, n44, n22, f1, f3, f4, f2, n13, n14, n12, n, n, n,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n,&n,&n);
    if (n22->index > n33->index)
       reorder( n11, n33, n22, n44, f1, f3, f2, f4, n13, n12, n14, n, n, n,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n,&n,&n);
    if (k)
      Treat_1face(f1,&f11,n22,n33,n44,pnflink);
    else
      treat_1face(f1,&f11,n22,n33,n44,faceNumber,pface,pfaceil,pfacebl,pnflink);
    treat_3face(f2,&f21,&f23,&f24,n13,n14,n11,n33,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f3,&f31,&f32,&f34,n12,n14,n11,n22,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f4,&f41,&f42,&f43,n12,n13,n11,n22,n33,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    k_inner_faces(3,fi,faceNumber,pface,pfaceil);  
    new_Fneighbour(n11,fi[1],pnflink);
    new_4Fneighbours(n12,fi[1],fi[2],f21,f24,pnflink);
    new_4Fneighbours(n13,fi[1],fi[2],f31,f34,pnflink);
    new_3SFneighbours(n14,fi[1],fi[2],f41,pnflink);
    new_2Fneighbours(n44,fi[1],fi[2],pnflink);
    elem2(theElement,0,1,n11,n12,n13,n14,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,1,1,n12,n13,n14,n44,f24,f34,fi[2],fi[1],pelemil,pelembl,
                                                                      pelement);
    pyramid(theElement,2,3,n12,n22,n33,n13,n44,f32,f11,f23,fi[2],fi[3],f42,f43,
                                              pnflink,pelemil,pelembl,pelement);
    nodes_and_facesToFaces(theElement,4,pfnlink,pflink);
 }  
                   
 void finer4A(theElement,n11,n22,n33,n44,n12,n13,n24,n34,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n24, *n34;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f12, *f13, *f14, *f21, *f23, *f24, *f31, *f32, *f34, 
         *f41, *f42, *f43, *fi[7]; 
    NODE *n;
    FLOAT c1, c2;
    
    if (n33->index < n22->index)
       reorder( n11, n33, n22, n44, f1, f3, f2, f4, n13, n12, n, n, n34, n24,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n,&n,&n24,&n34);
    if (n44->index < n11->index)
       reorder( n44, n22, n33, n11, f4, f2, f3, f1, n24, n34, n, n, n12, n13,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n,&n,&n24,&n34);
    treat_3face(f1,&f14,&f12,&f13,n24,n34,n44,n22,n33,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f2,&f23,&f21,&f24,n13,n34,n33,n11,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f3,&f32,&f31,&f34,n12,n24,n22,n11,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f4,&f41,&f42,&f43,n12,n13,n11,n22,n33,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    k_inner_faces(6,fi,faceNumber,pface,pfaceil);           
    new_Fneighbour(n11,fi[4],pnflink);
    new_Fneighbour(n22,fi[3],pnflink);   
    new_2SFneighbours(n33,fi[3],f32,pnflink);
    new_2SFneighbours(n44,fi[4],f41,pnflink);
    new_5Fneighbours(n12,fi[3],fi[4],fi[6],f12,f21,pnflink);
    new_3SFneighbours(n13,fi[4],fi[5],f31,pnflink);
    new_3SFneighbours(n24,fi[3],fi[6],f42,pnflink);
    new_Fneighbour(n34,fi[5],pnflink);
    c1 = crit1(n12,n13,n34) + crit1(n34,n24,n12);
    c2 = crit1(n24,n12,n13) + crit1(n13,n34,n24);
    if (c1 < c2){     /*   (n12,n34) is the diagonal   */
       pyramid(theElement,0,1,n12,n13,n34,n24,n33,f43,f23,f13,fi[3],fi[1],
                                  fi[5],fi[6],pnflink,pelemil,pelembl,pelement);
       pyramid(theElement,2,3,n12,n13,n34,n24,n44,fi[4],f24,f14,f34,fi[2],
                                  fi[5],fi[6],pnflink,pelemil,pelembl,pelement);
       edge(n12,n34,0,0,plink);
       new_Fneighbour(n12,fi[5],pnflink);
       new_Fneighbour(n34,fi[6],pnflink);
    }
    else{             /*   (n13,n24) is the diagonal   */
       pyramid(theElement,0,1,n13,n34,n24,n12,n33,f23,f13,fi[3],f43,fi[1],
                                  fi[5],fi[6],pnflink,pelemil,pelembl,pelement);
       pyramid(theElement,2,3,n13,n34,n24,n12,n44,f24,f14,f34,fi[4],fi[2],
                                  fi[5],fi[6],pnflink,pelemil,pelembl,pelement);
       edge(n13,n24,0,0,plink);
       new_Fneighbour(n13,fi[6],pnflink);
       new_Fneighbour(n24,fi[5],pnflink);
    }
    elem2(theElement,4,1,n11,n12,n13,n44,fi[4],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,5,1,n12,n22,n33,n24,f12,fi[3],f32,f42,pelemil,pelembl,
                                                                      pelement);
    nodes_and_facesToFaces(theElement,6,pfnlink,pflink);
 }  
      
 void finer4B(theElement,n11,n22,n33,n44,n12,n13,n14,n23,f1,f2,f3,f4,faceNumber,
    pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f12, *f13, *f21, *f23, *f24, *f31, *f32, *f34, 
         *f41, *f42, *f43, *f40, *fi[7]; 
    
    treat_2face(f1,&f12,&f13,n23,n22,n33,n44,faceNumber,pface,pfaceil,pfacebl,
                                                                 plink,pnflink);
    treat_4Face(f4,&f41,&f42,&f43,&f40,n12,n13,n23,n11,n22,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    if (n44->index > n22->index && n44->index > n33->index){
       treat_3face(f2,&f21,&f23,&f24,n13,n14,n11,n33,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
       treat_3face(f3,&f31,&f32,&f34,n12,n14,n11,n22,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
       k_inner_faces(4,fi,faceNumber,pface,pfaceil);  
       new_Fneighbour(n11,fi[1],pnflink);
       new_Fneighbour(n22,fi[2],pnflink);
       new_Fneighbour(n33,fi[3],pnflink);
       new_7Fneighbours(n44,fi[1],fi[2],fi[3],fi[4],f42,f43,f40,pnflink);
       new_7Fneighbours(n12,fi[1],fi[2],fi[3],fi[4],f12,f21,f24,pnflink);
       new_7Fneighbours(n13,fi[1],fi[2],fi[3],fi[4],f13,f31,f34,pnflink);
       new_5Fneighbours(n23,fi[2],fi[3],fi[4],f23,f32,pnflink);
       new_3SFneighbours(n14,fi[1],fi[4],f41,pnflink);
       elem2(theElement,0,1,n11,n12,n13,n14,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
       elem2(theElement,1,1,n12,n13,n14,n44,f24,f34,fi[4],fi[1],pelemil,pelembl,
                                                                      pelement);
       elem2(theElement,2,1,n12,n23,n13,n44,fi[3],fi[4],fi[2],f40,pelemil,
                                                              pelembl,pelement);
       elem2(theElement,3,1,n12,n22,n23,n44,f12,fi[2],f32,f42,pelemil,pelembl,
                                                                      pelement);
       elem2(theElement,4,1,n13,n23,n33,n44,f13,f23,fi[3],f43,pelemil,pelembl,
                                                                      pelement);
       nodes_and_facesToFaces(theElement,5,pfnlink,pflink);
    }
    else{
       edge(n14,n23,0,0,plink);
       k_inner_faces(6,fi,faceNumber,pface,pfaceil);
       if (n44->index > n22->index){
          treat_3face(f3,&f31,&f32,&f34,n12,n14,n11,n22,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
          pyramid(theElement,0,1,n12,n22,n44,n14,n23,f42,f12,fi[4],fi[2],fi[5],
                                      f32,f34,pnflink,pelemil,pelembl,pelement);
       }
       else{
          treat_3face(f3,&f31,&f34,&f32,n14,n12,n11,n44,n22,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
          pyramid(theElement,0,1,n14,n12,n22,n44,n23,fi[2],f42,f12,fi[4],fi[5],
                                      f32,f34,pnflink,pelemil,pelembl,pelement);
       }
       if (n44->index > n33->index){
          treat_3face(f2,&f21,&f23,&f24,n13,n14,n11,n33,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
          pyramid(theElement,2,3,n13,n33,n44,n14,n23,f43,f13,fi[4],fi[3],fi[6],
                                      f23,f24,pnflink,pelemil,pelembl,pelement);
       }
       else{
          treat_3face(f2,&f21,&f24,&f23,n14,n13,n11,n44,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
          pyramid(theElement,2,3,n14,n13,n33,n44,n23,fi[3],f43,f13,fi[4],fi[6],
                                      f23,f24,pnflink,pelemil,pelembl,pelement);
       } 
       new_Fneighbour(n11,fi[1],pnflink);
       new_Fneighbour(n44,fi[4],pnflink);    
       new_4Fneighbours(n12,fi[1],fi[2],fi[3],f21,pnflink);     
       new_4Fneighbours(n13,fi[1],fi[2],fi[3],f31,pnflink);     
       new_4Fneighbours(n23,fi[1],fi[2],fi[3],fi[4],pnflink);     
       new_3Fneighbours(n14,fi[1],fi[2],fi[3],pnflink);     
       new_3SFneighbours(n14,fi[4],f41,f40,pnflink);    
       elem2(theElement,4,1,n11,n12,n13,n14,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
       elem2(theElement,5,1,n12,n23,n13,n14,fi[3],fi[1],fi[2],f40,pelemil,
                                                              pelembl,pelement);
       nodes_and_facesToFaces(theElement,6,pfnlink,pflink);
    }
 }  
                   
 void finer5(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,f1,f2,f3,f4,
             faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                       pelemil,pelembl,pelement)
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24;
 LINK **plink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 FLINK **pflink;
 INT *faceNumber;
 {
    FACE *f12, *f13, *f14, *f21, *f23, *f24, *f31, *f32, *f34, *f30, 
         *f41, *f42, *f43, *f40, *fi[8]; 
    NODE *n;
    FLOAT c1, c2;
    
    if (n33->index > n44->index)
       reorder( n11, n22, n44, n33, f1, f2, f4, f3, n12, n14, n13, n24, n23, n,
               &n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n23,&n24,&n);
    treat_3face(f1,&f12,&f13,&f14,n23,n24,n22,n33,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_3face(f2,&f21,&f23,&f24,n13,n14,n11,n33,n44,faceNumber,pface,pfaceil,
                                                         pfacebl,plink,pnflink);
    treat_4Face(f3,&f31,&f32,&f34,&f30,n12,n14,n24,n11,n22,n44,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    treat_4Face(f4,&f41,&f42,&f43,&f40,n12,n13,n23,n11,n22,n33,faceNumber,pface,
                                                 pfaceil,pfacebl,plink,pnflink);
    k_inner_faces(7,fi,faceNumber,pface,pfaceil); 
    new_Fneighbour(n11,fi[1],pnflink);
    new_Fneighbour(n22,fi[2],pnflink);
    new_Fneighbour(n33,fi[3],pnflink);
    new_2SFneighbours(n44,fi[3],f43,pnflink);
    new_4Fneighbours(n12,fi[1],fi[2],f12,f21,pnflink);   
    new_3SFneighbours(n14,fi[1],fi[7],f41,pnflink);
    new_3SFneighbours(n24,fi[2],fi[7],f42,pnflink);
    new_5Fneighbours(n13,fi[1],fi[3],fi[6],f13,f31,pnflink);
    new_5Fneighbours(n23,fi[2],fi[3],fi[6],f23,f32,pnflink);
    c1 = crit1(n13,n24,n14) + crit1(n13,n23,n24);
    c2 = crit1(n14,n13,n23) + crit1(n14,n23,n24);
    if (c1 < c2){     /*   (n13,n24) is the diagonal   */
       pyramid(theElement,0,1,n13,n23,n24,n14,n12,f40,fi[2],f30,fi[1],fi[4],
                                  fi[6],fi[7],pnflink,pelemil,pelembl,pelement);
       pyramid(theElement,2,3,n13,n23,n24,n14,n44,fi[3],f14,f34,f24,fi[5],
                                  fi[6],fi[7],pnflink,pelemil,pelembl,pelement);
       edge(n13,n24,0,0,plink);
       new_Fneighbour(n13,fi[7],pnflink);
       new_Fneighbour(n24,fi[6],pnflink);
    }
    else{             /*   (n14,n23) is the diagonal   */
       pyramid(theElement,0,1,n14,n13,n23,n24,n12,fi[1],f40,fi[2],f30,fi[4],
                                  fi[6],fi[7],pnflink,pelemil,pelembl,pelement);
       pyramid(theElement,2,3,n14,n13,n23,n24,n44,f24,fi[3],f14,f34,fi[5],
                                  fi[6],fi[7],pnflink,pelemil,pelembl,pelement);
       edge(n14,n23,0,0,plink);
       new_Fneighbour(n14,fi[6],pnflink);
       new_Fneighbour(n23,fi[7],pnflink);
    }
    elem2(theElement,4,1,n11,n12,n13,n14,fi[1],f21,f31,f41,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,5,1,n12,n22,n23,n24,f12,fi[2],f32,f42,pelemil,pelembl,
                                                                      pelement);
    elem2(theElement,6,1,n13,n23,n33,n44,f13,f23,fi[3],f43,pelemil,pelembl,
                                                                      pelement);
    nodes_and_facesToFaces(theElement,7,pfnlink,pflink);
 }  
 
 void treat_faceEdge(mg,n12,n13,type,pnode,pnodeil,pnodebl,pvert,plink,
                                                                     nodeNumber)
 MULTIGRID *mg;                                 
 NODE *n12, *n13, **pnode, **pnodeil, **pnodebl;
 VERTEX **pvert;                               
 LINK **plink;                                 
 INT type, *nodeNumber;                             
 {
    NODE *m1, *m2;
    
    m1 = n12;
    m2 = n13;
    if (lower_edge(mg,&m1,&m2))
       nu_treat_edge(mg,m1,m2,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    else if(!get_link(n12,n13)){
       new_neighbour(n12,n13,0,type,plink);
       new_neighbour(n13,n12,0,type,plink);
    }
 }
  
void Treat_4face(mg,f,f1,f2,f3,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,plink,
  nodeNumber,faceNumber,pface,pfaceil,pfacebl) /* n11->index < n22->index ... */
MULTIGRID *mg;
FACE *f, **f1, **f2, **f3, **f0, **pface, **pfaceil, **pfacebl; 
NODE *n12, *n13, *n23, **pnode, **pnodeil, **pnodebl;
VERTEX **pvert;                               
LINK **plink;                                 
INT *nodeNumber, *faceNumber;
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
      treat_faceEdge(mg,n12,n13,IS_BF(f),pnode,pnodeil,pnodebl,pvert,plink,
                                                                    nodeNumber);
      treat_faceEdge(mg,n12,n23,IS_BF(f),pnode,pnodeil,pnodebl,pvert,plink,
                                                                    nodeNumber);
      treat_faceEdge(mg,n13,n23,IS_BF(f),pnode,pnodeil,pnodebl,pvert,plink,
                                                                    nodeNumber);
   }   
   *f1 = f->sons[1];
   *f2 = f->sons[2];
   *f3 = f->sons[3];
   *f0 = f->sons[0];
}
 
 void Treat_4Face(mg,f,f1,f2,f3,f0,n12,n13,n23,n11,n22,n33,pnode,pnodeil,
                pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl)
 MULTIGRID *mg;                                /* n11->index < n22->index ... */
 FACE *f, **f1, **f2, **f3, **f0, **pface, **pfaceil, **pfacebl; 
 NODE *n12, *n13, *n23, *n11, *n22, *n33, **pnode, **pnodeil, **pnodebl;
 VERTEX **pvert;                               
 LINK **plink;                                 
 INT *nodeNumber, *faceNumber;
 {
    if (n11->index < n22->index){
       if (n22->index < n33->index)        /*   n11 < n22 < n33   */
          Treat_4face(mg,f,f1,f2,f3,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
       else if (n11->index < n33->index)   /*   n11 < n33 < n22   */
          Treat_4face(mg,f,f1,f3,f2,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
       else                                /*   n33 < n11 < n22   */
          Treat_4face(mg,f,f3,f1,f2,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
    }
    else if (n33->index < n22->index)      /*   n33 < n22 < n11   */
       Treat_4face(mg,f,f3,f2,f1,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
    else if (n33->index < n11->index)      /*   n22 < n33 < n11   */
       Treat_4face(mg,f,f2,f3,f1,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
    else                                   /*   n22 < n11 < n33   */
       Treat_4face(mg,f,f2,f1,f3,f0,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,
                             plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
 }
                 
void Inner_elements(theElement,f11,f22,f33,f44,f10,f20,f30,f40,
      m1,m2,m3,m4,m5,m6,faceNumber,pface,pfaceil,plink,pelemil,pelembl,pelement)
FACE *f11,*f22,*f33,*f44,*f10,*f20,*f30,*f40,**pface,**pfaceil;
NODE *m1,*m2,*m3,*m4,*m5,*m6;
LINK **plink;
ELEMENT *theElement, **pelemil, **pelembl, **pelement; 
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
}

/* finer6 without neighbour relations.                                             */
void Finer6(mg,theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
            pnode,pnodeil,pnodebl,pvert,plink,nodeNumber,faceNumber,pface,
                                       pfaceil,pfacebl,pelemil,pelembl,pelement)
MULTIGRID *mg;
ELEMENT *theElement, **pelemil, **pelembl, **pelement;
FACE *f1, *f2, *f3, *f4, **pface, **pfaceil, **pfacebl;
NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34, **pnode, 
                                                           **pnodeil, **pnodebl;
VERTEX **pvert;                               
LINK **plink;                                 
INT *nodeNumber, *faceNumber;
{
 FACE *f11, *f12, *f13, *f14, *f10, *f21, *f22, *f23, *f24, *f20, 
      *f31, *f32, *f33, *f34, *f30, *f41, *f42, *f43, *f44, *f40;
 FLOAT c1, c2, c3;
 INT  i; 

 Treat_4Face(mg,f1,&f12,&f13,&f14,&f10,n23,n24,n34,n22,n33,n44,pnode,pnodeil,
              pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
 Treat_4Face(mg,f2,&f21,&f23,&f24,&f20,n13,n14,n34,n11,n33,n44,pnode,pnodeil,
              pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
 Treat_4Face(mg,f3,&f31,&f32,&f34,&f30,n12,n14,n24,n11,n22,n44,pnode,pnodeil,
              pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
 Treat_4Face(mg,f4,&f41,&f42,&f43,&f40,n12,n13,n23,n11,n22,n33,pnode,pnodeil,
              pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,pfacebl);
 f11 = *pface;
 f22 = *pface + 1;
 f33 = *pface + 2;
 f44 = *pface + 3;
 (*pfaceil)->succ = f11;
 *pfaceil = f44; 
 for (i=0; i<4; i++){
    fmake(*pface,*pface + 1,0,OR_SET,++(*faceNumber),(FACE*)NULL);
    (*pface)++;
 }
 ((*pface)-1)->succ = NULL;
 elem2(theElement,0,0,n11,n12,n13,n14,f11,f21,f31,f41,pelemil,pelembl,pelement);
 elem2(theElement,1,0,n22,n12,n23,n24,f22,f12,f32,f42,pelemil,pelembl,pelement);
 elem2(theElement,2,0,n33,n13,n23,n34,f33,f13,f23,f43,pelemil,pelembl,pelement);
 elem2(theElement,3,0,n44,n14,n24,n34,f44,f14,f24,f34,pelemil,pelembl,pelement);
 c1=crit1(n13,n24,n12)+crit1(n13,n24,n34)+crit1(n13,n24,n23)+crit1(n13,n24,n14);
 c2=crit1(n12,n34,n13)+crit1(n12,n34,n24)+crit1(n12,n34,n14)+crit1(n12,n34,n23);
 c3=crit1(n14,n23,n13)+crit1(n14,n23,n24)+crit1(n14,n23,n34)+crit1(n14,n23,n12);
 if ( c1<=c2 && c1<=c3 )  
    Inner_elements(theElement,f11,f44,f22,f33,f30,f20,f40,f10,n14,n23,n12,n34,
               n13,n24,faceNumber,pface,pfaceil,plink,pelemil,pelembl,pelement);
 else if ( c2<=c1 && c2<=c3 )
    Inner_elements(theElement,f11,f33,f44,f22,f20,f40,f30,f10,n13,n24,n14,n23,
               n12,n34,faceNumber,pface,pfaceil,plink,pelemil,pelembl,pelement);
 else
    Inner_elements(theElement,f44,f22,f33,f11,f10,f30,f20,f40,n24,n13,n34,n12,
               n14,n23,faceNumber,pface,pfaceil,plink,pelemil,pelembl,pelement);
}
 
 void nu_refine(k,theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,
      f1,f2,f3,f4,faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                       pelemil,pelembl,pelement)
 NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
 FACE *f1, *f2, *f3, *f4;
 ELEMENT *theElement, **pelemil, **pelembl, **pelement;
 FACE **pface, **pfaceil, **pfacebl;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT k, *faceNumber; 
 {
    permute(&n11,&n22,&n33,&n44,&f1,&f2,&f3,&f4,&n12,&n13,&n14,&n23,&n24,&n34);
    if (n12 == NULL)
       finer0(k,theElement,n11,n22,n33,n44,f1,f2,f3,f4,pnflink,pfnlink,pflink);
    else if (n13 == NULL){
       if (n34 == NULL)
          finer1(k,theElement,n11,n22,n33,n44,n12,f1,f2,f3,f4,faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
       else 
          finer2A(theElement,n11,n22,n33,n44,n12,n34,f1,f2,f3,f4,faceNumber,
                  pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,
                                                              pelembl,pelement);
    }
    else if (n14 == NULL){
       if (n23 == NULL){
          if (n24 == NULL){
             if (n34 == NULL)
                finer2B(k,theElement,n11,n22,n33,n44,n12,n13,f1,f2,f3,f4,
                        faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
             else 
                finer3A(theElement,n11,n33,n22,n44,n13,n12,n34,f1,f3,f2,f4,
                        faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
          }
          else if (n34==NULL)
             finer3A(theElement,n11,n22,n33,n44,n12,n13,n24,f1,f2,f3,f4,
                     faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
          else
             finer4A(theElement,n11,n22,n33,n44,n12,n13,n24,n34,f1,f2,f3,f4,
                     faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,
                                               pflink,pelemil,pelembl,pelement);
       }
       else
          finer3B(theElement,n11,n22,n33,n44,n12,n13,n23,f1,f2,f3,f4,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    }
    else if (n23 == NULL){
       if (n24 == NULL)
          finer3C(k,theElement,n11,n22,n33,n44,n12,n13,n14,f1,f2,f3,f4,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
       else 
          finer4B(theElement,n11,n22,n44,n33,n12,n14,n13,n24,f1,f2,f4,f3,
                  faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    }
    else if (n24 == NULL)
       finer4B(theElement,n11,n22,n33,n44,n12,n13,n14,n23,f1,f2,f3,f4,
               faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else if (n34 == NULL)
       finer5(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else
       finer6(theElement,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
 }  
 
 void treat_subface(ref,f,i1,j1,facs,fath,n11,n12,n13,ind,pface,pfaceil,pfacebl,
                                                             pnflink,faceNumber)
 FACE *f, *facs[4], *fath[4], **pface, **pfaceil, **pfacebl;
 NODE *n11, *n12, *n13;
 NFLINK **pnflink;
 INT ref, i1, j1, *faceNumber;
 INT ind;
 {
    if (facs[i1]){
       f->sons[j1] = facs[i1];
       if (ref == 3 && !facs[i1]->sons[0]) facs[i1]->sons[0] = facs[i1];
    }
    else if (ref == 3 && 
               (!get_link(n11,n12) || !get_link(n12,n13) || !get_link(n13,n11)))
       f->sons[j1] = fath[i1]; 
    else{
       if (NOT_FF(fath[i1]))
          *pfacebl = (*pfacebl)->succ = *pface;
       else
          *pfaceil = (*pfaceil)->succ = *pface;
       f->sons[j1] = *pface;
       fmake(f->sons[j1],++(*pface),FTYPE(fath[i1]),0,++(*faceNumber),fath[i1]);
       SET_ORIENT(f->sons[j1],norm_vect_or(n11,n12,n13,ind));
       fath[i1]->sons[0] = f->sons[j1];
       if (ref == 2 && IS_FF(f->sons[j1])){
          new_Fneighbour(n11,f->sons[j1],pnflink);
          new_Fneighbour(n12,f->sons[j1],pnflink);
          new_Fneighbour(n13,f->sons[j1],pnflink);
       }
    }
 }

void copy_faces(mg,pg,ref,fac,ind,ff,f,n1,n2,n3,n11,n22,n33,n12,n13,n23,pnode,
pnodeil,pnodebl,pvert,plink,pface,pfaceil,pfacebl,pnflink,nodeNumber,faceNumber)
MULTIGRID *mg;
GRID *pg;
NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23, **pnode, **pnodeil, 
                                                                      **pnodebl;
FACE *fac[4], *ff, **f, **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
NFLINK **pnflink;
INT ref, *faceNumber, *nodeNumber;
INT ind;
{
   FACE *facs[4], *fath[4];
   INT i, i1, i2, i3, j1, j2, j3;    
   
   *f = ff;
   (*f)->type = FTYPE(fac[0]);
   SET_ORIENT(*f,ind);
   (*f)->sons[1] = NULL;
   (*f)->sons[2] = NULL;
   (*f)->sons[3] = NULL;
   if (fac[3] == NULL){
      if(FACE_ON_NEW_LEVEL(fac[0],mg)){ 
         (*f)->sons[0] = fac[0];
         (*f)->sons[1] = fac[1];
         (*f)->sons[2] = fac[2];
         (*f)->father = *f; 
      }
      else{
         (*f)->sons[0] = NULL;
         (*f)->father = NULL;
      }
   }
   else{  /* fac[3] != NULL */
      if(FACE_ON_NEW_LEVEL(fac[0],mg) || fac[0]->index < 0){
         (*f)->sons[0] = fac[0];
         (*f)->sons[1] = fac[1];
         (*f)->sons[2] = fac[2];
         (*f)->sons[3] = fac[3];
         if (ref == 3)
            for (i = 0; i < 4; i++)
               if (!fac[i]->sons[0]) fac[i]->sons[0] = fac[i];
      }
      else{  /* the corresponding face of the refined "macroelement"  */   
         for (i = 0; i < 4; i++){         /*  consists of 4 subfaces  */
            facs[i] = fac[i];
            while (facs[i]->sons[0] && !facs[i]->sons[1]) 
               facs[i] = facs[i]->sons[0];
            if (!facs[i]->sons[0] && !FACE_ON_NEW_LEVEL(facs[i],mg)){
               fath[i] = facs[i];
               facs[i] = NULL;
            }
         }
         if (!(facs[0] && facs[1] && facs[2] && facs[3])){
            treat_faceEdge(mg,n12,n13,IS_BF(*f),pnode,pnodeil,pnodebl,pvert,
                                                              plink,nodeNumber);
            treat_faceEdge(mg,n12,n23,IS_BF(*f),pnode,pnodeil,pnodebl,pvert,
                                                              plink,nodeNumber);
            treat_faceEdge(mg,n13,n23,IS_BF(*f),pnode,pnodeil,pnodebl,pvert,
                                                              plink,nodeNumber);
         }
         while (pg->maxFaceIndex < fac[0]->index) pg = pg->finer;
         while (n1->index < pg->minNodeIndex) n1 = n1->son;
         while (n2->index < pg->minNodeIndex) n2 = n2->son;
         while (n3->index < pg->minNodeIndex) n3 = n3->son;
         find_order(n1,n2,n3,&i1,&i2,&i3);    /* fac[ik] lies at nk    */
         find_order(n11,n22,n33,&j1,&j2,&j3); /* nkk is the jk-th node */
         treat_subface(ref,*f,i1,j1,facs,fath,n11,n12,n13,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         treat_subface(ref,*f,i2,j2,facs,fath,n22,n23,n12,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         treat_subface(ref,*f,i3,j3,facs,fath,n33,n13,n23,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         treat_subface(ref,*f, 0, 0,facs,fath,n12,n23,n13,ind,pface,pfaceil,
                                                    pfacebl,pnflink,faceNumber);
         if (ref == 2){
            fac[0]->sons[3] = (*f)->sons[0]; 
            fac[1]->sons[3] = (*f)->sons[1]; 
            fac[2]->sons[3] = (*f)->sons[2]; 
            fac[3]->sons[3] = (*f)->sons[3]; 
         }
      }
      (*f)->father = *f; 
   }
}
 
 void fcopyS2(i1,i2,j1,j2,f,fac)
 INT i1, i2, j1, j2;
 FACE *f, *fac[4];
 {
    fac[i1]->sons[0] = f->sons[j1];
    fac[i2]->sons[0] = f->sons[j2];
    f->sons[j1]->father = fac[i1];
    f->sons[j2]->father = fac[i2];
 }
 
 void fcopy2(n1,n2,n11,n22,f,fac)
 NODE *n1, *n2, *n11, *n22;
 FACE *f, *fac[4];
 {
    if (n1->index < n2->index && n11->index < n22->index ||
        n1->index > n2->index && n11->index > n22->index )
       fcopyS2(0,1,0,1,f,fac);
    else
       fcopyS2(0,1,1,0,f,fac);
    fac[0]->sons[3] = f->sons[0];
    fac[1]->sons[3] = f->sons[1];   
 }
 
 void fcopy3(n2,n3,n22,n33,f,fac)
 NODE *n2, *n3, *n22, *n33;
 FACE *f, *fac[4];
 {
    fac[0]->sons[0] = f->sons[0];
    f->sons[0]->father = fac[0];
    if (n2->index < n3->index && n22->index < n33->index ||
        n2->index > n3->index && n22->index > n33->index )
       fcopyS2(1,2,1,2,f,fac);
    else
       fcopyS2(1,2,2,1,f,fac);
    fac[0]->sons[3] = f->sons[0];
    fac[1]->sons[3] = f->sons[1];   
    fac[2]->sons[3] = f->sons[2];   
 }
 
 void recopy_faces(pg,n1,n2,n3,n11,n22,n33,n12,n13,n23,f,fac)
 GRID *pg;
 NODE *n1, *n2, *n3, *n11, *n22, *n33, *n12, *n13, *n23;
 FACE *f, *fac[4];
 {
    INT i, j=0, k=0;
    
    if (f->father == NULL){
       set_norm_vect_or(n11,n22,n33,n12,n13,n23,f);
       while(fac[++j]);
       while(++k < 4 && f->sons[k]);
       if (j < k){
          for (i = 0; i < k; i++){
             fac[0]->sons[i] = f->sons[i];
             f->sons[i]->father = fac[0];
          }
          for (i = 1; i < j; i++)
             fac[i]->sons[0] = f->sons[0];
       }
       else   /*  j == k  */
          if (j < 2){ 
             fac[0]->sons[0] = f->sons[0];
             f->sons[0]->father = fac[0];
          }
          else{  
             while(pg->maxFaceIndex < fac[0]->index) pg = pg->finer;
             while(n1->index < pg->minNodeIndex) n1 = n1->son;
             while(n2->index < pg->minNodeIndex) n2 = n2->son;
             while(n3->index < pg->minNodeIndex) n3 = n3->son;
             if (j < 3)
                if (n12) 
                   fcopy2(n1,n2,n11,n22,f,fac);
                else if (n13)
                   fcopy2(n1,n3,n11,n33,f,fac);
                else 
                   fcopy2(n2,n3,n22,n33,f,fac);
             else if (j < 4)
                if (!n12)
                   fcopy3(n1,n2,n11,n22,f,fac);
                else if (!n13)
                   fcopy3(n1,n3,n11,n33,f,fac);
                else 
                   fcopy3(n2,n3,n22,n33,f,fac);
          }
    }
 }
 
 void set_ind1(n11,n22,n33,n12,n13,n23,ff,f)
 NODE *n11, *n22, *n33, *n12, *n13, *n23;
 FACE *ff, *f;
 {
    UNS type;
    
    type = FTYPE(f);
    FTYPE(f) = FTYPE(ff);
    set_norm_vect_or(n11,n22,n33,n12,n13,n23,f); 
    FTYPE(f) = type;
 }
 
/* refinement of an element with status=0 */
void refine1(mg,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,
          pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
MULTIGRID *mg;
ELEMENT *theElem, **pelemil, **pelembl, **pelement;
NODE **pnode, **pnodeil, **pnodebl;
FACE **pface, **pfaceil, **pfacebl;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *nodeNumber, *faceNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
    FACE *f1, *f2, *f3, *f4;
    
    if(FACE_ON_NEW_LEVEL(f1=top_face(theElem->f[0]),mg)) f1 = f1->father;
    if(FACE_ON_NEW_LEVEL(f2=top_face(theElem->f[1]),mg)) f2 = f2->father;
    if(FACE_ON_NEW_LEVEL(f3=top_face(theElem->f[2]),mg)) f3 = f3->father;
    if(FACE_ON_NEW_LEVEL(f4=top_face(theElem->f[3]),mg)) f4 = f4->father;
    treat_6edges(mg,theElem->n[0],theElem->n[1],theElem->n[2],theElem->n[3],
                 &n11,&n22,&n33,&n44,&n12,&n13,&n14,&n23,&n24,&n34,
                 pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    if (theElem->eflag>1)
       finer6(theElem,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    else
       nu_refine(0,theElem,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
                 faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    set_ind1(n11,n22,n33,n12,n13,n23,theElem->f[3],f4); 
    set_ind1(n11,n22,n44,n12,n14,n24,theElem->f[2],f3); 
    set_ind1(n11,n33,n44,n13,n14,n34,theElem->f[1],f2);  
    set_ind1(n22,n33,n44,n23,n24,n34,theElem->f[0],f1);    
    theElem->eflag = 0; 
 }
 
 /* single refinement of elements with status=1 */
 void refine2(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
              pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,
                                                          nodeNumber,faceNumber)
 MULTIGRID *mg;
 GRID *theGrid;
 ELEMENT *theElem, **pelemil, **pelembl, **pelement;
 NODE **pnode, **pnodeil, **pnodebl;
 FACE **pface, **pfaceil, **pfacebl;
 VERTEX **pvert;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT  *nodeNumber, *faceNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34,
         *no1, *no2, *no3, *no4;
    FACE ff1, ff2, ff3, ff4, *f1, *f2, *f3, *f4, *fac[4][4];
    ELEMENT *el[8];
    INT ind[4];
    INT i;
    
    data2_from_forefather(theGrid,theElem,&no1,&no2,&no3,&no4,el,fac,ind);
    treat_6edges(mg,no1,no2,no3,no4,&n11,&n22,&n33,&n44,&n12,&n13,&n14,
                   &n23,&n24,&n34,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    copy_faces(mg,theGrid,2,fac[0],ind[0],&ff1,&f1,no2,no3,no4,
                n22,n33,n44,n23,n24,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,2,fac[1],ind[1],&ff2,&f2,no1,no3,no4,
                n11,n33,n44,n13,n14,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,2,fac[2],ind[2],&ff3,&f3,no1,no2,no4,
                n11,n22,n44,n12,n14,n24,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,2,fac[3],ind[3],&ff4,&f4,no1,no2,no3,
                n11,n22,n33,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    nu_refine(0,el[0],n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,
              faceNumber,pface,pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,
                                                      pelemil,pelembl,pelement);
    for (i=1; el[i]; i++){
       el[i]->sons[0] = el[0]->sons[0];
       el[i]->eflag = 0;
    }
    el[0]->eflag = 0;
    recopy_faces(theGrid,no2,no3,no4,n22,n33,n44,n23,n24,n34,f1,fac[0]);
    recopy_faces(theGrid,no1,no3,no4,n11,n33,n44,n13,n14,n34,f2,fac[1]);
    recopy_faces(theGrid,no1,no2,no4,n11,n22,n44,n12,n14,n24,f3,fac[2]);
    recopy_faces(theGrid,no1,no2,no3,n11,n22,n33,n12,n13,n23,f4,fac[3]);
    for (i=0; el[i]; i++){
       if (el[i]->f[0]->sons[0] == NULL) el[i]->f[0]->sons[0] = el[i]->f[0];
       if (el[i]->f[1]->sons[0] == NULL) el[i]->f[1]->sons[0] = el[i]->f[1];
       if (el[i]->f[2]->sons[0] == NULL) el[i]->f[2]->sons[0] = el[i]->f[2];
       if (el[i]->f[3]->sons[0] == NULL) el[i]->f[3]->sons[0] = el[i]->f[3];
    }
 }
 
 void set_ind2(mg,n11,n22,n33,n12,n13,n23,f)
 MULTIGRID *mg;
 NODE *n11, *n22, *n33, *n12, *n13, *n23;
 FACE *f;
 {
    UNS type;
    
    type = FTYPE(f);
    if(f->sons[1] && f->index > 0)
       SET_ORIENT(f,correct_ind(TOP_GRID(mg),f,n11,n22,n33));
    set_norm_vect_or(n11,n22,n33,n12,n13,n23,f); 
    FTYPE(f) = type;
 }
 
/* refinement of products of Finer6 */
void second_refinement(mg,theElem,peil,pebl,pface,pfaceil,pfacebl,plink,pnflink,
                pfnlink,pflink,pelemil,pelembl,pelement,faceNumber,fictElNumber)
MULTIGRID *mg;
ELEMENT *theElem, *peil, *pebl, **pelemil, **pelembl, **pelement;
FACE **pface, **pfaceil, **pfacebl;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
INT  *faceNumber, *fictElNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34;
    
    NODES_OF_ELEMENT(n11,n22,n33,n44,theElem);
    n12 = middle(n11,n22); 
    n13 = middle(n11,n33);  
    n14 = middle(n11,n44);  
    n23 = middle(n22,n33);  
    n24 = middle(n22,n44);  
    n34 = middle(n33,n44);  
    nu_refine(1,theElem,n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,
       theElem->f[0],theElem->f[1],theElem->f[2],theElem->f[3],faceNumber,pface,
         pfaceil,pfacebl,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement);
    if (theElem->sons[0]) 
       treat_fictive_element(theElem,peil,pebl,fictElNumber);
    set_ind2(mg,n11,n22,n33,n12,n13,n23,theElem->f[3]); 
    set_ind2(mg,n11,n22,n44,n12,n14,n24,theElem->f[2]); 
    set_ind2(mg,n11,n33,n44,n13,n14,n34,theElem->f[1]); 
    set_ind2(mg,n22,n33,n44,n23,n24,n34,theElem->f[0]); 
 }

 /* double refinement of elements with status=1 */
 void refine3(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,
         pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,
                                             faceNumber,fictNumber,fictElNumber)
 MULTIGRID *mg;
 GRID *theGrid;
 ELEMENT *theElem, **pelemil, **pelembl, **pelement;
 NODE **pnode, **pnodeil, **pnodebl;
 FACE **pface, **pfaceil, **pfacebl;
 VERTEX **pvert;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 INT  *nodeNumber, *faceNumber, *fictNumber, *fictElNumber; 
 {
    NODE *n11, *n22, *n33, *n44, *n12, *n13, *n14, *n23, *n24, *n34,
         *no1, *no2, *no3, *no4;
    FACE *pfil, *pfbl, ff1, ff2, ff3, ff4, *f1, *f2, *f3, *f4, *fac[4][4];
    ELEMENT *peil, *pebl, *el[8];
    INT ind[4];
    INT i,j;
    
    pfil = *pfaceil;
    pfbl = *pfacebl;
    peil = *pelemil;
    pebl = *pelembl;
    data2_from_forefather(theGrid,theElem,&no1,&no2,&no3,&no4,el,fac,ind);
    treat_6edges(mg,no1,no2,no3,no4,&n11,&n22,&n33,&n44,&n12,&n13,&n14,
                   &n23,&n24,&n34,pnode,pnodeil,pnodebl,pvert,plink,nodeNumber);
    copy_faces(mg,theGrid,3,fac[0],ind[0],&ff1,&f1,no2,no3,no4,
               n22,n33,n44,n23,n24,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[1],ind[1],&ff2,&f2,no1,no3,no4,
               n11,n33,n44,n13,n14,n34,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[2],ind[2],&ff3,&f3,no1,no2,no4,
               n11,n22,n44,n12,n14,n24,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    copy_faces(mg,theGrid,3,fac[3],ind[3],&ff4,&f4,no1,no2,no3,
               n11,n22,n33,n12,n13,n23,pnode,pnodeil,pnodebl,pvert,plink,pface,
                                 pfaceil,pfacebl,pnflink,nodeNumber,faceNumber);
    Finer6(mg,el[0],n11,n22,n33,n44,n12,n13,n14,n23,n24,n34,f1,f2,f3,f4,pnode,
           pnodeil,pnodebl,pvert,plink,nodeNumber,faceNumber,pface,pfaceil,
                                              pfacebl,pelemil,pelembl,pelement);
    for (i=1; el[i]; i++){
       el[i]->sons[0] = el[0]->sons[0];
       el[i]->eflag = 0;
    }
    el[0]->eflag = 0;
    recopy_faces(theGrid,no2,no3,no4,n22,n33,n44,n23,n24,n34,f1,fac[0]);
    recopy_faces(theGrid,no1,no3,no4,n11,n33,n44,n13,n14,n34,f2,fac[1]);
    recopy_faces(theGrid,no1,no2,no4,n11,n22,n44,n12,n14,n24,f3,fac[2]);
    recopy_faces(theGrid,no1,no2,no3,n11,n22,n33,n12,n13,n23,f4,fac[3]);
    for (i=0; i<8; i++)
      second_refinement(mg,el[0]->sons[i],peil,pebl,pface,pfaceil,pfacebl,plink,
       pnflink,pfnlink,pflink,pelemil,pelembl,pelement,faceNumber,fictElNumber);
    for (i=0; i<8; i++)
       for (j=0; j<4; j++)
          if (el[0]->sons[i]->f[j]->sons[0] == el[0]->sons[i]->f[j])
             el[0]->sons[i]->f[j]->sons[0] = NULL;
          else if (el[0]->sons[i]->f[j]->sons[0] && 
                   el[0]->sons[i]->f[j]->index > TOP_GRID(mg)->maxFaceIndex)
             treat_fictive_face(el[0]->sons[i]->f[j],pfil,pfbl,fictNumber);
    for (i=0; el[i]; i++){
       if (el[i]->f[0]->sons[0] == NULL) el[i]->f[0]->sons[0] = el[i]->f[0];
       if (el[i]->f[1]->sons[0] == NULL) el[i]->f[1]->sons[0] = el[i]->f[1];
       if (el[i]->f[2]->sons[0] == NULL) el[i]->f[2]->sons[0] = el[i]->f[2];
       if (el[i]->f[3]->sons[0] == NULL) el[i]->f[3]->sons[0] = el[i]->f[3];
    }
 }
 
 void additional_neighbours(mg,n1,n2,n3,n4,f1,f2,f3,f4,plink,pnflink,pfnlink,
                                                                         pflink)
 MULTIGRID *mg;
 NODE *n1, *n2, *n3, *n4;
 FACE *f1, *f2, *f3, *f4;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 {
    FNLINK *pfnl;
    FLINK *pfl;
    
    if (NODE_ON_NEW_LEVEL(n1,mg) && IS_FN(n1)){
       if (IS_FN(n2) && !Get_link(n1,n2))
          new_neighbour(n1,n2,2,0,plink);
       if (IS_FN(n3) && !Get_link(n1,n3))
          new_neighbour(n1,n3,2,0,plink);
       if (IS_FN(n4) && !Get_link(n1,n4))
          new_neighbour(n1,n4,2,0,plink);
       if (IS_FF(f1))
          new_Fneighbour(n1,f1,pnflink);
       if (IS_FF(f2) && !FACE_ON_NEW_LEVEL(f2,mg) && !get_nflink(n1,f2))
          new_Fneighbour(n1,f2,pnflink);
       if (IS_FF(f3) && !FACE_ON_NEW_LEVEL(f3,mg) && !get_nflink(n1,f3))
          new_Fneighbour(n1,f3,pnflink);
       if (IS_FF(f4) && !FACE_ON_NEW_LEVEL(f4,mg) && !get_nflink(n1,f4))
          new_Fneighbour(n1,f4,pnflink);
    }
    if (IS_FF(f1) && FACE_ON_NEW_LEVEL(f1,mg)){
       if (IS_FN(n1)){
          if (f1->fnstart == NULL)
             f1->fnstart = *pfnlink;
          else{
             for (pfnl = f1->fnstart; pfnl->next != NULL; pfnl = pfnl->next);
             pfnl->next = *pfnlink;
          }
          (*pfnlink)->nbnode = n1;
          ((*pfnlink)++)->next = NULL;
       }
       if (IS_FF(f2) || IS_FF(f3) || IS_FF(f4)){
          if (f1->fstart == NULL)
             f1->fstart = *pflink;
          else{
             for (pfl = f1->fstart; pfl->next != NULL; pfl = pfl->next);
             pfl->next = *pflink;
          }
          if (IS_FF(f2)){
             (*pflink)->next = *pflink + 1;
             ((*pflink)++)->nbface = f2;
          } 
          if (IS_FF(f3)){
             (*pflink)->next = *pflink + 1;
             ((*pflink)++)->nbface = f3;
          } 
          if (IS_FF(f4)){
             (*pflink)->next = *pflink + 1;
             ((*pflink)++)->nbface = f4;
          } 
          (*pflink - 1)->next = NULL;
       }
    }
 }
          
 void add_neighbours_from_old_levels(mg,plink,pnflink,pfnlink,pflink)
 MULTIGRID *mg;
 LINK **plink;
 FLINK **pflink;
 NFLINK **pnflink;
 FNLINK **pfnlink;
 {         
    NODE *n1, *n2, *n3, *n4;
    FACE *f1, *f2, *f3, *f4;
    ELEMENT *pel;
    GRID *theGrid;
    
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
          if (IS_TOP_ELEMENT(pel)){
             TOPNODES_OF_ELEMENT(n1,n2,n3,n4,pel); 
             TOPFACES_OF_ELEMENT(f1,f2,f3,f4,pel);  
             additional_neighbours(mg,n1,n2,n3,n4,f1,f2,f3,f4,plink,pnflink,
                                                                pfnlink,pflink);
             additional_neighbours(mg,n2,n3,n4,n1,f2,f3,f4,f1,plink,pnflink,
                                                                pfnlink,pflink);
             additional_neighbours(mg,n3,n4,n1,n2,f3,f4,f1,f2,plink,pnflink,
                                                                pfnlink,pflink);
             additional_neighbours(mg,n4,n1,n2,n3,f4,f1,f2,f3,plink,pnflink,
                                                                pfnlink,pflink);
          }
 }

#else

void mark_primary_unif_refinement(mg)   
MULTIGRID *mg;
{  eprintf("Error: mark_primary_unif_refinement not available.\n");  }

INT mark_secondary_unif_refinement(mg,i1,i2) 
MULTIGRID *mg; INT i1,i2;
{  eprintf("Error: mark_secondary_unif_refinement not available.\n");  }

INT test_green_faces(mg,pnode,pnodeil,pnodebl,nodeNumber,i1)
MULTIGRID *mg; NODE **pnode, **pnodeil, **pnodebl; INT  *nodeNumber, i1;
{  eprintf("Error: test_green_faces not available.\n");  }

void mark_new_green_elements_and_double_refinement(mg,i1) 
MULTIGRID *mg; INT i1;
{  eprintf("Error: mark_new_green_elements_and_double_refinement not available.\n");  }

void refine1(mg,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
MULTIGRID *mg; ELEMENT *theElem, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber; 
{  eprintf("Error: refine1 not available.\n");  }

void refine2(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber)
MULTIGRID *mg; GRID *theGrid; ELEMENT *theElem, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber;
{  eprintf("Error: refine2 not available.\n");  }

void refine3(mg,theGrid,theElem,pnode,pnodeil,pnodebl,pface,pfaceil,pfacebl,pvert,plink,pnflink,pfnlink,pflink,pelemil,pelembl,pelement,nodeNumber,faceNumber,fictNumber,fictElNumber)
MULTIGRID *mg; GRID *theGrid; ELEMENT *theElem, **pelemil, **pelembl, **pelement; NODE **pnode, **pnodeil, **pnodebl; FACE **pface, **pfaceil, **pfacebl; VERTEX **pvert; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink; INT  *nodeNumber, *faceNumber, *fictNumber, *fictElNumber; 
{  eprintf("Error: refine3 not available.\n");  }

void add_neighbours_from_old_levels(mg,plink,pnflink,pfnlink,pflink)
MULTIGRID *mg; LINK **plink; FLINK **pflink; NFLINK **pnflink; FNLINK **pfnlink;
{  eprintf("Error: add_neighbours_from_old_levels not available.\n");  }

#endif

void perform_marking_of_edges(mg,pnode,pnodeil,pnodebl,nodeNumber)
MULTIGRID *mg;
NODE **pnode, **pnodeil, **pnodebl;
INT  *nodeNumber;
{
   INT i=1, i1=5, i2=6;
   
   prepare_links(mg);
   mark_primary_unif_refinement(mg);
   while( mark_secondary_unif_refinement(mg,(i1+=2),(i2+=2)) ) i++;
   printf("%i cycles for marking the secondary uniform refinement.\n",i);
   while( test_green_faces(mg,pnode,pnodeil,pnodebl,nodeNumber,i1) ){
      i=1;
      while( mark_secondary_unif_refinement(mg,(i1+=2),(i2+=2)) ) i++;
      printf("%i additional cycles ",i);
      printf("for marking the secondary uniform refinement.\n");
   }
   mark_new_green_elements_and_double_refinement(mg,i1);     
}   

void statistics(mg,elemNumber,fictNumber,fictElNumber)  
MULTIGRID *mg;
INT elemNumber,fictNumber,fictElNumber;
{
   NODE *pn;
   FACE *pf;
   ELEMENT *pe;
   GRID *pgrid;
   INT nodeNumber, faceNumber, elNumber, tfn, tnn;
   
   tnn = TOP_GRID(mg)->maxNodeIndex - TOP_GRID(mg)->minNodeIndex + 1;
   tfn = TOP_GRID(mg)->maxFaceIndex - TOP_GRID(mg)->minFaceIndex + 1-fictNumber;
   nodeNumber = tnn;
   faceNumber = tfn;
   elNumber   = elemNumber -= fictElNumber;
   for (pgrid = FIRSTGRID(mg); pgrid->finer != NULL; pgrid = pgrid->finer){
      for(pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
         if (IS_TOP_NODE(pn)) ++nodeNumber;
      for(pf = FIRSTF(pgrid); pf != NULL; pf = pf->succ)
         if (IS_TOP_FACE(pf)) ++faceNumber; 
      for(pe = FIRSTELEMENT(pgrid); pe != NULL; pe = pe->succ)
         if (IS_TOP_ELEMENT(pe)) ++elNumber; 
   }
   printf("Level %2i      | top level | top triang.| multigrid\n",mg->toplevel);
   printf("--------------------------------------------------\n");
   printf("node number   |%8i   |%9i   | %7i\n",tnn,nodeNumber,mg->nodeNumber);
   printf("face number   |%8i   |%9i   | %7i\n",tfn,faceNumber,
                                               mg->faceNumber - mg->fictNumber);
   printf("element number|%8i   |%9i   | %7i\n",elemNumber,elNumber,
                                             mg->elemNumber - mg->fictElNumber);
   printf("n. fict. fac. |%8i   |            | %7i\n",fictNumber,mg->fictNumber);
   printf("n. fict. el.  |%8i   |            | %7i\n",fictElNumber,
                                                              mg->fictElNumber);
   printf("tot. n. fac.  |%8i   |            | %7i\n",tfn + fictNumber,
                                                                mg->faceNumber);
   printf("tot. n. el.   |%8i   |            | %7i\n\n",
                                      elemNumber + fictElNumber,mg->elemNumber);
}

void nonuniform_refinement(mg,pnode,pface,pvert,plink,pnflink,pfnlink,pflink,
                                                                       pelement)
MULTIGRID *mg;
ELEMENT **pelement;
NODE **pnode;
FACE **pface;
VERTEX **pvert;
LINK **plink;
FLINK **pflink;
NFLINK **pnflink;
FNLINK **pfnlink;
{
   NODE lin, ldbn;
   FACE lif, ldbf, *pf;
   ELEMENT lie, ldbe, *theElem, *pelemil, *pelembl;
   GRID *pgrid, *theGrid;
   INT nodeNumber, faceNumber, elemNumber,fictNumber=0, fictElNumber=0;
   
   lin.succ = NULL;
   lif.succ = NULL;
   lie.succ = NULL;
   ldbn.succ = NULL;
   ldbf.succ = NULL;
   ldbe.succ = NULL;    
   pgrid = TOP_GRID(mg) + 1;
   pgrid->level = mg->toplevel + 1;
   pgrid->lastNode = &lin;
   pgrid->ldbn = &ldbn;
   pgrid->lastFace = &lif;
   pgrid->ldbf = &ldbf;
   pelemil = &lie;
   pelembl = &ldbe;
   pgrid->coarser = TOP_GRID(mg);
   pgrid->finer = NULL;
   pgrid->first_grid = TOP_GRID(mg)->first_grid;
   for (theGrid = pgrid; theGrid; theGrid = theGrid->coarser)
      theGrid->top_grid = pgrid;
   nodeNumber = mg->nodeNumber;
   faceNumber = mg->faceNumber;
   perform_marking_of_edges(mg,pnode,&pgrid->lastNode,&pgrid->ldbn,&nodeNumber);
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (theElem = FIRSTELEMENT(theGrid); theElem != NULL; theElem = 
                                                                  theElem->succ)
         if (IS_TOP_ELEMENT(theElem) && theElem->eflag >0)
            if (!theElem->status)
               refine1(mg,theElem,pnode,&pgrid->lastNode,&pgrid->ldbn,pface,
                     &pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,
                     pflink,&pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
            else if (theElem->eflag < 4)
               refine2(mg,theGrid,theElem,pnode,&pgrid->lastNode,&pgrid->ldbn,
                pface,&pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,pfnlink,
                pflink,&pelemil,&pelembl,pelement,&nodeNumber,&faceNumber);
            else
               refine3(mg,theGrid,theElem,pnode,&pgrid->lastNode,&pgrid->ldbn,
                       pface,&pgrid->lastFace,&pgrid->ldbf,pvert,plink,pnflink,
                       pfnlink,pflink,&pelemil,&pelembl,pelement,&nodeNumber,
                       &faceNumber,&fictNumber,&fictElNumber);
   add_neighbours_from_old_levels(mg,plink,pnflink,pfnlink,pflink);
   pgrid->firstNode = lin.succ;
   pgrid->fdbn = ldbn.succ;
   pgrid->firstFace = lif.succ;
   pgrid->fdbf = ldbf.succ;
   FIRSTELEMENT(pgrid) = lie.succ;
   FDBE(pgrid) = ldbe.succ;
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
   pelemil->succ = FDBE(pgrid);
   pgrid->minNodeIndex = (pgrid-1)->maxNodeIndex + 1;
   pgrid->maxNodeIndex = nodeNumber;
   pgrid->minFaceIndex = (pgrid-1)->maxFaceIndex + 1;    
   pgrid->maxFaceIndex = faceNumber;
   elemNumber = ((long)(*pelement) - 
      (long)(TMIN(FIRSTELEMENT(pgrid),FDBE(pgrid))))/sizeof(ELEMENT);
   mg->nodeNumber = nodeNumber;
   mg->faceNumber = faceNumber;
   mg->elemNumber += elemNumber;
   mg->fictNumber += fictNumber;
   mg->fictElNumber += fictElNumber;
   TOP_GRID(mg)->finer = pgrid;
   mg->toplevel++;
   TOP_GRID(mg) = pgrid;
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pf = FIRSTF(theGrid); pf != NULL; pf = pf->succ)
         if (IS_TOP_FACE(pf))
            SET_FTOPLEVEL(pf,mg->toplevel);
   check_new_level(mg);
   statistics(mg,elemNumber,fictNumber,fictElNumber);
}

