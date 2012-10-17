/******************************************************************************/
/*                                                                            */
/*                               triangulation                                */
/*                                                                            */
/******************************************************************************/

FLOAT barycentric_coordinates();

FACE *top_face(f)  /* can be used only for faces corresponding */
FACE *f;           /* to the highest triangulation             */
{
   while(f->sons[0]) f=f->sons[0];
   return(f);
}

void eprintf(text)
char text[];
{
   FILE *fp;

   NUMBER_OF_ERRORS++;
   printf(text);
   fp = fopen("detected_errors","a");
   fprintf(fp,text);
   fclose(fp);
}

void eeprintf(text,name,i,j)
char text[], name[], i, j;
{
   FILE *fp;

   NUMBER_OF_ERRORS++;
   printf(text,name,i,j);
   fp = fopen("detected_errors","a");
   fprintf(fp,text,name,i,j);
   fclose(fp);
}

void print_number_of_errors()
{
   FILE *fp;

   fp = fopen("detected_errors","a");
   if (NUMBER_OF_ERRORS == 1){
      printf("\n1 error detected.\n");
      fprintf(fp,"\n1 error detected.\n");
   }
   else{
      printf("\n%i errors detected.\n",NUMBER_OF_ERRORS);
      fprintf(fp,"\n%i errors detected.\n",NUMBER_OF_ERRORS);
   }
   fclose(fp);
}

int round0(x)
double x;
{
   int i=floor(x), j=ceil(x);  //  i <= x <= j

   if (x - i < j - x)
      return(i);
   else
      return(j);
}

void order2(i1,i2)
INT *i1, *i2;
{
   INT k;
  
   if (*i1 > *i2){
      k = *i1;
      *i1 = *i2;
      *i2 = k;
   }
}

void order3(i1,i2,i3)
INT *i1, *i2, *i3;
{
   order2(i1,i2);
   order2(i1,i3);
   order2(i2,i3);
}

void order4(i1,i2,i3,i4)
INT *i1, *i2, *i3, *i4;
{
   order2(i1,i2);
   order2(i1,i3);
   order2(i1,i4);
   order2(i2,i3);
   order2(i2,i4);
   order2(i3,i4);
}  

#if (DIM == 3 && ELEMENT_TYPE == SIMPLEX)

void elem(i1,i2,i3,i4,nodes,pelement)
INT i1,i2,i3,i4;
NODE *nodes;
ELEMENT **pelement;
{
   order4(&i1,&i2,&i3,&i4); 
   (*pelement)->succ = *pelement + 1;
   (*pelement)->status = 0;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = &nodes[i1];
   (*pelement)->n[1] = &nodes[i2];
   (*pelement)->n[2] = &nodes[i3];
   (*pelement)->n[3] = &nodes[i4];
   (*pelement)->father = NULL;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}

#elif (DIM == 2 && ELEMENT_TYPE == CUBE)

void elem(i1,i2,i3,i4,nodes,pelement)
INT i1,i2,i3,i4;
NODE *nodes;
ELEMENT **pelement;
{
   (*pelement)->succ = *pelement + 1;
   (*pelement)->status = 0;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = &nodes[i1];
   (*pelement)->n[1] = &nodes[i2];
   (*pelement)->n[2] = &nodes[i3];
   (*pelement)->n[3] = &nodes[i4];
   (*pelement)->father = NULL;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}

#elif (DIM == 2 && ELEMENT_TYPE == SIMPLEX)

void elem(i1,i2,i3,nodes,pelement)
INT i1,i2,i3;
NODE *nodes;
ELEMENT **pelement;
{
   order3(&i1,&i2,&i3);
   (*pelement)->succ = *pelement + 1;
   (*pelement)->status = 0;
   (*pelement)->eflag = 0;
   (*pelement)->n[0] = &nodes[i1];
   (*pelement)->n[1] = &nodes[i2];
   (*pelement)->n[2] = &nodes[i3];
   (*pelement)->father = NULL;
   (*pelement)->sons[0] = NULL;
   (*pelement)++;
}

#endif

#if DATA_STR & LG_DATA
#define PNODE_LGSTART  pnode->lgstart = NULL; 
#define PNODE_LGD      pnode->lgd = NULL;
#define PFACE_LGSTART  pface->lgstart = NULL;
#else
#define PNODE_LGSTART  
#define PNODE_LGD      
#define PFACE_LGSTART  
#endif

#if DATA_S & N_LINK_TO_FACES
#define PNODE_NFSTART  pnode->nfstart = NULL;   pnode->tnfstart = NULL;
#else
#define PNODE_NFSTART  
#endif

#if DATA_S & N_LINK_TO_ELEMENTS
#define PNODE_NESTART  pnode->nestart = NULL;
#else
#define PNODE_NESTART  
#endif

#if DATA_S & F_LINK_TO_NODES
#define PFACE_FNSTART  pface->fnstart = NULL;   pface->tfnstart = NULL;
#else
#define PFACE_FNSTART  
#endif

#if DATA_S & F_LINK_TO_FACES
#define PFACE_FSTART   pface->fstart = NULL;    pface->tfstart = NULL;
#else
#define PFACE_FSTART  
#endif

#if DATA_S & SPECIAL_NODES_AND_FACES
#define SET_SNODE      pnode->s_node = NULL;
#define SET_SFACE      pface->s_face = NULL;
#else
#define SET_SNODE      
#define SET_SFACE      
#endif

#if DIM == 3
#define PFACE_SONS  pface->sons[0] = NULL;  pface->sons[1] = NULL;             \
                    pface->sons[2] = NULL;  pface->sons[3] = NULL;
#else
#define PFACE_SONS  pface->sons[0] = NULL;  pface->sons[1] = NULL;
#endif

#if F_DATA & CURVED_FACE_MIDDLE
#define SET_F_CURVED_MIDDLE     pface->c_midpoint = NULL;
#else
#define SET_F_CURVED_MIDDLE
#endif

#if MOVING_BOUNDARY == YES
#define SET_NEW_VERTEX     pnode->newvertex = vertex;
#else
#define SET_NEW_VERTEX
#endif

#if N_DATA & NODE_IFLAG
#define SET_NODE_IFLAG     IFLAG(pnode) = 0; 
#else
#define SET_NODE_IFLAG
#endif

#if F_DATA & FACE_IFLAG
#define SET_FACE_IFLAG     IFLAG(pface) = 0; 
#else
#define SET_FACE_IFLAG
#endif

#if N_DATA & NUMBER_OF_N_NEL
#define SET_N_NEL          pnode->nel = 0;
#else
#define SET_N_NEL
#endif

void nmake(pnode,psucc,index,father,vertex)
NODE *pnode, *psucc, *father;
INT index;
VERTEX *vertex;
{
   pnode->succ = psucc;
   pnode->index = index;
   pnode->start = NULL;
   pnode->tstart = NULL;
   pnode->father = father;
   pnode->son = NULL;
   pnode->myvertex = vertex;
   SET_NEW_VERTEX
   vertex->topnode = pnode;
   SET_NODE_IFLAG
   SET_N_NEL
   SET_SNODE      
   PNODE_NFSTART
   PNODE_NESTART
   PNODE_LGSTART
   PNODE_LGD
}

void fmake(pface,psucc,type,or_set,index,father)
FACE *pface, *psucc, *father;
INT type, index;
UNS or_set;
{
   pface->succ = psucc;
   pface->type = ( ((UNS)type) & ~OR_SET ) | or_set;
   pface->index = index;
   pface->father = father;
   SET_FACE_IFLAG
   SET_SFACE      
   PFACE_SONS
   SET_F_CURVED_MIDDLE
   PFACE_FSTART
   PFACE_FNSTART
   PFACE_LGSTART
}

#if (DIM == 3 && ELEMENT_TYPE == SIMPLEX)

void cube(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, k, n, nf, il, nv, na; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv*nv;         /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) 
         for (k = 0; k < nv; k++) { 
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->x[2] = k/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[2]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[2]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   n = 0;
   nf = nv*nv;
   for (i = 0; i < nv-1; i++)
      for (j = -1; j < nv-2; j++)
         for (k = 0; k < nv-1; k++) {
            n++;
            il = n + j + (2*nv-1)*i;
            elem(il, il+nf+nv  , il+nv     , il+nv+1, nodes, pelement); 
            elem(il, il+nf+nv+1, il+nf+nv  , il+nv+1, nodes, pelement);    
            elem(il, il+1      , il+nf+nv+1, il+nv+1, nodes, pelement); 
            elem(il+nf+nv+1, il+1, il+nf+1 , il+nf, nodes, pelement);    
            elem(il+nf+nv+1, il+1, il      , il+nf, nodes, pelement); 
            elem(il+nf+nv+1, il  , il+nf+nv, il+nf, nodes, pelement);
   }
   (*pelement - 1)->succ = NULL;      
}

void read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, i4, type, nvert; 
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),
                                 &((*pvert)->x[2]),&type);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   while (i1 > -1){
      elem(i1,i2,i3,i4,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   }
   (*pelement - 1)->succ = NULL;
}

void Read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, i4, t1, t2, t3, nvert; 
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %lf %i %i %i",&((*pvert)->x[0]),&((*pvert)->x[1]),
                                       &((*pvert)->x[2]),&t1,&t2,&t3);
      (*pvert)->type = (UNS)t1;
      if (t2 == 1)
         (*pvert)->type = (*pvert)->type | GW;
      else if (t2 == 2)
         (*pvert)->type = (*pvert)->type | GLG;
      else if (t2 == 3)
         (*pvert)->type = (*pvert)->type | GLS;
      else if (t2 != 0)
         eprintf("Error in Read_vertexes_and_elements.\n");
      if (t3 == 1){
         (*pvert)->type = (*pvert)->type | GW;
         (*pvert)->type = (*pvert)->type | GLG;
      }
      else if (t3 == 2){
         (*pvert)->type = (*pvert)->type | GLG;
         (*pvert)->type = (*pvert)->type | GLS;
      }
      else if (t3 == 3)
         (*pvert)->type = (*pvert)->type | GWW;
      else if (t3 != 0)
         eprintf("Error in Read_vertexes_and_elements.\n");
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   while (i1 > -1){
      elem(i1,i2,i3,i4,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   }
   (*pelement - 1)->succ = NULL; 
}

#elif (DIM == 2 && ELEMENT_TYPE == SIMPLEX)

void cube(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na, type=CUBE_TYPE; /*  1, 2 ... Friedrichs-Keller  */
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1;
   }
   for (i = 0; i < na; i++)
      if (vertexes[i].type){
         if (SC_EXAMPLE != 9 && SC_EXAMPLE != 91)
            if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (SC_EXAMPLE != 80)
            if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (SC_EXAMPLE != 80)
            if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= NMASK;
/*
         if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA){
            if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 
                                                           NMASK|NMASK_FOR_SF;
            else if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 
                                                           NMASK|NMASK_FOR_SF;
            else vertexes[i].type |= NMASK_FOR_SF|SNMASK;
         }
*/
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            if (type == 1){
               elem(il, il+nv, il+nv+1, nodes, pelement); 
               elem(il, il+1,  il+nv+1, nodes, pelement);    
            }
            else if (type == 2){
               elem(il,   il+1,  il+nv,   nodes, pelement);    
               elem(il+1, il+nv, il+nv+1, nodes, pelement); 
            }
            else if (type == 3){
               if (2*((i+j)/2) == i+j){
                  elem(il, il+nv, il+nv+1, nodes, pelement); 
                  elem(il, il+1,  il+nv+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,  il+nv,   nodes, pelement);    
                  elem(il+1, il+nv, il+nv+1, nodes, pelement); 
               }
            }
            else if (type == 4){
               if (2*((i+j)/2) == i+j){
                  elem(il,   il+1,  il+nv,   nodes, pelement);    
                  elem(il+1, il+nv, il+nv+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nv, il+nv+1, nodes, pelement); 
                  elem(il, il+1,  il+nv+1, nodes, pelement);    
               }
            }
            else if (type == 5){
               if (2*(j/2) == j){
                  elem(il, il+nv, il+nv+1, nodes, pelement); 
                  elem(il, il+1,  il+nv+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,  il+nv,   nodes, pelement);    
                  elem(il+1, il+nv, il+nv+1, nodes, pelement); 
               }
            }
            else if (type == 6){
               if (2*(j/2) == j){
                  elem(il,   il+1,  il+nv,   nodes, pelement);    
                  elem(il+1, il+nv, il+nv+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nv, il+nv+1, nodes, pelement); 
                  elem(il, il+1,  il+nv+1, nodes, pelement);    
               }
            }
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_xy(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n=0, il, nvx=NVX, nvy=NVY, na=NVX*NVY, type=CUBE_TYPE;
   NODE *nodes;                        /*  type = 1, 2 ... Friedrichs-Keller  */
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           (*pvert)->x[0] = i/(nvx-1.);
           (*pvert)->x[1] = j/(nvy-1.);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1;
   }
   for (i = 0; i < na; i++)
      if (vertexes[i].type){
         if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++) {
            il = j + nvy*i;
            if (type == 1){
               elem(il, il+nvy, il+nvy+1, nodes, pelement); 
               elem(il, il+1,   il+nvy+1, nodes, pelement);    
            }
            else if (type == 2){
               elem(il,   il+1,   il+nvy,   nodes, pelement);    
               elem(il+1, il+nvy, il+nvy+1, nodes, pelement); 
            }
            else if (type == 3){
               if (2*((i+j)/2) == i+j){
                  elem(il, il+nvy, il+nvy+1, nodes, pelement); 
                  elem(il, il+1,   il+nvy+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nvy,   nodes, pelement);    
                  elem(il+1, il+nvy, il+nvy+1, nodes, pelement); 
               }
            }
            else if (type == 4){
               if (2*((i+j)/2) == i+j){
                  elem(il,   il+1,   il+nvy,   nodes, pelement);    
                  elem(il+1, il+nvy, il+nvy+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nvy, il+nvy+1, nodes, pelement); 
                  elem(il, il+1,   il+nvy+1, nodes, pelement);    
               }
            }
            else if (type == 5){
               if (2*(j/2) == j){
                  elem(il, il+nvy, il+nvy+1, nodes, pelement); 
                  elem(il, il+1,   il+nvy+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nvy,   nodes, pelement);    
                  elem(il+1, il+nvy, il+nvy+1, nodes, pelement); 
               }
            }
            else if (type == 6){
               if (2*(j/2) == j){
                  elem(il,   il+1,   il+nvy,   nodes, pelement);    
                  elem(il+1, il+nvy, il+nvy+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nvy, il+nvy+1, nodes, pelement); 
                  elem(il, il+1,   il+nvy+1, nodes, pelement);    
               }
            }
   }
   (*pelement - 1)->succ = NULL;      
}

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void cube_xy_surrounding(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, k=0, n=0, il, nvx=NVX, nvy=NVY, na=0, type=CUBE_TYPE, ind[100000],
       diff=0, i1, i2, i3, it, nvert;
   DOUBLE ax, bx, ay, by, x0, x1;
   NODE *nodes;
   VERTEX *vertexes;
   ELEMENT *pel;
   FILE *fp;
 
   ax = (NVX-1.)/(NV-1.);
   ay = (NVY-1.)/(NV-1.);
   bx = 0.5*(1.-ax);
   by = 0.5*(1.-ay);
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           x0 = ax*i/(nvx-1.) + bx;
           x1 = ay*j/(nvy-1.) + by;
           if (x0 > 0.999999 || x0 < 0.000001 ||
               x1 > 0.999999 || x1 < 0.000001){
              (*pvert)->x[0] = x0;
              (*pvert)->x[1] = x1;
              (*pvert)->type = 0; 
              nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
              if (x0 > -0.299999 && x0 < 1.299999 && 
                  x1 > -0.499999 && x1 < 1.499999)
                 ITYPE(*pnode) = 16;
              else
                 ITYPE(*pnode) = 32;
              (*pnode)++;
              na++;
              ind[k] = k-diff;
           }
           else{
              ind[k] = -3;
              diff++;
           }
           k++;
      } 
    
   for (i = 0; i < na; i++){
         if (fabs(vertexes[i].x[0]-bx) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[0]-ax-bx) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[1]-by) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[1]-ay-by) < EPSA) vertexes[i].type = 1|NMASK;
   }
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++){ 
         il = j + nvy*i;
         if (ind[il] > -1 && ind[il+nvy] > -1 && ind[il+nvy+1] > -1 &&
             ind[il+1] > -1){
            if (ITYPE(nodes + ind[il]) == 32 && 
                ITYPE(nodes + ind[il+nvy]) == 32 &&
                ITYPE(nodes + ind[il+nvy+1]) == 32)
               ITYPE(*pelement) = 32;
            else
               ITYPE(*pelement) = 16;
            elem(ind[il], ind[il+nvy], ind[il+nvy+1], nodes, pelement); 
            if (ITYPE(nodes + ind[il]) == 32 && 
                ITYPE(nodes + ind[il+1]) == 32 &&
                ITYPE(nodes + ind[il+nvy+1]) == 32)
               ITYPE(*pelement) = 32;
            else
               ITYPE(*pelement) = 16;
            elem(ind[il], ind[il+1],   ind[il+nvy+1], nodes, pelement);    
         }
      }
   (*pelement - 1)->succ = NULL;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
      if (ITYPE(pel) == 16)
         for (i = 0; i < NVERT; i++)
            if (ITYPE(pel->n[i]) == 32)
               ITYPE(pel->n[i]) = 48;
   (*pelement - 1)->succ = *pelement;
   
   fp = fopen("data_circular_magnet","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type,&it);
      (*pvert)->type = 0;
      nmake(*pnode,*pnode+1,n+i,(NODE*)NULL,(*pvert)++);
      ITYPE(*pnode) = it;
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n+nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   while (i1 > -1){
      ITYPE(*pelement) = it;
      elem(n+i1,n+i2,n+i3,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   }
   (*pelement - 1)->succ = NULL;
}

void cube_xy2_surrounding(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, k=0, n=0, il, nv=21, nvx=61, nvy=61, na=0, type=CUBE_TYPE, 
       diff=0, i1, i2, i3, it, nvert, ind[100000];
   DOUBLE ax, bx, ay, by, x0, x1;
   NODE *nodes;
   VERTEX *vertexes;
   ELEMENT *pel;
   FILE *fp;
 
   ax = (nvx-1.)/(nv-1.);
   ay = (nvy-1.)/(nv-1.);
   bx = 0.5*(1.-ax);
   by = 0.5*(1.-ay);
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           x0 = ax*i/(nvx-1.) + bx;
           x1 = ay*j/(nvy-1.) + by;
           if (x0 > 0.999999 || x0 < 0.000001 ||
               x1 > 0.999999 || x1 < 0.000001){
              (*pvert)->x[0] = x0;
              (*pvert)->x[1] = x1;
              (*pvert)->type = 0; 
              nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
              if (x0 > -0.299999 && x0 < 1.299999 && 
                  x1 > -0.499999 && x1 < 1.499999)
                 ITYPE(*pnode) = 16;
              else
                 ITYPE(*pnode) = 32;
              (*pnode)++;
              na++;
              ind[k] = k-diff;
           }
           else{
              ind[k] = -3;
              diff++;
           }
           k++;
      } 
    
   for (i = 0; i < na; i++){
         if (fabs(vertexes[i].x[0]-bx) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[0]-ax-bx) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[1]-by) < EPSA) vertexes[i].type = 1|NMASK;
         if (fabs(vertexes[i].x[1]-ay-by) < EPSA) vertexes[i].type = 1|NMASK;
   }
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++){ 
         il = j + nvy*i;
         if (ind[il] > -1 && ind[il+nvy] > -1 && ind[il+nvy+1] > -1 &&
             ind[il+1] > -1){
            if ((MYVERTEX(nodes + ind[il+nvy+1])->x[0] < 0.5000001 &&
                 MYVERTEX(nodes + ind[il+nvy+1])->x[1] < 0.5000001) ||
                (MYVERTEX(nodes + ind[il])->x[0] > 0.4999999 &&
                 MYVERTEX(nodes + ind[il])->x[1] > 0.4999999)){
               if (ITYPE(nodes + ind[il]) == 32 && 
                   ITYPE(nodes + ind[il+1]) == 32 &&
                   ITYPE(nodes + ind[il+nvy]) == 32)
                  ITYPE(*pelement) = 32;
               else
                  ITYPE(*pelement) = 16;
               elem(ind[il], ind[il+1], ind[il+nvy], nodes, pelement); 
               if (ITYPE(nodes + ind[il+1]) == 32 &&
                   ITYPE(nodes + ind[il+nvy]) == 32 &&
                   ITYPE(nodes + ind[il+nvy+1]) == 32)
                  ITYPE(*pelement) = 32;
               else
                  ITYPE(*pelement) = 16;
               elem(ind[il+1], ind[il+nvy],   ind[il+nvy+1], nodes, pelement);    
            }
            else{
               if (ITYPE(nodes + ind[il]) == 32 && 
                   ITYPE(nodes + ind[il+nvy]) == 32 &&
                   ITYPE(nodes + ind[il+nvy+1]) == 32)
                  ITYPE(*pelement) = 32;
               else
                  ITYPE(*pelement) = 16;
               elem(ind[il], ind[il+nvy], ind[il+nvy+1], nodes, pelement); 
               if (ITYPE(nodes + ind[il]) == 32 && 
                   ITYPE(nodes + ind[il+1]) == 32 &&
                   ITYPE(nodes + ind[il+nvy+1]) == 32)
                  ITYPE(*pelement) = 32;
               else
                  ITYPE(*pelement) = 16;
               elem(ind[il], ind[il+1],   ind[il+nvy+1], nodes, pelement);    
            }
         }
      }
   (*pelement - 1)->succ = NULL;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
      if (ITYPE(pel) == 16)
         for (i = 0; i < NVERT; i++)
            if (ITYPE(pel->n[i]) == 32)
               ITYPE(pel->n[i]) = 48;
   (*pelement - 1)->succ = *pelement;
   
   fp = fopen("data_circular_magnet","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type,&it);
      (*pvert)->type = 0;
      nmake(*pnode,*pnode+1,n+i,(NODE*)NULL,(*pvert)++);
      ITYPE(*pnode) = it;
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n+nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   while (i1 > -1){
      ITYPE(*pelement) = it;
      elem(n+i1,n+i2,n+i3,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   }
   (*pelement - 1)->succ = NULL;
}

#else

void cube_xy_surrounding(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: cube_xy_surrounding not available.\n"); }

void cube_xy2_surrounding(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: cube_xy2_surrounding not available.\n"); }

#endif

void step_xy(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n=0, il, nvx=NVX, nvy=NVY, nvx2, nvy2, nny, type=CUBE_TYPE;
   NODE *nodes;                        /*  type = 1, 2 ... Friedrichs-Keller  */
   VERTEX *vertexes;
 
   if (2*(nvx/2) == nvx)
      nvx += 1;
   if (2*(nvy/2) == nvy)
      nvy += 1;
   nvx2 = nvx/2;
   nvy2 = nvy/2;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) 
         if (i >= nvx2 || j >= nvy2){
            (*pvert)->x[0] = i/(nvx-1.);
            (*pvert)->x[1] = j/(nvy-1.);
            (*pvert)->type = 0; 
            nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
            (*pnode)++;
         }
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < n; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-0.5) < EPSA &&
               vertexes[i].x[1] < 0.5+EPSA)  vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.5) < EPSA &&
               vertexes[i].x[0] < 0.5+EPSA)  vertexes[i].type = 1;
   }
   for (i = 0; i < n; i++)
      if (vertexes[i].type){
         if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-0.5) < EPSA &&
                  vertexes[i].x[1] < 0.5+EPSA)  vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.5) < EPSA &&
                  vertexes[i].x[0] < 0.5+EPSA)  vertexes[i].type |= NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++)
         if (i >= nvx2 || j >= nvy2){
            if (i < nvx2)
               il = j - nvy2 + (nvy-nvy2)*i;
            else
               il = j + (nvy-nvy2)*nvx2 + nvy*(i-nvx2);
            if (i < nvx2-1)
               nny = nvy - nvy2;
            else
               nny = nvy;
            if (type == 1){
               elem(il, il+nny, il+nny+1, nodes, pelement); 
               elem(il, il+1,   il+nny+1, nodes, pelement);    
            }
            else if (type == 2){
               elem(il,   il+1,   il+nny,   nodes, pelement);    
               elem(il+1, il+nny, il+nny+1, nodes, pelement); 
            }
            else if (type == 3){
               if (2*((i+j)/2) == i+j){
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
            }
            else if (type == 4){
               if (2*((i+j)/2) == i+j){
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
            }
            else if (type == 5){
               if (2*(j/2) == j){
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
            }
            else if (type == 6){
               if (2*(j/2) == j){
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
            }
         }
   (*pelement - 1)->succ = NULL;      
}

void step2_xy(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n=0, il, nvx=NVX, nvy=NVY, nvx2, nvy2, nny, type=CUBE_TYPE;
   NODE *nodes;                        /*  type = 1, 2 ... Friedrichs-Keller  */
   VERTEX *vertexes;
 
   if (2*(nvx/2) == nvx)
      nvx += 1;
   if (2*(nvy/2) == nvy)
      nvy += 1;
   nvx2 = nvx/2;
   nvy2 = nvy/2;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) 
         if (i <= nvx2 || j >= nvy2){
            (*pvert)->x[0] = i/(nvx-1.);
            (*pvert)->x[1] = j/(nvy-1.);
            (*pvert)->type = 0; 
            nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
            (*pnode)++;
         }
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < n; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-0.5) < EPSA &&
               vertexes[i].x[1] < 0.5+EPSA)  vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.5) < EPSA &&
               vertexes[i].x[0] > 0.5-EPSA)  vertexes[i].type = 1;
   }
   for (i = 0; i < n; i++)
      if (vertexes[i].type){
         if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-0.5) < EPSA &&
                  vertexes[i].x[1] < 0.5+EPSA)  vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.5) < EPSA &&
                  vertexes[i].x[0] > 0.5-EPSA)  vertexes[i].type |= NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++)
         if (i < nvx2 || j >= nvy2){
            if (i <= nvx2)
               il = j + nvy*i;
            else
               il = j - nvy2 + nvy*(nvx2+1) + (nvy-nvy2)*(i-nvx2-1);
            if (i < nvx2)
               nny = nvy;
            else
               nny = nvy - nvy2;
            if (type == 1){
               elem(il, il+nny, il+nny+1, nodes, pelement); 
               elem(il, il+1,   il+nny+1, nodes, pelement);    
            }
            else if (type == 2){
               elem(il,   il+1,   il+nny,   nodes, pelement);    
               elem(il+1, il+nny, il+nny+1, nodes, pelement); 
            }
            else if (type == 3){
               if (2*((i+j)/2) == i+j){
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
            }
            else if (type == 4){
               if (2*((i+j)/2) == i+j){
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
            }
            else if (type == 5){
               if (2*(j/2) == j){
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
               else{
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
            }
            else if (type == 6){
               if (2*(j/2) == j){
                  elem(il,   il+1,   il+nny,   nodes, pelement);    
                  elem(il+1, il+nny, il+nny+1, nodes, pelement); 
               }
               else{
                  elem(il, il+nny, il+nny+1, nodes, pelement); 
                  elem(il, il+1,   il+nny+1, nodes, pelement);    
               }
            }
         }
   (*pelement - 1)->succ = NULL;      
}

void cube_for_donut(pnode,pvert,pelement,pgrid) /*  cube = square  */
VERTEX **pvert;                                 /*  NV has to be odd and > 4  */
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   ELEMENT *pel;
   NODE *pn, *nodes;
   INT i, n, nv; 
 
   cube(pnode,pvert,pelement,pgrid);
   for (nodes = FIRSTNODE(pgrid); nodes; nodes = nodes->succ)
      if (fabs(MYVERTEX(nodes)->x[0]-0.5) < 1.e-10 && 
               MYVERTEX(nodes)->x[1] < 0.500001)
         MYVERTEX(nodes)->type = 1|NMASK;
   for (nodes = FIRSTNODE(pgrid); nodes->succ; nodes = nodes->succ); 
   nodes = nodes->succ = *pnode; 
   nv = (NV-3)/2;
   n = pgrid->maxNodeIndex; 
   for (i = 1; i <= nv; i++){
      (*pvert)->x[0] = 0.5;
      (*pvert)->x[1] = i/(NV-1.);
      (*pvert)->type = 1; 
      nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   }
   (*pnode-1)->succ = NULL;
   pgrid->maxNodeIndex = n; 
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
      if(IS_IN_SQUARE(pel,0.500001,0.500001))
         for (i = 0; i < 3; i++)
            if (fabs(MYVERTEX(pel->n[i])->x[0]-0.5) < 1.e-10 && 
                     MYVERTEX(pel->n[i])->x[1] > 0.000001 &&
                     MYVERTEX(pel->n[i])->x[1] < 0.499999){
               for (pn = nodes; fabs(MYVERTEX(pn)->x[1]-
                            MYVERTEX(pel->n[i])->x[1]) > 1.e-10; pn = pn->succ);
               pel->n[i] = pn;
            }
}

void cube1(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n=0, il, nv=NV;
   DOUBLE d;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i <= nv; i++)
      for (j = 0; j < nv; j++) {
         if ((2*(j/2) != j) && i)
            d = 0.5;
         else
            d = 0.;
         if (i < nv){
            (*pvert)->x[0] = (i-d)/(nv-1.0);
            (*pvert)->x[1] = j/(nv-1.0);
            (*pvert)->type = 0; 
            nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
            (*pnode)++;
         }
         else if (d > 0.){
            (*pvert)->x[0] = 1.0;
            (*pvert)->x[1] = j/(nv-1.0);
            (*pvert)->type = 0; 
            nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
            (*pnode)++;
         }
      } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < n; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1;
   }
   for (i = 0; i < n; i++)
      if (vertexes[i].type){
         if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= NMASK;
         if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= NMASK;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++){
            il = j + nv*i;
            if (2*(j/2) == j){
               elem(il, il+nv, il+nv+1, nodes, pelement); 
               elem(il, il+1,  il+nv+1, nodes, pelement);    
            }
            else{
               elem(il,   il+1,  il+nv,   nodes, pelement);    
               elem(il+1, il+nv, il+nv+1, nodes, pelement); 
            }
      }
   il = (nv-1)*nv;
   i = nv*nv;
   n = 0;
   for (j = 0; j < nv-1; j++){
      elem(il, il+1, i, nodes, pelement);
      il++;
      if (n){
         i++;
         n = 0;
      }
      else
         n = 1;
   }
   (*pelement - 1)->succ = NULL;      
}

void cube2(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 13; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 13;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.5;
   vertexes[1].x[1] = 0.; 
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 1.;
   vertexes[3].x[1] = 0.5;
   vertexes[4].x[0] = 1.;
   vertexes[4].x[1] = 1.;
   vertexes[5].x[0] = 0.5;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.;
   vertexes[6].x[1] = 1.;
   vertexes[7].x[0] = 0.;
   vertexes[7].x[1] = 0.5;
   vertexes[8].x[0] = 0.5;
   vertexes[8].x[1] = 0.25;
   vertexes[9].x[0] = 0.75;
   vertexes[9].x[1] = 0.5;
   vertexes[10].x[0] = 0.5;
   vertexes[10].x[1] = 0.75;
   vertexes[11].x[0] = 0.25;
   vertexes[11].x[1] = 0.5;
   vertexes[12].x[0] = 0.5;
   vertexes[12].x[1] = 0.5;
   for (i = 0; i < 8; i++)
      vertexes[i].type = 1|NMASK; 
   for (i = 8; i < 13; i++)
      vertexes[i].type = 0; 
   for (i = 0; i < 13; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 8,nodes,pelement);
   elem( 1, 2, 8,nodes,pelement);
   elem( 2, 3, 9,nodes,pelement);
   elem( 3, 4, 9,nodes,pelement);
   elem( 4, 5,10,nodes,pelement);
   elem( 5, 6,10,nodes,pelement);
   elem( 6, 7,11,nodes,pelement);
   elem( 0, 7,11,nodes,pelement);
   elem( 0, 8,11,nodes,pelement);
   elem( 2, 8, 9,nodes,pelement);
   elem( 4, 9,10,nodes,pelement);
   elem( 6,10,11,nodes,pelement);
   elem( 8,11,12,nodes,pelement);
   elem( 8, 9,12,nodes,pelement);
   elem( 9,10,12,nodes,pelement);
   elem(10,11,12,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void cube3(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, im, nv, na; 
   FLOAT s;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   s = 0.5/(nv-1.0);
   for (i = 0; i < nv-1; i++) 
      for (j = 0; j < nv-1; j++) {
           (*pvert)->x[0] = s + i/(nv-1.0);
           (*pvert)->x[1] = s + j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            im = na + j + (nv-1)*i;
            elem(il,    il+nv,   im, nodes, pelement); 
            elem(il+nv, il+nv+1, im, nodes, pelement); 
            elem(il+1,  il+nv+1, im, nodes, pelement); 
            elem(il,    il+1,    im, nodes, pelement); 
   }
   (*pelement - 1)->succ = NULL;      
}

void cube4(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 9; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 9;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.5;
   vertexes[1].x[1] = 0.; 
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 1.;
   vertexes[3].x[1] = 0.5;
   vertexes[4].x[0] = 1.;
   vertexes[4].x[1] = 1.;
   vertexes[5].x[0] = 0.5;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.;
   vertexes[6].x[1] = 1.;
   vertexes[7].x[0] = 0.;
   vertexes[7].x[1] = 0.5;
   vertexes[8].x[0] = 0.5;
   vertexes[8].x[1] = 0.5;
   for (i = 0; i < 8; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[8].type = 0; 
   for (i = 0; i < 9; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 8,nodes,pelement);
   elem( 1, 2, 8,nodes,pelement);
   elem( 2, 3, 8,nodes,pelement);
   elem( 3, 4, 8,nodes,pelement);
   elem( 4, 5, 8,nodes,pelement);
   elem( 5, 6, 8,nodes,pelement);
   elem( 6, 7, 8,nodes,pelement);
   elem( 0, 7, 8,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void cube5(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 9; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 9;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.5;
   vertexes[1].x[1] = 0.; 
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 1.;
   vertexes[3].x[1] = 0.5;
   vertexes[4].x[0] = 1.;
   vertexes[4].x[1] = 1.;
   vertexes[5].x[0] = 0.5;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.;
   vertexes[6].x[1] = 1.;
   vertexes[7].x[0] = 0.;
   vertexes[7].x[1] = 0.5;
   vertexes[8].x[0] = 0.5;
   vertexes[8].x[1] = 0.5;
   for (i = 0; i < 8; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[8].type = 0; 
   for (i = 0; i < 9; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 1, 2, 3,nodes,pelement);
   elem( 3, 4, 5,nodes,pelement);
   elem( 5, 6, 7,nodes,pelement);
   elem( 0, 1, 7,nodes,pelement);
   elem( 1, 3, 8,nodes,pelement);
   elem( 3, 5, 8,nodes,pelement);
   elem( 5, 7, 8,nodes,pelement);
   elem( 1, 7, 8,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void cube6(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 5; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 5;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 1.;
   vertexes[1].x[1] = 0.;
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 1.;
   vertexes[3].x[0] = 0.;
   vertexes[3].x[1] = 1.;
   vertexes[4].x[0] = 0.5;
   vertexes[4].x[1] = 0.5;
   for (i = 0; i < 4; i++)
      vertexes[i].type = 1|NMASK;
   vertexes[4].type = 0; 
   for (i = 0; i < 5; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 4,nodes,pelement);
   elem( 1, 2, 4,nodes,pelement);
   elem( 2, 3, 4,nodes,pelement);
   elem( 0, 3, 4,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void rectangle_2x1(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nvx, nvy, na; 
   FLOAT hx, hy;
   NODE *nodes;
   VERTEX *vertexes;
 
   nvx = 2*NV-1;           /* number of vertices in x-direction */
   nvy = NV;               /* number of vertices in y-direction */
   na  = nvx*nvy;          /* number of all vertices */
   hx = 2./(nvx-1.);
   hy = 1./(nvy-1.);
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           (*pvert)->x[0] = i*hx;
           (*pvert)->x[1] = j*hy;
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-2.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++) {
            il = j + nvy*i;
            elem(il, il+nvy, il+nvy+1, nodes, pelement); 
            elem(il, il+1,   il+nvy+1, nodes, pelement);    
   }
   (*pelement - 1)->succ = NULL;      
}

void rectangle3_2x1(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, im, nvx, nvy, na; 
   FLOAT hx, hy, sx, sy;
   NODE *nodes;
   VERTEX *vertexes;
 
   nvx = 2*NV-1;           /* number of vertices in x-direction */
   nvy = NV;               /* number of vertices in y-direction */
   na  = nvx*nvy;
   hx = 2./(nvx-1.);
   hy = 1./(nvy-1.);
   sx = 0.5*hx;
   sy = 0.5*hy;
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           (*pvert)->x[0] = i*hx;
           (*pvert)->x[1] = j*hy;
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   for (i = 0; i < nvx-1; i++) 
      for (j = 0; j < nvy-1; j++) {
           (*pvert)->x[0] = sx + i*hx;
           (*pvert)->x[1] = sy + j*hy;
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-2.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++) {
            il = j + nvy*i;
            im = na + j + (nvy-1)*i;
            elem(il,     il+nvy,   im, nodes, pelement); 
            elem(il+nvy, il+nvy+1, im, nodes, pelement); 
            elem(il+1,   il+nvy+1, im, nodes, pelement); 
            elem(il,     il+1,     im, nodes, pelement); 
   }
   (*pelement - 1)->succ = NULL;      
}

void rectangle(pnode,pvert,pelement,pgrid) /*  rectangle dx X dy, lower left  */
VERTEX **pvert;                            /*  vertex is (lx,ly)              */
NODE **pnode;                              /*  NV point in both directions    */
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na; 
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT lx=LOWER_LEFT_X, ly=LOWER_LEFT_Y, dx=LENGTH_X, dy=LENGTH_Y;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) {
           (*pvert)->x[0] = dx*i/(nv-1.0) + lx;
           (*pvert)->x[1] = dy*j/(nv-1.0) + ly;
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-lx) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-ly) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-lx-dx) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-ly-dy) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            elem(il, il+nv, il+nv+1, nodes, pelement); 
            elem(il, il+1,  il+nv+1, nodes, pelement);    
   }
   (*pelement - 1)->succ = NULL;      
}

/* creates stair shape triangulation */
void stairs(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nvp, nvc, nv, na, dir; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   if (nv < 4) nv = 5;
   if (nv % 2 == 0) nv++;
   nvp = nv/2;
   nvc = nvp/2;
   na = nv*nv - nvp*nvp;  /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (j = 0; j < nv; j++) 
      for (i = 0; i < nv; i++) {
		 if (i > nvp-1 || j < nv-nvp) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
		 }
      } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1 | NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1 | NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1 | NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1 | NMASK_FOR_FB;
	  if (fabs(vertexes[i].x[1]-0.5) < EPSA && vertexes[i].x[0] < 0.5+EPSA) vertexes[i].type |= 1 | NMASK_FOR_FB;
	  if (fabs(vertexes[i].x[0]-0.5) < EPSA && vertexes[i].x[1] > 0.5-EPSA) vertexes[i].type = 1 | NMASK_FOR_FB;
   }
   for (i = 0; i < nvp; i++) {
	  dir = (i < nvc ? 0 : 1);
      for (j = 0; j < nv-1; j++) {
		  if (j == nvc || j == nvp || j == nvp+nvc) dir = (dir ? 0 : 1);
			il = j + nv*i;
			if (dir) {
				elem(il,   il+1,    il+nv, nodes, pelement); 
				elem(il+1, il+nv, il+nv+1, nodes, pelement);
			} else {
				elem(il, il+nv, il+nv+1, nodes, pelement); 
				elem(il, il+1,  il+nv+1, nodes, pelement);
			}
	  }
   }
   n = nvp*(nv+1);
   for (i = 0; i < nvp; i++) {
	   dir = (i < nvc ? 0 : 1);
	   for (j = 0; j < nvp; j++) {
			if (j == nvc) dir = (dir ? 0 : 1);
			il = n + j + (nvp+1)*i;
			if (dir) {
				elem(il,   il+1,     il+nvp+1, nodes, pelement);
				elem(il+1, il+nvp+1, il+nvp+2, nodes, pelement);
			} else {
				elem(il, il+nvp+1, il+nvp+2, nodes, pelement);
				elem(il, il+1,     il+nvp+2, nodes, pelement);
			}
	   }
   }
   (*pelement - 1)->succ = NULL;      
}

void stairs2(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nvp, nvc, nv, na, dir; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   if (nv < 4) nv = 5;
   if (nv % 2 == 0) nv++;
   nvp = nv/2;
   nvc = nvp/2;
   na = nv*nv - nvp*nvp;  /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (j = 0; j < nv; j++) 
      for (i = 0; i < nv; i++) {
		 if (i > nvp-1 || j < nv-nvp) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
		 }
      } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) {
          if (vertexes[i].x[1] < 0.25+EPSA) vertexes[i].type = 1 | NMASK;
          if (vertexes[i].x[1] > 0.25-EPSA) vertexes[i].type |= 1 | NMASK_ZN;
      }
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) {
          if (vertexes[i].x[1] < 0.75+EPSA) vertexes[i].type = 1 | NMASK;
          if (vertexes[i].x[1] > 0.75-EPSA) vertexes[i].type |= 1 | NMASK_ZN;
      }
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1 | NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1 | NMASK_FOR_FB;
	  if (fabs(vertexes[i].x[1]-0.5) < EPSA && vertexes[i].x[0] < 0.5+EPSA) vertexes[i].type |= 1 | NMASK_FOR_FB;
	  if (fabs(vertexes[i].x[0]-0.5) < EPSA && vertexes[i].x[1] > 0.5-EPSA) vertexes[i].type = 1 | NMASK_FOR_FB;
   }
   for (i = 0; i < nvp; i++) {
	  dir = (i < nvc ? 0 : 1);
      for (j = 0; j < nv-1; j++) {
		  if (j == nvc || j == nvp || j == nvp+nvc) dir = (dir ? 0 : 1);
			il = j + nv*i;
			if (dir) {
				elem(il,   il+1,    il+nv, nodes, pelement); 
				elem(il+1, il+nv, il+nv+1, nodes, pelement);
			} else {
				elem(il, il+nv, il+nv+1, nodes, pelement); 
				elem(il, il+1,  il+nv+1, nodes, pelement);
			}
	  }
   }
   n = nvp*(nv+1);
   for (i = 0; i < nvp; i++) {
	   dir = (i < nvc ? 0 : 1);
	   for (j = 0; j < nvp; j++) {
			if (j == nvc) dir = (dir ? 0 : 1);
			il = n + j + (nvp+1)*i;
			if (dir) {
				elem(il,   il+1,     il+nvp+1, nodes, pelement);
				elem(il+1, il+nvp+1, il+nvp+2, nodes, pelement);
			} else {
				elem(il, il+nvp+1, il+nvp+2, nodes, pelement);
				elem(il, il+1,     il+nvp+2, nodes, pelement);
			}
	   }
   }
   (*pelement - 1)->succ = NULL;      
}

void stairs3(pnode,pvert,pelement,pgrid)	//squere with cutted corner
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 8; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 8;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.5;
   vertexes[1].x[1] = 0.; 
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 1.;
   vertexes[3].x[1] = 0.5;
   vertexes[4].x[0] = 1.;
   vertexes[4].x[1] = 1.;
   vertexes[5].x[0] = 0.5;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.;
   vertexes[6].x[1] = 0.5;
   vertexes[7].x[0] = 0.5;
   vertexes[7].x[1] = 0.5;
   for (i = 0; i < 6; i++) {
      vertexes[i].type = 1;
	  if (i < 4)
		  vertexes[i].type |= NMASK;
   }
   for (i = 4; i < 7; i++) vertexes[i].type |= NMASK_FOR_FB;
   vertexes[0].type |= NMASK_ZN;
   vertexes[6].type |= NMASK_ZN;
   vertexes[3].type |= NMASK_ZN;
   vertexes[4].type |= NMASK_ZN;
   vertexes[7].type = 0; 
   for (i = 0; i < 8; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 7,nodes,pelement);
   elem( 1, 2, 7,nodes,pelement);
   elem( 2, 3, 7,nodes,pelement);
   elem( 3, 4, 7,nodes,pelement);
   elem( 4, 5, 7,nodes,pelement);
   elem( 5, 6, 7,nodes,pelement);
   elem( 0, 6, 7,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void cube_1(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = 5;               /* number of vertices in one direction */
   na = 28;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) {
      (*pvert)->x[0] = i/(nv-1.0);
      (*pvert)->x[1] = 0.;
      (*pvert)->type = 1|NMASK; 
      nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   for (j = 2; j < nv; j++)
	  for (i = 0; i < nv; i++) {
         (*pvert)->x[0] = i/(nv-1.0);
         (*pvert)->x[1] = j/(nv-1.0);
         (*pvert)->type = 0; 
         nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
         (*pnode)++;
      }
   for (j = 2; j < 4; j++)
	  for (i = 0; i < 4; i++) {
		 (*pvert)->x[0] = 1./(2*(nv-1.0)) + i/(nv-1.0);
         (*pvert)->x[1] = (1./(2*(nv-1.0)) + j/(nv-1.0));
         (*pvert)->type = 0; 
         nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
         (*pnode)++;
	  }

   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 5; i < 20; i++){
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK_FOR_FB;
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
   }   

   elem(0, 5, 6, nodes, pelement); 
   elem(0, 1, 6, nodes, pelement);
   elem(1, 6, 7, nodes, pelement); 
   elem(1, 2, 7, nodes, pelement);
   elem(2, 3, 7, nodes, pelement); 
   elem(3, 7, 8, nodes, pelement);
   elem(3, 4, 8, nodes, pelement); 
   elem(4, 8, 9, nodes, pelement);
   for (j = 0; j < 2; j++)
      for (i = 0; i < nv-1; i++) {
         il = 5 + i + 5*j;
         elem(il, il+1, il+15-j, nodes, pelement);
         elem(il, il+5, il+15-j, nodes, pelement);
 	     elem(il+1, il+6, il+15-j, nodes, pelement);
		 elem(il+5, il+6, il+15-j, nodes, pelement);
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_2(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na, na1;
   FLOAT y, s, r;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;							/* number of vertices in one direction */
   na1 = nv*nv;
   na = na1 + 2*(nv-1);         /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   y = 1.;
   for (i = 0; i < nv - 1; i++) y = y/2.;

   for (i = 0; i < nv; i++) {
	  s = 0.5; 
	  r = 0.;
	  for (j = 0; j < nv; j++) {
		(*pvert)->x[0] = i/(nv-1.0);
		if (j == 0) (*pvert)->x[1] = 0.;
		else
			if (j == nv-1) (*pvert)->x[1] = 1.;
			else {
				r = r + s;
				(*pvert)->x[1] = r;
				s = s/2.;
			}
	    (*pvert)->type = 0; 
		nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
		(*pnode)++;
	  } 
   }
   for (i = 0; i < nv-1; i++)
	  for (j = 0; j < 2; j++) {
		 (*pvert)->x[0] = 1./(2*(nv-1.0)) + i/(nv-1.0);
         (*pvert)->x[1] = 1. - y + j*y;
         (*pvert)->type = 0; 
         nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
         (*pnode)++;
	  }
/*   (*pvert)->x[0] = 0.;
   (*pvert)->x[1] = 1. - y;
   (*pvert)->type = 0; 
   nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
   (*pnode)++;
   (*pvert)->x[0] = 1.;
   (*pvert)->x[1] = 1. - y;
   (*pvert)->type = 0; 
   nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
   (*pnode)++;
*/
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
	  if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_FOR_FB;
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
	  printf("%i %e %e %i\n", i, vertexes[i].x[0], vertexes[i].x[1], vertexes[i].type);
   }   

   for (j = 0; j < nv-2; j++) {
       for (i = 0; i < nv/2; i++) {
		   il = nv*i + j;
		   elem(il,  il+1, il+nv+1, nodes, pelement);
		   elem(il, il+nv, il+nv+1, nodes, pelement);
	   }
	   for (i = nv/2; i < nv-1; i++) {
		   il = nv*i + j;
		   elem(  il,  il+1,   il+nv, nodes, pelement);
		   elem(il+1, il+nv, il+nv+1, nodes, pelement);
	   }
   }
   j = nv - 2;
   for (i = 0; i < nv-1; i++) {
      il = nv*i + j;
      elem(il, il+nv, na1+2*i, nodes, pelement);
      //if (i > 0)
		  elem(il, il+1, na1+2*i, nodes, pelement);
	  /*else {
		  elem(  il, na1, na-2, nodes, pelement);
		  elem(il+1, na1, na-2, nodes, pelement);
	  }
	  if (i == nv-2) {
		  elem(  il+nv, na1+2*i, na-1, nodes, pelement);
		  elem(il+nv+1, na1+2*i, na-1, nodes, pelement);
	  } else {*/
		  elem(il+nv, il+nv+1, na1+2*i, nodes, pelement);
	  //}
      elem(   il+1, na1+2*i, na1+2*i+1, nodes, pelement);
      elem(il+nv+1, na1+2*i, na1+2*i+1, nodes, pelement);
   }
   (*pelement - 1)->succ = NULL;      
}

INT distrib(k, n)
INT k, n;
{
	INT i, kf, nf;

	if (k == 0) return 0;
	if (k > n) k = k - n;
	if (k == 1) return 1;
	for (i = 2, kf = 1; i < k; i++) kf = i*kf;
	for (i = 2, nf = n-1; i < k; i++) nf = nf*(n-i);
	return nf/kf;
}

void cube_3(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, nv1, na, na1, p;
   FLOAT y, r;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;							/* number of vertices in one direction */
   nv1 = nv + nv/2;
   na1 = nv*nv1;
   na = na1 + nv1 - 1;         /* number of all vertices */
   p = nv1/2;

   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   y = 1.;
   for (i = 0; i < p; i++) y = y/2.;

   r = 0.;
   for (i = 0; i < nv1; i++) {
	  r += y*distrib(i, p);
	  for (j = 0; j < nv; j++) {
		(*pvert)->x[0] = r;
		(*pvert)->x[1] = j/(nv-1.0);
		(*pvert)->type = 0; 
		nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
		(*pnode)++;
	  } 
   }
 /*  for (i = 0; i < nv-1; i++)
	  for (j = 0; j < 2; j++) {
		 (*pvert)->x[0] = 1./(2*(nv-1.0)) + i/(nv-1.0);
         (*pvert)->x[1] = 1. - y + j*y;
         (*pvert)->type = 0; 
         nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
         (*pnode)++;
	  }
/*   (*pvert)->x[0] = 0.;
   (*pvert)->x[1] = 1. - y;
   (*pvert)->type = 0; 
   nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
   (*pnode)++;
   (*pvert)->x[0] = 1.;
   (*pvert)->x[1] = 1. - y;
   (*pvert)->type = 0; 
   nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
   (*pnode)++;
*/
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na1; i++){
	  if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_FOR_FB;
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
	  //printf("%i %e %e %i\n", i, vertexes[i].x[0], vertexes[i].x[1], vertexes[i].type);
   }   

   for (j = 0; j < nv-2; j++) {
       for (i = 0; i < p; i++) {
		   il = nv*i + j;
		   elem(il,  il+1, il+nv+1, nodes, pelement);
		   elem(il, il+nv, il+nv+1, nodes, pelement);
	   }
	   for (i = p; i < nv1-1; i++) {
		   il = nv*i + j;
		   elem(  il,  il+1,   il+nv, nodes, pelement);
		   elem(il+1, il+nv, il+nv+1, nodes, pelement);
	   }
   }
   j = nv-2;
   for (i = 0; i < p; i++) {
	   il = nv*i + j;
	   elem(  il,  il+1,   il+nv, nodes, pelement);
	   elem(il+1, il+nv, il+nv+1, nodes, pelement);
   }
   for (i = p; i < nv1-1; i++) {
	   il = nv*i + j;
	   elem(il,  il+1, il+nv+1, nodes, pelement);
	   elem(il, il+nv, il+nv+1, nodes, pelement);
   }
 /*  j = nv - 2;
   for (i = 0; i < nv-1; i++) {
      il = nv*i + j;
      elem(il, il+nv, na1+2*i, nodes, pelement);
		  elem(il, il+1, na1+2*i, nodes, pelement);
		  elem(il+nv, il+nv+1, na1+2*i, nodes, pelement);
      elem(   il+1, na1+2*i, na1+2*i+1, nodes, pelement);
      elem(il+nv+1, na1+2*i, na1+2*i+1, nodes, pelement);
   }*/
   (*pelement - 1)->succ = NULL;      
}

void cube_3a(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, nv1, na, na1, p;
   FLOAT x, y, r, s, t;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;							/* number of vertices in one direction */
   nv1 = nv + nv;
   if (nv1 % 2 == 0) nv1++;
   na1 = nv*nv1;
   na = na1 + nv1 - 1;         /* number of all vertices */
   p = nv1/2;

   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;

   x = 1.;
   for (i = 0; i < nv - 1; i++) x = x/2.;

   y = 1.;
   for (i = 0; i < p; i++) y = y/2.;

   r = 0.;
   for (i = 0; i < nv1; i++) {
	  r += y*distrib(i, p);
	  //r = i/(nv1-1.);
	  s = 0.5; t = 0.;
	  for (j = 0; j < nv; j++) {
		(*pvert)->x[0] = r;
		if (j == 0) (*pvert)->x[1] = 0.;
		else
			if (j == nv-1) (*pvert)->x[1] = SURFACE;
			else {
				t = t + s;
				(*pvert)->x[1] = SURFACE*t;
				s = s/2.;
			}
		(*pvert)->type = 0; 
		nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
		(*pnode)++;
	  } 
   }
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na1; i++){
	  if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= 1|NMASK;
      if (fabs(vertexes[i].x[1]-SURFACE) < EPSA) vertexes[i].type |= 1|NMASK_FOR_FB;
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
	  //printf("%i %e %e %i\n", i, vertexes[i].x[0], vertexes[i].x[1], vertexes[i].type);
   }   

   for (j = 0; j < nv-2; j++) {
       for (i = 0; i < p; i++) {
		   il = nv*i + j;
		   elem(il,  il+1, il+nv+1, nodes, pelement);
		   elem(il, il+nv, il+nv+1, nodes, pelement);
	   }
	   for (i = p; i < nv1-1; i++) {
		   il = nv*i + j;
		   elem(  il,  il+1,   il+nv, nodes, pelement);
		   elem(il+1, il+nv, il+nv+1, nodes, pelement);
	   }
   }
   j = nv-2;
   for (i = 0; i < p; i++) {
	   il = nv*i + j;
	   elem(  il,  il+1,   il+nv, nodes, pelement);
	   elem(il+1, il+nv, il+nv+1, nodes, pelement);
   }
   for (i = p; i < nv1-1; i++) {
	   il = nv*i + j;
	   elem(il,  il+1, il+nv+1, nodes, pelement);
	   elem(il, il+nv, il+nv+1, nodes, pelement);
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_3c(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, nv1, na, p;
   FLOAT y, r;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;							/* number of vertices in one direction */
   nv1 = nv + nv/2;
   nv = nv/2+1;
   na = nv*nv1;         /* number of all vertices */
   p = nv1/2;

   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   y = 1.;
   for (i = 0; i < p; i++) y = y/2.;

   r = 0.;
   for (i = 0; i < nv1; i++) {
	  r += y*distrib(i, p);
	  for (j = 0; j < nv; j++) {
		(*pvert)->x[0] = r;
		(*pvert)->x[1] = j/(nv-1.0);
		(*pvert)->type = 0; 
		nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
		(*pnode)++;
	  } 
   }

   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
	  if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_FOR_FB;
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
	  //printf("%i %e %e %i\n", i, vertexes[i].x[0], vertexes[i].x[1], vertexes[i].type);
   }   

   for (j = 0; j < nv-2; j++) {
       for (i = 0; i < p; i++) {
		   il = nv*i + j;
		   elem(il,  il+1, il+nv+1, nodes, pelement);
		   elem(il, il+nv, il+nv+1, nodes, pelement);
	   }
	   for (i = p; i < nv1-1; i++) {
		   il = nv*i + j;
		   elem(  il,  il+1,   il+nv, nodes, pelement);
		   elem(il+1, il+nv, il+nv+1, nodes, pelement);
	   }
   }
   j = nv-2;
   for (i = 0; i < p; i++) {
	   il = nv*i + j;
	   elem(  il,  il+1,   il+nv, nodes, pelement);
	   elem(il+1, il+nv, il+nv+1, nodes, pelement);
   }
   for (i = p; i < nv1-1; i++) {
	   il = nv*i + j;
	   elem(il,  il+1, il+nv+1, nodes, pelement);
	   elem(il, il+nv, il+nv+1, nodes, pelement);
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_4(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_ZN;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type |= 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type |= 1|NMASK_FOR_FB;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            elem(il, il+nv, il+nv+1, nodes, pelement); 
            elem(il, il+1,  il+nv+1, nodes, pelement);    
   }
   (*pelement - 1)->succ = NULL;      
}

void rectangle_with_hollow(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 11; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 11;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.5;
   vertexes[1].x[0] = 0.;
   vertexes[1].x[1] = 1.; 
   vertexes[2].x[0] = 0.5;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 0.5;
   vertexes[3].x[1] = 0.5;
   vertexes[4].x[0] = 0.5;
   vertexes[4].x[1] = 1.;
   vertexes[5].x[0] = 1.;
   vertexes[5].x[1] = 0.;
   vertexes[6].x[0] = 1.;
   vertexes[6].x[1] = 0.5;
   vertexes[7].x[0] = 1.;
   vertexes[7].x[1] = 1.;
   vertexes[8].x[0] = 1.5;
   vertexes[8].x[1] = 0.5;
   vertexes[9].x[0] = 1.5;
   vertexes[9].x[1] = 1.;
   vertexes[10].x[0] = 0.75;
   vertexes[10].x[1] = 0.25;
   for (i = 0; i < 10; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[10].type = 0; 
   for (i = 0; i < 11; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 4,nodes,pelement);
   elem( 0, 3, 4,nodes,pelement);
   elem( 3, 4, 7,nodes,pelement);
   elem( 3, 6, 7,nodes,pelement);
   elem( 6, 7, 9,nodes,pelement);
   elem( 6, 8, 9,nodes,pelement);
   elem( 2, 3, 10,nodes,pelement);
   elem( 2, 5, 10,nodes,pelement);
   elem( 3, 6, 10,nodes,pelement);
   elem( 5, 6, 10,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void cube_with_hole(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, hole_type;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 40; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 40;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.;
   vertexes[1].x[1] = 0.2;
   vertexes[2].x[0] = 0.;
   vertexes[2].x[1] = 0.4;
   vertexes[3].x[0] = 0.;
   vertexes[3].x[1] = 0.6;
   vertexes[4].x[0] = 0.;
   vertexes[4].x[1] = 0.8;
   vertexes[5].x[0] = 0.;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.2;
   vertexes[6].x[1] = 0.;
   vertexes[7].x[0] = 0.2;
   vertexes[7].x[1] = 0.2;
   vertexes[8].x[0] = 0.2;
   vertexes[8].x[1] = 0.4;
   vertexes[9].x[0] = 0.2;
   vertexes[9].x[1] = 0.6;
   vertexes[10].x[0] = 0.2;
   vertexes[10].x[1] = 0.8;
   vertexes[11].x[0] = 0.2;
   vertexes[11].x[1] = 1.;
   vertexes[12].x[0] = 0.4;
   vertexes[12].x[1] = 0.;
   vertexes[13].x[0] = 0.4;
   vertexes[13].x[1] = 0.2;
   vertexes[14].x[0] = 0.4;
   vertexes[14].x[1] = 0.4;
   vertexes[15].x[0] = 0.3;
   vertexes[15].x[1] = 0.5;
   vertexes[16].x[0] = 0.4;
   vertexes[16].x[1] = 0.6;
   vertexes[17].x[0] = 0.4;
   vertexes[17].x[1] = 0.8;
   vertexes[18].x[0] = 0.4;
   vertexes[18].x[1] = 1.;
   vertexes[19].x[0] = 0.5;
   vertexes[19].x[1] = 0.3;
   vertexes[20].x[0] = 0.5;
   vertexes[20].x[1] = 0.7;
   vertexes[21].x[0] = 0.6;
   vertexes[21].x[1] = 0.;
   vertexes[22].x[0] = 0.6;
   vertexes[22].x[1] = 0.2;
   vertexes[23].x[0] = 0.6;
   vertexes[23].x[1] = 0.4;
   vertexes[24].x[0] = 0.7;
   vertexes[24].x[1] = 0.5;
   vertexes[25].x[0] = 0.6;
   vertexes[25].x[1] = 0.6;
   vertexes[26].x[0] = 0.6;
   vertexes[26].x[1] = 0.8;
   vertexes[27].x[0] = 0.6;
   vertexes[27].x[1] = 1.;
   vertexes[28].x[0] = 0.8;
   vertexes[28].x[1] = 0.;
   vertexes[29].x[0] = 0.8;
   vertexes[29].x[1] = 0.2;
   vertexes[30].x[0] = 0.8;
   vertexes[30].x[1] = 0.4;
   vertexes[31].x[0] = 0.8;
   vertexes[31].x[1] = 0.6;
   vertexes[32].x[0] = 0.8;
   vertexes[32].x[1] = 0.8;
   vertexes[33].x[0] = 0.8;
   vertexes[33].x[1] = 1.;
   vertexes[34].x[0] = 1.;
   vertexes[34].x[1] = 0.;
   vertexes[35].x[0] = 1.;
   vertexes[35].x[1] = 0.2;
   vertexes[36].x[0] = 1.;
   vertexes[36].x[1] = 0.4;
   vertexes[37].x[0] = 1.;
   vertexes[37].x[1] = 0.6;
   vertexes[38].x[0] = 1.;
   vertexes[38].x[1] = 0.8;
   vertexes[39].x[0] = 1.;
   vertexes[39].x[1] = 1.;
   for (i = 7; i < 33; i++)
      vertexes[i].type = 0; 
   for (i = 0; i < 7; i++)
      vertexes[i].type = 1|NMASK; 
   for (i = 33; i < 40; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[11].type = 1|NMASK; 
   vertexes[12].type = 1|NMASK; 
   vertexes[18].type = 1|NMASK; 
   vertexes[21].type = 1|NMASK; 
   vertexes[27].type = 1|NMASK; 
   vertexes[28].type = 1|NMASK; 

   hole_type = 5;
   vertexes[14].type = hole_type|NMASK; 
   vertexes[15].type = hole_type|NMASK; 
   vertexes[16].type = hole_type|NMASK; 
   vertexes[19].type = hole_type|NMASK; 
   vertexes[20].type = hole_type|NMASK; 
   vertexes[23].type = hole_type|NMASK; 
   vertexes[24].type = hole_type|NMASK; 
   vertexes[25].type = hole_type|NMASK; 
   for (i = 0; i < 40; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 6,nodes,pelement);
   elem( 1, 6, 7,nodes,pelement);
   elem( 1, 2, 7,nodes,pelement);
   elem( 2, 7, 8,nodes,pelement);
   elem( 2, 3, 8,nodes,pelement);
   elem( 3, 8, 9,nodes,pelement);
   elem( 3, 4, 9,nodes,pelement);
   elem( 4, 9,10,nodes,pelement);
   elem( 4, 5,11,nodes,pelement);
   elem( 4,10,11,nodes,pelement);
   elem( 6, 7,13,nodes,pelement);
   elem( 6,12,13,nodes,pelement);
   elem( 7, 8,13,nodes,pelement);
   elem( 8,13,14,nodes,pelement);
   elem( 8,14,15,nodes,pelement);
   elem( 8, 9,15,nodes,pelement);
   elem( 9,15,16,nodes,pelement);
   elem( 9,10,17,nodes,pelement);
   elem( 9,16,17,nodes,pelement);
   elem(10,11,18,nodes,pelement);
   elem(10,17,18,nodes,pelement);
   elem(12,13,22,nodes,pelement);
   elem(12,21,22,nodes,pelement);
   elem(13,14,19,nodes,pelement);
   elem(13,19,22,nodes,pelement);
   elem(19,22,23,nodes,pelement);
   elem(16,17,20,nodes,pelement);
   elem(17,20,26,nodes,pelement);
   elem(20,25,26,nodes,pelement);
   elem(17,18,27,nodes,pelement);
   elem(17,26,27,nodes,pelement);
   elem(21,22,29,nodes,pelement);
   elem(21,28,29,nodes,pelement);
   elem(22,23,30,nodes,pelement);
   elem(22,29,30,nodes,pelement);
   elem(23,24,30,nodes,pelement);
   elem(24,30,31,nodes,pelement);
   elem(24,25,31,nodes,pelement);
   elem(25,26,31,nodes,pelement);
   elem(26,31,32,nodes,pelement);
   elem(26,27,33,nodes,pelement);
   elem(26,32,33,nodes,pelement);
   elem(28,29,35,nodes,pelement);
   elem(28,34,35,nodes,pelement);
   elem(29,30,35,nodes,pelement);
   elem(30,35,36,nodes,pelement);
   elem(30,31,36,nodes,pelement);
   elem(31,36,37,nodes,pelement);
   elem(31,32,37,nodes,pelement);
   elem(32,37,38,nodes,pelement);
   elem(32,33,38,nodes,pelement);
   elem(33,38,39,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void cube_with_magnet(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, hole_type;
   NODE *nodes, *pn;
   VERTEX *vertexes;
   ELEMENT *pel;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 41; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 41;
 
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.;
   vertexes[1].x[1] = 0.2;
   vertexes[2].x[0] = 0.;
   vertexes[2].x[1] = 0.4;
   vertexes[3].x[0] = 0.;
   vertexes[3].x[1] = 0.6;
   vertexes[4].x[0] = 0.;
   vertexes[4].x[1] = 0.8;
   vertexes[5].x[0] = 0.;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 0.2;
   vertexes[6].x[1] = 0.;
   vertexes[7].x[0] = 0.2;
   vertexes[7].x[1] = 0.2;
   vertexes[8].x[0] = 0.2;
   vertexes[8].x[1] = 0.4;
   vertexes[9].x[0] = 0.2;
   vertexes[9].x[1] = 0.6;
   vertexes[10].x[0] = 0.2;
   vertexes[10].x[1] = 0.8;
   vertexes[11].x[0] = 0.2;
   vertexes[11].x[1] = 1.;
   vertexes[12].x[0] = 0.4;
   vertexes[12].x[1] = 0.;
   vertexes[13].x[0] = 0.4;
   vertexes[13].x[1] = 0.2;
   vertexes[14].x[0] = 0.4;
   vertexes[14].x[1] = 0.4;
   vertexes[15].x[0] = 0.3;
   vertexes[15].x[1] = 0.5;
   vertexes[16].x[0] = 0.4;
   vertexes[16].x[1] = 0.6;
   vertexes[17].x[0] = 0.4;
   vertexes[17].x[1] = 0.8;
   vertexes[18].x[0] = 0.4;
   vertexes[18].x[1] = 1.;
   vertexes[19].x[0] = 0.5;
   vertexes[19].x[1] = 0.3;
   vertexes[20].x[0] = 0.5;
   vertexes[20].x[1] = 0.7;
   vertexes[21].x[0] = 0.6;
   vertexes[21].x[1] = 0.;
   vertexes[22].x[0] = 0.6;
   vertexes[22].x[1] = 0.2;
   vertexes[23].x[0] = 0.6;
   vertexes[23].x[1] = 0.4;
   vertexes[24].x[0] = 0.7;
   vertexes[24].x[1] = 0.5;
   vertexes[25].x[0] = 0.6;
   vertexes[25].x[1] = 0.6;
   vertexes[26].x[0] = 0.6;
   vertexes[26].x[1] = 0.8;
   vertexes[27].x[0] = 0.6;
   vertexes[27].x[1] = 1.;
   vertexes[28].x[0] = 0.8;
   vertexes[28].x[1] = 0.;
   vertexes[29].x[0] = 0.8;
   vertexes[29].x[1] = 0.2;
   vertexes[30].x[0] = 0.8;
   vertexes[30].x[1] = 0.4;
   vertexes[31].x[0] = 0.8;
   vertexes[31].x[1] = 0.6;
   vertexes[32].x[0] = 0.8;
   vertexes[32].x[1] = 0.8;
   vertexes[33].x[0] = 0.8;
   vertexes[33].x[1] = 1.;
   vertexes[34].x[0] = 1.;
   vertexes[34].x[1] = 0.;
   vertexes[35].x[0] = 1.;
   vertexes[35].x[1] = 0.2;
   vertexes[36].x[0] = 1.;
   vertexes[36].x[1] = 0.4;
   vertexes[37].x[0] = 1.;
   vertexes[37].x[1] = 0.6;
   vertexes[38].x[0] = 1.;
   vertexes[38].x[1] = 0.8;
   vertexes[39].x[0] = 1.;
   vertexes[39].x[1] = 1.;
   vertexes[40].x[0] = 0.5;
   vertexes[40].x[1] = 0.5;
   for (i = 7; i < 33; i++)
      vertexes[i].type = 0; 
   vertexes[40].type = 0;
   for (i = 0; i < 7; i++)
      vertexes[i].type = 1|NMASK; 
   for (i = 33; i < 40; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[11].type = 1|NMASK; 
   vertexes[12].type = 1|NMASK; 
   vertexes[18].type = 1|NMASK; 
   vertexes[21].type = 1|NMASK; 
   vertexes[27].type = 1|NMASK; 
   vertexes[28].type = 1|NMASK; 
   for (i = 0; i < 41; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 6,nodes,pelement);
   elem( 1, 6, 7,nodes,pelement);
   elem( 1, 2, 7,nodes,pelement);
   elem( 2, 7, 8,nodes,pelement);
   elem( 2, 3, 8,nodes,pelement);
   elem( 3, 8, 9,nodes,pelement);
   elem( 3, 4, 9,nodes,pelement);
   elem( 4, 9,10,nodes,pelement);
   elem( 4, 5,11,nodes,pelement);
   elem( 4,10,11,nodes,pelement);
   elem( 6, 7,13,nodes,pelement);
   elem( 6,12,13,nodes,pelement);
   elem( 7, 8,13,nodes,pelement);
   elem( 8,13,14,nodes,pelement);
   elem( 8,14,15,nodes,pelement);
   elem( 8, 9,15,nodes,pelement);
   elem( 9,15,16,nodes,pelement);
   elem( 9,10,17,nodes,pelement);
   elem( 9,16,17,nodes,pelement);
   elem(10,11,18,nodes,pelement);
   elem(10,17,18,nodes,pelement);
   elem(12,13,22,nodes,pelement);
   elem(12,21,22,nodes,pelement);
   elem(13,14,19,nodes,pelement);
   elem(13,19,22,nodes,pelement);
   elem(19,22,23,nodes,pelement);
   elem(16,17,20,nodes,pelement);
   elem(17,20,26,nodes,pelement);
   elem(20,25,26,nodes,pelement);
   elem(17,18,27,nodes,pelement);
   elem(17,26,27,nodes,pelement);
   elem(21,22,29,nodes,pelement);
   elem(21,28,29,nodes,pelement);
   elem(22,23,30,nodes,pelement);
   elem(22,29,30,nodes,pelement);
   elem(23,24,30,nodes,pelement);
   elem(24,30,31,nodes,pelement);
   elem(24,25,31,nodes,pelement);
   elem(25,26,31,nodes,pelement);
   elem(26,31,32,nodes,pelement);
   elem(26,27,33,nodes,pelement);
   elem(26,32,33,nodes,pelement);
   elem(28,29,35,nodes,pelement);
   elem(28,34,35,nodes,pelement);
   elem(29,30,35,nodes,pelement);
   elem(30,35,36,nodes,pelement);
   elem(30,31,36,nodes,pelement);
   elem(31,36,37,nodes,pelement);
   elem(31,32,37,nodes,pelement);
   elem(32,37,38,nodes,pelement);
   elem(32,33,38,nodes,pelement);
   elem(33,38,39,nodes,pelement);
   elem(14,19,40,nodes,pelement);
   elem(14,15,40,nodes,pelement);
   elem(15,16,40,nodes,pelement);
   elem(16,20,40,nodes,pelement);
   elem(20,25,40,nodes,pelement);
   elem(24,25,40,nodes,pelement);
   elem(23,24,40,nodes,pelement);
   elem(19,23,40,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
   for (pel = FIRSTELEMENT(pgrid), i = 0; pel; pel = pel->succ, i++)
      if (i < 52)
         ITYPE(pel) = 8;
      else
         ITYPE(pel) = 4;
   for (pn = FIRSTNODE(pgrid); pn; pn = pn->succ)
      ITYPE(pn) = 0;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         ITYPE(pel->n[i]) |= ITYPE(pel); 
}

#else

void cube_with_magnet(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: cube_with_magnet not available.\n"); }

#endif

void obdelnik3(pnode,pvert,pelement,pgrid)  
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT s,t=1.8;
   INT pocet=41;

   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = pocet; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + pocet;
  
   s=0.025*sqrt(3.);

   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.2;
   vertexes[1].x[1] = 0.; 
   vertexes[2].x[0] = 0.4;
   vertexes[2].x[1] = 0.;
   vertexes[3].x[0] = 0.825;
   vertexes[3].x[1] = 0.;
   vertexes[4].x[0] = 1.1;
   vertexes[4].x[1] = 0.;
   vertexes[5].x[0] = 1.375;
   vertexes[5].x[1] = 0.;
   vertexes[6].x[0] = 1.65;
   vertexes[6].x[1] = 0.;
   vertexes[7].x[0] = 1.925;
   vertexes[7].x[1] = 0.;
   vertexes[8].x[0] = 2.2;
   vertexes[8].x[1] = 0.;
  
   vertexes[9].x[0] = 0.;
   vertexes[9].x[1] = 0.205;
//   vertexes[10].x[0] = 0.275;
//   vertexes[10].x[1] = 0.205; 
   vertexes[11].x[0] = 0.4;
   vertexes[11].x[1] = 0.205;
   vertexes[12].x[0] = 0.825;
   vertexes[12].x[1] = 0.205;
   vertexes[13].x[0] = 1.1;
   vertexes[13].x[1] = 0.205;
   vertexes[14].x[0] = 1.375;
   vertexes[14].x[1] = 0.205;
   vertexes[15].x[0] = 1.65;
   vertexes[15].x[1] = 0.205;
   vertexes[16].x[0] = 1.925;
   vertexes[16].x[1] = 0.205;
   vertexes[17].x[0] = 2.2;
   vertexes[17].x[1] = 0.205;

   vertexes[18].x[0] = 0.;
   vertexes[18].x[1] = 0.41;
   vertexes[19].x[0] = 0.2;
   vertexes[19].x[1] = 0.41; 
   vertexes[20].x[0] = 0.4;
   vertexes[20].x[1] = 0.41;
   vertexes[21].x[0] = 0.825;
   vertexes[21].x[1] = 0.41;
   vertexes[22].x[0] = 1.1;
   vertexes[22].x[1] = 0.41;
   vertexes[23].x[0] = 1.375;
   vertexes[23].x[1] = 0.41;
   vertexes[24].x[0] = 1.65;
   vertexes[24].x[1] = 0.41;
   vertexes[25].x[0] = 1.925;
   vertexes[25].x[1] = 0.41;
   vertexes[26].x[0] = 2.2;
   vertexes[26].x[1] = 0.41;

   vertexes[10].x[0] = 0.175;
   vertexes[10].x[1] = 0.2+s; 
   vertexes[27].x[0] = 0.15;
   vertexes[27].x[1] = 0.2; 
   vertexes[28].x[0] = 0.25;
   vertexes[28].x[1] = 0.2;
   vertexes[29].x[0] = 0.175;
   vertexes[29].x[1] = 0.2-s;
   vertexes[30].x[0] = 0.225;
   vertexes[30].x[1] = 0.2-s;
   vertexes[31].x[0] = 0.225;
   vertexes[31].x[1] = 0.2+s;
   vertexes[32].x[0] = 0.1;
   vertexes[32].x[1] = 0.25;
   vertexes[33].x[0] = 0.2;
   vertexes[33].x[1] = 0.32;
   vertexes[34].x[0] = 0.3;
   vertexes[34].x[1] = 0.25;
   vertexes[35].x[0] = 0.3;
   vertexes[35].x[1] = 0.15;
   vertexes[36].x[0] = 0.2;
   vertexes[36].x[1] = 0.08;
   vertexes[37].x[0] = 0.1;
   vertexes[37].x[1] = 0.15;
   vertexes[38].x[0] = 0.612;
   vertexes[38].x[1] = 0.205;
   vertexes[39].x[0] = 0.612;
   vertexes[39].x[1] = 0.;
   vertexes[40].x[0] = 0.612;
   vertexes[40].x[1] = 0.41;


   for (i = 0; i < 10; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 17; i < 27; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 11; i < 17; i++)
      vertexes[i].type = 0;
   for (i=32;i<39;i++)
	   vertexes[i].type = 0;

   vertexes[10].type = 1|NMASK;
   vertexes[27].type = 1|NMASK;
   vertexes[28].type = 1|NMASK;
   vertexes[29].type = 1|NMASK;
   vertexes[30].type = 1|NMASK;
   vertexes[31].type = 1|NMASK;
   vertexes[39].type = 1|NMASK;
   vertexes[40].type = 1|NMASK;
   vertexes[17].type = 1|NMASK;  //=1 pro do nothing, je-li dirichlet tak =3
   
/* soupani */
   /*
   for(i=0;i<pocet;i++){
	   if(vertexes[i].x[0]<0.4 && vertexes[i].x[0]>0. && vertexes[i].x[1]>0. && vertexes[i].x[1]<0.41)
         vertexes[i].x[1]+=0.1*sin(t);
       if (vertexes[i].x[1]>0.41) vertexes[i].x[1]-=0.05*sin(t);
       if (vertexes[i].x[1]<0.) vertexes[i].x[1]+=0.05*sin(t);
   }
   /**/
/* dosoupano */      
   for (i = 0; i < pocet; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    


   elem( 12, 13, 22,nodes,pelement);
   elem( 13, 14, 23,nodes,pelement);
   elem( 14, 15, 24,nodes,pelement);
   elem( 15, 16, 25,nodes,pelement);
   elem( 16, 17, 26,nodes,pelement);

   elem( 12, 22, 21,nodes,pelement);
   elem( 13, 23, 22,nodes,pelement);
   elem( 14, 24, 23,nodes,pelement);
   elem( 15, 25, 24,nodes,pelement);
   elem( 16, 26, 25,nodes,pelement);

   elem( 12, 3, 4,nodes,pelement);
   elem( 13, 12, 4,nodes,pelement);
   elem( 13, 4, 5,nodes,pelement);
   elem( 14, 13, 5,nodes,pelement);
   elem( 14, 5, 6,nodes,pelement);
   elem( 15, 14, 6,nodes,pelement);
   elem( 15, 6, 7,nodes,pelement);
   elem( 16, 15, 7,nodes,pelement);
   elem( 16, 7, 8,nodes,pelement);
   elem( 17, 16, 8,nodes,pelement);
   elem( 32, 27, 10,nodes,pelement);
   elem( 18, 9, 32,nodes,pelement);
   elem( 33, 19, 18,nodes,pelement);
   elem( 33, 18, 32,nodes,pelement);
   elem( 33, 32, 10,nodes,pelement);
   elem( 33, 10, 31,nodes,pelement);
   elem( 33, 20, 19,nodes,pelement);
   elem( 34, 20, 33,nodes,pelement);
   elem( 34, 33, 31,nodes,pelement);
   elem( 34, 31, 28,nodes,pelement);
   elem( 34, 11, 20,nodes,pelement);
   elem( 35, 28, 30,nodes,pelement);
   elem( 35, 2, 11,nodes,pelement);
   elem( 35, 34, 28,nodes,pelement);
   elem( 35, 11, 34,nodes,pelement);
   elem( 36, 35, 30,nodes,pelement);
   elem( 36, 30, 29,nodes,pelement);
   elem( 36, 0, 1,nodes,pelement);
   elem( 36, 1, 2,nodes,pelement);
   elem( 36, 2, 35,nodes,pelement);
   elem( 37, 9, 0,nodes,pelement);
   elem( 37, 0, 36,nodes,pelement);
   elem( 37, 36, 29,nodes,pelement);
   elem( 37, 29, 27,nodes,pelement);
   elem( 37, 32, 9,nodes,pelement);
   elem( 37, 27, 32,nodes,pelement);
   elem( 38, 12, 21,nodes,pelement);
   elem( 38, 3, 12,nodes,pelement);
   elem( 39, 3, 38,nodes,pelement);
   elem( 39, 38, 11,nodes,pelement);
   elem( 39, 11, 2,nodes,pelement);
   elem( 40, 11, 38,nodes,pelement);
   elem( 40, 20, 11,nodes,pelement);
   elem( 40, 38, 21,nodes,pelement);
//   elem( , , ,nodes,pelement);
   
   (*pelement - 1)->succ = NULL;      
}

double hex(x,y,r)
FLOAT x,y,r;
{
FLOAT eps=1.e-4;

	if ((fabs(y-r*sqrt(3.)/2)<eps)&&(x<r/2.)) return(1.);
    else if ((x<r+eps)&&(x>=r/2.-eps)&&(fabs(sqrt(3.)*x+y-sqrt(3.)*r)<eps)) return(1.);
	     else return(0.);
}

void new_shift_nodes_of_hexagon2(mg)
MULTIGRID *mg;
{
   GRID *tGrid;
   NODE *pnode;
   FLOAT q, a0, a1, a, s, r, RR, q2,aa0, aa1;
   

   q=2.;r=0.05;RR=0.15;
   for (tGrid = FIRSTGRID(mg); tGrid; tGrid = tGrid->finer)
      for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode)){
            a0 = pnode->myvertex->x[0]-0.2;
            a1 = pnode->myvertex->x[1]-0.2;
            a  = sqrt(a0*a0 + a1*a1);
            if (RR>a){                  /*bod je blizko stredu kruznice => budeme ho soupat*/   
				aa0=fabs(a0);aa1=fabs(a1);
				if (hex(aa0,aa1,r)){    /* bod je na hexagonu, tak ho posuneme na kruznici*/
					q2=r/a;             /* o polomeru r*/
					//printf("hex %f %f",pnode->myvertex->x[0],pnode->myvertex->x[1]);
			    pnode->myvertex->x[0] =0.2+q2*a0;
                pnode->myvertex->x[1] =0.2+q2*a1;
				//printf("-> %f %f\n",pnode->myvertex->x[0],pnode->myvertex->x[1]);
				}
				else{          /* bod neni na hexagonu, ale musime ho posunout */
				//printf("    %f %f",pnode->myvertex->x[0],pnode->myvertex->x[1]);
		      	pnode->myvertex->x[0] =0.2+(1.+q*(RR-a))*a0;
                pnode->myvertex->x[1] =0.2+(1.+q*(RR-a))*a1;
				//printf("-> %f %f\n",pnode->myvertex->x[0],pnode->myvertex->x[1]);
				}
			}
			else {
			//printf("                     CHUDACI                         \n");
			//printf(" %f %f\n",pnode->myvertex->x[0],pnode->myvertex->x[1]);
			}
         }
}

void benchmark_channel(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT r, s, t, u, v, q=0.1, c0=0.2, c1=0.2, z=sqrt(3.);
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 49; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 49;
  
   r = 0.5*sqrt(3.)*q;
   s = 0.5*r;
   t = 0.5*q;
   u = 0.25*q;
   v = 0.75*q;
   vertexes[0].x[0] = c0 + t;
   vertexes[0].x[1] = c1;
   vertexes[1].x[0] = c0 + u;
   vertexes[1].x[1] = c1 + s;
   vertexes[2].x[0] = c0 - u;
   vertexes[2].x[1] = c1 + s;
   vertexes[3].x[0] = c0 - t;
   vertexes[3].x[1] = c1;
   vertexes[4].x[0] = c0 - u;
   vertexes[4].x[1] = c1 - s;
   vertexes[5].x[0] = c0 + u;
   vertexes[5].x[1] = c1 - s;
   vertexes[6].x[0] = c0 + q;
   vertexes[6].x[1] = c1;
   vertexes[7].x[0] = c0 + v;
   vertexes[7].x[1] = c1 + s;
   vertexes[8].x[0] = c0 + t;
   vertexes[8].x[1] = c1 + r;
   vertexes[9].x[0] = c0;
   vertexes[9].x[1] = c1 + r;
   vertexes[10].x[0] = c0 - t;
   vertexes[10].x[1] = c1 + r;
   vertexes[11].x[0] = c0 - v;
   vertexes[11].x[1] = c1 + s;
   vertexes[12].x[0] = c0 - q;
   vertexes[12].x[1] = c1;
   vertexes[13].x[0] = c0 - v;
   vertexes[13].x[1] = c1 - s;
   vertexes[14].x[0] = c0 - t;
   vertexes[14].x[1] = c1 - r;
   vertexes[15].x[0] = c0;
   vertexes[15].x[1] = c1 - r;
   vertexes[16].x[0] = c0 + t;
   vertexes[16].x[1] = c1 - r;
   vertexes[17].x[0] = c0 + v;
   vertexes[17].x[1] = c1 - s;

   vertexes[18].x[0] = 0.2+0.42/z;
   vertexes[18].x[1] = 0.2;
   vertexes[19].x[0] = 0.2+0.315/z;
   vertexes[19].x[1] = 0.305;
   vertexes[20].x[0] = 0.2+0.21/z;
   vertexes[20].x[1] = 0.41;
   vertexes[21].x[0] = 0.2;
   vertexes[21].x[1] = 0.41;
   vertexes[22].x[0] = 0.2-0.21/z;
   vertexes[22].x[1] = 0.41;
   vertexes[23].x[0] = 0.1-0.105/z;
   vertexes[23].x[1] = 0.305;
   vertexes[24].x[0] = 0.;
   vertexes[24].x[1] = 0.2;
   vertexes[25].x[0] = 0.1-0.1/z;
   vertexes[25].x[1] = 0.1;
   vertexes[26].x[0] = 0.2-0.2/z;
   vertexes[26].x[1] = 0.;
   vertexes[27].x[0] = 0.2;
   vertexes[27].x[1] = 0.;
   vertexes[28].x[0] = 0.2+0.2/z;
   vertexes[28].x[1] = 0.;
   vertexes[29].x[0] = 0.2+0.31/z;
   vertexes[29].x[1] = 0.1;

   vertexes[30].x[0] = 0.;
   vertexes[30].x[1] = 0.41;
   vertexes[31].x[0] = 0.;
   vertexes[31].x[1] = 0.305;
   vertexes[32].x[0] = 0.;
   vertexes[32].x[1] = 0.1;
   vertexes[33].x[0] = 0.;
   vertexes[33].x[1] = 0.;
   vertexes[34].x[0] = 0.765;
   vertexes[34].x[1] = 0.2;
   vertexes[35].x[0] = 1.175;
   vertexes[35].x[1] = 0.2;
   vertexes[36].x[0] = 1.585;
   vertexes[36].x[1] = 0.2;
   vertexes[37].x[0] = 1.995;
   vertexes[37].x[1] = 0.2;
   vertexes[38].x[0] = 0.56;
   vertexes[38].x[1] = 0.41;
   vertexes[39].x[0] = 0.97;
   vertexes[39].x[1] = 0.41;
   vertexes[40].x[0] = 1.38;
   vertexes[40].x[1] = 0.41;
   vertexes[41].x[0] = 1.79;
   vertexes[41].x[1] = 0.41;
   vertexes[42].x[0] = 2.2;
   vertexes[42].x[1] = 0.41;
   vertexes[43].x[0] = 0.56;
   vertexes[43].x[1] = 0.;
   vertexes[44].x[0] = 0.97;
   vertexes[44].x[1] = 0.;
   vertexes[45].x[0] = 1.38;
   vertexes[45].x[1] = 0.;
   vertexes[46].x[0] = 1.79;
   vertexes[46].x[1] = 0.;
   vertexes[47].x[0] = 2.2;
   vertexes[47].x[1] = 0.;
   vertexes[48].x[0] = 2.2;
   vertexes[48].x[1] = 0.2;

   for (i = 0; i < 49; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 6; i < 20; i++)
      vertexes[i].type = 0; 
   vertexes[23].type = 0; 
   vertexes[25].type = 0; 
   vertexes[29].type = 0; 
   for (i = 34; i < 38; i++)
      vertexes[i].type = 0; 

   for (i = 0; i < 49; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 6, 7,nodes,pelement);
   elem( 0, 1, 7,nodes,pelement);
   elem( 1, 7, 8,nodes,pelement);
   elem( 1, 8, 9,nodes,pelement);
   elem( 1, 2, 9,nodes,pelement);
   elem( 2, 9,10,nodes,pelement);
   elem( 2,10,11,nodes,pelement);
   elem( 2, 3,11,nodes,pelement);
   elem( 3,11,12,nodes,pelement);
   elem( 3,12,13,nodes,pelement);
   elem( 3, 4,13,nodes,pelement);
   elem( 4,13,14,nodes,pelement);
   elem( 4,14,15,nodes,pelement);
   elem( 4, 5,15,nodes,pelement);
   elem( 5,15,16,nodes,pelement);
   elem( 5,16,17,nodes,pelement);
   elem( 0, 5,17,nodes,pelement);
   elem( 0, 6,17,nodes,pelement);
   elem( 6,18,19,nodes,pelement);
   elem( 6, 7,19,nodes,pelement);
   elem( 7, 8,19,nodes,pelement);
   elem( 8,19,20,nodes,pelement);
   elem( 8,20,21,nodes,pelement);
   elem( 8, 9,21,nodes,pelement);
   elem( 9,10,21,nodes,pelement);
   elem(10,21,22,nodes,pelement);
   elem(10,22,23,nodes,pelement);
   elem(10,11,23,nodes,pelement);
   elem(11,12,23,nodes,pelement);
   elem(12,23,24,nodes,pelement);
   elem(12,24,25,nodes,pelement);
   elem(12,13,25,nodes,pelement);
   elem(13,14,25,nodes,pelement);
   elem(14,25,26,nodes,pelement);
   elem(14,26,27,nodes,pelement);
   elem(14,15,27,nodes,pelement);
   elem(15,16,27,nodes,pelement);
   elem(16,27,28,nodes,pelement);
   elem(16,28,29,nodes,pelement);
   elem(16,17,29,nodes,pelement);
   elem( 6,17,29,nodes,pelement);
   elem( 6,18,29,nodes,pelement);
   elem(22,23,30,nodes,pelement);
   elem(23,30,31,nodes,pelement);
   elem(23,24,31,nodes,pelement);
   elem(24,25,32,nodes,pelement);
   elem(25,32,33,nodes,pelement);
   elem(25,26,33,nodes,pelement);

   elem(19,20,38,nodes,pelement);
   elem(18,19,38,nodes,pelement);
   elem(18,34,38,nodes,pelement);
   elem(34,38,39,nodes,pelement);
   elem(34,35,39,nodes,pelement);
   elem(35,39,40,nodes,pelement);
   elem(35,36,40,nodes,pelement);
   elem(36,40,41,nodes,pelement);
   elem(36,37,41,nodes,pelement);
   elem(37,41,42,nodes,pelement);
   elem(37,42,48,nodes,pelement);
   elem(18,29,43,nodes,pelement);
   elem(28,29,43,nodes,pelement);
   elem(18,34,43,nodes,pelement);
   elem(34,43,44,nodes,pelement);
   elem(34,35,44,nodes,pelement);
   elem(35,44,45,nodes,pelement);
   elem(35,36,45,nodes,pelement);
   elem(36,45,46,nodes,pelement);
   elem(36,37,46,nodes,pelement);
   elem(37,46,47,nodes,pelement);
   elem(37,47,48,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

/* In the following, we consider a hexagon with vertices a1,...,a6
   connected to the point c inside the hexagon. On the lines connecting
   a_i with c we introduce points b_i in the distance r from c. Then we
   imagine that the points between the lines determined by a1,...,a6 and
   b1,...,b6 are shifted in such a way, that the line determined by b1,...,b6
   changes to a circle with radius r and the line determined by a1,...,a6
   does not change. The point x is shifted accordingly. */

INT is_in_triangle(x,a1,a2,a3)
FLOAT *x, *a1, *a2, *a3;
{
   FLOAT b[DIM2][DIM2];

   barycentric_coordinates(a1,a2,a3,b);
   if (DOT(b[0],x)+b[0][2] > -1.e-12 && DOT(b[1],x)+b[1][2] > -1.e-12 &&
       DOT(b[2],x)+b[2][2] > -1.e-12)
      return(1);
   else
      return(0);
}

void shift_point_in_triangle(x,a1,a2,c,r)
FLOAT *x, *a1, *a2, *c, r;
{
   FLOAT b1[DIM], b2[DIM], ax[DIM], bx[DIM], s, d0, d1, d2, e0, e1, e2;

   s = r/sqrt((a1[0]-c[0])*(a1[0]-c[0])+(a1[1]-c[1])*(a1[1]-c[1]));
   b1[0] = c[0] + s*(a1[0]-c[0]);
   b1[1] = c[1] + s*(a1[1]-c[1]);
   s = r/sqrt((a2[0]-c[0])*(a2[0]-c[0])+(a2[1]-c[1])*(a2[1]-c[1]));
   b2[0] = c[0] + s*(a2[0]-c[0]);
   b2[1] = c[1] + s*(a2[1]-c[1]);
   d0 = c[1] - x[1];
   d1 = x[0] - c[0];
   d2 = d0*x[0] + d1*x[1];
   e0 = a1[1] - a2[1];
   e1 = a2[0] - a1[0];
   e2 = e0*a2[0] + e1*a2[1];
   s = d0*e1 - d1*e0;
   ax[0] = (e1*d2 - d1*e2)/s;
   ax[1] = (d0*e2 - e0*d2)/s;
   e0 = b1[1] - b2[1];
   e1 = b2[0] - b1[0];
   e2 = e0*b2[0] + e1*b2[1];
   s = d0*e1 - d1*e0;
   bx[0] = (e1*d2 - d1*e2)/s;
   bx[1] = (d0*e2 - e0*d2)/s;
   e2 = sqrt((bx[0]- c[0])*(bx[0]- c[0])+(bx[1]- c[1])*(bx[1]- c[1]));
   s = sqrt( (( x[0]-ax[0])*( x[0]-ax[0])+( x[1]-ax[1])*( x[1]-ax[1]))/
             ((bx[0]-ax[0])*(bx[0]-ax[0])+(bx[1]-ax[1])*(bx[1]-ax[1])) )*(r-e2);
   e0 = (bx[0] - c[0])/e2;
   e1 = (bx[1] - c[1])/e2;
   s += sqrt((x[0]-c[0])*(x[0]-c[0])+(x[1]-c[1])*(x[1]-c[1]));
   x[0] = c[0] + s*e0;
   x[1] = c[1] + s*e1;
}

void shift_point_in_hexagon(x,a1,a2,a3,a4,a5,a6,c,r)
FLOAT *x, *a1, *a2, *a3, *a4, *a5, *a6, *c, r;
{
   if (is_in_triangle(x,a1,a2,c))
      shift_point_in_triangle(x,a1,a2,c,r);
   else if (is_in_triangle(x,a2,a3,c))
      shift_point_in_triangle(x,a2,a3,c,r);
   else if (is_in_triangle(x,a3,a4,c))
      shift_point_in_triangle(x,a3,a4,c,r);
   else if (is_in_triangle(x,a4,a5,c))
      shift_point_in_triangle(x,a4,a5,c,r);
   else if (is_in_triangle(x,a5,a6,c))
      shift_point_in_triangle(x,a5,a6,c,r);
   else if (is_in_triangle(x,a6,a1,c))
      shift_point_in_triangle(x,a6,a1,c,r);
}

void shift_points_in_benchmark_channel(mg)
MULTIGRID *mg;
{
   GRID *tGrid;
   NODE *pnode;
   FLOAT a1[DIM], a2[DIM], a3[DIM], a4[DIM], a5[DIM], a6[DIM], c[DIM], 
         z=sqrt(3.);
 
   a1[0] = 0.2+0.42/z;
   a1[1] = 0.2;
   a2[0] = 0.2+0.21/z;
   a2[1] = 0.41;
   a3[0] = 0.2-0.21/z;
   a3[1] = 0.41;
   a4[0] = 0.;
   a4[1] = 0.2;
   a5[0] = 0.2-0.2/z;
   a5[1] = 0.;
   a6[0] = 0.2+0.2/z;
   a6[1] = 0.;
   c[0] = 0.2;
   c[1] = 0.2;

   for (tGrid = FIRSTGRID(mg); tGrid; tGrid = tGrid->finer)
      for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode))
            shift_point_in_hexagon(pnode->myvertex->x,a1,a2,a3,a4,a5,a6,c,0.05);
}

void square_layer(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT r, s, q=2.;
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 12; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 12;
  
   r = 0.5*sqrt(2.)*q;
   s = 0.5*r;
   vertexes[0].x[0] =  s;
   vertexes[0].x[1] = -s;
   vertexes[1].x[0] =  s;
   vertexes[1].x[1] =  s;
   vertexes[2].x[0] = -s;
   vertexes[2].x[1] =  s;
   vertexes[3].x[0] = -s;
   vertexes[3].x[1] = -s;
   vertexes[4].x[0] =  r;
   vertexes[4].x[1] = -r;
   vertexes[5].x[0] =  r;
   vertexes[5].x[1] =  0.;
   vertexes[6].x[0] =  r;
   vertexes[6].x[1] =  r;
   vertexes[7].x[0] =  0.;
   vertexes[7].x[1] =  r;
   vertexes[8].x[0] = -r;
   vertexes[8].x[1] =  r;
   vertexes[9].x[0] = -r;
   vertexes[9].x[1] =  0.;
   vertexes[10].x[0] = -r;
   vertexes[10].x[1] = -r;
   vertexes[11].x[0] =  0.;
   vertexes[11].x[1] = -r;
   for (i = 0; i < 12; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 0; i < 12; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 5,nodes,pelement);
   elem( 1, 5, 6,nodes,pelement);
   elem( 1, 6, 7,nodes,pelement);
   elem( 1, 2, 7,nodes,pelement);
   elem( 2, 7, 8,nodes,pelement);
   elem( 2, 8, 9,nodes,pelement);
   elem( 2, 3, 9,nodes,pelement);
   elem( 3, 9,10,nodes,pelement);
   elem( 3,10,11,nodes,pelement);
   elem( 0, 3,11,nodes,pelement);
   elem( 0, 4,11,nodes,pelement);
   elem( 0, 4, 5,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void shift_nodes_of_square(mg)
MULTIGRID *mg;
{
   GRID *tGrid;
   NODE *pnode;
   FLOAT q=sqrt(2.), a0, a1, a, s;

   for (tGrid = FIRSTGRID(mg); tGrid; tGrid = tGrid->finer)
      for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode)){
            a0 = fabs(pnode->myvertex->x[0]);
            a1 = fabs(pnode->myvertex->x[1]);
            a  = sqrt(a0*a0 + a1*a1);
            s  = MAX(a0,a1)*q/a;
            pnode->myvertex->x[0] *= s;
            pnode->myvertex->x[1] *= s;
         }
}

void pentagon_layer(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT r, s, t, u, q=2.;
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 15; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 15;
  
   r = cos(PI*0.1)*q;
   s = sin(PI*0.1)*q;
   t = cos(PI*0.2)*q;
   u = sin(PI*0.2)*q;
   vertexes[0].x[0] =  r*0.5;
   vertexes[0].x[1] =  s*0.5;
   vertexes[1].x[0] =  0.;
   vertexes[1].x[1] =  q*0.5;
   vertexes[2].x[0] = -r*0.5;
   vertexes[2].x[1] =  s*0.5;
   vertexes[3].x[0] = -u*0.5;
   vertexes[3].x[1] = -t*0.5;
   vertexes[4].x[0] =  u*0.5;
   vertexes[4].x[1] = -t*0.5;
   vertexes[5].x[0] =  r;
   vertexes[5].x[1] =  s;
   vertexes[7].x[0] =  0.;
   vertexes[7].x[1] =  q;
   vertexes[9].x[0] = -r;
   vertexes[9].x[1] =  s;
   vertexes[11].x[0] = -u;
   vertexes[11].x[1] = -t;
   vertexes[13].x[0] =  u;
   vertexes[13].x[1] = -t;
   vertexes[6].x[0] = (vertexes[5].x[0] + vertexes[7].x[0])*0.5;
   vertexes[6].x[1] = (vertexes[5].x[1] + vertexes[7].x[1])*0.5;
   vertexes[8].x[0] = (vertexes[7].x[0] + vertexes[9].x[0])*0.5;
   vertexes[8].x[1] = (vertexes[7].x[1] + vertexes[9].x[1])*0.5;
   vertexes[10].x[0] = (vertexes[9].x[0] + vertexes[11].x[0])*0.5;
   vertexes[10].x[1] = (vertexes[9].x[1] + vertexes[11].x[1])*0.5;
   vertexes[12].x[0] = (vertexes[11].x[0] + vertexes[13].x[0])*0.5;
   vertexes[12].x[1] = (vertexes[11].x[1] + vertexes[13].x[1])*0.5;
   vertexes[14].x[0] = (vertexes[5].x[0] + vertexes[13].x[0])*0.5;
   vertexes[14].x[1] = (vertexes[5].x[1] + vertexes[13].x[1])*0.5;
   for (i = 0; i < 15; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 0; i < 15; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 5, 6,nodes,pelement);
   elem( 0, 1, 6,nodes,pelement);
   elem( 1, 6, 7,nodes,pelement);
   elem( 1, 7, 8,nodes,pelement);
   elem( 1, 2, 8,nodes,pelement);
   elem( 2, 8, 9,nodes,pelement);
   elem( 2, 9,10,nodes,pelement);
   elem( 2, 3,10,nodes,pelement);
   elem( 3,10,11,nodes,pelement);
   elem( 3,11,12,nodes,pelement);
   elem( 3, 4,12,nodes,pelement);
   elem( 4,12,13,nodes,pelement);
   elem( 4,13,14,nodes,pelement);
   elem( 0, 4,14,nodes,pelement);
   elem( 0, 5,14,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void shift_nodes_of_pentagon(mg)
MULTIGRID *mg;
{
   GRID *tGrid;
   NODE *pnode;
   FLOAT x0, x1, d, p, q, r, s, t, u, w, z;
 
   r = cos(PI*0.1);
   s = sin(PI*0.1);
   t = cos(PI*0.2);
   u = sin(PI*0.2);
   w = s/r;
   z = t/u;
   q = s*u + r*t;

   for (tGrid = FIRSTGRID(mg); tGrid; tGrid = tGrid->finer)
      for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode)){
            x0 = pnode->myvertex->x[0];
            x1 = pnode->myvertex->x[1];
            if (fabs(x0) < 1.e-14 && fabs(x1) < 1.e-14)
               d = 1.;
            else if (x1 < -1.e-14 && fabs(p=x0/x1) < u/t)
               d = 1./t/sqrt(1.+p*p);
            else if (fabs(x0) < 1.e-14)
               d = 1.;
            else{
               p = x1/x0;
               if (x1 > 1.e-14 && fabs(p) > w)
                  d = (fabs(p) + (1.-s)/r)/sqrt(1.+p*p);
               else if (fabs(p) < z && x0 > 0.)
                  d = (s+t-p*(r-u))/q/sqrt(1.+p*p);
               else if (fabs(p) < z && x0 < 0.)
                  d = (s+t+p*(r-u))/q/sqrt(1.+p*p);
               else 
                  d = 1.;
            }
            pnode->myvertex->x[0] *= d;
            pnode->myvertex->x[1] *= d;
         }
}

void hexagon(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT r;
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 7; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 7;
  
   r = 0.5*sqrt(3.);
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 1.;
   vertexes[1].x[1] = 0.;
   vertexes[2].x[0] = 0.5;
   vertexes[2].x[1] = r;
   vertexes[3].x[0] = -0.5;
   vertexes[3].x[1] = r;
   vertexes[4].x[0] = -1.;
   vertexes[4].x[1] = 0;
   vertexes[5].x[0] = -0.5;
   vertexes[5].x[1] = -r;
   vertexes[6].x[0] = 0.5;
   vertexes[6].x[1] = -r;
   for (i = 1; i < 7; i++)
      vertexes[i].type = 1|NMASK;
   vertexes[0].type = 0; 
   for (i = 0; i < 7; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 2,nodes,pelement);
   elem( 0, 2, 3,nodes,pelement);
   elem( 0, 3, 4,nodes,pelement);
   elem( 0, 4, 5,nodes,pelement);
   elem( 0, 5, 6,nodes,pelement);
   elem( 0, 1, 6,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void hexagon2(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 7; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 7;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 1.;
   vertexes[1].x[1] = 0.;
   vertexes[2].x[0] = 1.;
   vertexes[2].x[1] = 1.;
   vertexes[3].x[0] = 0.;
   vertexes[3].x[1] = 1.;
   vertexes[4].x[0] = -1.;
   vertexes[4].x[1] = 0;
   vertexes[5].x[0] = -1.;
   vertexes[5].x[1] = -1.;
   vertexes[6].x[0] = 0.;
   vertexes[6].x[1] = -1.;
   for (i = 1; i < 7; i++)
      vertexes[i].type = 1|NMASK;
   vertexes[0].type = 0; 
   for (i = 0; i < 7; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 2,nodes,pelement);
   elem( 0, 2, 3,nodes,pelement);
   elem( 0, 3, 4,nodes,pelement);
   elem( 0, 4, 5,nodes,pelement);
   elem( 0, 5, 6,nodes,pelement);
   elem( 0, 1, 6,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void hexagon_layer(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   NODE *nodes;
   VERTEX *vertexes;
   FLOAT r, s, t, u, v, q=2.;
   INT i;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 18; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 18;
  
   r = 0.5*sqrt(3.)*q;
   s = 0.5*r;
   t = 0.5*q;
   u = 0.25*q;
   v = 0.75*q;
   vertexes[0].x[0] = t;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = u;
   vertexes[1].x[1] = s;
   vertexes[2].x[0] = -u;
   vertexes[2].x[1] = s;
   vertexes[3].x[0] = -t;
   vertexes[3].x[1] = 0.;
   vertexes[4].x[0] = -u;
   vertexes[4].x[1] = -s;
   vertexes[5].x[0] = u;
   vertexes[5].x[1] = -s;
   vertexes[6].x[0] = q;
   vertexes[6].x[1] = 0.;
   vertexes[7].x[0] = v;
   vertexes[7].x[1] = s;
   vertexes[8].x[0] = t;
   vertexes[8].x[1] = r;
   vertexes[9].x[0] = 0.;
   vertexes[9].x[1] = r;
   vertexes[10].x[0] = -t;
   vertexes[10].x[1] = r;
   vertexes[11].x[0] = -v;
   vertexes[11].x[1] = s;
   vertexes[12].x[0] = -q;
   vertexes[12].x[1] = 0.;
   vertexes[13].x[0] = -v;
   vertexes[13].x[1] = -s;
   vertexes[14].x[0] = -t;
   vertexes[14].x[1] = -r;
   vertexes[15].x[0] = 0.;
   vertexes[15].x[1] = -r;
   vertexes[16].x[0] = t;
   vertexes[16].x[1] = -r;
   vertexes[17].x[0] = v;
   vertexes[17].x[1] = -s;
   for (i = 0; i < 18; i++)
      vertexes[i].type = 1|NMASK;
   for (i = 0; i < 18; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 6, 7,nodes,pelement);
   elem( 0, 1, 7,nodes,pelement);
   elem( 1, 7, 8,nodes,pelement);
   elem( 1, 8, 9,nodes,pelement);
   elem( 1, 2, 9,nodes,pelement);
   elem( 2, 9,10,nodes,pelement);
   elem( 2,10,11,nodes,pelement);
   elem( 2, 3,11,nodes,pelement);
   elem( 3,11,12,nodes,pelement);
   elem( 3,12,13,nodes,pelement);
   elem( 3, 4,13,nodes,pelement);
   elem( 4,13,14,nodes,pelement);
   elem( 4,14,15,nodes,pelement);
   elem( 4, 5,15,nodes,pelement);
   elem( 5,15,16,nodes,pelement);
   elem( 5,16,17,nodes,pelement);
   elem( 0, 5,17,nodes,pelement);
   elem( 0, 6,17,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void shift_nodes_of_hexagon(mg)
MULTIGRID *mg;
{
   GRID *tGrid;
   NODE *pnode;
   FLOAT q=sqrt(3.), a0, a1, a, s;

   for (tGrid = FIRSTGRID(mg); tGrid; tGrid = tGrid->finer)
      for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
         if (IS_TOP_NODE(pnode)){
            a0 = fabs(pnode->myvertex->x[0]);
            a1 = fabs(pnode->myvertex->x[1]);
            a  = sqrt(a0*a0 + a1*a1);
            if (a > 1.e-14){
               s = q*a;
               if ( (a0 < 1.e-14) || (a1/a0 > q) )
                  s = 2.*a1/s;
               else
                  s = (q*a0 + a1)/s; 
               pnode->myvertex->x[0] *= s;
               pnode->myvertex->x[1] *= s;
            }
         }
}

void rtriangle(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = 10; 
   nodes = *pnode;
   vertexes = *pvert;
   *pvert = *pvert + 10;
  
   vertexes[0].x[0] = 0.;
   vertexes[0].x[1] = 0.;
   vertexes[1].x[0] = 0.;
   vertexes[1].x[1] = 1.; 
   vertexes[2].x[0] = 0.;
   vertexes[2].x[1] = 2.;
   vertexes[3].x[0] = 0.;
   vertexes[3].x[1] = 3.;
   vertexes[4].x[0] = 1.;
   vertexes[4].x[1] = 0.;
   vertexes[5].x[0] = 1.;
   vertexes[5].x[1] = 1.;
   vertexes[6].x[0] = 1.;
   vertexes[6].x[1] = 2.;
   vertexes[7].x[0] = 2.;
   vertexes[7].x[1] = 0.;
   vertexes[8].x[0] = 2.;
   vertexes[8].x[1] = 1.;
   vertexes[9].x[0] = 3.;
   vertexes[9].x[1] = 0.;
   for (i = 0; i < 10; i++)
      vertexes[i].type = 1|NMASK; 
   vertexes[5].type = 0; 
   for (i = 0; i < 10; i++){
      nmake(*pnode,*pnode+1,i+1,(NODE*)NULL,vertexes+i);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
    
   elem( 0, 1, 4,nodes,pelement);
   elem( 1, 4, 5,nodes,pelement);
   elem( 1, 2, 5,nodes,pelement);
   elem( 2, 5, 6,nodes,pelement);
   elem( 2, 3, 6,nodes,pelement);
   elem( 4, 5, 7,nodes,pelement);
   elem( 5, 7, 8,nodes,pelement);
   elem( 5, 6, 8,nodes,pelement);
   elem( 7, 8, 9,nodes,pelement);
   (*pelement - 1)->succ = NULL;      
}

void smooth_triangulation(tGrid,n)
GRID *tGrid;
INT n;
{
   NODE *pnode, *pn;
   LINK *pl;
   INT i, m;
   FLOAT x0, x1;
   	
   for (i = 0; i < n; i++)
      for (pnode = FIRSTNODE(tGrid); pnode; pnode=SUCC(pnode))
         if (pnode->myvertex->type == 0){
            m = 0;
            x0 = 0.;
	    x1 = 0.;
            for (pl = TSTART(pnode); pl; pl = NEXT(pl)){
	       pn=NBNODE(pl);
	       x0 += pn->myvertex->x[0];
	       x1 += pn->myvertex->x[1];
               m++;
            }
	    pnode->myvertex->x[0] = x0/m;
	    pnode->myvertex->x[1] = x1/m;
	 }
}

#if N_DATA & NODE_ITYPE

void smooth_triangulation_for_itype(tGrid,n,itype)
GRID *tGrid;
INT n, itype;
{
   NODE *pnode, *pn;
   LINK *pl;
   INT i, m;
   FLOAT x0, x1;
   	
   for (i = 0; i < n; i++)
      for (pnode = FIRSTNODE(tGrid); pnode; pnode=SUCC(pnode))
         if (NOT_BN(pnode) && ITYPE(pnode) == itype){
            m = 0;
            x0 = 0.;
	    x1 = 0.;
            for (pl = TSTART(pnode); pl; pl = NEXT(pl)){
	       pn=NBNODE(pl);
	       x0 += pn->myvertex->x[0];
	       x1 += pn->myvertex->x[1];
               m++;
            }
	    pnode->myvertex->x[0] = x0/m;
	    pnode->myvertex->x[1] = x1/m;
	 }
}

#else

void smooth_triangulation_for_itype(tGrid,n,itype)
GRID *tGrid; INT n, itype;
{  eprintf("Error: smooth_triangulation_for_itype not available.\n"); }

#endif

#if F_DATA & CURVED_FACE_MIDDLE

void set_c_midpoint_on_all_faces(tGrid,pbpoint)
GRID *tGrid;
BPOINT **pbpoint;
{
   ELEMENT *pel;
   FACE *f0, *f1, *f2;
   FLOAT *x0, *x1, *x2;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      FACES_OF_ELEMENT(f0,f1,f2,pel);
      AVERAGE(x0,x1,(*pbpoint)->x);
      (*pbpoint)->type = 0; 
      f2->c_midpoint = (*pbpoint)++;
      AVERAGE(x0,x2,(*pbpoint)->x);
      (*pbpoint)->type = 0; 
      f1->c_midpoint = (*pbpoint)++;
      AVERAGE(x1,x2,(*pbpoint)->x);
      (*pbpoint)->type = 0; 
      f0->c_midpoint = (*pbpoint)++;
   }
}

void set_c_midpoint_on_boundary_faces(tGrid,pbpoint)
GRID *tGrid;
BPOINT **pbpoint;
{
   ELEMENT *pel;
   FACE *f0, *f1, *f2;
   FLOAT *x0, *x1, *x2;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      FACES_OF_ELEMENT(f0,f1,f2,pel);
      if (IS_BF(f2)){
         AVERAGE(x0,x1,(*pbpoint)->x);
         (*pbpoint)->type = 0; 
         f2->c_midpoint = (*pbpoint)++;
      }
      if (IS_BF(f1)){
         AVERAGE(x0,x2,(*pbpoint)->x);
         (*pbpoint)->type = 0; 
         f1->c_midpoint = (*pbpoint)++;
      }
      if (IS_BF(f0)){
         AVERAGE(x1,x2,(*pbpoint)->x);
         (*pbpoint)->type = 0; 
         f0->c_midpoint = (*pbpoint)++;
      }
   }
}

void shift_boundary_midpoints_on_unit_circle(mg)
MULTIGRID *mg;
{
   GRID *theGrid;
   FACE *pface;
   FLOAT d;

   for (theGrid = TOP_GRID(mg); theGrid; theGrid = theGrid->coarser)
      for (pface=FIRSTF(theGrid); pface; pface=pface->succ)
         if (IS_BF(pface) && pface->c_midpoint){
            d = sqrt(DOT(pface->c_midpoint->x,pface->c_midpoint->x));
            SET13(pface->c_midpoint->x,d)
         }
}

void shift_benchmark_boundary_midpoints_on_circle(mg)
MULTIGRID *mg;
{
   GRID *theGrid;
   FACE *pface;
   FLOAT xc[DIM], d, r=0.05;

   xc[0] = xc[1] = 0.2;
   for (theGrid = TOP_GRID(mg); theGrid; theGrid = theGrid->coarser)
      for (pface=FIRSTF(theGrid); pface; pface=pface->succ)
         if (IS_BF(pface) && pface->c_midpoint && 
             ((d=DISTANCE(pface->c_midpoint->x,xc)) < 0.1)){
            d = r/d;
            pface->c_midpoint->x[0] = xc[0] + d*(pface->c_midpoint->x[0]-xc[0]);
            pface->c_midpoint->x[1] = xc[1] + d*(pface->c_midpoint->x[1]-xc[1]);
         }
}

void check_c_midpoints(mg)
MULTIGRID *mg;
{
   GRID *theGrid;
   ELEMENT *pel;

   for (theGrid = TOP_GRID(mg); theGrid; theGrid = theGrid->coarser)
      for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
         if ( (pel->f[0]->c_midpoint && pel->f[1]->c_midpoint) ||
              (pel->f[0]->c_midpoint && pel->f[2]->c_midpoint) ||
              (pel->f[1]->c_midpoint && pel->f[2]->c_midpoint) )
            eprintf("Error: element with more than one curved edge.\n");
}

#else

void set_c_midpoint_on_all_faces(tGrid,pbpoint)
GRID *tGrid; BPOINT **pbpoint;
{  eprintf("Error: set_c_midpoint_on_all_faces not available.\n");  }

void set_c_midpoint_on_boundary_faces(tGrid,pbpoint)
GRID *tGrid; BPOINT **pbpoint;
{  eprintf("Error: set_c_midpoint_on_boundary_faces not available.\n");  }

void shift_boundary_midpoints_on_unit_circle(mg)
MULTIGRID *mg;
{  eprintf("Error: shift_boundary_midpoints_on_unit_circle not available.\n"); }
 
void shift_benchmark_boundary_midpoints_on_circle(mg)
MULTIGRID *mg;
{  eprintf("Error: shift_benchmark_boundary_midpoints_on_circle not available.\n"); }

void check_c_midpoints(mg)
MULTIGRID *mg;
{}

#endif

void read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, type, nvert; 
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
      (*pvert)->type = (UNS)type;
      if (SC_EXAMPLE != 9 && SC_EXAMPLE != 91)
         if (fabs((*pvert)->x[0]-0.0) < EPSA) (*pvert)->type |= NMASK;
      if (fabs((*pvert)->x[0]-1.0) < EPSA) (*pvert)->type |= NMASK;
      if (fabs((*pvert)->x[1]-0.0) < EPSA) (*pvert)->type |= NMASK;
      if (fabs((*pvert)->x[1]-1.0) < EPSA) (*pvert)->type |= NMASK;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   while (i1 > -1){
      elem(i1,i2,i3,nodes,pelement);
      fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   }
   (*pelement - 1)->succ = NULL;
}

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void read_vertexes_and_elements_with_itype(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, type, it, nvert; 
   DOUBLE x0, x1;
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type,&it);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
//    if (type & 1) it = 48;
      ITYPE(*pnode) = (UNS)it;
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   while (i1 > -1){
      ITYPE(*pelement) = it;
      elem(i1,i2,i3,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   }
   (*pelement - 1)->succ = NULL;
}

void read_vertexes_and_elements_from_two_files_with_itype(pnode,pvert,pelement,
                                                                          pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, type, it, nvert=0, n=0; 
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   nodes = *pnode;
   fp = fopen("data0","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type,&it);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      ITYPE(*pnode) = (UNS)it;
      (*pnode)++;
   } 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   while (i1 > -1){
      ITYPE(*pelement) = it;
      elem(i1,i2,i3,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&it);
   }
   fclose(fp);
   n = nvert;
   fp = fopen("data10","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
      if (type & 1)
         type = 1|NMASK;
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,n+i,(NODE*)NULL,(*pvert)++);
      ITYPE(*pnode) = 32;
      (*pnode)++;
   }
   fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   while (i1 > -1){
      ITYPE(*pelement) = 32;
      elem(n+i1,n+i2,n+i3,nodes,pelement);
      fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   }
   fclose(fp);
   (*pnode-1)->succ = NULL;
   (*pelement - 1)->succ = NULL;
   pgrid->maxNodeIndex = n+nvert;
}

#else

void read_vertexes_and_elements_with_itype(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: read_vertexes_and_elements_with_itype not available.\n"); }

void read_vertexes_and_elements_from_two_files_with_itype(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: read_vertexes_and_elements_from_two_files_with_itype not available.\n"); }

#endif

void read_vertexes_and_elements_mso(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, type, nvert; 
   NODE *nodes;
   FILE *fp;

   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
	  (*pvert)->type = (UNS)type;
	  if ( FBCRITERION((*pvert)->x[0],(*pvert)->x[1],type) )
		  (*pvert)->type |= NMASK_FOR_FB;
	  if ( ZNCRITERION((*pvert)->x[0],(*pvert)->x[1],type) )
		  (*pvert)->type |= NMASK_ZN;
      if ( DIRCRITERION((*pvert)->x[0],(*pvert)->x[1],type) ) 
		  (*pvert)->type |= NMASK;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   }
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   while (i1 > -1){
      elem(i1,i2,i3,nodes,pelement);
      fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   }
   (*pelement - 1)->succ = NULL;
}

void Read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: Read_vertexes_and_elements not available.\n");  }

void read_vertexes_and_elements_from_two_files(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, type, nvert=0, n=0;
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   pgrid->minNodeIndex = 1;
   nodes = *pnode;
   fp = fopen("data1","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   while (i1 > -1){
      elem(i1,i2,i3,nodes,pelement);
      fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   }
   fclose(fp);
   n = nvert;
   fp = fopen("data2","r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,n+i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   while (i1 > -1){
      elem(n+i1,n+i2,n+i3,nodes,pelement);
      fscanf(fp,"%i %i %i",&i1,&i2,&i3);
   }
   fclose(fp);
   (*pnode-1)->succ = NULL;
   (*pelement - 1)->succ = NULL;
   pgrid->maxNodeIndex = n+nvert; 
}

void remove_elements_and_nodes(tGrid,d0,d1,eps)
GRID *tGrid;
FLOAT d0, d1, eps;
{
   ELEMENT *pel, *prev_el;
   NODE *pn, *prev_n;
   FLOAT *x0, *x1, *x2;

   d0 -= eps*0.5;
   d1 += eps*0.5;
   prev_el = FIRSTELEMENT(tGrid);
   for (pel = FIRSTELEMENT(tGrid)->succ; pel; pel = pel->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      if (x0[0] > d0 && x1[0] > d0 && x2[0] > d0 &&
          x0[1] < d1 && x1[1] < d1 && x2[1] < d1)
         prev_el->succ = pel->succ;
      else
         prev_el = pel;
   }
   VERTICES_OF_ELEMENT(x0,x1,x2,FIRSTELEMENT(tGrid));
   if (x0[0] > d0 && x1[0] > d0 && x2[0] > d0 &&
       x0[1] < d1 && x1[1] < d1 && x2[1] < d1)
      FIRSTELEMENT(tGrid) = FIRSTELEMENT(tGrid)->succ;
   d0 += eps;
   d1 -= eps;
   prev_n = FIRSTN(tGrid);
   for (pn = FIRSTN(tGrid)->succ; pn; pn = pn->succ)
      if (pn->myvertex->x[0] > d0 && pn->myvertex->x[1] < d1)
         prev_n->succ = pn->succ;
      else
         prev_n = pn;
   if (FIRSTN(tGrid)->myvertex->x[0] > d0 && FIRSTN(tGrid)->myvertex->x[1] < d1)
      FIRSTN(tGrid) = FIRSTN(tGrid)->succ;
}

#if N_DATA & NODE_ITYPE
#define SET_ITYPE_FOR_DUPLICATE_NODES  ITYPE(pn) |= ITYPE(pn1);
#else
#define SET_ITYPE_FOR_DUPLICATE_NODES
#endif

void remove_duplicate_nodes(tGrid,eps)
GRID *tGrid;
FLOAT eps;
{
   ELEMENT *pel;
   NODE *pn, *pn1, *prev_n;
   INT i=0;

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ){
      prev_n = pn;
      for (pn1 = pn->succ; pn1; pn1 = pn1->succ)
         if (fabs(pn->myvertex->x[0]-pn1->myvertex->x[0]) < eps &&
             fabs(pn->myvertex->x[1]-pn1->myvertex->x[1]) < eps){
            prev_n->succ = pn1->succ;
            NTYPE(pn) &= NTYPE(pn1);
            SET_ITYPE_FOR_DUPLICATE_NODES
            for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
               if (pel->n[0] == pn1) pel->n[0] = pn;
               if (pel->n[1] == pn1) pel->n[1] = pn;
               if (pel->n[2] == pn1) pel->n[2] = pn;
            }
            i++;
         }
         else
            prev_n = pn1;
   }
   printf("%i nodes removed.\n",i);
}

void save_elements_and_nodes(tGrid)
GRID *tGrid;
{
   ELEMENT *pel;
   NODE *pn;
   INT i = 0;
   FILE *fp;

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      pn->index2 = i++;
   fp = fopen("data","w");
   fprintf(fp,"%ld \n",i);
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      fprintf(fp,"%e %e %i \n",pn->myvertex->x[0],pn->myvertex->x[1],NTYPE(pn));
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fprintf(fp,"%ld %ld %ld\n",pel->n[0]->index2,pel->n[1]->index2,
                                                   pel->n[2]->index2);
   fprintf(fp,"%i %i %i ",-3,-3,-3);
   fclose(fp);
}

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void save_elements_and_nodes_with_el_itype(tGrid)
GRID *tGrid;
{
   ELEMENT *pel;
   NODE *pn;
   INT i = 0;
   FILE *fp;

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      pn->index2 = i++;
   fp = fopen("data_circular_magnet_with_sur","w");
   fprintf(fp,"%ld \n",i);
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      fprintf(fp,"%e %e %i %i\n",pn->myvertex->x[0]-SHIFT,
                                 pn->myvertex->x[1]-SHIFT,NTYPE(pn),ITYPE(pn));
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fprintf(fp,"%ld %ld %ld %ld\n",pel->n[0]->index2,pel->n[1]->index2,
                                     pel->n[2]->index2,ITYPE(pel));
   fprintf(fp,"%i %i %i %i ",-3,-3,-3,-3);
   fclose(fp);
}

DOUBLE angle(r0,r1)
DOUBLE r0, r1;
{
   DOUBLE phi;

   if (fabs(r0) < 1.e-10){
      if (r1 > 0.)
         phi = 0.5*M_PI;
      else
         phi = 1.5*M_PI;
   }
   else{
      phi = atan(r1/r0);
      if (r0 < 0.)
         phi += M_PI;
      else if (r0 > 0. && r1 < 0.)
         phi += 2.*M_PI;
   }
   return(phi);
}

void save_elements_and_nodes_with_circ_surr(tGrid)
GRID *tGrid;
{
   ELEMENT *pel;
   VERTEX *bv[10000];
   NODE *pn;
   DOUBLE d, d0, d1, q, r=1.5, r2=2.25, pi2=2.*M_PI, phi, xr[200], yr[200],
          a[DIM], b[DIM], c[DIM], xm[DIM], xn[DIM], xz[DIM], ang[10000];
   INT i=0, j, k, n=0, m=50, ir[200][100];
   FILE *fp;

   for (pn = FIRSTN(tGrid); pn; pn = pn->succ){
      d0 = pn->myvertex->x[0] - 0.5;
      d1 = pn->myvertex->x[1] - 0.5;
      if (d0*d0 + d1*d1 > r2)
         pn->index2 = -1;
      else
         pn->index2 =  1;
      NTYPE(pn) = 0;
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (pel->n[0]->index2 == -1 || pel->n[1]->index2 == -1
                                  || pel->n[2]->index2 == -1)
         for (i = 0; i < 3; i++)
            if (pel->n[i]->index2 > 0)
               NTYPE(pel->n[i]) = 1|NMASK;
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      if (IS_BN(pn)){
         phi = angle(pn->myvertex->x[0]-0.5,pn->myvertex->x[1]-0.5);
         for (j = 0; j < n && phi > ang[j]; j++);
         for (k = n; k > j; k--){
            bv[k] = bv[k-1];
            ang[k] = ang[k-1];
         }
         bv[j] = pn->myvertex;
         ang[j] = phi;
         n++;
      }
   bv[n] = bv[0];
   ang[n] = pi2;
   q = pi2/n;
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      if (ITYPE(pn) == 32 && NOT_BN(pn)){
         b[0] = pn->myvertex->x[0]-0.5;
         b[1] = pn->myvertex->x[1]-0.5;
         phi = angle(b[0],b[1]);
         if (phi < 1.e-10 || fabs(phi-pi2) < 1.e-10)
            j = 1;
         else
            for (j = 0; j < n && phi > ang[j]; j++);
         SUBTR(bv[j-1]->x,bv[j]->x,a);
         c[0] = bv[j]->x[0]-0.5;
         c[1] = bv[j]->x[1]-0.5;
         d = (b[0]*c[1]-b[1]*c[0])/(a[0]*b[1]-a[1]*b[0]);
         xm[0] = d*bv[j-1]->x[0] + (1.-d)*bv[j]->x[0];
         xm[1] = d*bv[j-1]->x[1] + (1.-d)*bv[j]->x[1];
         xn[0] = 0.5 + r*(d*cos((j-1)*q) + (1.-d)*cos(j*q));
         xn[1] = 0.5 + r*(d*sin((j-1)*q) + (1.-d)*sin(j*q));
         if (pn->myvertex->x[0] > 1.299999 && fabs(b[1]/b[0]) < 1.250001){
            xz[0] = 1.3;
            xz[1] = 0.5 + 0.8*b[1]/b[0];
         }
         else if (pn->myvertex->x[0] < -0.299999 && fabs(b[1]/b[0]) < 1.250001){
            xz[0] = -0.3;
            xz[1] = 0.5 - 0.8*b[1]/b[0];
         }
         else if (pn->myvertex->x[1] > 1.499999 && fabs(b[0]/b[1]) < 0.800001){
            xz[0] = 0.5 + b[0]/b[1];
            xz[1] = 1.5;
         }
         else{
            xz[0] = 0.5 - b[0]/b[1];
            xz[1] = -0.5;
         }
         SUBTR(pn->myvertex->x,xm,a);
         SUBTR(xz,xm,b);
         if (fabs(b[0]) > fabs(b[1]))
            d = a[0]/b[0];
         else
            d = a[1]/b[1];
         pn->myvertex->x[0] = d*xz[0] + (1.-d)*xn[0];
         pn->myvertex->x[1] = d*xz[1] + (1.-d)*xn[1];
      }
   for (j = 0; j < n; j++){
      bv[j]->x[0] = 0.5 + r*cos(j*q);
      bv[j]->x[1] = 0.5 + r*sin(j*q);
   }
   i = 0;
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      if (pn->index2 > 0)
         pn->index2 = i++;
   if (m){
      for (j = 0; j < n; j++){
         xr[j] = cos(j*q);
         yr[j] = sin(j*q);
         ir[j][0] = bv[j]->topnode->index2;
         bv[j]->type = 0;
      }
      ir[n][0] = ir[0][0];
      for (k = 1; k <= m; k++){
         for (j = 0; j < n; j++)
            ir[j][k] = i++;
         ir[n][k] = ir[0][k];
      }
   }
   fp = fopen("data_circular_magnet_with_sur","w");
   fprintf(fp,"%ld \n",i);
   for (pn = FIRSTN(tGrid); pn; pn = pn->succ)
      if (pn->index2 >= 0)
         fprintf(fp,"%e %e %i %i\n",pn->myvertex->x[0]-SHIFT,
                                    pn->myvertex->x[1]-SHIFT,NTYPE(pn),
                                                             ITYPE(pn));
   if (m){
      d = r;
      q = 8./n;
      for (k = 1; k < m; k++){
         d *= 1. + q;
         for (j = 0; j < n; j++)
            fprintf(fp,"%e %e %i %i\n",0.5+d*xr[j]-SHIFT, 
                                       0.5+d*yr[j]-SHIFT, 0, 32);
      }
      d *= 1. + q;
      for (j = 0; j < n; j++)
         fprintf(fp,"%e %e %i %i\n",0.5+d*xr[j]-SHIFT, 
                                    0.5+d*yr[j]-SHIFT, 1|NMASK, 32);
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (pel->n[0]->index2 >= 0 && pel->n[1]->index2 >= 0
                                 && pel->n[2]->index2 >= 0)
      fprintf(fp,"%ld %ld %ld %ld\n",pel->n[0]->index2,pel->n[1]->index2,
                                     pel->n[2]->index2,ITYPE(pel));
   for (k = 1; k <= m; k++)
      for (j = 0; j < n; j++){
         fprintf(fp,"%ld %ld %ld %ld\n",ir[j][k-1],ir[j+1][k-1],ir[j][k],32);
         fprintf(fp,"%ld %ld %ld %ld\n",ir[j+1][k-1],ir[j][k],ir[j+1][k],32);
      }
   fprintf(fp,"%i %i %i %i ",-3,-3,-3,-3);
   fclose(fp);
}

#else

void save_elements_and_nodes_with_el_itype(tGrid)
GRID *tGrid;
{  eprintf("Error: save_elements_and_nodes_with_el_itype not available.\n");  }

void save_elements_and_nodes_with_circ_surr(tGrid)
GRID *tGrid;
{  eprintf("Error: save_elements_and_nodes_with_circ_surr not available.\n");  }

#endif

#elif (DIM == 2 && ELEMENT_TYPE == CUBE)

void cube(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na; 
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nv; i++) 
      for (j = 0; j < nv; j++) {
           (*pvert)->x[0] = i/(nv-1.0);
           (*pvert)->x[1] = j/(nv-1.0);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            elem(il, il+nv, il+nv+1, il+1, nodes, pelement); 
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_alt(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n, il, nv, na; 
   DOUBLE x, y, d;
   NODE *nodes;
   VERTEX *vertexes;
 
   nv = NV;               /* number of vertices in one direction */
   na = nv*nv;            /* number of all vertices */
 
   n = 0;
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   d=0.2/(nv-1.);
   for (i = 0; i < nv; i++){
      x = i/(nv-1.);
      if (i != 2*(i/2) && i < nv-1)
         x -= d;
      for (j = 0; j < nv; j++) {
           y = j/(nv-1.0);
           if (j != 2*(j/2) && j < nv-1)
              y -= d;
           (*pvert)->x[0] = x;
           (*pvert)->x[1] = y;
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   }
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nv-1; i++)
      for (j = 0; j < nv-1; j++) {
            il = j + nv*i;
            elem(il, il+nv, il+nv+1, il+1, nodes, pelement); 
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_xy(pnode,pvert,pelement,pgrid)  /*  cube = square  */
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, j, n=0, il, nvx=NVX, nvy=NVY, na=NVX*NVY;
   NODE *nodes;
   VERTEX *vertexes;
 
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   vertexes = *pvert;
  
   for (i = 0; i < nvx; i++) 
      for (j = 0; j < nvy; j++) {
           (*pvert)->x[0] = i/(nvx-1.);
           (*pvert)->x[1] = j/(nvy-1.);
           (*pvert)->type = 0; 
           nmake(*pnode,*pnode+1,++n,(NODE*)NULL,(*pvert)++);
           (*pnode)++;
         } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = n; 
    
   for (i = 0; i < na; i++){
      if (fabs(vertexes[i].x[0]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[0]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-0.0) < EPSA) vertexes[i].type = 1|NMASK;
      if (fabs(vertexes[i].x[1]-1.0) < EPSA) vertexes[i].type = 1|NMASK;
   }  
   for (i = 0; i < nvx-1; i++)
      for (j = 0; j < nvy-1; j++) {
            il = j + nvy*i;
            elem(il, il+nvy, il+nvy+1, il+1, nodes, pelement); 
   }
   (*pelement - 1)->succ = NULL;      
}

void cube_for_donut(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: cube_for_donut not available.\n");  }

void read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert;
NODE **pnode;
ELEMENT **pelement;
GRID *pgrid;
{
   INT i, i1, i2, i3, i4, type, nvert; 
   NODE *nodes;
   FILE *fp;
   
   FIRSTNODE(pgrid) = *pnode;
   FIRSTELEMENT(pgrid) = *pelement;
   nodes = *pnode;
   fp = fopen(COARSE_FILE,"r");
   fscanf(fp,"%i",&nvert);
   for(i = 1; i <= nvert; i++){
      fscanf(fp,"%lf %lf %i",&((*pvert)->x[0]),&((*pvert)->x[1]),&type);
      (*pvert)->type = (UNS)type;
      nmake(*pnode,*pnode+1,i,(NODE*)NULL,(*pvert)++);
      (*pnode)++;
   } 
   (*pnode-1)->succ = NULL;
   pgrid->minNodeIndex = 1;
   pgrid->maxNodeIndex = nvert; 
   fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   while (i1 > -1){
      elem(i1,i2,i3,i4,nodes,pelement);
      fscanf(fp,"%i %i %i %i",&i1,&i2,&i3,&i4);
   }
   (*pelement - 1)->succ = NULL;
}

void Read_vertexes_and_elements(pnode,pvert,pelement,pgrid)
VERTEX **pvert; NODE **pnode; ELEMENT **pelement; GRID *pgrid;
{  eprintf("Error: Read_vertexes_and_elements not available.\n");  }

#if F_DATA & CURVED_FACE_MIDDLE

void check_c_midpoints(mg)
MULTIGRID *mg;
{  eprintf("Error: check_c_midpoints not available.\n");  }

#else

void check_c_midpoints(mg)
MULTIGRID *mg;
{}

#endif

#endif

void new_node_order(pgrid,mask)  /*  place nodes with (type & mask) != 0     */
GRID *pgrid;                     /*  at the beginning of the list of nodes   */
INT mask;
{
   NODE *pn, *pnodeif, *pnodeil, *pnodebf, *pnodebl;
 
   pn = FIRSTNODE(pgrid);

   if (NTYPE(pn) & mask){
      pnodebf = pn;
      while (pn->succ && NTYPE(pn->succ) & mask) pn = pn->succ;
      pnodebl = pn;
      pnodeif = pnodeil = pn = pn->succ;
      if (pn)
         pnodebl->succ = pn->succ;
   }
   else{
      pnodeif = pn;
      while (pn->succ && !(NTYPE(pn->succ) & mask)) pn = pn->succ;
      pnodeil = pn;
      pnodebf = pnodebl = pn = pn->succ;
   }
   if (pn)
      pn = pn->succ;
   while (pn){
      if (NTYPE(pn) & mask) 
         pnodebl = pn;
      else{
         pnodeil = pnodeil->succ = pn;
         pnodebl->succ = pn->succ;
      }
      pn = pn->succ;
   }
   if (pnodebf){
      pgrid->firstN = pnodebf;
      pnodebl->succ = pnodeif;
   }
   else
      pgrid->firstN = pnodeif;
   if (pnodeil)
      pnodeil->succ = NULL;
   pgrid->firstNode = pnodeif;
   pgrid->lastNode = pnodeil;
   pgrid->fdbn = pnodebf;
   pgrid->ldbn = pnodebl;
}

#if (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void set_itypes_on_first_level(pgrid)
GRID *pgrid;
{
   ELEMENT *pel;
   NODE *pnode;
   INT i;
   
   for (pnode = FIRSTN(pgrid); pnode; pnode = pnode->succ)
      ITYPE(pnode) = 0;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         ITYPE(pel->n[i]) |= ITYPE(pel); 
}

void set_ibndry_face_bit(pgrid,itypes_pos)
GRID *pgrid;
ITYPES_POSITIONS *itypes_pos;
{
   ELEMENT *pel;
   FACE *pface;
   INT i, j, n=itypes_pos->n_e;
   
   for (pface = FIRSTF(pgrid); pface; pface = pface->succ)
      FTYPE(pface) = (FTYPE(pface) & ~IBNDRY_FACE_BIT) | AUX_FACE_BIT;
   for (i = 0; i < n; i++){
      for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
         if (ITYPE(pel) == itypes_pos->it_e[i])
            for (j = 0; j < SIDES; j++)
               if (FTYPE(pel->f[j]) & AUX_FACE_BIT){
                  if (FTYPE(pel->f[j]) & IBNDRY_FACE_BIT)
                     FTYPE(pel->f[j]) &= ~IBNDRY_FACE_BIT;
                  else
                     FTYPE(pel->f[j]) |= IBNDRY_FACE_BIT;
                     
               }
      for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ)
         if (ITYPE(pel) == itypes_pos->it_e[i])
            for (j = 0; j < SIDES; j++)
               FTYPE(pel->f[j]) &= ~AUX_FACE_BIT;
   }
}

void new_element_order2(pgrid,itypes_pos)
GRID *pgrid;
ITYPES_POSITIONS *itypes_pos;
{
   ELEMENT *pel, 
           *pe1[MAX_ITYPES]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}, 
           *pe2[MAX_ITYPES]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   INT it[MAX_ITYPES], jt[MAX_ITYPES]={0,0,0,0,0,0,0,0,0,0}, i, j, n;

   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ){
      i = ITYPE(pel);
      for (j = 0; j < itypes_pos->n && itypes_pos->it[j] != i; j++);
      if (j < itypes_pos->n)
         jt[j] = 1;
      else
         eprintf("Error: incompatible itype of element.\n");
   }
   for (i=0, n=0; i < itypes_pos->n; i++)
      if (jt[i])
         it[n++] = itypes_pos->it[i];
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ){
      GET_ITYPE_INDEX(i,ITYPE(pel),it);
      if (pe2[i])
         pe2[i] = (pe2[i])->succ = pel;
      else
         pe2[i] = pe1[i] = pel;
   }
   for (i = 1; i < n; i++)
      pe2[i-1]->succ = pe1[i];
   FIRSTELEMENT(pgrid) = pe1[0];
   pe2[n-1]->succ = NULL; 
   itypes_pos->n_e = n;
   for (i = 0; i < MAX_ITYPES; i++){
      itypes_pos->it_e[i] = it[i];
      itypes_pos->pe1[i] = pe1[i];
      itypes_pos->pe2[i] = pe2[i];
   }
}

void new_node_order2(pgrid,itypes_pos) /* orders the nodes according to itype */
GRID *pgrid;                    /*  itype of boundary nodes is first        */
ITYPES_POSITIONS *itypes_pos;   /*  order according to decreasing itype     */
{
   NODE *pnode,
        *pn1[MAX_ITYPES]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}, 
        *pn2[MAX_ITYPES]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   INT it[MAX_ITYPES]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, i, j, k, max, n=1;

   it[0] = ITYPE(FIRSTN(pgrid));
   for (pnode = FIRSTN(pgrid); pnode != FIRSTNODE(pgrid); pnode = pnode->succ)
      if (ITYPE(pnode) != it[0])
         eprintf("Error: the same itype assumed for all non-Dirichlet nodes.\n");
   for (pnode = FIRSTNODE(pgrid); pnode; pnode = pnode->succ){
      i = ITYPE(pnode);
      if (i!=it[0] && i!=it[1] && i!=it[2] && i!=it[3] && i!=it[4] &&
          i!=it[5] && i!=it[6] && i!=it[7] && i!=it[8] && i!=it[9]){
         if (n < MAX_ITYPES)
            it[n++] = i;
         else
            eprintf("Error: two many itypes of nodes.\n");
      }
   }
   for (i = 1; i < n; i++){
      max = it[i];
      k = i;
      for (j = i+1; j < n; j++)
         if (it[j] > max){
            max = it[j];
            k = j;
         }
      if (k != i){
         it[k] = it[i];
         it[i] = max;
      }
   }
   pgrid->ldbn = NULL;
   if (FIRSTN(pgrid) != FIRSTNODE(pgrid)){
      pn1[0] = FIRSTN(pgrid);
      for (pnode = FIRSTN(pgrid); pnode->succ != FIRSTNODE(pgrid); 
                                                         pnode = pnode->succ);
      pn2[0] = pgrid->ldbn = pnode;
   }
   for (pnode = FIRSTNODE(pgrid); pnode; pnode = pnode->succ){
      GET_ITYPE_INDEX(i,ITYPE(pnode),it);
      if (pn2[i])
         pn2[i] = (pn2[i])->succ = pnode;
      else
         pn2[i] = pn1[i] = pnode;
   }
   for (i = 1; i < n; i++)
      pn2[i-1]->succ = pn1[i];
   FIRSTN(pgrid) = pn1[0];
   pn2[n-1]->succ = NULL; 
   if (pgrid->ldbn)
      FIRSTNODE(pgrid) = pgrid->ldbn->succ;
   else
      FIRSTNODE(pgrid) = FIRSTN(pgrid);
   itypes_pos->n = n;
   for (i = 0; i < MAX_ITYPES; i++){
      itypes_pos->it[i] = it[i];
      itypes_pos->pn1[i] = pn1[i];
      itypes_pos->pn2[i] = pn2[i];
   }
}

void new_nnlink_order(pgrid,itypes_pos)
GRID *pgrid;
ITYPES_POSITIONS *itypes_pos;
{
   NODE *pnode;
   LINK *pli, *pl0, *pl1[MAX_ITYPES], *pl2[MAX_ITYPES];
   INT it[MAX_ITYPES], i, j, n=itypes_pos->n, nn, tstart_is_start;

   for (i = 0; i < n; i++)
      it[i] = itypes_pos->it[i];
   for (pnode = FIRSTN(pgrid); pnode; pnode = pnode->succ){
      pl1[0]=pl1[1]=pl1[2]=pl1[3]=pl1[4]=pl1[5]=pl1[6]=pl1[7]=pl1[8]=pl1[9]=NULL;
      if (pnode->start){
         if (pnode->tstart == pnode->start)
            tstart_is_start = 1;
         else{
            tstart_is_start = 0;
            pl1[0] = pnode->tstart;
            for (pli=pnode->tstart; pli->next != pnode->start; pli = pli->next);
            pl2[0] = pl0 = pli;
         }
         for (pli=pnode->start; pli; pli=pli->next){
            GET_ITYPE_INDEX(i,ITYPE(pli->nbnode),it);
            if (pl1[i])
               pl2[i] = (pl2[i])->next = pli;
            else
               pl2[i] = pl1[i] = pli;
         }
         if (tstart_is_start)
            for (i = 0; pl1[i] == NULL; i++);
         else
            i = 0;
         pnode->tstart = pl1[i];
         pli = pl2[i];
         for (i++; i < n; i++)
            if (pl1[i]){
               pli->next = pl1[i];
               pli = pl2[i];
            }
         pli->next = NULL;
         if (tstart_is_start)
            pnode->start = pnode->tstart;
         else
            pnode->start = pl0->next;
         if ((nn=itypes_pos->nn)+MAX_ITYPES+4 > N_NN_ITYPE)
            eprintf("Error: not enough entries in nn_link_itypes\n");
         nn_link_itypes[nn++] = pnode->tstart;
         nn_link_itypes[nn++] = pnode->start;
         if (!tstart_is_start)
            nn_link_itypes[nn++] = pl0;
         j = ITYPE(pnode->tstart->nbnode);
         for (pli=pnode->tstart; pli; pli=pli->next)
            if (pli->next)
               if(ITYPE(pli->next->nbnode) != j){
                  nn_link_itypes[nn++] = pli;
                  nn_link_itypes[nn++] = pli->next;
                  j = ITYPE(pli->next->nbnode);
               }
         if (!tstart_is_start || nn > itypes_pos->nn + 2){
            pnode->nn_i = itypes_pos->nn;
            pnode->nn_f = itypes_pos->nn = nn; 
         }
         else
            pnode->nn_f = 0;
      }
      else
         pnode->nn_f = 0;
   }
}

void first_and_last_nnlink_with_itype(pnode,itype,pl1,pl2)
NODE *pnode;
LINK **pl1, **pl2;
INT itype;
{
   LINK *pli;
   INT i, j;

   if (pnode->nn_f){
      if (pnode->tstart == pnode->start)
         i = j = pnode->nn_i + 2;
      else
         i = j = pnode->nn_i + 3;
      while (i < pnode->nn_f && ITYPE(nn_link_itypes[i]->nbnode) != itype)
         i += 2;
      if (i < pnode->nn_f){
         if (i == j)
            *pl1 = pnode->tstart;
         else
            *pl1 = nn_link_itypes[i-1];
         *pl2 = nn_link_itypes[i];
      }
      else if (ITYPE((*pl1=nn_link_itypes[--i])->nbnode) == itype){
         for (pli=*pl1; pli->next; pli=pli->next);
         *pl2 = pli; 
      }
      else
         *pl1 = NULL;
   }
   else if (ITYPE(pnode->tstart->nbnode) == itype){
      *pl1 = pnode->tstart;
      for (pli=*pl1; pli->next; pli=pli->next);
      *pl2 = pli;
   }
   else
      *pl1 = NULL;
}

void change_order_according_to_itypes(pgrid,itypes_pos)
GRID *pgrid;
ITYPES_POSITIONS *itypes_pos;
{
   ELEMENT *pel;
   NODE *pnode;

// set_itypes_on_first_level(pgrid);
   new_node_order2(pgrid,itypes_pos);
   new_element_order2(pgrid,itypes_pos);
   new_nnlink_order(pgrid,itypes_pos);
   set_ibndry_face_bit(pgrid,itypes_pos);

   itypes_pos->firstN = FIRSTN(pgrid);
   itypes_pos->firstNode = FIRSTNODE(pgrid);
   itypes_pos->ldbn = pgrid->ldbn;
}

void restrict_nodes_to_itypes(tGrid,nb,
              b_it1,b_it2,b_it3,b_it4,b_it5,b_it6,b_it7,b_it8,i_it1,itypes_pos)
GRID *tGrid;
ITYPES_POSITIONS *itypes_pos;
INT nb, b_it1, b_it2, b_it3, b_it4, b_it5, b_it6, b_it7, b_it8, i_it1;
{
   INT i, j, k, b_it[8]={b_it1,b_it2,b_it3,b_it4,b_it5,b_it6,b_it7,b_it8};

   GET_ITYPE_INDEX(i,b_it1,itypes_pos->it);
   FIRSTN(tGrid) = itypes_pos->pn1[i];
   for (k = 1; k < nb; k++){
      GET_ITYPE_INDEX(j,b_it[k],itypes_pos->it);
      itypes_pos->pn2[i]->succ = itypes_pos->pn1[j];
      i = j;
   }
   GET_ITYPE_INDEX(j,i_it1,itypes_pos->it);
   itypes_pos->pn2[i]->succ = FIRSTNODE(tGrid) = itypes_pos->pn1[j];
   itypes_pos->pn2[j]->succ = NULL;
}

void restrict_nnlinks_to_itypes(tGrid,nb,
                          b_it1,b_it2,b_it3,b_it4,b_it5,b_it6,b_it7,b_it8,i_it1)
GRID *tGrid;              /* DECREASING ORDER! */
INT nb, b_it1, b_it2, b_it3, b_it4, b_it5, b_it6, b_it7, b_it8, i_it1;
{
   NODE *pnode;
   LINK *pli1, *pli2, *pl1[8], *pl2[8];
   INT j, k, b_it[8]={b_it1,b_it2,b_it3,b_it4,b_it5,b_it6,b_it7,b_it8};

   for (pnode = FIRSTN(tGrid); pnode; pnode=SUCC(pnode)){
      first_and_last_nnlink_with_itype(pnode,i_it1,&pli1,&pli2);
      for (k = 0; k < nb; k++)
         first_and_last_nnlink_with_itype(pnode,b_it[k],&(pl1[k]),&(pl2[k]));
      for (k = 0; k < nb && pl1[k] == 0; k++);
      if (k < nb){
         pnode->tstart = pl1[k];
         j = k;
         for (k++; k < nb; k++)
            if (pl1[k]){
               pl2[j]->next = pl1[k];
               j = k;
            }
         pl2[j]->next = pli1;
      }
      else
         pnode->tstart = pli1;
      pnode->start = pli1;
      if (pli1)
         pli2->next = NULL;
   }
}

void restrict_elements_to_itype(tGrid,itype,itypes_pos)
GRID *tGrid;
ITYPES_POSITIONS *itypes_pos;
INT itype;
{
   INT i;

   GET_ITYPE_INDEX(i,itype,itypes_pos->it_e);
   FIRSTELEMENT(tGrid) = itypes_pos->pe1[i];
   itypes_pos->pe2[i]->succ = NULL;
}

void remove_restrictions_to_itypes(pgrid,itypes_pos)
GRID *pgrid;
ITYPES_POSITIONS *itypes_pos;
{
   NODE *pnode;
   LINK *pli;
   INT i, n;

   FIRSTN(pgrid) = itypes_pos->firstN;
   FIRSTNODE(pgrid) = itypes_pos->firstNode;
   for (i = 1; i < itypes_pos->n; i++)
      itypes_pos->pn2[i-1]->succ = itypes_pos->pn1[i];
   itypes_pos->pn2[itypes_pos->n - 1]->succ = NULL;
   if (itypes_pos->ldbn)
      itypes_pos->ldbn->succ = FIRSTNODE(pgrid);
   FIRSTELEMENT(pgrid) = itypes_pos->pe1[0];
   for (i = 1; i < itypes_pos->n_e; i++)
      itypes_pos->pe2[i-1]->succ = itypes_pos->pe1[i];
   itypes_pos->pe2[itypes_pos->n_e - 1]->succ = NULL;

   for (pnode = FIRSTN(pgrid); pnode; pnode = pnode->succ)
      if (pnode->nn_f){
         n = pnode->nn_i;
         pnode->tstart = nn_link_itypes[n++];
         pnode->start = nn_link_itypes[n++];
         if (pnode->tstart != pnode->start)
            nn_link_itypes[n++]->next = pnode->start;
         while (n < pnode->nn_f){
            nn_link_itypes[n]->next = nn_link_itypes[n+1];
            n += 2;
         }
      }
}

#else

void change_order_according_to_itypes(pgrid,itypes_pos)
GRID *pgrid; ITYPES_POSITIONS *itypes_pos;
{}

#endif

void new_face_order(pgrid,mask) /* place faces with (type & mask) != 0        */
GRID *pgrid;                    /* at the beginning of list of faces          */
INT mask;
{
   FACE *pf, *pfaceif, *pfaceil, *pfacebf, *pfacebl;

   pf = FIRSTFACE(pgrid);
   if (FTYPE(pf) & mask){
      pfacebf = pf;
      while (pf->succ && FTYPE(pf->succ) & mask) pf = pf->succ;
      pfacebl = pf;
      pfaceif = pfaceil = pf = pf->succ;
      if (pf)
         pfacebl->succ = pf->succ;
   }
   else{
      pfaceif = pf;
      while (pf->succ && !(FTYPE(pf->succ) & mask)) pf = pf->succ;
      pfaceil = pf;
      pfacebf = pfacebl = pf = pf->succ;
   }
   if (pf)
      pf = pf->succ;
   while (pf){
      if (FTYPE(pf) & mask) 
         pfacebl = pf;
      else{
         pfaceil = pfaceil->succ = pf;
         pfacebl->succ = pf->succ;
      }
      pf = pf->succ;
   }
   if (pfacebf){
      pgrid->firstF = pfacebf;
      pfacebl->succ = pfaceif;
   }
   else
      pgrid->firstF = pfaceif;
   pfaceil->succ = NULL;
   pgrid->firstFace = pfaceif;
   pgrid->lastFace = pfaceil;
   pgrid->fdbf = pfacebf;
   pgrid->ldbf = pfacebl;
} 

void new_element_order(pgrid,mask)/*  elements having a face with           */
GRID *pgrid;                      /*  (type & mask) != 0 are placed at the end*/
INT mask;                         /*  of the list of elements                 */
{
 ELEMENT *pe, *pelemif, *pelemil, *pelembf, *pelembl;
 
 pe = FIRSTELEMENT(pgrid);

 if (IS_FTYPE(pe,mask)){
   pelembf = pe;
   while (pe->succ && IS_FTYPE(pe->succ,mask)) 
      pe = pe->succ;
   pelembl = pe;
   pelemif = pelemil = pe = pe->succ;
   if (pe)
      pelembl->succ = pe->succ;
 }
 else{
   pelemif = pe;
   while (pe->succ && !IS_FTYPE(pe->succ,mask)) 
      pe = pe->succ;
   pelemil = pe;
   pelembf = pelembl = pe = pe->succ;
 }
 if (pe)
    pe = pe->succ;
 while (pe){
   if (IS_FTYPE(pe,mask))
     pelembl = pe;
   else{
     pelemil = pelemil->succ = pe;
     pelembl->succ = pe->succ;
   }
   pe = pe->succ;
 }
 if (pelemif){
    pelemil->succ = pelembf;
    FIRSTELEMENT(pgrid) = pelemif;
 }
 else
    FIRSTELEMENT(pgrid) = pelembf;
 pelembl->succ = NULL;
 FDBE(pgrid) = pelembf;
}

void el_neighbour(n,i,plink)
NODE *n;
LINK **plink;
INT i;
{
   LINK *pli;
   
   if (n->start == NULL) n->start = *plink;
   else {
      for (pli = n->start; pli->next != NULL; pli = pli->next);
      pli->next = *plink;
   }
   (*plink)->flag = i;
   ((*plink)++)->next = NULL;
}

#if DIM == 3

#if ELEMENT_TYPE == SIMPLEX

INT find3(i1,i2,i3,n)  /*   i1 < i2 < i3                              */
INT i1,i2,i3;          /*    (i1,i2,i3) == (j1,j2,j3) -> 0            */
NODE *n[NVERT];        /*                  (j0,j2,j3) -> 1            */
{                      /*                  (j0,j1,j3) -> 2            */
   INT j0,j1,j2,j3;    /*                  (j0,j1,j2) -> 3            */
 
   j0 = n[0]->index;
   j1 = n[1]->index;
   j2 = n[2]->index;
   j3 = n[3]->index;
 
   if (i1==j0){
      if (i2==j1){
         if (i3==j2)
            return(3);
         else if (i3==j3)
            return(2);
         else
            return(-1);
      }
      else if (i2==j2 && i3==j3)
         return(1);
      else
         return(-1);
   } 
   else if (i1==j1 && i2==j2 && i3==j3)
      return(0);
   else
      return(-1);
} 
 
void findel(i1,i2,i3,l,k,faceNumber,pelement,pface,firstelement)
INT i1,i2,i3,l,k,*faceNumber;/* l satisfies pelement->n[l]->index \in {i1,i2,i3} */
ELEMENT *pelement, *firstelement;  /*  pelement->f[k] is the face (i1,i2,i3)  */
FACE **pface;
{
   ELEMENT *pelem;
   LINK *pli;
   INT indf, ipel;
 
   ipel = (INT)(pelement - firstelement);
   for (pli = pelement->n[l]->start; pli != NULL &&
       ((indf = find3(i1,i2,i3,(pelem=firstelement + pli->flag)->n)) < 0 || 
                                       pli->flag == ipel); pli = pli->next);
   if (pli && pelem < pelement){
      pelement->f[k] = pelem->f[indf];
      FTYPE(pelement->f[k]) = 0;
   }
   else{
      fmake(*pface,*pface + 1,BNDRY_FACE_BIT,OR_SET,++(*faceNumber),(FACE*)NULL);
      pelement->f[k] = (*pface)++;
   }
}

void set_type_of_boundary_faces(pgrid,n_mask,f_mask)
GRID *pgrid;
INT n_mask, f_mask;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *f0, *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n0,n1,n2,n3,pel);
      FACES_OF_ELEMENT(f0,f1,f2,f3,pel);
      if (IS_BF(f0) && NTYPE(n1) & n_mask && NTYPE(n2) & n_mask && 
                       NTYPE(n3) & n_mask) FTYPE(f0) |= f_mask;
      if (IS_BF(f1) && NTYPE(n0) & n_mask && NTYPE(n2) & n_mask && 
                       NTYPE(n3) & n_mask) FTYPE(f1) |= f_mask;
      if (IS_BF(f2) && NTYPE(n0) & n_mask && NTYPE(n1) & n_mask && 
                       NTYPE(n3) & n_mask) FTYPE(f2) |= f_mask;
      if (IS_BF(f3) && NTYPE(n0) & n_mask && NTYPE(n1) & n_mask && 
                       NTYPE(n2) & n_mask) FTYPE(f3) |= f_mask;
   }
} 

void make_faces(pface,plink,pgrid,n_mask,f_mask)
FACE **pface;
LINK **plink;
GRID *pgrid;
INT n_mask, f_mask;
{
   INT i0, i1, i2, i3, nface=0, ipel; 
   ELEMENT *pel, *firstelement;
   NODE *pn;
   LINK *pl;

   firstelement = FIRSTELEMENT(pgrid);
   pgrid->firstFace = *pface;
   pl = *plink;

   for (pel = firstelement; pel != NULL; pel = pel->succ){
      ipel = (INT)(pel - firstelement);
      el_neighbour(pel->n[0],ipel,plink);
      el_neighbour(pel->n[1],ipel,plink);
      el_neighbour(pel->n[2],ipel,plink);
      el_neighbour(pel->n[3],ipel,plink);
   }
   for (pel = firstelement; pel != NULL; pel = pel->succ){
      i0 = pel->n[0]->index;
      i1 = pel->n[1]->index; 
      i2 = pel->n[2]->index;
      i3 = pel->n[3]->index;
      findel(i0,i1,i2,0,3,&nface,pel,pface,firstelement);
      findel(i0,i1,i3,0,2,&nface,pel,pface,firstelement);
      findel(i0,i2,i3,0,1,&nface,pel,pface,firstelement);
      findel(i1,i2,i3,1,0,&nface,pel,pface,firstelement);
   }
   (*pface-1)->succ = NULL;
   pgrid->minFaceIndex = 1;
   pgrid->maxFaceIndex = nface;
   for (pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
      pn->start = NULL;
   *plink = pl;
   set_type_of_boundary_faces(pgrid,n_mask,f_mask);
   new_face_order(pgrid,f_mask);
} 

#endif

#else   /*  if (DIM == 2)  */

INT find3();

void findel(i1,i2,l,k,faceNumber,pelement,pface,firstelement)
INT i1,i2,l,k,*faceNumber; /*  l satisfies pelement->n[l]->index \in {i1,i2} */
ELEMENT *pelement, *firstelement; /*  pelement->f[k] is the face (i1,i2)     */
FACE **pface;
{
   ELEMENT *pelem;
   LINK *pli;
   INT indf, ipel;
 
   ipel = (INT)(pelement - firstelement);
   for (pli = pelement->n[l]->start; pli != NULL &&
     ((indf = find3(i1,i2,(pelem=firstelement + pli->flag)->n)) < 0 || 
                                           pli->flag == ipel); pli = pli->next);
   if (pli && pelem < pelement){
     pelement->f[k] = pelem->f[indf];
     FTYPE(pelement->f[k]) = 0;
   }
   else{
     fmake(*pface,*pface + 1,BNDRY_FACE_BIT,OR_SET,++(*faceNumber),(FACE*)NULL);
     pelement->f[k] = (*pface)++;
   }
}

#if ELEMENT_TYPE == CUBE

INT find3(i1,i2,n)
INT i1,i2;          /*    (i1,i2) == (j0,j1) -> 0            */
NODE *n[NVERT];     /*               (j1,j2) -> 1            */
{                   /*               (j2,j3) -> 2            */
   INT j0,j1,j2,j3; /*               (j3,j0) -> 3            */
 
   j0 = n[0]->index;
   j1 = n[1]->index;
   j2 = n[2]->index;
   j3 = n[3]->index;
 
   if (SAME_COUPLES(i1,i2,j0,j1)) return(0);
   else if (SAME_COUPLES(i1,i2,j1,j2)) return(1);
   else if (SAME_COUPLES(i1,i2,j2,j3)) return(2);
   else if (SAME_COUPLES(i1,i2,j3,j0)) return(3);
   else return(-1);
} 
 
void set_type_of_boundary_faces(pgrid,n_mask,f_mask)
GRID *pgrid;
INT n_mask, f_mask;
{
   NODE *n0, *n1, *n2, *n3;
   FACE *f0, *f1, *f2, *f3;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      NODES_OF_4ELEMENT(n0,n1,n2,n3,pel);
      FACES_OF_4ELEMENT(f0,f1,f2,f3,pel);
      if (IS_BF(f0) && NTYPE(n0) & n_mask && NTYPE(n1) & n_mask)
                                                  FTYPE(f0) |= f_mask;
      if (IS_BF(f1) && NTYPE(n1) & n_mask && NTYPE(n2) & n_mask)
                                                  FTYPE(f1) |= f_mask;
      if (IS_BF(f2) && NTYPE(n2) & n_mask && NTYPE(n3) & n_mask)
                                                  FTYPE(f2) |= f_mask;
      if (IS_BF(f3) && NTYPE(n3) & n_mask && NTYPE(n0) & n_mask)
                                                  FTYPE(f3) |= f_mask;
   }
} 

void make_faces(pface,plink,pgrid,n_mask,f_mask)
FACE **pface;
LINK **plink;
GRID *pgrid;
INT n_mask, f_mask;
{
   INT i0, i1, i2, i3, nface=0, ipel; 
   ELEMENT *pel, *firstelement;
   NODE *pn;
   LINK *pl;

   firstelement = FIRSTELEMENT(pgrid);
   pgrid->firstFace = *pface;
   pl = *plink;

   for (pel = firstelement; pel != NULL; pel = pel->succ){
      ipel = (INT)(pel - firstelement);
      el_neighbour(pel->n[0],ipel,plink);
      el_neighbour(pel->n[1],ipel,plink);
      el_neighbour(pel->n[2],ipel,plink);
      el_neighbour(pel->n[3],ipel,plink);
   }
   for (pel = firstelement; pel != NULL; pel = pel->succ){
      i0 = pel->n[0]->index;
      i1 = pel->n[1]->index; 
      i2 = pel->n[2]->index;
      i3 = pel->n[3]->index;
      findel(i0,i1,0,0,&nface,pel,pface,firstelement);
      findel(i1,i2,1,1,&nface,pel,pface,firstelement);
      findel(i2,i3,2,2,&nface,pel,pface,firstelement);
      findel(i3,i0,3,3,&nface,pel,pface,firstelement);
   }
   (*pface-1)->succ = NULL;
   pgrid->minFaceIndex = 1;
   pgrid->maxFaceIndex = nface;
   for (pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
      pn->start = NULL;
   *plink = pl;
   set_type_of_boundary_faces(pgrid,n_mask,f_mask);
   new_face_order(pgrid,f_mask);
}

#elif ELEMENT_TYPE == SIMPLEX

INT find3(i1,i2,n)  /*   i1 < i2                             */
INT i1,i2;          /*    (i1,i2) == (j1,j2) -> 0            */
NODE *n[NVERT];     /*               (j0,j2) -> 1            */
{                   /*               (j0,j1) -> 2            */
   INT j0,j1,j2;
 
   j0 = n[0]->index;
   j1 = n[1]->index;
   j2 = n[2]->index;
 
   if (i1==j0){
      if (i2==j1)
         return(2);
      else if (i2==j2)
         return(1);
      else
         return(-1);
   } 
   else if (i1==j1 && i2==j2)
      return(0);
   else
      return(-1);
} 
 
void set_type_of_boundary_faces(pgrid,n_mask,f_mask)
GRID *pgrid;
INT n_mask, f_mask;
{
   NODE *n0, *n1, *n2;
   FACE *f0, *f1, *f2;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(f0,f1,f2,pel);
      if (IS_BF(f0) && NTYPE(n1) & n_mask && NTYPE(n2) & n_mask)
                                                  FTYPE(f0) |= f_mask;
      if (IS_BF(f1) && NTYPE(n0) & n_mask && NTYPE(n2) & n_mask)
                                                  FTYPE(f1) |= f_mask;
      if (IS_BF(f2) && NTYPE(n0) & n_mask && NTYPE(n1) & n_mask)
                                                  FTYPE(f2) |= f_mask;
   }
} 

void set_type_of_boundary_faces2(pgrid,n_mask1,n_mask2,f_mask)
GRID *pgrid;
INT n_mask1, n_mask2, f_mask;
{
   NODE *n0, *n1, *n2;
   FACE *f0, *f1, *f2;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(f0,f1,f2,pel);
      if (IS_BF(f0) && NTYPE(n1) & n_mask1 && NTYPE(n2) & n_mask2)
													FTYPE(f0) |= f_mask;
      if (IS_BF(f1) && NTYPE(n0) & n_mask1 && NTYPE(n2) & n_mask2)
													FTYPE(f1) |= f_mask;
      if (IS_BF(f2) && NTYPE(n0) & n_mask1 && NTYPE(n1) & n_mask2)
													FTYPE(f2) |= f_mask;
	  if (IS_BF(f0) && NTYPE(n1) & n_mask2 && NTYPE(n2) & n_mask1)
													FTYPE(f0) |= f_mask;
      if (IS_BF(f1) && NTYPE(n0) & n_mask2 && NTYPE(n2) & n_mask1)
													FTYPE(f1) |= f_mask;
      if (IS_BF(f2) && NTYPE(n0) & n_mask2 && NTYPE(n1) & n_mask1)
													FTYPE(f2) |= f_mask;
   }
}

void make_faces(pface,plink,pgrid,n_mask,f_mask)
FACE **pface;
LINK **plink;
GRID *pgrid;
INT n_mask, f_mask;
{
   INT i0, i1, i2, nface=0, ipel; 
   ELEMENT *pel, *firstelement;
   NODE *pn;
   LINK *pl;

   firstelement = FIRSTELEMENT(pgrid);
   pgrid->firstFace = *pface;
   pl = *plink;

   for (pel = firstelement; pel != NULL; pel = pel->succ){
      ipel = (INT)(pel - firstelement);
      el_neighbour(pel->n[0],ipel,plink);
      el_neighbour(pel->n[1],ipel,plink);
      el_neighbour(pel->n[2],ipel,plink);
   }
   for (pel = firstelement; pel != NULL; pel = pel->succ){
      i0 = pel->n[0]->index;
      i1 = pel->n[1]->index; 
      i2 = pel->n[2]->index;
      findel(i0,i1,0,2,&nface,pel,pface,firstelement);
      findel(i0,i2,0,1,&nface,pel,pface,firstelement);
      findel(i1,i2,1,0,&nface,pel,pface,firstelement);
   }
   (*pface-1)->succ = NULL;
   pgrid->minFaceIndex = 1;
   pgrid->maxFaceIndex = nface;
   for (pn = FIRSTN(pgrid); pn != NULL; pn = pn->succ)
      pn->start = NULL;
   *plink = pl;
   set_type_of_boundary_faces(pgrid,n_mask,f_mask);
   new_face_order(pgrid,f_mask);
}

#endif  /*  ELEMENT_TYPE == SIMPLEX  */

#endif  /*  DIM == 2  */
 
void nodesToNodes(plink,pelem,n_mask)
LINK **plink;
ELEMENT *pelem;
INT n_mask;
{
   LINK *pl, *plif, *plil, *plbf, *plbl;
   INT i, j;

   for (i = 0; i < NVERT; i++){
      plif = plbf = NULL;
      for (j = 0; j < NVERT; j++)
         if (j != i){
            if (pelem->n[i]->tstart)
               for (pl = pelem->n[i]->tstart; pl->nbnode != pelem->n[j] &&
                                                       pl->next; pl = pl->next);
            else
               pl = NULL;
            if ((pl == NULL) || pl->nbnode != pelem->n[j]){
               if (NTYPE(pelem->n[j]) & n_mask){
                  if (plbf == NULL)
                     plbf = plbl = *plink;
                  else
                     plbl = plbl->next = *plink;
               }
               else{
                  if (plif == NULL)
                     plif = plil = *plink;
                  else
                     plil = plil->next = *plink;
               }
               (*plink)->flag = 0;
               ((*plink)++)->nbnode = pelem->n[j];
            }
         }
      if (pelem->n[i]->tstart == NULL){
         if (plbf == NULL)
            plbf = plif;
         else
            plbl->next = plif;
         if (plif) 
            plil->next = NULL;
         pelem->n[i]->tstart = plbf;
         pelem->n[i]->start = plif;
      }
      else{
         if (plif){
            if (pelem->n[i]->tstart == pelem->n[i]->start)
               pelem->n[i]->tstart = plif;
            else{
               for (pl = pelem->n[i]->tstart; pl->next != pelem->n[i]->start; 
                                                                 pl = pl->next);
               pl->next = plif;
            }
            plil->next = pelem->n[i]->start;
            pelem->n[i]->start = plif;
         }
         if (plbf){
            plbl->next = pelem->n[i]->tstart;
            pelem->n[i]->tstart = plbf;
         }
      }
   }
}

#if DATA_S & F_LINK_TO_NODES

#if ELEMENT_TYPE == SIMPLEX

void nodesToFaces(pfnlink,pelem,n_mask,f_mask)
FNLINK **pfnlink;
ELEMENT *pelem;
INT n_mask, f_mask;
{
   FNLINK *pfnl, *pfnlif, *pfnlil, *pfnlbf, *pfnlbl;
   INT i, j;

   for (i = 0; i < SIDES; i++)
      if (pelem->f[i]->tfnstart == NULL){
         pfnlif = pfnlbf = NULL;
         for (j = 0; j < NVERT; j++){
            if (NTYPE(pelem->n[j]) & n_mask){
               if (pfnlbf == NULL)
                  pfnlbf = pfnlbl = *pfnlink;
               else
                  pfnlbl = pfnlbl->next = *pfnlink;
            }
            else{
               if (pfnlif == NULL)
                  pfnlif = pfnlil = *pfnlink;
               else
                  pfnlil = pfnlil->next = *pfnlink;
            }
            ((*pfnlink)++)->nbnode = pelem->n[j];
         }
         if (pfnlbf == NULL)
            pfnlbf = pfnlif;
         else
            pfnlbl->next = pfnlif;
         if (pfnlif) 
            pfnlil->next = NULL;
         pelem->f[i]->tfnstart = pfnlbf;
         pelem->f[i]->fnstart = pfnlif;
      }
      else{
         if (NTYPE(pelem->n[i]) & n_mask){
            (*pfnlink)->next = pelem->f[i]->tfnstart;
            pelem->f[i]->tfnstart = *pfnlink;
         }
         else{
            if (pelem->f[i]->tfnstart == pelem->f[i]->fnstart)
               pelem->f[i]->tfnstart = *pfnlink;
            else{
               for (pfnl = pelem->f[i]->tfnstart; 
                         pfnl->next != pelem->f[i]->fnstart; pfnl = pfnl->next);
               pfnl->next = *pfnlink;
            }
            (*pfnlink)->next = pelem->f[i]->fnstart;
            pelem->f[i]->fnstart = *pfnlink;
         }
         ((*pfnlink)++)->nbnode = pelem->n[i];
      }
}
 
#else  /*  if ELEMENT_TYPE != SIMPLEX  */

void nodesToFaces(pfnlink,pelem,n_mask,f_mask)
FNLINK **pfnlink;
ELEMENT *pelem;
INT n_mask, f_mask;
{
   FNLINK *pfnl, *pfnlif, *pfnlil, *pfnlbf, *pfnlbl;
   INT i, j;

   for (i = 0; i < SIDES; i++){
      pfnlif = pfnlbf = NULL;
      for (j = 0; j < NVERT; j++){
         if (pelem->f[i]->tfnstart)
            for (pfnl = pelem->f[i]->tfnstart; pfnl->nbnode != pelem->n[j] &&
                                                 pfnl->next; pfnl = pfnl->next);
         else
            pfnl = NULL;
         if ((pfnl == NULL) || pfnl->nbnode != pelem->n[j]){
            if (NTYPE(pelem->n[j]) & n_mask){
               if (pfnlbf == NULL)
                  pfnlbf = pfnlbl = *pfnlink;
               else
                  pfnlbl = pfnlbl->next = *pfnlink;
            }
            else{
               if (pfnlif == NULL)
                  pfnlif = pfnlil = *pfnlink;
               else
                  pfnlil = pfnlil->next = *pfnlink;
            }
            ((*pfnlink)++)->nbnode = pelem->n[j];
         }
      }
      if (pelem->f[i]->tfnstart == NULL){
         if (pfnlbf == NULL)
            pfnlbf = pfnlif;
         else
            pfnlbl->next = pfnlif;
         if (pfnlif) 
            pfnlil->next = NULL;
         pelem->f[i]->tfnstart = pfnlbf;
         pelem->f[i]->fnstart = pfnlif;
      }
      else{
         if (pfnlif){
            if (pelem->f[i]->tfnstart == pelem->f[i]->fnstart)
               pelem->f[i]->tfnstart = pfnlif;
            else{
               for (pfnl = pelem->f[i]->tfnstart; 
                         pfnl->next != pelem->f[i]->fnstart; pfnl = pfnl->next);
               pfnl->next = pfnlif;
            }
            pfnlil->next = pelem->f[i]->fnstart;
            pelem->f[i]->fnstart = pfnlif;
         }
         if (pfnlbf){
            pfnlbl->next = pelem->f[i]->tfnstart;
            pelem->f[i]->tfnstart = pfnlbf;
         }
      }
   }
}
 
#endif

#else  /*  if !(DATA_S & F_LINK_TO_NODES)  */

void nodesToFaces(pfnlink,pelem,n_mask,f_mask)
FNLINK **pfnlink; ELEMENT *pelem; INT n_mask, f_mask;
{}

#endif
 
#if DATA_S & N_LINK_TO_FACES

void facesToNodes(pnflink,pelem,n_mask,f_mask)
NFLINK **pnflink;
ELEMENT *pelem;
INT n_mask, f_mask;
{
   NFLINK *pnfl, *pnflif, *pnflil, *pnflbf, *pnflbl;
   INT i, j;

   for (i = 0; i < NVERT; i++){
      pnflif = pnflbf = NULL;
      for (j = 0; j < SIDES; j++){
         if (pelem->n[i]->tnfstart)
            for (pnfl = pelem->n[i]->tnfstart; pnfl->nbface != pelem->f[j] &&
                                                 pnfl->next; pnfl = pnfl->next);
         else
            pnfl = NULL;
         if ((pnfl == NULL) || pnfl->nbface != pelem->f[j]){
            if (FTYPE(pelem->f[j]) & f_mask){
               if (pnflbf == NULL)
                  pnflbf = pnflbl = *pnflink;
               else
                  pnflbl = pnflbl->next = *pnflink;
            }
            else{
               if (pnflif == NULL)
                  pnflif = pnflil = *pnflink;
               else
                  pnflil = pnflil->next = *pnflink;
            }
            ((*pnflink)++)->nbface = pelem->f[j];
         }
      }
      if (pelem->n[i]->tnfstart == NULL){
         if (pnflbf == NULL)
            pnflbf = pnflif;
         else
            pnflbl->next = pnflif;
         if (pnflif) 
            pnflil->next = NULL;
         pelem->n[i]->tnfstart = pnflbf;
         pelem->n[i]->nfstart = pnflif;
      }
      else{
         if (pnflif){
            if (pelem->n[i]->tnfstart == pelem->n[i]->nfstart)
               pelem->n[i]->tnfstart = pnflif;
            else{
               for (pnfl = pelem->n[i]->tnfstart;
                         pnfl->next != pelem->n[i]->nfstart; pnfl = pnfl->next);
               pnfl->next = pnflif;
            }
            pnflil->next = pelem->n[i]->nfstart;
            pelem->n[i]->nfstart = pnflif;
         }
         if (pnflbf){
            pnflbl->next = pelem->n[i]->tnfstart;
            pelem->n[i]->tnfstart = pnflbf;
         }
      }
   }
}

#else

void facesToNodes(pnflink,pel,n_mask,f_mask)
NFLINK **pnflink; ELEMENT *pel; INT n_mask, f_mask;
{}

#endif
 
#if DATA_S & F_LINK_TO_FACES

void facesToFaces(pflink,pelem,f_mask)
FLINK **pflink;
ELEMENT *pelem;
INT f_mask;
{
   FLINK *pfl, *pflif, *pflil, *pflbf, *pflbl;
   INT i, j;

   for (i = 0; i < SIDES; i++){
      pflif = pflbf = NULL;
      for (j = 0; j < SIDES; j++)
         if (j != i){
            if (FTYPE(pelem->f[j]) & f_mask){
               if (pflbf == NULL)
                  pflbf = pflbl = *pflink;
               else
                  pflbl = pflbl->next = *pflink;
            }
            else{
               if (pflif == NULL)
                  pflif = pflil = *pflink;
               else
                  pflil = pflil->next = *pflink;
            }
            ((*pflink)++)->nbface = pelem->f[j];
         }
      if (pelem->f[i]->tfstart == NULL){
         if (pflbf == NULL)
            pflbf = pflif;
         else
            pflbl->next = pflif;
         if (pflif) 
            pflil->next = NULL;
         pelem->f[i]->tfstart = pflbf;
         pelem->f[i]->fstart = pflif;
      }
      else{
         if (pflif){
            if (pelem->f[i]->tfstart == pelem->f[i]->fstart)
               pelem->f[i]->tfstart = pflif;
            else{
               for (pfl = pelem->f[i]->tfstart;
                             pfl->next != pelem->f[i]->fstart; pfl = pfl->next);
               pfl->next = pflif;
            }
            pflil->next = pelem->f[i]->fstart;
            pelem->f[i]->fstart = pflif;
         }
         if (pflbf){
            pflbl->next = pelem->f[i]->tfstart;
            pelem->f[i]->tfstart = pflbf;
         }
      }
   }
}

#else  /*  if !(DATA_S & F_LINK_TO_FACES)  */

void facesToFaces(pflink,pel,f_mask)
FLINK **pflink; ELEMENT *pel; INT f_mask;
{}  

#endif

#if DATA_S & N_LINK_TO_ELEMENTS

void elementsToNodes(pnelink,pelem)
NELINK **pnelink;
ELEMENT *pelem;
{
   NELINK *pnel;
   INT i;

   for (i = 0; i < NVERT; i++){
      if (pelem->n[i]->nestart){
         for (pnel = pelem->n[i]->nestart; pnel->next; pnel = pnel->next);
         pnel->next = *pnelink;
      }
      else
         pelem->n[i]->nestart = *pnelink;
      (*pnelink)->nbel = pelem;
      ((*pnelink)++)->next = NULL;
   }
}

#else

void elementsToNodes(pnelink,pelem)
NELINK **pnelink; ELEMENT *pelem;
{}

#endif

#if DIM == 3

#if ELEMENT_TYPE == SIMPLEX

FLOAT normal_vector(x1,x2,x3,f,n)  /* the ordering of the vertices x1, x2, x3 */
FLOAT *x1, *x2, *x3, *n;           /* of the face f has to correspond to the  */
FACE *f;                           /* increasing node index                   */
{
   FLOAT a[3], b[3], d;
  
   SUBTR(x2,x1,a);
   SUBTR(x3,x2,b);
   n[0] = a[1]*b[2]-a[2]*b[1];
   n[1] = a[2]*b[0]-a[0]*b[2];
   n[2] = a[0]*b[1]-a[1]*b[0];
   d = sqrt(DOT(n,n));
   if (FTYPE(f) & ORIENT){ 
      n[0] = -n[0]/d;
      n[1] = -n[1]/d;
      n[2] = -n[2]/d;
   }
   else{
      n[0] = n[0]/d;
      n[1] = n[1]/d;
      n[2] = n[2]/d;
   }
   return(d/2.0); /* d/2.0 is the area of triangle determined by x1, x2, x3  */
}

/* nn is normal vector to the face (x1,x2,x3) of tetrahedron (x0,x1,x2,x3).  */
void test_normal_vector_direction(nn,x0,x1,area)
FLOAT  *nn, *x0, *x1, *area;
{
   if ( (x0[0]-x1[0])*nn[0]+(x0[1]-x1[1])*nn[1]+(x0[2]-x1[2])*nn[2] > 0)
      *area = - *area;     /*   nn directed to the interior of (x1,x2,x3,x0) */
}     

FLOAT test_dir(i,pel)  /*  test_dir > 0 ... orientation of normal vector to  */
ELEMENT *pel;          /*                       pel->f[i] has to be changed  */
INT i;
{
   FLOAT *x0, *x1, *x2, *x3, n[3], a = -1.0;
   
   VERTICES_OF_ELEMENT(x0,x1,x2,x3,pel);
   switch (i){
      case 0: normal_vector(x1,x2,x3,pel->f[i],n);
              break;
      case 1: normal_vector(x2,x3,x0,pel->f[i],n);
              break;
      case 2: normal_vector(x3,x0,x1,pel->f[i],n);
              break;
      case 3: normal_vector(x0,x1,x2,pel->f[i],n);
   }
   test_normal_vector_direction(n,pel->n[i]->myvertex->x,
                                  pel->n[MAX(1-i,0)]->myvertex->x,&a); 
   return(a); 
}

void mark_edges_of_boundary_face(n1,n2,n3)
NODE *n1, *n2, *n3;                   /*  n1->index < n2->index < n3->index  */
{
   LINK *pli;
   
   for (pli = n1->tstart; pli->nbnode != n2; pli = pli->next); 
   MARK_BND_EDGE(pli);
   for (pli = n1->tstart; pli->nbnode != n3; pli = pli->next); 
   MARK_BND_EDGE(pli);
   for (pli = n2->tstart; pli->nbnode != n3; pli = pli->next); 
   MARK_BND_EDGE(pli);
}
 
void mark_boundary_edges(theGrid) 
GRID *theGrid;
{
   ELEMENT *pelem;
   
   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
      if (IS_BF(pelem->f[0]))
         mark_edges_of_boundary_face(pelem->n[1],pelem->n[2],pelem->n[3]);
      if (IS_BF(pelem->f[1]))
         mark_edges_of_boundary_face(pelem->n[0],pelem->n[2],pelem->n[3]);
      if (IS_BF(pelem->f[2]))
         mark_edges_of_boundary_face(pelem->n[0],pelem->n[1],pelem->n[3]);
      if (IS_BF(pelem->f[3]))
         mark_edges_of_boundary_face(pelem->n[0],pelem->n[1],pelem->n[2]);
   }
}         

#endif  /*  ELEMENT_TYPE == SIMPLEX  */

#else  /* if DIM == 2 */

FLOAT normal_vector(x1,x2,f,n)  /*  the ordering of the vertices x1, x2 of  */
FLOAT *x1, *x2, *n;             /*  the face f has to correspond to the     */
FACE *f;                        /*  increasing node index                   */
{
   FLOAT a[DIM], d;
  
   SUBTR(x2,x1,a);
   d = sqrt(DOT(a,a));
   if (FTYPE(f) & ORIENT){ 
      n[0] =  a[1]/d;
      n[1] = -a[0]/d;
   }
   else{
      n[0] = -a[1]/d;
      n[1] =  a[0]/d;
   }
   return(d); /* d is the length of f */
}

/* nn is normal vector to the face (x1,x2) of the triangle (x0,x1,x2).  */
void test_normal_vector_direction(nn,x0,x1,area)
FLOAT  *nn, *x0, *x1, *area;
{
   if ( (x0[0]-x1[0])*nn[0]+(x0[1]-x1[1])*nn[1] > 0)
      *area = - *area;   /*  nn directed to the interior of (x1,x2,x0)    */
}     

void outer_nonunit_normal_vector_to_element_edge(n0,n1,n2,nn)
NODE *n0, *n1, *n2; /* nodes of an element */
DOUBLE *nn;         /* outer normal vector to n0,n1; |nn|=|n0,n1|  */
{
   DOUBLE ss[DIM];

   nn[0] = n1->myvertex->x[1] - n0->myvertex->x[1];
   nn[1] = n0->myvertex->x[0] - n1->myvertex->x[0];
   ss[0] = n0->myvertex->x[0] - n2->myvertex->x[0];
   ss[1] = n0->myvertex->x[1] - n2->myvertex->x[1];
   if (DOT(nn,ss) < 0.)
      SET8(nn,nn)
}

void mark_edges_of_boundary_face(n1,n2)  /*  n1->index < n2->index  */
NODE *n1, *n2;
{
   LINK *pli;

   for (pli = n1->tstart; pli->nbnode != n2; pli = pli->next);
   MARK_BND_EDGE(pli);
}

#if ELEMENT_TYPE == SIMPLEX

FLOAT test_dir(i,pel)  /*  test_dir > 0 ... orientation of normal vector to  */
ELEMENT *pel;          /*                       pel->f[i] has to be changed  */
INT i;
{
   FLOAT *x0, *x1, *x2, n[DIM], a = -1.0;
   
   VERTICES_OF_ELEMENT(x0,x1,x2,pel);
   switch (i){
      case 0: normal_vector(x1,x2,pel->f[i],n);
              break;
      case 1: normal_vector(x0,x2,pel->f[i],n);
              break;
      case 2: normal_vector(x0,x1,pel->f[i],n);
   }
   test_normal_vector_direction(n,pel->n[i]->myvertex->x,
                                  pel->n[MAX(1-i,0)]->myvertex->x,&a); 
   return(a); 
}

void mark_boundary_edges(theGrid)
GRID *theGrid;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
      if (IS_BF(pelem->f[0]))
         mark_edges_of_boundary_face(pelem->n[1],pelem->n[2]);
      if (IS_BF(pelem->f[1]))
         mark_edges_of_boundary_face(pelem->n[0],pelem->n[2]);
      if (IS_BF(pelem->f[2]))
         mark_edges_of_boundary_face(pelem->n[0],pelem->n[1]);
   }
}

#elif ELEMENT_TYPE == CUBE

FLOAT test_dir(i,pel)  /*  test_dir > 0 ... orientation of normal vector to  */
ELEMENT *pel;          /*                       pel->f[i] has to be changed  */
INT i;
{
   NODE *n0, *n1;
   FLOAT n[DIM], a = -1.0;
   INT j;
   
   n0 = pel->n[i];
   if (i < 3)
      n1 = pel->n[i+1];
   else
      n1 = pel->n[0];
   if (INDEX(n0) > INDEX(n1)){
      n0 = n1;
      n1 = pel->n[i];
   }
   normal_vector(n0->myvertex->x,n1->myvertex->x,pel->f[i],n);
   if (i)
      j = i-1;
   else
      j = 3;
   test_normal_vector_direction(n,pel->n[j]->myvertex->x,
                                  pel->n[i]->myvertex->x,&a); 
   return(a); 
}

void mark_boundary_edges(theGrid)
GRID *theGrid;
{
   ELEMENT *pelem;

   for (pelem = FIRSTELEMENT(theGrid); pelem != NULL; pelem = pelem->succ){
      if (IS_BF(pelem->f[0]))
         if (INDEX(pelem->n[0]) < INDEX(pelem->n[1]))
            mark_edges_of_boundary_face(pelem->n[0],pelem->n[1]);
         else
            mark_edges_of_boundary_face(pelem->n[1],pelem->n[0]);
      if (IS_BF(pelem->f[1]))
         if (INDEX(pelem->n[1]) < INDEX(pelem->n[2]))
            mark_edges_of_boundary_face(pelem->n[1],pelem->n[2]);
         else
            mark_edges_of_boundary_face(pelem->n[2],pelem->n[1]);
      if (IS_BF(pelem->f[2]))
         if (INDEX(pelem->n[2]) < INDEX(pelem->n[3]))
            mark_edges_of_boundary_face(pelem->n[2],pelem->n[3]);
         else
            mark_edges_of_boundary_face(pelem->n[3],pelem->n[2]);
      if (IS_BF(pelem->f[3]))
         if (INDEX(pelem->n[3]) < INDEX(pelem->n[0]))
            mark_edges_of_boundary_face(pelem->n[3],pelem->n[0]);
         else
            mark_edges_of_boundary_face(pelem->n[0],pelem->n[3]);
   }
}

#endif  /*  ELEMENT_TYPE == CUBE  */

#endif

void set_outer_normal_vector(theGrid)
GRID *theGrid;
{
   ELEMENT *pel;
   INT i;
   
   for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i < SIDES; i++)
         if (IS_BF(pel->f[i]))
            if (test_dir(i,pel) > 0.){
               if (FTYPE(pel->f[i]) & ORIENT)
                  SET_ORIENT(pel->f[i],0);
               else
                  SET_ORIENT(pel->f[i],ORIENT);
            }
}

#if DIM == 2

FLOAT barycentric_coordinates(x0,x1,x2,b)
FLOAT *x0, *x1, *x2, b[DIM2][DIM2];
{
   FLOAT deta, a[2][2];
  
   a[0][0] = x0[0] - x2[0];
   a[0][1] = x1[0] - x2[0];
   a[1][0] = x0[1] - x2[1];
   a[1][1] = x1[1] - x2[1];
  
   deta = a[0][0]*a[1][1]-a[0][1]*a[1][0];
  
   if ( ( (deta>0)?deta:-deta ) < 1.e-15 ){
      eprintf("ZERO DETERMINANT\n");
      return(1.);
   }
  
   b[0][0] =  a[1][1]/deta;
   b[0][1] = -a[0][1]/deta;
   b[1][0] = -a[1][0]/deta;
   b[1][1] =  a[0][0]/deta;
   b[2][0] = -b[0][0]-b[1][0];
   b[2][1] = -b[0][1]-b[1][1];
   b[0][2] = -DOT(b[0],x2);
   b[1][2] = -DOT(b[1],x2);
   b[2][2] = -DOT(b[2],x0);
   
   return(fabs(deta)/2.0);  /* area of the triangle */
}

/* The reference mapping is   a[0][0]*x[0] + a[0][1]*x[1] + c[0]
                              a[1][0]*x[0] + a[1][1]*x[1] + c[1]  */
FLOAT P1_reference_mapping0(n0,n1,n2,a,c)
NODE *n0, *n1, *n2;
FLOAT a[DIM][DIM], c[DIM];
{
   FLOAT *x0, *x1, *x2, vol;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   c[0] = x0[0];
   c[1] = x0[1];
   vol = fabs(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
   if (vol < 1.e-10)
      eprintf("Error in P1_reference_mapping0\n");
   return(vol*0.5); /* area of the triangle */
}

FLOAT P1_reference_mapping0_with_inverse(n0,n1,n2,b,a,c)  /*  b =  (a)^{-1}  */
NODE *n0, *n1, *n2;
FLOAT b[DIM][DIM], a[DIM][DIM], c[DIM];
{
   FLOAT *x0, *x1, *x2, vol;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   c[0] = x0[0];
   c[1] = x0[1];
   vol = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   if (fabs(vol) < 1.e-10)
      eprintf("Error in P1_reference_mapping0_with_inverse\n");
   else{
      b[0][0] =   a[1][1]/vol;
      b[0][1] =  -a[0][1]/vol;
      b[1][0] =  -a[1][0]/vol;
      b[1][1] =   a[0][0]/vol;
   }
   return(fabs(vol)*0.5); /* area of the triangle */
}

#define QTEST_ZERO_JACOBIAN(jac,s0,s1,s2,s3)                                   \
   s0 = -jac[0] - jac[1] + jac[2];                                             \
   s1 =  jac[0] - jac[1] + jac[2];                                             \
   s2 =  jac[0] + jac[1] + jac[2];                                             \
   s3 = -jac[0] + jac[1] + jac[2];                                             \
   if ( fabs(s0) < 1.e-7 || fabs(s1) < 1.e-7 || fabs(s2) < 1.e-7 ||            \
        fabs(s3) < 1.e-7 || s0*s1 < 1.e-14 || s1*s2 < 1.e-14 ||                \
        s2*s3 < 1.e-14 || s3*s0 < 1.e-14)                                      \
      eprintf("ZERO JACOBIAN\n");

/* The reference mapping is 
                     a[0][0]*x[0] + a[0][1]*x[1] + c[0] + alpha[0]*x[0]*x[1]
                     a[1][0]*x[0] + a[1][1]*x[1] + c[1] + alpha[1]*x[0]*x[1]  */
/* The reference element is the square crossing the axes at +1 and -1.        */
/* jac = jac[0]*x + jac[1]*y + jac[2] is the Jacobian;
   (b[i][j][0]*x + b[i][j][1]*y + b[i][j][2])/jac is the inverse of the 
   Jacobi matrix of the reference mapping */
void Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac)
NODE *n0, *n1, *n2, *n3;
FLOAT b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, *x3, s0, s1, s2, s3;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   a[0][0]  = (-x0[0] + x1[0] + x2[0] - x3[0])*0.25;
   a[1][0]  = (-x0[1] + x1[1] + x2[1] - x3[1])*0.25;
   a[0][1]  = (-x0[0] - x1[0] + x2[0] + x3[0])*0.25;
   a[1][1]  = (-x0[1] - x1[1] + x2[1] + x3[1])*0.25;
   alpha[0] = ( x0[0] - x1[0] + x2[0] - x3[0])*0.25;
   alpha[1] = ( x0[1] - x1[1] + x2[1] - x3[1])*0.25;
   c[0]     = ( x0[0] + x1[0] + x2[0] + x3[0])*0.25;
   c[1]     = ( x0[1] + x1[1] + x2[1] + x3[1])*0.25;

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   QTEST_ZERO_JACOBIAN(jac,s0,s1,s2,s3)
  
   b[0][0][0] =  alpha[1];
   b[0][0][1] =  0.;
   b[0][0][2] =  a[1][1];
   b[0][1][0] =  -alpha[0];
   b[0][1][1] =  0.;
   b[0][1][2] =  -a[0][1];
   b[1][0][0] =  0.;
   b[1][0][1] =  -alpha[1];
   b[1][0][2] =  -a[1][0];
   b[1][1][0] =  0.;
   b[1][1][1] =  alpha[0];
   b[1][1][2] =  a[0][0];
}

void inverse_of_Q1_reference_mapping(n0,n1,n2,n3,b,jac)
NODE *n0, *n1, *n2, *n3;
FLOAT b[2][2][DIM2], jac[DIM2];
{
   FLOAT a[DIM][DIM], c[DIM], alpha[DIM];

   Q1_reference_mapping_with_inverse(n0,n1,n2,n3,b,a,c,alpha,jac);
}

void Q1_reference_mapping(n0,n1,n2,n3,a,alpha,c,jac)
NODE *n0, *n1, *n2, *n3;
FLOAT a[DIM][DIM], alpha[DIM], c[DIM], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, *x3, s0, s1, s2, s3;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   x3 = n3->myvertex->x;
   a[0][0]  = (-x0[0] + x1[0] + x2[0] - x3[0])*0.25;
   a[1][0]  = (-x0[1] + x1[1] + x2[1] - x3[1])*0.25;
   a[0][1]  = (-x0[0] - x1[0] + x2[0] + x3[0])*0.25;
   a[1][1]  = (-x0[1] - x1[1] + x2[1] + x3[1])*0.25;
   alpha[0] = ( x0[0] - x1[0] + x2[0] - x3[0])*0.25;
   alpha[1] = ( x0[1] - x1[1] + x2[1] - x3[1])*0.25;
   c[0]     = ( x0[0] + x1[0] + x2[0] + x3[0])*0.25;
   c[1]     = ( x0[1] + x1[1] + x2[1] + x3[1])*0.25;

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   QTEST_ZERO_JACOBIAN(jac,s0,s1,s2,s3)
}

#else  /*  if DIM == 3  */

FLOAT barycentric_coordinates(x0,x1,x2,x3,b)
FLOAT *x0, *x1, *x2, *x3, b[DIM2][DIM2];
{
   FLOAT deta, a[DIM][DIM];
  
   a[0][0] = x0[0] - x3[0];
   a[0][1] = x1[0] - x3[0];
   a[0][2] = x2[0] - x3[0];
   a[1][0] = x0[1] - x3[1];
   a[1][1] = x1[1] - x3[1];
   a[1][2] = x2[1] - x3[1];
   a[2][0] = x0[2] - x3[2];
   a[2][1] = x1[2] - x3[2];
   a[2][2] = x2[2] - x3[2];
  
   deta = a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2]) +
          a[1][0]*(a[2][1]*a[0][2]-a[0][1]*a[2][2]) +
          a[2][0]*(a[0][1]*a[1][2]-a[0][2]*a[1][1]);
  
   if ( ( (deta>0)?deta:-deta ) < 1.e-15 ){
      eprintf("ZERO DETERMINANT\n");
      return(1.);
   }
  
   b[0][0] = (a[1][1]*a[2][2]-a[2][1]*a[1][2])/deta;
   b[0][1] = (a[0][2]*a[2][1]-a[0][1]*a[2][2])/deta;
   b[0][2] = (a[0][1]*a[1][2]-a[1][1]*a[0][2])/deta;
   b[1][0] = (a[1][2]*a[2][0]-a[1][0]*a[2][2])/deta;
   b[1][1] = (a[0][0]*a[2][2]-a[2][0]*a[0][2])/deta;
   b[1][2] = (a[0][2]*a[1][0]-a[0][0]*a[1][2])/deta;   
   b[2][0] = (a[1][0]*a[2][1]-a[2][0]*a[1][1])/deta;
   b[2][1] = (a[2][0]*a[0][1]-a[0][0]*a[2][1])/deta; 
   b[2][2] = (a[0][0]*a[1][1]-a[1][0]*a[0][1])/deta;
   b[3][0] = -(b[0][0]+b[1][0]+b[2][0]);
   b[3][1] = -(b[0][1]+b[1][1]+b[2][1]);
   b[3][2] = -(b[0][2]+b[1][2]+b[2][2]);
   b[0][3] = -DOT(b[0],x3);
   b[1][3] = -DOT(b[1],x3);
   b[2][3] = -DOT(b[2],x3);
   b[3][3] = -DOT(b[3],x0);
   
   return(fabs(deta)/6.0);  /* volume of the tetrahedron */
}

#endif  /*  DIM == 3  */

#if (F_DATA & CURVED_FACE_MIDDLE) && (DIM == 2)

#define TEST_ZERO_JACOBIAN(jac,s0,s1,s2)                                       \
   s0 = jac[2];                                                                \
   s1 = jac[0]+jac[2];                                                         \
   s2 = jac[1]+jac[2];                                                         \
   if ( fabs(s0) < 1.e-7 || fabs(s1) < 1.e-7 || fabs(s2) < 1.e-7 ||            \
        s0*s1 < 1.e-14 || s0*s2 < 1.e-14 || s1*s2 < 1.e-14)                    \
      eprintf("ZERO JACOBIAN\n");


/* jac = jac[0]*x + jac[1]*y + jac[2] is the Jacobian;
   (b[i][j][0]*x + b[i][j][1]*y + b[i][j][2])/jac is the inverse of the 
   Jacobi matrix of the reference mapping */
/* The curved edge n1n2 is mapped on the longest edge of the ref. element */
void inverse_of_P2_reference_mapping0(n0,n1,n2,fa0,b,jac)
NODE *n0, *n1, *n2;
FACE *fa0;
FLOAT b[2][2][DIM2], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa0->c_midpoint->x[0] - 2.*(x1[0] + x2[0]);
   alpha[1] = 4.*fa0->c_midpoint->x[1] - 2.*(x1[1] + x2[1]);

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
  
   b[0][0][0] =  alpha[1];
   b[0][0][1] =  0.;
   b[0][0][2] =  a[1][1];
   b[0][1][0] =  -alpha[0];
   b[0][1][1] =  0.;
   b[0][1][2] =  -a[0][1];
   b[1][0][0] =  0.;
   b[1][0][1] =  -alpha[1];
   b[1][0][2] =  -a[1][0];
   b[1][1][0] =  0.;
   b[1][1][1] =  alpha[0];
   b[1][1][2] =  a[0][0];
}
 
/* The curved edge n0n2 is mapped on the edge between (0,0) and (0,1). */
void inverse_of_P2_reference_mapping1(n0,n1,n2,fa1,b,jac)
NODE *n0, *n1, *n2;
FACE *fa1;
FLOAT b[2][2][DIM2], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa1->c_midpoint->x[0] - 2.*(x0[0] + x2[0]);
   alpha[1] = 4.*fa1->c_midpoint->x[1] - 2.*(x0[1] + x2[1]);

   jac[0] = -a[0][0]*alpha[1]+a[1][0]*alpha[0];
   jac[1] = 2.*(a[1][0]*alpha[0] - a[0][0]*alpha[1]) -
            alpha[0]*a[1][1]+alpha[1]*a[0][1];
   jac[2] = a[0][0]*a[1][1]+a[0][0]*alpha[1]-a[1][0]*alpha[0]-a[1][0]*a[0][1];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
  
   b[0][0][0] =  -alpha[1];
   b[0][0][1] =  -2.*alpha[1];
   b[0][0][2] =  a[1][1]+alpha[1];
   b[0][1][0] =  alpha[0];
   b[0][1][1] =  2.*alpha[0];
   b[0][1][2] =  -a[0][1]-alpha[0];
   b[1][0][0] =  0.;
   b[1][0][1] =  alpha[1];
   b[1][0][2] =  -a[1][0];
   b[1][1][0] =  0.;
   b[1][1][1] =  -alpha[0];
   b[1][1][2] =  a[0][0];
}
 
/* The curved edge n0n1 is mapped on the edge between (0,0) and (1,0). */
void inverse_of_P2_reference_mapping2(n0,n1,n2,fa2,b,jac)
NODE *n0, *n1, *n2;
FACE *fa2;
FLOAT b[2][2][DIM2], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa2->c_midpoint->x[0] - 2.*(x0[0] + x1[0]);
   alpha[1] = 4.*fa2->c_midpoint->x[1] - 2.*(x0[1] + x1[1]);

   jac[0] = 2.*(a[0][1]*alpha[1]-alpha[0]*a[1][1])+
            alpha[0]*a[1][0]-a[0][0]*alpha[1];
   jac[1] = a[0][1]*alpha[1]-alpha[0]*a[1][1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0]+alpha[0]*a[1][1]-a[0][1]*alpha[1];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
  
   b[0][0][0] =  -alpha[1];
   b[0][0][1] =  0.;
   b[0][0][2] =  a[1][1];
   b[0][1][0] =  alpha[0];
   b[0][1][1] =  0.;
   b[0][1][2] =  -a[0][1];
   b[1][0][0] =  2.*alpha[1];
   b[1][0][1] =  alpha[1];
   b[1][0][2] =  -a[1][0]-alpha[1];
   b[1][1][0] =  -2.*alpha[0];
   b[1][1][1] =  -alpha[0];
   b[1][1][2] =  a[0][0]+alpha[0];
}

/* the reference mapping always maps (0,0)->n0, (1,0)->n1, (0,1)->n2  */
void inverse_of_P2_reference_mapping(n0,n1,n2,fa0,fa1,fa2,b,jac)
NODE *n0, *n1, *n2; 
FACE *fa0, *fa1, *fa2;
FLOAT b[2][2][DIM2], jac[DIM2];
{
   if (fa0->c_midpoint)
      inverse_of_P2_reference_mapping0(n0,n1,n2,fa0,b,jac);
   else if (fa1->c_midpoint)
      inverse_of_P2_reference_mapping1(n0,n1,n2,fa1,b,jac);
   else if (fa2->c_midpoint)
      inverse_of_P2_reference_mapping2(n0,n1,n2,fa2,b,jac);
   else
      eprintf("Error: element has no curved edge.\n");
/*
   if (fa0->c_midpoint)
      inverse_of_P2_reference_mapping0(n0,n1,n2,fa0,b,jac);
   else{ 
      if (fa1->c_midpoint){
         EXCHANGE(n1,n2,n)
         EXCHANGE(fa1,fa2,f)
      }
      inverse_of_P2_reference_mapping2(n0,n1,n2,fa2,b,jac);
   }
*/
}

/* The curved edge n1n2 is mapped on the longest edge of the ref. element */
/* If n1n2 is not curved, then all edges of the element are straight.     */
/* The reference mapping is 
                     a[0][0]*x[0] + a[0][1]*x[1] + c[0] + alpha[0]*x[0]*x[1]
                     a[1][0]*x[0] + a[1][1]*x[1] + c[1] + alpha[1]*x[0]*x[1]  */
void P2_reference_mapping0(n0,n1,n2,fa0,a,c,alpha,jac)
NODE *n0, *n1, *n2;
FACE *fa0;
FLOAT a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, s0, s1, s2;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   if (fa0->c_midpoint){
      alpha[0] = 4.*fa0->c_midpoint->x[0] - 2.*(x1[0] + x2[0]);
      alpha[1] = 4.*fa0->c_midpoint->x[1] - 2.*(x1[1] + x2[1]);
   }
   else
      alpha[0] = alpha[1] = 0.;      
   c[0] = x0[0];
   c[1] = x0[1];

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
}

void P2_reference_mapping0_with_inverse(n0,n1,n2,fa0,b,a,c,alpha,jac)
NODE *n0, *n1, *n2;
FACE *fa0;
FLOAT b[2][2][DIM2], a[DIM][DIM], c[DIM], alpha[DIM], jac[DIM2];
{
   FLOAT *x0, *x1, *x2, s0, s1, s2;
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   if (fa0->c_midpoint){
      alpha[0] = 4.*fa0->c_midpoint->x[0] - 2.*(x1[0] + x2[0]);
      alpha[1] = 4.*fa0->c_midpoint->x[1] - 2.*(x1[1] + x2[1]);
   }
   else
      alpha[0] = alpha[1] = 0.;      
   c[0] = x0[0];
   c[1] = x0[1];

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
  
   b[0][0][0] =  alpha[1];
   b[0][0][1] =  0.;
   b[0][0][2] =  a[1][1];
   b[0][1][0] =  -alpha[0];
   b[0][1][1] =  0.;
   b[0][1][2] =  -a[0][1];
   b[1][0][0] =  0.;
   b[1][0][1] =  -alpha[1];
   b[1][0][2] =  -a[1][0];
   b[1][1][0] =  0.;
   b[1][1][1] =  alpha[0];
   b[1][1][2] =  a[0][0];
}
 
/* The curved edge n1n2 is mapped on the longest edge of the ref. element */
void jacobian0(n0,n1,n2,fa0,jac)
NODE *n0, *n1, *n2;
FACE *fa0;
FLOAT jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa0->c_midpoint->x[0] - 2.*(x1[0] + x2[0]);
   alpha[1] = 4.*fa0->c_midpoint->x[1] - 2.*(x1[1] + x2[1]);

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
}
 
/* The curved edge n0n2 is mapped on the edge between (0,0) and (0,1). */
void jacobian1(n0,n1,n2,fa1,jac)
NODE *n0, *n1, *n2;
FACE *fa1;
FLOAT jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa1->c_midpoint->x[0] - 2.*(x0[0] + x2[0]);
   alpha[1] = 4.*fa1->c_midpoint->x[1] - 2.*(x0[1] + x2[1]);

   jac[0] = -a[0][0]*alpha[1]+a[1][0]*alpha[0];
   jac[1] = 2.*(a[1][0]*alpha[0] - a[0][0]*alpha[1]) -
            alpha[0]*a[1][1]+alpha[1]*a[0][1];
   jac[2] = a[0][0]*a[1][1]+a[0][0]*alpha[1]-a[1][0]*alpha[0]-a[1][0]*a[0][1];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
}
 
/* The curved edge n0n1 is mapped on the edge between (0,0) and (1,0). */
void jacobian2(n0,n1,n2,fa2,jac)
NODE *n0, *n1, *n2;
FACE *fa2;
FLOAT jac[DIM2];
{
   FLOAT *x0, *x1, *x2, a[2][2], s0, s1, s2, alpha[2];
  
   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   x2 = n2->myvertex->x;
   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];
   alpha[0] = 4.*fa2->c_midpoint->x[0] - 2.*(x0[0] + x1[0]);
   alpha[1] = 4.*fa2->c_midpoint->x[1] - 2.*(x0[1] + x1[1]);

   jac[0] = 2.*(a[0][1]*alpha[1]-alpha[0]*a[1][1])+
            alpha[0]*a[1][0]-a[0][0]*alpha[1];
   jac[1] = a[0][1]*alpha[1]-alpha[0]*a[1][1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0]+alpha[0]*a[1][1]-a[0][1]*alpha[1];
   TEST_ZERO_JACOBIAN(jac,s0,s1,s2)
}

/* the reference mapping always maps (0,0)->n0, (1,0)->n1, (0,1)->n2  */
void jacobian_c(n0,n1,n2,fa0,fa1,fa2,jac)
NODE *n0, *n1, *n2; 
FACE *fa0, *fa1, *fa2;
FLOAT jac[DIM2];
{
   if (fa0->c_midpoint)
      jacobian0(n0,n1,n2,fa0,jac);
   else if (fa1->c_midpoint)
      jacobian1(n0,n1,n2,fa1,jac);
   else if (fa2->c_midpoint)
      jacobian2(n0,n1,n2,fa2,jac);
   else
      eprintf("Error: element has no curved edge.\n");
}

#else

void inverse_of_P2_reference_mapping(n0,n1,n2,fa0,fa1,fa2,b,jac)
NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; FLOAT b[2][2][DIM2], jac[DIM2];
{  eprintf("Error: inverse_of_P2_reference_mapping not available.\n");  }

void jacobian_c(n0,n1,n2,fa0,fa1,fa2,jac)
NODE *n0, *n1, *n2; FACE *fa0, *fa1, *fa2; FLOAT jac[DIM2];
{  eprintf("Error: jacobian_c not available.\n");  }

#endif

#if DIM == 2 && ELEMENT_TYPE == SIMPLEX

void reference_mapping(pel,type,ref_map)
ELEMENT *pel;
INT type;
REF_MAPPING *ref_map;
{
   NODE *n0, *n1, *n2;

   ref_map->type = type;
   ref_map->domain = SIMPLEX;
   NODES_OF_ELEMENT(n0,n1,n2,pel);
   if (type == P1_REF_MAP){
      ref_map->area = P1_reference_mapping0(n0,n1,n2,
                                                   ref_map->p1_a,ref_map->p1_c);
      ref_map->p1_jac = ref_map->p1_a[0][0]*ref_map->p1_a[1][1]-
                        ref_map->p1_a[0][1]*ref_map->p1_a[1][0];
   }
   else
      eprintf("Error: reference_mapping not available.\n");
}

void reference_mapping_with_inverse(pel,type,ref_map)
ELEMENT *pel;
INT type;
REF_MAPPING *ref_map;
{
   NODE *n0, *n1, *n2;

   ref_map->type = type;
   ref_map->domain = SIMPLEX;
   NODES_OF_ELEMENT(n0,n1,n2,pel);
   if (type == P1_REF_MAP){
      ref_map->area = P1_reference_mapping0_with_inverse(n0,n1,n2,
                                     ref_map->p1_b,ref_map->p1_a,ref_map->p1_c);
      ref_map->p1_jac = ref_map->p1_a[0][0]*ref_map->p1_a[1][1]-
                        ref_map->p1_a[0][1]*ref_map->p1_a[1][0];
   }
   else
      eprintf("Error: reference_mapping_with_inverse not available.\n");
}

#elif DIM == 2 && ELEMENT_TYPE == CUBE

void reference_mapping(pel,type,ref_map)
ELEMENT *pel;
INT type;
REF_MAPPING *ref_map;
{
   NODE *n0, *n1, *n2, *n3;

   ref_map->type = type;
   ref_map->domain = CUBE;
   NODES_OF_4ELEMENT(n0,n1,n2,n3,pel);
   if (type == Q1_REF_MAP)
      Q1_reference_mapping(n0,n1,n2,n3,
                 ref_map->q1_a,ref_map->q1_alpha,ref_map->q1_c,ref_map->q1_jac);
   else
      eprintf("Error: reference_mapping not available.\n");
}

void reference_mapping_with_inverse(pel,type,ref_map)
ELEMENT *pel;
INT type;
REF_MAPPING *ref_map;
{
   NODE *n0, *n1, *n2, *n3;

   ref_map->type = type;
   ref_map->domain = CUBE;
   NODES_OF_4ELEMENT(n0,n1,n2,n3,pel);
   if (type == Q1_REF_MAP)
      Q1_reference_mapping_with_inverse(n0,n1,n2,n3,ref_map->q1_b,
                 ref_map->q1_a,ref_map->q1_c,ref_map->q1_alpha,ref_map->q1_jac);
   else
      eprintf("Error: reference_mapping_with_inverse not available.\n");
}

#else

void reference_mapping(pel,type,ref_map)
ELEMENT *pel; INT type; REF_MAPPING *ref_map;
{  eprintf("Error: reference_mapping not available.\n");  }

void reference_mapping_with_inverse(pel,type,ref_map)
ELEMENT *pel; INT type; REF_MAPPING *ref_map;
{  eprintf("Error: reference_mapping_with_inverse not available.\n");  }

#endif

FLOAT jacobian(x,ref_map)
FLOAT *x;
REF_MAPPING *ref_map;
{
   if (ref_map->type == Q1_REF_MAP)
      return(LINV(ref_map->q1_jac,x));
   else if (ref_map->type == P1_REF_MAP)
      return(ref_map->p1_jac);
   else{
      eprintf("Error: jacobian not available.\n");
      return(0.);
   }
}

void inv_of_jac_matr_times_jacobian(x,b,ref_map)
FLOAT *x, b[2][2];
REF_MAPPING *ref_map;
{
   if (ref_map->type == Q1_REF_MAP){
      b[0][0] = LINV(ref_map->q1_b[0][0],x);
      b[0][1] = LINV(ref_map->q1_b[0][1],x);
      b[1][0] = LINV(ref_map->q1_b[1][0],x);
      b[1][1] = LINV(ref_map->q1_b[1][1],x);
   }
   else if (ref_map->type == P1_REF_MAP){
      b[0][0] =  ref_map->p1_a[1][1];
      b[0][1] = -ref_map->p1_a[0][1];
      b[1][0] = -ref_map->p1_a[1][0];
      b[1][1] =  ref_map->p1_a[0][0];
   }
   else
      eprintf("Error: inv_of_jac_matr_times_jacobian not available.\n");
}

#if N_DATA & NUMBER_OF_N_NEL

#if (DIM == 2 && ELEMENT_TYPE == SIMPLEX)

void compute_number_of_neighbouring_elements_for_nodes(theGrid)
GRID *theGrid;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ){
      (pel->n[0]->nel)++;
      (pel->n[1]->nel)++;
      (pel->n[2]->nel)++;
   }
}

#elif (DIM == 3 && ELEMENT_TYPE == SIMPLEX) || (DIM == 2 && ELEMENT_TYPE == CUBE)

void compute_number_of_neighbouring_elements_for_nodes(theGrid)
GRID *theGrid;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ){
      (pel->n[0]->nel)++;
      (pel->n[1]->nel)++;
      (pel->n[2]->nel)++;
      (pel->n[3]->nel)++;
   }
}

#endif

#else

void compute_number_of_neighbouring_elements_for_nodes(theGrid)
GRID *theGrid;
{}

#endif

#if (DATA_STR & LG_DATA) && (DIM == 3)

void lg_neighbours(tGrid,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata,
                                                                          x,y,z)
GRID *tGrid;
NLGLINK **pnlglink;
LGNLINK **plgnlink;
LGFLINK **plgflink;
LGLGLINK **plglglink;
FLGLINK **pflglink;
LGDATA **plgdata;
INT x, y, z;  /* auxiliary face variables */
{
   ELEMENT *pel;
   NODE *pnode;
   NLGLINK *pnlg;
   LGNLINK *plgn;
   LGFLINK *plgf;
   LGLGLINK *plglg;
   FLGLINK *pflg;
   FLOAT ar, d, n[4];
   INT i, j, in[3];

   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (IS_IN_GLG(pnode)){
         pnode->lgd = (*plgdata)++;
         pnode->lgd->nstart = NULL;
         pnode->lgd->fstart = NULL;
         pnode->lgd->lgstart = NULL;
         pnode->lgd->n[0] = pnode->lgd->n[1] = pnode->lgd->n[2] = 0.0;
      }
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i < 4; i++){
         if (pel->n[i]->lgd)
            for (j = 0; j < 4; j++){
               if (j != i){
                  if (IS_FN(pel->n[j])){
                     if (pel->n[i]->lgd->nstart == NULL){
                        pel->n[i]->lgd->nstart = *plgnlink;
                        (*plgnlink)->nbnode = pel->n[j];
                        ((*plgnlink)++)->next = NULL;
                     }
                     else{
                        for (plgn = pel->n[i]->lgd->nstart; plgn->next != NULL 
                               && plgn->nbnode != pel->n[j]; plgn = plgn->next);
                        if (plgn->nbnode != pel->n[j]){
                           plgn->next = *plgnlink;
                           (*plgnlink)->nbnode = pel->n[j];
                           ((*plgnlink)++)->next = NULL;
                        }
                     }
                     if (pel->n[j]->lgstart == NULL){
                        pel->n[j]->lgstart = *pnlglink;
                        (*pnlglink)->nbnode = pel->n[i];
                        ((*pnlglink)++)->next = NULL;
                     }
                     else{
                        for (pnlg = pel->n[j]->lgstart; pnlg->next != NULL && 
                                  pnlg->nbnode != pel->n[i]; pnlg = pnlg->next);
                        if (pnlg->nbnode != pel->n[i]){
                           pnlg->next = *pnlglink;
                           (*pnlglink)->nbnode = pel->n[i];
                           ((*pnlglink)++)->next = NULL;
                        }
                     }
                  }
                  else if (pel->n[j]->lgd)
                     if (pel->n[i]->lgd->lgstart == NULL){
                        pel->n[i]->lgd->lgstart = *plglglink;
                        (*plglglink)->nbnode = pel->n[j];
                        ((*plglglink)++)->next = NULL;
                     }
                     else{
                        for (plglg = pel->n[i]->lgd->lgstart; plglg->next !=NULL
                            && plglg->nbnode != pel->n[j]; plglg = plglg->next);
                        if (plglg->nbnode != pel->n[j]){
                           plglg->next = *plglglink;
                           (*plglglink)->nbnode = pel->n[j];
                           ((*plglglink)++)->next = NULL;
                        }
                     }
               }
               if (IS_FF(pel->f[j])){
                  if (pel->n[i]->lgd->fstart == NULL){
                     pel->n[i]->lgd->fstart = *plgflink;
                     (*plgflink)->nbface = pel->f[j];
                     ((*plgflink)++)->next = NULL;
                  }
                  else{
                     for (plgf = pel->n[i]->lgd->fstart; plgf->next != NULL && 
                                  plgf->nbface != pel->f[j]; plgf = plgf->next);
                     if (plgf->nbface != pel->f[j]){
                        plgf->next = *plgflink;
                        (*plgflink)->nbface = pel->f[j];
                        ((*plgflink)++)->next = NULL;
                     }
                  }
                  if (pel->f[j]->lgstart == NULL){
                     pel->f[j]->lgstart = *pflglink;
                     (*pflglink)->nbnode = pel->n[i];
                     ((*pflglink)++)->next = NULL;
                  }
                  else{
                     for (pflg = pel->f[j]->lgstart; pflg->next != NULL && 
                                  pflg->nbnode != pel->n[i]; pflg = pflg->next);
                     if (pflg->nbnode != pel->n[i]){
                        pflg->next = *pflglink;
                        (*pflglink)->nbnode = pel->n[i];
                        ((*pflglink)++)->next = NULL;
                     }
                  }
               }
            }
         if (NOT_FF(pel->f[i]) && 
             (pel->n[0]->lgd && i != 0 || pel->n[1]->lgd && i != 1 || 
              pel->n[2]->lgd && i != 2 || pel->n[3]->lgd && i != 3)){
            j = 0;
            if (i != 0) in[j++] = 0;
            if (i != 1) in[j++] = 1;
            if (i != 2) in[j++] = 2;
            if (i != 3) in[j++] = 3;
            ar = normal_vector(MYVERTEX(pel->n[in[0]])->x,
                               MYVERTEX(pel->n[in[1]])->x,
                               MYVERTEX(pel->n[in[2]])->x,pel->f[i],n);
            FD(pel->f[i],x) = n[0]*ar;
            FD(pel->f[i],y) = n[1]*ar;
            FD(pel->f[i],z) = n[2]*ar;
         }
      }        
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      for (i = 0; i < 4; i++)
         if (NOT_FF(pel->f[i]))
            for (j = 0; j < 4; j++)
               if (pel->n[j]->lgd && j != i){
                  pel->n[j]->lgd->n[0] += FD(pel->f[i],x);
                  pel->n[j]->lgd->n[1] += FD(pel->f[i],y);
                  pel->n[j]->lgd->n[2] += FD(pel->f[i],z);
               }   
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (pnode->lgd){
         d = sqrt(DOT(pnode->lgd->n,pnode->lgd->n));
         n[1] = pnode->lgd->n[0] /= d;
         n[2] = pnode->lgd->n[1] /= d;
         n[3] = pnode->lgd->n[2] /= d;
         d = sqrt(n[2]*n[2] + n[3]*n[3]);
         pnode->lgd->t1[0] = 0.;
         pnode->lgd->t1[1] = n[3]/d;
         pnode->lgd->t1[2] = -n[2]/d;
         d = sqrt((n[1]*n[1]-1.)*(n[1]*n[1]-1.) + n[1]*n[1]*n[2]*n[2]
                                                + n[1]*n[1]*n[3]*n[3]);
         pnode->lgd->t2[0] = (n[1]*n[1]-1.)/d;
         pnode->lgd->t2[1] = n[1]*n[2]/d;
         pnode->lgd->t2[2] = n[1]*n[3]/d;
      }
}

#else

void lg_neighbours(tGrid,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata,
                                                                          x,y,z)
GRID *tGrid;
NLGLINK **pnlglink;
LGNLINK **plgnlink;
LGFLINK **plgflink;
LGLGLINK **plglglink;
FLGLINK **pflglink;
LGDATA **plgdata;
INT x, y, z;  /* auxiliary face variables */
{}

#endif

#if E_DATA & E_E_NEIGHBOURS

void new_elem_neighbour(pelem,pel,pelink)
ELEMENT *pelem, *pel;
ELINK **pelink;
{
   ELINK *peli;
   
   if (pelem->estart == NULL)
      pelem->estart = *pelink;
   else{
      for (peli = pelem->estart; peli->next != NULL; peli = peli->next);
      peli->next = *pelink;
   }
   (*pelink)->nbel = pel;
   ((*pelink)++)->next = NULL;
}

void element_neighbours_on_first_grid(tGrid,pelink)
GRID *tGrid;
ELINK **pelink;
{
   ELEMENT *pelem, *pel;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ)
      pelem->estart = NULL;
   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ)
      for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
         if (pelem != pel && ARE_NEIGHBOURS(pelem,pel))
            new_elem_neighbour(pelem,pel,pelink);
}

#else

void element_neighbours_on_first_grid(tGrid,pelink)
GRID *tGrid;
ELINK **pelink;
{}

#endif

#if E_DATA & E_E_FNEIGHBOURS

void add_elem_fneighbour(pelem,pel)
ELEMENT *pelem, *pel;
{
   INT i, j=-1;
   
   for (i = 0; i < SIDES; i++)
      if (FACE_CONTAINED(pelem->f[i],pel)) j = i;
   if (j > -1)
      NB_EL(pelem,j) = pel;
}
               
void element_fneighbours_on_first_grid(tGrid)
GRID *tGrid;
{
   ELEMENT *pelem, *pel;
   INT i;
   
   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
         if (pelem != pel)
            add_elem_fneighbour(pelem,pel);
      for (i = 0; i < SIDES; i++)
         if (IS_BF(pelem->f[i]))
            NB_EL(pelem,i) = NULL;
   }
}
               
#else

void element_fneighbours_on_first_grid(tGrid)
GRID *tGrid;
{}

#endif

#if DATA_S & PREVIOUS_NODE

void make_previous_nodes(tGrid)
GRID *tGrid;
{
   NODE *theNode;

   PREV(FIRSTN(tGrid)) = NULL;
   for (theNode=FIRSTN(tGrid); SUCC(theNode); theNode=SUCC(theNode))
      PREV(SUCC(theNode)) = theNode;
}

#else 

void make_previous_nodes(tGrid)
GRID *tGrid;
{}

#endif

#if DATA_S & PREVIOUS_FACE

void make_previous_faces(tGrid)
GRID *tGrid;
{
   FACE *theFace;

   PREV(FIRSTF(tGrid)) = NULL;
   for (theFace=FIRSTF(tGrid); SUCC(theFace); theFace=SUCC(theFace))
      PREV(SUCC(theFace)) = theFace;
}

#else

void make_previous_faces(tGrid)
GRID *tGrid;
{}

#endif

#if (DATA_S & SPECIAL_NODES_AND_FACES) && (DIM == 2)

void fill_snode(sn,n,f,pel)
SNODE *sn;
NODE *n;
FACE *f;
ELEMENT *pel;
{
   if (sn->n1 && sn->n1 != n){
      sn->n2 = n;
      sn->f2 = f;
      sn->pel2 = pel;
   }
   else{
      sn->n1 = n;
      sn->f1 = f;
      sn->pel1 = pel;
   }
}

void make_special_nodes_and_faces(pgrid,psnode,psface,nmask,fmask)
GRID *pgrid; 
SNODE **psnode;
SFACE **psface;
INT nmask, fmask;
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   SNODE *lsnode, *fsnode, *psn;
   SFACE *lsface, *fsface;

   lsnode = NULL;
   fsnode = *psnode;
   for (pnode = FIRSTN(pgrid); pnode; pnode = pnode->succ)
      if (NTYPE(pnode) & nmask){
         pnode->s_node = *psnode;         
         (*psnode)->n0 = pnode;
         (*psnode)->n1 = (*psnode)->n2 = NULL;
         if (lsnode)
            lsnode = lsnode->succ = (*psnode)++;
         else
            lsnode = (*psnode)++;
      }
   if (lsnode){
      lsnode->succ = NULL;
      FIRSTSN(pgrid) = fsnode;
   }
   else
      FIRSTSN(pgrid) = NULL;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ){
      if (   (NTYPE(pel->n[0]) & nmask) && (FTYPE(pel->f[1]) & fmask))
         fill_snode(pel->n[0]->s_node,  pel->n[2],pel->f[1],pel);
      if (   (NTYPE(pel->n[0]) & nmask) && (FTYPE(pel->f[2]) & fmask))
         fill_snode(pel->n[0]->s_node,  pel->n[1],pel->f[2],pel);
      if (   (NTYPE(pel->n[1]) & nmask) && (FTYPE(pel->f[0]) & fmask))
         fill_snode(pel->n[1]->s_node,  pel->n[2],pel->f[0],pel);
      if (   (NTYPE(pel->n[1]) & nmask) && (FTYPE(pel->f[2]) & fmask))
         fill_snode(pel->n[1]->s_node,  pel->n[0],pel->f[2],pel);
      if (   (NTYPE(pel->n[2]) & nmask) && (FTYPE(pel->f[0]) & fmask))
         fill_snode(pel->n[2]->s_node,  pel->n[1],pel->f[0],pel);
      if (   (NTYPE(pel->n[2]) & nmask) && (FTYPE(pel->f[1]) & fmask))
         fill_snode(pel->n[2]->s_node,  pel->n[0],pel->f[1],pel);
   }
   for (psn = FIRSTSN(pgrid); psn; psn = psn->succ)
      if (psn->n1 == NULL)
         eprintf("Error in make_special_nodes_and_faces.\n");
      else if (psn->n2 == NULL){
         psn->n2 = psn->n1;
         psn->f2 = psn->f1;
         psn->pel2 = NULL;
      }

   lsface = NULL;
   fsface = *psface;
   for (pface = FIRSTF(pgrid); pface; pface = pface->succ)
      if (FTYPE(pface) & fmask){
         pface->s_face = *psface;         
         (*psface)->f = pface;
         if (lsface)
            lsface = lsface->succ = (*psface)++;
         else
            lsface = (*psface)++;
      }
   if (lsface){
      lsface->succ = NULL;
      FIRSTSF(pgrid) = fsface;
   }
   else
      FIRSTSF(pgrid) = NULL;
   for (pel = FIRSTELEMENT(pgrid); pel; pel = pel->succ){
      if (FTYPE(pel->f[0]) & fmask){
         pel->f[0]->s_face->n1  = pel->n[1];
         pel->f[0]->s_face->n2  = pel->n[2];
         pel->f[0]->s_face->pel = pel;
      }
      if (FTYPE(pel->f[1]) & fmask){
         pel->f[1]->s_face->n1  = pel->n[0];
         pel->f[1]->s_face->n2  = pel->n[2];
         pel->f[1]->s_face->pel = pel;
      }
      if (FTYPE(pel->f[2]) & fmask){
         pel->f[2]->s_face->n1  = pel->n[0];
         pel->f[2]->s_face->n2  = pel->n[1];
         pel->f[2]->s_face->pel = pel;
      }
   }
}

#else

void make_special_nodes_and_faces(pgrid,psnode,psface,nmask,fmask)
GRID *pgrid; SNODE **psnode; SFACE **psface; INT nmask, fmask;
{
   FIRSTSN(pgrid) = NULL;
   FIRSTSF(pgrid) = NULL;
}

#endif

#if (DATA_S & PERIODIC_NODES) && (DIM == 2)

void make_periodic_nodes(pgrid,pp_node,nmask,k)
GRID *pgrid; 
P_NODE **pp_node;
INT nmask, k;
{
   NODE *pnode, *pn;
   INT i=1;

   FIRSTPN(pgrid) = *pp_node;
   for (pnode = FIRSTNODE(pgrid); pnode; pnode = pnode->succ)
      if (NTYPE(pnode) & nmask){
         for (pn = FIRSTNODE(pgrid); pn && (NOT_BN(pn) || pn == pnode ||
                     fabs(pn->myvertex->x[k] - pnode->myvertex->x[k]) > EPS); 
                                                                pn = pn->succ);
         if (pn == NULL)
            eprintf("Error in make_periodic_nodes.\n");
         else{
            (*pp_node)->n1 = pnode;
            (*pp_node)->n2 = pn;
            (*pp_node)->succ = (*pp_node) + 1;
            (*pp_node)++;
            i = 0;
         }
      }
   if (i){
      eprintf("Error: no nodes for periodic b.c.\n");
      FIRSTPN(pgrid) = NULL;
   }
   else
      (*pp_node-1)->succ = NULL;
}

#else

void make_periodic_nodes(pgrid,pp_node,nmask,k)
GRID *pgrid; P_NODE **pp_node; INT nmask, k;
{
   FIRSTPN(pgrid) = NULL;
}

#endif

FLOAT edge_length_square();

FLOAT update_max_and_min_length(n1,n2,max,min,sum,m)
NODE *n1, *n2;
FLOAT *max, *min, *sum;
INT *m;
{
   FLOAT x;
  
   if (*max < (x=sqrt(edge_length_square(n1,n2)))) *max = x;
   else if (*min > x) *min = x;
   *sum += x;
   if (IS_BN(n1) && IS_BN(n2)){
      *sum += x;
      (*m)++;  
   }
   return(x);
}

#if DIM == 3

FLOAT edge_length_square(n1,n2)
NODE *n1, *n2;
{
   FLOAT x, y, z;
   
   x = MYVERTEX(n1)->x[0] - MYVERTEX(n2)->x[0];
   y = MYVERTEX(n1)->x[1] - MYVERTEX(n2)->x[1];
   z = MYVERTEX(n1)->x[2] - MYVERTEX(n2)->x[2];
   return(x*x+y*y+z*z);
}

#if ELEMENT_TYPE == SIMPLEX

FLOAT volume(x0,x1,x2,x3)
FLOAT *x0, *x1, *x2, *x3;
{
   FLOAT deta, a[4][4];
  
   a[1][1] = x1[0] - x0[0];
   a[1][2] = x2[0] - x0[0];
   a[1][3] = x3[0] - x0[0];
   a[2][1] = x1[1] - x0[1];
   a[2][2] = x2[1] - x0[1];
   a[2][3] = x3[1] - x0[1];
   a[3][1] = x1[2] - x0[2];
   a[3][2] = x2[2] - x0[2];
   a[3][3] = x3[2] - x0[2];
  
   deta = a[1][1]*(a[2][2]*a[3][3]-a[3][2]*a[2][3])+
          a[2][1]*(a[3][2]*a[1][3]-a[1][2]*a[3][3])+
          a[3][1]*(a[1][2]*a[2][3]-a[1][3]*a[2][2]);
  
   if ( ( (deta>0)?deta:-deta ) < 1.e-15 ) {
      eprintf("ZERO DETERMINANT: %e\n",deta);
      return(-1.);
   }
   
   return(fabs(deta)/6.0);  /* volume of the tetrahedron (x0,x1,x2,x3) */
}
 
FLOAT max_edge_length(pel)
ELEMENT *pel;
{
   FLOAT m=0., x;
   
   if (m < (x=edge_length_square(pel->n[0],pel->n[1]))) m = x;
   if (m < (x=edge_length_square(pel->n[0],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[0],pel->n[3]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[3]))) m = x;
   if (m < (x=edge_length_square(pel->n[2],pel->n[3]))) m = x;
   return(sqrt(m));
}

void max_and_min_edge_length(pel,max,min,sum,m,sigma)
ELEMENT *pel;
FLOAT *max, *min, *sum, *sigma;
INT *m;
{
   update_max_and_min_length(pel->n[0],pel->n[1],max,min,sum,m);
   update_max_and_min_length(pel->n[0],pel->n[2],max,min,sum,m);
   update_max_and_min_length(pel->n[0],pel->n[3],max,min,sum,m);
   update_max_and_min_length(pel->n[1],pel->n[2],max,min,sum,m);
   update_max_and_min_length(pel->n[1],pel->n[3],max,min,sum,m);
   update_max_and_min_length(pel->n[2],pel->n[3],max,min,sum,m);
   *sigma = -1.;
}

FLOAT ratio(pel,lmax)
ELEMENT *pel;
FLOAT *lmax;
{
   FLOAT maxvol, m;
    
   m = max_edge_length(pel);
   if (m > *lmax) *lmax = m;
   maxvol = 0.11785113*m*m*m;
   return(volume(MYVERTEX(pel->n[0])->x,MYVERTEX(pel->n[1])->x,
                        MYVERTEX(pel->n[2])->x,MYVERTEX(pel->n[3])->x)/maxvol);
}

void check_ratio(pgrid)
GRID *pgrid;
{
   ELEMENT *pel;
   FLOAT x, min, max, lmax;

   lmax = 0.;
   pel = FIRSTELEMENT(pgrid);
   min = max = ratio(pel,&lmax);
   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ)
      if (max < (x=ratio(pel,&lmax)))
         max = x;
      else if (x < min)
         min = x;
   printf("Ratio: min = %f, max = %f\n",min,max);
   printf("Max. edge length: %f\n",lmax);
}

FLOAT diameter(pel)
ELEMENT *pel;
{
   return(max_edge_length(pel));
}

FLOAT sdiameter(pel,b,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel; FLOAT b[DIM2][DIM2], s0, s1;
{  eprintf("Error: sdiameter not available.\n"); return(0.);  }

FLOAT gsdiameter(pel,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel; FLOAT s0, s1;
{  eprintf("Error: gsdiameter not available.\n"); return(0.);  }

#endif

#else  /*  DIM == 2  */

FLOAT edge_length_square(n1,n2)
NODE *n1, *n2;
{
   FLOAT x, y;
   
   x = MYVERTEX(n1)->x[0] - MYVERTEX(n2)->x[0];
   y = MYVERTEX(n1)->x[1] - MYVERTEX(n2)->x[1];
   return(x*x+y*y);
}

#if ELEMENT_TYPE == SIMPLEX

FLOAT volume(x1,x2,x3)
FLOAT *x1, *x2, *x3;
{
   FLOAT a[DIM], b[DIM];

   SUBTR(x2,x1,a);
   SUBTR(x3,x1,b);
   return(fabs(a[0]*b[1]-a[1]*b[0])/2.0);
}

FLOAT max_edge_length(pel)
ELEMENT *pel;
{
   FLOAT m=0., x;
   
   if (m < (x=edge_length_square(pel->n[0],pel->n[1]))) m = x;
   if (m < (x=edge_length_square(pel->n[0],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[2]))) m = x;
   return(sqrt(m));
}

DOUBLE sigma_squared(n0,n1,n2)
NODE *n0, *n1, *n2;
{
   FLOAT h0, h1, h2, s;

   h2 = sqrt(edge_length_square(n0,n1));
   h1 = sqrt(edge_length_square(n0,n2));
   h0 = sqrt(edge_length_square(n1,n2));
   s = (h0 + h1 + h2)/2.;
   s = 4.*(s-h0)*(s-h1)*(s-h2)/s;
   if (h0 > h2) h2 = h0;
   if (h1 > h2) h2 = h1;
   return(h2*h2/s);
}

void max_and_min_edge_length(pel,max,min,sum,m,sigma)
ELEMENT *pel;
FLOAT *max, *min, *sum, *sigma;
INT *m;
{
   FLOAT h0, h1, h2, s;

   h2 = update_max_and_min_length(pel->n[0],pel->n[1],max,min,sum,m);
   h1 = update_max_and_min_length(pel->n[0],pel->n[2],max,min,sum,m);
   h0 = update_max_and_min_length(pel->n[1],pel->n[2],max,min,sum,m);
   s = (h0 + h1 + h2)/2.;
   s = 4.*(s-h0)*(s-h1)*(s-h2)/s;
   if (h0 > h2) h2 = h0;
   if (h1 > h2) h2 = h1;
   s = h2*h2/s;
   if (s > *sigma) *sigma = s;
}

FLOAT diameter(pel)
ELEMENT *pel;
{
   return(max_edge_length(pel));
}

FLOAT sdiameter(pel,b,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel;
FLOAT b[DIM2][DIM2], s0, s1;
{ 
   FLOAT z=sqrt(s0*s0+s1*s1);

   if (z < 1.e-10)
      return(diameter(pel));
   else
      return(2*z/(fabs(b[0][0]*s0+b[0][1]*s1)+
                  fabs(b[1][0]*s0+b[1][1]*s1)+
                  fabs(b[2][0]*s0+b[2][1]*s1)));
}

FLOAT gsdiameter(pel,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel;
FLOAT s0, s1;
{ 
   FLOAT bar[DIM2][DIM2];

   barycentric_coordinates(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                           pel->n[2]->myvertex->x,bar);
   return(sdiameter(pel,bar,s0,s1));
}

#elif ELEMENT_TYPE == CUBE

FLOAT volume_q(x0,x1,x2,x3)  /*  volume of a quadrilateral  */
FLOAT *x0, *x1, *x2, *x3;
{
   return(fabs(
      (-x0[0] + x1[0] + x2[0] - x3[0])*(-x0[1] - x1[1] + x2[1] + x3[1]) -
      (-x0[1] + x1[1] + x2[1] - x3[1])*(-x0[0] - x1[0] + x2[0] + x3[0]) )*0.25);
}

FLOAT max_edge_length(pel)
ELEMENT *pel;
{
   FLOAT m=0., x;
   
   if (m < (x=edge_length_square(pel->n[0],pel->n[1]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[2],pel->n[3]))) m = x;
   if (m < (x=edge_length_square(pel->n[3],pel->n[0]))) m = x;
   return(sqrt(m));
}

void max_and_min_edge_length(pel,max,min,sum,m,sigma)
ELEMENT *pel;
FLOAT *max, *min, *sum, *sigma;
INT *m;
{
   FLOAT h0, h1, h2, h3, s;

   h0 = update_max_and_min_length(pel->n[0],pel->n[1],max,min,sum,m);
   h1 = update_max_and_min_length(pel->n[1],pel->n[2],max,min,sum,m);
   h2 = update_max_and_min_length(pel->n[2],pel->n[3],max,min,sum,m);
   h3 = update_max_and_min_length(pel->n[3],pel->n[0],max,min,sum,m);
/*
   s = (h0 + h1 + h2)/2.;
   s = 2.*sqrt((s-h0)*(s-h1)*(s-h2)/s);
   if (h0 > h2) h2 = h0;
   if (h1 > h2) h2 = h1;
   s = h2/s;
   if (s > *sigma) *sigma = s;
*/
   *sum = 0.;
   *sigma = -1.;
}

FLOAT diameter(pel)
ELEMENT *pel;
{
   FLOAT m=0., x;
   
   if (m < (x=edge_length_square(pel->n[0],pel->n[1]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[2],pel->n[3]))) m = x;
   if (m < (x=edge_length_square(pel->n[3],pel->n[0]))) m = x;
   if (m < (x=edge_length_square(pel->n[0],pel->n[2]))) m = x;
   if (m < (x=edge_length_square(pel->n[1],pel->n[3]))) m = x;
   return(sqrt(m));
}

FLOAT sdiameter(pel,b,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel; FLOAT b[DIM2][DIM2], s0, s1;
{  eprintf("Error: sdiameter not available.\n"); return(0.);  }

FLOAT r_q1_0_0();
FLOAT r_q1_0_1();
FLOAT r_q1_1_0();
FLOAT r_q1_1_1();
FLOAT r_q1_2_0();
FLOAT r_q1_2_1();
FLOAT r_q1_3_0();
FLOAT r_q1_3_1();

FLOAT gsdiameter(pel,s0,s1)  /* diameter of pel in the direction of (s0,s1) */
ELEMENT *pel;
FLOAT s0, s1;
{
   FLOAT x[DIM]={0.,0.}, b[DIM][DIM], z=sqrt(s0*s0+s1*s1);
   REF_MAPPING ref_map;

   reference_mapping_with_inverse(pel,REF_MAP,&ref_map);
   inv_of_jac_matr_times_jacobian(x,b,&ref_map);
   if (z < 1.e-10)
      return(diameter(pel));
   else
      return(2*z*jacobian(x,&ref_map)/
                 (fabs(s0*(r_q1_0_0(x)*b[0][0] + r_q1_0_1(x)*b[1][0]) +
                       s1*(r_q1_0_0(x)*b[0][1] + r_q1_0_1(x)*b[1][1])) +
                  fabs(s0*(r_q1_1_0(x)*b[0][0] + r_q1_1_1(x)*b[1][0]) +
                       s1*(r_q1_1_0(x)*b[0][1] + r_q1_1_1(x)*b[1][1])) +
                  fabs(s0*(r_q1_2_0(x)*b[0][0] + r_q1_2_1(x)*b[1][0]) +
                       s1*(r_q1_2_0(x)*b[0][1] + r_q1_2_1(x)*b[1][1])) +
                  fabs(s0*(r_q1_3_0(x)*b[0][0] + r_q1_3_1(x)*b[1][0]) +
                       s1*(r_q1_3_0(x)*b[0][1] + r_q1_3_1(x)*b[1][1]))));
}

#endif

FLOAT max_diameter_in_grid();

void check_ratio(pgrid)
GRID *pgrid;
{
   printf("Ratio: min = ?, max = ?\n");
   printf("Max. diameter: %f\n",max_diameter_in_grid(pgrid));
}

#endif

FLOAT max_diameter_in_grid(tGrid)
GRID *tGrid;
{
   ELEMENT *pelem;
   FLOAT hT, h=0.;

   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      hT = diameter(pelem);
      if (hT > h)
         h = hT;
   } 
   return(h);
}

void max_and_min_edge_in_grid(tGrid)
GRID *tGrid;
{
   ELEMENT *pelem;
   FLOAT max=0., min=1.e40, aver=0., sigma=0.;
   INT m=0, n=0;

   for (pelem = FIRSTELEMENT(tGrid); pelem != NULL; pelem = pelem->succ){
      max_and_min_edge_length(pelem,&max,&min,&aver,&m,&sigma);
      n++;
   } 
   aver /= n*3*(DIM-1) + m;
   if (sigma > 0) sigma = sqrt(sigma);
   printf("Edge lengths: max = %e, min = %e, average = %e;   sigma = %e\n",
           max,min,aver,sigma);
   if (tGrid->level == 1){
      if (sigma > MAX_SIGMA)
         MAX_SIGMA = 1.1*sigma;
   }
   else
      if (sigma > MAX_SIGMA && REF_TYPE != RED_GREEN)
          eprintf("Error: sigma > MAX_SIGMA.\n");
}

/***********  new level check ***********/

#if DATA_S & N_LINK_TO_FACES

void nf_check(pel,k)
ELEMENT *pel;
INT k;
{
   INT i, j, ind=0;
   NFLINK *pnf;
   
   for(i=0; i < SIDES; i++)
      for(j=0; j < NVERT; j++){
         for (pnf=pel->n[j]->tnfstart; pnf && pnf->nbface != pel->f[i]; 
                                                                 pnf=pnf->next);
         if (pnf==NULL) {
            if (ind==0) printf("\n%i: ",k);
            eprintf("NF_ERROR!!!   ");
            ind = 1;
         }
      }
}

#else

void nf_check(pel,k)
ELEMENT *pel; INT k;
{}

#endif

#if DATA_S & F_LINK_TO_NODES

void fn_check(pel,k)
ELEMENT *pel;
INT k;
{
   INT i, j, ind=0;
   FNLINK *pfn;
   
   for(i=0; i < SIDES; i++)
      for(j=0; j < NVERT; j++){
         for (pfn=pel->f[i]->tfnstart; pfn && pfn->nbnode != pel->n[j]; 
                                                                 pfn=pfn->next);
         if (pfn==NULL) {
            if (ind==0) printf("\n%i:  ",k);
            eprintf("FN_ERROR!!!   ");
            ind = 1;
         }
      }
}

#else

void fn_check(pel,k)
ELEMENT *pel; INT k;
{}

#endif

#if DATA_S & F_LINK_TO_FACES

void ff_check(pel,k)
ELEMENT *pel;
INT k;
{
   INT i, j, ind=0;
   FLINK *pfl;
   
   for(i=0; i < SIDES; i++)
      for(j=0; j < SIDES; j++){
         if (i!=j){
            for (pfl= pel->f[i]->tfstart; pfl && pfl->nbface != pel->f[j];
                                                                 pfl=pfl->next);
            if (pfl==NULL){
               if (ind==0) printf("\n%i:  ",k);
               eprintf("F_ERROR!!!   ");
               ind = 1;
            } 
         }        
      }
}

#else

void ff_check(pel,k)
ELEMENT *pel; INT k;
{}

#endif
 
void rel_check(pel,k)
ELEMENT *pel;
INT k;
{
   INT i, j, ind=0;
   LINK *pli;
   
   nf_check(pel,k);
   fn_check(pel,k);
   ff_check(pel,k);
   for(i=0; i < NVERT; i++)
      for(j=0; j < NVERT; j++){
         if (i!=j){
            for (pli = pel->n[i]->tstart; pli!=NULL && pli->nbnode!=pel->n[j];
                                                                 pli=pli->next);
            if (pli==NULL){
               if (ind==0) printf("\n%i:  ",k);
               eprintf("N_ERROR!!!   ");
               ind = 1;
            } 
            if (NOT_FF(pel->f[i]) && IS_FN(pel->n[j]) && ELEMENT_TYPE == SIMPLEX)
               eprintf("T_ERROR!!!\n");
         }        
      }
}
 
#if DATA_S & N_LINK_TO_FACES

 void node_face_check(mg)
 MULTIGRID *mg;
 {
    NODE *pn;
    NFLINK *pnf, *pnf2;
    
    for (pn = FIRSTN(TOP_GRID(mg)); pn != NULL; pn = pn->succ){
       for (pnf = pn->nfstart; pnf != NULL; pnf = pnf->next){
          if (NOT_FF(pnf->nbface)) 
             eprintf("T_ERROR!\n");       
          for (pnf2 = pnf->next; pnf2 != NULL; pnf2 = pnf2->next)
             if (pnf->nbface == pnf2->nbface)
                eprintf("D_ERROR!\n");
       }
       for (pnf = pn->tnfstart; pnf != pn->nfstart; pnf = pnf->next){
          if (IS_FF(pnf->nbface)) 
             eprintf("T_ERROR!\n");       
          for (pnf2 = pnf->next; pnf2 != NULL; pnf2 = pnf2->next)
             if (pnf->nbface == pnf2->nbface)
                eprintf("D_ERROR!\n");
       }
    }
 }
 
#else

 void node_face_check(mg)
 MULTIGRID *mg;
 {}

#endif

#if DATA_S & F_LINK_TO_NODES

 void face_node_check(mg)
 MULTIGRID *mg;
 {
    FACE *pf;
    FNLINK *pfn, *pfn2;
    
    for (pf = FIRSTF(TOP_GRID(mg)); pf != NULL; pf = pf->succ){
       for (pfn = pf->fnstart; pfn != NULL; pfn = pfn->next){
          if (NOT_FN(pfn->nbnode)) 
             eprintf("T_ERROR!\n");
          for (pfn2 = pfn->next; pfn2 != NULL; pfn2 = pfn2->next)
             if (pfn->nbnode == pfn2->nbnode)
                eprintf("D_ERROR!\n");
       }
       for (pfn = pf->tfnstart; pfn != pf->fnstart; pfn = pfn->next){
          if (IS_FN(pfn->nbnode)) 
             eprintf("T_ERROR!\n");
          for (pfn2 = pfn->next; pfn2 != NULL; pfn2 = pfn2->next)
             if (pfn->nbnode == pfn2->nbnode)
                eprintf("D_ERROR!\n");
       }
    }
 }
 
#else

 void face_node_check(mg)
 MULTIGRID *mg;
 {}

#endif

#if DATA_S & F_LINK_TO_FACES

 void face_face_check(mg)
 MULTIGRID *mg;
 {
    FACE *pf;
    FLINK *pfl, *pfl2;
    
    for (pf = FIRSTF(TOP_GRID(mg)); pf != NULL; pf = pf->succ){
       for (pfl = pf->fstart; pfl != NULL; pfl = pfl->next){
          if (NOT_FF(pfl->nbface)) 
              eprintf("T_ERROR!\n");  
          for(pfl2 = pfl->next; pfl2 != NULL; pfl2 = pfl2->next)
             if (pfl->nbface == pfl2->nbface)
                eprintf("D_ERROR!\n");
       }
       for (pfl = pf->tfstart; pfl != pf->fstart; pfl = pfl->next){
          if (IS_FF(pfl->nbface)) 
              eprintf("T_ERROR!\n");  
          for(pfl2 = pfl->next; pfl2 != NULL; pfl2 = pfl2->next)
             if (pfl->nbface == pfl2->nbface)
                eprintf("D_ERROR!\n");
       }
    }           
 }
 
#else

 void face_face_check(mg)
 MULTIGRID *mg;
 {}

#endif

 void check_new_level(mg)
 MULTIGRID *mg;
 {
    ELEMENT *pel;
    NODE *pn;
    LINK *pli, *pli2;
    INT j=0;
    
    printf("Checking the new level...\n");
    for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel != NULL; pel = pel->succ)
       rel_check(pel,++j);
    for (pn = FIRSTN(TOP_GRID(mg)); pn != NULL; pn = pn->succ){
       for(pli = pn->start; pli != NULL; pli = pli->next){
          if (NOT_FN(pli->nbnode)) 
             eprintf("T_ERROR!\n");
          for (pli2 = pli->next; pli2 != NULL; pli2 = pli2->next)
             if (pli->nbnode == pli2->nbnode)
                eprintf("D_ERROR!\n");
       }
       for(pli = pn->tstart; pli != pn->start; pli = pli->next){
          if (IS_FN(pli->nbnode)) 
             eprintf("T_ERROR!\n");
          for (pli2 = pli->next; pli2 != NULL; pli2 = pli2->next)
             if (pli->nbnode == pli2->nbnode)
                eprintf("D_ERROR!\n");
       }
    }
    node_face_check(mg);
    face_node_check(mg);
    face_face_check(mg);
    printf("New level checked.\n"); 
 }

/***********  end of new level check ***********/

void triangulation(vertexes_and_elements,pgrid,pnode,pvert,pelement,pface,plink,
                  pelink,pnelink,pfnlink,pnflink,pflink,psnode,psface,pp_node,
                  pp_face,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata)
GRID *pgrid;
NODE **pnode;
VERTEX **pvert;
ELEMENT **pelement;
FACE **pface;
LINK **plink;
ELINK **pelink;
NELINK **pnelink;
FNLINK **pfnlink;
NFLINK **pnflink;
FLINK **pflink;
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
void (*vertexes_and_elements)();
{
   FACE *pf;
   ELEMENT *pelem;
   
   (*vertexes_and_elements)(pnode,pvert,pelement,pgrid); 
   new_node_order(pgrid,NMASK);
   make_faces(pface,plink,pgrid,NMASK,FMASK);
   new_element_order(pgrid,FMASK);
   for (pelem = FIRSTELEMENT(pgrid); pelem != NULL; pelem = pelem->succ){
      nodesToNodes(plink,pelem,NMASK);
      nodesToFaces(pfnlink,pelem,NMASK,FMASK);
      facesToNodes(pnflink,pelem,NMASK,FMASK);
      facesToFaces(pflink,pelem,FMASK); 
      elementsToNodes(pnelink,pelem);
   }
   change_order_according_to_itypes(pgrid,&itypes_pos);
   set_outer_normal_vector(pgrid);
   mark_boundary_edges(pgrid);
   compute_number_of_neighbouring_elements_for_nodes(pgrid);
   lg_neighbours(pgrid,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata,
                                                                         U,F,D);
   element_neighbours_on_first_grid(pgrid,pelink);
   element_fneighbours_on_first_grid(pgrid);
   make_previous_nodes(pgrid);
   make_previous_faces(pgrid);
   for (pf = FIRSTF(pgrid); pf != NULL; pf = pf->succ)
      SET_FTOPLEVEL(pf,1);
   if (DATA_S & SPECIAL_NODES_AND_FACES)
      set_type_of_boundary_faces(pgrid,NMASK_FOR_SF,SFMASK);
   make_special_nodes_and_faces(pgrid,psnode,psface,SNMASK,SFMASK);
   make_periodic_nodes(pgrid,pp_node,PNMASK,0);
   FIRSTPF(pgrid) = NULL;
   pgrid->level = 1;
   pgrid->finer = NULL;
   pgrid->coarser = NULL;
   pgrid->first_grid = pgrid->top_grid = pgrid;
}
 
void first_level(vertexes_and_elements,mg,pnode,pvert,pbpoint,pelement,pface,
                 plink,pelink,pnelink,pfnlink,pnflink,pflink,psnode,psface,
                 pp_node,pp_face,pnlglink,plgnlink,plgflink,plglglink,pflglink,
                 plgdata)
MULTIGRID *mg;
NODE **pnode;
VERTEX **pvert;
BPOINT *pbpoint;
ELEMENT **pelement;
FACE **pface;
LINK **plink;
ELINK **pelink;
NELINK **pnelink;
FNLINK **pfnlink;
NFLINK **pnflink;
FLINK **pflink;
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
void (*vertexes_and_elements)();
{ 
   triangulation(vertexes_and_elements,FIRSTGRID(mg),pnode,pvert,pelement,pface,
              plink,pelink,pnelink,pfnlink,pnflink,pflink,psnode,psface,pp_node,
              pp_face,pnlglink,plgnlink,plgflink,plglglink,pflglink,plgdata);
   mg->toplevel = 1;
   mg->nodeNumber = FIRSTGRID(mg)->maxNodeIndex;
   mg->faceNumber = FIRSTGRID(mg)->maxFaceIndex;
   mg->elemNumber = ((long)(*pelement) - (long)(TMIN(FIRSTELEMENT(FIRSTGRID(mg)),
                                         FDBE(FIRSTGRID(mg)))))/sizeof(ELEMENT);
   mg->fictNumber = 0;
   mg->fictElNumber = 0;
   printf("First level: node number = %i, face number = %i,\n",mg->nodeNumber,
                                                               mg->faceNumber);
   printf("             element number = %i\n",mg->elemNumber);
   max_and_min_edge_in_grid(TOP_GRID(mg));
   check_ratio(TOP_GRID(mg));
   boundary_points();
   check_new_level(mg);
   printf("\n");
}

void size_of_data(mg,nodes,snodes,vertexes,bpoints,elements,faces,sfaces,links,
                  elinks,nelinks,flinks,nflinks,fnlinks,pnode,psnode,pvert,
                  pbpoint,pelement,pface,psface,plink,pelink,pnelink,pflink,
                  pnflink,pfnlink,pnlglink,nlglinks,plgnlink,lgnlinks,plgflink,
                  lgflinks,plglglink,lglglinks,pflglink,flglinks,plgdata,
                  lgdatas)
MULTIGRID *mg;
NODE *nodes, *pnode;
SNODE *snodes, *psnode;
VERTEX *vertexes, *pvert;
BPOINT *bpoints, *pbpoint;
ELEMENT *elements, *pelement;
FACE *faces, *pface;
SFACE *sfaces, *psface;
LINK *links, *plink;
ELINK *elinks, *pelink;
NELINK *nelinks, *pnelink;
FNLINK *fnlinks, *pfnlink;
NFLINK *nflinks, *pnflink;
FLINK *flinks, *pflink;
NLGLINK *nlglinks, *pnlglink;
LGNLINK *lgnlinks, *plgnlink;
LGFLINK *lgflinks, *plgflink;
LGLGLINK *lglglinks, *plglglink;
FLGLINK *flglinks, *pflglink;
LGDATA *lgdatas, *plgdata;
{
   printf("NODES:%5.1f %%, VERTEX:%5.1f %%, ELEMENT:%5.1f %%, FACE: %5.1f %%,\n",
     ((long)(pnode)-(long)(nodes))/sizeof(NODE)*100.0/(MAXNODE),
     ((long)(pvert)-(long)(vertexes))/sizeof(VERTEX)*100.0/(MAXVERT),
     ((long)(pelement)-(long)(elements))/sizeof(ELEMENT)*100.0/(MAXELEM),
     ((long)(pface)-(long)(faces))/sizeof(FACE)*100.0/(MAXFACE));
   printf("LINK: %5.1f %%, FNLINK:%5.1f %%, NFLINK: %5.1f %%, FLINK:%5.1f %%, ",
     ((long)(plink)-(long)(links))/sizeof(LINK)*100.0/(MAXLINK),
     ((long)(pfnlink)-(long)(fnlinks))/sizeof(FNLINK)*100.0/(MAXFNLINK),
     ((long)(pnflink)-(long)(nflinks))/sizeof(NFLINK)*100.0/(MAXNFLINK),
     ((long)(pflink)-(long)(flinks))/sizeof(FLINK)*100.0/(MAXFLINK));
   printf("GRID:%5.1f %%\n",((INT)mg->toplevel + 1)*100.0/(MAXLEVEL));
   printf("NLGLINK:  %5.1f %%, LGNLINK:%5.1f %%, LGFLINK: %5.1f %%\n",
     ((long)(pnlglink)-(long)(nlglinks))/sizeof(NLGLINK)*100.0/(MAXNLGLINK),
     ((long)(plgnlink)-(long)(lgnlinks))/sizeof(LGNLINK)*100.0/(MAXLGNLINK),
     ((long)(plgflink)-(long)(lgflinks))/sizeof(LGFLINK)*100.0/(MAXLGFLINK));
   printf("LGLGLINK: %5.1f %%, FLGLINK:%5.1f %%, LGDATA:  %5.1f %%\n",
     ((long)(plglglink)-(long)(lglglinks))/sizeof(LGLGLINK)*100.0/(MAXLGLGLINK),
     ((long)(pflglink)-(long)(flglinks))/sizeof(FLGLINK)*100.0/(MAXFLGLINK),
     ((long)(plgdata)-(long)(lgdatas))/sizeof(LGDATA)*100.0/(MAXLGDATA));
   printf("ELINK: %4.1f %%, NELINK: %4.1f %%, BPOINT: %4.1f %%, SNODE: %4.1f %%, SFACE: %4.1f %%\n",
            ((long)(pelink)-(long)(elinks))/sizeof(ELINK)*100.0/(MAXELINK),
            ((long)(pnelink)-(long)(nelinks))/sizeof(NELINK)*100.0/(MAXNELINK),
            ((long)(pbpoint)-(long)(bpoints))/sizeof(BPOINT)*100.0/(MAXBPOINT),
            ((long)(psnode)-(long)(snodes))/sizeof(SNODE)*100.0/(MAXSNODE),
            ((long)(psface)-(long)(sfaces))/sizeof(SFACE)*100.0/(MAXSFACE));

/*
   printf("NODES: %4i, VERTEX: %4i, ELEMENT: %4i, FACE:  %4i,\n",
     ((long)(pnode)-(long)(nodes))/sizeof(NODE),
     ((long)(pvert)-(long)(vertexes))/sizeof(VERTEX),
     ((long)(pelement)-(long)(elements))/sizeof(ELEMENT),
     ((long)(pface)-(long)(faces))/sizeof(FACE));
   printf("LINK:  %4i, FNLINK: %4i, NFLINK:  %4i, FLINK: %4i, ",
     ((long)(plink)-(long)(links))/sizeof(LINK),
     ((long)(pfnlink)-(long)(fnlinks))/sizeof(FNLINK),
     ((long)(pnflink)-(long)(nflinks))/sizeof(NFLINK),
     ((long)(pflink)-(long)(flinks))/sizeof(FLINK));
   printf("GRID: %4i\n",((INT)mg->toplevel + 1));
   printf("NLGLINK:  %4i, LGNLINK: %4i, LGFLINK: %4i\n",
     ((long)(pnlglink)-(long)(nlglinks))/sizeof(NLGLINK),
     ((long)(plgnlink)-(long)(lgnlinks))/sizeof(LGNLINK),
     ((long)(plgflink)-(long)(lgflinks))/sizeof(LGFLINK));
   printf("LGLGLINK: %4i, FLGLINK: %4i, LGDATA:  %4i\n",
     ((long)(plglglink)-(long)(lglglinks))/sizeof(LGLGLINK),
     ((long)(pflglink)-(long)(flglinks))/sizeof(FLGLINK),
     ((long)(plgdata)-(long)(lgdatas))/sizeof(LGDATA));
   printf("ELINK: %4i\n",((long)(pelink)-(long)(elinks))/sizeof(ELINK));
*/
 
   if ( pnode > nodes + MAXNODE )
                          eprintf("Error: maximum number of nodes exceeded!\n");
   if ( psnode > snodes + MAXSNODE )
                         eprintf("Error: maximum number of snodes exceeded!\n");
   if ( pvert > vertexes + MAXVERT )
                       eprintf("Error: maximum number of vertexes exceeded!\n");
   if ( pbpoint > bpoints + MAXBPOINT )
                       eprintf("Error: maximum number of bpoints exceeded!\n");
   if ( pelement > elements + MAXELEM )
                       eprintf("Error: maximum number of elements exceeded!\n");
   if ( pface > faces + MAXFACE )
                          eprintf("Error: maximum number of faces exceeded!\n");
   if ( psface > sfaces + MAXSFACE )
                         eprintf("Error: maximum number of sfaces exceeded!\n");
   if ( plink > links + MAXLINK )
                          eprintf("Error: maximum number of links exceeded!\n");
   if ( pelink > elinks + MAXELINK )
                         eprintf("Error: maximum number of elinks exceeded!\n");
   if ( pnelink > nelinks + MAXNELINK )
                        eprintf("Error: maximum number of nelinks exceeded!\n");
   if ( pflink > flinks + MAXFLINK )
                         eprintf("Error: maximum number of flinks exceeded!\n");
   if ( pnflink > nflinks + MAXNFLINK )
                        eprintf("Error: maximum number of nflinks exceeded!\n");
   if ( pfnlink > fnlinks + MAXFNLINK )
                        eprintf("Error: maximum number of fnlinks exceeded!\n");
   if ( pnlglink > nlglinks + MAXNLGLINK )
                       eprintf("Error: maximum number of nlglinks exceeded!\n");
   if ( plgnlink > lgnlinks + MAXLGNLINK )
                       eprintf("Error: maximum number of lgnlinks exceeded!\n");
   if ( plgflink > lgflinks + MAXLGFLINK )
                       eprintf("Error: maximum number of lgflinks exceeded!\n");
   if ( plglglink > lglglinks + MAXLGLGLINK )
                      eprintf("Error: maximum number of lglglinks exceeded!\n");
   if ( pflglink > flglinks + MAXFLGLINK )
                       eprintf("Error: maximum number of flglinks exceeded!\n");
   if ( plgdata > lgdatas + MAXLGDATA )
                        eprintf("Error: maximum number of lgdatas exceeded!\n");
   if ( (INT)mg->toplevel >= MAXLEVEL )
                          eprintf("Error: maximum number of grids exceeded!\n");
}
