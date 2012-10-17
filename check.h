/******************************************************************************/
/*                                                                            */
/*                 comparison with data obtained on one grid                  */
/*                                                                            */
/******************************************************************************/

/* usage:  1. perform a computation with grid refinement and save the
              finest triangulation.
           2. perform a computation on one grid obtained by reading the finest
              triangulation and save all necessary data by "Save_results".
              Now, a second computation on one grid with the saved matrix can
              be performed. The saved data can be read by "Read_results".
           3. perform the computation with grid refinement again and compare
              obtained data with data saved in 2 ("compare_results").
              It can be also performed a computation with data structures 
              obtained by the nonuniform refinement and with matrix and rhs 
              saved in 2. The data are read by "read_matrix".                 */

#if !(DATA_STR & KORN_MATRIX) && !(DATA_STR & REDUCED) && (DIM == 3)

 void Save_results(mg,vertexes,pvert)
 MULTIGRID *mg;
 VERTEX *vertexes, *pvert;
 {
    ELEMENT *pel;
    VERTEX *pv;
    NODE *pnode;
    FACE *pface;
    INT sv;
    LINK *pli;
    FLINK *pfl;
    NFLINK *pnf;
    FNLINK *pfn;
    FILE *fp;
    
    fp = fopen("data.res","w");
    sv = sizeof(VERTEX);
    pv = vertexes;
    fprintf(fp,"%ld \n",((long)(pvert) - (long)(vertexes))/sv);
    while(pv < pvert){
       fprintf(fp,"%e %e %e %e %e %e\n",pv->x[0],pv->x[1],pv->x[2],
                    ND(pv->topnode,U,0),ND(pv->topnode,U,1),ND(pv->topnode,U,2));
       pv++;
    }
    for (pel = FIRSTELEMENT(FIRSTGRID(mg)); pel != NULL; pel = pel->succ)
       fprintf(fp,"%ld %ld %ld %ld %i %i %i %i %e %e %e %e\n",
          ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
          ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
          ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv,
          ((long)(pel->n[3]->myvertex)  - (long)(vertexes))/sv,
          pel->f[0]->index,pel->f[1]->index,
          pel->f[2]->index,pel->f[3]->index,
          FD(pel->f[0],U),FD(pel->f[1],U),
          FD(pel->f[2],U),FD(pel->f[3],U));
    fprintf(fp,"%i %i %i %i %i %i %i %i %e %e %e %e\n",
                                    -3,-3,-3,-3,-3,-3,-3,-3,-3.0,-3.0,-3.0,-3.0);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ){
       fprintf(fp,"%ld %e  ",
           ((long)(pnode->myvertex)  - (long)(vertexes))/sv,COEFFN(pnode,A));
       for (pli = pnode->start; pli != NULL; pli = pli->next)
          fprintf(fp,"%ld %e  ",
           ((long)(pli->nbnode->myvertex)  - (long)(vertexes))/sv,COEFFL(pli,A));
       fprintf(fp,"%i\n",-1);
    }
    fprintf(fp,"%i\n",-1);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ){
       fprintf(fp,"%i %e  ",pface->index,COEFF_FF(pface,A));
       for (pfl = pface->fstart; pfl != NULL; pfl = pfl->next)
          fprintf(fp,"%i %e  ",pfl->nbface->index,COEFF_FL(pfl,A));
       fprintf(fp,"%i\n",-1);
    }
    fprintf(fp,"%i\n",-1);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ){
       fprintf(fp,"%ld ",((long)(pnode->myvertex)  - (long)(vertexes))/sv);
       for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next)
          fprintf(fp,"%i %e %e %e  ",pnf->nbface->index,COEFF_NF(pnf,A,0),
                                            COEFF_NF(pnf,A,1),COEFF_NF(pnf,A,2));
       fprintf(fp,"%i\n",-1);
    }
    fprintf(fp,"%i\n",-1);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ){
       fprintf(fp,"%i  ",pface->index);
       for (pfn = pface->fnstart; pfn != NULL; pfn = pfn->next)
          fprintf(fp,"%ld %e %e %e  ",
                          ((long)(pfn->nbnode->myvertex)  - (long)(vertexes))/sv,
                          COEFF_FN(pfn,A,0),COEFF_FN(pfn,A,1),COEFF_FN(pfn,A,2));
       fprintf(fp,"%i\n",-1);
    }
    fprintf(fp,"%i\n",-1);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ)
       fprintf(fp,"%ld %e %e %e\n",((long)(pnode->myvertex) - 
                 (long)(vertexes))/sv,ND(pnode,F,0),ND(pnode,F,1),ND(pnode,F,2));
    fprintf(fp,"%i %e %e %e\n",-1,1.0,1.0,1.0);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ)
       fprintf(fp,"%i %e\n",pface->index,FD(pface,F));
    fprintf(fp,"%i %e\n",-1,1.0);
    
    fclose(fp);
 }
 
 void Read_results(mg,vertexes,pvert)
 MULTIGRID *mg;
 VERTEX *vertexes, *pvert;
 {
    ELEMENT *pel;
    VERTEX *pv;
    NODE *pnode;
    FACE *pface;
    LINK *pli;
    FLINK *pfl;
    NFLINK *pnf;
    FNLINK *pfn;
    FILE *fp;
    FLOAT xx, *x;
    INT kk, *k;
    
    x = &xx;
    k = &kk;
    
    fp = fopen("data.res","r");
    pv = vertexes;
    fscanf(fp,"%i",k);
    while(pv < pvert){
       fscanf(fp,"%le %le %le %le %le %le",x,x,x,x,x,x);
       pv++;
    }
    for (pel = FIRSTELEMENT(FIRSTGRID(mg)); pel != NULL; pel = pel->succ)
       fscanf(fp,"%i %i %i %i %i %i %i %i %le %le %le %le",
                                                        k,k,k,k,k,k,k,k,x,x,x,x);
    fscanf(fp,"%i %i %i %i %i %i %i %i %le %le %le %le",k,k,k,k,k,k,k,k,x,x,x,x);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ){
       fscanf(fp,"%i %le",k,&(COEFFN(pnode,A)));
       for (pli = pnode->start; pli != NULL; pli = pli->next)
          fscanf(fp,"%i %le",k,&(COEFFL(pli,A)));
       fscanf(fp,"%i",k);
    }
    fscanf(fp,"%i",k);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ){
       fscanf(fp,"%i %le",k,&(COEFF_FF(pface,A)));
       for (pfl = pface->fstart; pfl != NULL; pfl = pfl->next)
          fscanf(fp,"%i %le",k,&(COEFF_FL(pfl,A)));
       fscanf(fp,"%i",k);
    }
    fscanf(fp,"%i",k);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ){
       fscanf(fp,"%i",k);
       for (pnf = pnode->nfstart; pnf != NULL; pnf = pnf->next)
          fscanf(fp,"%i %le %le %le",
               k,&(COEFF_NF(pnf,A,0)),&(COEFF_NF(pnf,A,1)),&(COEFF_NF(pnf,A,2)));
       fscanf(fp,"%i",k);
    }
    fscanf(fp,"%i",k);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ){
       fscanf(fp,"%i",k);
       for (pfn = pface->fnstart; pfn != NULL; pfn = pfn->next)
          fscanf(fp,"%i %le %le %le",
               k,&(COEFF_FN(pfn,A,0)),&(COEFF_FN(pfn,A,1)),&(COEFF_FN(pfn,A,2)));
       fscanf(fp,"%i",k);
    }
    fscanf(fp,"%i",k);
    
    for (pnode = FIRSTNODE(FIRSTGRID(mg)); pnode != NULL; pnode = pnode->succ)
       fscanf(fp,"%i %le %le %le",
                           k,&(ND(pnode,F,0)),&(ND(pnode,F,1)),&(ND(pnode,F,2)));
    fscanf(fp,"%i %le %le %le",k,x,x,x);
    
    for (pface = FIRSTFACE(FIRSTGRID(mg)); pface != NULL; pface = pface->succ)
       fscanf(fp,"%i %le",k,&(FD(pface,F)));
    fscanf(fp,"%i %le",k,x);
    
    fclose(fp);
 }

 void set_y(j1,j2,j3,j4,l1,y1,pel,yf,pf,ind)
 INT j1,j2,j3,j4,l1,yf,*ind;
 FLOAT y1;
 ELEMENT *pel;
 FACE **pf;
 {
    INT k,l;
    
    if (j1==pel->n[0]->index) k=0;
    else if (j1==pel->n[1]->index) k=1;
    else if (j1==pel->n[2]->index) k=2;
    else if (j1==pel->n[3]->index) k=3;
    else eprintf("ERROR3 in compare_results.\n");
    if (j2 < j3)
       if (j3 < j4 || j4 < j2) l = 1;
       else l = -1;
    else if (j3 < j4 && j4 < j2) l = 1;
    else l = -1;
    if (l > 0 && (FTYPE(pel->f[k]) & ORIENT) == 0 || 
        l < 0 && (FTYPE(pel->f[k]) & ORIENT)){
       FD(top_face(pel->f[k]),yf) = y1;
       ind[l1] = 1;
    }
    else{
       FD(top_face(pel->f[k]),yf) = -y1;
       ind[l1] = -1;
    }
    pf[l1] = top_face(pel->f[k]);
 }
 
 void test_err(x,y,amax,rmax,min,max,xr,yr)
 FLOAT x, y, *amax, *rmax, *min, *max, *xr, *yr;
 {
    FLOAT r, r1, r2;
    
    if (fabs(x) > 0.0 && fabs(x) < *min) *min = fabs(x);
    if (fabs(x) > *max) *max = fabs(x);
    if ((r=fabs(y-x)) > *amax) *amax = r;
    if (fabs(x) > 0.0 && (r1=r/fabs(x)) > *rmax){
       *rmax = r1;
       *xr = x;
       *yr = y;
    }
    if (fabs(y) > 0.0 && (r2=r/fabs(y)) > *rmax){
       *rmax = r2;
       *xr = x;
       *yr = y;
    }
/*    if ((fabs(x) > 0.0 && r1 > 0.1) || (fabs(y) > 0.0 && r2 > 0.1))
       printf("Big relative error between   %e  and  %e.\n",x,y);     */
 }
  
 void print_err(amax,rmax,min,max,xr,yr)
 FLOAT *amax, *rmax, *min, *max, xr, yr;
 {
    printf("---------------------------------------\n");
    printf("Maximum absolute error: %f\n",*amax);
    printf("Maximum relative error: %f   (between %e and %e)\n",*rmax,xr,yr);
    printf("min. nonzero element: %f\n",*min);
    printf("max. element: %f\n",*max);
    *amax = -1.0;
    *rmax = -1.0;
    *min = 100000000.0;
    *max = 0.0;
 }
 
 void compare_results(mg,vertexes,pvert,y,yf,o)
 MULTIGRID *mg;
 VERTEX *vertexes, *pvert;
 INT y, yf,o;
 {
    ELEMENT *pelem, *pel;
    GRID *theGrid, *pg;
    VERTEX *pv;
    INT i, j, n;  
    FILE *fp;
    INT k, i1, i2, i3, i4, j1, j2, j3, j4, l1, l2, l3, l4, *ind;
    NODE *pnode, *n1, *n2, *n3, *n4, *pnj, **pn;
    FACE *pface, *pfj, **pf;
    LINK *pli;
    FLINK *pfl;
    NFLINK *pnf;
    FNLINK *pfn;
    FLOAT x, x1, x2, x3, y1, y2, y3, y4, amax=0.0, rmax=-1.0, max=0.0, 
          min=100000000.0, xr, yr;
    
    ind = (INT *) calloc(mg->faceNumber,sizeof(INT));
    pn  = (NODE **) calloc(mg->nodeNumber,sizeof(NODE *));
    pf  = (FACE **) calloc(mg->faceNumber,sizeof(FACE *));
    
    fp = fopen("data.res","r");
    fscanf(fp,"%i",&n);
    for (i=0; i<n; i++){
       fscanf(fp,"%lf %lf %lf %lf %lf %lf",&x1,&x2,&x3,&y1,&y2,&y3);
       pv = vertexes;
       while(pv < pvert && 
         (fabs(pv->x[0]-x1)>EPS || fabs(pv->x[1]-x2)>EPS || 
                                   fabs(pv->x[2]-x3)>EPS)) pv++;
       if (pv < pvert){
          pn[i] = pv->topnode;
          ND(pv->topnode,y,0) = y1;
          ND(pv->topnode,y,1) = y2;
          ND(pv->topnode,y,2) = y3;
       }
       else
          eprintf("ERROR1 in compare_results.\n");
    }
    fscanf(fp,"%i %i %i %i %i %i %i %i %lf %lf %lf %lf",
                                &i1,&i2,&i3,&i4,&l1,&l2,&l3,&l4,&y1,&y2,&y3,&y4);
    while(i1 > -1){
     k = 1;
     for (theGrid = FIRSTGRID(mg); theGrid != NULL && k; theGrid = theGrid->finer){
       n1 = pn[i1];
       n2 = pn[i2];
       n3 = pn[i3];
       n4 = pn[i4];
       while (n1->index > theGrid->maxNodeIndex) n1 = n1->father;
       while (n2->index > theGrid->maxNodeIndex) n2 = n2->father;
       while (n3->index > theGrid->maxNodeIndex) n3 = n3->father;
       while (n4->index > theGrid->maxNodeIndex) n4 = n4->father;
       j1 = n1->index;
       j2 = n2->index;
       j3 = n3->index;
       j4 = n4->index;
       order4(&j1,&j2,&j3,&j4);
       for (pelem = FIRSTELEMENT(theGrid); pelem != NULL && k; pelem = pelem->succ)
          if (IS_TOP_ELEMENT(pelem) && pelem->n[0]->index == j1 && 
              pelem->n[1]->index == j2 && pelem->n[2]->index == j3 && 
              pelem->n[3]->index == j4){
                    pel = pelem;
                    k = 0;
          }
     }
     if (k==1){
        eprintf("ERROR2 in compare_results.\n"); 
        printf("i:   %i %i %i %i\n",i1,i2,i3,i4);
        printf("j:   %i %i %i %i\n",j1,j2,j3,j4);
     }      
     else{
        j1 = n1->index;
        j2 = n2->index;
        j3 = n3->index;
        j4 = n4->index;
        set_y(j1,j2,j3,j4,l1,y1,pel,yf,pf,ind);
        set_y(j2,j3,j4,j1,l2,y2,pel,yf,pf,ind);
        set_y(j3,j4,j1,j2,l3,y3,pel,yf,pf,ind);
        set_y(j4,j1,j2,j3,l4,y4,pel,yf,pf,ind);
     }
     fscanf(fp,"%i %i %i %i %i %i %i %i %lf %lf %lf %lf",
                                &i1,&i2,&i3,&i4,&l1,&l2,&l3,&l4,&y1,&y2,&y3,&y4);
    } 
    
    printf("Comparison with a solution obtained on one grid:\n");
    printf("------------------------------------------------\n");
    printf("Nodes:\n");
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pnode = FIRSTNODE(theGrid); pnode != NULL; pnode = pnode->succ)
          if (IS_TOP_NODE(pnode)){
             if (o) printf("%f  %f\n",ND(pnode,U,0)-ND(pnode,y,0),
                                    (ND(pnode,U,0)-ND(pnode,y,0))/ND(pnode,y,0));
             if (o) printf("%f  %f\n",ND(pnode,U,1)-ND(pnode,y,1),
                                    (ND(pnode,U,1)-ND(pnode,y,1))/ND(pnode,y,1));
             if (o) printf("%f  %f\n",ND(pnode,U,2)-ND(pnode,y,2),
                                    (ND(pnode,U,2)-ND(pnode,y,2))/ND(pnode,y,2));
             test_err(ND(pnode,U,0),ND(pnode,y,0),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(ND(pnode,U,1),ND(pnode,y,1),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(ND(pnode,U,2),ND(pnode,y,2),&amax,&rmax,&min,&max,&xr,&yr);
          }  
    print_err(&amax,&rmax,&min,&max,xr,yr);      
    printf("Faces:\n");
    for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
       for (pface = FIRSTFACE(theGrid); pface != NULL; pface = pface->succ)
          if (IS_TOP_FACE(pface)){
             if (o) printf("%f  %f\n",FD(pface,U)-FD(pface,yf),
                                        (FD(pface,U)-FD(pface,yf))/FD(pface,yf));
             test_err(FD(pface,U),FD(pface,yf),&amax,&rmax,&min,&max,&xr,&yr);
          }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    
    printf("Errors in the matrix:\n");
    printf("---------------------\n");
    printf("nodes-nodes:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%lf",&x);
       if(o) printf("\n%f  %f\n",COEFFN(pn[i],A)-x,(COEFFN(pn[i],A)-x)/x);
       test_err(x,COEFFN(pn[i],A),&amax,&rmax,&min,&max,&xr,&yr);
       pg = TOP_GRID(mg);
       if(!ON_FINEST_GRID(pn[i],mg))
          while (pn[i]->index < pg->minNodeIndex) pg = pg->coarser;
       fscanf(fp,"%i",&j);
       while(j >= 0){
          fscanf(fp,"%lf",&x);
          pnj = pn[j];
          while (pnj->index > pg->maxNodeIndex) pnj = pnj->father;
          for(pli = pn[i]->start; pli!=NULL && pli->nbnode != pnj; pli=pli->next);
          if (pli){
             if(o) printf("%f  %f\n",COEFFL(pli,A)-x,(COEFFL(pli,A)-x)/x);
             test_err(x,COEFFL(pli,A),&amax,&rmax,&min,&max,&xr,&yr);
          }
          else
             eprintf("ERROR5 in compare_results.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    printf("\nfaces-faces:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%lf",&x);
       if(o) printf("\n%f  %f\n",COEFF_FF(pf[i],A)-x,(COEFF_FF(pf[i],A)-x)/x);
       test_err(x,COEFF_FF(pf[i],A),&amax,&rmax,&min,&max,&xr,&yr);
       pg = TOP_GRID(mg);
       if(!F_ON_FINEST_GRID(pf[i],mg))
          while (pf[i]->index < pg->minFaceIndex) pg = pg->coarser;
       fscanf(fp,"%i",&j);
       while(j >= 0){
          fscanf(fp,"%lf",&x);
          x = x*ind[i]*ind[j];
          pfj = pf[j];
          while (pfj->index > pg->maxFaceIndex) pfj = pfj->father;
          for(pfl = pf[i]->fstart; pfl!=NULL && pfl->nbface != pfj; pfl=pfl->next);
          if (pfl){
             if(o) printf("%f  %f\n",COEFF_FL(pfl,A)-x,(COEFF_FL(pfl,A)-x)/x);
             test_err(x,COEFF_FL(pfl,A),&amax,&rmax,&min,&max,&xr,&yr);
          }
          else
             eprintf("ERROR5 in compare_results.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    printf("\nnodes-faces:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%i",&j);
       pg = TOP_GRID(mg);
       if(!ON_FINEST_GRID(pn[i],mg))
          while (pn[i]->index < pg->minNodeIndex) pg = pg->coarser;
       while(j >= 0){
          fscanf(fp,"%lf %lf %lf",&x1,&x2,&x3);
          x1 = x1*ind[j];
          x2 = x2*ind[j];
          x3 = x3*ind[j];
          pfj = pf[j];
          while (pfj->index > pg->maxFaceIndex) pfj = pfj->father;          
          for(pnf = pn[i]->nfstart; pnf!=NULL && pnf->nbface != pfj; pnf=pnf->next);
          if (pnf){
             if(o) printf("%f  %f  %f      %f  %f  %f\n",
                  COEFF_NF(pnf,A,0)-x1,COEFF_NF(pnf,A,1)-x2,COEFF_NF(pnf,A,2)-x3,
                  (COEFF_NF(pnf,A,0)-x1)/x1,(COEFF_NF(pnf,A,1)-x2)/x2,
                                                      (COEFF_NF(pnf,A,2)-x3)/x3);
             test_err(x1,COEFF_NF(pnf,A,0),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(x2,COEFF_NF(pnf,A,1),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(x3,COEFF_NF(pnf,A,2),&amax,&rmax,&min,&max,&xr,&yr);
          }
          else
             eprintf("ERROR5 in compare_results.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
       if (o) printf("\n");
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    printf("\nfaces-nodes:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%i",&j);
       pg = TOP_GRID(mg);
       if(!F_ON_FINEST_GRID(pf[i],mg))
          while (pf[i]->index < pg->minFaceIndex) pg = pg->coarser;
       while(j >= 0){
          fscanf(fp,"%lf %lf %lf",&x1,&x2,&x3); 
          x1 = x1*ind[i];
          x2 = x2*ind[i];
          x3 = x3*ind[i];
          pnj = pn[j];
          while (pnj->index > pg->maxNodeIndex) pnj = pnj->father;          
          for(pfn = pf[i]->fnstart; pfn!=NULL && pfn->nbnode != pnj; pfn=pfn->next);
          if (pfl){
             if(o) printf("%f %f %f      %f %f %f\n",
              COEFF_FN(pfn,A,0)-x1,COEFF_FN(pfn,A,1)-x2,COEFF_FN(pfn,A,2)-x3,
              (COEFF_FN(pfn,A,0)-x1)/x1,(COEFF_FN(pfn,A,1)-x2)/x2,
                                                      (COEFF_FN(pfn,A,2)-x3)/x3);
             test_err(x1,COEFF_FN(pfn,A,0),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(x2,COEFF_FN(pfn,A,1),&amax,&rmax,&min,&max,&xr,&yr);
             test_err(x3,COEFF_FN(pfn,A,2),&amax,&rmax,&min,&max,&xr,&yr);
          }
          else
             eprintf("ERROR5 in compare_results.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
       if (o) printf("\n");
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    printf("\nrhs:\n");
    printf("------\n");
    printf("Nodes:\n");
    fscanf(fp,"%i %lf %lf %lf",&i,&x1,&x2,&x3);
    while(i >= 0){
       if (o) printf("%f %f %f      %f %f %f\n",ND(pn[i],F,0)-x1,
                         ND(pn[i],F,1)-x2,ND(pn[i],F,2)-x3,(ND(pn[i],F,0)-x1)/x1,
                                    (ND(pn[i],F,1)-x2)/x2,(ND(pn[i],F,2)-x3)/x3);
       test_err(x1,ND(pn[i],F,0),&amax,&rmax,&min,&max,&xr,&yr);
       test_err(x2,ND(pn[i],F,1),&amax,&rmax,&min,&max,&xr,&yr);
       test_err(x3,ND(pn[i],F,2),&amax,&rmax,&min,&max,&xr,&yr);
       fscanf(fp,"%i %lf %lf %lf",&i,&x1,&x2,&x3);
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    printf("Faces:\n");
    fscanf(fp,"%i %lf",&i,&x);
    while(i >= 0){
       x = x*ind[i];
       if (o) printf("%f %f\n",FD(pf[i],F)-x,(FD(pf[i],F)-x)/x);
       test_err(x,FD(pf[i],F),&amax,&rmax,&min,&max,&xr,&yr);
       fscanf(fp,"%i %lf",&i,&x);
    }
    print_err(&amax,&rmax,&min,&max,xr,yr);
    fclose(fp);  
 }
 
 void read_matrix(mg,vertexes,pvert,yf)
 MULTIGRID *mg;
 VERTEX *vertexes, *pvert;
 INT yf;
 {
    ELEMENT *theElem, *pel;
    GRID *theGrid, *pg;
    VERTEX *pv;
    INT i, j, n;  
    FILE *fp;
    INT k, i1, i2, i3, i4, j1, j2, j3, j4, l1, l2, l3, l4, ind[10000];
    NODE *n1, *n2, *n3, *n4, *pnj, *pn[10000];
    FACE *pfj, *pf[10000];
    LINK *pli;
    FLINK *pfl;
    NFLINK *pnf;
    FNLINK *pfn;
    FLOAT x, x1, x2, x3, y1, y2, y3, y4;
    
    fp = fopen("data.res","r");
    fscanf(fp,"%i",&n);
    for (i=0; i<n; i++){
       fscanf(fp,"%lf %lf %lf %lf %lf %lf",&x1,&x2,&x3,&y1,&y2,&y3);
       pv = vertexes;
       while(pv < pvert && 
         (fabs(pv->x[0]-x1)>EPS || fabs(pv->x[1]-x2)>EPS || 
                                   fabs(pv->x[2]-x3)>EPS)) pv++;
       if (pv < pvert)
          pn[i] = pv->topnode;
       else
          eprintf("ERROR1 in read_matrix.\n");
    }
    fscanf(fp,"%i %i %i %i %i %i %i %i %lf %lf %lf %lf",
                                &i1,&i2,&i3,&i4,&l1,&l2,&l3,&l4,&y1,&y2,&y3,&y4);
    while(i1 > -1){
     k = 1;
     for (theGrid = FIRSTGRID(mg); theGrid != NULL && k; theGrid = theGrid->finer){
       n1 = pn[i1];
       n2 = pn[i2];
       n3 = pn[i3];
       n4 = pn[i4];
       while (n1->index > theGrid->maxNodeIndex) n1 = n1->father;
       while (n2->index > theGrid->maxNodeIndex) n2 = n2->father;
       while (n3->index > theGrid->maxNodeIndex) n3 = n3->father;
       while (n4->index > theGrid->maxNodeIndex) n4 = n4->father;
       j1 = n1->index;
       j2 = n2->index;
       j3 = n3->index;
       j4 = n4->index;
       order4(&j1,&j2,&j3,&j4);
       for (theElem = FIRSTELEMENT(theGrid); theElem != NULL && k; 
                                                             theElem = theElem->succ)
          if (IS_TOP_ELEMENT(theElem) && theElem->n[0]->index == j1 && 
              theElem->n[1]->index == j2 && theElem->n[2]->index == j3 && 
                                            theElem->n[3]->index == j4){
                    pel = theElem;
                    k = 0;
          }
     }
     if (k==1){
        eprintf("ERROR2 in read_matrix.\n"); 
        printf("i:   %i %i %i %i\n",i1,i2,i3,i4);
        printf("j:   %i %i %i %i\n",j1,j2,j3,j4);
     }      
     else{
        j1 = n1->index;
        j2 = n2->index;
        j3 = n3->index;
        j4 = n4->index;
        set_y(j1,j2,j3,j4,l1,y1,pel,yf,pf,ind);
        set_y(j2,j3,j4,j1,l2,y2,pel,yf,pf,ind);
        set_y(j3,j4,j1,j2,l3,y3,pel,yf,pf,ind);
        set_y(j4,j1,j2,j3,l4,y4,pel,yf,pf,ind);
     }
     fscanf(fp,"%i %i %i %i %i %i %i %i %lf %lf %lf %lf",
                                &i1,&i2,&i3,&i4,&l1,&l2,&l3,&l4,&y1,&y2,&y3,&y4);
    } 
    
    printf("read_matrix:\n");
    printf("nodes-nodes:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%lf",&x);
       COEFFN(pn[i],A)=x;
       pg = TOP_GRID(mg);
       if(!ON_FINEST_GRID(pn[i],mg))
          while (pn[i]->index < pg->minNodeIndex) pg = pg->coarser;
       fscanf(fp,"%i",&j);
       while(j >= 0){
          fscanf(fp,"%lf",&x);
          pnj = pn[j];
          while (pnj->index > pg->maxNodeIndex) pnj = pnj->father;
          for(pli = pn[i]->start; pli!=NULL && pli->nbnode != pnj; pli=pli->next);
          if (pli)
             COEFFL(pli,A)=x;
          else
             eprintf("ERROR5 in read_matrix.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    printf("\nfaces-faces:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%lf",&x);
       COEFF_FF(pf[i],A)=x;
       pg = TOP_GRID(mg);
       if(!F_ON_FINEST_GRID(pf[i],mg))
          while (pf[i]->index < pg->minFaceIndex) pg = pg->coarser;
       fscanf(fp,"%i",&j);
       while(j >= 0){
          fscanf(fp,"%lf",&x);
          x = x*ind[i]*ind[j];
          pfj = pf[j];
          while (pfj->index > pg->maxFaceIndex) pfj = pfj->father;
          for(pfl = pf[i]->fstart; pfl!=NULL && pfl->nbface != pfj; pfl=pfl->next);
          if (pfl)
             COEFF_FL(pfl,A)=x;
          else
             eprintf("ERROR5 in read_matrix.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    printf("\nnodes-faces:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%i",&j);
       pg = TOP_GRID(mg);
       if(!ON_FINEST_GRID(pn[i],mg))
          while (pn[i]->index < pg->minNodeIndex) pg = pg->coarser;
       while(j >= 0){
          fscanf(fp,"%lf %lf %lf",&x1,&x2,&x3);
          x1 = x1*ind[j];
          x2 = x2*ind[j];
          x3 = x3*ind[j];
          pfj = pf[j];
          while (pfj->index > pg->maxFaceIndex) pfj = pfj->father;          
          for(pnf = pn[i]->nfstart; pnf!=NULL && pnf->nbface != pfj; pnf=pnf->next);
          if (pnf){
             COEFF_NF(pnf,A,0)=x1;
             COEFF_NF(pnf,A,1)=x2;
             COEFF_NF(pnf,A,2)=x3;
          }
          else
             eprintf("ERROR5 in read_matrix.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    printf("\nfaces-nodes:\n");
    fscanf(fp,"%i",&i);
    while (i >= 0){
       fscanf(fp,"%i",&j);
       pg = TOP_GRID(mg);
       if(!F_ON_FINEST_GRID(pf[i],mg))
          while (pf[i]->index < pg->minFaceIndex) pg = pg->coarser;
       while(j >= 0){
          fscanf(fp,"%lf %lf %lf",&x1,&x2,&x3); 
          x1 = x1*ind[i];
          x2 = x2*ind[i];
          x3 = x3*ind[i];
          pnj = pn[j];
          while (pnj->index > pg->maxNodeIndex) pnj = pnj->father;          
          for(pfn = pf[i]->fnstart; pfn!=NULL && pfn->nbnode != pnj; pfn=pfn->next);
          if (pfl){
             COEFF_FN(pfn,A,0)=x1;
             COEFF_FN(pfn,A,1)=x2;
             COEFF_FN(pfn,A,2)=x3;
          }
          else
             eprintf("ERROR5 in read_matrix.\n");
          fscanf(fp,"%i",&j);
       }
       fscanf(fp,"%i",&i);
    }
    printf("\nrhs:\n");
    printf("------\n");
    printf("Nodes:\n");
    fscanf(fp,"%i %lf %lf %lf",&i,&x1,&x2,&x3);
    while(i >= 0){
       ND(pn[i],F,0)=x1;
       ND(pn[i],F,1)=x2;
       ND(pn[i],F,2)=x3;
       fscanf(fp,"%i %lf %lf %lf",&i,&x1,&x2,&x3);
    }
    printf("Faces:\n");
    fscanf(fp,"%i %lf",&i,&x);
    while(i >= 0){
       x = x*ind[i];
       FD(pf[i],F)=x;
       fscanf(fp,"%i %lf",&i,&x);
    }
    fclose(fp);  
 }

void tests(mg,nodes,snodes,p_nodes,vertexes,elements,faces,sfaces,p_faces,links,
           elinks,flinks,nflinks,fnlinks,pnode,psnode,pp_node,pvert,pelement,
           pface,psface,pp_face,plink,pelink,pflink,pnflink,pfnlink,pnlglink,
           nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,plglglink,lglglinks,
           pflglink,flglinks,plgdata,lgdatas)
MULTIGRID *mg;
NODE *nodes, *pnode;
NODE *snodes, *psnode;
P_NODE *p_nodes, *pp_node;
VERTEX *vertexes, *pvert;
ELEMENT *elements, *pelement;
FACE *faces, *pface;
FACE *sfaces, *psface;
P_FACE *p_faces, *pp_face;
LINK *links, *plink;
ELINK *elinks, *pelink;
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
   INT i, j;  
   FLOAT err; 
   
   printf("1 ... computation with refinement and saving the\n");
   printf("        finest triangulation.\n");
   printf("2 ... computation on one grid obtained by reading\n");
   printf("        the finest triangulation; saving the res-\n");
   printf("        ults, matrix, rhs, and other necessary data.\n");
   printf("3 ... computation on one grid with matrix and rhs\n");
   printf("        saved in 2; saving the results, matrix,\n"); 
   printf("        rhs, and other necessary data.\n");
   printf("4 ... computation with nonuniform refinement and\n");
   printf("        comparing the results with results\n"); 
   printf("        obtained in 2 or 3.\n");
   printf("5 ... computation with data sructures constructed\n");
   printf("        by nonuniform refinement and with matrix\n");
   printf("        and rhs saved in 2 or 3. Comparisson with\n");
   printf("        results saved in 2 or 3.\n");
   do{
      printf("Input: ");
      scanf("%i",&j);
   } while(j < 1 || j > 5);

   if (j == 2 || j == 3)
      first_level(read_vertexes_and_elements,
                      mg,&pnode,&pvert,&pbpoint,&pelement,&pface,&plink,&pelink,
                      &pfnlink,&pnflink,&pflink,&psnode,&psface,&pp_node,
                      &pp_face,&pnlglink,
                      &plgnlink,&plgflink,&plglglink,&pflglink,&plgdata);
   else
      first_level(cube,
                      mg,&pnode,&pvert,&pbpoint,&pelement,&pface,&plink,&pelink,
                      &pfnlink,&pnflink,&pflink,&psnode,&psface,&pp_node,
                      &pp_face,&pnlglink,
                      &plgnlink,&plgflink,&plglglink,&pflglink,&plgdata);
              
   size_of_data(mg,nodes,snodes,vertexes,bpoints,elements,faces,
                  sfaces,links,elinks,flinks,nflinks,fnlinks,pnode,psnode,pvert,
                  pbpoint,pelement,pface,psface,plink,pelink,pflink,pnflink,
                  pfnlink,pnlglink,nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,
                  plglglink,lglglinks,pflglink,flglinks,plgdata,lgdatas);

   bound_val(mg,U,D,Q);
   if (j == 3)
      Read_results(mg,vertexes,pvert);
   else
      stiff_matr(mg,TOP_GRID(mg),1.0,1.0,A,F,U,EF);
   test4(mg);
   L2test(TOP_GRID(mg),1,U,&err);
   H10errors(TOP_GRID(mg),U,&err,&err,&err);
   
   if (j == 2 || j == 3)
      Save_results(mg,vertexes,pvert);
   else{
      for (i=1; i<11; i++){
         printf("\n\n");
         printf("Refinement:\n\n");
         
         estimate1(mg);
         nonuniform_refinement(mg,&pnode,&pface,&pvert,&plink,&pnflink,&pfnlink,
                                                             &pflink,&pelement);
         size_of_data(mg,nodes,snodes,vertexes,bpoints,elements,faces,
                  sfaces,links,elinks,flinks,nflinks,fnlinks,pnode,psnode,pvert,
                  pbpoint,pelement,pface,psface,plink,pelink,pflink,pnflink,
                  pfnlink,pnlglink,nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,
                  plglglink,lglglinks,pflglink,flglinks,plgdata,lgdatas);
         bound_val(mg,U,D,Q);
         stiff_matr(mg,TOP_GRID(mg),1.0,1.0,A,F,U,EF); 
         test4(mg);
         L2test(TOP_GRID(mg),1,U,&err);
         H10errors(TOP_GRID(mg),U,&err,&err,&err);
      }
      if (j == 1)
         save_finest_triangulation(mg,vertexes,pvert);
      else{
         if (j == 5){
            read_matrix(mg,vertexes,pvert,D);
            test4(mg);
            L2test(TOP_GRID(mg),1,U,&err);
            H10errors(TOP_GRID(mg),U,&err,&err,&err); 
         }         
         compare_results(mg,vertexes,pvert,D,D,0);          
      } 
   }
}

#else

void tests(mg,nodes,snodes,p_nodes,vertexes,elements,faces,sfaces,p_faces,links,
           elinks,flinks,nflinks,fnlinks,pnode,psnode,pp_node,pvert,pelement,
           pface,psface,pp_face,plink,pelink,pflink,pnflink,pfnlink,pnlglink,
           nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,plglglink,lglglinks,
           pflglink,flglinks,plgdata,lgdatas)
MULTIGRID *mg; NODE *nodes, *pnode; NODE *snodes, *psnode; P_NODE *p_nodes, *pp_node; VERTEX *vertexes, *pvert; ELEMENT *elements, *pelement; FACE *faces, *pface; FACE *sfaces, *psface; P_FACE *p_faces, *pp_face; LINK *links, *plink; ELINK *elinks, *pelink; FNLINK *fnlinks, *pfnlink; NFLINK *nflinks, *pnflink; FLINK *flinks, *pflink; NLGLINK *nlglinks, *pnlglink; LGNLINK *lgnlinks, *plgnlink; LGFLINK *lgflinks, *plgflink; LGLGLINK *lglglinks, *plglglink; FLGLINK *flglinks, *pflglink; LGDATA *lgdatas, *plgdata;
{
   eprintf("No tests available.\n\n");
}

#endif
