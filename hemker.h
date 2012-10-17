/******************************************************************************/
/*                                                                            */
/*                     triangulation  for hemker test case                    */
/*                                                                            */
/******************************************************************************/

void stretch_and_shift(vertexes,pvert,stretch_x,stretch_y,shift_x,shift_y)
VERTEX *vertexes, *pvert;
FLOAT stretch_x, stretch_y, shift_x, shift_y;
{
   VERTEX *pv;

   pv = vertexes;
   while(pv < pvert){
      pv->x[0] = pv->x[0]*stretch_x + shift_x;
      pv->x[1] = pv->x[1]*stretch_y + shift_y;
      pv++;
   }
}

void mark_vertices_for_hemker(vertexes,pvert)           /*  in om  ...  0,    */
VERTEX *vertexes, *pvert;     /*  on boundary  ...  1, outside om  ...  4     */
{
   VERTEX *pv;
   FLOAT r2;

   pv = vertexes;
   while(pv < pvert){
      r2 = pv->x[0]*pv->x[0] + pv->x[1]*pv->x[1];
      if (fabs(r2-1.) < EPSA)
         pv->type = 1;
      else if (r2 < 1.)
         pv->type = 4;
      NDS(pv->topnode,U) = 20.;
      pv++;
   }
}

INT find_intersection(v1,v2)
VERTEX *v1, *v2;
{
   FLOAT a0, a1, b0, b1, a, b, c, d, s;
   INT i=0;

   if (v1->type + v2->type == 4){
      i = 1;
      a0 = v1->x[0];
      a1 = v1->x[1];
      b0 = v2->x[0] - v1->x[0];
      b1 = v2->x[1] - v1->x[1];
      a = b0*b0 + b1*b1;
      b = a0*b0 + a1*b1;
      c = a0*a0 + a1*a1 -1.;
      d = sqrt(b*b - a*c);
      s = -(b+d)/a;
      if (s < EPSA || s > 1.-EPSA)
         s = (-b+d)/a;
      if (s < EPSA || s > 1.-EPSA)
         eprintf("Error in find_intersection.\n");
      a0 += s*b0;
      a1 += s*b1;     /*  (a0,a1)  is the intersection with the boundary  */
      if (s < 0.5){   /*  (a0,a1) is nearer to v1 than to v2  */ 
         a = s*sqrt(a);
         if (a < NDS(v1->topnode,U)){
            NDS(v1->topnode,U) = a;
            NDS(v1->topnode,F) = v2->topnode->index;
            ND(v1->topnode,U,0) = a0;
            ND(v1->topnode,U,1) = a1;
         }
      }
      else{
         a = (1.-s)*sqrt(a);
         if (a < NDS(v2->topnode,U)){
            NDS(v2->topnode,U) = a;
            NDS(v2->topnode,F) = v1->topnode->index;
            ND(v2->topnode,U,0) = a0;
            ND(v2->topnode,U,1) = a1;
         }
      }
   }
   return(i);
}

INT compute_intersections(pgrid)
GRID *pgrid;
{
   VERTEX *v0, *v1, *v2;
   ELEMENT *pel;
   INT i=0;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      v0 = pel->n[0]->myvertex;
      v1 = pel->n[1]->myvertex;
      v2 = pel->n[2]->myvertex;
      if (find_intersection(v0,v1)) i = 1;
      if (find_intersection(v0,v2)) i = 1;
      if (find_intersection(v1,v2)) i = 1;
   }
   return(i);
}

void check_new_coordinates(pgrid)
GRID *pgrid;
{
   VERTEX *v0, *v1, *v2;
   ELEMENT *pel;
   FLOAT d0[DIM], d1[DIM], d2[DIM], l0, l1, l2, c;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      v0 = pel->n[0]->myvertex;
      v1 = pel->n[1]->myvertex;
      v2 = pel->n[2]->myvertex;
      if (NDS(v0->topnode,U) < 15. || NDS(v1->topnode,U) < 15. || 
          NDS(v2->topnode,U) < 15.){
         SUBTR(v0->x,v1->x,d2);
         SUBTR(v1->x,v2->x,d0);
         SUBTR(v2->x,v0->x,d1);
         l0 = DOT(d0,d0);
         l1 = DOT(d1,d1);
         l2 = DOT(d2,d2);
         if (l2 > l0 - EPSA && l2 > l1 - EPSA)
            c = (l0 + l1)/l2;
         else if (l1 > l0 - EPSA && l1 > l2 - EPSA)
            c = (l0 + l2)/l1;
         else
            c = (l1 + l2)/l0;
         printf("%f  %f  %f  %f\n",l0,l1,l2,c);
         if (c < 1.2){
            if (NDS(v0->topnode,U) < 15.) NDS(v0->topnode,U) = -1.;
            if (NDS(v1->topnode,U) < 15.) NDS(v1->topnode,U) = -1.;
            if (NDS(v2->topnode,U) < 15.) NDS(v2->topnode,U) = -1.;
         }
      }
   }
}
      
void print_new_coordinates(pgrid)
GRID *pgrid;
{
   VERTEX *pv, *v0, *v1, *v2;
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(pgrid); pel != NULL; pel = pel->succ){
      v0 = pel->n[0]->myvertex;
      v1 = pel->n[1]->myvertex;
      v2 = pel->n[2]->myvertex;
      if (NDS(v0->topnode,U) < 15. || NDS(v1->topnode,U) < 15. || 
          NDS(v2->topnode,U) < 15.){
          printf("----------------------------------------\n");
          if (NDS(v0->topnode,U) < 15.)
             printf("(%7.2f, %7.2f)  ->  (%7.2f, %7.2f)      %7.2f    %7.2f\n",
                     v0->x[0],v0->x[1],
                     ND(v0->topnode,U,0),ND(v0->topnode,U,1),
                     sqrt(ND(v0->topnode,U,0)*ND(v0->topnode,U,0)+
                          ND(v0->topnode,U,1)*ND(v0->topnode,U,1)),
                     NDS(v0->topnode,U));
          else
             printf("(%7.2f, %7.2f)\n",v0->x[0],v0->x[1]);
          if (NDS(v1->topnode,U) < 15.)
             printf("(%7.2f, %7.2f)  ->  (%7.2f, %7.2f)      %7.2f    %7.2f\n",
                     v1->x[0],v1->x[1],
                     ND(v1->topnode,U,0),ND(v1->topnode,U,1),
                     sqrt(ND(v1->topnode,U,0)*ND(v1->topnode,U,0)+
                          ND(v1->topnode,U,1)*ND(v1->topnode,U,1)),
                     NDS(v1->topnode,U));
          else
             printf("(%7.2f, %7.2f)\n",v1->x[0],v1->x[1]);
          if (NDS(v2->topnode,U) < 15.)
             printf("(%7.2f, %7.2f)  ->  (%7.2f, %7.2f)      %7.2f    %7.2f\n",
                     v2->x[0],v2->x[1],
                     ND(v2->topnode,U,0),ND(v2->topnode,U,1),
                     sqrt(ND(v2->topnode,U,0)*ND(v2->topnode,U,0)+
                          ND(v2->topnode,U,1)*ND(v2->topnode,U,1)),
                     NDS(v2->topnode,U));
          else
             printf("(%7.2f, %7.2f)\n",v2->x[0],v2->x[1]);
      }
   }
}

void shitf_vertices_to_boundary(pgrid,vertexes,pvert)
GRID *pgrid;
VERTEX *vertexes, *pvert;
{
   VERTEX *pv, *v0, *v1, *v2;
   LINK *pl;

   while (compute_intersections(pgrid)){
      pv = vertexes;
      while(pv < pvert){
         if (NDS(pv->topnode,U) < 15.){
            pv->x[0] = ND(pv->topnode,U,0);
            pv->x[1] = ND(pv->topnode,U,1);
            pv->type = 1;
            NDS(pv->topnode,U) = 20.;
            for (pl = TSTART(pv->topnode); pl != NULL; pl = NEXT(pl))
               NDS(NBNODE(pl),U) = 20.;
         }
         pv++;
      }
   }
}

void save_inner_elements_for_gnuplot(mg,name)
MULTIGRID *mg;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   FILE *fp;
   FLOAT *x0, *x1, *x2;
  
   fp = fopen(name,"w");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel) && NTYPE(pel->n[0]) != 4 &&
             NTYPE(pel->n[1]) != 4 &&  NTYPE(pel->n[2]) != 4){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            fprintf(fp,"%e %e\n",  x0[0],x0[1]);
            fprintf(fp,"%e %e\n",  x1[0],x1[1]);
            fprintf(fp,"%e %e\n",  x2[0],x2[1]);
            fprintf(fp,"%e %e\n\n",x0[0],x0[1]);
         }
   fclose(fp);
}

void save_inner_triangulation(mg,vertexes,pvert,name)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv;
   INT n=0;
   FILE *fp;

   fp = fopen(name,"w");
   pv = vertexes;
   while(pv < pvert){
      if (pv->type != 4)
         pv->topnode->index = n++;
      pv++;
   }
   fprintf(fp,"%ld \n",n);
   pv = vertexes;
   while(pv < pvert){
      if (pv->type != 4)
         fprintf(fp,"%e %e %i \n",pv->x[0],pv->x[1],pv->type);
      pv++;
   }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel) && NTYPE(pel->n[0]) != 4 &&
             NTYPE(pel->n[1]) != 4 &&  NTYPE(pel->n[2]) != 4)
            fprintf(fp,"%ld %ld %ld \n",pel->n[0]->index,pel->n[1]->index,
                                                         pel->n[2]->index);
   fprintf(fp,"%i %i %i ",-3,-3,-3);
   fclose(fp);
}
