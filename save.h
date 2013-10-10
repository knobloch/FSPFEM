#if N_DATA & SCALAR_NODE_DATA

void save_data_for_tecplot(tGrid,name)
GRID *tGrid;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   INT i=0;
   FILE *fp;

   for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
      pnode->index2 = ++i;
   fp = fopen(name,"w");
   fprintf(fp,"VARIABLES = \"X\", \"Y\", \"U\", \"U1\", \"U2\"\n");
   fprintf(fp,"ZONE N=%i",i);
   i = 0;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      i++;
   fprintf(fp,", E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",i);
   for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
      fprintf(fp,"%e %e %+e %+e %+e\n",
         pnode->myvertex->x[0],pnode->myvertex->x[1],NDS(pnode,U),
         u01(pnode->myvertex->x),u02(pnode->myvertex->x));
   fprintf(fp,"\n");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fprintf(fp,"%i %i %i\n",
         pel->n[0]->index2,pel->n[1]->index2,pel->n[2]->index2);
   fclose(fp);
}

#else

void save_data_for_tecplot(tGrid,name)
GRID *tGrid; char name[];
{ eprintf("Error: save_data_for_tecplot not available.\n"); }

#endif

#if DIM == 2 && ELEMENT_TYPE == SIMPLEX

INT is_in_triang(x,pel,eps,x_ref)
ELEMENT *pel;
FLOAT *x, eps, *x_ref;
{
   FLOAT *x0, *x1, *x2, b[DIM2][DIM2], 
         y0=x[0]+eps, z0=x[0]-eps, y1=x[1]+eps, z1=x[1]-eps;

   VERTICES_OF_ELEMENT(x0,x1,x2,pel);
   if ((y0 > x0[0] && y0 > x1[0] && y0 > x2[0]) ||
       (z0 < x0[0] && z0 < x1[0] && z0 < x2[0]) ||
       (y1 > x0[1] && y1 > x1[1] && y1 > x2[1]) ||
       (z1 < x0[1] && z1 < x1[1] && z1 < x2[1]))
      return(0.); 
   barycentric_coordinates(x0,x1,x2,b);
   if (IS_IN_SIMPLEX(b,x,eps)){
      x_ref[0] = LINV(b[1],x);
      x_ref[1] = LINV(b[2],x);
      return(1);
   }
   else
      return(0);
}

INT is_in_quadrilateral(x,pel,eps,x_ref)
ELEMENT *pel; FLOAT *x, eps, *x_ref;
{ eprintf("Error: is_in_quadrilateral not available.\n"); }

#elif DIM == 2 && ELEMENT_TYPE == CUBE

INT is_in_triang(x,pel,eps,x_ref)
ELEMENT *pel; FLOAT *x, eps, *x_ref;
{ eprintf("Error: is_in_triang not available.\n"); }

/* ONLY FOR PARALLELOGRAMS !!!!!!!!!!!!!!!!! */
INT is_in_quadrilateral(x,pel,eps,x_ref)
ELEMENT *pel;
FLOAT *x, eps, *x_ref;
{
   FLOAT *x0, *x1, *x2, *x3, b[DIM2][DIM2], bb[DIM2][DIM2];
   INT i;

   VERTICES_OF_4ELEMENT(x0,x1,x2,x3,pel);
   barycentric_coordinates(x0,x1,x3,b);
   i = IS_IN_SIMPLEX(b,x,eps);
   if (!i){
      barycentric_coordinates(x1,x2,x3,bb);
      i = IS_IN_SIMPLEX(bb,x,eps);
   }
   if (i){
      if (fabs(x0[0]-x1[0]+x2[0]-x3[0])+fabs(x0[1]-x1[1]+x2[1]-x3[1]) <
          1.e-10*(fabs(x0[0]-x1[0])+fabs(x1[0]-x2[0])
                 +fabs(x2[0]-x3[0])+fabs(x3[0]-x0[0])
                 +fabs(x0[1]-x1[1])+fabs(x1[1]-x2[1])
                 +fabs(x2[1]-x3[1])+fabs(x3[1]-x0[1]))){
         x_ref[0] = 2.*LINV(b[1],x) - 1.;
         x_ref[1] = 2.*LINV(b[2],x) - 1.;
      }
      else
         eprintf("is_in_quadrilateral implemented for parallelograms only.\n");
   }
   return(i);
}

#else

INT is_in_triang(x,pel,eps,x_ref)
ELEMENT *pel; FLOAT *x, eps, *x_ref;
{ eprintf("Error: is_in_triang not available.\n"); }

INT is_in_quadrilateral(x,pel,eps,x_ref)
ELEMENT *pel; FLOAT *x, eps, *x_ref;
{ eprintf("Error: is_in_quadrilateral not available.\n"); }

#endif

INT is_in_element(x,pel,eps,x_ref)
ELEMENT *pel;
FLOAT *x, eps, *x_ref;
{
   if (DIM == 2 && ELEMENT_TYPE == SIMPLEX)
      return(is_in_triang(x,pel,eps,x_ref));
   else if (DIM == 2 && ELEMENT_TYPE == CUBE)
      return(is_in_quadrilateral(x,pel,eps,x_ref));
   else
      eprintf("Error: is_in_element not available.\n");
}

INT intersection(x,v,y,w,a,b,u)  /* intersection u of x+a*v and y+b*w */
DOUBLE x[DIM], v[DIM], y[DIM], w[DIM], *a, *b, u[DIM]; /* v, w, x-y != 0 */
{
   DOUBLE z[DIM], d;

   SUBTR(y,x,z);  /* z := y - x */
   d = v[0]*w[1] - v[1]*w[0];
   if (fabs(d) < 1.e-20*(fabs(v[0])+fabs(v[1]))*(fabs(w[0])+fabs(w[1]))){
      /* v, w are parallel */
      if (fabs(z[0]*w[1] - z[1]*w[0]) < 
          1.e-20*(fabs(z[0])+fabs(z[1]))*(fabs(w[0])+fabs(w[1])))
         return(-1);   /* the lines coincide */ 
      else
         return(0);    /* no intersection */
   }
   else{
      *a = (z[0]*w[1] - z[1]*w[0])/d;
      *b = (z[0]*v[1] - z[1]*v[0])/d;
      u[0] = x[0] + (*a)*v[0];
      u[1] = x[1] + (*a)*v[1];
      return(1);   /* one intersection */
   }
}

/* intersection u of x+a*v and edge x0, x1 defining line x0 + b*(x1-x0) */
INT intersection_with_edge(x,v,x0,x1,a,b,u)
DOUBLE x[DIM], v[DIM], x0[DIM], x1[DIM], *a, *b, u[DIM]; 
{                                                     /* v, x-x0, x1-x0 != 0 */
   DOUBLE w[DIM];
   INT i;
 
   SUBTR(x1,x0,w);  /* w := x1 - x0 */
   i = intersection(x,v,x0,w,a,b,u);
   if (i == 1 && (*b < -1.e-20 || *b-1. > 1.e-20))
      i = 0;
   return(i);
}

INT is_between_points(x,y,z)  /* z is on the line through x, y between x, y */
DOUBLE x[DIM], y[DIM], z[DIM];  /* x != y */
{
   DOUBLE u[DIM], v[DIM], a;

   SUBTR(y,x,u);  /* u := y - x */
   SUBTR(z,x,v);  /* v := z - x */
   if (fabs(v[0])+fabs(v[1]) < 
       1.e-20*(fabs(x[0])+fabs(x[1])+fabs(z[0])+fabs(z[1])))
      return(1);
   else if (fabs(u[0]*v[1] - u[1]*v[0]) < 
            1.e-20*(fabs(u[0])+fabs(u[1]))*(fabs(v[0])+fabs(v[1]))){
      if (fabs(u[0]) > fabs(u[1]))
         a = v[0]/u[0];
      else
         a = v[1]/u[1];
      if (a > 0. && a-1. < 1.e-20)
         return(1);
      else
         return(0);
   }
   else
      return(0);
}

ELEMENT *element_containing_boundary_point(tGrid,z)
GRID *tGrid;
DOUBLE z[DIM];
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      if ((IS_BN(pel->n[0]) || IS_BN(pel->n[1]) || IS_BN(pel->n[2])) &&
          (is_between_points(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,z) 
        || is_between_points(pel->n[1]->myvertex->x,pel->n[2]->myvertex->x,z) 
        || is_between_points(pel->n[2]->myvertex->x,pel->n[0]->myvertex->x,z)))
         return(pel);
   eprintf("Error in element_containing_boundary_point.\n");
   return(NULL);
}

FLOAT svalue_at_point(tGrid,x,u,space)
GRID *tGrid;
FLOAT *x;
INT u, space;
{
   ELEMENT *pel;
   FLOAT x_ref[DIM], y[DIM], eps=0., val=0.;
   INT i=1;
   REF_MAPPING ref_map;

   for (pel = FIRSTELEMENT(tGrid); pel && i; pel = pel->succ)
      if (is_in_element(x,pel,eps,x_ref)){
         val = sfunc_value_ref(pel,x_ref,u,space);
         i = 0;
      }
   if (i)
      eprintf("Error in svalue_at_point.\n"); 
   return(val);
}

FLOAT gsvalue_at_point(tGrid,x,u,rmtype,finite_el)
GRID *tGrid;
FLOAT *x;
INT u, rmtype;
FINITE_ELEMENT finite_el;
{
   ELEMENT *pel;
   FLOAT x_ref[DIM], y[DIM], eps=0., val=0.;
   INT i=1;
   REF_MAPPING ref_map;

   for (pel = FIRSTELEMENT(tGrid); pel && i; pel = pel->succ)
      if (is_in_element(x,pel,eps,x_ref)){
/*
         reference_mapping(pel,rmtype,&ref_map);
         ref_map_point(x_ref,y,&ref_map);
         if (fabs(x[0]-y[0]) < 1.e-12 && fabs(x[1]-y[1]) < 1.e-12)
            printf("POINT O.K.\n");
         else
            printf("WRONG POINT.\n");
*/
         val = gfunc_value_ref(tGrid,pel,u,x_ref,finite_el);
         i = 0;
      }
   if (i)
      eprintf("Error in gsvalue_at_point.\n"); 
   return(val);
}

#if (N_DATA & SCALAR_NODE_DATA) && (DIM == 2)

void sn_grid_data(theGrid,xmin,xmax,ymin,ymax,nx,ny,u,name)
GRID *theGrid;
FLOAT xmin, xmax, ymin, ymax;
INT nx, ny, u;
char name[];
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FLOAT *x0, *x1, *x2, l0, l1, l2, b[DIM2][DIM2], x[DIM], eps = -1.e-10;
   FLOAT dx = (xmax-xmin)/nx, dy = (ymax-ymin)/ny;
   INT i, j, k;
   FILE *fp;
   
   fp = fopen(name,"w");
   for (i = 0; i <= nx; i++){
      x[0] = xmin + i*dx;
      for (j = 0; j <= ny; j++){
         x[1] = ymin + j*dy;
         k = 1;
         for (pel = FIRSTELEMENT(theGrid); k && pel; pel = pel->succ){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            barycentric_coordinates(x0,x1,x2,b); /* DOT(b[i],x)+b[i][2] is the
                             barycentric coordinate of x with respect to x_i */
            l0 = DOT(b[0],x)+b[0][2];
            l1 = DOT(b[1],x)+b[1][2];
            l2 = DOT(b[2],x)+b[2][2];
            if(l0 >= eps && l1 >= eps && l2 >= eps){
               NODES_OF_ELEMENT(n0,n1,n2,pel);
               if (fabs(x0[0]-x[0]) < 1.e-14 && fabs(x0[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],NDS(n0,u));
               else if (fabs(x1[0]-x[0]) < 1.e-14 && fabs(x1[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],NDS(n1,u));
               else if (fabs(x2[0]-x[0]) < 1.e-14 && fabs(x2[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],NDS(n2,u));
               else
                  fprintf(fp,"%e %e %e\n",x[0],x[1],
                                   NDS(n0,u)*l0 + NDS(n1,u)*l1 + NDS(n2,u)*l2);
               k = 0;
            }
         }
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

#else

void sn_grid_data(theGrid,xmin,xmax,ymin,ymax,nx,ny,u,name)
GRID *theGrid; FLOAT xmin, xmax, ymin, ymax; INT nx, ny, u; char name[];
{  eprintf("Error: sn_grid_data not available.\n");  }

#endif



#if (E_DATA & SCALAR_ELEMENT_DATA) && (DIM == 2)

void sn_grid_data2(mg,xmin,xmax,ymin,ymax,nx,ny,u,name)
MULTIGRID *mg;
FLOAT xmin, xmax, ymin, ymax;
INT nx, ny, u;
char name[];
{
   GRID *theGrid=TOP_GRID(mg);
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FLOAT *x0, *x1, *x2, l0, l1, l2, b[DIM2][DIM2], x[DIM], eps = -1.e-10;
   FLOAT dx = (xmax-xmin)/nx, dy = (ymax-ymin)/ny;
   INT i, j, k;
   FILE *fp;
   
   fp = fopen(name,"w");
   for (i = 0; i <= nx; i++){
      x[0] = xmin + i*dx;
      for (j = 0; j <= ny; j++){
         x[1] = ymin + j*dy;
         k = 1;
         for (pel = FIRSTELEMENT(theGrid); k && pel; pel = pel->succ){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            barycentric_coordinates(x0,x1,x2,b); /* DOT(b[i],x)+b[i][2] is the
                             barycentric coordinate of x with respect to x_i */
            l0 = DOT(b[0],x)+b[0][2];
            l1 = DOT(b[1],x)+b[1][2];
            l2 = DOT(b[2],x)+b[2][2];
            if(l0 >= eps && l1 >= eps && l2 >= eps){
               fprintf(fp,"%e %e %e\n",x[0],x[1],ED(pel,u));
               k = 0;
            }
         }
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

#else

void sn_grid_data2(mg,xmin,xmax,ymin,ymax,nx,ny,u,name)
MULTIGRID *mg; FLOAT xmin, xmax, ymin, ymax; INT nx, ny, u; char name[];
{  eprintf("Error: sn_grid_data2 not available.\n");  }

#endif

#if (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (DIM == 2)

void sn_grid_data3(mg,xmin,xmax,ymin,ymax,nx,ny,u,name)
MULTIGRID *mg;
FLOAT xmin, xmax, ymin, ymax;
INT nx, ny, u;
char name[];
{
   GRID *theGrid=TOP_GRID(mg);
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FLOAT *x0, *x1, *x2, l0, l1, l2, b[DIM2][DIM2], x[DIM], eps = -1.e-10;
   FLOAT dx = (xmax-xmin)/nx, dy = (ymax-ymin)/ny;
   INT i, j, k;
   FILE *fp;
   
   fp = fopen(name,"w");
   for (i = 0; i <= nx; i++){
      x[0] = xmin + i*dx;
      for (j = 0; j <= ny; j++){
         x[1] = ymin + j*dy;
         k = 1;
         for (pel = FIRSTELEMENT(theGrid); k && pel; pel = pel->succ){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            barycentric_coordinates(x0,x1,x2,b);
            l0 = DOT(b[0],x)+b[0][2];
            l1 = DOT(b[1],x)+b[1][2];
            l2 = DOT(b[2],x)+b[2][2];
            if(l0 >= eps && l1 >= eps && l2 >= eps){
               NODES_OF_ELEMENT(n0,n1,n2,pel);
               if (fabs(x0[0]-x[0]) < 1.e-14 && fabs(x0[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],EDSN(pel,u,0));
               else if (fabs(x1[0]-x[0]) < 1.e-14 && fabs(x1[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],EDSN(pel,u,1));
               else if (fabs(x2[0]-x[0]) < 1.e-14 && fabs(x2[1]-x[1]) < 1.e-14)
                  fprintf(fp,"%e %e %e\n",x[0],x[1],EDSN(pel,u,2));
               else
                  fprintf(fp,"%e %e %e\n",x[0],x[1],l0*EDSN(pel,u,0)
                                                   +l1*EDSN(pel,u,1)
                                                   +l2*EDSN(pel,u,2));
               k = 0;
            }
         }
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

#else

void sn_grid_data3(mg,xmin,xmax,ymin,ymax,nx,ny,u,name)
MULTIGRID *mg; FLOAT xmin, xmax, ymin, ymax; INT nx, ny, u; char name[];
{  eprintf("Error: sn_grid_data3 not available.\n");  }

#endif







#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void gsn_grid_data(theGrid,xmin,xmax,ymin,ymax,nx,ny,u,name)
GRID *theGrid;
FLOAT xmin, xmax, ymin, ymax;
INT nx, ny, u;
char name[];
{
   NODE *pnode;
   FLOAT x[DIM];
   FLOAT dx = (xmax-xmin)/nx, dy = (ymax-ymin)/ny;
   INT i, j, k;
   FILE *fp;
   
   fp = fopen(name,"w");
   for (i = 0; i <= nx; i++){
      x[0] = xmin + i*dx;
      for (j = 0; j <= ny; j++){
         x[1] = ymin + j*dy;
         k = 1;
         pnode = FIRSTN(theGrid);
         while (pnode && k)
            if (fabs(pnode->myvertex->x[0]-x[0]) > 1.e-10 ||
                fabs(pnode->myvertex->x[1]-x[1]) > 1.e-10)
               pnode = pnode->succ;
            else
               k = 0;
         if (k)
            eprintf("Error: node not found in gsn_grid_data.\n");
         else
            fprintf(fp,"%e %e %e\n",x[0],x[1],NDMV(pnode,u,0));
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

#else

void gsn_grid_data(theGrid,xmin,xmax,ymin,ymax,nx,ny,u,name)
GRID *theGrid; FLOAT xmin, xmax, ymin, ymax; INT nx, ny, u; char name[];
{  eprintf("Error: gsn_grid_data not available.\n");  }

#endif

/******************************************************************************/
/*                                                                            */
/*                         saving and reading of data                         */
/*                                                                            */
/******************************************************************************/

char *cfn3(i)
INT i;
{
        static char pom[18];

//        memset(&pom[0],0,sizeof(char)*18);
        sprintf(&pom[0],"profile%.4i.gnu",i);
        return &pom[0];
}

void make_string_from_number(n,s)
int n;
char s[];
{
   int i;

   if (0 <= n && n < 1000){
      i = n/100;
      s[0] = '0' + i;
      n -= i*100;
      i = n/10;
      s[1] = '0' + i;
      n -= i*10;
      s[2] = '0' + n;
      s[3] = '\0';
   }
   else
      printf("Error: bad value of n in make_string_from_number.\n");
} 

void connect_strings(s1,s2,s)
char s1[], s2[], s[];
{
   int i=-1, j=0;

   while (s1[++i] != '\0') s[i] = s1[i];
   while (s2[j] != '\0') s[i++] = s2[j++];
   s[i] = '\0';
}

void open_pipe_for_gnuplot()
{
   if (PIPE_FOR_GNUPLOT_OPENED == YES)
      printf("Pipe for gnuplot have been opened already.\n");
   else{
      if (access(FIFO_NAME,F_OK) == -1) 
         if (mkfifo(FIFO_NAME,07777)){
            fprintf(stderr,"Error: Cannot create the pipe  %s\n",FIFO_NAME);
            exit(EXIT_FAILURE);
         }
      if ((PIPE_FD_FOR_GNUPLOT=open(FIFO_NAME, O_WRONLY)) == -1)
         exit(EXIT_FAILURE);
      else
         PIPE_FOR_GNUPLOT_OPENED = YES;
   }
}

void close_pipe_for_gnuplot()
{
   if (PIPE_FOR_GNUPLOT_OPENED == YES){
      (void)close(PIPE_FD_FOR_GNUPLOT);
      PIPE_FOR_GNUPLOT_OPENED = NO;
   }
}

int send_command(int pipe_fd, const char *command)
{
   if (write(pipe_fd,command,strlen(command)) == -1){
      fprintf(stderr,"Error when writing in the pipe.\n");
      return(EXIT_FAILURE);
   }
}

void send_command_to_gnuplot(const char *command)
{
   if (PIPE_FOR_GNUPLOT_OPENED == NO)
      open_pipe_for_gnuplot();
   send_command(PIPE_FD_FOR_GNUPLOT, command);
}

void plot_with_gnuplot(filename,options)
char *filename, *options;
{
   char s[1000], command[1000];

   if (PIPE_FOR_GNUPLOT_OPENED == NO)
      open_pipe_for_gnuplot();
   connect_strings("\nplot '",filename,s);
   connect_strings(s,"' ",command);
   connect_strings(command,options,s);
   connect_strings(s,"\n",command);
   send_command(PIPE_FD_FOR_GNUPLOT, command);
   send_command(PIPE_FD_FOR_GNUPLOT,"\nreplot\n");
}

/*
   send_command_to_gnuplot("\n set xrange [8:15]; set yrange [10:12]\n");
   plot_with_gnuplot("datafile","");
   sleep(2);
   plot_with_gnuplot("datafile","with lines");
   sleep(2);
   plot_with_gnuplot("datafile1","with lines");
   sleep(2);
*/

void splot_with_gnuplot(filename,options)
char *filename, *options;
{
   char s[1000], command[1000];

   if (PIPE_FOR_GNUPLOT_OPENED == NO)
      open_pipe_for_gnuplot();
   connect_strings("\nsplot '",filename,s);
   connect_strings(s,"' ",command);
   connect_strings(command,options,s);
   connect_strings(s,"\n",command);
   send_command(PIPE_FD_FOR_GNUPLOT, command);
   send_command(PIPE_FD_FOR_GNUPLOT,"\nreplot\n");
}

#if DIM == 3

void save_finest_triangulation(mg,vertexes,pvert)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv;
   INT sv;  
   FILE *fp;
   
   fp = fopen("data","w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%ld \n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e %e %i \n",pv->x[0],pv->x[1],pv->x[2],pv->type);
      pv++;
   }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld %ld \n",
              ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[3]->myvertex)  - (long)(vertexes))/sv);
   fprintf(fp,"%i %i %i %i ",-3,-3,-3,-3);
   fclose(fp);
}

void save_triangulation_for_gnu2gr(mg)
MULTIGRID *mg;
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv1, *pv2, *pv3, *pv4;
   FILE *fp;
   
   fp = fopen("Tdata","w");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            pv1 = pel->n[0]->myvertex;
            pv2 = pel->n[1]->myvertex;
            pv3 = pel->n[2]->myvertex;
            pv4 = pel->n[3]->myvertex;
            fprintf(fp,"%e %e %e\n",pv1->x[0],pv1->x[1],pv1->x[2]);
            fprintf(fp,"%e %e %e\n",pv2->x[0],pv2->x[1],pv2->x[2]);
            fprintf(fp,"%e %e %e\n",pv3->x[0],pv3->x[1],pv3->x[2]);
            fprintf(fp,"%e %e %e\n\n",pv4->x[0],pv4->x[1],pv4->x[2]);
         }
   fclose(fp);
}

#if F_DATA & SCALAR_FACE_DATA

void save_triangulation_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   FACE *pface;
   VERTEX *pv;
   INT i=0, j, sv;  
   FILE *fp;
    
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pface = FIRSTF(theGrid); pface != NULL; pface = pface->succ)
         FD(pface,p) = FD(pface,q) = -1.;
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            pel->index2 = i++;
            for (j=0; j < 4; j++)
               if (FD(pel->f[j],p) < -0.5)
                  FD(pel->f[j],p) = pel->index2;
               else
                  FD(pel->f[j],q) = pel->index2;
          } 
   fp = fopen(name,"w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%i\n",i);
   fprintf(fp,"%ld\n\n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e %e\n",pv->x[0],pv->x[1],pv->x[2]);
      pv++;
   }
   fprintf(fp,"\n");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld %ld \n",
               ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
               ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
               ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv,
               ((long)(pel->n[3]->myvertex)  - (long)(vertexes))/sv);
   fprintf(fp,"\n");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            for (j=0; j < 4; j++){
               if ((int)(FD(pel->f[j],q)) == pel->index2)
                  i = (int)(FD(pel->f[j],p));
               else
                  i = (int)(FD(pel->f[j],q));
               fprintf(fp,"%i ",i);
            }
            fprintf(fp,"\n");
         }
   fclose(fp);
}

#else

void save_triangulation_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg; VERTEX *vertexes, *pvert; INT p, q; char name[];
{  eprintf("Error: save_triangulation_for_grape not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void save_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
{
   VERTEX *pv;
   FILE *fp;

   save_triangulation_for_grape(mg,vertexes,pvert,p,q,name);
   fp = fopen(name,"a");
   fprintf(fp,"\n3\n");   /* number of values in each row */
   pv = vertexes;
   while(pv < pvert){
/*      fprintf(fp,"%e %e %e %e\n",ND(pv->topnode,U,0),ND(pv->topnode,U,1),
                                    ND(pv->topnode,U,2),NDS(pv->topnode,U)); */
      fprintf(fp,"%e %e %e\n",ND(pv->topnode,U,0),ND(pv->topnode,U,1),
                                                       ND(pv->topnode,U,2));
      pv++;
   }
   fprintf(fp,"\n");
   fclose(fp);
}

#else

void save_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg; VERTEX *vertexes, *pvert; INT p, q; char name[];
{  eprintf("Error: save_for_grape not available.\n");  }

#endif

void save_exact_solution_for_grape(mg,vertexes,pvert,p,q,name,u1,u2,u3)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
FLOAT (*u1)(), (*u2)(), (*u3)();
{
   VERTEX *pv;
   FILE *fp;

   save_triangulation_for_grape(mg,vertexes,pvert,p,q,name);
   fp = fopen(name,"a");
   fprintf(fp,"\n3\n");   /* number of values in each row */
   pv = vertexes;
   while(pv < pvert){
      fprintf(fp,"%e %e %e\n",u1(pv->x),u2(pv->x),u3(pv->x));
      pv++;
   }
   fprintf(fp,"\n");
   fclose(fp);
}

void save_triangulation_for_gnuplot(mg,name)
MULTIGRID *mg;
char name[];
{
   eprintf("Error: save_triangulation_for_gnuplot not available.\n");
}

void save_boundary_for_gnuplot(mg,name)
MULTIGRID *mg;
char name[];
{
   eprintf("Error: save_boundary_for_gnuplot not available.\n");
}

void save_arrows_for_gnuplot(mg,u,name)
MULTIGRID *mg;
INT u;
char name[];
{
   eprintf("Error: save_arrows_for_gnuplot not available.\n");
}

void save_exact_arrows_for_gnuplot(mg,name,u1,u2)
MULTIGRID *mg; char name[]; FLOAT (*u1)(), (*u2)();
{
   eprintf("Error: save_exact_arrows_for_gnuplot not available.\n");
}

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void save_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;
  
   fp = fopen(name,"w");
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fprintf(fp,"%e %e %e\n",ND(pnode,u,0),ND(pnode,u,1),ND(pnode,u,2));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e\n",FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,ue));
   fclose(fp);
}

void read_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"r");
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fscanf(fp,"%lf %lf %lf\n",&ND(pnode,u,0),&ND(pnode,u,1),&ND(pnode,u,2));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf\n",&FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf\n",&ED(pel,ue));
   fclose(fp);
}

#else

void save_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   eprintf("Error: save_res not available.\n");
}

void read_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   eprintf("Error: read_res not available.\n");
}

#endif

#else  /* if DIM == 2 */

#if ELEMENT_TYPE == SIMPLEX

void save_finest_triangulation(mg,vertexes,pvert)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv;
   INT sv;  
   FILE *fp;
   
   fp = fopen("data","w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%ld \n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e %i \n",pv->x[0],pv->x[1],pv->type);
      pv++;
   }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld \n",
              ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv);
   fprintf(fp,"%i %i %i ",-3,-3,-3);
   fclose(fp);
}

#else /* if ELEMENT_TYPE == CUBE */

void save_finest_triangulation(mg,vertexes,pvert)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv;
   INT sv;  
   FILE *fp;
   
   fp = fopen("data","w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%ld \n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e %i \n",pv->x[0],pv->x[1],pv->type);
      pv++;
   }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld %ld \n",
              ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[3]->myvertex)  - (long)(vertexes))/sv);
   fprintf(fp,"%i %i %i %i ",-3,-3,-3,-3);
   fclose(fp);
}

#endif

#if (ELEMENT_TYPE == SIMPLEX) && (N_DATA & NODE_ITYPE) && (E_DATA & ELEM_ITYPE)

void save_finest_triangulation_with_itype(mg,vertexes,pvert)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
{
   ELEMENT *pel;
   GRID *theGrid;
   VERTEX *pv;
   INT sv;  
   FILE *fp;
   
   fp = fopen("data","w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%ld \n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e %i %i\n",pv->x[0],pv->x[1],pv->type,ITYPE(pv->topnode));
      pv++;
   }
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld %i\n",
              ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
              ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv,
              ITYPE(pel));
   fprintf(fp,"%i %i %i %i",-3,-3,-3,-3);
   fclose(fp);
}

#else

void save_finest_triangulation_with_itype(mg,vertexes,pvert)
MULTIGRID *mg; VERTEX *vertexes, *pvert;
{  eprintf("Error: save_finest_triangulation_with_itype not available.\n");  }

#endif

void save_triangulation_for_gnu2gr(mg)
MULTIGRID *mg;
{
   eprintf("Error: save_triangulation_for_gnu2gr not available.\n");
}
   
#if (ELEMENT_TYPE == SIMPLEX) && (F_DATA & SCALAR_FACE_DATA)

void save_triangulation_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   FACE *pface;
   VERTEX *pv;
   INT i=0, j, sv;  
   FILE *fp;
    
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pface = FIRSTF(theGrid); pface != NULL; pface = pface->succ)
         FD(pface,p) = FD(pface,q) = -1.;
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            pel->index2 = i++;
            for (j=0; j < SIDES; j++)
               if (FD(pel->f[j],p) < -0.5)
                  FD(pel->f[j],p) = pel->index2;
               else
                  FD(pel->f[j],q) = pel->index2;
          } 
   fp = fopen(name,"w");
   sv = sizeof(VERTEX);
   pv = vertexes;
   fprintf(fp,"%i\n",i);
   fprintf(fp,"%ld\n\n",((long)(pvert) - (long)(vertexes))/sv);
   while(pv < pvert){
      fprintf(fp,"%e %e\n",pv->x[0],pv->x[1]);
      pv++;
   }
   fprintf(fp,"\n");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel))
            fprintf(fp,"%ld %ld %ld\n",
               ((long)(pel->n[0]->myvertex)  - (long)(vertexes))/sv,
               ((long)(pel->n[1]->myvertex)  - (long)(vertexes))/sv,
               ((long)(pel->n[2]->myvertex)  - (long)(vertexes))/sv);
   fprintf(fp,"\n");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            for (j=0; j < SIDES; j++){
               if ((int)(FD(pel->f[j],q)) == pel->index2)
                  i = (int)(FD(pel->f[j],p));
               else
                  i = (int)(FD(pel->f[j],q));
               fprintf(fp,"%i ",i);
            }
            fprintf(fp,"\n");
         }
   fclose(fp);
}

#else

void save_triangulation_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg; VERTEX *vertexes, *pvert; INT p, q; char name[];
{  eprintf("Error: save_triangulation_for_grape not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void save_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
{
   VERTEX *pv;
   FILE *fp;

   save_triangulation_for_grape(mg,vertexes,pvert,p,q,name);
   fp = fopen(name,"a");
   fprintf(fp,"\n2\n");   /* number of values in each row */
   pv = vertexes;
   while(pv < pvert){
      fprintf(fp,"%e %e\n",ND(pv->topnode,U,0),ND(pv->topnode,U,1));
      pv++;
   }
   fprintf(fp,"\n");
   fclose(fp);
}

#else

void save_for_grape(mg,vertexes,pvert,p,q,name)
MULTIGRID *mg; VERTEX *vertexes, *pvert; INT p, q; char name[];
{  eprintf("Error: save_for_grape not available.\n");  }

#endif

void save_exact_solution_for_grape(mg,vertexes,pvert,p,q,name,u1,u2,u3)
MULTIGRID *mg;
VERTEX *vertexes, *pvert;
INT p, q;
char name[];
FLOAT (*u1)(), (*u2)(), (*u3)();
{
   VERTEX *pv;
   FILE *fp;

   save_triangulation_for_grape(mg,vertexes,pvert,p,q,name);
   fp = fopen(name,"a");
   fprintf(fp,"\n2\n");   /* number of values in each row */
   pv = vertexes;
   while(pv < pvert){
      fprintf(fp,"%e %e\n",u1(pv->x),u2(pv->x));
      pv++;
   }
   fprintf(fp,"\n");
   fclose(fp);
}

#if ELEMENT_TYPE == SIMPLEX

void save_triangulation_for_gnuplot(mg,name)
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
         if (IS_TOP_ELEMENT(pel)){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            fprintf(fp,"%e %e\n",  x0[0],x0[1]);
            fprintf(fp,"%e %e\n",  x1[0],x1[1]);
            fprintf(fp,"%e %e\n",  x2[0],x2[1]);
            fprintf(fp,"%e %e\n\n",x0[0],x0[1]);
         }
   fclose(fp);
}

#elif ELEMENT_TYPE == CUBE

void save_triangulation_for_gnuplot(mg,name)
MULTIGRID *mg;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   FILE *fp;
   FLOAT *x0, *x1, *x2, *x3;
   
   fp = fopen(name,"w");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            VERTICES_OF_4ELEMENT(x0,x1,x2,x3,pel);
            fprintf(fp,"%e %e\n",  x0[0],x0[1]);
            fprintf(fp,"%e %e\n",  x1[0],x1[1]);
            fprintf(fp,"%e %e\n",  x2[0],x2[1]);
            fprintf(fp,"%e %e\n",  x3[0],x3[1]);
            fprintf(fp,"%e %e\n\n",x0[0],x0[1]);
         }
   fclose(fp);
}

#else

void save_triangulation_for_gnuplot(mg,name)
MULTIGRID *mg; char name[];
{  eprintf("Error: save_triangulation_for_gnuplot not available.\n");  }

#endif

void save_boundary_for_gnuplot(mg,name)
MULTIGRID *mg;
char name[];
{
   ELEMENT *pel;
   FILE *fp;
   FLOAT *x0, *x1, *x2;
   
   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel != NULL; pel = pel->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      if (IS_BF(pel->f[2])) 
         fprintf(fp,"%e %e\n%e %e\n\n",x0[0],x0[1],x1[0],x1[1]);
      if (IS_BF(pel->f[0]))
         fprintf(fp,"%e %e\n%e %e\n\n",x1[0],x1[1],x2[0],x2[1]);
      if (IS_BF(pel->f[1]))
         fprintf(fp,"%e %e\n%e %e\n\n",x2[0],x2[1],x0[0],x0[1]);
   }
   fclose(fp);
}

void new_coord(ax,ay,bx,by,ymax)
FLOAT *ax, *ay, *bx, *by, ymax;
{
   FLOAT a;

   a = ((*ay)-ymax)/((*ay) - (*by));
   *ax = (1.-a)*(*ax) + a*(*bx);
   *ay = ymax;
}

void draw_line(ax,ay,bx,by,xmin,xmax,ymin,ymax,fp)
FLOAT ax,ay,bx,by,xmin,xmax,ymin,ymax;
FILE *fp;
{
   if (xmin <= ax && ax <= xmax && ymin <= ay && ay <= ymax &&
       xmin <= bx && bx <= xmax && ymin <= by && by <= ymax)
      fprintf(fp,"%e %e\n%e %e\n\n",ax,ay,bx,by);
   else if (!( (ay >= ymax && by >= ymax) || (ay <= ymin && by <= ymin) || 
               (ax >= xmax && bx >= ymax) || (ax <= xmin && bx <= xmin) )){
      if (ay > ymax) new_coord(&ax,&ay,&bx,&by,ymax);
      else if (by > ymax) new_coord(&bx,&by,&ax,&ay,ymax);
      if (ay < ymin) new_coord(&ax,&ay,&bx,&by,ymin);
      else if (by < ymin) new_coord(&bx,&by,&ax,&ay,ymin);
      if (!( (ax >= xmax && bx >= ymax) || (ax <= xmin && bx <= xmin) )){
         if (ax > xmax) new_coord(&ay,&ax,&by,&bx,xmax);
         else if (bx > xmax) new_coord(&by,&bx,&ay,&ax,xmax);
         if (ax < xmin) new_coord(&ay,&ax,&by,&bx,xmin);
         else if (bx < xmin) new_coord(&by,&bx,&ay,&ax,xmin);
         fprintf(fp,"%e %e\n%e %e\n\n",ax,ay,bx,by);
      }
   }
}

void draw_arrow(x0,x1,u0,u1,scale,xmin,xmax,ymin,ymax,fp)
FLOAT x0,x1,u0,u1,scale,xmin,xmax,ymin,ymax;
FILE *fp;
{
   FLOAT y0, y1, r0, r1, s0, s1;

   u0 *= scale;
   u1 *= scale;
   y0 = x0 + u0;
   y1 = x1 + u1;
   u0 /= 4.;
   u1 /= 4.;
   r0 = y0 - u0 - u1;
   r1 = y1 + u0 - u1;
   s0 = y0 - u0 + u1;
   s1 = y1 - u0 - u1;
   draw_line(x0,x1,y0,y1,xmin,xmax,ymin,ymax,fp);
   draw_line(y0,y1,r0,r1,xmin,xmax,ymin,ymax,fp);
   draw_line(y0,y1,s0,s1,xmin,xmax,ymin,ymax,fp);
/*   fprintf(fp,"%e %e\n%e %e\n\n",x0,x1,y0,y1);
   fprintf(fp,"%e %e\n%e %e\n\n",y0,y1,r0,r1);
   fprintf(fp,"%e %e\n%e %e\n\n",y0,y1,s0,s1); */
}

#define NPOINTS       17
#define NNPOINTS     289  /*  NPOINTS*NPOINTS  */

void points_for_visualization(mg,xmin,xmax,ymin,ymax,x,f,n)
MULTIGRID *mg;
FLOAT *xmin, *xmax, *ymin, *ymax, x[NNPOINTS][2];
INT f[NNPOINTS], *n;
{
   NODE *pnode;
   FLOAT xp, yp, d;
   
   *xmin = *xmax = FIRSTN(FIRSTGRID(mg))->myvertex->x[0];
   *ymin = *ymax = FIRSTN(FIRSTGRID(mg))->myvertex->x[1];
   for (pnode = FIRSTN(FIRSTGRID(mg)); pnode != NULL; pnode=pnode->succ){
      if (pnode->myvertex->x[0] > *xmax) 
         *xmax = pnode->myvertex->x[0];
      else if (pnode->myvertex->x[0] < *xmin)
         *xmin = pnode->myvertex->x[0];
      if (pnode->myvertex->x[1] > *ymax) 
         *ymax = pnode->myvertex->x[1];
      else if (pnode->myvertex->x[1] < *ymin)
         *ymin = pnode->myvertex->x[1];
   }
   printf("xmin = %6.2e, xmax = %6.2e, ymin = %6.2e, ymax = %6.2e\n",
           *xmin,*xmax,*ymin,*ymax);
   *n = 0;
   d = MAX(*xmax - *xmin,*ymax - *ymin)/(NPOINTS-1);
   xp = *xmin;
   while (xp <= *xmax){
      yp = *ymin;
      while (yp <= *ymax){
         x[*n][0] = xp;
         x[*n][1] = yp;
         f[*n] = 1;
         (*n)++;
         yp += d;
      }
      xp += d;
   }
   printf("number of points = %i\n",*n);
}

#if N_DATA & VECTOR_NODE_DATA

void find_elements(mg,xmin,xmax,ymin,ymax,x,f,n,u,scale,eps,name)
MULTIGRID *mg;
FLOAT xmin, xmax, ymin, ymax, x[NNPOINTS][2], scale,eps;
INT *f, n, u;
char name[];
{
   ELEMENT *pel;
   GRID *theGrid;
   NODE *n0, *n1, *n2;
   FLOAT *x0, *x1, *x2, l0, l1, l2, u0, u1, b[DIM2][DIM2];
   INT i;
   FILE *fp;
   
   fp = fopen(name,"a");
   for (theGrid = FIRSTGRID(mg); theGrid != NULL; theGrid = theGrid->finer)
      for (pel = FIRSTELEMENT(theGrid); pel != NULL; pel = pel->succ)
         if (IS_TOP_ELEMENT(pel)){
            VERTICES_OF_ELEMENT(x0,x1,x2,pel);
            barycentric_coordinates(x0,x1,x2,b); /* DOT(b[i],x)+b[i][2] is the
                             barycentric coordinate of x with respect to x_i */
            for (i = 0; i < n; i++){
               l0 = DOT(b[0],x[i])+b[0][2];
               l1 = DOT(b[1],x[i])+b[1][2];
               l2 = DOT(b[2],x[i])+b[2][2];
               if(f[i] && l0 >= eps && l1 >= eps && l2 >= eps){
                  NODES_OF_ELEMENT(n0,n1,n2,pel);
                  u0 = ND(n0,u,0)*l0 + ND(n1,u,0)*l1 + ND(n2,u,0)*l2;
                  u1 = ND(n0,u,1)*l0 + ND(n1,u,1)*l1 + ND(n2,u,1)*l2;
                  draw_arrow(x[i][0],x[i][1],u0,u1,scale,xmin,xmax,ymin,ymax,fp);
                  f[i] = 0;
               }
            }
         }
   fclose(fp);
}

#else

void find_elements(mg,xmin,xmax,ymin,ymax,x,f,n,u,scale,eps,name)
MULTIGRID *mg; FLOAT xmin, xmax, ymin, ymax, x[NNPOINTS][2], scale,eps;
INT *f, n, u; char name[];
{  eprintf("Error: find_elements not available.\n");  }

#endif

void save_arrows_for_gnuplot(mg,u,name)
MULTIGRID *mg;
INT u;
char name[];
{
   FLOAT xmin, xmax, ymin, ymax, d, x[NNPOINTS][2], eps=0.;
   FLOAT scale=4.;
   INT i, j, n, f[NNPOINTS];
   
   points_for_visualization(mg,&xmin,&xmax,&ymin,&ymax,x,f,&n);
   save_boundary_for_gnuplot(mg,name);
   d = MIN(xmax-xmin,ymax-ymin)/5.;
   xmin -= d;
   ymin -= d;
   xmax += d;
   ymax += d;
   find_elements(mg,xmin,xmax,ymin,ymax,x,f,n,u,scale,eps,name);
   i = j = 0;
   while (i < n && !j)
      if (f[i++]) j = 1;
   eps = -1.e-15;
   while (j && eps > 1.e-8){
      find_elements(mg,xmin,xmax,ymin,ymax,x,f,n,u,scale,eps,name);
      i = j = 0;
      while (i < n && !j)
         if (f[i++]) j = 1;
      eps *= 10.;
   }
}

void save_exact_arrows_for_gnuplot(mg,name,u1,u2)
MULTIGRID *mg;
char name[];
FLOAT (*u1)(), (*u2)();
{
   FLOAT xmin, xmax, ymin, ymax, d, x[NNPOINTS][2];
   FLOAT scale=4.;
   INT i, n, f[NNPOINTS];
   FILE *fp;
   
   points_for_visualization(mg,&xmin,&xmax,&ymin,&ymax,x,f,&n);
   save_boundary_for_gnuplot(mg,name);
   d = MIN(xmax-xmin,ymax-ymin)/5.;
   xmin -= d;
   ymin -= d;
   xmax += d;
   ymax += d;
   fp = fopen(name,"a");
   for (i = 0; i < n; i++)
      draw_arrow(x[i][0],x[i][1],u1(x[i]),u2(x[i]),
                 scale,xmin,xmax,ymin,ymax,fp);
   fclose(fp);
}

/******************************************************************************/

void draw_graph_on_triangle(fp,x1,x2,x3,e1,e2,e3,eps,i)
FLOAT *x1, *x2, *x3, e1, e2, e3, eps;
INT *i;
FILE *fp;
{
   if (fabs(e1) > eps || fabs(e2) > eps || fabs(e3) > eps){
      fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],e1);
      fprintf(fp,"%e  %e  %e\n",  x2[0],x2[1],e2);
      fprintf(fp,"%e  %e  %e\n",  x3[0],x3[1],e3);
      if (*i){
         fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],e1);
         *i = 0;
      }
      fprintf(fp,"%e  %e  %e\n\n",x1[0],x1[1],e1);
   }
}

void draw_graph_on_quadrilateral(fp,x1,x2,x3,x4,e1,e2,e3,e4,eps,i)
FLOAT *x1, *x2, *x3, *x4, e1, e2, e3, e4, eps;
INT *i;
FILE *fp;
{
   if (fabs(e1) > eps || fabs(e2) > eps || fabs(e3) > eps || fabs(e4) > eps){
      fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],e1);
      fprintf(fp,"%e  %e  %e\n",  x2[0],x2[1],e2);
      fprintf(fp,"%e  %e  %e\n",  x3[0],x3[1],e3);
      fprintf(fp,"%e  %e  %e\n",  x4[0],x4[1],e4);
      if (*i){
         fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],e1);
         *i = 0;
      }
      fprintf(fp,"%e  %e  %e\n\n",x1[0],x1[1],e1);
   }
}

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void gsnode_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      draw_graph_on_triangle(fp,x1,x2,x3,NDMV(n1,u,0),NDMV(n2,u,0),NDMV(n3,u,0),-1.,&i);
   }
}

void gsnode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, e1, e2, e3;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      e1 = k*(NDMV(n1,u,0) - (*u0)(x1));
      e2 = k*(NDMV(n2,u,0) - (*u0)(x2));
      e3 = k*(NDMV(n3,u,0) - (*u0)(x3));
      draw_graph_on_triangle(fp,x1,x2,x3,e1,e2,e3,eps,&i);
   }
}

void gp1x3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, xc[DIM], uc;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      POINT3(x1,x2,x3,xc);
      uc = (NDMV(n1,u,0) + NDMV(n2,u,0) + NDMV(n3,u,0))/3. + EDMV(pel,u,0);
      draw_graph_on_triangle(fp,x1,x2,xc,NDMV(n1,u,0),NDMV(n2,u,0),uc,-1.,&i);
      draw_graph_on_triangle(fp,x2,x3,xc,NDMV(n2,u,0),NDMV(n3,u,0),uc,-1.,&i);
      draw_graph_on_triangle(fp,x3,x1,xc,NDMV(n3,u,0),NDMV(n1,u,0),uc,-1.,&i);
   }
}

void gp2c_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], u01, u02, u12;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      u01 = 0.5*(NDMV(n0,u,0)+NDMV(n1,u,0)) + 0.25*FDMV(fa2,u,0);
      u02 = 0.5*(NDMV(n0,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa1,u,0);
      u12 = 0.5*(NDMV(n1,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa0,u,0);
      draw_graph_on_triangle(fp,x0,x01,x02,NDMV(n0,u,0),u01,u02,-1.,&i);
      draw_graph_on_triangle(fp,x1,x12,x01,NDMV(n1,u,0),u12,u01,-1.,&i);
      draw_graph_on_triangle(fp,x2,x02,x12,NDMV(n2,u,0),u02,u12,-1.,&i);
   }
}

void gp2c_difference_graph_for_gnuplot(fp,tGrid,u,v)
GRID *tGrid;
INT u, v;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], u01, u02, u12,
                                                      v01, v02, v12;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      u01 = 0.5*(NDMV(n0,u,0)+NDMV(n1,u,0)) + 0.25*FDMV(fa2,u,0);
      u02 = 0.5*(NDMV(n0,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa1,u,0);
      u12 = 0.5*(NDMV(n1,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa0,u,0);
      v01 = 0.5*(NDMV(n0,v,0)+NDMV(n1,v,0)) + 0.25*FDMV(fa2,v,0);
      v02 = 0.5*(NDMV(n0,v,0)+NDMV(n2,v,0)) + 0.25*FDMV(fa1,v,0);
      v12 = 0.5*(NDMV(n1,v,0)+NDMV(n2,v,0)) + 0.25*FDMV(fa0,v,0);
      draw_graph_on_triangle(fp,x0,x01,x02,
                             NDMV(n0,u,0)-NDMV(n0,v,0),u01-v01,u02-v02,-1.,&i);
      draw_graph_on_triangle(fp,x1,x12,x01,
                             NDMV(n1,u,0)-NDMV(n1,v,0),u12-v12,u01-v01,-1.,&i);
      draw_graph_on_triangle(fp,x2,x02,x12,
                             NDMV(n2,u,0)-NDMV(n2,v,0),u02-v02,u12-v12,-1.,&i);
   }
}

void gp2x3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, xc[DIM], x01[DIM], x02[DIM], x12[DIM], 
         x0c[DIM], x1c[DIM], x2c[DIM], u01, u02, u12, u0c, u1c, u2c, uc, r, s;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      POINT3(x0,x1,x2,xc);
      AVERAGE(x0,xc,x0c)
      AVERAGE(x1,xc,x1c)
      AVERAGE(x2,xc,x2c)
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      r = NDMV(n0,u,0) + NDMV(n1,u,0) + NDMV(n2,u,0);
      s = FDMV(fa0,u,0) + FDMV(fa1,u,0) + FDMV(fa2,u,0);
      uc = (3.*r + s)/9. + EDMV(pel,u,3);
      s = (1.5*r + s)/9. + 0.5*EDMV(pel,u,3);
      u0c = 0.5*NDMV(n0,u,0) - FDMV(fa0,u,0)/12. + 0.25*EDMV(pel,u,0) + s;
      u1c = 0.5*NDMV(n1,u,0) - FDMV(fa1,u,0)/12. + 0.25*EDMV(pel,u,1) + s;
      u2c = 0.5*NDMV(n2,u,0) - FDMV(fa2,u,0)/12. + 0.25*EDMV(pel,u,2) + s;
      u01 = 0.5*(NDMV(n0,u,0)+NDMV(n1,u,0)) + 0.25*FDMV(fa2,u,0);
      u02 = 0.5*(NDMV(n0,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa1,u,0);
      u12 = 0.5*(NDMV(n1,u,0)+NDMV(n2,u,0)) + 0.25*FDMV(fa0,u,0);
      draw_graph_on_triangle(fp,x0,x01,x0c,NDMV(n0,u,0),u01,u0c,-1.,&i);
      draw_graph_on_triangle(fp,x0,x02,x0c,NDMV(n0,u,0),u02,u0c,-1.,&i);
      draw_graph_on_triangle(fp,x1,x12,x1c,NDMV(n1,u,0),u12,u1c,-1.,&i);
      draw_graph_on_triangle(fp,x1,x01,x1c,NDMV(n1,u,0),u01,u1c,-1.,&i);
      draw_graph_on_triangle(fp,x2,x02,x2c,NDMV(n2,u,0),u02,u2c,-1.,&i);
      draw_graph_on_triangle(fp,x2,x12,x2c,NDMV(n2,u,0),u12,u2c,-1.,&i);
      draw_graph_on_triangle(fp,x0c,x1c,xc,u0c,u1c,uc,-1.,&i);
      draw_graph_on_triangle(fp,x1c,x2c,xc,u1c,u2c,uc,-1.,&i);
      draw_graph_on_triangle(fp,x2c,x0c,xc,u2c,u0c,uc,-1.,&i);
   }
}

void gp2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], e0, e1, e2, e01, e02, e12;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      e0  = k*(NDMV(n0,u,0) - (*u0)(x0));
      e1  = k*(NDMV(n1,u,0) - (*u0)(x1));
      e2  = k*(NDMV(n2,u,0) - (*u0)(x2));
      e01 = k*(0.5*(NDMV(n0,u,0)+NDMV(n1,u,0))+0.25*FDMV(fa2,u,0) - (*u0)(x01));
      e02 = k*(0.5*(NDMV(n0,u,0)+NDMV(n2,u,0))+0.25*FDMV(fa1,u,0) - (*u0)(x02));
      e12 = k*(0.5*(NDMV(n1,u,0)+NDMV(n2,u,0))+0.25*FDMV(fa0,u,0) - (*u0)(x12));
      draw_graph_on_triangle(fp,x0,x01,x02,e0,e01,e02,eps,&i);
      draw_graph_on_triangle(fp,x1,x12,x01,e1,e12,e01,eps,&i);
      draw_graph_on_triangle(fp,x2,x02,x12,e2,e02,e12,eps,&i);
   }
}

void gq1_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FLOAT *x1, *x2, *x3, *x4;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      draw_graph_on_quadrilateral(fp,x1,x2,x3,x4,
                    NDMV(n1,u,0),NDMV(n2,u,0),NDMV(n3,u,0),NDMV(n4,u,0),-1.,&i);
   }
}

void gq1x4_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FACE *fa1, *fa2, *fa3, *fa4;
   FLOAT *x1, *x2, *x3, *x4, x12[DIM], x23[DIM], x34[DIM], x41[DIM], x1234[DIM];
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      FACES_OF_4ELEMENT(fa1,fa2,fa3,fa4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      AVERAGE(x1,x2,x12)
      AVERAGE(x2,x3,x23)
      AVERAGE(x3,x4,x34)
      AVERAGE(x4,x1,x41)
      POINT4(x1,x2,x3,x4,x1234);
      draw_graph_on_quadrilateral(fp,x1,x12,x1234,x41,
                 NDMV(n1,u,0),FDMV(fa1,u,0),EDMV(pel,u,0),FDMV(fa4,u,0),-1.,&i);
      draw_graph_on_quadrilateral(fp,x2,x23,x1234,x12,
                 NDMV(n2,u,0),FDMV(fa2,u,0),EDMV(pel,u,0),FDMV(fa1,u,0),-1.,&i);
      draw_graph_on_quadrilateral(fp,x3,x34,x1234,x23,
                 NDMV(n3,u,0),FDMV(fa3,u,0),EDMV(pel,u,0),FDMV(fa2,u,0),-1.,&i);
      draw_graph_on_quadrilateral(fp,x4,x41,x1234,x34,
                 NDMV(n4,u,0),FDMV(fa4,u,0),EDMV(pel,u,0),FDMV(fa3,u,0),-1.,&i);
   }
}

void gq2_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FACE *fa1, *fa2, *fa3, *fa4;
   FLOAT *x1, *x2, *x3, *x4, x12[DIM], x23[DIM], x34[DIM], x41[DIM], x1234[DIM],
         u12, u23, u34, u41, u1234;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      FACES_OF_4ELEMENT(fa1,fa2,fa3,fa4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      AVERAGE(x1,x2,x12)
      AVERAGE(x2,x3,x23)
      AVERAGE(x3,x4,x34)
      AVERAGE(x4,x1,x41)
      POINT4(x1,x2,x3,x4,x1234);
      u12 = 0.5*(NDMV(n1,u,0)+NDMV(n2,u,0)) + FDMV(fa1,u,0);
      u23 = 0.5*(NDMV(n2,u,0)+NDMV(n3,u,0)) + FDMV(fa2,u,0);
      u34 = 0.5*(NDMV(n3,u,0)+NDMV(n4,u,0)) + FDMV(fa3,u,0);
      u41 = 0.5*(NDMV(n4,u,0)+NDMV(n1,u,0)) + FDMV(fa4,u,0);
      u1234 = 0.25*(NDMV(n1,u,0)+NDMV(n2,u,0)+NDMV(n3,u,0)+NDMV(n4,u,0))
          +0.5*(FDMV(fa1,u,0)+FDMV(fa2,u,0)+FDMV(fa3,u,0)+FDMV(fa4,u,0))
          +EDMV(pel,u,0);
      draw_graph_on_quadrilateral(fp,x1,x12,x1234,x41,
                           NDMV(n1,u,0),u12,u1234,u41,-1.,&i);
      draw_graph_on_quadrilateral(fp,x2,x23,x1234,x12,
                           NDMV(n2,u,0),u23,u1234,u12,-1.,&i);
      draw_graph_on_quadrilateral(fp,x3,x34,x1234,x23,
                           NDMV(n3,u,0),u34,u1234,u23,-1.,&i);
      draw_graph_on_quadrilateral(fp,x4,x41,x1234,x34,
                           NDMV(n4,u,0),u41,u1234,u34,-1.,&i);
   }
}

void gq2x4_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   FLOAT x[5][5][DIM], s[5][5], f[4][2];
   INT i, j, k=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      SET1(x[0][0],pel->n[0]->myvertex->x)
      SET1(x[4][0],pel->n[1]->myvertex->x)
      SET1(x[4][4],pel->n[2]->myvertex->x)
      SET1(x[0][4],pel->n[3]->myvertex->x)
      AVERAGE(x[0][0],x[4][0],x[2][0])
      AVERAGE(x[0][0],x[2][0],x[1][0])
      AVERAGE(x[2][0],x[4][0],x[3][0])
      AVERAGE(x[0][4],x[4][4],x[2][4])
      AVERAGE(x[0][4],x[2][4],x[1][4])
      AVERAGE(x[2][4],x[4][4],x[3][4])
      for (i = 0; i < 5; i++){
         AVERAGE(x[i][0],x[i][4],x[i][2])
         AVERAGE(x[i][0],x[i][2],x[i][1])
         AVERAGE(x[i][2],x[i][4],x[i][3])
      }
      for (i = 0; i < 4; i++){
         j = i + 1;
         if (j == 4)
            j = 0;
         if (pel->n[i]->index < pel->n[j]->index){
            f[i][0] = FDMV(pel->f[i],u,0);
            f[i][1] = FDMV(pel->f[i],u,2);
         }
         else{
            f[i][0] = FDMV(pel->f[i],u,2);
            f[i][1] = FDMV(pel->f[i],u,0);
         }
      }
      s[0][0] = NDMV(pel->n[0],u,0);
      s[4][0] = NDMV(pel->n[1],u,0);
      s[4][4] = NDMV(pel->n[2],u,0);
      s[0][4] = NDMV(pel->n[3],u,0);
      s[2][0] = FDMV(pel->f[0],u,1);
      s[4][2] = FDMV(pel->f[1],u,1);
      s[2][4] = FDMV(pel->f[2],u,1);
      s[0][2] = FDMV(pel->f[3],u,1);
      s[2][2] = EDMV(pel,u,8);
      s[1][0] = f[0][0] + 0.5*(s[0][0] + s[2][0]);
      s[3][0] = f[0][1] + 0.5*(s[2][0] + s[4][0]);
      s[4][1] = f[1][0] + 0.5*(s[4][0] + s[4][2]);
      s[4][3] = f[1][1] + 0.5*(s[4][2] + s[4][4]);
      s[3][4] = f[2][0] + 0.5*(s[4][4] + s[2][4]);
      s[1][4] = f[2][1] + 0.5*(s[2][4] + s[0][4]);
      s[0][3] = f[3][0] + 0.5*(s[0][4] + s[0][2]);
      s[0][1] = f[3][1] + 0.5*(s[0][2] + s[0][0]);
      s[2][1] = EDMV(pel,u,0) + 0.5*(s[2][0] + s[2][2]);
      s[3][2] = EDMV(pel,u,1) + 0.5*(s[4][2] + s[2][2]);
      s[2][3] = EDMV(pel,u,2) + 0.5*(s[2][4] + s[2][2]);
      s[1][2] = EDMV(pel,u,3) + 0.5*(s[0][2] + s[2][2]);
      s[1][1] = EDMV(pel,u,4) - 0.25*(s[0][0]+s[2][0]+s[2][2]+s[0][2])
                              +  0.5*(s[1][0]+s[2][1]+s[1][2]+s[0][1]);
      s[3][1] = EDMV(pel,u,5) - 0.25*(s[4][0]+s[4][2]+s[2][2]+s[2][0])
                              +  0.5*(s[4][1]+s[3][2]+s[2][1]+s[3][0]);
      s[3][3] = EDMV(pel,u,6) - 0.25*(s[4][4]+s[2][4]+s[2][2]+s[4][2])
                              +  0.5*(s[3][4]+s[2][3]+s[3][2]+s[4][3]);
      s[1][3] = EDMV(pel,u,7) - 0.25*(s[0][4]+s[0][2]+s[2][2]+s[2][4])
                              +  0.5*(s[0][3]+s[1][2]+s[2][3]+s[1][4]);
      for (i = 0; i < 4; i++)
         for (j = 0; j < 4; j++)
            draw_graph_on_quadrilateral(fp,
                                x[i][j],x[i+1][j],x[i+1][j+1],x[i][j+1],
                                s[i][j],s[i+1][j],s[i+1][j+1],s[i][j+1],-1.,&k);
   }
}

void gq2b3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FACE *fa1, *fa2, *fa3, *fa4;
   FLOAT *x1, *x2, *x3, *x4, x12[DIM], x23[DIM], x34[DIM], x41[DIM], x1234[DIM],
         u12, u23, u34, u41, u1234;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      FACES_OF_4ELEMENT(fa1,fa2,fa3,fa4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      AVERAGE(x1,x2,x12)
      AVERAGE(x2,x3,x23)
      AVERAGE(x3,x4,x34)
      AVERAGE(x4,x1,x41)
      POINT4(x1,x2,x3,x4,x1234);
      u12 = 0.5*(NDMV(n1,u,0)+NDMV(n2,u,0)) + FDMV(fa1,u,0);
      u23 = 0.5*(NDMV(n2,u,0)+NDMV(n3,u,0)) + FDMV(fa2,u,0);
      u34 = 0.5*(NDMV(n3,u,0)+NDMV(n4,u,0)) + FDMV(fa3,u,0);
      u41 = 0.5*(NDMV(n4,u,0)+NDMV(n1,u,0)) + FDMV(fa4,u,0);
      u1234 = 0.25*(NDMV(n1,u,0)+NDMV(n2,u,0)+NDMV(n3,u,0)+NDMV(n4,u,0))
          +0.5*(FDMV(fa1,u,0)+FDMV(fa2,u,0)+FDMV(fa3,u,0)+FDMV(fa4,u,0))
          +EDMV(pel,u,0);
      if (IS_BN(n1) || IS_BN(n2) || IS_BN(n3) || IS_BN(n4)){
         u12 += 0.5625*FDMV(fa1,u,1);
         u23 += 0.5625*FDMV(fa2,u,1);
         u34 += 0.5625*FDMV(fa3,u,1);
         u41 += 0.5625*FDMV(fa4,u,1);
         u1234 += 81./256*(EDMV(pel,u,1)+EDMV(pel,u,2)+EDMV(pel,u,3))
            - 9./256.*(FDMV(fa1,u,1)+FDMV(fa2,u,1)+FDMV(fa3,u,1)+FDMV(fa4,u,1));
      }
      draw_graph_on_quadrilateral(fp,x1,x12,x1234,x41,
                           NDMV(n1,u,0),u12,u1234,u41,-1.,&i);
      draw_graph_on_quadrilateral(fp,x2,x23,x1234,x12,
                           NDMV(n2,u,0),u23,u1234,u12,-1.,&i);
      draw_graph_on_quadrilateral(fp,x3,x34,x1234,x23,
                           NDMV(n3,u,0),u34,u1234,u23,-1.,&i);
      draw_graph_on_quadrilateral(fp,x4,x41,x1234,x34,
                           NDMV(n4,u,0),u41,u1234,u34,-1.,&i);
   }
}

#else

void gsnode_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gsnode_solution_graph_for_gnuplot not available.\n");  } 

void gsnode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: gsnode_error_graph_for_gnuplot not available.\n");  } 

void gp1x3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gp1x3_solution_graph_for_gnuplot not available.\n");  }

void gp2c_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gp2c_solution_graph_for_gnuplot not available.\n");  }

void gp2x3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gp2x3_solution_graph_for_gnuplot not available.\n");  }

void gp2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: gp2c_error_graph_for_gnuplot not available.\n");  }

void gq1_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gq1_solution_graph_for_gnuplot not available.\n");  } 

void gq1x4_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gq1x4_solution_graph_for_gnuplot not available.\n");  } 

void gq2_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gq2_solution_graph_for_gnuplot not available.\n");  } 

void gq2x4_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gq2x4_solution_graph_for_gnuplot not available.\n");  } 

void gq2b3_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: gq2b3_solution_graph_for_gnuplot not available.\n");  } 

#endif

#if E_DATA & SCALAR_ELEMENT_DATA

void sel_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   FLOAT *x1, *x2, *x3;
   DOUBLE e;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      e = ED(pel,u);
      draw_graph_on_triangle(fp,x1,x2,x3,e,e,e,-1.,&i);
   }
}

#else

void sel_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{ eprintf("Error: sel_solution_graph_for_gnuplot not available.\n"); }

#endif

#if E_DATA & VECTOR_ELEMENT_DATA

void vel_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid;
INT u, j;
FILE *fp;
{
   ELEMENT *pel;
   FLOAT *x1, *x2, *x3;
   DOUBLE e;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      e = EDV(pel,u,j);
      draw_graph_on_triangle(fp,x1,x2,x3,e,e,e,-1.,&i);
   }
}

#else

void vel_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid; INT u, j; FILE *fp;
{ eprintf("Error: vel_solution_graph_for_gnuplot not available.\n"); }

#endif

#if N_DATA & SCALAR_NODE_DATA

void snode_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      draw_graph_on_triangle(fp,x1,x2,x3,NDS(n1,u),NDS(n2,u),NDS(n3,u),-1.,&i);
   }
}

void snode_solution_graph_for_gnuplot0(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, tol=0.005;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      if (fabs(NDS(n1,u)) > tol || fabs(NDS(n2,u)) > tol || 
          fabs(NDS(n3,u)) > tol)
      draw_graph_on_triangle(fp,x1,x2,x3,NDS(n1,u),NDS(n2,u),NDS(n3,u),-1.,&i);
   }
   tol *= 1.5;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      for (i = 0; i < 3; i++)
         if (IS_BF(pel->f[i])){
            n1 = pel->n[i];
            n2 = pel->n[(i+1)%3];
            n3 = pel->n[(i+2)%3];
            if (fabs(NDS(n1,u)) < tol && fabs(NDS(n2,u)) < tol && 
                fabs(NDS(n3,u)) < tol){
               fprintf(fp,"%e  %e  %e\n",  n2->myvertex->x[0],
                                           n2->myvertex->x[1],NDS(n2,u));
               fprintf(fp,"%e  %e  %e\n\n",n3->myvertex->x[0],
                                           n3->myvertex->x[1],NDS(n3,u));
            }
         }
}

void snode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, e1, e2, e3;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      e1 = k*(NDS(n1,u) - (*u0)(x1));
      e2 = k*(NDS(n2,u) - (*u0)(x2));
      e3 = k*(NDS(n3,u) - (*u0)(x3));
      draw_graph_on_triangle(fp,x1,x2,x3,e1,e2,e3,eps,&i);
   }
}

void q1_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FLOAT *x1, *x2, *x3, *x4;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      draw_graph_on_quadrilateral(fp,x1,x2,x3,x4,
                                NDS(n1,u),NDS(n2,u),NDS(n3,u),NDS(n4,u),-1.,&i);
   }
}

void q1_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3, *n4;
   FLOAT *x1, *x2, *x3, *x4, e1, e2, e3, e4;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
      VERTICES_OF_4ELEMENT(x1,x2,x3,x4,pel);
      e1 = k*(NDS(n1,u) - (*u0)(x1));
      e2 = k*(NDS(n2,u) - (*u0)(x2));
      e3 = k*(NDS(n3,u) - (*u0)(x3));
      e4 = k*(NDS(n4,u) - (*u0)(x4));
      draw_graph_on_quadrilateral(fp,x1,x2,x3,x4,e1,e2,e3,e4,eps,&i);
   }
}

void draw_cut_for_gnuplot_p1c(tGrid,u,name,x,v) /* cut along the line x+a*v */
GRID *tGrid;
INT u;
char name[];
DOUBLE x[DIM], v[DIM];
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3;
   DOUBLE u1[DIM], u2[DIM], u3[DIM], a, b1, b2, b3;
   INT i=1, int1, int2, int3;
   FILE *fp;

   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      int1 = intersection_with_edge(x,v,x2,x3,&a,&b1,u1);
      int2 = intersection_with_edge(x,v,x3,x1,&a,&b2,u2);
      int3 = intersection_with_edge(x,v,x1,x2,&a,&b3,u3);
      if (int1 == -1){
         fprintf(fp,"%e  %e  %e\n",  x2[0],x2[1],NDS(n2,u));
         fprintf(fp,"%e  %e  %e\n\n",x3[0],x3[1],NDS(n3,u));
      }
      else if (int2 == -1){
         fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],NDS(n1,u));
         fprintf(fp,"%e  %e  %e\n\n",x3[0],x3[1],NDS(n3,u));
      }
      else if (int3 == -1){
         fprintf(fp,"%e  %e  %e\n",  x1[0],x1[1],NDS(n1,u));
         fprintf(fp,"%e  %e  %e\n\n",x2[0],x2[1],NDS(n2,u));
      }
      else if (int1 + int2 + int3 > 1){
         if (int1 == 1 && int2 == 1 && 
             fabs(u1[0]-u2[0]) + fabs(u1[1]-u2[1]) < 1.e-20)
            int2 = 0;
         if (int2 == 1 && int3 == 1 && 
             fabs(u2[0]-u3[0]) + fabs(u2[1]-u3[1]) < 1.e-20)
            int3 = 0;
         if (int3 == 1 && int1 == 1 && 
             fabs(u3[0]-u1[0]) + fabs(u3[1]-u1[1]) < 1.e-20)
            int1 = 0;
         if (int1 == 1)
            fprintf(fp,"%e  %e  %e\n",u1[0],u1[1],
                                      (1.-b1)*NDS(n2,u)+b1*NDS(n3,u));
         if (int2 == 1)
            fprintf(fp,"%e  %e  %e\n",u2[0],u2[1],
                                      (1.-b2)*NDS(n3,u)+b2*NDS(n1,u));
         if (int3 == 1)
            fprintf(fp,"%e  %e  %e\n",u3[0],u3[1],
                                      (1.-b3)*NDS(n1,u)+b3*NDS(n2,u));
         if (i){
            if (int3 == 0)
               fprintf(fp,"%e  %e  %e\n",u2[0],u2[1],
                                         (1.-b2)*NDS(n3,u)+b2*NDS(n1,u));
            else
               fprintf(fp,"%e  %e  %e\n",u3[0],u3[1],
                                         (1.-b3)*NDS(n1,u)+b3*NDS(n2,u));
            i = 0;
         }
         fprintf(fp,"\n");
      }
   }
   fclose(fp);
}

#else

void snode_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: snode_solution_graph_for_gnuplot not available.\n");  } 

void snode_solution_graph_for_gnuplot0(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: snode_solution_graph_for_gnuplot0 not available.\n");  } 

void snode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: snode_error_graph_for_gnuplot not available.\n");  } 

void q1_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: q1_solution_graph_for_gnuplot not available.\n");  } 

void q1_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: q1_error_graph_for_gnuplot not available.\n");  } 

void draw_cut_for_gnuplot_p1c(tGrid,u,name,x,v) /* cut along the line x+a*v */
GRID *tGrid; INT u; char name[]; DOUBLE x[DIM], v[DIM];
{  eprintf("Error: draw_cut_for_gnuplot_p1cn not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (E_DATA & VECTOR_ELEMENT_DATA)

INT find_oscillation_regions_p1c(tGrid,d,v,w)
GRID *tGrid;
DOUBLE d;
INT v, w; /* v ... scalar node variable, w ... vector element variable */
{
   ELEMENT *pel;
   NODE *pn;
   DOUBLE a[10000], r;
   INT i, j, m=0, n=0;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      EDV(pel,w,1) = 0.;
   a[0] = -1.;
   for (pn = FIRSTNODE(tGrid); pn; pn = pn->succ)
      if ((r=fabs(NDS(pn,v))) > d){
         n++;
         for (i = 0; r < a[i]; i++);
         for (j = n; j > i; j--)
            a[j] = a[j-1];
         a[i] = r;
      }
      else
         m++;
   r = a[n/2];
   printf("oscillations: %i yes, %i no; r = %f\n",n,m,r);
   for (pn = FIRSTNODE(tGrid); pn; pn = pn->succ)
      if (fabs(NDS(pn,v)) > d)
         NDS(pn,v) = fabs(NDS(pn,v))/r;
      else
         NDS(pn,v) = 0.;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      r = MAX(NDS(pel->n[0],v),NDS(pel->n[1],v));
      EDV(pel,w,1) = MAX(r,NDS(pel->n[2],v));
   }
   if (n == 0){
      eprintf("There are no oscillations.\n");
      return(0);
   }
   else{
      if (n > m/4)
         eprintf("Too much nodes with oscillations.\n");
      return(1);
   }
}

#else

void find_oscillation_regions_p1c(tGrid,d,v,w)
GRID *tGrid; DOUBLE d; INT v, w;
{  eprintf("Error: find_oscillation_regions_p1c not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (DATA_S & N_LINK_TO_ELEMENTS)

void averaging_in_crosswind_direction_p1c(tGrid,u,u_aver)
GRID *tGrid;
INT u, u_aver;
{
   NELINK *pnel;
   ELEMENT *pel;
   NODE *pn;
   DOUBLE *x0, *x1, *x2, v[DIM], z[DIM], a, b, am, um, ap, up;
   INT im, ip, i, j, k;

   for (pn = FIRSTNODE(tGrid); pn; pn = pn->succ)
      if (NOT_BN(pn)){
         v[0] = -bb1(pn->myvertex->x);
         v[1] =  bb0(pn->myvertex->x);
         im = ip = 1;
         for (pnel = NESTART(pn); pnel && (im || ip); pnel = pnel->next){
            pel = pnel->nbel;
            if (pn == pel->n[0])
               i = 0;
            else if (pn == pel->n[1])
               i = 1;
            else
               i = 2;
            j = (i+1) % SIDES;
            k = (i+2) % SIDES;
            x0 = pel->n[i]->myvertex->x;
            x1 = pel->n[j]->myvertex->x;
            x2 = pel->n[k]->myvertex->x;
            if (intersection_with_edge(x0,v,x1,x2,&a,&b,z) == 1){
               if (a < 0. && im){
                  im = 0;
                  am = -a;
                  um = (1.-b)*NDS(pel->n[j],u) + b*NDS(pel->n[k],u);
               }
               if (a > 0. && ip){
                  ip = 0;
                  ap = a;
                  up = (1.-b)*NDS(pel->n[j],u) + b*NDS(pel->n[k],u);
               }
            }
         }
         NDS(pn,u_aver) = (ap*um + am*up)/(ap + am);
      }
}

#else

void averaging_in_crosswind_direction_p1c(tGrid,u,u_aver)
GRID *tGrid; INT u, u_aver;
{  eprintf("Error: averaging_in_crosswind_direction_p1c not available.\n");  }

#endif

#if N_DATA & VECTOR_NODE_DATA

void vnode_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid;
INT u, j;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      draw_graph_on_triangle(fp,x1,x2,x3,
                             ND(n1,u,j),ND(n2,u,j),ND(n3,u,j),-1.,&i);
   }
}

void vnode_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps)
GRID *tGrid;
INT u, j;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n1, *n2, *n3;
   FLOAT *x1, *x2, *x3, e1, e2, e3;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n1,n2,n3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      e1 = k*(ND(n1,u,j) - (*u0)(x1));
      e2 = k*(ND(n2,u,j) - (*u0)(x2));
      e3 = k*(ND(n3,u,j) - (*u0)(x3));
      draw_graph_on_triangle(fp,x1,x2,x3,e1,e2,e3,eps,&i);
   }
}

#else

void vnode_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid; INT u, j; FILE *fp;
{  eprintf("Error: vnode_solution_graph_for_gnuplot not available.\n");  } 

void vnode_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps)
GRID *tGrid; INT u, j; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: vnode_error_graph_for_gnuplot not available.\n");  } 

#endif

#if F_DATA & SCALAR_FACE_DATA

void f_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM];
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
      draw_graph_on_triangle(fp,xc1,xc2,xc3,FD(f1,u),FD(f2,u),FD(f3,u),-1.,&i);
   }
}

void f_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM], e1, e2, e3;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
      e1 = k*(FD(f1,u) - (*u0)(xc1));
      e2 = k*(FD(f2,u) - (*u0)(xc2));
      e3 = k*(FD(f3,u) - (*u0)(xc3));
      draw_graph_on_triangle(fp,xc1,xc2,xc3,e1,e2,e3,eps,&i);
   }
}

#else  /*  !(F_DATA & SCALAR_FACE_DATA)  */

void f_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: f_solution_graph_for_gnuplot not available.\n");  }

void f_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: f_error_graph_for_gnuplot not available.\n");  }

#endif

#if F_DATA & VECTOR_FACE_DATA

void vf_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid;
INT u, j;
FILE *fp;
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM];
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
      draw_graph_on_triangle(fp,xc1,xc2,xc3,
                             FDV(f1,u,j),FDV(f2,u,j),FDV(f3,u,j),-1.,&i);
   }
}

void vf_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps)
GRID *tGrid;
INT u, j;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   FACE *f1, *f2, *f3;
   FLOAT *x1, *x2, *x3, xc1[DIM], xc2[DIM], xc3[DIM], e1, e2, e3;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      FACES_OF_ELEMENT(f1,f2,f3,pel);
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      MIDPOINTS(x1,x2,x3,xc3,xc2,xc1)
      e1 = k*(FDV(f1,u,j) - (*u0)(xc1));
      e2 = k*(FDV(f2,u,j) - (*u0)(xc2));
      e3 = k*(FDV(f3,u,j) - (*u0)(xc3));
      draw_graph_on_triangle(fp,xc1,xc2,xc3,e1,e2,e3,eps,&i);
   }
}

#else  /*  !(F_DATA & VECTOR_FACE_DATA)  */

void vf_solution_graph_for_gnuplot(fp,tGrid,u,j)
GRID *tGrid; INT u, j; FILE *fp;
{  eprintf("Error: vf_solution_graph_for_gnuplot not available.\n");  } 

void vf_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps)
GRID *tGrid; INT u, j; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: vf_error_graph_for_gnuplot not available.\n");  } 

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void p2c_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid;
INT u;
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], u01, u02, u12;
   INT i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      u01 = 0.5*(NDS(n0,u)+NDS(n1,u)) + 0.25*FD(fa2,u);
      u02 = 0.5*(NDS(n0,u)+NDS(n2,u)) + 0.25*FD(fa1,u);
      u12 = 0.5*(NDS(n1,u)+NDS(n2,u)) + 0.25*FD(fa0,u);
      draw_graph_on_triangle(fp,x0,x01,x02,NDS(n0,u),u01,u02,-1.,&i);
      draw_graph_on_triangle(fp,x1,x12,x01,NDS(n1,u),u12,u01,-1.,&i);
      draw_graph_on_triangle(fp,x2,x02,x12,NDS(n2,u),u02,u12,-1.,&i);
   }
}

void p2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid;
INT u;
FLOAT eps, (*u0)();
FILE *fp;
{
   ELEMENT *pel;
   NODE *n0, *n1, *n2;
   FACE *fa0, *fa1, *fa2;
   FLOAT *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], e0, e1, e2, e01, e02, e12;
   INT k=-1, i=1;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){ 
      NODES_OF_ELEMENT(n0,n1,n2,pel);
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      e0  = k*(NDS(n0,u) - (*u0)(x0));
      e1  = k*(NDS(n1,u) - (*u0)(x1));
      e2  = k*(NDS(n2,u) - (*u0)(x2));
      e01 = k*(0.5*(NDS(n0,u)+NDS(n1,u)) + 0.25*FD(fa2,u) - (*u0)(x01));
      e02 = k*(0.5*(NDS(n0,u)+NDS(n2,u)) + 0.25*FD(fa1,u) - (*u0)(x02));
      e12 = k*(0.5*(NDS(n1,u)+NDS(n2,u)) + 0.25*FD(fa0,u) - (*u0)(x12));
      draw_graph_on_triangle(fp,x0,x01,x02,e0,e01,e02,eps,&i);
      draw_graph_on_triangle(fp,x1,x12,x01,e1,e12,e01,eps,&i);
      draw_graph_on_triangle(fp,x2,x02,x12,e2,e02,e12,eps,&i);
   }
}

#else

void p2c_solution_graph_for_gnuplot(fp,tGrid,u)
GRID *tGrid; INT u; FILE *fp;
{  eprintf("Error: p2c_solution_graph_for_gnuplot not available.\n");  }

void p2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps)
GRID *tGrid; INT u; FLOAT eps, (*u0)(); FILE *fp;
{  eprintf("Error: p2c_error_graph_for_gnuplot not available.\n");  }

#endif

void exact_solution_graph_for_gnuplot(tGrid,u,name)
GRID *tGrid;
FLOAT (*u)();
char name[];
{
   ELEMENT *pel;
   FLOAT *x1, *x2, *x3;
   INT i=1;
   FILE *fp;

   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      VERTICES_OF_ELEMENT(x1,x2,x3,pel);
      draw_graph_on_triangle(fp,x1,x2,x3,u(x1),u(x2),u(x3),-1.,&i);
   }
   fclose(fp);
}

void solution_graph_for_gnuplot(tGrid,u,j,name,space,structure)
GRID *tGrid;
INT u, j, space, structure;
char name[];
{
   FILE *fp;

   fp = fopen(name,"w");
   switch(space){
   case P0:    if (structure == SCALAR)
                  sel_solution_graph_for_gnuplot(fp,tGrid,u);
               else if (structure == VECTOR)
                  vel_solution_graph_for_gnuplot(fp,tGrid,u,j);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case P1_NC: if (structure == SCALAR)
                  f_solution_graph_for_gnuplot(fp,tGrid,u);
               else if (structure == VECTOR)
                  vf_solution_graph_for_gnuplot(fp,tGrid,u,j);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case P1_MOD: if (structure == SCALAR)
                  vf_solution_graph_for_gnuplot(fp,tGrid,u,0);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case P1C:   if (structure == SCALAR)
                  snode_solution_graph_for_gnuplot(fp,tGrid,u);
               else if (structure == VECTOR)
                  vnode_solution_graph_for_gnuplot(fp,tGrid,u,j);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GP1C:
   case GP1C_ELBUB: if (structure == SCALAR)
                  gsnode_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GP1X3C: if (structure == SCALAR)
                  gp1x3_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case P1C_ELBUB:   if (structure == SCALAR)
                  snode_solution_graph_for_gnuplot(fp,tGrid,u);
               else if (structure == VECTOR)
                  vnode_solution_graph_for_gnuplot(fp,tGrid,u,j);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case IP1C:   if (structure == SCALAR)
                  snode_solution_graph_for_gnuplot(fp,tGrid,u);
               else if (structure == VECTOR)
                  vnode_solution_graph_for_gnuplot(fp,tGrid,u,j);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case Q1C:   if (structure == SCALAR)
                  q1_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GQ1C:
   case GQ1C_ELBUB:  if (structure == SCALAR)
                  gq1_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GQ1X4C: if (structure == SCALAR)
                  gq1x4_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case P2C:   if (structure == SCALAR)
                  p2c_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GP2C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:  if (structure == SCALAR)
                  gp2c_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GP2X3C: if (structure == SCALAR)
                  gp2x3_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GQ2C:
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB:  if (structure == SCALAR)
                  gq2_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GQ2X4C: if (structure == SCALAR)
                  gq2x4_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   case GQ2B3C: if (structure == SCALAR)
                  gq2b3_solution_graph_for_gnuplot(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   default:
        eprintf("Error: solution_graph_for_gnuplot not available.\n");
        break;
   }
   fclose(fp);
}

void solution_graph_for_gnuplot0(tGrid,u,j,name,space,structure)
GRID *tGrid;
INT u, j, space, structure;
char name[];
{
   FILE *fp;

   fp = fopen(name,"w");
   switch(space){
   case P1C:   if (structure == SCALAR)
                  snode_solution_graph_for_gnuplot0(fp,tGrid,u);
               else
                  eprintf("Error: solution_graph_for_gnuplot0 not available.\n");
        break;
   default:
        eprintf("Error: solution_graph_for_gnuplot0 not available.\n");
        break;
   }
   fclose(fp);
}

void error_graph_for_gnuplot(tGrid,u0,u,j,eps,name,space,structure)
GRID *tGrid;
FLOAT eps, (*u0)();
INT u, j, space, structure;
char name[];
{
   FILE *fp;

   fp = fopen(name,"w");
   switch(space){
   case P1_NC: if (structure == SCALAR)
                  f_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else if (structure == VECTOR)
                  vf_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case P1_MOD: if (structure == SCALAR)
                  vf_error_graph_for_gnuplot(fp,tGrid,u0,u,0,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case P1C:   if (structure == SCALAR)
                  snode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else if (structure == VECTOR)
                  vnode_error_graph_for_gnuplot(fp,tGrid,u0,u,j,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case GP1C:  
   case GP1C_ELBUB: if (structure == SCALAR)
                  gsnode_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case Q1C:   if (structure == SCALAR)
                  q1_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case P2C:   if (structure == SCALAR)
                  p2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case GP2C:  if (structure == SCALAR)
                  gp2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   case GP2C_3ELBUB:
   case GP2C_6ELBUB:  if (structure == SCALAR)
                  gp2c_error_graph_for_gnuplot(fp,tGrid,u0,u,eps);
               else
                  eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   default:
        eprintf("Error: error_graph_for_gnuplot not available.\n");
        break;
   }
   fclose(fp);
}

void p1_isolines_on_triangle(fp,x0,x1,x2,u0,u1,u2,iso_val,delta)
DOUBLE *x0, *x1, *x2, u0, u1, u2, iso_val, delta;
FILE *fp;
{
   DOUBLE min=u0, max=u0, iso, start[DIM], end[DIM], p, q, eps=1.e-15*delta;

   if (u1 > u0)
      max = u1;
   else
      min = u1; 
   if (u2 > max)
      max = u2;
   else if (u2 < min)
      min = u2;
   iso = iso_val + delta*ceil((min-iso_val)/delta);
   if (max - min > eps)
      while (iso <= max){
         if (fabs(u0-iso) < eps && fabs(u1-iso) < eps){
            SET1(start,x0)
            SET1(end,x1)
         }
         else if (fabs(u0-iso) < eps && fabs(u2-iso) < eps){
            SET1(start,x0)
            SET1(end,x2)
         }
         else if (fabs(u1-iso) < eps && fabs(u2-iso) < eps){
            SET1(start,x1)
            SET1(end,x2)
         }
         else if (fabs(u0-iso) < eps && min < u0 && u0 < max){
            SET1(start,x0)
            p = (iso-u2)/(u1-u2);
            q = 1.-p;
            SET20(end,x1,p,x2,q)
         }
         else if (fabs(u1-iso) < eps && min < u1 && u1 < max){
            SET1(start,x1)
            p = (iso-u2)/(u0-u2);
            q = 1.-p;
            SET20(end,x0,p,x2,q)
         }
         else if (fabs(u2-iso) < eps && min < u2 && u2 < max){
            SET1(start,x2)
            p = (iso-u1)/(u0-u1);
            q = 1.-p;
            SET20(end,x0,p,x1,q)
         }
         else{
            if ((u0 <= iso && iso <= u1) || (u1 <= iso && iso <= u0)){
               p = (iso-u1)/(u0-u1);
               q = 1.-p;
               SET20(start,x0,p,x1,q)
            }
            else{
               p = (iso-u2)/(u0-u2);
               q = 1.-p;
               SET20(start,x0,p,x2,q)
            }
            if ((u1 <= iso && iso <= u2) || (u2 <= iso && iso <= u1)){
               p = (iso-u2)/(u1-u2);
               q = 1.-p;
               SET20(end,x1,p,x2,q)
            }
            else{
               p = (iso-u2)/(u0-u2);
               q = 1.-p;
               SET20(end,x0,p,x2,q)
            }
         }
         fprintf(fp,"%e %e\n",start[0],start[1]);
         fprintf(fp,"%e %e\n\n",end[0],end[1]);
         iso += delta;
      }
}

#if N_DATA & SCALAR_NODE_DATA

void p1_isolines(tGrid,u,name,n)
GRID *tGrid;
char name[];
INT u, n;
{
   ELEMENT *pel;
   NODE *pnode;
   DOUBLE min=NDS(FIRSTN(tGrid),u), max=NDS(FIRSTN(tGrid),u), delta;
   FILE *fp;

   for (pnode=FIRSTN(tGrid); pnode; pnode=SUCC(pnode))
      if (NDS(pnode,u) > max)
         max = NDS(pnode,u);
      else if (NDS(pnode,u) < min)
         min = NDS(pnode,u);
   if (max-min > 1.e-20){
      fp = fopen(name,"w");
      delta = (max-min)/n;
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
         p1_isolines_on_triangle(fp,pel->n[0]->myvertex->x,
             pel->n[1]->myvertex->x,pel->n[2]->myvertex->x,
             NDS(pel->n[0],u),NDS(pel->n[1],u),NDS(pel->n[2],u),0.,delta);
      fclose(fp);
   }
}

#else

void p1_isolines(tGrid,u,name,n)
GRID *tGrid; char name[]; INT u, n;
{ eprintf("Error: p1_isolines not available.\n"); }

#endif

/******************************************************************************/

#if (N_DATA & VECTOR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void save_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"w");
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fprintf(fp,"%e %e\n",ND(pnode,u,0),ND(pnode,u,1));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e\n",FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,ue));
   fclose(fp);
}

void read_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"r");
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fscanf(fp,"%lf %lf\n",&ND(pnode,u,0),&ND(pnode,u,1));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf\n",&FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf\n",&ED(pel,ue));
   fclose(fp);
}

#else

void save_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   eprintf("Error: save_res not available.\n");
}

void read_res(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   eprintf("Error: read_res not available.\n");
}

#endif

#endif  /* DIM == 2 */


#if (N_DATA & SCALAR_NODE_DATA) && (DATA_STR & LG_DATA) && (DIM == 3)
 
void save_results(tGrid,u,ue,ut,t,name)
GRID *tGrid;
INT u, ue, ut, t;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;
   
   fp = fopen(name,"w");
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (pnode->lgd)
         fprintf(fp,"%e %e %e\n",ND(pnode,u,0),ND(pnode,u,1),ND(pnode,u,2));
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fprintf(fp,"%e %e %e\n",ND(pnode,u,0),ND(pnode,u,1),ND(pnode,u,2));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e\n",FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,ue));
   if (t){
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
         if (!IS_DIR(pnode,t))
            fprintf(fp,"%e\n",NDS(pnode,ut));
      for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
         fprintf(fp,"%e\n",NDS(pnode,ut));
   }
   fclose(fp);
} 

void read_results(tGrid,u,ue,ut,t,name)
GRID *tGrid;
INT u, ue, ut, t;
char name[];
{
   ELEMENT *pel;
   NODE *pnode;
   FACE *pface;
   FILE *fp;
   
   fp = fopen(name,"r");
   for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
      if (pnode->lgd){
        fscanf(fp,"%lf %lf %lf\n",&ND(pnode,u,0),&ND(pnode,u,1),&ND(pnode,u,2));
        NDLG(pnode,u,0) = DOT(NDD(pnode,u),pnode->lgd->t1);
        NDLG(pnode,u,1) = DOT(NDD(pnode,u),pnode->lgd->t2);
      }
   for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
      fscanf(fp,"%lf %lf %lf\n",&ND(pnode,u,0),&ND(pnode,u,1),&ND(pnode,u,2));
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf\n",&FD(pface,u));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf\n",&ED(pel,ue));
   if (t){
      for (pnode = FDBN(tGrid); pnode != FIRSTNODE(tGrid); pnode = pnode->succ)
         if (!IS_DIR(pnode,t))
            fscanf(fp,"%lf\n",&NDS(pnode,ut));
      for (pnode = FIRSTNODE(tGrid); pnode != NULL; pnode = pnode->succ)
         fscanf(fp,"%lf\n",&NDS(pnode,ut));
   }
   fclose(fp);
}

#endif

#if (F_DATA & VECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (DIM==2)

void save_res_p1nc_p0(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;
 
   fp = fopen(name,"w");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e %e\n",FDV(pface,u,0),FDV(pface,u,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,ue));
   fclose(fp);
}

void read_res_p1nc_p0(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"r");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf %lf\n",&FDV(pface,u,0),&FDV(pface,u,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf\n",&ED(pel,ue));
   fclose(fp);
}

#else  /*  !((F_DATA & VECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)
                                         && (DIM==2))                         */
void save_res_p1nc_p0(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: save_res_p1nc_p0 not available.\n");  }

void read_res_p1nc_p0(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: read_res_p1nc_p0 not available.\n");  }

#endif

#if (F_DATA & DVECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA) && (DIM==2)

void save_res_p1mod_p0(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;
 
   fp = fopen(name,"w");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e %e %e %e\n",FDDV(pface,u,0,0),FDDV(pface,u,0,1),
                                 FDDV(pface,u,1,0),FDDV(pface,u,1,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,ue));
   fclose(fp);
}

void read_res_p1mod_p0(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"r");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf %lf %lf %lf\n",&FDDV(pface,u,0,0),&FDDV(pface,u,0,1),
                                    &FDDV(pface,u,1,0),&FDDV(pface,u,1,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf\n",&ED(pel,ue));
   fclose(fp);
}

#else  /*  !((F_DATA & DVECTOR_FACE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)
                                          && (DIM==2))                        */
void save_res_p1mod_p0(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: save_res_p1mod_p0 not available.\n");  }

void read_res_p1mod_p0(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: read_res_p1mod_p0 not available.\n");  }

#endif

#if (F_DATA & DVECTOR_FACE_DATA) && (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (DIM==2)

void save_res_p1mod_p1disc(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;
 
   fp = fopen(name,"w");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fprintf(fp,"%e %e %e %e\n",FDDV(pface,u,0,0),FDDV(pface,u,0,1),
                                 FDDV(pface,u,1,0),FDDV(pface,u,1,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fprintf(fp,"%e %e\n",EDSN(pel,ue,0),EDSN(pel,ue,1));
   fclose(fp);
}

void read_res_p1mod_p1disc(tGrid,u,ue,name)
GRID *tGrid;
INT u, ue;
char name[];
{
   ELEMENT *pel;
   FACE *pface;
   FILE *fp;

   fp = fopen(name,"r");
   for (pface = FIRSTFACE(tGrid); pface != NULL; pface = pface->succ)
      fscanf(fp,"%lf %lf %lf %lf\n",&FDDV(pface,u,0,0),&FDDV(pface,u,0,1),
                                    &FDDV(pface,u,1,0),&FDDV(pface,u,1,1));
   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
      fscanf(fp,"%lf %lf\n",&EDSN(pel,ue,0),&EDSN(pel,ue,1));
   fclose(fp);
}

#else  /*  !((F_DATA & DVECTOR_FACE_DATA) && 
             (E_DATA & SCALAR_DATA_IN_ELEMENT_NODES) && (DIM==2))             */

void save_res_p1mod_p1disc(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: save_res_p1mod_p1disc not available.\n");  }

void read_res_p1mod_p1disc(tGrid,u,ue,name)
GRID *tGrid; INT u, ue; char name[];
{  eprintf("Error: read_res_p1mod_p1disc not available.\n");  }

#endif

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) && (F_DATA & SCALAR_FACE_DATA) && (DIM == 2)

void save_p1nc_matrix_and_vectors(tGrid,Z,u,f)
GRID *tGrid;
INT Z, u, f;
{
   FACE *theFace;
   FLINK *pfl;
   INT ind_j[5], i=1, j, l=1, m, max, ind;
   DOUBLE aij[5], a;
   FILE *fp_ira, *fp_dia, *fp_jsl, *fp_u, *fp_f, *fp_m;

   fp_ira = fopen("DIA_IRA_JSL/ira","w");
   fp_dia = fopen("DIA_IRA_JSL/dia","w");
   fp_jsl = fopen("DIA_IRA_JSL/jsl","w");
   fp_u = fopen("DIA_IRA_JSL/vector_u","w");
   fp_f = fopen("DIA_IRA_JSL/vector_f","w");
   fp_m = fopen("DIA_IRA_JSL/aij","w");
   fprintf(fp_ira,"%i\n",1);
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      theFace->index2 = i++;
      fprintf(fp_u,"%e\n",FD(theFace,u));
      fprintf(fp_f,"%e\n",FD(theFace,f));
   }
   fclose(fp_u);
   fclose(fp_f);

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      fprintf(fp_m,"%i %i %e\n",
                           theFace->index2,theFace->index2,COEFF_FF(theFace,Z));
      aij[m=0] = COEFF_FF(theFace,Z); 
      ind_j[m] = theFace->index2; 
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl)){
         fprintf(fp_m,"%i %i %e\n",
                           theFace->index2,NBFACE(pfl)->index2,COEFF_FL(pfl,Z));
         aij[++m] = COEFF_FL(pfl,Z); 
         ind_j[m] = NBFACE(pfl)->index2;
      } 
      m++; 
      for (i = m; i > 1; i--){
         max = 0; 
         for (j = 0; j < i; j++)
            if (ind_j[j] > max){
               max = ind_j[j]; 
               ind = j;
            }
         EXCHANGE(ind_j[i-1],ind_j[ind],j) 
         EXCHANGE(aij[i-1],aij[ind],a)
      }
      fprintf(fp_ira,"%i\n",(l+=m));
      for (i = 0; i < m; i++){
         fprintf(fp_dia,"%e\n",aij[i]);
         fprintf(fp_jsl,"%i\n",ind_j[i]);
      }
   }
   fclose(fp_ira);
   fclose(fp_dia);
   fclose(fp_jsl);
}

#else

void save_p1nc_matrix_and_vectors(tGrid,Z,u,f)
GRID *tGrid; INT Z, u, f;
{ eprintf("Error: save_p1nc_matrix_and_vectors not available.\n"); }

#endif

void new_node_indices(tGrid,fn,ni,n)
INT *ni, *n;
GRID *tGrid;
{
   INT i=fn;
   NODE *theNode, *stop=FIRSTNODE(tGrid);

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
      theNode->index2 = i++;
   *ni = i-1;
   for (theNode=FIRSTN(tGrid); theNode!=stop; theNode=SUCC(theNode))
      theNode->index2 = i++;
   *n = i-1;
}

void new_face_indices(tGrid,ni,n)
INT *ni, *n;
GRID *tGrid;
{
   INT i=0;
   FACE *theFace, *stop=FIRSTFACE(tGrid);

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
      theFace->index2 = ++i;
   *ni = i;
   for (theFace=FIRSTF(tGrid); theFace!=stop; theFace=SUCC(theFace))
      theFace->index2 = ++i;
   *n = i;
}

void new_element_indices(tGrid,nf,n)
INT nf, *n;
GRID *tGrid;
{
   INT i=nf;
   ELEMENT *theElement;

   for (theElement=FIRSTELEMENT(tGrid); theElement; theElement=SUCC(theElement))
      theElement->index2 = ++i;
   *n = i;
}

#if (F_DATA & ONE_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) 

void save_f_matrix(tGrid,Z,m1,m2,fp)
GRID *tGrid;
INT Z, m1, m2;
FILE *fp;
{
   FACE *theFace;
   FLINK *pfl;

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      fprintf(fp,"%i %i %e\n",
                    m1+theFace->index2,m2+theFace->index2,COEFF_FF(theFace,Z));
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",
                    m1+theFace->index2,m2+NBFACE(pfl)->index2,COEFF_FL(pfl,Z));
   }
}

#else 

void save_f_matrix(tGrid,Z,m1,m2,fp)
GRID *tGrid; INT Z, m1, m2; FILE *fp;
{  eprintf("Error: save_f_matrix not available.\n");  }

#endif

#if (F_DATA & DxD_FACE_MATR) && (DATA_S & F_LINK_TO_FACES) 

void save_vf_matrix(tGrid,Z,m1,m2,nfi,fp)
GRID *tGrid;
INT Z, m1, m2, nfi;
FILE *fp;
{
   FACE *theFace;
   FLINK *pfl;
   INT n1=m1+nfi, n2=m2+nfi;

   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      fprintf(fp,"%i %i %e\n",
                  m1+theFace->index2,m2+theFace->index2,COEFFNN(theFace,Z,0,0));
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",
                  m1+theFace->index2,m2+NBFACE(pfl)->index2,COEFFNN(pfl,Z,0,0));
      fprintf(fp,"%i %i %e\n",
                  m1+theFace->index2,n2+theFace->index2,COEFFNN(theFace,Z,0,1));
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",
                  m1+theFace->index2,n2+NBFACE(pfl)->index2,COEFFNN(pfl,Z,0,1));
   }
   for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace)){
      fprintf(fp,"%i %i %e\n",
                  n1+theFace->index2,m2+theFace->index2,COEFFNN(theFace,Z,1,0));
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",
                  n1+theFace->index2,m2+NBFACE(pfl)->index2,COEFFNN(pfl,Z,1,0));
      fprintf(fp,"%i %i %e\n",
                  n1+theFace->index2,n2+theFace->index2,COEFFNN(theFace,Z,1,1));
      for (pfl=FSTART(theFace); pfl; pfl=NEXT(pfl))
         fprintf(fp,"%i %i %e\n",
                  n1+theFace->index2,n2+NBFACE(pfl)->index2,COEFFNN(pfl,Z,1,1));
   }
}

#else 

void save_vf_matrix(tGrid,Z,m1,m2,nfi,fp)
GRID *tGrid; INT Z, m1, m2, nfi; FILE *fp;
{  eprintf("Error: save_vf_matrix not available.\n");  }

#endif

#if (E_DATA & ExDF_MATR) && (ELEMENT_TYPE == SIMPLEX)

void save_e_X_vf_matrix(tGrid,Z,m1,m2,nfi,fp)
GRID *tGrid;
INT Z, m1, m2, nfi;
FILE *fp;
{
   ELEMENT *pel;
   INT k, n2=m2+nfi;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      k = m1 + pel->index2;
      if (IS_FF(pel->f[0]))
         fprintf(fp,"%i %i %e\n",k,m2+pel->f[0]->index2,COEFF_BDF(pel,Z,0,0));
      if (IS_FF(pel->f[1]))
         fprintf(fp,"%i %i %e\n",k,m2+pel->f[1]->index2,COEFF_BDF(pel,Z,1,0));
      if (IS_FF(pel->f[2]))
         fprintf(fp,"%i %i %e\n",k,m2+pel->f[2]->index2,COEFF_BDF(pel,Z,2,0));
      if (IS_FF(pel->f[0]))
         fprintf(fp,"%i %i %e\n",k,n2+pel->f[0]->index2,COEFF_BDF(pel,Z,0,1));
      if (IS_FF(pel->f[1]))
         fprintf(fp,"%i %i %e\n",k,n2+pel->f[1]->index2,COEFF_BDF(pel,Z,1,1));
      if (IS_FF(pel->f[2]))
         fprintf(fp,"%i %i %e\n",k,n2+pel->f[2]->index2,COEFF_BDF(pel,Z,2,1));
   }
}

#else

void save_e_X_vf_matrix(tGrid,Z,m1,m2,nfi,fp)
GRID *tGrid; INT Z, m1, m2, nfi; FILE *fp;
{  eprintf("Error: save_e_X_vf_matrix not available.\n");  }

#endif

void save_Stokes_matrices_and_vectors(tGrid,ZA,ZB,f,fe,u,ue,save_vectors,
                                      u_type,p_type,A_struct,B_struct)
GRID *tGrid;
INT ZA, ZB, f, fe, u, ue, u_type, p_type, A_struct, B_struct, save_vectors;
{
   INT nni, nn, nfi, nf, ne;
   FILE *fp_a, *fp_b;

   new_node_indices(tGrid,1,&nni,&nn);
   new_face_indices(tGrid,&nfi,&nf);
   new_element_indices(tGrid,0,&ne);
   if (save_vectors){
      save_vector(tGrid,u ,0,nn,nni,nf,nfi,ne,0,u_type,"MATRICES/vector_u");   
      save_vector(tGrid,ue,0,nn,nni,nf,nfi,ne,0,p_type,"MATRICES/vector_p");   
      save_vector(tGrid,f ,0,nn,nni,nf,nfi,ne,0,u_type,"MATRICES/vector_f");   
      save_vector(tGrid,fe,0,nn,nni,nf,nfi,ne,0,p_type,"MATRICES/vector_g");   
   }
   fp_a = fopen("MATRICES/matrix_a","w");
   fp_b = fopen("MATRICES/matrix_b","w");
   switch(u_type){
   case Q_VF: 
        if (p_type == Q_SE && (A_struct & (Q_FULL | Q_FEBDIAG))){
           fprintf(fp_a,"%i\n",2*nfi);
           fprintf(fp_b,"%i %i\n",ne,2*nfi);
           if (A_struct & Q_FULL)
              save_vf_matrix(tGrid,ZA,0,0,nfi,fp_a);
           else if (A_struct & Q_FEBDIAG){
              save_f_matrix(tGrid,ZA,0,0,fp_a);
              save_f_matrix(tGrid,ZA,nfi,nfi,fp_a);
           }
           save_e_X_vf_matrix(tGrid,ZB,0,0,nfi,fp_b);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
/*
   case Q_DVF: 
        if ((p_type == Q_SE || p_type == Q_SNE) && (A_struct & Q_DVFEBDIAG)){
           vmultA_dvf(tGrid,Z,x,y,t_row);
           if (p_type == Q_SE && (B_struct & Q_FIRST_DV))
              multA_dvf_X_e(tGrid,Z,x,y,q,t_row);
           else if (p_type == Q_SNE)
              multA_dvf_X_sef(tGrid,Z,x,y,t_row);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
   case Q_VNSF: 
        if ((p_type == Q_SE) && (A_struct & (Q_FULL | Q_NEBDIAG))){
           if (A_struct & Q_FULL)
              vmultA_vn_sf(tGrid,Z,x,y);
           else if (structure & Q_NEBDIAG)
              smultA_vn_sf(tGrid,Z,x,y,t_row);
           multA_vn_sf_X_e(tGrid,Z,x,y,q,t_row);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
   case Q_VNVF: 
        if ((p_type == Q_SE || p_type == Q_SN) && 
            (A_struct & (Q_FULL | Q_EBDIAG))){
           if (A_struct & Q_FULL)
              vmultA_vn_vf(tGrid,Z,x,y,t_row);
           else if (A_struct & Q_EBDIAG)
              smultA_vn_vf(tGrid,Z,x,y,t_row);
           if (p_type == Q_SN)
              multA_vn_vf_X_sn(tGrid,Z,x,y);
           else if (p_type == Q_SE)
              multA_vn_vf_X_e(tGrid,Z,x,y,q,t_row);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
   case Q_VNVE:
        if (p_type == Q_SN && (A_struct & Q_NEBDIAG_EDIAG){
           smultA_vn(tGrid,Z,x,y,t_row);
           multA_vn_X_sn(tGrid,Z,x,y,t_row);
           multA_ve_X_sn(tGrid,Z,x,y,t_row);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
   case Q_VNVFVE: 
        if ((p_type == Q_SE || p_type == Q_SN || p_type == Q_SNE) &&
            (A_struct & Q_EBDIAG)){
           smultA_vn_vf_ve(tGrid,Z,x,y,q,t_row);
           if (p_type == Q_SE){
              multA_vn_vf_X_e(tGrid,Z,x,y,q,t_row);
              set_value(tGrid,0.,y,t_column,Q_VE);
           }
           else if (p_type == Q_SN){
              multA_vn_X_sn(tGrid,Z,x,y,t_row);
              multA_ve_X_sn(tGrid,Z,x,y,t_row);
           }
           else if (p_type == Q_SNE)
              multA_vn_vf_ve_X_sef(tGrid,Z,x,y,t_row);
        }
        else
           eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
*/
   default:
        eprintf("Error: save_Stokes_matrices_and_vectors not available.\n");
        break;
   }
   fprintf(fp_a,"%i %i %e\n",-1,-1,-1.);
   fprintf(fp_b,"%i %i %e\n",-1,-1,-1.);
   fclose(fp_a);
   fclose(fp_b);
}

void umf_mult_A(Ap,Ai,Ax,x,y,nj)
INT *Ap, *Ai, nj;
DOUBLE *Ax, *x, *y;
{
   NODE *pni;
   LINK *plij;
   INT i, j, k, m;

   for (i = 0; i < nj; i++)
      y[i] = 0.;
   for (j = 0; j < nj; j++){
      m = Ap[j+1];
      for (k = Ap[j]; k < m; k++)
         y[Ai[k]] += Ax[k]*x[j];
   }
}

void put_index_into_Ai_list(ind,j,n,Ap,Ai,max)
INT ind, j, *n, *Ap, *Ai, max;
{
   INT i, k, l;

   i = (*n)++;
   while (i >= Ap[j] && ind < Ai[i]) i--;
   i++;
   if (*n >= max)
      eprintf("Error: maximum number of entries in Ai exceeded.\n");
   k = *n;
   while ((l=k) > i)
      Ai[l] = Ai[--k];
   Ai[i] = ind;
}

#if DATA_S & N_LINK_TO_NODES 

void fill_Ap_and_Ai_for_one_node_matr(tGrid,Ap,Ai,nj,maxj,max)
GRID *tGrid;
INT *Ap, *Ai, *nj, maxj, max;
{
   NODE *pnj;
   LINK *plji;
   INT i, j, k, l, n;

   new_node_indices(tGrid,0,nj,&n);
   (*nj)++;
   if (*nj >= maxj)
      eprintf("Error: Ap too short.\n");
   Ap[0] = 0;
   for (pnj = FIRSTNODE(tGrid); pnj; pnj = pnj->succ){
      j = pnj->index2;
      n = Ap[j];
      Ai[n] = j;
      for (plji = START(pnj); plji; plji = plji->next)
         put_index_into_Ai_list(NBNODE(plji)->index2,j,&n,Ap,Ai,max);
      Ap[j+1] = n + 1;
   }
/*
   printf("nj = %i\n",*nj);
   printf("Ap:\n");
   for (j=0; j <= *nj; j++)
      printf("Ap[%i] = %i, ",j,Ap[j]);
   printf("\nAi:");
   for (j=0; j < *nj; j++){
      printf("\nj = %i:  ",j);
      for (k = Ap[j]; k < Ap[j+1]; k++)
         printf("%i ",Ai[k]);
   }
   printf("\n====================\n");
   for (pnj = FIRSTNODE(tGrid); pnj; pnj = pnj->succ){
      j = pnj->index2;
      printf("(%i,%i) ",j,j);
      for (plji = START(pnj); plji; plji = plji->next)
         printf("(%i,%i) ",j,NBNODE(plji)->index2);
      printf("\n");
   }
*/
}

void fill_Ap_and_Ai_for_sn_se(tGrid,Ap,Ai,nj,maxj,max)
GRID *tGrid;
INT *Ap, *Ai, *nj, maxj, max;
{
   ELEMENT *pel;
   NODE *pnj;
   LINK *plji;
   INT i, j, k, l, n;

   new_node_indices(tGrid,0,&j,&n);
   new_element_indices(tGrid,j,nj);
   (*nj)++;
   if (*nj >= maxj)
      eprintf("Error: Ap too short.\n");
   Ap[0] = 0;
   for (pnj = FIRSTNODE(tGrid); pnj; pnj = pnj->succ){
      j = pnj->index2;
      n = Ap[j];
      Ai[n] = j;
      for (plji = START(pnj); plji; plji = plji->next)
         put_index_into_Ai_list(NBNODE(plji)->index2,j,&n,Ap,Ai,max);
      for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
         for (k = 0; k < SIDES; k++)
            if (pel->n[k]->index2 == j)
               put_index_into_Ai_list(pel->index2,j,&n,Ap,Ai,max);
      Ap[j+1] = n + 1;
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      j = pel->index2;
      n = Ap[j];
      Ai[n] = j;
      for (k = 0; k < SIDES; k++)
         if (IS_FN(pel->n[k]))
            put_index_into_Ai_list(pel->n[k]->index2,j,&n,Ap,Ai,max);
      Ap[j+1] = n + 1;
   }
}

#else

void fill_Ap_and_Ai_for_one_node_matr(tGrid,Ap,Ai,nj,maxj,max)
GRID *tGrid; INT *Ap, *Ai, *nj, maxj, max;
{ eprintf("Error: fill_Ap_and_Ai_for_one_node_matr not available.\n"); }

void fill_Ap_and_Ai_for_sn_se(tGrid,Ap,Ai,nj,maxj,max)
GRID *tGrid; INT *Ap, *Ai, *nj, maxj, max;
{ eprintf("Error: fill_Ap_and_Ai_for_sn_se not available.\n"); }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) 

void fill_Ax_for_one_node_matr(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid;
INT Z, *Ap, *Ai;
DOUBLE *Ax;
{
   NODE *pni;
   LINK *plij;
   INT i, j, k;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
      i = pni->index2;
      for (k = Ap[i]; Ai[k] != i; k++);
      Ax[k] = COEFFN(pni,Z);
      for (plij = START(pni); plij; plij = plij->next){
         j = NBNODE(plij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFFL(plij,Z);
      }
   }
}

void test_Ax_for_one_node_matr(tGrid,Z,Ap,Ai,Ax,nj)
GRID *tGrid;
INT Z, *Ap, *Ai, nj;
DOUBLE *Ax;
{
   NODE *pni;
   LINK *plij;
   INT i, j, k;

   for (j=0; j < nj; j++){
      for (k = Ap[j]; k < Ap[j+1]; k++)
         printf("(%i,%i) %e, ",Ai[k],j,Ax[k]);
      printf("\n");
   }
   printf("\n====================\n");
   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
      i = pni->index2;
      printf("(%i,%i) %e, ",i,i,COEFFN(pni,Z));
      for (plij = START(pni); plij; plij = plij->next)
         printf("(%i,%i) %e, ",i,NBNODE(plij)->index2,COEFFL(plij,Z));
      printf("\n");
   }
}

#else

void fill_Ax_for_one_node_matr(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid; INT Z, *Ap, *Ai; DOUBLE *Ax;
{ eprintf("Error: fill_Ax_for_one_node_matr not available.\n"); }

void test_Ax_for_one_node_matr(tGrid,Z,Ap,Ai,Ax,nj)
GRID *tGrid; INT Z, *Ap, *Ai, nj; DOUBLE *Ax;
{ eprintf("Error: test_Ax_for_one_node_matr not available.\n"); }

#endif

#if (N_DATA & ONE_NODE_MATR) && (N_DATA & ONE_NODE_FACE_MATR) && (F_DATA & ONE_FACE_MATR) && (F_DATA & ONE_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) 

void fill_Ax_for_sn_sf(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid;
INT Z, *Ap, *Ai;
DOUBLE *Ax;
{
   NODE *pni;
   FACE *pfi;
   LINK *plij;
   NFLINK *pnflij;
   FNLINK *pfnlij;
   FLINK *pflij;
   INT i, j, k;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
      i = pni->index2;
      for (k = Ap[i]; Ai[k] != i; k++);
      Ax[k] = COEFFN(pni,Z);
      for (plij = START(pni); plij; plij = plij->next){
         j = NBNODE(plij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFFL(plij,Z);
      }
      for (pnflij = NFSTART(pni); pnflij; pnflij = pnflij->next){
         j = NBFACE(pnflij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFFL(pnflij,Z);
      }
   }
   for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ){
      i = pfi->index2;
      for (k = Ap[i]; Ai[k] != i; k++);
      Ax[k] = COEFF_FF(pfi,Z);
      for (pflij = FSTART(pfi); pflij; pflij = pflij->next){
         j = NBFACE(pflij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFF_FL(pflij,Z);
      }
      for (pfnlij = FNSTART(pfi); pfnlij; pfnlij = pfnlij->next){
         j = NBNODE(pfnlij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFFL(pfnlij,Z);
      }
   }
}

#else

void fill_Ax_for_sn_sf(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid; INT Z, *Ap, *Ai; DOUBLE *Ax;
{  eprintf("Error: fill_Ax_for_sn_sf not available.\n");  }

#endif

#if (N_DATA & ONE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (E_DATA & ExN_MATR) && (E_DATA & NxE_MATR) && (E_DATA & ExE_MATR)

void fill_Ax_for_sn_se_matr(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid;
INT Z, *Ap, *Ai;
DOUBLE *Ax;
{
   ELEMENT *pel;
   NODE *pni;
   LINK *plij;
   INT i, j, k, l;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
      i = pni->index2;
      for (k = Ap[i]; Ai[k] != i; k++);
      Ax[k] = COEFFN(pni,Z);
      for (plij = START(pni); plij; plij = plij->next){
         j = NBNODE(plij)->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         Ax[k] = COEFFL(plij,Z);
      }
   }
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ){
      for (l = 0; l < SIDES; l++)
         if (IS_FN(pel->n[l])){
            i = pel->n[l]->index2;
            for (k = Ap[pel->index2]; Ai[k] != i; k++);
            Ax[k] = COEFF_NE(pel,Z,l);
            i = pel->index2;
            for (k = Ap[pel->n[l]->index2]; Ai[k] != i; k++);
            Ax[k] = COEFF_EN(pel,Z,l);
         }
      i = pel->index2;
      for (k = Ap[i]; Ai[k] != i; k++);
      Ax[k] = COEFF_EE(pel,Z);
   }
}

#else

void fill_Ax_for_sn_se_matr(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid; INT Z, *Ap, *Ai; DOUBLE *Ax;
{  eprintf("Error: fill_Ax_for_sn_se_matr not available.\n");  }

#endif

#if N_DATA & SCALAR_NODE_DATA

void make_vector_from_sn(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   NODE *pni;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      x[pni->index2] = NDS(pni,u);
}

void make_sn_from_vector(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   NODE *pni;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      NDS(pni,u) = x[pni->index2];
}

#else

void make_vector_from_sn(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_vector_from_sn not available.\n");  }

void make_sn_from_vector(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_sn_from_vector not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (F_DATA & SCALAR_FACE_DATA)

void make_vector_from_grid_data_sn_sf(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;
   FACE *pfi;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      x[pni->index2] = NDS(pni,u);
   for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ)
      x[pfi->index2] = FD(pfi,u);
}

void make_grid_data_from_vector_sn_sf(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;
   FACE *pfi;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      NDS(pni,u) = x[pni->index2];
   for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ)
      FD(pfi,u) = x[pfi->index2];
}

#else

void make_vector_from_grid_data_sn_sf(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_vector_from_grid_data_sn_sf not available.\n");  }

void make_grid_data_from_vector_sn_sf(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_grid_data_from_vector_sn_sf not available.\n");  }

#endif

#if (N_DATA & SCALAR_NODE_DATA) && (E_DATA & SCALAR_ELEMENT_DATA)

void make_vector_from_sn_se(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      x[pni->index2] = NDS(pni,u);
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      x[pel->index2] = ED(pel,u);
}


void make_sn_se_from_vector(tGrid,u,x)
GRID *tGrid;
INT u;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;

   for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
      NDS(pni,u) = x[pni->index2];
   for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
      ED(pel,u) = x[pel->index2];
}

#else

void make_vector_from_sn_se(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_vector_from_sn_se not available.\n");  }

void make_sn_se_from_vector(tGrid,u,x)
GRID *tGrid; INT u; DOUBLE *x;
{  eprintf("Error: make_sn_se_from_vector not available.\n");  }

#endif

#if (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (DIM == 2)

void put_indices_into_Ai_list(ind,kk,j,n,Ap,Ai,max)
INT ind, kk, j, *n, *Ap, *Ai, max;
{
   INT i=(*n), k=(*n), l;

   l = *n += kk;
   while (i >= Ap[j] && ind < Ai[i]) i--;
   i++;
   if (*n >= max)
      eprintf("Error: maximum number of entries in Ai exceeded.\n");
   else{
      while (k >= i)
         Ai[l--] = Ai[k--];
      Ai[i] = ind;
      for (l = 0; l < kk; l++)
         Ai[i+l] = ind + l;
   }
}

void copy_columns(j,kk,Ap,Ai,max)
INT j, kk, *Ap, *Ai, max;
{
   INT i, k=Ap[j+1]-Ap[j], l, m;

   if (Ap[j] + kk*k >= max)
      eprintf("Error: maximum number of entries in Ai exceeded.\n");
   else{
      for (l = 1; l < kk; l++){
         m = j + l;
         Ap[m+1] = Ap[m] + k;
         m = l*k;
         for (i = Ap[j]; i < Ap[j+1]; i++)
            Ai[i+m] = Ai[i];
      }
   }
}

void fill_Ap_and_Ai_general(tGrid,Ap,Ai,nj,maxj,max,kk,nn,mm)
GRID *tGrid;
INT *Ap, *Ai, *nj, maxj, max, kk, nn, mm;
{
   NODE *theNode, *pnj;
   FACE *theFace, *pfj;
   ELEMENT *pel;
   LINK *plji;
   NFLINK *pnflji;
   FNLINK *pfnlji;
   FLINK *pflji;
   INT i, j, k, l, n;

   i = -kk;
   if (kk)
      for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode))
         theNode->index2 = i += kk;
   i += kk - nn;
   if (nn)
      for (theFace=FIRSTFACE(tGrid); theFace; theFace=SUCC(theFace))
         theFace->index2 = i += nn;
   i += nn - mm;
   if (mm)
      for (pel=FIRSTELEMENT(tGrid); pel; pel=SUCC(pel))
         pel->index2 = i += mm;
   (*nj) = i += mm;
   if (*nj >= maxj)
      eprintf("Error: Ap too short.\n");
   Ap[0] = 0;
   if (kk)
      for (pnj = FIRSTNODE(tGrid); pnj; pnj = pnj->succ){
         j = pnj->index2;
         n = Ap[j];
         for (l = 0; l < kk; l++)
            Ai[n+l] = j + l;
         n += kk - 1;
         for (plji = START(pnj); plji; plji = plji->next)
            put_indices_into_Ai_list(NBNODE(plji)->index2,kk,j,&n,Ap,Ai,max);
         if (nn)
            for (pnflji = NFSTART(pnj); pnflji; pnflji = pnflji->next)
               put_indices_into_Ai_list(NBFACE(pnflji)->index2,nn,j,&n,Ap,Ai,max);
         if (mm)
            for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
               for (k = 0; k < NVERT; k++)
                  if (IS_FN(pel->n[k]) && pel->n[k]->index2 == j)
                     put_indices_into_Ai_list(pel->index2,mm,j,&n,Ap,Ai,max);
         Ap[j+1] = n + 1;
         copy_columns(j,kk,Ap,Ai,max);
      }
   if (nn)
      for (pfj = FIRSTFACE(tGrid); pfj; pfj = pfj->succ){
         j = pfj->index2;
         n = Ap[j];
         for (l = 0; l < nn; l++)
            Ai[n+l] = j + l;
         n += nn - 1;
         if (kk)
            for (pfnlji = FNSTART(pfj); pfnlji; pfnlji = pfnlji->next)
               put_indices_into_Ai_list(NBNODE(pfnlji)->index2,kk,j,&n,Ap,Ai,max);
         for (pflji = FSTART(pfj); pflji; pflji = pflji->next)
            put_indices_into_Ai_list(NBFACE(pflji)->index2,nn,j,&n,Ap,Ai,max);
         if (mm)
            for (pel = FIRSTELEMENT(tGrid); pel; pel=pel->succ)
               for (k = 0; k < SIDES; k++)
                  if (IS_FF(pel->f[k]) && pel->f[k]->index2 == j)
                     put_indices_into_Ai_list(pel->index2,mm,j,&n,Ap,Ai,max);
         Ap[j+1] = n + 1;
         copy_columns(j,nn,Ap,Ai,max);
      }
   if (mm)
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         j = pel->index2;
         n = Ap[j];
         for (l = 0; l < mm; l++)
            Ai[n+l] = j + l;
         n += mm - 1;
         if (kk)
            for (k = 0; k < NVERT; k++)
               if (IS_FN(pel->n[k]))
                  put_indices_into_Ai_list(pel->n[k]->index2,kk,j,&n,Ap,Ai,max);
         if (nn)
            for (k = 0; k < SIDES; k++)
               if (IS_FF(pel->f[k]))
                  put_indices_into_Ai_list(pel->f[k]->index2,nn,j,&n,Ap,Ai,max);
         Ap[j+1] = n + 1;
         copy_columns(j,mm,Ap,Ai,max);
      }
}

#else

void fill_Ap_and_Ai_general(tGrid,Ap,Ai,nj,maxj,max,kk,nn,mm)
GRID *tGrid; INT *Ap, *Ai, *nj, maxj, max, kk, nn, mm;
{ eprintf("Error: fill_Ap_and_Ai_general not available.\n"); }

#endif

#if (N_DATA & NxN_NODE_MATR) && (N_DATA & NxM_NODE_FACE_MATR) && (F_DATA & NxN_FACE_MATR) && (F_DATA & MxN_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & MxM_E_E_MATR) && (E_DATA & MxN_E_N_MATR) && (E_DATA & NxM_N_E_MATR) && (E_DATA & MxN_E_F_MATR) && (E_DATA & NxM_F_E_MATR) && (DIM == 2)

/*
void fill_Ax_for_general_matr(tGrid,Z,Ap,Ai,Ax)
GRID *tGrid;
INT Z, *Ap, *Ai;
DOUBLE *Ax;
{
   NODE *pni;
   LINK *plij;
   INT i, j, k, kk=N_OF_NODE_FUNC, nn=N_OF_FACE_FUNC, mm=N_OF_ELEM_FUNC;

   if (kk > 1)
      eprintf("Error: kk > 1 not implemented in general_stiff_matr.\n");
   else if (nn > 0)
      eprintf("Error: nn > 0 not implemented in general_stiff_matr.\n");
   else if (mm > 0)
      eprintf("Error: mm > 0 not implemented in general_stiff_matr.\n");
   else
      for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
         i = pni->index2;
         for (k = Ap[i]; Ai[k] != i; k++);
         Ax[k] = COEFF_NN(pni,Z,0,0);
         for (plij = START(pni); plij; plij = plij->next){
            j = NBNODE(plij)->index2;
            for (k = Ap[j]; Ai[k] != i; k++);
            Ax[k] = COEFF_NN(plij,Z,0,0);
         }
      }
}
*/

void fill_Ax_for_general_matr(tGrid,Z,Ap,Ai,Ax,kk,nn,mm)
GRID *tGrid;
INT Z, *Ap, *Ai, kk, nn, mm;
DOUBLE *Ax;
{
   ELEMENT *pel;
   NODE *pni;
   FACE *pfi;
   LINK *plij;
   NFLINK *pnflij;
   FNLINK *pfnlij;
   FLINK *pflij;
   INT i, j, k, l, li, lj, m, r;

   if (kk)
      for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ){
         i = j = pni->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         m = Ap[j+1] - Ap[j];
         for (lj = 0; lj < kk; lj++){
            r = k+lj*m;
            for (li = 0; li < kk; li++)
               Ax[r+li] = COEFF_NN(pni,Z,li,lj);
         }
         for (plij = START(pni); plij; plij = plij->next){
            j = NBNODE(plij)->index2;
            for (k = Ap[j]; Ai[k] != i; k++);
            m = Ap[j+1] - Ap[j];
            for (lj = 0; lj < kk; lj++){
               r = k+lj*m;
               for (li = 0; li < kk; li++)
                  Ax[r+li] = COEFF_NN(plij,Z,li,lj);
            }
         }
         if (nn)
            for (pnflij = NFSTART(pni); pnflij; pnflij = pnflij->next){
               j = NBFACE(pnflij)->index2;
               for (k = Ap[j]; Ai[k] != i; k++);
               m = Ap[j+1] - Ap[j];
               for (lj = 0; lj < nn; lj++){
                  r = k+lj*m;
                  for (li = 0; li < kk; li++)
                     Ax[r+li] = COEFF_NN(pnflij,Z,li,lj);
               }
            }
      }
   if (nn)
      for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ){
         i = j = pfi->index2;
         for (k = Ap[j]; Ai[k] != i; k++);
         m = Ap[j+1] - Ap[j];
         for (lj = 0; lj < nn; lj++){
            r = k+lj*m;
            for (li = 0; li < nn; li++)
               Ax[r+li] = COEFF_NN(pfi,Z,li,lj);
         }
         for (pflij = FSTART(pfi); pflij; pflij = pflij->next){
            j = NBFACE(pflij)->index2;
            for (k = Ap[j]; Ai[k] != i; k++);
            m = Ap[j+1] - Ap[j];
            for (lj = 0; lj < nn; lj++){
               r = k+lj*m;
               for (li = 0; li < nn; li++)
                  Ax[r+li] = COEFF_NN(pflij,Z,li,lj);
            }
         }
         if (kk)
            for (pfnlij = FNSTART(pfi); pfnlij; pfnlij = pfnlij->next){
               j = NBNODE(pfnlij)->index2;
               for (k = Ap[j]; Ai[k] != i; k++);
               m = Ap[j+1] - Ap[j];
               for (lj = 0; lj < kk; lj++){
                  r = k+lj*m;
                  for (li = 0; li < nn; li++)
                     Ax[r+li] = COEFF_NN(pfnlij,Z,li,lj);
               }
            }
      }
   if (mm)
      for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
         i = j = pel->index2;
         m = Ap[j+1] - Ap[j];
         for (k = Ap[j]; Ai[k] != i; k++);
         for (lj = 0; lj < mm; lj++){
            r = k+lj*m;
            for (li = 0; li < mm; li++)
               Ax[r+li] = COEFF_EE_MM(pel,Z,li,lj);
         }
         if(kk)
            for (l = 0; l < NVERT; l++)
               if (IS_FN(pel->n[l])){
                  i = pel->n[l]->index2;
                  for (k = Ap[j]; Ai[k] != i; k++);
                  for (lj = 0; lj < mm; lj++){
                     r = k+lj*m;
                     for (li = 0; li < kk; li++)
                        Ax[r+li] = COEFF_NE_NM(pel,Z,l,li,lj);
                  }
               }
         if(nn)
            for (l = 0; l < SIDES; l++)
               if (IS_FF(pel->f[l])){
                  i = pel->f[l]->index2;
                  for (k = Ap[j]; Ai[k] != i; k++);
                  for (lj = 0; lj < mm; lj++){
                     r = k+lj*m;
                     for (li = 0; li < nn; li++)
                        Ax[r+li] = COEFF_FE_NM(pel,Z,l,li,lj);
                  }
               }
         i = pel->index2;
         if(kk)
            for (l = 0; l < NVERT; l++)
               if (IS_FN(pel->n[l])){
                  j = pel->n[l]->index2;
                  for (k = Ap[j]; Ai[k] != i; k++);
                  m = Ap[j+1] - Ap[j];
                  for (lj = 0; lj < kk; lj++){
                     r = k+lj*m;
                     for (li = 0; li < mm; li++)
                        Ax[r+li] = COEFF_EN_MN(pel,Z,l,li,lj);
                  }
               }
         if(nn)
            for (l = 0; l < SIDES; l++)
               if (IS_FF(pel->f[l])){
                  j = pel->f[l]->index2;
                  for (k = Ap[j]; Ai[k] != i; k++);
                  m = Ap[j+1] - Ap[j];
                  for (lj = 0; lj < nn; lj++){
                     r = k+lj*m;
                     for (li = 0; li < mm; li++)
                        Ax[r+li] = COEFF_EF_MN(pel,Z,l,li,lj);
                  }
               }
      }
}

#else

void fill_Ax_for_general_matr(tGrid,Z,Ap,Ai,Ax,kk,nn,mm)
GRID *tGrid; INT Z, *Ap, *Ai, kk, nn, mm; DOUBLE *Ax;
{ eprintf("Error: fill_Ax_for_general_matr not available.\n"); }

#endif

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA)

void make_vector_from_grid_data_general(tGrid,u,x,kk,nn,mm)
GRID *tGrid;
INT u, kk, nn, mm;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;
   FACE *pfi;
   INT l;

   if (kk)
      for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
         for (l = 0; l < kk; l++)
            x[pni->index2 + l] = NDMV(pni,u,l);
   if (nn)
      for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ)
         for (l = 0; l < nn; l++)
            x[pfi->index2 + l] = FDMV(pfi,u,l);
   if (mm)
      for (pel = FIRSTELEMENT(tGrid); pel; pel=SUCC(pel))
         for (l = 0; l < mm; l++)
            x[pel->index2 + l] = EDMV(pel,u,l);
}

void make_grid_data_from_vector_general(tGrid,u,x,kk,nn,mm)
GRID *tGrid;
INT u, kk, nn, mm;
DOUBLE *x;
{
   ELEMENT *pel;
   NODE *pni;
   FACE *pfi;
   INT l;

   if (kk)
      for (pni = FIRSTNODE(tGrid); pni; pni = pni->succ)
         for (l = 0; l < kk; l++)
            NDMV(pni,u,l) = x[pni->index2 + l];
   if (nn)
      for (pfi = FIRSTFACE(tGrid); pfi; pfi = pfi->succ)
         for (l = 0; l < nn; l++)
            FDMV(pfi,u,l) = x[pfi->index2 + l];
   if (mm)
      for (pel = FIRSTELEMENT(tGrid); pel; pel=SUCC(pel))
         for (l = 0; l < mm; l++)
            EDMV(pel,u,l) = x[pel->index2 + l];
}

#else

void make_vector_from_grid_data_general(tGrid,u,x,kk,nn,mm)
GRID *tGrid; INT u, kk, nn, mm; DOUBLE *x;
{  eprintf("Error: make_vector_from_grid_data_general not available.\n");  }

void make_grid_data_from_vector_general(tGrid,u,x,kk,nn,mm)
GRID *tGrid; INT u, kk, nn, mm; DOUBLE *x;
{  eprintf("Error: make_grid_data_from_vector_general not available.\n");  }

#endif

void fill_Ap_and_Ai(tGrid,Ap,Ai,nj,maxj,max,type,space)
GRID *tGrid;
INT *Ap, *Ai, *nj, maxj, max, type, space;
{
   switch(type){
   case Q_SN:    fill_Ap_and_Ai_for_one_node_matr(tGrid,Ap,Ai,nj,maxj,max);
        break;
   case Q_SNSF:  fill_Ap_and_Ai_general(tGrid,Ap,Ai,nj,maxj,max,1,1,0);
        break;
   case Q_SNSE:  fill_Ap_and_Ai_for_sn_se(tGrid,Ap,Ai,nj,maxj,max);
        break;
   case Q_GENERAL: fill_Ap_and_Ai_general(tGrid,Ap,Ai,nj,maxj,max,
                             N_OF_NODE_FUNC,N_OF_FACE_FUNC,N_OF_ELEM_FUNC);
        break;
   default:
        eprintf("Error: fill_Ap_and_Ai not available.\n");
        break;
   }
}

void fill_Ax(tGrid,Z,Ap,Ai,Ax,type,space)
GRID *tGrid;
INT Z, *Ap, *Ai, type, space;
DOUBLE *Ax;
{
   switch(type){
   case Q_SN:    fill_Ax_for_one_node_matr(tGrid,Z,Ap,Ai,Ax);
        break;
   case Q_SNSF:  fill_Ax_for_sn_sf(tGrid,Z,Ap,Ai,Ax);
        break;
   case Q_SNSE:  fill_Ax_for_sn_se_matr(tGrid,Z,Ap,Ai,Ax);
        break;
   case Q_GENERAL: fill_Ax_for_general_matr(tGrid,Z,Ap,Ai,Ax,
                             N_OF_NODE_FUNC,N_OF_FACE_FUNC,N_OF_ELEM_FUNC);
        break;
   default:
        eprintf("Error: fill_Ax not available.\n");
        break;
   }
}

void make_vector_from_grid_data(tGrid,u,x,type,space)
GRID *tGrid;
INT u, type, space;
DOUBLE *x;
{
   switch(type){
   case Q_SN:    make_vector_from_sn(tGrid,u,x);
        break;
   case Q_SNSF:  make_vector_from_grid_data_sn_sf(tGrid,u,x);
        break;
   case Q_SNSE:  make_vector_from_sn_se(tGrid,u,x);
        break;
   case Q_GENERAL: make_vector_from_grid_data_general(tGrid,u,x,
                                  N_OF_NODE_FUNC,N_OF_FACE_FUNC,N_OF_ELEM_FUNC);
        break;
   default:
        eprintf("Error: make_vector_from_grid_data not available.\n");
        break;
   }
}

void make_grid_data_from_vector(tGrid,u,x,type,space)
GRID *tGrid;
INT u, type, space;
DOUBLE *x;
{
   switch(type){
   case Q_SN:    make_sn_from_vector(tGrid,u,x);
        break;
   case Q_SNSF:  make_grid_data_from_vector_sn_sf(tGrid,u,x);
        break;
   case Q_SNSE:  make_sn_se_from_vector(tGrid,u,x);
        break;
   case Q_GENERAL: make_grid_data_from_vector_general(tGrid,u,x,
                                  N_OF_NODE_FUNC,N_OF_FACE_FUNC,N_OF_ELEM_FUNC);
        break;
   default:
        eprintf("Error: make_grid_data_from_vector not available.\n");
        break;
   }
}

void set_min_and_max(min,max,e)  /*  assume *min <= *max  */
FLOAT *min, *max, e;
{
   if (e > *max)
      *max = e;
   else if (e < *min)
      *min = e;
}

#if DIM == 2

#if N_DATA & SCALAR_NODE_DATA

void check_maximum_principle_for_p1c(tGrid,u,tmin,tmax,rmin,rmax,l2tmin,l2tmax,
                                                                 l2rmin,l2rmax)
GRID *tGrid;
INT u;
FLOAT *tmin, *tmax, *rmin, *rmax, *l2tmin, *l2tmax, *l2rmin, *l2rmax;
{
   NODE *pnode;
   LINK *pli;
   FLOAT min, max, val;

   *tmin = *tmax = *rmin = *rmax = *l2tmin = *l2tmax = *l2rmin = *l2rmax = 0.;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ){
      max = min = NDS(pnode->tstart->nbnode,u) - pnode->tstart->nbnode->myvertex->x[0];
      for (pli=pnode->tstart->next; pli; pli=pli->next){
         val = NDS(pli->nbnode,u) - pli->nbnode->myvertex->x[0];
         min = MIN(min,val);
         max = MAX(max,val);
      }
      val = NDS(pnode,u) - pnode->myvertex->x[0];
      min = MAX(min-val,0.);
      max = MAX(val-max,0.);
      *tmin = MAX(*tmin,min);
      *tmax = MAX(*tmax,max);
      *l2tmin += min*min;
      *l2tmax += max*max;
      val = fabs(val);
      if (val > 1.e-12){
         min /= val;
         max /= val;
      }
      *rmin = MAX(*rmin,min);
      *rmax = MAX(*rmax,max);
      *l2rmin += min*min;
      *l2rmax += max*max;
   }
}

void p1c_L_inf_errors(tGrid,u0,u,err0p,err1p,err2p,err3p,
                                 err0m,err1m,err2m,err3m,errl2m,errl2ma)
GRID *tGrid;
FLOAT (*u0)(), *err0p, *err1p, *err2p, *err3p, *err0m, *err1m, *err2m, *err3m,
      *errl2m, *errl2ma;
INT u;
{
   GRID *theGrid;
   NODE *theNode;
   FLOAT d, *x;
   INT n = 0;

   *err0p = *err1p = *err2p = *err3p = *err0m = *err1m = *err2m = *err3m = 0.;
   *errl2m = 0.;
   for (theNode=FIRSTN(tGrid); theNode; theNode=SUCC(theNode)){
      x = MYVERTEX(theNode)->x;
      d = NDS(theNode,u) - (*u0)(x);
      if (d < 0.){
         *errl2m += d*d;
         n++;
      }
      if (x[0] < 0.9){
         if (x[1] > 0.1 && x[1] < 0.9)
            set_min_and_max(err0m,err0p,d);
         else
            set_min_and_max(err1m,err1p,d);
      }
      else{
         if (x[1] > 0.1 && x[1] < 0.9)
            set_min_and_max(err3m,err3p,d);
         else
            set_min_and_max(err2m,err2p,d);
      }
   }
   if (n)
      *errl2ma = sqrt(*errl2m/n);
   *errl2m = sqrt(*errl2m);
   printf("+ L-infinity errors: err0 = %9.4e  err1 = %9.4e  err2 = %9.4e  err3 = %9.4e\n",*err0p, *err1p, *err2p, *err3p);
   printf("- L-infinity errors: err0 = %9.4e  err1 = %9.4e  err2 = %9.4e  err3 = %9.4e\n",*err0m, *err1m, *err2m, *err3m);
   printf("- l2-errors: errl2m = %9.4e  errl2ma = %9.4e\n",*errl2m,*errl2ma);
}

void p1c_smear_and_osc(tGrid,u,smear_2,smear_inf,in_osc_2,in_osc_inf,
                                                  b_osc_2, b_osc_inf)
GRID *tGrid;
FLOAT *smear_2, *smear_inf, *in_osc_2, *in_osc_inf, *b_osc_2, *b_osc_inf;
INT u;
{
   GRID *theGrid;
   NODE *theNode;
   FLOAT val, *x, a0, a1, a2, max0=0., max1=0.;

   *smear_2 = *smear_inf = *b_osc_2 = *b_osc_inf= *in_osc_2 = 0.;
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      x = theNode->myvertex->x;
      val = NDS(theNode,u);
      a0 = -MIN(val,0.);
      a1 = MAX(val-1.,0.);
      a2 = MAX(1.-val,0.);
      if (x[0] >= 0.7){
         *smear_2 += a2*a2;
         *smear_inf = MAX(*smear_inf,a2);
//         if (x[0] >= 0.9){
         if (x[0] >= 0.7){
            *b_osc_2 += a1*a1;
            *b_osc_inf = MAX(*b_osc_inf,a1);
         }
      }
//      else if (x[0] <= 0.5){
      else if (x[0] <= 0.5 && x[1] >= 0.1){
         *in_osc_2 += a0*a0 + a1*a1;
         max0 = MAX(max0,a0);
         max1 = MAX(max1,a1);
      }
   }
   *in_osc_inf = max0 + max1;
   *smear_2 = sqrt(*smear_2);
   *b_osc_2 = sqrt(*b_osc_2);
   *in_osc_2 = sqrt(*in_osc_2);
   printf("smear:   l2: %e   max: %e\n",*smear_2,*smear_inf);
   printf("in osc:  l2: %e   max: %e\n",*in_osc_2,*in_osc_inf);
   printf("b osc:   l2: %e   max: %e\n",*b_osc_2,*b_osc_inf);
}

void min_and_max_on_cut_p1c(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{
   NODE *pnode;
   DOUBLE middle;

   *min =  1.e10;
   *max = -1.e10;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
      if (fabs(pnode->myvertex->x[0]-0.5) < 1.e-10){
         *min = MIN(*min,NDS(pnode,U));
         *max = MAX(*max,NDS(pnode,U));
         if (fabs(pnode->myvertex->x[1]-0.5) < 1.e-10)
            middle = NDS(pnode,U);
      }
   *min -= middle;
   *max -= middle;
}

void compute_last_point_err_p1c(mg,u,err)
MULTIGRID *mg;
FLOAT *err;
INT u;
{
   GRID *theGrid;
   NODE *pnode;
   FLOAT x=-1., value;

   for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
      for (pnode = FIRSTNODE(theGrid); pnode; pnode=pnode->succ)
         if (IS_TOP_NODE(pnode))
            if ((pnode->myvertex->x[0] > x) &&
                (fabs(pnode->myvertex->x[1] - 0.5) < 1.e-14)){
               x = pnode->myvertex->x[0];
               value = NDS(pnode,u);
            }
   *err = value - x;
}

void min_and_max_for_p1c(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{
   NODE *pnode;

   *min = *max = NDS(FIRSTNODE(tGrid),U);
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      *min = MIN(*min,NDS(pnode,U));
      *max = MAX(*max,NDS(pnode,U));
   }
}

void min_and_max_for_p1c_and_hill(tGrid,min,diff)
GRID *tGrid;
FLOAT *min, *diff;
{
   NODE *pnode;
   FLOAT *x, min1, max1;

   *min = min1 = 100.;
   max1 = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      if (x[0] >= 0.8){
         min1 = MIN(min1,NDS(pnode,U));
         max1 = MAX(max1,NDS(pnode,U));
      }
      else if (x[0] >= 0.4 && x[0] <= 0.6)
         *min = MIN(*min,NDS(pnode,U));
   }
   *min = fabs(*min);
   *diff = max1 - min1;
}

void min_and_max_for_p1c_and_hemker(tGrid,min1,max1,min2,max2)
GRID *tGrid;
FLOAT *min1, *max1, *min2, *max2;
{
   NODE *pnode;
   FLOAT *x;

   *min1 = *min2 =  100.;
   *max1 = *max2 = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      if (x[0] > -1.5 && x[0] < 1.5){
         *min1 = MIN(*min1,NDS(pnode,U));
         *max1 = MAX(*max1,NDS(pnode,U));
      }
      else if (x[0] > 2.){
         *min2 = MIN(*min2,NDS(pnode,U));
         *max2 = MAX(*max2,NDS(pnode,U));
      }
   }
}

void sn_min_and_max_for_KLR(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{
   NODE *pnode;
   FLOAT *x;

   *min =  100.;
   *max = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      *min = MIN(*min,NDS(pnode,U));
      *max = MAX(*max,NDS(pnode,U));
   }
}

#else

void check_maximum_principle_for_p1c(tGrid,u,tmin,tmax,rmin,rmax,l2tmin,l2tmax,l2rmin,l2rmax)
GRID *tGrid; INT u; FLOAT *tmin, *tmax, *rmin, *rmax, *l2tmin, *l2tmax, *l2rmin, *l2rmax;
{  eprintf("Error: check_maximum_principle_for_p1c not available.\n");  }

void p1c_L_inf_errors(tGrid,u0,u,err0p,err1p,err2p,err3p,err0m,err1m,err2m,err3m,errl2m,errl2ma)
GRID *tGrid; FLOAT (*u0)(), *err0p, *err1p, *err2p, *err3p, *err0m, *err1m, *err2m, *err3m, *errl2m, *errl2ma; INT u;
{  eprintf("Error: p1c_L_inf_errors not available.\n");  }

void p1c_smear_and_osc(tGrid,u,smear_2,smear_inf,in_osc_2,in_osc_inf,b_osc_2, b_osc_inf)
GRID *tGrid; FLOAT *smear_2, *smear_inf, *in_osc_2, *in_osc_inf, *b_osc_2, *b_osc_inf; INT u;
{  eprintf("Error: p1c_smear_and_osc not available.\n");  }

void min_and_max_on_cut_p1c(tGrid,min,max)
GRID *tGrid; FLOAT *min, *max;
{  eprintf("Error: min_and_max_on_cut_p1c not available.\n");  }

void compute_last_point_err_p1c(mg,u,err)
MULTIGRID *mg; FLOAT *err; INT u;
{  eprintf("Error: compute_last_point_err_p1c not available.\n");  }

void min_and_max_for_p1c(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{  eprintf("Error: min_and_max_for_p1c not available.\n");  }

void min_and_max_for_p1c_and_hill(tGrid,min,diff)
GRID *tGrid; FLOAT *min, *diff;
{  eprintf("Error: min_and_max_for_p1c_and_hill not available.\n");  }

void min_and_max_for_p1c_and_hemker(tGrid,min1,max1,min2,max2)
GRID *tGrid; FLOAT *min1, *max1, *min2, *max2;
{  eprintf("Error: min_and_max_for_p1c_and_hemker not available.\n");  }

void sn_min_and_max_for_KLR(tGrid,min,max)
GRID *tGrid; FLOAT *min, *max;
{  eprintf("Error: sn_min_and_max_for_KLR not available.\n");  }

#endif

void gmin_and_max_for_KLR2(tGrid,min,max,min2,max2,min_out,max_out,mid_val,space)
GRID *tGrid;
FLOAT *min, *max, *min2, *max2, *min_out, *max_out, *mid_val;
INT space;
{
   NODE *pnode;
   FLOAT *x, val, s, diff=1.;

   *min =  100.;
   *max = -100.;
   *min_out =  100.;
   *max_out = -100.;
   *mid_val = -100.;
   *min2 = 0.;
   *max2 = 0.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      val = sn_node_value(pnode,U,space);
      *min = MIN(*min,val);
      *max = MAX(*max,val);
      s = MIN(val,0.);
      *min2 += s*s;
      s = MAX(val-1.,0.);
      *max2 += s*s;
      if (x[0] < 1.e-7){
         *min_out = MIN(*min_out,val);
         *max_out = MAX(*max_out,val);
         if (fabs(x[1]-0.5) < diff){
            *mid_val = val;
            diff = fabs(x[1]-0.5);
         }
      }
   }
   *min2 = sqrt(*min2);
   *max2 = sqrt(*max2);
}

void evaluate_profile(y,v,n,maxo)
DOUBLE *y, *v, *maxo;
INT n;
{
   DOUBLE s;
   int i=0, j, imin, imax, imid=0;

   while (v[++i] < 0.1); 
   while (v[i] < v[i+1]) i++; 
   imin = i;
   i = n; 
   while (v[--i] < 0.1); 
   while (v[i] < v[i-1]) i--; 
   imax = i;
   s = 0.5*(y[imin] + y[imax]);
   for (i = imin; i <= imax; i++)
      if (fabs(y[i]-s) < fabs(y[imid]-s))
         imid = i;
   s = 0.;
   for (i = imin; i <= imid; i++)
      for (j = i; j <= imid; j++)
         if (v[i] - v[j] > s)
            s = v[i] - v[j];
   for (i = imid; i <= imax; i++)
      for (j = imid; j <= i; j++)
         if (v[i] - v[j] > s)
            s = v[i] - v[j];
/*
      printf("zacatek: %e %e\n",y[imin],v[imin]);
      printf("stred: %e %e\n",y[imid],v[imid]);
      printf("konec: %e %e\n",y[imax],v[imax]);
      printf("maxo = %e\n",s);
*/
   *maxo = s;
}

void add_point_to_1d_graph();

void compute_maxo_for_klr(tGrid,u,space,maxo)
GRID *tGrid;
DOUBLE *maxo;
INT u, space;
{
   NODE *pnode;
   DOUBLE y[1000], v[1000], y1[1000], v1[1000], min, max;
   INT i, j, m, n=0;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      add_point_to_1d_graph(pnode->myvertex->x[1],pnode->myvertex->x[0],
                            sn_node_value(pnode,u,space),0.,y1,v1,&n);
   max = y1[0];
   for (i = 1; i < n; i++)
      if (y1[i] > max)
         max = y1[i];
   max = 10.*max + 10.;
   for (i = 0; i < n; i++){
      min = y1[0];
      m = 0;
      for (j = 1; j < n; j++)
         if (y1[j] < min){
            min = y1[j];
            m = j;
         }
      y[i] = y1[m];
      v[i] = v1[m];
      y1[m] = max;
   }
   evaluate_profile(y,v,n,maxo);
}

#if (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void min_and_max_on_cut_gp1c(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{
   NODE *pnode;
   DOUBLE middle;

   *min =  1.e10;
   *max = -1.e10;
   for (pnode = FIRSTNODE(tGrid); pnode; pnode = pnode->succ)
      if (fabs(pnode->myvertex->x[0]-0.5) < 1.e-10){
         *min = MIN(*min,NDMV(pnode,U,0));
         *max = MAX(*max,NDMV(pnode,U,0));
         if (fabs(pnode->myvertex->x[1]-0.5) < 1.e-10)
            middle = NDMV(pnode,U,0);
      }
   *min -= middle;
   *max -= middle;
}

void gmin_and_max_for_KLR(tGrid,min,max)
GRID *tGrid;
FLOAT *min, *max;
{
   NODE *pnode;
   FLOAT *x;

   *min =  100.;
   *max = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      *min = MIN(*min,NDMV(pnode,U,0));
      *max = MAX(*max,NDMV(pnode,U,0));
   }
}

void min_and_max_for_gp1c_and_hill(tGrid,min,diff)
GRID *tGrid;
FLOAT *min, *diff;
{
   NODE *pnode;
   FLOAT *x, min1, max1;

   *min = min1 = 100.;
   max1 = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      if (x[0] >= 0.8){
         min1 = MIN(min1,NDMV(pnode,U,0));
         max1 = MAX(max1,NDMV(pnode,U,0));
      }
      else if (x[0] >= 0.4 && x[0] <= 0.6)
         *min = MIN(*min,NDMV(pnode,U,0));
   }
   *min = fabs(*min);
   *diff = max1 - min1;
}

void min_and_max_for_gp1c_and_hemker(tGrid,min1,max1,min2,max2)
GRID *tGrid;
FLOAT *min1, *max1, *min2, *max2;
{
   NODE *pnode;
   FLOAT *x;

   *min1 = *min2 =  100.;
   *max1 = *max2 = -100.;
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ){
      x = pnode->myvertex->x;
      if (x[0] > -1.5 && x[0] < 1.5){
         *min1 = MIN(*min1,NDMV(pnode,U,0));
         *max1 = MAX(*max1,NDMV(pnode,U,0));
      }
      else if (x[0] > 2.){
         *min2 = MIN(*min2,NDMV(pnode,U,0));
         *max2 = MAX(*max2,NDMV(pnode,U,0));
      }
   }
}

void nodes_of_edges();
void save_1d_graph();
INT not_among_points();

void save_gq1x4c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid;
FLOAT d;
INT k, l, u, space;
char name[];
{
   ELEMENT *pel;
   NODE *pnode, *n1[4], *n2[4];
   FACE *f[4];
   DOUBLE x[DIM], y[1000], v[1000];
   INT i, n=0;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      add_point_to_1d_graph(pnode->myvertex->x[l],
                            pnode->myvertex->x[k],NDMV(pnode,u,0),d,y,v,&n);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      nodes_of_edges(pel,f,n1,n2);
      for (i = 0; i < EDGES; i++){
         AVERAGE(n1[i]->myvertex->x,n2[i]->myvertex->x,x)
         if (not_among_points(x[l],x[k],d,y,n))
            add_point_to_1d_graph(x[l],x[k],FDMV(f[i],u,0),d,y,v,&n);
      }
   }
   save_1d_graph(y,v,n,name);
}

void save_gq2x4c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid;
FLOAT d;
INT k, l, u, space;
char name[];
{
   ELEMENT *pel;
   NODE *pnode, *n1[4], *n2[4];
   FACE *f[4];
   DOUBLE x[DIM], x1[DIM], x2[DIM], v1, v2, y[1000], v[1000];
   INT i, n=0;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      add_point_to_1d_graph(pnode->myvertex->x[l],
                            pnode->myvertex->x[k],NDMV(pnode,u,0),d,y,v,&n);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      nodes_of_edges(pel,f,n1,n2);
      for (i = 0; i < EDGES; i++){
         AVERAGE(n1[i]->myvertex->x,n2[i]->myvertex->x,x)
         if (not_among_points(x[l],x[k],d,y,n)){
            add_point_to_1d_graph(x[l],x[k],FDMV(f[i],u,1),d,y,v,&n);
            AVERAGE(n1[i]->myvertex->x,x,x1)
            AVERAGE(n2[i]->myvertex->x,x,x2)
            if (n1[i]->index < n2[i]->index){
               v1 = FDMV(pel->f[i],u,0) + 0.5*(NDMV(n1[i],u,0) + FDMV(pel->f[i],u,1));
               v2 = FDMV(pel->f[i],u,2) + 0.5*(NDMV(n2[i],u,0) + FDMV(pel->f[i],u,1));
            }
            else{
               v1 = FDMV(pel->f[i],u,2) + 0.5*(NDMV(n1[i],u,0) + FDMV(pel->f[i],u,1));
               v2 = FDMV(pel->f[i],u,0) + 0.5*(NDMV(n2[i],u,0) + FDMV(pel->f[i],u,1));
            }
            add_point_to_1d_graph(x1[l],x1[k],v1,d,y,v,&n);
            add_point_to_1d_graph(x2[l],x2[k],v2,d,y,v,&n);
         }
      }
   }
   save_1d_graph(y,v,n,name);
}
#else

void min_and_max_on_cut_gp1c(tGrid,min,max)
GRID *tGrid; FLOAT *min, *max;
{  eprintf("Error: min_and_max_on_cut_gp1c not available.\n");  }

void gmin_and_max_for_KLR(tGrid,min,max)
GRID *tGrid; FLOAT *min, *max;
{  eprintf("Error: gmin_and_max_for_KLR not available.\n");  }

void min_and_max_for_gp1c_and_hill(tGrid,min,diff)
GRID *tGrid; FLOAT *min, *diff;
{  eprintf("Error: min_and_max_for_gp1c_and_hill not available.\n");  }

void min_and_max_for_gp1c_and_hemker(tGrid,min1,max1,min2,max2)
GRID *tGrid; FLOAT *min1, *max1, *min2, *max2;
{  eprintf("Error: min_and_max_for_gp1c_and_hemker not available.\n");  }

void save_gq1x4c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid; FLOAT d; INT k, l, u, space; char name[];
{  eprintf("Error: save_gq1x4c_profile not available.\n");  }

void save_gq2x4c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid; FLOAT d; INT k, l, u, space; char name[];
{  eprintf("Error: save_gq2x4c_profile not available.\n");  }

#endif

void p1c_smear_and_osc2(tGrid,u,space,smear_2,in_osc_2,b_osc_2)
GRID *tGrid;
FLOAT *smear_2, *in_osc_2, *b_osc_2;
INT u, space;
{
   GRID *theGrid;
   NODE *theNode; 
   FLOAT val, *x, a0, a1, a2, max0=0., max1=0.;

   *smear_2 = *b_osc_2 = *in_osc_2 = 0.;
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      x = theNode->myvertex->x;
      val = sn_node_value(theNode,u,space);
      a0 = -MIN(val,0.);
      a1 = MAX(val-1.,0.);
      a2 = MAX(1.-val,0.);
      if (x[0] >= 0.7){
         *smear_2 += a2*a2;
         if (x[0] >= 0.7){
            *b_osc_2 += a1*a1;
         }
      }
      else if (x[0] <= 0.5 && x[1] >= 0.1){
         *in_osc_2 += a0*a0 + a1*a1;
         max0 = MAX(max0,a0);
         max1 = MAX(max1,a1);
      }
   }
   *smear_2 = sqrt(*smear_2);
   *b_osc_2 = sqrt(*b_osc_2);
   *in_osc_2 = sqrt(*in_osc_2);
   printf("in osc:  l2: %e\n",*in_osc_2);
   printf("b osc:   l2: %e\n",*b_osc_2);
   printf("smear:   l2: %e\n",*smear_2);
}

void p1c_smear_and_osc2_ref(tGrid,u,space,smear_2b,smear_2r,in_osc_2,b_osc_2b,b_osc_2r) 
GRID *tGrid;
FLOAT *smear_2b, *smear_2r, *in_osc_2, *b_osc_2b, *b_osc_2r;
INT u, space;
{
   GRID *theGrid;
   NODE *theNode; 
   FLOAT val, *x, a0, a1, a2, max0=0., max1=0.;

   *smear_2b = *smear_2r = *b_osc_2b = *b_osc_2r = *in_osc_2 = 0.;
   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      x = theNode->myvertex->x;
      val = sn_node_value(theNode,u,space);
      a0 = -MIN(val,0.);
      a1 = MAX(val-1.,0.);
      a2 = MAX(1.-val,0.);
      if (x[0] >= 0.7){
         if (x[0] < 0.85 && x[1] < 0.15){
            *smear_2b += a2*a2;
            *b_osc_2b += a1*a1;
         }
         else if (x[0] > 0.85 && x[1] > 0.15){
            *smear_2r += a2*a2;
            *b_osc_2r += a1*a1;
         }
      }
      else if (x[0] <= 0.5 && x[1] >= 0.1){
         *in_osc_2 += a0*a0 + a1*a1;
         max0 = MAX(max0,a0);
         max1 = MAX(max1,a1);
      }
   }
   *smear_2b = sqrt(*smear_2b);
   *smear_2r = sqrt(*smear_2r);
   *b_osc_2b = sqrt(*b_osc_2b);
   *b_osc_2r = sqrt(*b_osc_2r);
   *in_osc_2 = sqrt(*in_osc_2);
   printf("in osc:  l2: %e\n",*in_osc_2);
   printf("b osc b: l2: %e\n",*b_osc_2b);
   printf("b osc r: l2: %e\n",*b_osc_2r);
   printf("smear b: l2: %e\n",*smear_2b);
   printf("smear r: l2: %e\n",*smear_2r);
}

void compute_derivatives(tGrid,u,space,d,
                         derx_max0,derx_min0,dery_max0,dery_min0,
                         derx_max1,derx_min1,dery_max1,dery_min1,
                         derx_max2,derx_min2,dery_max2,dery_min2)
GRID *tGrid;
INT u, space;
FLOAT d, *derx_max0, *derx_min0, *dery_max0, *dery_min0,
         *derx_max1, *derx_min1, *dery_max1, *dery_min1,
         *derx_max2, *derx_min2, *dery_max2, *dery_min2;
{
   ELEMENT *pelem;
   INT d0=1, d1=1, d2=1;
   FLOAT *x0, *x1, *x2, b[DIM2][DIM2], xc[DIM], grad_u[DIM];

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
      barycentric_coordinates(x0,x1,x2,b);
      coord_of_barycentre(pelem,xc);
      sgrad_value(pelem,xc,b,u,space,grad_u);
      if (xc[0] < 1.-d){
         if (xc[1] < d){
            if (d1){
               *derx_min1 = *derx_max1 = grad_u[0];
               *dery_min1 = *dery_max1 = grad_u[1];
               d1 = 0;
            } 
            set_min_and_max(derx_min1,derx_max1,grad_u[0]);
            set_min_and_max(dery_min1,dery_max1,grad_u[1]);
         }
         else if (xc[1] > 1.-d){
            if (d1){
               *derx_min1 = *derx_max1 = grad_u[0];
               *dery_min1 = *dery_max1 = -grad_u[1];
               d1 = 0;
            } 
            set_min_and_max(derx_min1,derx_max1,grad_u[0]);
            set_min_and_max(dery_min1,dery_max1,-grad_u[1]);
         }
         else{
            if (d0){
               *derx_min0 = *derx_max0 = grad_u[0];
               *dery_min0 = *dery_max0 = grad_u[1];
               d0 = 0;
            } 
            set_min_and_max(derx_min0,derx_max0,grad_u[0]);
            set_min_and_max(dery_min0,dery_max0,grad_u[1]);
         }
      }
      else if (xc[1] > d && xc[1] < 1.-d){
         if (d2){
            *derx_min2 = *derx_max2 = grad_u[0];
            *dery_min2 = *dery_max2 = grad_u[1];
            d2 = 0;
         } 
         set_min_and_max(derx_min2,derx_max2,grad_u[0]);
         set_min_and_max(dery_min2,dery_max2,grad_u[1]);
      }
   }
}

INT non_zero_array(f,n)
INT *f, n;
{
   INT i=0, j=0;

   while (i < n && !j)
      if (f[i++]) j = 1;
   return(j);
}

#define N_OF_POINTS 100020

void new_operate_points(mg,x,n,u,space,new_operation,data_i,data_f,fp)
MULTIGRID *mg;
FLOAT x[N_OF_POINTS][2], *data_f;
INT n, u, space, *data_i;
FILE *fp;
void (*new_operation)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT value, eps=0., x_ref[2];
   INT i, j=1, f[N_OF_POINTS];
   
   if (n > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in operate_points.\n");
   else{   
      for (i = 0; i < n; i++)
         f[i] = 1;
      while (non_zero_array(f,n)){
         if (eps < -1.e-8)
            eprintf("Error: too large eps in operate_points.\n");
         for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
            for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
               if (IS_TOP_ELEMENT(pel)){
                  for (i = 0; i < n; i++)
                     if (f[i] && is_in_element(x[i],pel,eps,x_ref)){
                        new_operation(theGrid,pel,x[i],x_ref,u,space,data_i,data_f,fp);
                        f[i] = 0;
                     }
               }
         if (j){
            eps = -1.e-15;
            j = 0;
         }
         else
            eps*= 10.;
      }
   }
}

void operate_points(mg,x,n,u,space,operation,data_i,data_f,fp)
MULTIGRID *mg;
FLOAT x[N_OF_POINTS][2], *data_f;
INT n, u, space, *data_i;
FILE *fp;
void (*operation)();
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT *x0, *x1, *x2, value, b[DIM2][DIM2], eps=0.;
   INT i, j=1, f[N_OF_POINTS];
   
   if (n > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in operate_points.\n");
   else{   
      for (i = 0; i < n; i++)
         f[i] = 1;
      while (non_zero_array(f,n)){
         if (eps < -1.e-8)
            eprintf("Error: too large eps in operate_points.\n");
         for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
            for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
               if (IS_TOP_ELEMENT(pel)){
                  VERTICES_OF_ELEMENT(x0,x1,x2,pel);
                  barycentric_coordinates(x0,x1,x2,b);
                  for (i = 0; i < n; i++)
                     if (f[i] && IS_IN_SIMPLEX(b,x[i],eps)){
                        operation(pel,x[i],b,u,space,data_i,data_f,fp);
                        f[i] = 0;
                     }
               }
         if (j){
            eps = -1.e-15;
            j = 0;
         }
         else
            eps*= 10.;
      }
   }
}

void min_max_and_value(pel,x,b,u,space,data_i,data_f,fp)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], *data_f;
INT u, space, *data_i;
FILE *fp;
{
   FLOAT value=sfunc_value(pel,x,b,u,space);

   data_f[0] = MIN(data_f[0],value);
   data_f[1] = MAX(data_f[1],value);
   if ((fabs(x[0]-data_f[3]) < 1.e-14) && (fabs(x[1]-data_f[4]) < 1.e-14))
      data_f[2] = value;
}

void compute_min_max_and_value(mg,u,space,npoints,min,max,value)
MULTIGRID *mg;
FLOAT *min, *max, *value;
INT u, space, npoints;
{
   INT i, n=0;
   FLOAT d=1./(NV_POINTS-1.)-1.e-6, h=1./(npoints-1.), y, x[N_OF_POINTS][2], data_f[5];

   if (npoints > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in compute_min_max_and_value.\n");
   else{
      for (i = 0; i < npoints; i++){
         y = i*h;
         if (y >= d && y <= 1.-d){
            x[n][0] = 0.5;
            x[n++][1] = y;
         }
      }
      data_f[0] =  1.e10;
      data_f[1] = -1.e10;
      data_f[2] =  1.e10;
      data_f[3] = 0.5;
      data_f[4] = 0.5;
      operate_points(mg,x,n,u,space,min_max_and_value,NULL,data_f,NULL);
      if (data_f[2] > 9.e9)
         eprintf("Error: value probably not found.\n");
      *min = data_f[0] - data_f[2];
      *max = data_f[1] - data_f[2];
      *value = data_f[2] - 0.5;
   }
}

void smearing_diff(pel,x,b,u,space,data_i,data_f,fp)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], *data_f;
INT u, space, *data_i;
FILE *fp;
{
   FLOAT value=sfunc_value(pel,x,b,u,space);

   if (x[0] < data_f[0] && value >= 0.1)
      data_f[0] = x[0];
   if (x[0] < data_f[1] && value > 0.9)
      data_f[1] = x[0];
   if (x[0] > data_f[2] && value >= 0.9)
      data_f[2] = x[0];
}

void new_smearing_diff(tGrid,pel,x,x_ref,u,space,data_i,data_f,fp)
GRID *tGrid;
ELEMENT *pel;
FLOAT *x, *x_ref, *data_f;
INT u, space, *data_i;
FILE *fp;
{
   FLOAT value=gfunc_value_ref(tGrid,pel,u,x_ref,ELEM);

   if (x[0] < data_f[0] && value >= 0.1)
      data_f[0] = x[0];
   if (x[0] < data_f[1] && value > 0.9)
      data_f[1] = x[0];
   if (x[0] > data_f[2] && value >= 0.9)
      data_f[2] = x[0];
}

void smearing_diff2(pel,x,b,u,space,data_i,data_f,fp)
ELEMENT *pel;
FLOAT *x, b[DIM2][DIM2], *data_f;
INT u, space, *data_i;
FILE *fp;
{
   FLOAT value=sfunc_value(pel,x,b,u,space);

   if (x[1] < data_f[0] && value >= 0.1)
      data_f[0] = x[1];
   if (x[1] < data_f[1] && value > 0.9)
      data_f[1] = x[1];
   if (x[1] > data_f[2] && value > 0.9)
      data_f[2] = x[1];
   if (x[1] > data_f[3] && value >= 0.1)
      data_f[3] = x[1];
   if (fabs(x[1]-0.5) < data_f[4]){
      data_f[4] = fabs(x[1]-0.5);
      data_f[5] = value;
   }
}

void new_smearing_diff2(tGrid,pel,x,x_ref,u,space,data_i,data_f,fp)
GRID *tGrid;
ELEMENT *pel;
FLOAT *x, *x_ref, *data_f;
INT u, space, *data_i;
FILE *fp;
{
   FLOAT value=gfunc_value_ref(tGrid,pel,u,x_ref,ELEM);

   if (x[1] < data_f[0] && value >= 0.1)
      data_f[0] = x[1];
   if (x[1] < data_f[1] && value > 0.9)
      data_f[1] = x[1];
   if (x[1] > data_f[2] && value > 0.9)
      data_f[2] = x[1];
   if (x[1] > data_f[3] && value >= 0.1)
      data_f[3] = x[1];
   if (fabs(x[1]-0.5) < data_f[4]){
      data_f[4] = fabs(x[1]-0.5);
      data_f[5] = value;
   }
}

void evaluate_smearing(mg,u,space,npoints,diff,point)
MULTIGRID *mg;
FLOAT *diff, *point;
INT u, space, npoints;
{
   INT i;
   FLOAT h=1./(npoints-1.), x[N_OF_POINTS][2], data_f[3];

   if (npoints > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in evaluate_smearing.\n");
   else{
      for (i = 0; i < npoints; i++){
         x[i][0] = i*h;
         x[i][1] = 0.25;
      }
      data_f[0] = 1.;
      data_f[1] = 1.;
      data_f[2] = 0.;
      if (space == P1C || space == Q1C)
         operate_points(mg,x,npoints,u,space,smearing_diff,NULL,data_f,NULL);
      else if (space == GP1C || space == GQ1C)
         new_operate_points(mg,x,npoints,u,space,new_smearing_diff,NULL,data_f,NULL);
      *diff = data_f[1] - data_f[0] - h;
      *point = data_f[2];
printf("points: %e %e\n",data_f[0],data_f[1] - h);
   }
}

void evaluate_smearing1(mg,u,space,npoints,diff,width,mid_val)
MULTIGRID *mg;
FLOAT *diff, *width, *mid_val;
INT u, space, npoints;
{
   INT i;
   FLOAT h=1./(npoints-1.), x[N_OF_POINTS][2], data_f[6];

   if (npoints > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in evaluate_smearing.\n");
   else{
      for (i = 0; i < npoints; i++){
         x[i][0] = 0.;
         x[i][1] = i*h;
      }
      data_f[0] = 1.;
      data_f[1] = 1.;
      data_f[2] = 0.;
      data_f[3] = 0.;
      data_f[4] = 1.;
      data_f[5] = -100.;
      if (space == P1C || space == Q1C)
         operate_points(mg,x,npoints,u,space,smearing_diff2,NULL,data_f,NULL);
      else if (space == GP1C || space == GQ1C)
         new_operate_points(mg,x,npoints,u,space,new_smearing_diff2,NULL,data_f,NULL);
      *diff = data_f[3] - data_f[2] - h + data_f[1] - data_f[0] - h;
      *width = data_f[3] - data_f[0];
      *mid_val = data_f[5];
//printf("%e %e %e %e\n",data_f[0],data_f[1]-h,data_f[2]+h,data_f[3]);
   }
}

void evaluate_smearing2(mg,u,space,npoints,diff,point)
MULTIGRID *mg;
FLOAT *diff, *point;
INT u, space, npoints;
{
   INT i;
   FLOAT h=1./(npoints-1.), x[N_OF_POINTS][2], data_f[4];

   if (npoints > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in evaluate_smearing.\n");
   else{
      for (i = 0; i < npoints; i++){
         x[i][0] = 0.;
         x[i][1] = i*h;
      }
      data_f[0] = 1.;
      data_f[1] = 1.;
      data_f[2] = 0.;
      data_f[3] = 0.;
      if (space == P1C || space == Q1C)
         operate_points(mg,x,npoints,u,space,smearing_diff2,NULL,data_f,NULL);
      else if (space == GP1C || space == GQ1C)
         new_operate_points(mg,x,npoints,u,space,new_smearing_diff2,NULL,data_f,NULL);
      *diff = data_f[3] - data_f[2] - h + data_f[1] - data_f[0] - h;
//printf("%e %e %e %e\n",data_f[0],data_f[1]-h,data_f[2]+h,data_f[3]);
   }
}

void compute_values(mg,x,v,n,u,space)
MULTIGRID *mg;
FLOAT x[N_OF_POINTS][2], *v;
INT n, u, space;
{
   GRID *theGrid;
   ELEMENT *pel;
   FLOAT *x0, *x1, *x2, value, b[DIM2][DIM2], eps=0.;
   INT i, j=1, f[N_OF_POINTS];
   
   if (n > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in compute_value.\n");
   else{
      for (i = 0; i < n; i++)
         f[i] = 1;
      while (non_zero_array(f,n)){
         if (eps < -1.e-8)
            eprintf("Error: two large eps in compute_values.\n");
         for (theGrid = FIRSTGRID(mg); theGrid; theGrid = theGrid->finer)
            for (pel = FIRSTELEMENT(theGrid); pel; pel = pel->succ)
               if (IS_TOP_ELEMENT(pel)){
                  VERTICES_OF_ELEMENT(x0,x1,x2,pel);
                  barycentric_coordinates(x0,x1,x2,b);
                  for (i = 0; i < n; i++)
                     if (f[i] && IS_IN_SIMPLEX(b,x[i],eps)){
                        v[i] = sfunc_value(pel,x[i],b,u,space);
                        f[i] = 0;
                     }
               }
         if (j){
            eps = -1.e-15;
            j = 0;
         }
         else
            eps*= 10.;
      }
   }
}

void integrals_of_errors(mg,u,space,npoints,err1,err2)
MULTIGRID *mg;
FLOAT *err1, *err2;
INT u, space, npoints;
{
   INT i, n;
   FLOAT h=1./(npoints-1.), x[N_OF_POINTS][2], v[N_OF_POINTS], 
         b = 1.-1./(NV_POINTS-1.), s, z=0., u1, u2, v1, v2;

   if (npoints > 4*N_OF_POINTS-8)
      eprintf("Error: value of N_OF_POINTS too small in integrals_of_errors.\n");
   else{
      x[0][0] = 0.75;
      x[0][1] = 0.25;
      for (n = 0; (z=n*h) < 0.75; n++);
      x[1][0] = z;
      x[1][1] = 0.25;
      z -= h;
      for (n = 2; (s=z + n*h) < b; n++){
         x[n][0] = s;
         x[n][1] = 0.25;
      }
      x[n][0] = b;
      x[n++][1] = 0.25;
      compute_values(mg,x,v,n,u,space);
      n--;
      *err1 = *err2 = 0.;
      for (i = 0; i < n; i++){
         u2 = 1. - v[i];
         v2 = 1. - v[i+1];
         u1 = MAX(0,-u2);
         v1 = MAX(0,-v2);
         u2 = MAX(0,u2);
         v2 = MAX(0,v2);
         h = x[i+1][0]-x[i][0];
         *err1 += (u1 + v1)*h;
         *err2 += (u2 + v2)*h;
      }
      *err1 *= 0.5;
      *err2 *= 0.5;
   }
}

void nodes_of_edges(pel,f,n1,n2)
ELEMENT *pel;
NODE **n1, **n2;
FACE **f;
{
   f[0] = pel->f[0];
   f[1] = pel->f[1];
   f[2] = pel->f[2];
   if (ELEMENT_TYPE == SIMPLEX){
      n1[0] = pel->n[1];
      n1[1] = pel->n[2];
      n1[2] = pel->n[0];
      n2[0] = pel->n[2];
      n2[1] = pel->n[0];
      n2[2] = pel->n[1];
   }
   else if (ELEMENT_TYPE == CUBE){
      f[3] = pel->f[3];
      n1[0] = pel->n[0];
      n1[1] = pel->n[1];
      n1[2] = pel->n[2];
      n1[3] = pel->n[3];
      n2[0] = pel->n[1];
      n2[1] = pel->n[2];
      n2[2] = pel->n[3];
      n2[3] = pel->n[0];
   }
}

void graph_on_cut(mg,k,d,u,name,space,npoints)
MULTIGRID *mg;
FLOAT d;
INT k, u, space, npoints;
char name[];
{
   INT i, l=1;
   FLOAT h=1./(npoints-1.), x[N_OF_POINTS][2], v[N_OF_POINTS];
   FILE *fp;

   if (npoints > N_OF_POINTS)
      eprintf("Error: value of N_OF_POINTS too small in graph_on_cut.\n");
   else{
      if (k) l = 0;
      for (i = 0; i < npoints; i++){
         x[i][k] = d;
         x[i][l] = i*h;
      }
      compute_values(mg,x,v,npoints,u,space);
      fp = fopen(name,"w");
      for (i = 0; i < npoints; i++)
         fprintf(fp,"%e %e\n",x[i][l],v[i]);
      fclose(fp);
   }
}

void save_exact_diag_profile(tGrid,u,name)
GRID *tGrid;
FLOAT (*u)();
char name[];
{
   DOUBLE x[DIM], d, h = 0.001;
   FILE *fp;

   fp = fopen(name,"w");
   for (d = 0.; d < 1.; d += h){
      x[0] = x[1] = d;
      fprintf(fp,"%e %e\n",d,u(x));
   }
   x[0] = x[1] = 1.;
   fprintf(fp,"%e %e\n",d,u(x));
   fclose(fp);
}

void add_point_to_1d_graph(x0,x1,value,d,y,v,n)
FLOAT x0, x1, value, d, *y, *v;
INT *n;
{
   if (fabs(x1-d) < 1.e-10){
      y[*n] = x0;
      v[*n] = value;
      (*n)++;
   }
}

INT not_among_points(x0,x1,d,y,n)
FLOAT x0, x1, d, *y;
INT n;
{
   INT i=0;

   if (fabs(x1-d) >= 1.e-10)
      return(1);
   else{
      n--;
      while ((fabs(y[i]-x0) > 1.e-10) && i < n) i++;
      if (fabs(y[i]-x0) > 1.e-10)
         return(1);
      else
         return(0);
   }
}

void save_1d_graph(y,v,n,name)
FLOAT *y, *v;
INT n;
char name[];
{
   DOUBLE miny;
   INT i, mini=n;
   FILE *fp;

   fp = fopen(name,"w");
   do{
      miny = 200.;
      for (i=0; i < n; i++)
         if (y[i] < miny){
            miny = y[i];
            mini = i;
         }
      if (miny < 2.)
         fprintf(fp,"%e %e\n",miny,v[mini]);
      y[mini] = 10.;
   }while (miny < 2. && mini < n); 
   fclose(fp);
}

void save_p1c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid;
FLOAT d;
INT k, l, u, space;
char name[];
{
   NODE *pnode;
   DOUBLE y[1000], v[1000];
   INT n=0;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
/*
   if (IS_FN(pnode) || pnode->myvertex->x[0] < EPSC || 
                       pnode->myvertex->x[1] < EPSC ||
                       fabs(pnode->myvertex->x[0]-1.) < EPSC || 
                       fabs(pnode->myvertex->x[1]-1.) < EPSC ||
                      (fabs(pnode->myvertex->x[0]-0.5) < EPSC && 
                       fabs(pnode->myvertex->x[1]-0.5) < EPSC))
*/
      add_point_to_1d_graph(pnode->myvertex->x[l],
                            pnode->myvertex->x[k],
                            sn_node_value(pnode,u,space),d,y,v,&n);
   save_1d_graph(y,v,n,name);
}

void save_p2c_profile(tGrid,k,l,d,u,name,space)
GRID *tGrid;
FLOAT d;
INT k, l, u, space;
char name[];
{
   ELEMENT *pel;
   NODE *pnode, *n1[4], *n2[4];
   FACE *f[4];
   DOUBLE x[DIM], y[1000], v[1000];
   INT i, n=0;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      add_point_to_1d_graph(pnode->myvertex->x[l],
                            pnode->myvertex->x[k],
                            sn_node_value(pnode,u,space),d,y,v,&n);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      nodes_of_edges(pel,f,n1,n2);
      for (i = 0; i < EDGES; i++){
         AVERAGE(n1[i]->myvertex->x,n2[i]->myvertex->x,x)
         if (not_among_points(x[l],x[k],d,y,n)){
            if (space == P2C || space == GP2C)
               add_point_to_1d_graph(x[l],x[k],0.25*sf_face_value(f[i],u,space) 
                + 0.5*(sn_node_value(n1[i],u,space)+sn_node_value(n2[i],u,space)),d,y,v,&n);
            else if (space == Q2C || space == GQ2C)
               add_point_to_1d_graph(x[l],x[k],sf_face_value(f[i],u,space) 
                + 0.5*(sn_node_value(n1[i],u,space)+sn_node_value(n2[i],u,space)),d,y,v,&n);
         }
      }
   }
   save_1d_graph(y,v,n,name);
}

#if N_DATA & SCALAR_NODE_DATA

void save_p1c_diag_profile(tGrid,u,name)
GRID *tGrid;
INT u;
char name[];
{
   NODE *pnode;
   DOUBLE y[1000], v[1000], miny;
   INT i, n=0, mini;
   FILE *fp;

   fp = fopen(name,"w");
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      if (fabs(pnode->myvertex->x[0]-pnode->myvertex->x[1]) < 1.e-10){
         y[n] = pnode->myvertex->x[0];
         v[n] = NDS(pnode,u);
         n++;
      }
   do{
      miny = 200.;
      for (i=0; i < n; i++)
         if (y[i] < miny){
            miny = y[i];
            mini = i;
         }
      if (miny < 2.)
         fprintf(fp,"%e %e\n",miny,v[mini]);
      y[mini] = 10.;
   }while (miny < 2.); 
   fclose(fp);
}

void save_p1c_diag2_profile(tGrid,u,name)
GRID *tGrid;
INT u;
char name[];
{
   NODE *pnode;
   DOUBLE y[1000], v[1000], miny;
   INT i, n=0, mini;
   FILE *fp;

   fp = fopen(name,"w");
   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      if (fabs(pnode->myvertex->x[0]+pnode->myvertex->x[1]-1.) < 1.e-10){
         y[n] = pnode->myvertex->x[0];
         v[n] = NDS(pnode,u);
         n++;
      }
   do{
      miny = 200.;
      for (i=0; i < n; i++)
         if (y[i] < miny){
            miny = y[i];
            mini = i;
         }
      if (miny < 2.)
         fprintf(fp,"%e %e\n",miny,v[mini]);
      y[mini] = 10.;
   }while (miny < 2.); 
   fclose(fp);
}

void evaluate_parabolic_layers_for_unstructured_grid(mg,u0,u,space,d)
MULTIGRID *mg;
INT u, space;
FLOAT (*u0)(), d;
{
   GRID *tGrid=TOP_GRID(mg);
   ELEMENT *pelem;
   NODE *theNode;
   FLOAT osc_para2, osc_exp, osc_para1, smear;
   FLOAT *x, *x0, *x1, *x2, b[DIM2][DIM2], xc[DIM], grad_u[DIM], s;
   INT d1=1, d2=1;

   osc_para1 = smear = 0.;
   for (theNode=FIRSTN(tGrid); theNode; theNode=SUCC(theNode)){
      x = MYVERTEX(theNode)->x;
      s = NDS(theNode,u) - (*u0)(x);
      if (s < 0.)
         smear += s*s;
      if (x[0] < 0.9 && (x[1] <= 0.1 || x[1] >= 0.9) && s > osc_para1)
         osc_para1 = s;
   }
   smear = sqrt(smear);

   for (pelem = FIRSTELEMENT(tGrid); pelem; pelem = pelem->succ){
      VERTICES_OF_ELEMENT(x0,x1,x2,pelem);
      barycentric_coordinates(x0,x1,x2,b);
      coord_of_barycentre(pelem,xc);
      sgrad_value(pelem,xc,b,u,space,grad_u);
      if (xc[0] < 1.-d){
         if (xc[1] < d){
            if (d1){
               osc_para2 = grad_u[1];
               d1 = 0;
            } 
            else if (osc_para2 > grad_u[1])
               osc_para2 = grad_u[1];
         }
         else if (xc[1] > 1.-d){
            if (d1){
               osc_para2 = -grad_u[1];
               d1 = 0;
            } 
            else if (osc_para2 > -grad_u[1])
               osc_para2 = -grad_u[1];
         }
      }
      else if (xc[1] > d && xc[1] < 1.-d){
         if (d2){
            osc_exp = grad_u[0];
            d2 = 0;
         } 
         else if (osc_exp < grad_u[0])
            osc_exp = grad_u[0];
      }
   }
   osc_para2 *= -1.;
   printf("osc_para1  %9.3e, osc_para2  %9.3e, osc_exp  %9.3e, smear  %9.3e\n",
           osc_para1,osc_para2,osc_exp,smear);
}

#else

void save_p1c_diag_profile(tGrid,u,name)
GRID *tGrid; INT u; char name[];
{ eprintf("Error: save_p1c_diag_profile not available.\n"); }

void save_p1c_diag2_profile(tGrid,u,name)
GRID *tGrid; INT u; char name[];
{ eprintf("Error: save_p1c_diag2_profile not available.\n"); }

void evaluate_parabolic_layers_for_unstructured_grid(mg,u0,u,space,d)
MULTIGRID *mg; INT u, space; FLOAT (*u0)(), d;
{ eprintf("Error: evaluate_parabolic_layers_for_unstructured_grid not available.\n"); }

#endif

#if F_DATA & SCALAR_FACE_DATA

void save_p1nc_profile(tGrid,k,l,d,u,name)
GRID *tGrid;
FLOAT d;
INT k, l, u;
char name[];
{
   ELEMENT *pel;
   FACE *fa0, *fa1, *fa2;
   DOUBLE *x0, *x1, *x2, x01[DIM], x02[DIM], x12[DIM], y[1000], v[1000];
   INT n=0;

   for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ){
      FACES_OF_ELEMENT(fa0,fa1,fa2,pel);
      VERTICES_OF_ELEMENT(x0,x1,x2,pel);
      MIDPOINTS(x0,x1,x2,x01,x02,x12)
      add_point_to_1d_graph(x01[l],x01[k],FD(fa2,u),d,y,v,&n);
      add_point_to_1d_graph(x02[l],x02[k],FD(fa1,u),d,y,v,&n);
      add_point_to_1d_graph(x12[l],x12[k],FD(fa0,u),d,y,v,&n);
   }
   save_1d_graph(y,v,n,name);
}

#else

void save_p1nc_profile(tGrid,k,l,d,u,name)
GRID *tGrid; FLOAT d; INT k, l, u; char name[];
{ eprintf("Error: save_p1nc_profile not available.\n"); }

#endif

void save_profile(tGrid,k,d,u,name,space)
GRID *tGrid;
FLOAT d;
INT k, u, space;
char name[];
{
   INT l=1;

   if (k) l = 0;
   switch(space){
   case P1_NC: save_p1nc_profile(tGrid,k,l,d,u,name);
        break;
   case P1C:
   case Q1C:
   case GP1C:
   case GQ1C:  save_p1c_profile(tGrid,k,l,d,u,name,space);
        break;
   case P2C:
   case Q2C:
   case GP2C:
   case GQ2C:  save_p2c_profile(tGrid,k,l,d,u,name,space);
        break;
   case GP1C_ELBUB: save_p1c_profile(tGrid,k,l,d,u,name,GP1C);
        break;
   case GQ1C_ELBUB: save_p1c_profile(tGrid,k,l,d,u,name,GQ1C);
        break;
   case GQ1X4C: save_gq1x4c_profile(tGrid,k,l,d,u,name,space);
        break;
   case GP2X3C:
   case GP2C_3ELBUB:
   case GP2C_6ELBUB: save_p2c_profile(tGrid,k,l,d,u,name,GP2C);
        break;
   case GQ2C_2ELBUB:
   case GQ2C_3ELBUB: save_p2c_profile(tGrid,k,l,d,u,name,GQ2C);
        break;
   case GQ2X4C: save_gq2x4c_profile(tGrid,k,l,d,u,name,space);
        break;
   default:
        eprintf("Error: save_profile not available.\n");
        break;
   }
}

INT intersections_with_line(pel,a,b,c,z1,z2)
ELEMENT *pel;   /* intersection of pel with the line a*x[0] + b*x[1] + c = 0 */
FLOAT a, b, c, *z1, *z2;
{
   FLOAT *x0, *x1, *x2, c0, c1, c2, min, mid, max, 
         xmin[DIM], xmid[DIM], xmax[DIM], p, q, eps=1.e-15;

   VERTICES_OF_ELEMENT(x0,x1,x2,pel);
   c0 = -a*x0[0] - b*x0[1];   
   c1 = -a*x1[0] - b*x1[1];   
   c2 = -a*x2[0] - b*x2[1];   
   if (c0 < c1){
      min = c0;
      max = c1;
      SET1(xmin,x0)
      SET1(xmax,x1)
   }
   else{
      min = c1;
      max = c0;
      SET1(xmin,x1)
      SET1(xmax,x0)
   }
   if (c2 < min){
      mid = min;
      min = c2;
      SET1(xmid,xmin)
      SET1(xmin,x2)
   } else if (c2 > max){
      mid = max;
      max = c2;
      SET1(xmid,xmax)
      SET1(xmax,x2)
   } else{
      mid = c2;
      SET1(xmid,x2)
   }
   if (min-eps < c && c < max+eps){
      p = (max-c)/(max-min);
      q = 1.-p;
      SET20(z1,xmin,p,xmax,q)
      if (c < mid){
         p = (mid-c)/(mid-min);
         q = 1.-p;
         SET20(z2,xmin,p,xmid,q)
      } else{
         p = (max-c)/(max-mid);
         q = 1.-p;
         SET20(z2,xmid,p,xmax,q)
      }
      return(1);
   }
   else
      return(0);
}

#if N_DATA & SCALAR_NODE_DATA

INT value_at_point_on_edge(n0,n1,x,u,val)  /*  edge through n0 and n1  */
NODE *n0, *n1;
FLOAT *x, *val;
INT u;
{
   FLOAT *x0, *x1, d[DIM], e[DIM], f[DIM], p, eps=1.e-15, eps2=1.e-10;
   INT i=0;

   x0 = n0->myvertex->x;
   x1 = n1->myvertex->x;
   if (fabs(x0[0]-x[0]) < eps && fabs(x0[1]-x[1]) < eps){
      *val = NDS(n0,u);
      i = 1;
   }
   else if (fabs(x1[0]-x[0]) < eps && fabs(x1[1]-x[1]) < eps){
      *val = NDS(n1,u);
      i = 1;
   }
   else {
      SUBTR(x, x0,d);
      SUBTR(x1,x0,e);
      ORT_VECT(f,e)
      p = DOT(e,e);
      if (fabs(DOT(d,f)) < eps2*p){
         p = DOT(d,e)/p;
         if (p > -eps && p < 1.+eps){
            *val = p*NDS(n1,u) + (1.-p)*NDS(n0,u);
            i = 1;
         }
      }
   }
   if (i)
      return(1);
   else
      return(0);
}

#else

INT value_at_point_on_edge(n0,n1,x,u,val)
NODE *n0, *n1; FLOAT *x, *val; INT u;
{ eprintf("Error: value_at_point_on_edge not available.\n"); }

#endif

void value_at_point_on_element_boundary(pel,x,u,val)
ELEMENT *pel;
FLOAT *x, *val;
INT u;
{
   if (!value_at_point_on_edge(pel->n[0],pel->n[1],x,u,val))
      if (!value_at_point_on_edge(pel->n[0],pel->n[2],x,u,val))
         if (!value_at_point_on_edge(pel->n[1],pel->n[2],x,u,val))
            eprintf("Error in value_at_point_on_element_boundary.\n");
}

void sn_sol_graph_on_line(tGrid,a,b,c,u,name)
GRID *tGrid;    /*  graph of u along the line a*x[0] + b*x[1] + c = 0  */
FLOAT a, b, c;
INT u;
char name[];
{
   ELEMENT *pel;
   FLOAT z1[DIM], z2[DIM], val1, val2;
   INT i=1;
   FILE *fp;

   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (intersections_with_line(pel,a,b,c,z1,z2)){
         value_at_point_on_element_boundary(pel,z1,u,&val1);
         value_at_point_on_element_boundary(pel,z2,u,&val2);
         fprintf(fp,"%e %e %e\n",  z1[0],z1[1],val1);
         if (i){
            fprintf(fp,"%e %e %e\n",z2[0],z2[1],val2);
            i = 0;
         }
         fprintf(fp,"%e %e %e\n\n",z2[0],z2[1],val2);
      } 
   fclose(fp);
}

void evaluate_parabolic_bdry_layers_for_p1c(mg,tGrid,u,r,err_sq,u0,
                                                               t,type,structure)
MULTIGRID *mg;
GRID *tGrid;
FLOAT err_sq, (*u0)();
INT u, r, t, type, structure;
{
   INT space=P1C;
   FLOAT err1, err2, err3, err5, err6, err8, errl2m, errl2ma,
         err0p, err1p, err2p, err3p, err0m, err1m, err2m, err3m,
         derx_max0, derx_min0, dery_max0, dery_min0,
         derx_max1, derx_min1, dery_max1, dery_min1,
         derx_max2, derx_min2, dery_max2, dery_min2;
   FILE *fp;

   min_and_max_on_cut_p1c(tGrid,&err1,&err2);
   printf("min = %e, max = %e\n",err1,err2);
   compute_min_max_and_value(mg,u,space,1025,&err1,&err2,&err8);
   printf("min = %e, max = %e, mid = %e\n",err1,err2,err8);
   compute_last_point_err_p1c(mg,u,&err3);
   printf("last point err = %e\n",err3);
   printf("L infinity err = %e\n",err6);
   initialize(mg,r,u0,u0,u0,space,structure);
   subtr(tGrid,u,r,r,t,type);
   err5 = sqrt(dot(tGrid,r,r,t,type));
   printf("l2 error = %e\n",err5);
   printf("L2 error on square = %e\n",err_sq);
   p1c_L_inf_errors(tGrid,u0,u,&err0p,&err1p,&err2p,&err3p,
                    &err0m,&err1m,&err2m,&err3m,&errl2m,&errl2ma);
   check_maximum_principle_for_p1c(tGrid,u,&err0m,&err0p,&err1m,&err1p,
                                           &err2m,&err2p,&err3m,&err3p);
   printf("tmax = %e, rmax = %e, l2tmax = %e, l2rmax = %e\n",
           err0p,err1p,err2p,err3p);
   printf("tmin = %e, rmin = %e, l2tmin = %e, l2rmin = %e\n",
           err0m,err1m,err2m,err3m);
   compute_derivatives(tGrid,u,space,0.1,
                          &derx_max0,&derx_min0,&dery_max0,&dery_min0,
                          &derx_max1,&derx_min1,&dery_max1,&dery_min1,
                          &derx_max2,&derx_min2,&dery_max2,&dery_min2);
   printf("domain 0: u_x in (%8.1e,%8.1e), u_y in (%8.1e,%8.1e)\n",
           derx_min0,derx_max0,dery_min0,dery_max0);
   printf("domain 1: u_x in (%8.1e,%8.1e), u_y in (%8.1e,%8.1e)\n",
           derx_min1,derx_max1,dery_min1,dery_max1);
   printf("domain 2: u_x in (%8.1e,%8.1e), u_y in (%8.1e,%8.1e)\n",
           derx_min2,derx_max2,dery_min2,dery_max2);
   fp = fopen("bex_errors","w");
   fprintf(fp,"                 %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
           err0p,err1p,err2p,err3p,err0m,err1m,err2m,err3m,errl2m,errl2ma);
   fclose(fp);
   fp = fopen("ex_errors","w");
   fprintf(fp,"%8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e\n",
           err2,err1,err8,err3,err6,err5,err_sq);
   fprintf(fp,"                 %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
           err0p,err1p,err2p,err3p,err0m,err1m,err2m,err3m);
   fprintf(fp,"                 %9.4e  %9.4e  %9.4e  %9.4e\n",
           err0m,err1m,err2m,err3m);
   fprintf(fp,"(%8.1e,%9.2e) (%8.1e,%8.1e)  (%8.1e,%8.1e) (%8.1e,%8.1e)  (%8.1e,%8.1e) (%8.1e,%8.1e)\n",
           derx_min0,derx_max0,dery_min0,dery_max0,
           derx_min1,derx_max1,dery_min1,dery_max1,
           derx_min2,derx_max2,dery_min2,dery_max2);
   fprintf(fp,"         & %9.3e &  & %9.3e &  & %9.3e &  & %9.3e & &  \\\\\n",
           err1p,errl2m,dery_min1,derx_max2);
   fclose(fp);
}

void evaluate_parabolic_layers2(mg,u,it,def)
MULTIGRID *mg;
FLOAT def;
INT u, it;
{
   INT s2, s1, scm=SC_TYPE;
   FLOAT err1, err2;
   FILE *fp;

   if (U_SPACE == P1C || U_SPACE == Q1C)
      min_and_max_on_cut_p1c(TOP_GRID(mg),&err1,&err2);
   else if (U_SPACE == GP1C || U_SPACE == GQ1C || 
            U_SPACE == GP1C_ELBUB || U_SPACE == GQ1C_ELBUB)
      min_and_max_on_cut_gp1c(TOP_GRID(mg),&err1,&err2);
   if      (err1 > -1.e-5) s1 = 2;
   else if (err1 > -1.e-3) s1 = 1;
   else if (err1 > -1.e-1) s1 = 0;
   else                    s1 = -2;
   if      (err2 < 1.e-3) s2 = 4;
   else if (err2 < 1.e-2) s2 = 2;
   else if (err2 < 1.e-1) s2 = 0;
   else                   s2 = -4;
   if (COMPUTE_SC==NO) scm = 0;
   printf(
   "SCM=%2i C=%4.2f  osc=%9.3e (%2i) smear=%10.3e (%2i) sc=%2i damp=%4.2f %3iit.  res/rhs=%7.1e\n",
   scm,SC_BETA,err2,s2,err1,s1,s1+s2,DAMPING_PARAM,it,def);
   fp=fopen("ex1_res","a");
   fprintf(fp,
   "SCM=%2i C=%4.2f  osc=%9.3e (%2i) smear=%10.3e (%2i) sc=%2i damp=%4.2f %3iit.  res/rhs=%7.1e\n",
   scm,SC_BETA,err2,s2,err1,s1,s1+s2,DAMPING_PARAM,it,def);
   fclose(fp);
}

void evaluate_skew_convection_for_p1c(mg,tGrid,u)
MULTIGRID *mg;
GRID *tGrid;
INT u;
{
   INT space=P1C;
   FLOAT err1, err2, err3, err5, err6, err8,
         err0p, err1p, err2p, err3p, err0m, err1m;
   FILE *fp;

   min_and_max_for_p1c(tGrid,&err1,&err2);
   printf("min = %e, max-1 = %e, max-min = %e\n",err1,err2-1.,err2-err1);
   evaluate_smearing(mg,u,space,N_FOR_EVAL,&err3,&err8);
   printf("smearing = %e,  %e\n",err3,err8);
/*
   integrals_of_errors(mg,u,space,N_FOR_EVAL,&err5,&err6);
   printf("integrals: %e,  %e\n",err5,err6);
   fp = fopen("ex_errors","w");
   fprintf(fp,"         %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
           err1,err2-1.,err2-err1,err3,err8,err5,err6);
   fclose(fp);
*/
   p1c_smear_and_osc(tGrid,u,&err0p,&err1p,&err2p,&err3p,&err0m,&err1m);
   fp = fopen("bex_errors","w");
   fprintf(fp,"         %9.3e   %9.3e   %9.3e   %9.3e   %9.3e   %9.3e   %9.3e\n",
           err2p,err3p,err0m,err1m,err0p,err1p,err3);
   fclose(fp);
}

void evaluate_skew_convection2(mg,tGrid,u,space,it,def,ref)
MULTIGRID *mg;
GRID *tGrid;
FLOAT def;
INT u, space, it, ref;
{
   INT s1, s2, s3, s4, scm=SC_TYPE;
   FLOAT err3=-1., err0p, err1p, err2p;
   FLOAT smear_2b,smear_2r,in_osc_2,b_osc_2b,b_osc_2r;
   FILE *fp;

   evaluate_smearing(mg,u,space,N_FOR_EVAL,&err3,&err0p);
   printf("smearing = %e\n",err3);
   if (ref == YES){
      p1c_smear_and_osc2_ref(tGrid,u,space,&smear_2b,&smear_2r,&in_osc_2,&b_osc_2b,&b_osc_2r);
      printf(
      "SCM=%2i C=%4.2f iosc=%9.3e eoscb=%9.3e eoscr=%9.3e ismear=%10.3e esmearb=%10.3e esmearr=%10.3e damp=%4.2f %3iit. res/rhs=%7.1e\n",
      scm,SC_BETA,in_osc_2,b_osc_2b,b_osc_2r,err3,smear_2b,smear_2r,DAMPING_PARAM,it,def);
      fp=fopen("ex2_res","a");
      fprintf(fp,
      "SCM=%2i C=%4.2f iosc=%9.3e eoscb=%9.3e eoscr=%9.3e ismear=%10.3e esmearb=%10.3e esmearr=%10.3e damp=%4.2f %3iit. res/rhs=%7.1e\n",
      scm,SC_BETA,in_osc_2,b_osc_2b,b_osc_2r,err3,smear_2b,smear_2r,DAMPING_PARAM,it,def);
      fclose(fp);
   }
   else{
   p1c_smear_and_osc2(tGrid,u,space,&err0p,&err1p,&err2p);
   fp = fopen("bex_errors","w");
   if      (err1p < 1.e-4) s1 = 4;
   else if (err1p < 1.e-2) s1 = 2;
   else if (err1p < 1.e-1) s1 = 0;
   else                    s1 = -4;
   if      (err2p < 1.e-5) s2 = 4;
   else if (err2p < 1.e-4) s2 = 3;
   else if (err2p < 1.e-3) s2 = 2;
   else if (err2p < 1.e-1) s2 = 0;
   else                    s2 = -4;
   if      (err3 < 4.e-2) s3 = 2;
   else if (err3 < 6.e-2) s3 = 1;
   else if (err3 < 8.e-2) s3 = 0;
   else                   s3 = -2;
   if      (err0p < 1.e-4) s4 = 2;
   else if (err0p < 1.e-2) s4 = 1;
   else if (err0p < 5.e-1) s4 = 0;
   else                    s4 = -2;
   if (COMPUTE_SC==NO) scm = 0;
   printf(
   "SCM=%2i C=%4.2f  iosc=%9.3e (%2i) eosc=%9.3e (%2i) ismear=%10.3e (%2i) esmear=%10.3e (%2i)  sc=%2i damp=%4.2f %3iit.  res/rhs=%7.1e\n",
   scm,SC_BETA,err1p,s1,err2p,s2,err3,s3,err0p,s4,s1+s2+s3+s4,DAMPING_PARAM,it,def);
   fp=fopen("ex2_res","a");
   fprintf(fp,
   "SCM=%2i C=%4.2f  iosc=%9.3e (%2i) eosc=%9.3e (%2i) ismear=%10.3e (%2i) esmear=%10.3e (%2i)  sc=%2i damp=%4.2f %3iit.  res/rhs=%7.1e\n",
   scm,SC_BETA,err1p,s1,err2p,s2,err3,s3,err0p,s4,s1+s2+s3+s4,DAMPING_PARAM,it,def);
   fclose(fp);
   }
}

DOUBLE p1c_in_osc(tGrid,u,space)
GRID *tGrid;
INT u, space;
{
   NODE *theNode; 
   FLOAT val, *x, imax=1., imin=0., inosc;

   for (theNode=FIRSTNODE(tGrid); theNode; theNode=SUCC(theNode)){
      x = theNode->myvertex->x;
      val = sn_node_value(theNode,u,space);
      if (x[0] <= 0.5 && x[1] >= 0.25){
         imin = MIN(imin,val);
         imax = MAX(imax,val);
      }
   }
   inosc = MAX(imax-1.,fabs(imin));
   printf("imin = %e, imax = %e, inosc = %e\n",imin,imax,inosc);
   return(inosc);
}

void evaluate_klr(mg,u,space,it,def)
MULTIGRID *mg;
INT u, space, it;
FLOAT def;
{
   INT scm=SC_TYPE;
   FLOAT err1, err2, err4=-1., err5;
   FILE *fp;

   if (space == P1C || space == Q1C || space == GP1C || space == GQ1C)
      evaluate_smearing2(mg,U,U_SPACE,N_FOR_EVAL,&err4,&err1);
   if (space == P1C || space == Q1C)
      sn_min_and_max_for_KLR(TOP_GRID(mg),&err1,&err2);
   else if (space == GP1C || space == GQ1C)
      gmin_and_max_for_KLR(TOP_GRID(mg),&err1,&err2);
   if (COMPUTE_SC==NO) scm = 0;
   printf("%3i  C=%8.2e  min=%9.3e  max=%9.3e  smear=%9.3e  damp=%5.0e %3i it.  res/rhs=%7.1e\n",
   scm,SC_BETA,fabs(err1),err2-1.,err4,DAMPING_PARAM,it,def);
   fp=fopen("KLR_res","a");
   fprintf(fp,"%3i  C=%8.2e  min=%9.3e  max=%9.3e  smear=%9.3e  damp=%5.0e %3i it.  res/rhs=%7.1e\n",
   scm,SC_BETA,fabs(err1),err2-1.,err4,DAMPING_PARAM,it,def);
   fclose(fp);
}

char *sc_name(sc_type)
INT sc_type;
{
     switch(sc_type) {
      case 0: if (SDFEM==NO) return("GAL ");
              else return("SUPG ");
      break;
      case SCM_C: return("C93  ");
      break;
      case SCM_MH: return("MH85 ");
      break;
      case SCM_GC2: return("dCG91");
      break;
      case SCM_K4: return("BE02 ");
      break;
      case SCM_BH: return("BH04 ");
      break;
      case SCM_BE2: return("BE51 ");
      break;
      case SCM_BE3: return("BE05 ");
      break;
      default:
         eprintf("Error: name not defined for sc_type used.\n");
   }
}

void evaluate_klr2(mg,u,space,it,def)
MULTIGRID *mg;
INT u, space, it;
FLOAT def;
{
   INT scm=SC_TYPE;
   FLOAT err1, err2, err4=-1., err5, 
         min, max, min2, max2, min_out, max_out, mid_val, smear, width;
   FILE *fp;

   if (space == P1C || space == Q1C || space == GP1C || space == GQ1C){
      evaluate_smearing1(mg,U,U_SPACE,N_FOR_EVAL,&smear,&width,&mid_val);
      gmin_and_max_for_KLR2(TOP_GRID(mg),&min,&max,&min2,&max2,&min_out,
                            &max_out,&mid_val,U_SPACE);
      compute_maxo_for_klr(TOP_GRID(mg),U,U_SPACE,&max_out);
   }
   min = fabs(min);
   max = max-1.;
   min_out = fabs(min_out);
//   max_out = max_out-mid_val;
   printf("Grid - absolut: smear = %f, width = %f\n",smear,width);
   printf("Grid - relativ: smear/0.0898437 - 1. = %f\n",smear/0.0898437 - 1.);
   printf("Grid - relativ: width/0.3788512 - 1. = %f\n",width/0.3788512 - 1.);
   width = width/0.385697 - 1.;
   smear = smear/0.114518 - 1.;
   if (COMPUTE_SC==NO) scm = 0;
   printf("%3i  C=%8.2e  min=%9.3e  max=%9.3e  min2=%9.3e  max2=%9.3e mino=%9.3e  maxo=%9.3e  smear=%9.3e  width=%9.3e  midv=%9.3e  damp=%5.0e %3i it.  res/rhs=%7.1e\n",
   scm,SC_BETA,min,max,min2,max2,min_out,max_out,smear,width,1.-mid_val,DAMPING_PARAM,it,def);
   fp=fopen("KLR_res","a");
/*
   fprintf(fp,"%3i  C=%8.2e  min=%9.3e  max=%9.3e  min2=%9.3e  max2=%9.3e mino=%9.3e  maxo=%9.3e  smear=%9.3e  width=%9.3e  midv=%9.3e  damp=%5.0e %3i it.  res/rhs=%7.1e\n",
   scm,SC_BETA,min,max,min2,max2,min_out,max_out,smear,width,1.-mid_val,DAMPING_PARAM,it,def);
*/
printf("ha: %s\n",sc_name(scm));
   fprintf(fp,"%5s min=%9.3e max=%9.3e min2=%9.3e max2=%9.3e mino=%9.3e maxo=%9.3e smear=%9.3e width=%9.3e midv=%9.3e res/f=%7.1e\n",
   sc_name(scm),min,max,min2,max2,min_out,max_out,smear,width,1.-mid_val,def);
   fclose(fp);
   fp=fopen("KLR_tab","a");
   fprintf(fp,"%5s & %9.3e & %9.3e & %9.3e & %9.3e & %9.3e & %9.3e & %9.3e & %9.3e \\\\\n",
   sc_name(scm),min,max,min2,max2,min_out,max_out,smear,width);
   fclose(fp);
}

void evaluate_hill(mg,u,space,it,def)
MULTIGRID *mg;
INT u, space, it;
FLOAT def;
{
   FLOAT err1, err2, err4, err5;
   FILE *fp;

   min_and_max_for_gp1c_and_hill(TOP_GRID(mg),&err1,&err2);
   printf("min=%9.3e  diff=%9.3e\n",err1,err2);
}

void evaluate_hemker(mg,u,space,it,def)
MULTIGRID *mg;
INT u, space, it;
FLOAT def;
{
   FLOAT err1, err2, err4, err5;
   FILE *fp;

   if (space == P1C || space == Q1C)
      min_and_max_for_p1c_and_hemker(TOP_GRID(mg),&err1,&err2,&err4,&err5);
   else if (space == GP1C || space == GQ1C)
      min_and_max_for_gp1c_and_hemker(TOP_GRID(mg),&err1,&err2,&err4,&err5);
   printf(
   "%3i  C=%8.2e  min1=%9.3e  max1=%9.3e  min2=%9.3e  max2=%9.3e  damp=%5.0e %3i it.  def=%7.1e\n",
   SC_TYPE,SC_BETA,fabs(err1),err2-1.,fabs(err4),err5-1.,DAMPING_PARAM,it,def);
   fp=fopen("hemker_res","a");
   fprintf(fp,
   "%3i  C=%8.2e  min1=%9.3e  max1=%9.3e  min2=%9.3e  max2=%9.3e  damp=%5.0e %3i it.  def=%7.1e\n",
   SC_TYPE,SC_BETA,fabs(err1),err2-1.,fabs(err4),err5-1.,DAMPING_PARAM,it,def);
   fclose(fp);
}

void g_conv_diff_ref_values(tGrid,u0,u01,u02,u03,u,d,r,
                                             u_space,u_structure,t_for_u,u_type)
GRID *tGrid;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*u03)();
INT u, d, r, u_space, u_structure, t_for_u, u_type;
{
   FLOAT err1, err2, x1[2], x2[2], x3[2], val1, val2, val3, 
         val1n, val2n, val3n, val1f, val2f, val3f, val1e, val2e, val3e;

   sL2_error(tGrid,u,r,&err1,u0,u_space,u_structure,t_for_u,u_type);
   sH10_error(tGrid,u,r,&err2,u01,u02,u03,u_space);
   if (ELEMENT_TYPE == SIMPLEX){
      x1[0] = x2[0] = x3[0] = 0.5 + 1./128.;
      x1[1] = 0.25 - 1./128.;
      x2[1] = 0.5  - 1./128.;
      x3[1] = 0.75 - 1./128.;
   }
   else if (ELEMENT_TYPE == CUBE){
      x1[0] = x2[0] = x3[0] = 0.5 + 1./64.;
      x1[1] = 0.25 + 1./64.;
      x2[1] = 0.5  + 1./64.;
      x3[1] = 0.75 + 1./64.;
   }
   val1 = gsvalue_at_point(tGrid,x1,u,REF_MAP,ELEM);     
   val2 = gsvalue_at_point(tGrid,x2,u,REF_MAP,ELEM);     
   val3 = gsvalue_at_point(tGrid,x3,u,REF_MAP,ELEM);     
/*
   val1 = svalue_at_point(tGrid,x1,u,u_space);
   val2 = svalue_at_point(tGrid,x2,u,u_space);
   val3 = svalue_at_point(tGrid,x3,u,u_space);
*/
   printf("L2: %12.6e\n",err1);
   printf("H1-semi: %12.6e\n",err2);
   printf("error in (%12.6e,%12.6e): %12.6e\n",x1[0],x1[1],fabs(val1-u0(x1)));
   printf("error in (%12.6e,%12.6e): %12.6e\n",x2[0],x2[1],fabs(val2-u0(x2)));
   printf("error in (%12.6e,%12.6e): %12.6e\n",x3[0],x3[1],fabs(val3-u0(x3)));
/*
   copy(tGrid,u,d,0,u_type);
   mv_set_value_f(tGrid,0.,d,0);
   mv_set_value_e(tGrid,0.,d);
   val1n = gsvalue_at_point(tGrid,x1,d,REF_MAP,ELEM);     
   val2n = gsvalue_at_point(tGrid,x2,d,REF_MAP,ELEM);     
   val3n = gsvalue_at_point(tGrid,x3,d,REF_MAP,ELEM);     
   copy(tGrid,u,d,0,u_type);
   mv_set_value_n(tGrid,0.,d,0);
   mv_set_value_e(tGrid,0.,d);
   val1f = gsvalue_at_point(tGrid,x1,d,REF_MAP,ELEM);     
   val2f = gsvalue_at_point(tGrid,x2,d,REF_MAP,ELEM);     
   val3f = gsvalue_at_point(tGrid,x3,d,REF_MAP,ELEM);     
   copy(tGrid,u,d,0,u_type);
   mv_set_value_n(tGrid,0.,d,0);
   mv_set_value_f(tGrid,0.,d,0);
   val1e = gsvalue_at_point(tGrid,x1,d,REF_MAP,ELEM);     
   val2e = gsvalue_at_point(tGrid,x2,d,REF_MAP,ELEM);     
   val3e = gsvalue_at_point(tGrid,x3,d,REF_MAP,ELEM);     
   printf("1: n=%12.6e f=%12.6e e=%12.6e sum-val=%12.6e\n",
                                val1n,val1f,val1e,fabs(val1n+val1f+val1e-val1));
   printf("2: n=%12.6e f=%12.6e e=%12.6e sum-val=%12.6e\n",
                                val2n,val2f,val2e,fabs(val2n+val2f+val2e-val2));
   printf("3: n=%12.6e f=%12.6e e=%12.6e sum-val=%12.6e\n",
                                val3n,val3f,val3e,fabs(val3n+val3f+val3e-val3));

   printf("Without bubble:\n");
   copy(tGrid,u,d,0,u_type);
   mv_set_value_e(tGrid,0.,d);
   sL2_error(tGrid,d,r,&err1,u0,u_space,u_structure,t_for_u,u_type);
   sH10_error(tGrid,d,r,&err2,u01,u02,u03,u_space);
   val1 = gsvalue_at_point(tGrid,x1,d,REF_MAP,ELEM);     
   val2 = gsvalue_at_point(tGrid,x2,d,REF_MAP,ELEM);     
   val3 = gsvalue_at_point(tGrid,x3,d,REF_MAP,ELEM);     
   printf("L2: %12.6e\n",err1);
   printf("H1-semi: %12.6e\n",err2);
   printf("error in (%12.6e,%12.6e): %12.6e\n",x1[0],x1[1],fabs(val1-u0(x1)));
   printf("error in (%12.6e,%12.6e): %12.6e\n",x2[0],x2[1],fabs(val2-u0(x2)));
   printf("error in (%12.6e,%12.6e): %12.6e\n",x3[0],x3[1],fabs(val3-u0(x3)));
*/
}

DOUBLE compute_consistency_error_for_conv_lp(mg,u0,u01,u02,rhs)
MULTIGRID *mg;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*rhs)();
{
   GRID *topGrid=TOP_GRID(mg);
   FLOAT err=1., err1, err2, err3, err4, err44, err55, err6, err7;
   INT i=1;

   TNU = VALUE_OF_EPS;
   boundary_values(mg,u0,u0,u0,U,D,Q,U_SPACE,U_STRUCTURE);
   set_value(topGrid,0.,U,T_FOR_U,U_TYPE);

   set_mat_value(topGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
   if (LOC_PROJECT == YES)
      add_lp_term_matr(topGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (CONV_LOC_PROJECT == YES)
      add_conv_lp_term_matr(topGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
/*
   set_value(topGrid,0.,UU,0,U_TYPE);
   boundary_values(mg,u0,u0,u0,UU,D,Q,U_SPACE,U_STRUCTURE);
   mult_A(topGrid,A,UU,R,Q,STOP_IS_FIRST_INNER,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   initialize(mg,UU,u0,u0,u0,GQ2C,U_STRUCTURE);
   mult_A(topGrid,A,UU,F,Q,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
          0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   add(topGrid,R,F,F,T_FOR_U,U_TYPE);
*/
   set_value(topGrid,0.,F,0,U_TYPE);
   general_consistency_err_for_conv_lp(topGrid,F,u01,u02,TNU,3,bb0,bb1,
                           one_fcn,r_l1,r_l2,NULL,NULL,NULL,NULL,NULL,
                           g_consist_err_conv_local_proj_ref,
                           ELEM,REF_MAP,REACT_Q_RULE);

   if (METHOD == UMFPACK){
      fill_Ax(topGrid,A,Ap,Ai,Ax,U_TYPE,U_SPACE);
      make_vector_from_grid_data(topGrid,F,Rhs,U_TYPE,U_SPACE);
      solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
      make_grid_data_from_vector(topGrid,U,Sol,U_TYPE,U_SPACE);
   }
   defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
   printf("LP:    def/rhs= %e,  |z_h| = %e\n",
          sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/
               dot(topGrid,F,F,T_FOR_U,U_TYPE)),
          sqrt(dot(topGrid,F,U,T_FOR_U,U_TYPE)));

//-----------------------------------------------------
   set_mat_value(topGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
   add_react_term_matr(topGrid,A,react,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (METHOD == UMFPACK){
      fill_Ax(topGrid,A,Ap,Ai,Ax,U_TYPE,U_SPACE);
      make_vector_from_grid_data(topGrid,F,Rhs,U_TYPE,U_SPACE);
      solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
      make_grid_data_from_vector(topGrid,U,Sol,U_TYPE,U_SPACE);
   }
   defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
   printf("L2:    def/rhs= %e,  |z_h| = %e\n",
          sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/
               dot(topGrid,F,F,T_FOR_U,U_TYPE)),
          sqrt(dot(topGrid,F,U,T_FOR_U,U_TYPE)));
//-----------------------------------------------------
   set_mat_value(topGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
   add_Laplace_matr(topGrid,TNU,A,U_SPACE,A_STRUCT,U_STRUCTURE,KORN_LAPLACE);
   add_react_term_matr(topGrid,A,react,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (LOC_PROJECT == YES)
      add_lp_term_matr(topGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (CONV_LOC_PROJECT == YES)
      add_conv_lp_term_matr(topGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (METHOD == UMFPACK){
      fill_Ax(topGrid,A,Ap,Ai,Ax,U_TYPE,U_SPACE);
      make_vector_from_grid_data(topGrid,F,Rhs,U_TYPE,U_SPACE);
      solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
      make_grid_data_from_vector(topGrid,U,Sol,U_TYPE,U_SPACE);
   }
   defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
   printf("all:   def/rhs= %e,  |z_h| = %e\n",
          sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/
               dot(topGrid,F,F,T_FOR_U,U_TYPE)),
          sqrt(dot(topGrid,F,U,T_FOR_U,U_TYPE)));

//    mult_A(topGrid,A,U,R,Q,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
//        0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
//    printf("|z_h| = %e\n",sqrt(dot(topGrid,R,U,T_FOR_U,U_TYPE)));
//    sL2_error(topGrid,U,R,&err1,p0,U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
//    sH10_error(topGrid,U,R,&err2,p0,p0,p0,U_SPACE);
/*
   set_value(topGrid,0.,UU,0,U_TYPE);
   initialize(mg,UU,u0,u0,u0,GQ2C,U_STRUCTURE);
   if (LOC_PROJECT == YES)
      s_lp_error(topGrid,UU,&err3,TNU,bb0,bb1,bb1,p0,p0,p0,U_SPACE);
   if (CONV_LOC_PROJECT == YES)
      s_conv_lp_error(topGrid,UU,&err3,TNU,bb0,bb1,bb1,p0,p0,p0,U_SPACE);
*/
   set_value(topGrid,0.,D,0,U_TYPE);
   if (LOC_PROJECT == YES)
      s_lp_error(topGrid,D,&err3,TNU,bb0,bb1,bb1,u01,u02,p0,U_SPACE);
   if (CONV_LOC_PROJECT == YES)
      s_conv_lp_error(topGrid,D,&err3,TNU,bb0,bb1,bb1,u01,u02,p0,U_SPACE);
}

void evaluate_ex_55(tGrid,y,d,tau0,u,ie,name_graph,name_err,name_sum)
GRID *tGrid;
FLOAT y, d, tau0;
INT u, ie;
char name_graph[], name_err[], name_sum[];
{
   DOUBLE v0, v1, z[2]={0.,y}, s, sum=0., err=0.;
   INT i=1, n=0;
   FILE *fp;

   fp = fopen(name_graph,"w");
   v0 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
   fprintf(fp,"%e %e\n",z[0],v0);
   while ((z[0]+=d) < 0.9999){
      v1 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
      fprintf(fp,"%e %e\n",z[0],v1);
      s = z[0] - v1;
      if (ie){
         if (i){
            err += s*s;
            i = 0;
         }
         else
            i = 1; 
      }
      else
         err += s*s;
      sum += fabs(v0-v1);
      v0 = v1;
   }
   fprintf(fp,"%e %e\n",1.,0.);
   fclose(fp);
   err = sqrt(err);
   printf("err = %12.6e, sum = %12.6e\n",err,sum);
   fp = fopen(name_err,"a");
   fprintf(fp,"%e %e\n",tau0,err);
   fclose(fp);
   fp = fopen(name_sum,"a");
   fprintf(fp,"%e %e\n",tau0,sum);
   fclose(fp);
}

void evaluate_ex_55_half(tGrid,y,d,tau0,u,u0,name_graph,name_err,name_sum)
GRID *tGrid;
FLOAT y, d, tau0, (*u0)();
INT u;
char name_graph[], name_err[], name_sum[];
{
   DOUBLE v0, v1, z[2]={0.,y}, s, sum=0., err=0., max;
   FILE *fp;

   fp = fopen(name_graph,"w");
   v0 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
   fprintf(fp,"%e %e\n",z[0],v0);
   max = fabs(v0);
   while ((z[0]+=d) < 0.5000001){
      v1 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
      fprintf(fp,"%e %e\n",z[0],v1);
      max = MAX(max,fabs(v1));
      s = u0(z) - v1;
      err += s*s;
      sum += fabs(v0-v1);
      v0 = v1;
   }
   fclose(fp);
   err = sqrt(err);
   sum /= max;
   printf("err = %12.6e, sum = %12.6e\n",err,sum);
   fp = fopen(name_err,"a");
   fprintf(fp,"%e %e\n",tau0,err);
   fclose(fp);
   fp = fopen(name_sum,"a");
   fprintf(fp,"%e %e\n",tau0,sum);
   fclose(fp);
}

void evaluate_ex_for_gert(tGrid,x,d,tau0,u0,u,name_graph,name_err)
GRID *tGrid;
FLOAT x, d, tau0, (*u0)();
INT u;
char name_graph[], name_err[];
{
   DOUBLE v0, v1, z[2]={x,0.}, err, sum1=0., sum2=0., max=0., tv;
   FILE *fp;

   fp = fopen(name_graph,"w");
   v0 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
   tv = -v0;
   fprintf(fp,"%e %e\n",z[1],v0);
   while ((z[1]+=d) < 0.9999){
      v1 = gsvalue_at_point(tGrid,z,u,REF_MAP,ELEM);
      fprintf(fp,"%e %e\n",z[1],v1);
      err = fabs(u0(z) - v1);
      sum1 += err;
      sum2 += err*err;
      max = MAX(max,err);
      tv += fabs(v0-v1);
      v0 = v1;
   }
   fprintf(fp,"%e %e\n",1.,0.);
   fclose(fp);
   sum2 = sqrt(sum2);
   printf("l1 = %12.6e, l2 = %12.6e, linf = %12.6e, tv = %12.6e\n",
          sum1,sum2,max,tv);
   fp = fopen(name_err,"a");
   fprintf(fp,"%e %e %e %e %e\n",tau0,sum1,sum2,max,tv);
   fclose(fp);
}

void compute_errors_for_gert(tGrid,u,r,tau0,eps,u0,u01,u02,u03,
                        t_for_u,u_type,u_space,u_structure,finite_el,name,name1)
GRID *tGrid;
FLOAT tau0, eps, (*u0)(), (*u01)(), (*u02)(), (*u03)();
INT u, r, t_for_u, u_type, u_space, u_structure;
FINITE_ELEMENT finite_el;
char name[], name1[];
{
   FLOAT err1, err2, err3, err4, err5, err6,
         err11, err22, err33, err44, err55, err66, err77;
   FILE *fp;

/* ------------------------------------------------------------------------ */
   printf("\ncomplete solution, with belems:\n");
   general_sL2test(tGrid,u0,u,L2_Q_RULE,REF_MAP,finite_el,&err1,general_sL2);
   general_sH10_test(tGrid,u01,u02,u03,u,
                     H10_Q_RULE,REF_MAP,finite_el,&err2,general_sH10);
   if (LOC_PROJECT == YES)
      s_lp_error(tGrid,u,&err3,eps,bb0,bb1,bb1,u01,u02,u03,u_space);
   if (CONV_LOC_PROJECT == YES)
      s_conv_lp_error(tGrid,u,&err3,eps,bb0,bb1,bb1,u01,u02,u03,u_space);
   if (LOC_PROJECT == YES || CONV_LOC_PROJECT == YES)
      printf("lp err= %e\n",err3=sqrt(err3*err3+eps*err2*err2));
   else
      err3=0.;
/* ------------------------------------------------------------------------ */
   printf("\ncomplete solution, without belems:\n");
   general_sL2test_without_belems(tGrid,u0,u,L2_Q_RULE,REF_MAP,finite_el,
                                  &err11,general_sL2);
   general_sH10_test_without_belems(tGrid,u01,u02,u03,u,
                              H10_Q_RULE,REF_MAP,finite_el,&err22,general_sH10);
   if (LOC_PROJECT == YES)
      s_lp_error_without_belems(tGrid,u,&err33,eps,
                                bb0,bb1,bb1,u01,u02,u03,u_space);
   if (CONV_LOC_PROJECT == YES)
      s_conv_lp_error_without_belems(tGrid,u,&err33,eps,
                                     bb0,bb1,bb1,u01,u02,u03,u_space);
   if (LOC_PROJECT == YES || CONV_LOC_PROJECT == YES)
      printf("lp err= %e\n",err33=sqrt(err33*err33+eps*err22*err22));
   else
      err33=0.;
/*========================================================================= */
   if (u_space == GP1C_ELBUB  || u_space == GQ1C_ELBUB ||
       u_space == GP2C_3ELBUB || u_space == GP2C_6ELBUB || 
       u_space == GQ2C_2ELBUB || u_space == GQ2C_3ELBUB){
      if (u_space == GQ2C_2ELBUB || u_space == GQ2C_3ELBUB || 
          u_space == GP2C_3ELBUB || u_space == GP2C_6ELBUB)
         mv_set_value_es(tGrid,0.,u);
      else
         mv_set_value_e(tGrid,0.,u);
/* ------------------------------------------------------------------------ */
      printf("\nsolution without bubbles, with belems:\n");
      general_sL2test(tGrid,u0,u,L2_Q_RULE,REF_MAP,finite_el,&err4,general_sL2);
      general_sH10_test(tGrid,u01,u02,u03,u,
                        H10_Q_RULE,REF_MAP,finite_el,&err5,general_sH10);
      if (LOC_PROJECT == YES)
         s_lp_error(tGrid,u,&err6,eps,bb0,bb1,bb1,u01,u02,u03,u_space);
      if (CONV_LOC_PROJECT == YES)
         s_conv_lp_error(tGrid,u,&err6,eps,bb0,bb1,bb1,u01,u02,u03,u_space);
      printf("lp err= %e\n",err6=sqrt(err6*err6+eps*err5*err5));
/* ------------------------------------------------------------------------ */
      printf("\nsolution without bubbles, without belems:\n");
      general_sL2test_without_belems(tGrid,u0,u,L2_Q_RULE,REF_MAP,finite_el,
                                     &err44,general_sL2);
      general_sH10_test_without_belems(tGrid,u01,u02,u03,u,H10_Q_RULE,REF_MAP,
                                       finite_el,&err55,general_sH10);
      if (LOC_PROJECT == YES)
         s_lp_error_without_belems(tGrid,u,&err66,eps,
                                   bb0,bb1,bb1,u01,u02,u03,u_space);
      if (CONV_LOC_PROJECT == YES)
         s_conv_lp_error_without_belems(tGrid,u,&err66,eps,
                                        bb0,bb1,bb1,u01,u02,u03,u_space);
      printf("lp err= %e\n",err66=sqrt(err66*err66+eps*err55*err55));
/* ------------------------------------------------------------------------ */
   }
   sL_infinity_error_without_belems(tGrid,u,r,&err77,u0,u_space);
   if (u_space == GP1C_ELBUB  || u_space == GQ1C_ELBUB ||
       u_space == GP2C_3ELBUB || u_space == GP2C_6ELBUB || 
       u_space == GQ2C_2ELBUB || u_space == GQ2C_3ELBUB){
      fp = fopen(name,"a");
      fprintf(fp,"%e %e %e %e %e %e %e\n",tau0,err1,err2,err3,err4,err5,err6);
      fclose(fp);
      fp = fopen(name1,"a");
      fprintf(fp,"%e %e %e %e %e %e %e %e\n",
                                tau0,err11,err22,err33,err44,err55,err66,err77);
      fclose(fp);
   }
   else{
      fp = fopen(name,"a");
      fprintf(fp,"%e %e %e %e %e %e %e %e\n",
                                   tau0,err1,err2,err3,err11,err22,err33,err77);
      fclose(fp);
   }
}

#if N_DATA & MVECTOR_NODE_DATA

void check_symmetry(tGrid,u)
GRID *tGrid;
INT u;
{
   NODE *pn1, *pn2;
   FLOAT sum=0.; 

   for (pn1 = FIRSTNODE(tGrid); pn1; pn1 = pn1->succ)
   if (fabs(pn1->myvertex->x[1]-0.5) > 1.e-6){
      for (pn2 = FIRSTNODE(tGrid); pn2 == pn1 || 
         fabs(pn2->myvertex->x[0] - pn1->myvertex->x[0]) > 1.e-6 ||
         fabs(pn2->myvertex->x[1] + pn1->myvertex->x[1] -1.) > 1.e-6;
                                                        pn2 = pn2->succ);
/*
      printf("(%f,%f) (%f,%f) (%e,%e) %e %e %e\n",
          pn1->myvertex->x[0],pn1->myvertex->x[1],pn2->myvertex->x[0],
          pn2->myvertex->x[1],fabs(pn2->myvertex->x[0] - pn1->myvertex->x[0]),
          fabs(pn2->myvertex->x[1] + pn1->myvertex->x[1] -1.),
          fabs(NDMV(pn1,u,0)-NDMV(pn2,u,0)),NDMV(pn1,u,0),NDMV(pn2,u,0));
*/
      sum += fabs(NDMV(pn1,u,0)-NDMV(pn2,u,0));
   }
   printf("measure of symmetry: %e\n",sum);
}

#else

void check_symmetry(tGrid,u)
GRID *tGrid; INT u;
{  eprintf("Error: check_symmetry not available.\n");  }

#endif

#if (N_DATA & NxN_NODE_MATR) && (N_DATA & NxM_NODE_FACE_MATR) && (F_DATA & NxN_FACE_MATR) && (F_DATA & MxN_FACE_NODE_MATR) && (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES) && (E_DATA & MxM_E_E_MATR) && (E_DATA & MxN_E_N_MATR) && (E_DATA & NxM_N_E_MATR) && (E_DATA & MxN_E_F_MATR) && (E_DATA & NxM_F_E_MATR) && (N_DATA & MVECTOR_NODE_DATA) && (F_DATA & MVECTOR_FACE_DATA) && (E_DATA & MVECTOR_ELEMENT_DATA) && (DIM == 2)

void store_q1_on_fine_el(Z,n41,n11,n12,ni,f4,n1,f1,pel,i)
ELEMENT *pel;
NODE *n41, *n11, *n12, *ni, *n1;
FACE *f4, *f1;
INT Z, i;
{
   LINK *pli;
   NFLINK *pnf;
   FLINK *pfl;
   FNLINK *pfn;

   COEFF_NN(n1,Z,0,0) = COEFF_NN(n11,Z,0,0);
   COEFF_NN(f1,Z,0,0) = COEFF_NN(n12,Z,0,0);
   for (pli = n11->tstart; pli->nbnode != n12; pli = pli->next);
   for (pnf = n1->tnfstart; pnf->nbface != f1; pnf = pnf->next);
   COEFF_NN(pnf,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n11->tstart; pli->nbnode != n41; pli = pli->next);
   for (pnf = n1->tnfstart; pnf->nbface != f4; pnf = pnf->next);
   COEFF_NN(pnf,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n41->tstart; pli->nbnode != n11; pli = pli->next);
   for (pfn = f4->tfnstart; pfn->nbnode != n1; pfn = pfn->next);
   COEFF_NN(pfn,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n41->tstart; pli->nbnode != n12; pli = pli->next);
   for (pfl = f4->tfstart; pfl->nbface != f1; pfl = pfl->next);
   COEFF_NN(pfl,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n12->tstart; pli->nbnode != n11; pli = pli->next);
   for (pfn = f1->tfnstart; pfn->nbnode != n1; pfn = pfn->next);
   COEFF_NN(pfn,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n12->tstart; pli->nbnode != n41; pli = pli->next);
   for (pfl = f1->tfstart; pfl->nbface != f4; pfl = pfl->next);
   COEFF_NN(pfl,Z,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n11->tstart; pli->nbnode != ni; pli = pli->next);
   COEFF_NE_NM(pel,Z,i,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = ni->tstart; pli->nbnode != n11; pli = pli->next);
   COEFF_EN_MN(pel,Z,i,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = n12->tstart; pli->nbnode != ni; pli = pli->next);
   COEFF_FE_NM(pel,Z,i,0,0) = COEFF_NN(pli,Z,0,0);
   for (pli = ni->tstart; pli->nbnode != n12; pli = pli->next);
   COEFF_EF_MN(pel,Z,i,0,0) = COEFF_NN(pli,Z,0,0);
}

void store_q1_on_four_el(pel,Z)
ELEMENT *pel;
INT Z;
{
   NODE *n1, *n2, *n3, *n4, *n11, *n22, *n33, *n44, *n12, *n23, *n34, *n41, *ni;
   FACE *f1, *f2, *f3, *f4;
 
   NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
   FACES_OF_4ELEMENT(f1,f2,f3,f4,pel);
   NODES_OF_4ELEMENT(n41,n11,n12,ni,pel->sons[0]);
   NODES_OF_4ELEMENT(n23,n33,n34,ni,pel->sons[2]);
   n22 = n2->son;
   n44 = n4->son;
   store_q1_on_fine_el(Z,n41,n11,n12,ni,f4,n1,f1,pel,0);
   store_q1_on_fine_el(Z,n12,n22,n23,ni,f1,n2,f2,pel,1);
   store_q1_on_fine_el(Z,n23,n33,n34,ni,f2,n3,f3,pel,2);
   store_q1_on_fine_el(Z,n34,n44,n41,ni,f3,n4,f4,pel,3);
   COEFF_EE_MM(pel,Z,0,0) = COEFF_NN(ni,Z,0,0);
}

void store_q1_matrix_on_coarser_grid(tGrid,Z) /* tGrid is the coarse grid */
GRID *tGrid;
INT Z;
{
   ELEMENT *pel;

   set_mat_value(tGrid,Z,0.,0,0,U_TYPE,U_TYPE,A_STRUCT);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      store_q1_on_four_el(pel,Z);
}

void store_q1_vector_on_coarser_grid(tGrid,u) /* tGrid is the coarse grid */
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   NODE *pnode;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDMV(pnode,u,0) = NDMV(pnode->son,u,0); 
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      FDMV(pel->f[0],u,0) = NDMV(pel->sons[0]->n[2],u,0);
      FDMV(pel->f[1],u,0) = NDMV(pel->sons[1]->n[2],u,0);
      FDMV(pel->f[2],u,0) = NDMV(pel->sons[2]->n[2],u,0);
      FDMV(pel->f[3],u,0) = NDMV(pel->sons[3]->n[2],u,0);
      EDMV(pel,u,0) = NDMV(pel->sons[0]->n[3],u,0);
   }
}

void store_q1_vector_on_finer_grid(tGrid,u) /* tGrid is the coarse grid */
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   NODE *pnode;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDMV(pnode->son,u,0) = NDMV(pnode,u,0);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      NDMV(pel->sons[0]->n[2],u,0) = FDMV(pel->f[0],u,0);
      NDMV(pel->sons[1]->n[2],u,0) = FDMV(pel->f[1],u,0);
      NDMV(pel->sons[2]->n[2],u,0) = FDMV(pel->f[2],u,0);
      NDMV(pel->sons[3]->n[2],u,0) = FDMV(pel->f[3],u,0);
      NDMV(pel->sons[0]->n[3],u,0) = EDMV(pel,u,0);
   }
}

void store_q2_on_fine_el(Z,n41,n11,n12,ni,f411,f112,f12i,f41i,
                         f4,n1,f1,pel,i,i12,i14)
ELEMENT *pel;
NODE *n41, *n11, *n12, *ni, *n1;
FACE *f411, *f112, *f12i, *f41i, *f4, *f1;
INT Z, i, i12, i14;
{
   NFLINK *pnff, *pnfc;
   FLINK *pflf, *pflc;
   FNLINK *pfnf, *pfnc;
   INT i1=i+1, i4=i, i5=i+5, im=i-1;

   if (i == 0){
      i4 = 4;
      im = 3;
   }

   store_q1_on_fine_el(Z,n41,n11,n12,ni,f4,n1,f1,pel,i);

   for (pnff = n11->tnfstart; pnff->nbface != f112; pnff = pnff->next);   
   for (pnfc = n1->tnfstart; pnfc->nbface != f1; pnfc = pnfc->next);
   COEFF_NN(pnfc,Z,0,i12) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n11->tnfstart; pnff->nbface != f411; pnff = pnff->next);   
   for (pnfc = n1->tnfstart; pnfc->nbface != f4; pnfc = pnfc->next);
   COEFF_NN(pnfc,Z,0,i14) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n11->tnfstart; pnff->nbface != f12i; pnff = pnff->next);   
   COEFF_NE_NM(pel,Z,i,0,i1) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n11->tnfstart; pnff->nbface != f41i; pnff = pnff->next);   
   COEFF_NE_NM(pel,Z,i,0,i4) = COEFF_NN(pnff,Z,0,0);
   COEFF_NE_NM(pel,Z,i,0,i5) = COEFF_NE_NM(pel->sons[i],Z,1,0,0);

   for (pnff = n12->tnfstart; pnff->nbface != f112; pnff = pnff->next);   
   COEFF_NN(f1,Z,0,i12) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n12->tnfstart; pnff->nbface != f411; pnff = pnff->next);   
   for (pflc = f1->tfstart; pflc->nbface != f4; pflc = pflc->next);
   COEFF_NN(pflc,Z,0,i14) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n12->tnfstart; pnff->nbface != f12i; pnff = pnff->next);   
   COEFF_FE_NM(pel,Z,i,0,i1) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n12->tnfstart; pnff->nbface != f41i; pnff = pnff->next);   
   COEFF_FE_NM(pel,Z,i,0,i4) = COEFF_NN(pnff,Z,0,0);
   COEFF_FE_NM(pel,Z,i,0,i5) = COEFF_NE_NM(pel->sons[i],Z,2,0,0);

   for (pnff = n41->tnfstart; pnff->nbface != f411; pnff = pnff->next);   
   COEFF_NN(f4,Z,0,i14) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n41->tnfstart; pnff->nbface != f112; pnff = pnff->next);   
   for (pflc = f4->tfstart; pflc->nbface != f1; pflc = pflc->next);
   COEFF_NN(pflc,Z,0,i12) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n41->tnfstart; pnff->nbface != f12i; pnff = pnff->next);   
   COEFF_FE_NM(pel,Z,im,0,i1) = COEFF_NN(pnff,Z,0,0);
   for (pnff = n41->tnfstart; pnff->nbface != f41i; pnff = pnff->next);   
   COEFF_FE_NM(pel,Z,im,0,i4) = COEFF_NN(pnff,Z,0,0);
   COEFF_FE_NM(pel,Z,im,0,i5) = COEFF_NE_NM(pel->sons[i],Z,0,0,0);

   for (pnff = ni->tnfstart; pnff->nbface != f112; pnff = pnff->next);   
   COEFF_EF_MN(pel,Z,i,0,i12) = COEFF_NN(pnff,Z,0,0);
   for (pnff = ni->tnfstart; pnff->nbface != f411; pnff = pnff->next);   
   COEFF_EF_MN(pel,Z,im,0,i14) = COEFF_NN(pnff,Z,0,0);
   for (pnff = ni->tnfstart; pnff->nbface != f12i; pnff = pnff->next);   
   COEFF_EE_MM(pel,Z,0,i1) = COEFF_NN(pnff,Z,0,0);
   for (pnff = ni->tnfstart; pnff->nbface != f41i; pnff = pnff->next);   
   COEFF_EE_MM(pel,Z,0,i4) = COEFF_NN(pnff,Z,0,0);
   COEFF_EE_MM(pel,Z,0,i5) = COEFF_NE_NM(pel->sons[i],Z,3,0,0);

   for (pfnf = f112->tfnstart; pfnf->nbnode != n11; pfnf = pfnf->next);
   for (pfnc = f1->tfnstart; pfnc->nbnode != n1; pfnc = pfnc->next);
   COEFF_NN(pfnc,Z,i12,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f112->tfnstart; pfnf->nbnode != n41; pfnf = pfnf->next);
   for (pflc = f1->tfstart; pflc->nbface != f4; pflc = pflc->next);
   COEFF_NN(pflc,Z,i12,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f112->tfnstart; pfnf->nbnode != n12; pfnf = pfnf->next);
   COEFF_NN(f1,Z,i12,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f112->tfnstart; pfnf->nbnode != ni; pfnf = pfnf->next);
   COEFF_FE_NM(pel,Z,i,i12,0) = COEFF_NN(pfnf,Z,0,0);
   COEFF_NN(f1,Z,i12,i12) = COEFF_NN(f112,Z,0,0);
   for (pflf = f112->tfstart; pflf->nbface != f411; pflf = pflf->next);
   for (pflc = f1->tfstart; pflc->nbface != f4; pflc = pflc->next);
   COEFF_NN(pflc,Z,i12,i14) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f112->tfstart; pflf->nbface != f12i; pflf = pflf->next);
   COEFF_FE_NM(pel,Z,i,i12,i1) = COEFF_NN(pflf,Z,0,0);   
   for (pflf = f112->tfstart; pflf->nbface != f41i; pflf = pflf->next);
   COEFF_FE_NM(pel,Z,i,i12,i4) = COEFF_NN(pflf,Z,0,0);   
   COEFF_FE_NM(pel,Z,i,i12,i5) = COEFF_FE_NM(pel->sons[i],Z,1,0,0);

   for (pfnf = f411->tfnstart; pfnf->nbnode != n11; pfnf = pfnf->next);
   for (pfnc = f4->tfnstart; pfnc->nbnode != n1; pfnc = pfnc->next);
   COEFF_NN(pfnc,Z,i14,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f411->tfnstart; pfnf->nbnode != n12; pfnf = pfnf->next);
   for (pflc = f4->tfstart; pflc->nbface != f1; pflc = pflc->next);
   COEFF_NN(pflc,Z,i14,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f411->tfnstart; pfnf->nbnode != n41; pfnf = pfnf->next);
   COEFF_NN(f4,Z,i14,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f411->tfnstart; pfnf->nbnode != ni; pfnf = pfnf->next);
   COEFF_FE_NM(pel,Z,im,i14,0) = COEFF_NN(pfnf,Z,0,0);
   COEFF_NN(f4,Z,i14,i14) = COEFF_NN(f411,Z,0,0);
   for (pflf = f411->tfstart; pflf->nbface != f112; pflf = pflf->next);
   for (pflc = f4->tfstart; pflc->nbface != f1; pflc = pflc->next);
   COEFF_NN(pflc,Z,i14,i12) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f411->tfstart; pflf->nbface != f12i; pflf = pflf->next);
   COEFF_FE_NM(pel,Z,im,i14,i1) = COEFF_NN(pflf,Z,0,0);   
   for (pflf = f411->tfstart; pflf->nbface != f41i; pflf = pflf->next);
   COEFF_FE_NM(pel,Z,im,i14,i4) = COEFF_NN(pflf,Z,0,0);   
   COEFF_FE_NM(pel,Z,im,i14,i5) = COEFF_FE_NM(pel->sons[i],Z,0,0,0);

   for (pfnf = f41i->tfnstart; pfnf->nbnode != n11; pfnf = pfnf->next);
   COEFF_EN_MN(pel,Z,i,i4,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f41i->tfnstart; pfnf->nbnode != n41; pfnf = pfnf->next);
   COEFF_EF_MN(pel,Z,im,i4,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f41i->tfnstart; pfnf->nbnode != n12; pfnf = pfnf->next);
   COEFF_EF_MN(pel,Z,i,i4,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f41i->tfnstart; pfnf->nbnode != ni; pfnf = pfnf->next);
   COEFF_EE_MM(pel,Z,i4,0) = COEFF_NN(pfnf,Z,0,0);
   COEFF_EE_MM(pel,Z,i4,i4) = COEFF_NN(f41i,Z,0,0);
   for (pflf = f41i->tfstart; pflf->nbface != f112; pflf = pflf->next);
   COEFF_EF_MN(pel,Z,i,i4,i12) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f41i->tfstart; pflf->nbface != f411; pflf = pflf->next);
   COEFF_EF_MN(pel,Z,im,i4,i14) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f41i->tfstart; pflf->nbface != f12i; pflf = pflf->next);
   COEFF_EE_MM(pel,Z,i4,i1) = COEFF_NN(pflf,Z,0,0);
   COEFF_EE_MM(pel,Z,i4,i5) = COEFF_FE_NM(pel->sons[i],Z,3,0,0);

   for (pfnf = f12i->tfnstart; pfnf->nbnode != n11; pfnf = pfnf->next);
   COEFF_EN_MN(pel,Z,i,i1,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f12i->tfnstart; pfnf->nbnode != n41; pfnf = pfnf->next);
   COEFF_EF_MN(pel,Z,im,i1,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f12i->tfnstart; pfnf->nbnode != n12; pfnf = pfnf->next);
   COEFF_EF_MN(pel,Z,i,i1,0) = COEFF_NN(pfnf,Z,0,0);
   for (pfnf = f12i->tfnstart; pfnf->nbnode != ni; pfnf = pfnf->next);
   COEFF_EE_MM(pel,Z,i1,0) = COEFF_NN(pfnf,Z,0,0);
   COEFF_EE_MM(pel,Z,i1,i1) = COEFF_NN(f12i,Z,0,0);
   for (pflf = f12i->tfstart; pflf->nbface != f112; pflf = pflf->next);
   COEFF_EF_MN(pel,Z,i,i1,i12) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f12i->tfstart; pflf->nbface != f411; pflf = pflf->next);
   COEFF_EF_MN(pel,Z,im,i1,i14) = COEFF_NN(pflf,Z,0,0);
   for (pflf = f12i->tfstart; pflf->nbface != f41i; pflf = pflf->next);
   COEFF_EE_MM(pel,Z,i1,i4) = COEFF_NN(pflf,Z,0,0);
   COEFF_EE_MM(pel,Z,i1,i5) = COEFF_FE_NM(pel->sons[i],Z,2,0,0);

   COEFF_EF_MN(pel,Z,im,i5,0) = COEFF_EN_MN(pel->sons[i],Z,0,0,0);
   COEFF_EN_MN(pel,Z,i,i5,0)  = COEFF_EN_MN(pel->sons[i],Z,1,0,0);
   COEFF_EF_MN(pel,Z,i,i5,0)  = COEFF_EN_MN(pel->sons[i],Z,2,0,0);
   COEFF_EE_MM(pel,Z,i5,0)    = COEFF_EN_MN(pel->sons[i],Z,3,0,0);
   COEFF_EF_MN(pel,Z,im,i5,i14) = COEFF_EF_MN(pel->sons[i],Z,0,0,0);
   COEFF_EF_MN(pel,Z,i,i5,i12)  = COEFF_EF_MN(pel->sons[i],Z,1,0,0);
   COEFF_EE_MM(pel,Z,i5,i1)     = COEFF_EF_MN(pel->sons[i],Z,2,0,0);
   COEFF_EE_MM(pel,Z,i5,i4)     = COEFF_EF_MN(pel->sons[i],Z,3,0,0);
   COEFF_EE_MM(pel,Z,i5,i5)     = COEFF_EE_MM(pel->sons[i],Z,0,0);
}

void compute_indices(n11,n22,n33,n44,i12,i21,i23,i32,i34,i43,i14,i41)
NODE *n11, *n22, *n33, *n44;
INT *i12, *i21, *i23, *i32, *i34, *i43, *i14, *i41;
{
   *i12 = *i21 = *i23 = *i32 = *i34 = *i43 = *i14 = *i41 = 1;
   if (n11->index < n22->index)
      *i21 = 2;
   else
      *i12 = 2;
   if (n22->index < n33->index)
      *i32 = 2;
   else
      *i23 = 2;
   if (n33->index < n44->index)
      *i43 = 2;
   else
      *i34 = 2;
   if (n44->index < n11->index)
      *i14 = 2;
   else
      *i41 = 2;
}

void store_q2_on_four_el(pel,Z)
ELEMENT *pel;
INT Z;
{
   NODE *n1, *n2, *n3, *n4, *n11, *n22, *n33, *n44, *n12, *n23, *n34, *n41, *ni;
   FACE *f112, *f122, *f223, *f233, *f334, *f344, *f441, *f411, 
        *f12i, *f23i, *f34i, *f41i, *f1, *f2, *f3, *f4;
   INT i12, i21, i23, i32, i34, i43, i14, i41;
 
   NODES_OF_4ELEMENT(n1,n2,n3,n4,pel);
   FACES_OF_4ELEMENT(f1,f2,f3,f4,pel);
   FACES_OF_4ELEMENT(f411,f112,f12i,f41i,pel->sons[0]);
   FACES_OF_4ELEMENT(f122,f223,f23i,f12i,pel->sons[1]);
   FACES_OF_4ELEMENT(f233,f334,f34i,f23i,pel->sons[2]);
   FACES_OF_4ELEMENT(f344,f441,f41i,f34i,pel->sons[3]);
   NODES_OF_4ELEMENT(n41,n11,n12,ni,pel->sons[0]);
   NODES_OF_4ELEMENT(n23,n33,n34,ni,pel->sons[2]);
   n22 = n2->son;
   n44 = n4->son;
   compute_indices(n11,n22,n33,n44,&i12,&i21,&i23,&i32,&i34,&i43,&i14,&i41);
   store_q2_on_fine_el(Z,n41,n11,n12,ni,f411,f112,f12i,f41i,
                       f4,n1,f1,pel,0,i12,i14);
   store_q2_on_fine_el(Z,n12,n22,n23,ni,f122,f223,f23i,f12i,
                       f1,n2,f2,pel,1,i23,i21);
   store_q2_on_fine_el(Z,n23,n33,n34,ni,f233,f334,f34i,f23i,
                       f2,n3,f3,pel,2,i34,i32);
   store_q2_on_fine_el(Z,n34,n44,n41,ni,f344,f441,f41i,f34i,
                       f3,n4,f4,pel,3,i41,i43);
   COEFF_EE_MM(pel,Z,0,0) = COEFF_NN(ni,Z,0,0);
}

void store_q2_matrix_on_coarser_grid(tGrid,Z) /* tGrid is the coarse grid */
GRID *tGrid;
INT Z;
{
   ELEMENT *pel;

   set_mat_value(tGrid,Z,0.,0,0,U_TYPE,U_TYPE,A_STRUCT);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      store_q2_on_four_el(pel,Z);
}

void store_q2_vector_on_coarser_grid(tGrid,u) /* tGrid is the coarse grid */
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   NODE *pnode;
   INT i12, i21, i23, i32, i34, i43, i14, i41;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDMV(pnode,u,0) = NDMV(pnode->son,u,0); 
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      compute_indices(pel->sons[0]->n[1],pel->sons[1]->n[1],
                      pel->sons[2]->n[1],pel->sons[3]->n[1],
                      &i12,&i21,&i23,&i32,&i34,&i43,&i14,&i41);
      FDMV(pel->f[0],u,0) = NDMV(pel->sons[0]->n[2],u,0);
      FDMV(pel->f[1],u,0) = NDMV(pel->sons[1]->n[2],u,0);
      FDMV(pel->f[2],u,0) = NDMV(pel->sons[2]->n[2],u,0);
      FDMV(pel->f[3],u,0) = NDMV(pel->sons[3]->n[2],u,0);
      FDMV(pel->f[3],u,i14) = FDMV(pel->sons[0]->f[0],u,0);
      FDMV(pel->f[0],u,i21) = FDMV(pel->sons[1]->f[0],u,0);
      FDMV(pel->f[1],u,i32) = FDMV(pel->sons[2]->f[0],u,0);
      FDMV(pel->f[2],u,i43) = FDMV(pel->sons[3]->f[0],u,0);
      FDMV(pel->f[0],u,i12) = FDMV(pel->sons[0]->f[1],u,0);
      FDMV(pel->f[1],u,i23) = FDMV(pel->sons[1]->f[1],u,0);
      FDMV(pel->f[2],u,i34) = FDMV(pel->sons[2]->f[1],u,0);
      FDMV(pel->f[3],u,i41) = FDMV(pel->sons[3]->f[1],u,0);
      EDMV(pel,u,0) = NDMV(pel->sons[0]->n[3],u,0);
      EDMV(pel,u,1) = FDMV(pel->sons[0]->f[2],u,0);
      EDMV(pel,u,2) = FDMV(pel->sons[1]->f[2],u,0);
      EDMV(pel,u,3) = FDMV(pel->sons[2]->f[2],u,0);
      EDMV(pel,u,4) = FDMV(pel->sons[3]->f[2],u,0);
      EDMV(pel,u,5) = EDMV(pel->sons[0],u,0);
      EDMV(pel,u,6) = EDMV(pel->sons[1],u,0);
      EDMV(pel,u,7) = EDMV(pel->sons[2],u,0);
      EDMV(pel,u,8) = EDMV(pel->sons[3],u,0);
   }
}

void store_q2_vector_on_finer_grid(tGrid,u) /* tGrid is the coarse grid */
GRID *tGrid;
INT u;
{
   ELEMENT *pel;
   NODE *pnode;
   INT i12, i21, i23, i32, i34, i43, i14, i41;

   for (pnode = FIRSTN(tGrid); pnode; pnode = pnode->succ)
      NDMV(pnode->son,u,0) = NDMV(pnode,u,0);
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      compute_indices(pel->sons[0]->n[1],pel->sons[1]->n[1],
                      pel->sons[2]->n[1],pel->sons[3]->n[1],
                      &i12,&i21,&i23,&i32,&i34,&i43,&i14,&i41);
      NDMV(pel->sons[0]->n[2],u,0) = FDMV(pel->f[0],u,0);
      NDMV(pel->sons[1]->n[2],u,0) = FDMV(pel->f[1],u,0);
      NDMV(pel->sons[2]->n[2],u,0) = FDMV(pel->f[2],u,0);
      NDMV(pel->sons[3]->n[2],u,0) = FDMV(pel->f[3],u,0);
      FDMV(pel->sons[0]->f[0],u,0) = FDMV(pel->f[3],u,i14);
      FDMV(pel->sons[1]->f[0],u,0) = FDMV(pel->f[0],u,i21);
      FDMV(pel->sons[2]->f[0],u,0) = FDMV(pel->f[1],u,i32);
      FDMV(pel->sons[3]->f[0],u,0) = FDMV(pel->f[2],u,i43);
      FDMV(pel->sons[0]->f[1],u,0) = FDMV(pel->f[0],u,i12);
      FDMV(pel->sons[1]->f[1],u,0) = FDMV(pel->f[1],u,i23);
      FDMV(pel->sons[2]->f[1],u,0) = FDMV(pel->f[2],u,i34);
      FDMV(pel->sons[3]->f[1],u,0) = FDMV(pel->f[3],u,i41);
      NDMV(pel->sons[0]->n[3],u,0) = EDMV(pel,u,0);
      FDMV(pel->sons[0]->f[2],u,0) = EDMV(pel,u,1);
      FDMV(pel->sons[1]->f[2],u,0) = EDMV(pel,u,2);
      FDMV(pel->sons[2]->f[2],u,0) = EDMV(pel,u,3);
      FDMV(pel->sons[3]->f[2],u,0) = EDMV(pel,u,4);
      EDMV(pel->sons[0],u,0)       = EDMV(pel,u,5);
      EDMV(pel->sons[1],u,0)       = EDMV(pel,u,6);
      EDMV(pel->sons[2],u,0)       = EDMV(pel,u,7);
      EDMV(pel->sons[3],u,0)       = EDMV(pel,u,8);
   }
}

void solve_q1_on_coarser_grid(tGrid) /* tGrid is the coarse grid */
GRID *tGrid;
{
   INT nj;

   store_q1_matrix_on_coarser_grid(tGrid,A);
   store_q1_vector_on_coarser_grid(tGrid,F);
   fill_Ap_and_Ai_general(tGrid,Ap,Ai,&nj,MAX_ROW,MAX_ENT,1,1,1);
   fill_Ax_for_general_matr(tGrid,A,Ap,Ai,Ax,1,1,1);
   make_vector_from_grid_data_general(tGrid,F,Rhs,1,1,1);
   solve_system_using_umfpack(nj,Ap,Ai,Ax,Rhs,Sol);
   make_grid_data_from_vector_general(tGrid,R,Sol,1,1,1);
   store_q1_vector_on_finer_grid(tGrid,R);
   copy(tGrid->finer,R,U,T_FOR_U,U_TYPE);
}

void solve_q2_on_coarser_grid(tGrid) /* tGrid is the coarse grid */
GRID *tGrid;
{
   INT nj;

   store_q2_matrix_on_coarser_grid(tGrid,A);
   store_q2_vector_on_coarser_grid(tGrid,F);
   fill_Ap_and_Ai_general(tGrid,Ap,Ai,&nj,MAX_ROW,MAX_ENT,1,3,9);
   fill_Ax_for_general_matr(tGrid,A,Ap,Ai,Ax,1,3,9);
   make_vector_from_grid_data_general(tGrid,F,Rhs,1,3,9);
   solve_system_using_umfpack(nj,Ap,Ai,Ax,Rhs,Sol);
   make_grid_data_from_vector_general(tGrid,R,Sol,1,3,9);
   store_q2_vector_on_finer_grid(tGrid,R);
   copy(tGrid->finer,R,U,T_FOR_U,U_TYPE);
}

#endif

#endif /*  DIM == 2  */

