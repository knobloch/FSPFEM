#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#include "umfpack.h"
#include "macros1.h"
#include "def-cd-opt.h"
#include "def_data.h"
#include "data_sizes.h"
#include "macros2.h"
#include "ifndef.h"
#include "data_str.h"
#include "quadr_rules.h"
#include "boundary_approx.h"
#include "triang.h"
#include "refinement.h"
#include "basis_functions.h"
#include "b_cond.h"
#include "init.h"
#include "estimators.h"
//#include "nonunif_ref.h"
#include "vector_op.h"
#include "b_copy.h"
#include "matrix_op.h"
#include "ILU.h"
#include "CG_and_GMRES.h"
#include "matrices3D.h"
#include "matrices2D.h"
#include "smoother.h"
#include "new_restrict.h"
#include "restrict.h"
#include "solver.h"
#include "errors.h"
#include "postprocessing.h"
#include "save.h"
#include "check.h"
#include "mtests.h"
#include "nc_mtests.h"
#include "nc_prolong.h"
#include "tests-cd-opt.h"

main()
{
   INT i,j;
   NODE *nodes, *pnode;
   FACE *faces, *pface;
   VERTEX *vertexes, *pvert;
   BPOINT *bpoints, *pbpoint;
   ELEMENT *elements, *pelement;
   LINK *links, *plink;
   ELINK *elinks, *pelink;
   NELINK *nelinks, *pnelink;
   FNLINK *fnlinks, *pfnlink;
   NFLINK *nflinks, *pnflink;
   FLINK *flinks, *pflink;
   SNODE *snodes, *psnode;
   SFACE *sfaces, *psface;
   P_NODE *p_nodes, *pp_node;
   P_FACE *p_faces, *pp_face;
   NLGLINK *nlglinks, *pnlglink;
   LGNLINK *lgnlinks, *plgnlink;
   LGFLINK *lgflinks, *plgflink;
   LGLGLINK *lglglinks, *plglglink;
   FLGLINK *flglinks, *pflglink;
   LGDATA *lgdatas, *plgdata;
   GRID *grids, *pgrid;
   MULTIGRID multigrid, *mg;
   clock_t start_time, stop_time;

   start_time = clock();
   pnode = nodes = (NODE *) calloc(MAXNODE, sizeof(NODE));
   pface = faces = (FACE *) calloc(MAXFACE, sizeof(FACE));
   pvert = vertexes = (VERTEX *) calloc(MAXVERT, sizeof(VERTEX));
   pbpoint = bpoints = (BPOINT *) calloc(MAXBPOINT, sizeof(BPOINT));
   pelement = elements = (ELEMENT *) calloc(MAXELEM, sizeof(ELEMENT));
   plink = links = (LINK *) calloc(MAXLINK, sizeof(LINK));
   pelink = elinks = (ELINK *) calloc(MAXELINK, sizeof(ELINK));
   pnelink = nelinks = (NELINK *) calloc(MAXNELINK, sizeof(NELINK));
   pfnlink = fnlinks = (FNLINK *) calloc(MAXFNLINK, sizeof(FNLINK)); 
   pnflink = nflinks = (NFLINK *) calloc(MAXNFLINK, sizeof(NFLINK));
   pflink = flinks = (FLINK *) calloc(MAXFLINK, sizeof(FLINK));
   psnode = snodes = (SNODE *) calloc(MAXSNODE, sizeof(SNODE));
   psface = sfaces = (SFACE *) calloc(MAXSFACE, sizeof(SFACE));
   pp_node = p_nodes = (P_NODE *) calloc(NS, sizeof(P_NODE));
   pp_face = p_faces = (P_FACE *) calloc(NS, sizeof(P_FACE));
   pnlglink = nlglinks = (NLGLINK *) calloc(MAXNLGLINK, sizeof(NLGLINK));
   plgnlink = lgnlinks = (LGNLINK *) calloc(MAXLGNLINK, sizeof(LGNLINK));
   plgflink = lgflinks = (LGFLINK *) calloc(MAXLGFLINK, sizeof(LGFLINK));
   plglglink = lglglinks = (LGLGLINK *) calloc(MAXLGLGLINK, sizeof(LGLGLINK));
   pflglink = flglinks = (FLGLINK *) calloc(MAXFLGLINK, sizeof(FLGLINK));
   plgdata = lgdatas = (LGDATA *) calloc(MAXLGDATA, sizeof(LGDATA));
   pgrid = grids = (GRID *) calloc(MAXLEVEL, sizeof(GRID));
 
   creat("detected_errors",0644);
   if ( pnode==NULL || pface==NULL || pvert==NULL || pbpoint == NULL ||
        pelement==NULL || psnode==NULL || psface==NULL || pp_node==NULL ||
        pp_face==NULL || plink==NULL || pelink==NULL || pnelink==NULL || 
        pfnlink==NULL || pnflink==NULL || 
        pflink==NULL || pnlglink==NULL || plgnlink==NULL || plgflink==NULL || 
        plglglink==NULL || pflglink==NULL || plgdata==NULL || pgrid==NULL ) 
      eprintf("Unsufficient memory!\n");
   else{
      i = ALLOCATED_MEMORY;
      printf("%i bytes (%3.2f MB) allocated.\n",i,i/1048576.);
      mg = &multigrid;
      FIRSTGRID(mg) = grids;
      
      for (NV = 3; NV < 36; NV++)
         if(NV == NV_POINTS){
            NVX = NVY = NV;
//            NVX = 21;
//            NVY = 21;
            pnode = nodes;
            pface = faces;
            pvert = vertexes;
            pbpoint = bpoints;
            pelement = elements;
            plink = links;
            pelink = elinks;
            pnelink = nelinks;
            pfnlink = fnlinks; 
            pnflink = nflinks;
            pflink = flinks;
            psnode = snodes;
            psface = sfaces;
            pp_node = p_nodes;
            pp_face = p_faces;
            pnlglink = nlglinks;
            plgnlink = lgnlinks;
            plgflink = lgflinks;
            plglglink = lglglinks;
            pflglink = flglinks;
            plgdata = lgdatas;
            pgrid = grids;

            if (USE_COARSE_F == YES)
               first_level(read_vertexes_and_elements,
                      mg,&pnode,&pvert,&pbpoint,&pelement,&pface,&plink,&pelink,
                      &pnelink,&pfnlink,&pnflink,&pflink,
                      &psnode,&psface,&pp_node,&pp_face,&pnlglink,
                      &plgnlink,&plgflink,&plglglink,&pflglink,&plgdata);
            else
               first_level(COARSE_GRID,
                      mg,&pnode,&pvert,&pbpoint,&pelement,&pface,&plink,&pelink,
                      &pnelink,&pfnlink,&pnflink,&pflink,
                      &psnode,&psface,&pp_node,&pp_face,&pnlglink,
                      &plgnlink,&plgflink,&plglglink,&pflglink,&plgdata);
            size_of_data(mg,nodes,snodes,vertexes,bpoints,elements,faces,
                  sfaces,links,elinks,nelinks,flinks,nflinks,fnlinks,pnode,
                  psnode,pvert,pbpoint,pelement,pface,psface,
                  plink,pelink,pnelink,pflink,pnflink,
                  pfnlink,pnlglink,nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,
                  plglglink,lglglinks,pflglink,flglinks,plgdata,lgdatas);
                                
            for (i=1; i < REFINEMENTS + 1; i++){
               printf("\n\n");
               printf("Refinement:\n\n");
                              
               estimate1(mg);
               refinement(mg,&pnode,&pface,&pvert,&plink,&pelink,&pnelink,
                  &pnflink,&pfnlink,&pflink,&psnode,&psface,&pp_node,&pp_face,
                  &pelement,&pnlglink,
                  &plgnlink,&plgflink,&plglglink,&pflglink,&plgdata,REF_TYPE);
               size_of_data(mg,nodes,snodes,vertexes,bpoints,elements,faces,
                  sfaces,links,elinks,nelinks,flinks,nflinks,fnlinks,pnode,
                  psnode,pvert,pbpoint,pelement,pface,psface,
                  plink,pelink,pnelink,pflink,pnflink,
                  pfnlink,pnlglink,nlglinks,plgnlink,lgnlinks,plgflink,lgflinks,
                  plglglink,lglglinks,pflglink,flglinks,plgdata,lgdatas);
            }
            check_c_midpoints(mg);
            i = 1;
            if (INCLUDE_UMFPACK == YES)
               fill_Ap_and_Ai(TOP_GRID(mg),Ap,Ai,&Nj,
                                                MAX_ROW,MAX_ENT,U_TYPE,U_SPACE);
//          if (DELTA_TYPE == OUTFLOW_D)
//             set_supg_tau(TOP_GRID(mg),TNU,bb0,bb1,OUT_MASK);
//            solve_conv_diff(mg,t0,t01,t02,t03,ft0,YES,YES,YES,&i,1,1.e-8);
//            iterate_SDFEM(mg,t0,t01,t02,ft0);
            optimize_stab_par(mg,USE_BEL);
//            zkus_to(TOP_GRID(mg));
//            test_edge_stabilization(mg,ft0);
         }
   }
//   save_finest_triangulation(mg,vertexes,pvert);
   if (SAVE_TRIANG == YES)
      save_triangulation_for_gnuplot(mg,"gnutr");
/*   save_boundary_for_gnuplot(mg,"gnub"); */
/*   exact_solution_graph_for_gnuplot(TOP_GRID(mg),t0,"exact.gnu");  */
/*   save_exact_diag_profile(TOP_GRID(mg),t0,"exact_diag.gnu");  */
   stop_time = clock();
   printf("start_time = %i, stop_time = %i, CLOCKS_PER_SEC = %i\n",
           start_time,stop_time,CLOCKS_PER_SEC);
   printf("time taken = %i secs.\n",(stop_time - start_time)/CLOCKS_PER_SEC);
   close_pipe_for_gnuplot();
   print_number_of_errors();
   return(0);
}
