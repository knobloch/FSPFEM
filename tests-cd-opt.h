/******************************************************************************/
/*                                                                            */
/*                                   tests                                    */
/*                                                                            */
/******************************************************************************/

void supg_based_sold_tau(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], bbo[DIM], max=0., res;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (pel->eflag == 1)
         EDV(pel,x,1) = 0.;
      else{
         barycentric_coordinates(pel->n[0]->myvertex->x,
                                 pel->n[1]->myvertex->x,
                                 pel->n[2]->myvertex->x,bar);
         coord_of_barycentre(pel,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         ORT_VECT(bbo,bb)
         sp1c_grad(pel->n[0],pel->n[1],pel->n[2],bar,U,grad);
         res = DOT(bb,grad) +
               react(xc)*(NDS(pel->n[0],U)+NDS(pel->n[1],U)
                                          +NDS(pel->n[2],U))/3. - ft0(xc);
         EDV(pel,x,1) = 
          fabs(EDV(pel,x,0)-0.1*EDV(pel,W,0))*(fabs(res)+fabs(DOT(bbo,grad)));
//        fabs(res)+fabs(DOT(bbo,grad));
//        fabs(res)*fabs(DOT(bbo,grad));
         max = MAX(EDV(pel,x,1),max);
      }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      EDV(pel,x,1) *= 0.3/max;
}

void mult1(tGrid,r,x)
GRID *tGrid;
DOUBLE r;
INT x;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      EDV(pel,x,1) *= r;
}

void mult1s(tGrid,r,x,y)
GRID *tGrid;
DOUBLE r;
INT x, y;
{
   ELEMENT *pel;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], bbo[DIM], max=0., res, q;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (pel->eflag != 1){
         barycentric_coordinates(pel->n[0]->myvertex->x,
                                 pel->n[1]->myvertex->x,
                                 pel->n[2]->myvertex->x,bar);
         coord_of_barycentre(pel,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         ORT_VECT(bbo,bb)
         sp1c_grad(pel->n[0],pel->n[1],pel->n[2],bar,U,grad);
         res = DOT(bb,grad) +
               react(xc)*(NDS(pel->n[0],U)+NDS(pel->n[1],U)
                                          +NDS(pel->n[2],U))/3. - ft0(xc);
         EDV(pel,y,1) = fabs(res)+fabs(DOT(bbo,grad));
         max = MAX(EDV(pel,y,1),max);
      }
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (pel->eflag != 1){
         q = r*EDV(pel,y,1)/max;
         if (q > 1.)
            EDV(pel,x,1) *= q;
      }
}

void norm_of_param(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;
   FLOAT sum0=0., sum1=0.;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      sum0 += EDV(pel,x,0)*EDV(pel,x,0);
      sum1 += EDV(pel,x,1)*EDV(pel,x,1);
   }
   printf("norm0 = %17.11e, norm1 = %17.11e\n",sqrt(sum0), sqrt(sum1));
}

void compute_scaling_par();

DOUBLE solve_conv_diff_lp(tGrid,u0,u01,u02,u03,rhs,
                                  out_file,out_gnu,out_err,itercnt,max_iter,tol)
GRID *tGrid;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*u03)(), (*rhs)(), tol;
INT out_file, out_gnu, out_err, *itercnt, max_iter;
{
   FLOAT err=1., err1, err2, err3, err4, err6, err7, def;
   FLOAT u_rhs, area;
   INT i=1;
   FILE *fp;

   set_mat_value(tGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
   add_Laplace_matr(tGrid,TNU,A,U_SPACE,A_STRUCT,U_STRUCTURE,KORN_LAPLACE);
   add_convective_term_matr(tGrid,A,bb0,bb1,div_bb,CONV_DISCR,
                                                  U_SPACE,A_STRUCT,U_STRUCTURE);
   add_react_term_matr(tGrid,A,react,U_SPACE,A_STRUCT,U_STRUCTURE);
   add_macro_conv_stab_term_matr(tGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
   if (U_SPACE == GQ2B3C)
      partial_set_mat_to_zero_for_q2b3(tGrid,A);
   integrate_rhs(tGrid,F,rhs,rhs,rhs,
                                 T_FOR_U,U_TYPE,U_SPACE,U_STRUCTURE,RHS_INTEGR);
   subtract_Acontribution_from_bc(tGrid,A,U,F,Q,R,D,T_FOR_U,U_TYPE,A_STRUCT);
   set_value(tGrid,0.,D,T_FOR_U,U_TYPE);
   subtract_local_proj(tGrid,TNU,bb0,bb1,U,D,U_SPACE);
   subtr(tGrid,F,D,F,T_FOR_U,U_TYPE);
   if (U_SPACE == GQ2B3C){
      partial_set_vect_to_zero_for_q2b3(tGrid,U);
      partial_set_vect_to_zero_for_q2b3(tGrid,F);
   }
{

/*
   if (METHOD == BASIC_GMRES)
      GMRES(tGrid,U,F,3,0,100000,R,15,1.e-8,1.e-8,NONZERO_INIT,
      mult_A_lp,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,U_SPACE,0.,0.,0.);
   else if (METHOD == ILU_GMRES){
*/
      make_ILU(tGrid,A,1,0,ILU_TYPE,ILU_STRUCT);
      PGMRES(tGrid,U,F,D,Q,3,0,100000,R,15,1.e-10,1.e-10,NONZERO_INIT,
       mult_A_lp,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,U_SPACE,0.,0.,0.,
       1,A_STRUCT,0,0,ILU_PR);
/*
   }
   else if (METHOD == UMFPACK){
         fill_Ax(tGrid,A,Ap,Ai,Ax,U_TYPE,U_SPACE);
         make_vector_from_grid_data(tGrid,F,Rhs,U_TYPE,U_SPACE);
         solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
         make_grid_data_from_vector(tGrid,U,Sol,U_TYPE,U_SPACE);
   }
*/
      mult_A_lp(tGrid,A,U,D,F,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
                0,1,2,3,4,5,6,7,8,U_SPACE,0.,0.,0.);
      subtr(tGrid,D,F,D,T_FOR_U,U_TYPE);
      if (U_SPACE == GQ2B3C)
         partial_set_vect_to_zero_for_q2b3(tGrid,D);
      err4 = sqrt(dot(tGrid,D,D,T_FOR_U,U_TYPE)/
                  dot(tGrid,F,F,T_FOR_U,U_TYPE));
      printf("def/rhs: %e\n",err4);

      (*itercnt)++;
   }

/*
      if (out_err == YES){
         sL2_error(tGrid,U,R,&err1,u0,U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
         sH10_error(tGrid,U,R,&err2,u01,u02,u03,U_SPACE);
         sL_infinity_error(tGrid,U,R,&err6,u0,U_SPACE);
         g_conv_diff_ref_values(tGrid,u0,u01,u02,u03,U,D,R,
                                            U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
      }
*/
   return(def);
}

void iterate_SDFEM_lp(mg,u0,u01,u02,rhs)
MULTIGRID *mg;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*rhs)();
{
   GRID *tGrid, *topGrid;
   FLOAT tol=1.e-8, def=0.;
   INT itercnt=1, max_iter=15000, it=1;

//    tGrid = FIRSTGRID(mg);
      tGrid = TOP_GRID(mg);
      topGrid = TOP_GRID(mg);
      TNU = VALUE_OF_EPS;
      boundary_values(mg,u0,u0,u0,U,D,Q,U_SPACE,U_STRUCTURE);
      boundary_values(mg,u0,u0,u0,UU,D,Q,U_SPACE,U_STRUCTURE);
      if (INITIAL == ZERO)
         set_value(tGrid,0.,U,T_FOR_U,U_TYPE);
      else if (INITIAL == SOLUTION)
         initialize(mg,U,u0,u0,u0,U_SPACE,U_STRUCTURE);
      set_value(tGrid,0.,UU,T_FOR_U,U_TYPE);

   solve_conv_diff_lp(tGrid,u0,u01,u02,u02,rhs,NO,YES,YES,&itercnt,1,tol);
   solution_graph_for_gnuplot(tGrid,U,0,"sol_graph.gnu",U_SPACE,U_STRUCTURE);
// evaluate_ex_55_half(tGrid,0.5,1./32.,TAU0,U,u0,
//                     "half_cut.gnu","err.gnu","tot_var.gnu");
// check_symmetry(tGrid,U);
/*
      sn_grid_data(tGrid,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph.grid_gnu");
*/
}

DOUBLE solve_conv_diff(mg,u0,u01,u02,u03,rhs,
                                  out_file,out_gnu,out_err,itercnt,max_iter,tol)
MULTIGRID *mg;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*u03)(), (*rhs)(), tol;
INT out_file, out_gnu, out_err, *itercnt, max_iter;
{
   GRID *topGrid=TOP_GRID(mg), *tGrid;
   FLOAT err=1., err1, err2, err3, err4, err44, err55, err6, err7, def;
   FLOAT u_rhs, area, max_h=0.;
   INT i=1;
   FILE *fp;

   if (out_file == YES){
      fp = fopen(F_NAME,"w");
      fprintf(fp,"bc2sd.h (u = x*x + 0.1). Cube with refinement. Convection-difusion equation. \n\n");
      if (U_SPACE == P1_MOD)
         fprintf(fp,"Modified nonconforming P1 element.\n\n");
      else if (U_SPACE == P1_NC)
         fprintf(fp,"Nonconforming P1 element.\n\n");
      else if (U_SPACE == P1C)
         fprintf(fp,"Conforming P1 element.\n\n");
      else if (U_SPACE == P2C)
         fprintf(fp,"Conforming P2 element.\n\n");
      fprintf(fp," h        eps       L2-err    H10-err   L2conv    SD-err    L-infty   def/rhs\n");
      fprintf(fp,"------------------------------------------------------");
      fprintf(fp,"------------------------\n");
      fclose(fp);
   }

   if (*itercnt == 1){
      max_h = max_diameter_in_grid(topGrid);
      TNU = VALUE_OF_EPS;
//    printf("TNU = %e,  max. diameter = %e\n",TNU,max_h);
      boundary_values(mg,u0,u0,u0,U,D,Q,U_SPACE,U_STRUCTURE);
      boundary_values(mg,u0,u0,u0,UU,D,Q,U_SPACE,U_STRUCTURE);
      if (INITIAL == ZERO)
         set_value(topGrid,0.,U,T_FOR_U,U_TYPE);
      else if (INITIAL == SOLUTION)
         initialize(mg,U,u0,u0,u0,U_SPACE,U_STRUCTURE);
      set_value(topGrid,0.,UU,T_FOR_U,U_TYPE);
   }
   if (PERIODIC_BC == YES)
      copy_1st_to_2nd_periodic_boundary(topGrid,U,U_TYPE);

   if (SOLVE_BY_IMPL_EULER == YES){
      set_mat_value(topGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
      add_mass_matr(topGrid,IE_TAU,A,U_SPACE,A_STRUCT,U_STRUCTURE);
      mult_A(topGrid,A,U,B,B,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT,
             0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   }

   if (COMPUTE_SC == YES && SC_TYPE== SCM_MH)
      set_value(topGrid,0.,D,T_FOR_U,U_TYPE);
   for (tGrid = topGrid; tGrid && (METHOD == MULTI_GRID || i); 
                                                          tGrid=tGrid->coarser){
      set_mat_value(tGrid,A,0.,T_FOR_U,T_FOR_U,U_TYPE,U_TYPE,A_STRUCT);
      add_Laplace_matr(tGrid,TNU,A,U_SPACE,A_STRUCT,U_STRUCTURE,KORN_LAPLACE);
      if (COMPUTE_SC == YES && SC_TYPE == SCM_MH){
         if (MAKE_AVERAGES == YES){
            make_averages(topGrid,U,UU);
            add_MH_conv_term_matr(tGrid,A,UU,D,bb0,bb1,react,rhs,U_SPACE);
         }
         else
            add_MH_conv_term_matr(tGrid,A,U,D,bb0,bb1,react,rhs,U_SPACE);
// P1C_tabata_upwind_convective_term(tGrid,A,Q,R,bb0,bb1);
// P1C_upwind_convective_term(tGrid,TNU,A,bb0,bb1,barycentric_fluxes,heaviside);
// P1C_upwind_convective_term(tGrid,TNU,A,bb0,bb1,circumcentric_fluxes,heaviside);
// P1C_upwind_t_convective_term(tGrid,TNU,A,U,UU,bb0,bb1,barycentric_fluxes,heaviside);
// P1C_upwind_t_convective_term(tGrid,TNU,A,U,UU,bb0,bb1,circumcentric_fluxes,heaviside);
      }
      else{
         add_convective_term_matr(tGrid,A,bb0,bb1,div_bb,CONV_DISCR,
                                                  U_SPACE,A_STRUCT,U_STRUCTURE);
         add_react_term_matr(tGrid,A,react,U_SPACE,A_STRUCT,U_STRUCTURE);
      }
      if (SDFEM == YES)
         add_sd_term_matr(tGrid,TNU,A,bb0,bb1,react,
                                                  U_SPACE,A_STRUCT,U_STRUCTURE);
      if (LOC_PROJECT == YES)
         add_lp_term_matr(tGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
      if (CONV_LOC_PROJECT == YES)
         add_conv_lp_term_matr(tGrid,TNU,A,bb0,bb1,U_SPACE,A_STRUCT,U_STRUCTURE);
      if (COMPUTE_SC == YES && SC_TYPE != SCM_MH && SC_TYPE != SCM_LP &&
                                      (*itercnt > 1 || is_linear_sc(SC_TYPE))){
         add_stab_matr(tGrid,TNU,A,U,bb0,bb1,react,rhs,U_SPACE,SC_TYPE,SC_BETA,
                                                                      SC_BETA2);
         if (SC_TYPE == SCM_CS)
            add_stab_matr(tGrid,TNU,A,U,bb0,bb1,react,rhs,U_SPACE,SCM_CS1,
                                                              SC_BETA,SC_BETA2);
      }
      i = 0;
   }
   integrate_rhs(topGrid,F,rhs,rhs,rhs,
                                 T_FOR_U,U_TYPE,U_SPACE,U_STRUCTURE,RHS_INTEGR);
   if (SDFEM == YES)
      add_sd_rhs(topGrid,F,TNU,rhs,rhs,rhs,bb0,bb1,bb1,
                                 T_FOR_U,U_TYPE,U_SPACE,U_STRUCTURE,RHS_INTEGR);
   if (COMPUTE_SC == YES && SC_TYPE== SCM_LP && 
                                        (*itercnt > 1 || is_linear_sc(SC_TYPE)))
      add_penalty_term_to_rhs(topGrid,F,U,U_SPACE,SC_BETA);
   if (COMPUTE_SC == YES && SC_TYPE== SCM_MH)
      add(topGrid,D,F,F,T_FOR_U,U_TYPE);
   subtract_Acontribution_from_bc(topGrid,A,U,F,Q,R,D,T_FOR_U,U_TYPE,A_STRUCT);
   if (PERIODIC_BC == YES){
      add_2nd_to_1st_periodic_boundary(topGrid,F,U_TYPE);
      set_zero_on_2nd_periodic_boundary(topGrid,F,U_TYPE);
      set_zero_on_2nd_periodic_boundary(topGrid,U,U_TYPE);
   }

   if ((COMPUTE_SC == YES || SOLVE_BY_IMPL_EULER == YES) && *itercnt > 1) {
      defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
      err = sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/
                 dot(topGrid,F,F,T_FOR_U,U_TYPE));
      printf("residuum/rhs: %10.5e\n", err);
      def = err;
   }

   if (*itercnt <= max_iter && err > tol){
      if (max_iter > 1)
         printf("\npass %d:\n",*itercnt);

      if (SOLVE_BY_IMPL_EULER == YES){
         add(topGrid,B,F,F,T_FOR_U,U_TYPE);
         for (tGrid = topGrid, i = 1; tGrid && (METHOD == MULTI_GRID || i); 
                                                          tGrid=tGrid->coarser){
            add_mass_matr(topGrid,IE_TAU,A,U_SPACE,A_STRUCT,U_STRUCTURE);
            i = 0;
         }
      }

      copy(topGrid,U,UU,T_FOR_U,U_TYPE);
      if (METHOD == BASIC_GMRES)
         GMRES(topGrid,U,F,3,0,100000,R,15,1.e-8,1.e-8,NONZERO_INIT,
               mult_A,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
      else if (METHOD == DIAG_GMRES){
         copy_diag_to_vector(topGrid,A,R,T_FOR_U,U_TYPE,A_STRUCT);
         PGMRES(topGrid,U,F,D,Q,3,0,100000,B,15,1.e-8,1.e-8,NONZERO_INIT,
                mult_A,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
                R,0,0,0,DIAG_PR);
      }
      else if (METHOD == ILU_GMRES){
         make_ILU(topGrid,A,1,0,ILU_TYPE,ILU_STRUCT);
         PGMRES(topGrid,U,F,D,Q,3,0,100000,R,15,1.e-8,1.e-8,NONZERO_INIT,
                mult_A,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
                1,A_STRUCT,0,0,ILU_PR);
      }
      else if (METHOD == MULTI_GRID){
         u_rhs = dot(topGrid,F,F,T_FOR_U,U_TYPE);
         err4 = 1.;
         for (i=0; i<200000  && (err4 > 1.e-8); i++){
            printf("%i\n",i);
            mgm(topGrid,1,W_CYCLE,5,5,YES,1.,A,1,F,U,D,Q,
                T_FOR_U,U_TYPE,U_SPACE,A_STRUCT,SCALAR,SOR_F);
            defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
            err4 = sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/u_rhs);
            printf("def/rhs: %e\n",err4);
         }
      }
      else if (METHOD == UMFPACK){
         fill_Ax(topGrid,A,Ap,Ai,Ax,U_TYPE,U_SPACE);
         make_vector_from_grid_data(topGrid,F,Rhs,U_TYPE,U_SPACE);
         solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
         make_grid_data_from_vector(topGrid,U,Sol,U_TYPE,U_SPACE);
      }
      if (COMPUTE_SC == YES && !is_linear_sc(SC_TYPE) && 
          SOLVE_BY_IMPL_EULER == NO)
         damp(topGrid,DAMPING_PARAM,UU,U,U,T_FOR_U,U_TYPE);

      defect(topGrid,A,F,U,D,F,T_FOR_U,U_TYPE,A_STRUCT);
      err4 = sqrt(dot(topGrid,D,D,T_FOR_U,U_TYPE)/
                  dot(topGrid,F,F,T_FOR_U,U_TYPE));
//    printf("def/rhs: %e\n",err4);

      if (PERIODIC_BC == YES)
         copy_1st_to_2nd_periodic_boundary(topGrid,U,U_TYPE);
      if (out_err == YES){
         sL2_error(topGrid,U,R,&err1,u0,U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
         sH10_error(topGrid,U,R,&err2,u01,u02,u03,U_SPACE);
         s_conv_error(topGrid,U,R,&err44,bb0,bb1,bb1,u01,u02,u03,U_SPACE);
         s_sd_error(topGrid,U,R,&err55,TNU,bb0,bb1,bb1,u0,u01,u02,u03,U_SPACE);
         sL_infinity_error(topGrid,U,R,&err6,u0,U_SPACE);
         if (LOC_PROJECT == YES)
            s_lp_error(topGrid,U,&err3,TNU,bb0,bb1,bb1,u01,u02,u03,U_SPACE);
         if (CONV_LOC_PROJECT == YES)
            s_conv_lp_error(topGrid,U,&err3,TNU,bb0,bb1,bb1,u01,u02,u03,U_SPACE);
         printf("lp err= %e\n",err3=sqrt(err3*err3+TNU*err2*err2));
/*
         mv_set_value_e(topGrid,0.,U);
         sL2_error(topGrid,U,R,&err4,u0,U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
         sH10_error(topGrid,U,R,&err6,u01,u02,u03,U_SPACE);
         if (LOC_PROJECT == YES)
            s_lp_error(topGrid,U,&err7,TNU,bb0,bb1,bb1,u01,u02,u03,U_SPACE);
         if (CONV_LOC_PROJECT == YES)
            s_conv_lp_error(topGrid,U,&err7,TNU,bb0,bb1,bb1,u01,u02,u03,U_SPACE);
         printf("lp err= %e\n",err7=sqrt(err7*err7+TNU*err6*err6));
         fp = fopen("errors.gnu","a");
         fprintf(fp,"%e %e %e %e %e %e %e\n",TAU0,err1,err2,err3,err4,err6,err7);
         fclose(fp);
*/
/*
         compute_errors_for_gert(topGrid,U,R,TAU0,TNU,u0,u01,u02,u03,
          T_FOR_U,U_TYPE,U_SPACE,U_STRUCTURE,ELEM,"errors.gnu","errors_we.gnu");
         g_conv_diff_ref_values(topGrid,u0,u01,u02,u03,U,D,R,
                                            U_SPACE,U_STRUCTURE,T_FOR_U,U_TYPE);
*/
      }
      if (out_file == YES){
         fp = fopen(F_NAME,"a");
         fprintf(fp,"%8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  \n",
                                     max_h,TNU,err1,err2,err44,err55,err6,err4);
         fprintf(fp,"%11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e\n",
                                     max_h,TNU,err1,err2,err44,err55,err6,err4);
      }
      if (out_err == YES && ERR_ON_SQUARE == YES){
         errors_on_square(topGrid,X_SQUARE,Y_SQUARE,TNU,bb0,bb1,bb1,
              u0,u01,u02,u03,U,R,&err7,&err2,&err44,&err55,&err3,&area,U_SPACE);
         if (out_file == YES){
            fprintf(fp,"Error on square (0,%4.2f)x(0,%4.2f)",X_SQUARE,Y_SQUARE);
            fprintf(fp,"   (computed area = %f)\n",area);
            fprintf(fp,"%8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  \n",
                                     max_h,TNU,err7,err2,err44,err55,err3,err4);
            fprintf(fp,"%11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e\n",
                                     max_h,TNU,err7,err2,err44,err55,err3,err4);
         }
      }
//      save_profile(topGrid,1,0.5,U,cfn3(*itercnt),U_SPACE);
      if (out_file == YES)
         fclose(fp);
      (*itercnt)++;
   }
   else
      *itercnt = max_iter + 2;
   if (PERIODIC_BC == YES)
      copy_1st_to_2nd_periodic_boundary(topGrid,U,U_TYPE);
   if (out_gnu == YES && (max_iter == 1 || *itercnt == max_iter + 2)){
      solution_graph_for_gnuplot(topGrid,U,0,"sol_graph.gnu",U_SPACE,
                                                                   U_STRUCTURE);
/*
      error_graph_for_gnuplot(topGrid,u0,U,0,0.1,"err_graph.gnu",U_SPACE,
                                                                   U_STRUCTURE);
      sn_grid_data(topGrid,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph.grid_gnu");
      save_profile(topGrid,0,1./(NV_POINTS-1.),U,"cut1_y.gnu",U_SPACE);
      save_profile(topGrid,0,1.-1./(NV_POINTS-1.),U,"cut2_y.gnu",U_SPACE);
      evaluate_ex_for_gert(topGrid,1./32.,1./64.,
                           TAU0,u0,U,"prof1.gnu","err1.gnu");
      evaluate_ex_for_gert(topGrid,1.-1./32.,1./64.,
                           TAU0,u0,U,"prof2.gnu","err2.gnu");
*/
      if ((U_SPACE == P1_MOD || U_SPACE == P1_NC)){
         if (U_SPACE == P1_MOD)
            vector_to_scalar_f(topGrid,U,U,0);
         nc_conf_averaging(topGrid,U,U,R,u0);
         solution_graph_for_gnuplot(topGrid,U,0,"n_sol_graph.gnu",P1C,
                                                                   U_STRUCTURE);
         error_graph_for_gnuplot(topGrid,u0,U,0,0.1,"n_err_graph.gnu",P1C,
                                                                   U_STRUCTURE);
      }
/*
      save_profile(topGrid,1,0.,U,"profile_x.gnu",U_SPACE);
      save_profile(topGrid,0,0.,U,"profile_y.gnu",U_SPACE);
      save_profile(topGrid,1,0.5,U,"cut_x.gnu",U_SPACE);
      save_profile(topGrid,0,0.5,U,"cut_y.gnu",U_SPACE);
      save_profile(topGrid,1,0.25,U,"cut25_x.gnu",U_SPACE);
      save_profile(topGrid,0,0.75,U,"cut75_y.gnu",U_SPACE);
      graph_on_cut(mg,1,0.25,U,"cut25_x_new.gnu",U_SPACE,1001);
      if (U_SPACE == P1C){
         save_p1c_diag_profile(topGrid,U,"diag_profile.gnu");
         save_p1c_diag2_profile(topGrid,U,"diag2_profile.gnu");
      }
      evaluate_parabolic_bdry_layers_for_p1c(mg,topGrid,U,R,err7,u0,
                                             T_FOR_U,U_TYPE,U_STRUCTURE);
      evaluate_skew_convection_for_p1c(mg,topGrid,U);
*/
   }
   return(def);
}

/*   
   CG(topGrid,F,U,D,Q,R,0,100,0.0001,EPS,0,mult_A,A,
       R,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   PCG(topGrid,F,U,D,Q,R,0,1000,0.1,EPS,0,mult_A,A,
       R,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
       1,ILU_STRUCT,0,0,ILU_PR);
*/

void iterate_SDFEM(mg,u0,u01,u02,rhs)
MULTIGRID *mg;
FLOAT (*u0)(), (*u01)(), (*u02)(), (*rhs)();
{
   FLOAT tol=1.e-8, def=0.;
   INT itercnt=1, max_iter=15000, it=1;

   if ((COMPUTE_SC == NO || is_linear_sc(SC_TYPE)) &&
       (SOLVE_BY_IMPL_EULER == NO))
      solve_conv_diff(mg,u0,u01,u02,u02,rhs,NO,NO,NO,&itercnt,1,tol);
   else
     while(itercnt <= max_iter+1){
        it  = itercnt-1;
        def = solve_conv_diff(mg,u0,u01,u02,u02,rhs,NO,YES,YES,&itercnt,
                                                                  max_iter,tol);
     }
/*
   check_tabata(TOP_GRID(mg),U,0.5,-0.5*sqrt(3.));
   check_baba_tabata(TOP_GRID(mg),U,0.5,-0.5*sqrt(3.));
   check_kanayama(TOP_GRID(mg),U,0.5,-0.5*sqrt(3.));
   if (SC_EXAMPLE == 20)
      evaluate_hill(mg,U,U_SPACE,it,def);
   if (SC_EXAMPLE == 25)
      evaluate_hemker(mg,U,U_SPACE,it,def);
   else if (SC_EXAMPLE == 55 && USE_COARSE_F == NO)
      evaluate_parabolic_layers2(mg,U,it,def);
   else if (SC_EXAMPLE == 55 && USE_COARSE_F == YES)
     evaluate_parabolic_layers_for_unstructured_grid(mg,u0,U,U_SPACE,0.1);
   else if (SC_EXAMPLE == 88)
      evaluate_skew_convection2(mg,TOP_GRID(mg),U,U_SPACE,it,def,NO);
   else if (SC_EXAMPLE == 93)
      evaluate_klr(mg,U,U_SPACE,it,def);
*/
}

#if E_DATA & SCALAR_ELEMENT_DATA

void save_stab_par(tGrid,name,u)
GRID *tGrid;
char name[];
INT u;
{
   ELEMENT *pel;
   FILE *fp;

   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fprintf(fp,"%e\n",ED(pel,u));
   fclose(fp);
}

void read_stab_par(tGrid,name,u)
GRID *tGrid;
char name[];
INT u;
{
   ELEMENT *pel;
   FILE *fp;

   fp = fopen(name,"r");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fscanf(fp,"%lf",&ED(pel,u));
   fclose(fp);
}

#elif E_DATA & VECTOR_ELEMENT_DATA

void save_stab_par(tGrid,name,u)
GRID *tGrid;
char name[];
INT u;
{
   ELEMENT *pel;
   FILE *fp;

   fp = fopen(name,"w");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fprintf(fp,"%e %e\n",EDV(pel,u,0),EDV(pel,u,1));
   fclose(fp);
}

void read_stab_par(tGrid,name,u)
GRID *tGrid;
char name[];
INT u;
{
   ELEMENT *pel;
   FILE *fp;

   fp = fopen(name,"r");
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      fscanf(fp,"%lf %lf",&EDV(pel,u,0),&EDV(pel,u,1));
   fclose(fp);
}

#else

void save_stab_par(tGrid,name,u)
GRID *tGrid; char name[]; INT u;
{  eprintf("Error: save_stab_par not available.\n");  }

void read_stab_par(tGrid,name,u)
GRID *tGrid; char name[]; INT u;
{  eprintf("Error: read_stab_par not available.\n");  }

#endif

void compute_solution_of_adjoint_problem(tGrid,Z,Z1,
                       u,der,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;                      /* g is the solution of the adjoint problem */
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT Z, Z1, u, der, g, v, w, tot_der, use_bel; 
{
   derivatives_of_residual_based_error_estimator(tGrid,u,der,g,v,w,
                                            eps,bb0,bb1,react,rhs,beta,use_bel);
   smakeAT(tGrid,Z,Z1,1);
   if (METHOD == ILU_GMRES){
      make_ILU(tGrid,Z1,2,0,ILU_TYPE,ILU_STRUCT);
      PGMRES(tGrid,g,der,R,B,3,0,100000,W,15,1.e-8,1.e-8,0,
             mult_A,Z1,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
             2,A_STRUCT,0,0,ILU_PR);
   }
   else if (METHOD == UMFPACK){
      fill_Ax(tGrid,Z1,Ap,Ai,Ax,U_TYPE,U_SPACE);
      make_vector_from_grid_data(tGrid,der,Rhs,U_TYPE,U_SPACE);
      solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
      make_grid_data_from_vector(tGrid,g,Sol,U_TYPE,U_SPACE);
   }
   set_value(tGrid,0.,g,STOP_IS_FIRST_INNER,U_TYPE);
}

#if (E_DATA & SCALAR_ELEMENT_DATA) && (PAR_TYPE == Q_SE)

void total_derivative_of_residual_based_error_estimator(tGrid,Z,Z1,
                        u,der,d,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT Z, Z1, u, der, d, g, v, w, tot_der, use_bel; 
                /*  u ... solution at nodes, der ... derivatives at nodes  */
                /*  d, g ... auxiliary scalar node variable  */
                /*  v, w ... auxiliary scalar face variables  */
                /*  tot_der ... scalar element variable  */
{
   ELEMENT *pel;
   FLOAT b[DIM2][DIM2];

   compute_solution_of_adjoint_problem(tGrid,Z,Z1,
                       u,der,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel);
                                 /* g is the solution of the adjoint problem */
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      barycentric_coordinates(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                              pel->n[2]->myvertex->x,b);
      ED(pel,tot_der) = -cd_res_times_sd(pel,b,u,g,eps,bb0,bb1,react,rhs);
   }
}

#elif (E_DATA & VECTOR_ELEMENT_DATA) && (PAR_TYPE == Q_VE)

void total_derivative_of_residual_based_error_estimator(tGrid,Z,Z1,
                        u,der,d,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT Z, Z1, u, der, d, g, v, w, tot_der, use_bel; 
                /*  u ... solution at nodes, der ... derivatives at nodes  */
                /*  d, g ... auxiliary scalar node variable  */
                /*  v, w ... auxiliary scalar face variables  */
                /*  tot_der ... scalar element variable  */
{
   ELEMENT *pel;
   FLOAT b[DIM2][DIM2], bb[DIM], bbo[DIM], grad[DIM], xc[DIM], rdetB, q, res, tau;

   compute_solution_of_adjoint_problem(tGrid,Z,Z1,
                       u,der,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel);
                                 /* g is the solution of the adjoint problem */
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      rdetB = barycentric_coordinates(pel->n[0]->myvertex->x,
                                      pel->n[1]->myvertex->x,
                                      pel->n[2]->myvertex->x,b);
      EDV(pel,tot_der,0) = -cd_res_times_sd(pel,b,u,g,eps,bb0,bb1,react,rhs);
      if (COMPUTE_SC == YES){
         coord_of_barycentre(pel,xc);
         bb[0] = bb0(xc);
         bb[1] = bb1(xc);
         bbo[0] = -bb1(xc);
         bbo[1] =  bb0(xc);
         sp1c_grad(pel->n[0],pel->n[1],pel->n[2],b,u,grad);
         if ((q=DOT(grad,grad)) < 1.e-60)
            EDV(pel,tot_der,1) = 0.;
         else{
            res = DOT(bb,grad) +
                  react(xc)*(NDS(pel->n[0],u)+NDS(pel->n[1],u)
                                             +NDS(pel->n[2],u))/3. - ft0(xc);
            EDV(pel,tot_der,1) = 0.5*diameter(pel)*fabs(res)/sqrt(q);
            EDV(pel,tot_der,1) *= -DOT(bbo,grad)*rdetB/DOT(bbo,bbo);
            sp1c_grad(pel->n[0],pel->n[1],pel->n[2],b,g,grad);
            EDV(pel,tot_der,1) *= DOT(bbo,grad);
         }
//EDV(pel,tot_der,0) = 0.;
      }
      else
         EDV(pel,tot_der,1) = 0.;
      if (RESTRICT_DER == YES && (pel->eflag == 0 || pel->eflag == 2))
         EDV(pel,tot_der,0) = EDV(pel,tot_der,1) = 0.;
   }
}

#else

void total_derivative_of_residual_based_error_estimator(tGrid,Z,Z1,u,der,d,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid; FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta; INT Z, Z1, u, der, d, g, v, w, tot_der, use_bel; 
{  eprintf("Error: total_derivative_of_residual_based_error_estimator not available.\n");  }

#endif

/*
DOUBLE e_fcn_for_minim(mg,x,t,type)  // the parameters are in ED(pel,x)
MULTIGRID *mg;
INT x, t, type;
{
   ELEMENT *pel;
   DOUBLE cmin=1.e-5, a=1.e16, sum=0.;
// DOUBLE cmin=1.e-4, a=0., sum=0.;

   copy_ed_to_tau(TOP_GRID(mg),x);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel; pel = pel->succ)
      if (ED(pel,x) < cmin)
         sum += a*(ED(pel,x) - cmin)*(ED(pel,x) - cmin);
   return(sum + residual_based_error_estimator(TOP_GRID(mg),U,U,UU,TNU,bb0,bb1,
                                               react,ft0,0.,USE_BEL));
}

void e_grad_fcn_for_minim(mg,x,y,t,type)
MULTIGRID *mg;
INT x, y, t, type;
{
   ELEMENT *pel;
   DOUBLE cmin=1.e-5, a=1.e16, sum=0.;
// DOUBLE cmin=1.e-4, a=0., sum=0.;

   copy_ed_to_tau(TOP_GRID(mg),x);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   total_derivative_of_residual_based_error_estimator(TOP_GRID(mg),
                           A,1,U,D,R,Q,U,UU,y,TNU,bb0,bb1,react,ft0,0.,USE_BEL);
   for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel; pel = pel->succ)
      if (ED(pel,x) < cmin)
         ED(pel,y) += 2.*a*(ED(pel,x) - cmin);
}
*/

DOUBLE e_fcn_for_minim(mg,x,t,type)  /* the parameters are in ED(pel,x) */
MULTIGRID *mg;
INT x, t, type;
{
   copy_ed_to_tau(TOP_GRID(mg),x);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   return(residual_based_error_estimator(TOP_GRID(mg),U,U,UU,TNU,bb0,bb1,
                                                react,ft0,0.,USE_BEL));
}

void e_grad_fcn_for_minim(mg,x,y,t,type)
MULTIGRID *mg;
INT x, y, t, type;
{
   copy_ed_to_tau(TOP_GRID(mg),x);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   total_derivative_of_residual_based_error_estimator(TOP_GRID(mg),
                           A,1,U,D,R,Q,U,UU,y,TNU,bb0,bb1,react,ft0,0.,USE_BEL);
}

DOUBLE e_line_phi(mg,a,x,p,z,t,type)
MULTIGRID *mg;
DOUBLE a;
INT x, p, z, t, type;
{
   GRID *tGrid=TOP_GRID(mg);

   mult_and_add(tGrid,a,p,x,z,t,type);
   return(e_fcn_for_minim(mg,z,t,type));
}

DOUBLE e_der_line_phi(mg,a,x,p,y,z,t,type)
MULTIGRID *mg;
DOUBLE a;
INT x, p, y, z, t, type;
{
   GRID *tGrid=TOP_GRID(mg);

   mult_and_add(tGrid,a,p,x,z,t,type);
   e_grad_fcn_for_minim(mg,z,y,t,type);
   return(dot(tGrid,y,p,t,type));
}

DOUBLE e_zoom(mg,tGrid,zoom_it,a1,a2,phi1,phi2,dphi1,phi0,dphi0,phi3,dphi3,
            c1,c2,x,p,y,z,t,type)
MULTIGRID *mg;
GRID *tGrid;
DOUBLE a1, a2, phi1, phi2, dphi1, phi0, dphi0, *phi3, *dphi3, c1, c2;
INT zoom_it, x, p, y, z, t, type;
{
   DOUBLE c, a3;
   INT i=0, k=1;

   while (k && i < zoom_it){
//    a3 = argmin q(x) where q(x) = phi1 + dphi1*(x-a1) + c*(x-a1)*(x-a1) 
//    with q(a2)=phi2 
      c = (phi2 - phi1 - dphi1*(a2-a1))/((a2-a1)*(a2-a1));
      a3 = a1 - 0.5*dphi1/c;
      if (!( (a1 <= a3 && a3 <= a2) || (a2 <= a3 && a3 <= a1) ))
         a3 = 0.5*(a1 + a2);
      *phi3 =  e_line_phi(mg,a3,x,p,z,t,type);
      if (*phi3 > phi0 + c1*a3*dphi0 || *phi3 >= phi1){
         a2 = a3;
         phi2 = *phi3;
      }
      else{
         *dphi3 = e_der_line_phi(mg,a3,x,p,y,z,t,type);
         if (fabs(*dphi3) <= -c2*dphi0)
            k = 0;
         if ((*dphi3)*(a2-a1) >= 0.){
            a2 = a1;
            phi2 = phi1;
         }
         a1 = a3;
         phi1 = *phi3;
         dphi1 = *dphi3;
      }
      i++;
   }
   printf("%i iterations in zoom\n",i);
   return(a3);
}

INT e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,phi0_new,step,
                  c1,c2,x,p,y,z,t,type)
MULTIGRID *mg;
GRID *tGrid;
DOUBLE phi0, dphi0, *phi0_new, *step, c1, c2;
INT line_it, zoom_it, x, p, y, z, t, type;
{
   DOUBLE a1=0., a2=*step, phi1=phi0, phi2, dphi1=dphi0, dphi2;
   INT i=0, k=1;

   if (dphi0 >= 0.)
      eprintf("Error in e_line_search: not a descent direction.\n");
   while (k && i < line_it){
      phi2 = e_line_phi(mg,a2,x,p,z,t,type);
      if (phi2 > phi0 + c1*a2*dphi0 || (i && phi2 >= phi1)){
         a2 = e_zoom(mg,tGrid,zoom_it,a1,a2,phi1,phi2,dphi1,phi0,dphi0,
                   &phi2,&dphi2,c1,c2,x,p,y,z,t,type);
         k = 0;
      }
      else{
         dphi2 = e_der_line_phi(mg,a2,x,p,y,z,t,type);
         if (fabs(dphi2) <= -c2*dphi0)
            k = 0;
         else if (dphi2 >= 0){
               a2 = e_zoom(mg,tGrid,zoom_it,a2,a1,phi2,phi1,dphi2,phi0,dphi0,
                         &phi2,&dphi2,c1,c2,x,p,y,z,t,type);
               k = 0;
         }
         else{
            phi1 = phi2;
            dphi1 = dphi2;
            a1 = a2;
            a2 *= 2.;
         }
      }
      i++;
   }
printf("%i iterations in e_line_search\n",i);
   if (phi2 <= phi0 + c1*a2*dphi0 && fabs(dphi2) <= -c2*dphi0){
      *step = a2;
      *phi0_new = phi2;
      return(1);
   }
   else{
      eprintf("Error in e_line_search: parameter not found.\n");
      return(0);
   }
}

#if 0
DOUBLE max_factor_for_nonnegativity(tGrid,x,p)
GRID *tGrid;  /* components of x are nonnegative */
INT x, p;     /* x + factor*p should be nonnegative as well */
{
   ELEMENT *pel;
   DOUBLE factor=1.e150, a;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      if (ED(pel,p) < 0. && (a=-ED(pel,x)/ED(pel,p)) < factor)
         factor = a;
   return(factor); 
}
#endif

#if (E_DATA & SCALAR_ELEMENT_DATA) && (PAR_TYPE == Q_SE)

void take_positive_part(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      if (ED(pel,x) < 0.)
         ED(pel,x) = 0.;
/*
      if (ED(pel,x) > ED(pel,W))
         ED(pel,x) = ED(pel,W);
*/
   }
}

#elif (E_DATA & VECTOR_ELEMENT_DATA) && (PAR_TYPE == Q_VE)

void take_positive_part(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      if (EDV(pel,x,0) < 0.)
         EDV(pel,x,0) = 0.;
      if (EDV(pel,x,1) < 0.)
         EDV(pel,x,1) = 0.;
      if (EDV(pel,x,0) > EDV(pel,W,0))
         EDV(pel,x,0) = EDV(pel,W,0);
      if (EDV(pel,x,1) > 1.)
         EDV(pel,x,1) = 1.;
   }
}

void compute_scaling_par(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], q, res, tau;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      barycentric_coordinates(pel->n[0]->myvertex->x,
                              pel->n[1]->myvertex->x,
                              pel->n[2]->myvertex->x,bar);
      coord_of_barycentre(pel,xc);
      bb[0] = bb0(xc);
      bb[1] = bb1(xc);
      sp1c_grad(pel->n[0],pel->n[1],pel->n[2],bar,U,grad);
      q = DOT(grad,grad);
      if ((q=DOT(grad,grad)) < 1.e-60)
         EDV(pel,x,1) = 0.;
      else{
         res = DOT(bb,grad) +
               react(xc)*(NDS(pel->n[0],U)+NDS(pel->n[1],U)
                                          +NDS(pel->n[2],U))/3. - ft0(xc);
         EDV(pel,x,1) = 0.5*diameter(pel)*fabs(res)/sqrt(q);
      }
   }
}

void adjust_sold_par_to_h2(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;
   FLOAT q;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      q = diameter(pel)*diameter(pel);
      if (EDV(pel,x,1) > q)
         EDV(pel,x,1) = 1.;
      else
         EDV(pel,x,1) /= q;
   }
}

void adjust_sold_par_to_1(tGrid,x)
GRID *tGrid;
INT x;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
      EDV(pel,x,1) *= diameter(pel)*diameter(pel);
}

void check_scaling_par(tGrid,x,y)
GRID *tGrid;
INT x, y;
{
   ELEMENT *pel;
   FLOAT bar[DIM2][DIM2], grad[DIM], xc[DIM], bb[DIM], q, res, tau;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      barycentric_coordinates(pel->n[0]->myvertex->x,
                              pel->n[1]->myvertex->x,
                              pel->n[2]->myvertex->x,bar);
      coord_of_barycentre(pel,xc);
      bb[0] = bb0(xc);
      bb[1] = bb1(xc);
      sp1c_grad(pel->n[0],pel->n[1],pel->n[2],bar,U,grad);
      q = DOT(grad,grad);
      if ((q=DOT(grad,grad)) > 1.e-60){
         res = DOT(bb,grad) +
               react(xc)*(NDS(pel->n[0],U)+NDS(pel->n[1],U)
                                          +NDS(pel->n[2],U))/3. - ft0(xc);
         q = 0.5*diameter(pel)*fabs(res)/sqrt(q);
      }
//q = diameter(pel)*diameter(pel);
      if (EDV(pel,x,1) > q)
         EDV(pel,x,1) = q;
   }
}

#else

void take_positive_part(tGrid,x)
GRID *tGrid; INT x;
{  eprintf("Error: take_positive_part not available.\n");  }

#endif

void simple_opt(mg,x,p,t,type,use_bel)
MULTIGRID *mg;
INT x, p, t, type, use_bel;
{
   GRID *tGrid=TOP_GRID(mg);
   DOUBLE a=1.e-6, a0, v0, v1, max;
   INT k=0, l=0, ll=0, m=1;

   copy_ed_to_tau(tGrid,UU);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   v0 = residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,react,ft0,0.,use_bel);
   total_derivative_of_residual_based_error_estimator(TOP_GRID(mg),
                           A,1,U,D,R,Q,U,UU,p,TNU,bb0,bb1,react,ft0,0.,use_bel);
   inv(tGrid,p,p,t,type);
// max = max_factor_for_nonnegativity(tGrid,UU,p);
// a = MIN(a,max);
// if (a < 1.e-10) a = 1.e-8;
   if (a > 1.e-20){
   while (m){
      mult_and_add(tGrid,a,p,UU,x,t,type);
      take_positive_part(tGrid,x);
      copy_ed_to_tau(tGrid,x);
      iterate_SDFEM(mg,t0,t01,t02,ft0);
      v1 = residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,react,ft0,0.,use_bel);
      if (v1 < v0){
         v0 = v1;
         a0 = a;
         k = 1;
         if (l && ll)
            m = 0;
         else{
            a *= 2.;
/*
            if (a > max){
               a = max;
               ll = 1.; 
            }
*/
         }
      }
      else if (k)
         m = 0;
      else{
         a *= 0.5;
         l = 1;
      }
   }
   printf("error est = %e, a = %e\n\n",sqrt(v0),a0);
   mult_and_add(tGrid,a0,p,UU,UU,t,type);
   take_positive_part(tGrid,UU);
   }
}

void simple_opt2(mg,x,p,t,type,use_bel,fcn,grad_fcn)
MULTIGRID *mg;
INT x, p, t, type, use_bel;
DOUBLE (*fcn)(), (*grad_fcn)();
{
   GRID *tGrid=TOP_GRID(mg);
   DOUBLE a=1.e-6, a0, v0, v1, max;
   INT k=0, l=0, m=1;

   copy(tGrid,UU,x,t,type);
   v0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,p,t,type);
   inv(tGrid,p,p,t,type);
   if (a > 1.e-20){
   while (m){
      mult_and_add(tGrid,a,p,UU,x,t,type);
      take_positive_part(tGrid,x);
      v1 = fcn(mg,x,t,type);
      if (v1 < v0){
         v0 = v1;
         a0 = a;
         k = 1;
         if (l)
            m = 0;
         else{
            a *= 2.;
         }
      }
      else if (k)
         m = 0;
      else{
         a *= 0.5;
         l = 1;
      }
   }
   printf("error est = %e, a = %e\n\n",sqrt(v0),a0);
   mult_and_add(tGrid,a0,p,UU,UU,t,type);
   take_positive_part(tGrid,UU);
   }
}

DOUBLE simple_line_search(mg,tGrid,x,p,v0,new_v0,a,t,type,fcn)
MULTIGRID *mg;
GRID *tGrid;
INT x, p, t, type;
DOUBLE a, v0, *new_v0, (*fcn)();
{
   DOUBLE a0, v1, amin=1.e-12;
   INT i=0, k=0, l=0, m=1, max_it=1000;

   while (m){
      mult_and_add(tGrid,a,p,UU,x,t,type);
      take_positive_part(tGrid,x);
      v1 = fcn(mg,x,t,type);
      if (a < amin || i == max_it){
         if (v1 < v0){
            a0 = a;
            v0 = v1;
         }
         else if (!k)
            a0 = 0.;
         m = 0;
      }
      else if (v1 < v0){
         v0 = v1;
         a0 = a;
         k = 1;
         if (l)
            m = 0;
         else
            a *= 2.;
      }
      else if (k)
         m = 0;
      else{
         a *= 0.5;
         l = 1;
      }
      i++;
   }
   *new_v0 = v0;
/*
   if (USE_GRAD_IND == YES)
      printf("error est = %e, a = %e\n",v0,a0);
   else
      printf("error est = %e, a = %e\n",sqrt(v0),a0);
*/
   return(a0);
}

void vmult_e_2d_s(tGrid,p,q,r0,r1)
GRID *tGrid;
INT p, q;
DOUBLE r0, r1;
{
   ELEMENT *pel;

   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      EDV(pel,q,0) = r0*EDV(pel,p,0);
      EDV(pel,q,1) = r1*EDV(pel,p,1);
   }
}

DOUBLE multiple_line_search(mg,tGrid,x,p,q,v0,new_v0,a,a0_old,a1_old,t,type,fcn)
MULTIGRID *mg;
GRID *tGrid;
INT x, p, q, t, type;
DOUBLE a, *a0_old, *a1_old, v0, *new_v0, (*fcn)();
{
   DOUBLE phi0, phi1, a0=MAX(*a0_old,1.e-8), a1=MAX(*a1_old,1.e-8), 
          s0=0., s1=1., z;
   INT i=1;

   vmult_e_2d_s(tGrid,p,q,s0,1.-s0);
   a0 = *a0_old = simple_line_search(mg,tGrid,x,q,v0,&phi0,a0,t,type,fcn);
   while (fabs(s0-s1) > 0.2){
      vmult_e_2d_s(tGrid,p,q,s1,1.-s1);
      a1 = simple_line_search(mg,tGrid,x,q,v0,&phi1,a1,t,type,fcn);
      if (i){
         *a1_old = a1;
         i = 0;
      }
      if (phi0 > phi1){
         EXCHANGE(phi0,phi1,z)
         EXCHANGE(a0,a1,z)
         EXCHANGE(s0,s1,z)
      }
      s1 = 0.5*(s0+s1);
   }
   *new_v0 = phi0;
   z = MAX(s0,1.-s0);
   vmult_e_2d_s(tGrid,p,p,s0/z,(1.-s0)/z);
   return(a0*z);
}

void steepest_descent_ar(mg,tGrid,max_it,line_it,zoom_it,eps,x,p,y,z,v,w,t,type,
                  fcn,grad_fcn)
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE c1=1.e-4, c2=0.1, q, r, sum, phi0, dphi0, dphi0_old;
   DOUBLE a, a0, v0, v1, max;
   INT i=0, k1=1, k2=1, k, l, m;

   copy(tGrid,UU,x,t,type);
   phi0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,v,t,type);
   inv(tGrid,v,p,t,type);
   dphi0_old = -dot(tGrid,v,v,t,type);
   while (k2 && 
          max_abs_value(tGrid,v,t,type) > eps*(1.+fabs(phi0)) && i < max_it){
      a=1.e-6;
      copy(tGrid,UU,x,t,type);
      v0 = fcn(mg,x,t,type);
      grad_fcn(mg,x,v,t,type);
      inv(tGrid,v,p,t,type);
      dphi0 = dot(tGrid,v,p,t,type);
      a *= dphi0_old/dphi0;
      if (1 || !e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,
                         &a,c1,c2,x,p,y,z,t,type))
         a = simple_line_search(mg,tGrid,x,p,v0,&v0,a,t,type,fcn);
      mult_and_add(tGrid,a,p,UU,UU,t,type);
      take_positive_part(tGrid,UU);
      dphi0_old = dphi0;
      i++;
   }
}

void steepest_descent(mg,tGrid,max_it,line_it,zoom_it,eps,x,p,y,z,v,w,t,type,
                  fcn,grad_fcn)
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE c1=1.e-4, c2=0.99, q, r, sum, phi0, dphi0, dphi1, dphi0_old;
   DOUBLE a=1.e-6, a0, v0, v1, max;
   INT i=0, k, l, m;

   copy(tGrid,UU,x,t,type);
   phi0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,v,t,type);
   inv(tGrid,v,p,t,type);
   dphi0_old = -dot(tGrid,v,v,t,type);
   while ((q=max_abs_value(tGrid,v,t,type))>eps*(1.+fabs(phi0)) && i < max_it){
      printf("%i %e\n",i+1,q);
      copy(tGrid,UU,x,t,type);
      v0 = fcn(mg,x,t,type);
      grad_fcn(mg,x,v,t,type);
      inv(tGrid,v,p,t,type);
      dphi0 = dot(tGrid,v,p,t,type);
      a *= dphi0_old/dphi0;
//    if(!e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,
//                      &a,c1,c2,x,p,y,z,t,type)){

         if (i == 0)
            a = 1.e-6;
         a = simple_line_search(mg,tGrid,x,p,v0,&v0,a,t,type,fcn);
//    }
      mult_and_add(tGrid,a,p,UU,UU,t,type);
      take_positive_part(tGrid,UU);
      dphi0_old = dphi0;
      i++;

/*
copy(tGrid,UU,x,t,type);
v1 = fcn(mg,x,t,type);
grad_fcn(mg,x,v,t,type);
dphi1 = dot(tGrid,v,p,t,type);

if (v1 <= v0 + c1*a*dphi0 && fabs(dphi1) <= -c2*dphi0)
printf("O.K.\n");
else
printf("NOT NOT NOT %e %e    %e %e\n",v1,v0 + c1*a*dphi0,fabs(dphi1),-c2*dphi0);
*/
   }
}

void multiply_by_H(tGrid,p,g,mm,n,oldest,ind_s,ind_y,rho,t,type)
GRID *tGrid;
INT p, mm, n, ind_s, ind_y, t, type;
DOUBLE g, *rho;
{
   DOUBLE a[1000], b;
   INT i, ii;

   for (i = n-1; i >=0; i--){
      ii = (oldest + i) % mm;
      a[ii] = rho[ii]*dot(tGrid,ind_s+ii,p,t,type);
      mult_and_add(tGrid,-a[ii],ind_y+ii,p,p,t,type);
   }
   mult(tGrid,g,p,p,t,type);
   for (i = 0; i < n; i++){
      ii = (oldest + i) % mm;
      b = rho[ii]*dot(tGrid,ind_y+ii,p,t,type);
      mult_and_add(tGrid,a[ii]-b,ind_s+ii,p,p,t,type);
   }
}

void LBFGS(mg,tGrid,max_it,steep_it,line_it,zoom_it,eps,m,x,p,y,z,v,w,ind_s,
           t,type,fcn,grad_fcn)  /*  memory for up to ind_s+2*m-1 needed  */
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, m, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE c1=1.e-4, c2=0.9, q, r, phi0=1.+eps, dphi0, dphi0_old, a, g, 
          rho[1000], previous_fcn[10], fcn_var, a0=1.e-6, a1=1.e-6;
   INT i, ii, j, k, l, ll=0, ind_y=ind_s+m, n, oldest, itot=0, use_mls=0;

   for (i = 0; i < 10; i++)
      previous_fcn[i] = 1.e100;
   while (itot < max_it && phi0 > eps){
      printf("\nRestart\n");
      i = n = oldest = 0;
      l = 1;
      copy(tGrid,UU,x,t,type);
      phi0 = fcn(mg,x,t,type);
      grad_fcn(mg,x,v,t,type);
      g = 1.;
      while (l && i < max_it-itot && phi0 > eps){
         printf("\n%i\n",itot+i+1);
         inv(tGrid,v,p,t,type);
//       printf("|p| = %e, g = %e\n",sqrt(dot(tGrid,p,p,t,type)),g);
//norm_of_param(tGrid,p);
         if (!(itot == 0 && i < steep_it) && i > 0)
            multiply_by_H(tGrid,p,g,m,n,oldest,ind_s,ind_y,rho,t,type);
//       printf("|Hp| = %e\n",sqrt(dot(tGrid,p,p,t,type)));
//norm_of_param(tGrid,p);
         dphi0 = dot(tGrid,v,p,t,type);
         if (itot == 0 && i == 0)
            a = 1.e-6;
         else if (fabs(dphi0) > 1.e-30)
            a *= dphi0_old/dphi0;
         else
            a = 1.e-6;
         if (a < 1.e-12)
            a = 1.e-6;
         if (!(itot == 0 && i < steep_it) && i > 0)
            a = MIN(1.,a);
if (i==0) a = 1.e-6;
//       if(!e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,
//                         &a,c1,c2,x,p,y,z,t,type)) a = 0.;
         if (use_mls && COMPUTE_SC == YES)
            a = multiple_line_search(mg,tGrid,x,p,269,phi0,&phi0,a,&a0,&a1,
                                     t,type,fcn);
         else
            a = simple_line_search(mg,tGrid,x,p,phi0,&phi0,a,t,type,fcn);
//       if (!(itot == 0 && i < steep_it) && i > 0 && a < 1.e-12)
//          a = 0.;
         dphi0_old = dphi0;
         copy(tGrid,UU,x,t,type);
         mult_and_add(tGrid,a,p,UU,UU,t,type);
         take_positive_part(tGrid,UU);
//if (10*((itot+i+1)/10) == itot+i+1)
if (itot+i+1 == 80 || itot+i+1 == 100 || itot+i+1 == 120){
//mult1(tGrid,20.,UU);
//mult1s(tGrid,10.,UU,270);
take_positive_part(tGrid,UU);
}
//printf("UU:\n");
//norm_of_param(tGrid,UU);
if (itot+i+1 == -80){
supg_based_sold_tau(tGrid,UU);
take_positive_part(tGrid,UU);
solution_graph_for_gnuplot(tGrid,UU,1,"delta_sold_supg_based_graph.gnu",P0,VECTOR);
COMPUTE_SC=YES;
}
compute_scaling_par(TOP_GRID(mg),W);
phi0 = fcn(mg,UU,t,type);
printf("phi0 = %e, a = %e\n",phi0,a);
         fcn_var = fabs(previous_fcn[0] - phi0)/previous_fcn[0];
         if (!use_mls && fcn_var <= 1.e-2){
            fcn_var = use_mls = 1;
            if (SC_EXAMPLE == 55 || SC_EXAMPLE == 81)
               l = 0;
         }
if (l){
         if (!(itot == 0 && i < steep_it) && i > 0 && a < 1.e-6)
            l = 0;
         else if (fcn_var <= 1.e-4){
            if (ll)
               i = max_it;
            else{
               ll = 1;
               l = 0;
            }
         }
         else{
            ll = 0;
            for (j = 0; j < 9; j++)
               previous_fcn[j] = previous_fcn[j+1];
            previous_fcn[9] = phi0;
            if (itot == 0 && i < steep_it){
               copy(tGrid,UU,x,t,type);
               grad_fcn(mg,x,v,t,type);
            }
            else{
               if (n < m)
                  ii = n;
               else{
                  ii = oldest;
                  oldest++;
                  if (oldest == m)
                     oldest = 0;
               }
               subtr(tGrid,UU,x,ind_s+ii,t,type);
               copy(tGrid,UU,x,t,type);
               copy(tGrid,v,z,t,type);
               grad_fcn(mg,x,v,t,type);
               subtr(tGrid,v,z,ind_y+ii,t,type);
               r = dot(tGrid,ind_y+ii,ind_s+ii,t,type);
               if (r > 1.e-30)
                  rho[ii] = 1./r;
               else{
                  l = 0;
                  printf("y*s too small\n");
               }
               if (i == 0 || i >= m){
                  q = dot(tGrid,ind_y+ii,ind_y+ii,t,type);
                 if (q > 1.e-30)
                     g = r/q;
                  else{
                     l = 0;
                     printf("y*y too small\n");
                  }
               }
               if (n < m)
                  n++;
            }
         }
}
         i++;
         if (itot == 0 && i == steep_it)
            l = 0;
      }
      if (l)
         i = max_it;
      itot += i; 
   }
   grad_fcn(mg,UU,v,t,type);
      solution_graph_for_gnuplot(tGrid,v,0,"der_graph0.gnu",P0,VECTOR);
      solution_graph_for_gnuplot(tGrid,v,1,"der_graph1.gnu",P0,VECTOR);
      solution_graph_for_gnuplot(tGrid,ERR_VAR,0,"ind_graph0.gnu",P0,VECTOR);
      solution_graph_for_gnuplot(tGrid,ERR_VAR,1,"ind_graph1.gnu",P0,VECTOR);
}

INT LBFGSi_old(mg,tGrid,itot,max_it,line_it,zoom_it,eps,m,x,p,y,z,v,w,ind_s,
           t,type,fcn,grad_fcn,previous_fcn)
MULTIGRID *mg;                   /*  memory for up to ind_s+2*m-1 needed  */
GRID *tGrid;
INT itot, max_it, line_it, zoom_it, m, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)(), *previous_fcn;
{
   DOUBLE c1=1.e-4, c2=0.9, q, r, phi0, dphi0, dphi0_old, a=1.e-6, g, rho[1000];
   INT i=0, ii, j, k, l=1, ind_y=ind_s+m, n=0, oldest=0;

   copy(tGrid,UU,x,t,type);
   phi0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,v,t,type);
   inv(tGrid,v,p,t,type);
   while (l && i < max_it){
      printf("\n%i\n",itot+i+1);
printf("|p| = %e, g = %e\n",sqrt(dot(tGrid,p,p,t,type)),g);
      if (i > 0)
         multiply_by_H(tGrid,p,g,m,n,oldest,ind_s,ind_y,rho,t,type);
printf("|Hp| = %e\n",sqrt(dot(tGrid,p,p,t,type)));
      dphi0 = dot(tGrid,v,p,t,type);
      if (i > 0){
         if (fabs(dphi0) > 1.e-30)
            a *= dphi0_old/dphi0;
         a = MIN(1.,a);
      }
//    if(!e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,
//                      &a,c1,c2,x,p,y,z,t,type)) a = 0.;
      a = simple_line_search(mg,tGrid,x,p,phi0,&phi0,a,t,type,fcn);
      dphi0_old = dphi0;
      copy(tGrid,UU,x,t,type);
      mult_and_add(tGrid,a,p,UU,UU,t,type);
      take_positive_part(tGrid,UU);
      if (i > 0 && a < 1.e-6)
         l = 0;
      else if ((previous_fcn[0] - phi0)/previous_fcn[0] <= 1.e-4)
         i = max_it;
      else{
         for (j = 0; j < 9; j++)
            previous_fcn[j] = previous_fcn[j+1];
         previous_fcn[9] = phi0;
         if (n < m)
            ii = n;
         else{
            ii = oldest;
            oldest++;
            if (oldest == m)
               oldest = 0;
         }
         subtr(tGrid,UU,x,ind_s+ii,t,type);
         copy(tGrid,UU,x,t,type);
         copy(tGrid,v,z,t,type);
         grad_fcn(mg,x,v,t,type);
         subtr(tGrid,v,z,ind_y+ii,t,type);
         r = dot(tGrid,ind_y+ii,ind_s+ii,t,type);
         if (r > 1.e-30)
            rho[ii] = 1./r;
         else{
            l = 0;
            eprintf("Error: y*s too small.\n");
         }
         inv(tGrid,v,p,t,type);
         if (i == 0 || i >= m){
            q = dot(tGrid,ind_y+ii,ind_y+ii,t,type);
            if (q > 1.e-30)
               g = r/q;
            else{
               l = 0;
               eprintf("Error: y*y too small.\n");
            }
         }
         if (n < m)
            n++;
      }
      i++;
   }
   if (l)
      i = max_it;
   return(i);
}

/*
copy(tGrid,UU,x,t,type);
v1 = fcn(mg,x,t,type);
grad_fcn(mg,x,v,t,type);
dphi1 = dot(tGrid,v,p,t,type);

if (v1 <= v0 + c1*a*dphi0 && fabs(dphi1) <= -c2*dphi0)
printf("O.K.\n");
else
printf("NOT NOT NOT %e %e    %e %e\n",v1,v0 + c1*a*dphi0,fabs(dphi1),-c2*dphi0);
*/

void LBFGS_old(mg,tGrid,max_it,line_it,zoom_it,eps,m,x,p,y,z,v,w,ind_s,t,type,
           fcn,grad_fcn)  /*  memory for up to ind_s+2*m-1 needed  */
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, m, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE a=1.e-6, v0, previous_fcn[10];
   INT i, itot=0;

   for (i = 0; i < 10; i++)
      previous_fcn[i] = 1.e100;
   while (itot < max_it){
printf("Restart\n");
      i = LBFGSi_old(mg,tGrid,itot,max_it-itot,line_it,zoom_it,eps,
                 m,x,p,y,z,v,w,ind_s,t,type,fcn,grad_fcn,previous_fcn);
      if (i == 0){
         copy(tGrid,UU,x,t,type);
         v0 = fcn(mg,x,t,type);
         grad_fcn(mg,x,p,t,type);
         inv(tGrid,p,p,t,type);
         a = simple_line_search(mg,tGrid,x,p,v0,&v0,a,t,type,fcn);
         mult_and_add(tGrid,a,p,UU,UU,t,type);
         take_positive_part(tGrid,UU);
         printf("One step of steepest descent method performed.\n");
      }
      else
         itot += i;
   }
}

void LBFGS_ar(mg,tGrid,max_it,line_it,zoom_it,eps,mm,x,p,y,z,v,w,ind_s,t,type,
           fcn,grad_fcn)  /*  memory for up to ind_s+2*mm-1 needed  */
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, mm, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE c1=1.e-4, c2=0.99, q, r, sum, phi0, dphi0, dphi1, dphi0_old;
   DOUBLE a=1.e-6, a0, v0, v1, g, rho[1000];
   INT i=0, k, l, m, ind_y=ind_s+mm, n=0;

   copy(tGrid,UU,x,t,type);
   phi0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,v,t,type);
   inv(tGrid,v,p,t,type);
   dphi0_old = -dot(tGrid,v,v,t,type);
   while (max_abs_value(tGrid,v,t,type) > eps*(1.+fabs(phi0)) && i < max_it){
      if (i > 0)
         multiply_by_H(tGrid,p,v,g,mm,n,ind_s,ind_y,rho,t,type);
      dphi0 = dot(tGrid,v,p,t,type);
      a *= dphi0_old/dphi0;
      if(!e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,
                        &a,c1,c2,x,p,y,z,t,type)){
         if (i == 0)
            a = 1.e-6;
         copy(tGrid,UU,x,t,type);
         v0 = fcn(mg,x,t,type);
         a = simple_line_search(mg,tGrid,x,p,v0,&v0,a,t,type,fcn);
         phi0 = e_line_phi(mg,a,x,p,z,t,type);
      }
      dphi0_old = dphi0;
      copy(tGrid,UU,x,t,type);
      mult_and_add(tGrid,a,p,UU,UU,t,type);
      take_positive_part(tGrid,UU);
      subtr(tGrid,UU,x,ind_s+n,t,type);
      copy(tGrid,UU,x,t,type);
      copy(tGrid,v,z,t,type);
      grad_fcn(mg,x,v,t,type);
      subtr(tGrid,v,z,ind_y+n,t,type);
      rho[n] = 1./dot(tGrid,ind_y+n,ind_s+n,t,type);
      inv(tGrid,v,p,t,type);
      if (i == 0)
         g = 1./(rho[n]*dot(tGrid,ind_y+n,ind_y+n,t,type));
      i++;
      n++;
/*
copy(tGrid,UU,x,t,type);
v1 = fcn(mg,x,t,type);
grad_fcn(mg,x,v,t,type);
dphi1 = dot(tGrid,v,p,t,type);

if (v1 <= v0 + c1*a*dphi0 && fabs(dphi1) <= -c2*dphi0)
printf("O.K.\n");
else
printf("NOT NOT NOT %e %e    %e %e\n",v1,v0 + c1*a*dphi0,fabs(dphi1),-c2*dphi0);
*/
   }
}

void e_nonlinear_CG(mg,tGrid,max_it,line_it,zoom_it,eps,x,p,y,z,v,w,t,type,
                  fcn,grad_fcn)
MULTIGRID *mg;
GRID *tGrid;
INT max_it, line_it, zoom_it, x, p, y, z, v, w, t, type;
DOUBLE eps, (*fcn)(), (*grad_fcn)();
{
   DOUBLE a=1., c1=1.e-4, c2=0.1, q, r, sum, phi0, dphi0, dphi0_old;
   INT i=0, k1=1, k2=1;

// set_value(tGrid,0.,x,t,type);
   copy(tGrid,UU,x,t,type);
   phi0 = fcn(mg,x,t,type);
   grad_fcn(mg,x,v,t,type);
   copy(tGrid,v,w,t,type);
   inv(tGrid,v,p,t,type);
   r = dot(tGrid,v,v,t,type);
   dphi0_old = -r;
   while (k2 && 
          max_abs_value(tGrid,v,t,type) > eps*(1.+fabs(phi0)) && i < max_it){
      dphi0 = dot(tGrid,v,p,t,type);
      a *= dphi0_old/dphi0;
      if (e_line_search(mg,tGrid,line_it,zoom_it,phi0,dphi0,&phi0,&a,c1,c2,
                      x,p,y,z,t,type)){
         mult_and_add(tGrid,a,p,x,x,t,type);
         take_positive_part(tGrid,x);
         grad_fcn(mg,x,v,t,type);
         sum = dot(tGrid,v,v,t,type);
         q = dot(tGrid,v,w,t,type);
         if (fabs(q) > 0.1*r || q > sum)
            inv(tGrid,v,p,t,type);
         else
            mult_and_subtr(tGrid,(sum-q)/r,p,v,p,t,type);
         r = sum;
         dphi0_old = dphi0;
         copy(tGrid,v,w,t,type);
         k1 = 1;
      }
      else if (k1){
         k1 = 0;
         inv(tGrid,v,p,t,type);
         a = 1.;
      }
      else
         k2 = 0;
      i++;
   }
   copy(tGrid,x,UU,t,type);
   printf("%i iterations in nonlinear CG.\n",i);
}

#if 0
void compare_der_and_fd(pn,v0,v1,dx,d)
NODE *pn;
DOUBLE v0, v1, dx;
INT d;
{
   printf("err = %e, der = %+e, FD = %+e  ind = %i (%i)\n",
                 fabs(1.-(v1-v0)/dx/NDS(pn,d)), NDS(pn,d),(v1-v0)/dx,pn->index,
                 NTYPE(pn));
}

void ecompare_der_and_fd(pel,v0,v1,dx,d)
ELEMENT *pel;
DOUBLE v0, v1, dx;
INT d;
{
   DOUBLE xc[DIM], *x0, *x1, *x2, fd=(v1-v0)/dx, err;

/*
   coord_of_barycentre(pel,xc);
   printf("err = %e, der = %+e, FD = %+e (%5.3f,%5.3f)\n",
                 fabs(1.-(v1-v0)/dx/ED(pel,d)), ED(pel,d),(v1-v0)/dx,xc[0],xc[1]);
*/
   VERTICES_OF_ELEMENT(x0,x1,x2,pel);
   if (fabs(ED(pel,d))*1.e10 > fabs(fd))
      err = fabs(1.-fd/ED(pel,d));
   else
      err = 1.;
   printf(
     "err = %e, der = %+e, FD = %+e (%4.2f,%4.2f) (%4.2f,%4.2f) (%4.2f,%4.2f)\n",
     err, ED(pel,d),fd,x0[0],x0[1],x1[0],x1[1],x2[0],x2[1]);
}
#endif

void optimize_stab_par(mg,use_bel)
MULTIGRID *mg;
INT use_bel;
{
   GRID *tGrid=TOP_GRID(mg);
   ELEMENT *pel;
   NODE *pn;
   DOUBLE eps, v0, v1, dx, xc[DIM], bo[DIM], p, r, tau;
   INT i;

   mark_elements_for_error_ind(tGrid,bb0,bb1,Q);
   set_tau(tGrid,TNU,bb0,bb1,UU);
   mult(tGrid,10.,UU,W,0,PAR_TYPE);
for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ)
EDV(pel,W,1) = EDV(pel,CROSS_FCN_W,1) = 1.;
//read_stab_par(tGrid,"stab_par33_constr_a",UU);
   copy_ed_to_tau(tGrid,UU);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   if (SCALED_IND == YES){
      MAX_CROSS_DER = max_crosswind_derivative(tGrid,U,bb0,bb1);
      printf("max crosswind derivative: %e\n",MAX_CROSS_DER);
      p = CROSS_DER_WEIGHT;
      CROSS_DER_WEIGHT = 0.;
      printf("fcn = %e\n",
      r=res_and_der_based_error_indicator(tGrid,U,TNU,bb0,bb1,react,ft0));
      RES_NORM_WEIGHT = 1./r;
      CROSS_DER_WEIGHT = p;
   }
/*
   solution_graph_for_gnuplot(tGrid,U,0,"sol_graph.gnu",U_SPACE,U_STRUCTURE);
   sn_grid_data(tGrid,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph.grid_gnu");
   derivatives_of_residual_based_error_estimator(tGrid,
                                   U,D,Q,U,UU,TNU,bb0,bb1,react,ft0,0.,use_bel);
   v0 = residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                          react,ft0,0.,use_bel);
   if (USE_GRAD_IND == YES)
      printf("error est = %e\n",v0);
   else
      printf("error est = %e\n",sqrt(v0));
   total_derivative_of_residual_based_error_estimator(tGrid,
                           A,1,U,D,R,Q,U,UU,D,TNU,bb0,bb1,react,ft0,0.,use_bel);
   for (i=0; i < -100; i++)
      simple_opt(mg,D,Q,0,PAR_TYPE,use_bel);
*/
/*
   steepest_descent(mg,tGrid,10,1000,100,1.e-38,R,B,Q,P,X,Y,
                    0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
*/
   eps = 0.0001*residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                          react,ft0,0.,use_bel);
   LBFGS(mg,tGrid,5000,10,10,10,eps,100,U,F,D,Q,R,B,X,
         0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
/*
   printf("*** fcn = %e\n",e_fcn_for_minim(mg,UU,0,PAR_TYPE));
   LBFGS(mg,tGrid,500,0,10,10,1.e-38,100,U,F,D,Q,R,B,X,
         0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
   printf("*** fcn = %e\n",e_fcn_for_minim(mg,UU,0,PAR_TYPE));
   LBFGS(mg,tGrid,500,0,10,10,1.e-38,100,U,F,D,Q,R,B,X,
         0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
   printf("*** fcn = %e\n",e_fcn_for_minim(mg,UU,0,PAR_TYPE));
   LBFGS(mg,tGrid,500,0,10,10,1.e-38,100,U,F,D,Q,R,B,X,
         0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
*/
// e_nonlinear_CG(mg,tGrid,1000,1000,100,1.e-8,R,B,Q,P,X,Y,
//                0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
// save_stab_par(tGrid,"stab_par",UU);
/*
   if (SC_EXAMPLE == 81)
      sL_infinity_error(tGrid,U,R,&tau,t0,U_SPACE);
   if (PAR_TYPE == Q_SE)
      solution_graph_for_gnuplot(tGrid,UU,0,"delta_graph.gnu",P0,SCALAR);
   else if (PAR_TYPE == Q_VE){
      solution_graph_for_gnuplot(tGrid,UU,0,"delta_graph.gnu",P0,VECTOR);
      solution_graph_for_gnuplot(tGrid,UU,1,"delta_sold_graph.gnu",P0,VECTOR);
   }
   copy_ed_to_tau(tGrid,UU);
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   v0 = residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                          react,ft0,0.,use_bel);
   if (USE_GRAD_IND == YES)
      printf("error est = %e\n",v0);
   else
      printf("error est = %e\n",sqrt(v0));
*/
/*
   for (i = 1; i < 13; i++)
      if (i < 3 || i == 4 || i == 8 || i == 12){
         P_IN_LP_NORM = i;
         printf("L%2i est = %e\n",i,
            residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                         react,ft0,0.,use_bel));
      }
*/
   if (SC_EXAMPLE == 9 || SC_EXAMPLE == 91)
      save_profile(tGrid,0,0.,U,"profile_y.gnu",U_SPACE);
   solution_graph_for_gnuplot(tGrid,U,0,"sol_graph_opt.gnu",U_SPACE,U_STRUCTURE);
   sn_grid_data(tGrid,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph_opt.grid_gnu");
/*
   compute_res(tGrid,U,D,TNU,bb0,bb1,react,ft0);
   solution_graph_for_gnuplot(tGrid,D,0,"res_graph.gnu",P0,VECTOR);
   compute_norm_of_res(tGrid,D,L1_NORM,1.);
   compute_norm_of_res(tGrid,D,L2_NORM,1.);
   compute_norm_of_res(tGrid,D,LP_NORM,1.);
   compute_norm_of_res(tGrid,D,LP_NORM,2.);
   compute_norm_of_res(tGrid,D,LP_NORM,4.);
   compute_norm_of_res(tGrid,D,LP_NORM,8.);
   compute_norm_of_res(tGrid,D,LP_NORM,12.);
   compute_norm_of_res(tGrid,D,RES_FCN,1.);
   printf("jump_indicator: %e\n",jump_indicator(tGrid,U,R,Q,TNU));
   res_err_to_p1c_prolongation(tGrid,D,D,Q);
   sn_grid_data(tGrid,0.,1.,0.,1.,NVX-1,NVY-1,D,"res_graph.grid_gnu");
   bo[0] = -bb1(xc);
   bo[1] =  bb0(xc);
   xc[0] = -0.33*bo[0]/bo[1];
   xc[1] = 0.;
   draw_cut_for_gnuplot_p1c(tGrid,U,"cut.gnu",xc,bo);
*/
/*
   v0 = residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                          react,ft0,0.,use_bel);
   averaging_in_crosswind_direction_p1c(tGrid,U,R);
   subtr(tGrid,U,R,D,T_FOR_U,U_TYPE);
   set_zero_along_outflow_boundary(tGrid,D);
   solution_graph_for_gnuplot(tGrid,D,0,"sol_graph_diff.gnu",
                                                           U_SPACE,U_STRUCTURE);
   if(find_oscillation_regions_p1c(tGrid,0.01,D,CROSS_FCN_W) && v0 > eps){
      RESTRICT_FCN = NO;
      RESTRICT_DER = NO;
      COMPUTE_SC = YES;
      if (SCALED_IND == YES){
         MAX_CROSS_DER = max_crosswind_derivative(tGrid,U,bb0,bb1);
         RES_NORM_WEIGHT = 1.;
         CROSS_DER_WEIGHT = 0.;
         p = res_and_der_based_error_indicator(tGrid,U,TNU,bb0,bb1,react,ft0);
         RES_NORM_WEIGHT = 0.;
         CROSS_DER_WEIGHT = 1.;
         r = res_and_der_based_error_indicator(tGrid,U,TNU,bb0,bb1,react,ft0);
         RES_NORM_WEIGHT = 1./p;
         CROSS_DER_WEIGHT = 5./r;
      }
      else
         CROSS_DER_WEIGHT = 1.;
      eps = 0.0001*residual_based_error_estimator(tGrid,U,U,UU,TNU,bb0,bb1,
                                                          react,ft0,0.,use_bel);
      LBFGS(mg,tGrid,5000,10,10,10,eps,100,U,F,D,Q,R,B,X,
            0,PAR_TYPE,e_fcn_for_minim,e_grad_fcn_for_minim);
   }
   solution_graph_for_gnuplot(tGrid,U,0,"sol_graph_opt.gnu",U_SPACE,U_STRUCTURE);
   sn_grid_data(tGrid,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph_opt.grid_gnu");
*/
/*
   dx = 1.e-7;
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      tau = pel->tau;
      pel->tau += dx;
      iterate_SDFEM(mg,t0,t01,t02,ft0);
      v1 = residual_based_error_estimator(tGrid,
                                       U,U,UU,TNU,bb0,bb1,react,ft0,0.,use_bel);
      ecompare_der_and_fd(pel,v0,v1,dx,D);
      pel->tau = tau;
   }
   dx = 1.e-7;
   for (pn = FIRSTNODE(tGrid); pn; pn = pn->succ){
      NDS(pn,U) += dx;
      v1 = residual_based_error_estimator(tGrid,
                                       U,U,UU,TNU,bb0,bb1,react,ft0,0.,use_bel);
      compare_der_and_fd(pn,v0,v1,dx,D);
      NDS(pn,U) -= dx;
   }
*/
}
