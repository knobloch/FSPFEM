/******************************************************************************/
/*                                                                            */
/*                                   tests                                    */
/*                                                                            */
/******************************************************************************/

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
         GMRES(topGrid,U,F,3,0,10,R,15,1.e-8,1.e-8,NONZERO_INIT,
               mult_A,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
      else if (METHOD == DIAG_GMRES){
         copy_diag_to_vector(topGrid,A,R,T_FOR_U,U_TYPE,A_STRUCT);
         PGMRES(topGrid,U,F,D,Q,3,0,100000,B,15,1.e-8,1.e-8,NONZERO_INIT,
                mult_A,A,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
                R,0,0,0,DIAG_PR);
      }
      else if (METHOD == ILU_GMRES){
         make_ILU(topGrid,A,1,0,ILU_TYPE,ILU_STRUCT);
         PGMRES(topGrid,U,F,D,Q,3,0,500,R,15,1.e-8,1.e-8,NONZERO_INIT, //POZOR
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
      if ((U_SPACE == P1_MOD || U_SPACE == P1_NC)){
         if (U_SPACE == P1_MOD)
            vector_to_scalar_f(topGrid,U,U,0);
         nc_conf_averaging(topGrid,U,U,R,u0);
         solution_graph_for_gnuplot(topGrid,U,0,"n_sol_graph.gnu",P1C,
                                                                   U_STRUCTURE);
         error_graph_for_gnuplot(topGrid,u0,U,0,0.1,"n_err_graph.gnu",P1C,
                                                                   U_STRUCTURE);
      }
   }
   return(def);
}

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
}

void compute_solution_of_adjoint_problem(tGrid,Z,Z1,
                       u,der,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;                      /* g is the solution of the adjoint problem */
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT Z, Z1, u, der, g, v, w, tot_der, use_bel; 
{
   derivatives_of_residual_based_error_estimator(tGrid,u,der,g,v,w,
                                            eps,bb0,bb1,react,rhs,beta,use_bel);
   make_AT(tGrid,Z,Z1,U_TYPE,U_TYPE,A_STRUCT);
   if (METHOD == ILU_GMRES){
      make_ILU(tGrid,Z1,2,0,ILU_TYPE,ILU_STRUCT);
      PGMRES(tGrid,g,der,R,B,3,0,100000,W,15,1.e-8,1.e-8,0,
             mult_A,Z1,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,
             2,A_STRUCT,0,0,ILU_PR);
   }
   else if (METHOD == BASIC_GMRES){
      GMRES(tGrid,g,der,3,0,10,R,15,1.e-8,1.e-8,0,
            mult_A,Z1,F,T_FOR_U,U_TYPE,A_STRUCT,0,1,2,3,4,5,6,7,8,9,0.,0.,0.);
   }
   else if (METHOD == UMFPACK){
      fill_Ax(tGrid,Z1,Ap,Ai,Ax,U_TYPE,U_SPACE);
      make_vector_from_grid_data(tGrid,der,Rhs,U_TYPE,U_SPACE);
      solve_system_using_umfpack(Nj,Ap,Ai,Ax,Rhs,Sol);
      make_grid_data_from_vector(tGrid,g,Sol,U_TYPE,U_SPACE);
   }
   set_value(tGrid,0.,g,STOP_IS_FIRST_INNER,U_TYPE);
}

void total_derivative_of_residual_based_error_estimator(tGrid,Z,Z1,
                        u,der,d,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel)
GRID *tGrid;
FLOAT eps, (*bb0)(), (*bb1)(), (*react)(), (*rhs)(), beta;
INT Z, Z1, u, der, d, g, v, w, tot_der, use_bel; 
                /*  u ... solution at nodes, der ... derivatives at nodes  */
                /*  d, g ... auxiliary scalar node variable  */
                /*  v, w ... auxiliary scalar face variables  */
                /*  tot_der  */
{
   ELEMENT *pel;
   FLOAT b[DIM2][DIM2];

   compute_solution_of_adjoint_problem(tGrid,Z,Z1,
                       u,der,g,v,w,tot_der,eps,bb0,bb1,react,rhs,beta,use_bel);
                                 /* g is the solution of the adjoint problem */
   for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ){
      barycentric_coordinates(pel->n[0]->myvertex->x,pel->n[1]->myvertex->x,
                              pel->n[2]->myvertex->x,b);
      cd_res_times_sd(pel,tot_der,b,u,g,eps,bb0,bb1,react,rhs);
   }
}

DOUBLE e_fcn_for_minim(mg,x,t,type)
MULTIGRID *mg;
INT x, t, type;
{
   TAU_VARIABLE = x;
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   return(residual_based_error_estimator(TOP_GRID(mg),U,K,L,TNU,bb0,bb1,react,ft0,0.,USE_BEL));
}

void e_grad_fcn_for_minim(mg,x,y,t,type)
MULTIGRID *mg;
INT x, y, t, type;
{
   TAU_VARIABLE = x;
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   total_derivative_of_residual_based_error_estimator(TOP_GRID(mg),A,1,
                               U,P,G,J,K,L,y,TNU,bb0,bb1,react,ft0,0.,USE_BEL);
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
      mult_and_add(tGrid,a,p,N,x,t,type);
      take_positive_part(tGrid,x,W,t,type);
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
   return(a0);
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
   DOUBLE c1=1.e-4, c2=0.9, q, r, phi0, dphi0, dphi0_old, a, g, rho[1000],
          previous_fcn[10], fcn_var, a0=1.e-6, a1=1.e-6,
          total_time, seconds; // total_time a seconds jsou k mereni casu
   INT i, ii, j, k, l, ll=0, ind_y=ind_s+m, n, oldest, itot=0;
   clock_t start_time, stop_time; // k mereni casu
   FILE *fp;                      // k mereni casu

   fp = fopen("graph_of_minim.txt", "wb"); // k mereni casu
   start_time = clock(); // k mereni casu
   for (i = 0; i < 10; i++)
      previous_fcn[i] = 1.e100;
   while (itot < max_it){
      printf("\nRestart\n");
      i = n = oldest = 0;
      l = 1;
      copy(tGrid,N,x,t,type);
      phi0 = fcn(mg,x,t,type);
      grad_fcn(mg,x,v,t,type);
      g = 1.;
      while (l && i < max_it-itot){
         printf("\n%i\n",itot+i+1);
         inv(tGrid,v,p,t,type);
         if (!(itot == 0 && i < steep_it) && i > 0)
            multiply_by_H(tGrid,p,g,m,n,oldest,ind_s,ind_y,rho,t,type);
         dphi0 = dot(tGrid,v,p,t,type);
         printf("dphi0 = %e\n",dphi0);
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
         if (i==0)
            a = 1.e-6;
         a = simple_line_search(mg,tGrid,x,p,phi0,&phi0,a,t,type,fcn);
         dphi0_old = dphi0;
         copy(tGrid,N,x,t,type);
         mult_and_add(tGrid,a,p,N,N,t,type);
         take_positive_part(tGrid,N,W,t,type);
         printf("phi0 = %e, a = %e\n",phi0,a);
         fcn_var = fabs(previous_fcn[0] - phi0)/previous_fcn[0];
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
            // k mereni casu - zacatek
                 stop_time = clock();
                 total_time = (stop_time - start_time);
                 seconds = total_time/CLOCKS_PER_SEC;
                 fprintf(fp,"%e  %e \n",seconds,phi0);
            // k mereni casu - konec
            for (j = 0; j < 9; j++)
               previous_fcn[j] = previous_fcn[j+1];
            previous_fcn[9] = phi0;
            if (itot == 0 && i < steep_it){
               copy(tGrid,N,x,t,type);
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
               subtr(tGrid,N,x,ind_s+ii,t,type);
               copy(tGrid,N,x,t,type);
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
         i++;
         if (itot == 0 && i == steep_it)
            l = 0;
      }
      if (l){
         i = max_it;
         fclose(fp); // k mereni casu
      }
      itot += i; 
   }
}

void optimize_stab_par(mg,use_bel)
MULTIGRID *mg;
INT use_bel;
{
   GRID *tGrid=TOP_GRID(mg);
   DOUBLE v0, v1, p, r, variable;
   INT i;

   set_tau(tGrid,TNU,bb0,bb1,N);
   //sn_grid_data3(mg,0.,1.,0.,1.,100,100,N,"sol_graph_opt_P1_NC.grid_gnu");
   mult(tGrid,100.,N,W,0,Q_SE);
   TAU_VARIABLE = N;
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   LBFGS(mg,tGrid,20,10,10,10,1.e-38,100,16,17,18,19,20,21,22,
         0,Q_SE,e_fcn_for_minim,e_grad_fcn_for_minim);
   TAU_VARIABLE = N;
   iterate_SDFEM(mg,t0,t01,t02,ft0);
   
   variable = residual_based_error_estimator(tGrid,U,K,L,TNU,bb0,bb1,react,ft0,0.,USE_BEL);
   printf("%e\n",variable);
   
   /*
   derivatives_of_residual_based_error_estimator(tGrid,U,UU,J,K,L,
                                            TNU,bb0,bb1,react,ft0,0.,USE_BEL);
   solution_graph_for_gnuplot(tGrid,UU,0,"graph_der.gnu",U_SPACE,U_STRUCTURE);
   */


   /*  visualization of solution  */
   solution_graph_for_gnuplot(tGrid,U,0,"sol_graph_opt.gnu",U_SPACE,U_STRUCTURE);
   //save_profile(tGrid,0,0. ,U,"profile_y.gnu",U_SPACE);  // save profile
   //save_profile(tGrid,1,0.1,U,"profile_x.gnu",U_SPACE);  // save profile
   //sn_grid_data(mg,0.,1.,0.,1.,NVX-1,NVY-1,U,"sol_graph_opt.grid_gnu");

   /*  visualization of tau  */
   //solution_graph_for_gnuplot(tGrid,N,0,"sol_graph_opt_P1C.gnu",U_SPACE,U_STRUCTURE);
   //sn_grid_data(mg,0.,1.,0.,1.,600,600,N,"sol_graph_opt_P1C.grid_gnu");
   //sn_grid_data2(mg,0.,1.,0.,1.,600,600,N,"sol_graph_opt_P0.grid_gnu");
   //sn_grid_data3(mg,0.,1.,0.,1.,600,600,N,"sol_graph_opt_P1_NC.grid_gnu");
}
