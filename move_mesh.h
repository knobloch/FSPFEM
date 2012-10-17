/******************************************************************************/
/*                                                                            */
/*                    functions for moving triangulation                      */
/*                          (by Martin Sodomka)                               */
/******************************************************************************/

#if DIM == 2

#define QPI180		0.017453292519943295769236907684886
#define QPI2		1.5707963267948966192313216916398

/* moves nodes to adjust angle between wall and surface */
void surface_angle_old(mg, angle)
MULTIGRID *mg;
FLOAT angle;	/* in degrees */
{
	GRID *tGrid = TOP_GRID(mg);
	NODE *n, *nb;
	FACE *f;
	LINK *pli;
	NFLINK *pnf;
	FLOAT h, h0, d, dqx, x0, x1, y0, y1, z0, z1, p0, p1, q0, q1;
	FLOAT theta = angle*QPI180 - QPI2;

	printf("\nchanging angles..\n");
	for (n = FIRSTNODE(tGrid); n; n = n->succ)
		if (IS_ZNN(n) && IS_FBN(n) && IS_TOP_NODE(n)) {
			x0 = n->newvertex->x[0];
			x1 = n->newvertex->x[1];
			for (pli = TSTART(n); pli; pli = pli->next) {
				nb = NBNODE(pli);
				if (IS_ZNN(nb) && IS_TOP_NODE(nb)) {
					q0 = nb->newvertex->x[0];
					q1 = nb->newvertex->x[1];
				} else
				if (IS_FBN(nb) && IS_TOP_NODE(nb)) {
					z0 = nb->newvertex->x[0];
					z1 = nb->newvertex->x[1];
				} else 
				if (IS_TOP_NODE(nb)) {
					p0 = nb->newvertex->x[0];
					p1 = nb->newvertex->x[1];
				}
			}
			for (pnf = NFSTART(n); pnf; pnf = pnf->next)
				if (IS_CURVED_EDGE(f=NBFACE(pnf)) && IS_TOP_FACE(f))
					break;
			y0 = f->newvertex->x[0];
			y1 = f->newvertex->x[1];
			dqx = q0 - x0;
			if (fabs(dqx) > 1E-14) {
				eprintf("Neni kolma stena!");
			}
			h0 = z0 - x0;
			h = (y0 - (z0 + 3.*p0)/4.)/h0;
			if (h0 > 0) angle = theta;
				else	angle = -theta;
			d = 4.*y0 - 3.*x0 - z0;
			n->newvertex->x[1] = (4.*(y1 + h*x1) - z1 - d*tan(angle))/(4.*h + 3.);
			f->newvertex->x[1] = y1 + h*(x1 - n->newvertex->x[1]);
		}
}

/* moves nodes to adjust angle between wall and surface */
void surface_angle(mg, angle)
MULTIGRID *mg;
FLOAT angle;	/* in degrees */
{
	GRID *tGrid = TOP_GRID(mg);
	NODE *n, *nb;
	FACE *f;
	LINK *pli;
	NFLINK *pnf;
	FLOAT h, h0, d, dqx, x0, x1, y0, y1, z0, z1, p0, p1, q0, q1;
	FLOAT beta, theta = angle*QPI180 - QPI2;

	printf("changing angles..\n");
	for (n = FIRSTNODE(tGrid); n; n = n->succ)
		if (IS_ZNN(n) && IS_FBN(n) && IS_TOP_NODE(n)) {
			x0 = n->newvertex->x[0];
			x1 = n->newvertex->x[1];
			for (pli = TSTART(n); pli; pli = pli->next) {
				nb = NBNODE(pli);
				if (IS_ZNN(nb) && IS_TOP_NODE(nb)) {
					q0 = nb->newvertex->x[0];
					q1 = nb->newvertex->x[1];
				} else
				if (IS_FBN(nb) && IS_TOP_NODE(nb)) {
					z0 = nb->newvertex->x[0];
					z1 = nb->newvertex->x[1];
				} else 
				if (IS_TOP_NODE(nb)) {
					p0 = nb->newvertex->x[0];
					p1 = nb->newvertex->x[1];
				}
			}
			for (pnf = NFSTART(n); pnf; pnf = pnf->next)
				if (IS_CURVED_EDGE(f=NBFACE(pnf)) && IS_TOP_FACE(f))
					break;
			y0 = f->newvertex->x[0];
			y1 = f->newvertex->x[1];
			dqx = q0 - x0;
			beta = 0.;
			if (fabs(dqx) > 1E-15) {
				//eprintf("Neni kolma stena!");
				if (fabs(q1-x1) > 1E-15) beta = atan(dqx/(q1 - x1));
				else beta = (dqx > 0 ? -QPI2 : QPI2);
			}
			h0 = z0 - x0;
			if (h0 > 0) {
				angle = theta - beta;
				h = (y0 - (z0 + 3.*q0)/4.)/h0;
			} else {
				angle = -theta - beta;
				h = (y0 - (z0 + 3.*q0)/4.)/h0;
			}
			d = 4.*y0 - 3.*x0 - z0;
			n->newvertex->x[1] = (4.*(y1 + h*x1) - z1 - d*tan(angle))/(4.*h + 3.);
			f->newvertex->x[1] = y1 + h*(x1 - n->newvertex->x[1]);
		}
}

FLOAT c_volume_newvertex(f0,n0,n1,n2)
FACE *f0;
NODE *n0, *n1, *n2;
{
	FLOAT a, b, c, d, e, g;

	a = n1->newvertex->x[0] - n0->newvertex->x[0];
	e = n1->newvertex->x[1] - n0->newvertex->x[1];
	g = n2->newvertex->x[0] - n0->newvertex->x[0];
	c = n2->newvertex->x[1] - n0->newvertex->x[1];
	b = n2->newvertex->x[1]*(f0->newvertex->x[0] - n1->newvertex->x[0]) 
		- n1->newvertex->x[1]*(f0->newvertex->x[0] - n2->newvertex->x[0]);
	d = f0->newvertex->x[1]*(a-g);
	return ( 0.5*fabs(a*c - e*g + 4.*(b + d)/3.) );
}

void print_volumes(mg)
MULTIGRID *mg;
{
	ELEMENT *pel;
	NODE *nx, *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	FLOAT vol;

	printf("\n");
	for (nx = FIRSTNODE(TOP_GRID(mg)); nx; nx = SUCC(nx)) {
		if (IS_ZNN(nx) && IS_FBN(nx) && IS_TOP_NODE(nx)) {
			vol = 0.;
			for (pel = FIRSTELEMENT(TOP_GRID(mg)); pel; pel = SUCC(pel)) 
				if (IS_TOP_ELEMENT(pel)) {
					NODES_OF_ELEMENT(n0,n1,n2,pel);
					FACES_OF_ELEMENT(f0,f1,f2,pel);
					if (n0 == nx || n1 == nx || n2 == nx)
						if (IS_CURVED_EDGE(f0))
							vol += c_volume_newvertex(f0,n0,n1,n2);
						else if (IS_CURVED_EDGE(f1))
							vol += c_volume_newvertex(f1,n1,n0,n2);
						else if (IS_CURVED_EDGE(f2))
							vol += c_volume_newvertex(f2,n2,n0,n1);
						else
							vol += volume(pel->n[0]->newvertex->x,pel->n[1]->newvertex->x,
																  pel->n[2]->newvertex->x);
				}
			printf("volume = %g\n", vol);
		}
	}
}

void print_angles(mg)
MULTIGRID *mg;
{
	NODE *nx, *nz;
	FACE *fy;
	LINK *pl;
	NFLINK *pnf;
	FLOAT a, d;

	for (nx = FIRSTNODE(TOP_GRID(mg)); nx; nx = SUCC(nx)) {
		if (IS_ZNN(nx) && IS_FBN(nx) && IS_TOP_NODE(nx)) {
			for (pl = START(nx); NOT_FBN(NBNODE(pl)); pl = NEXT(pl));
			nz = NBNODE(pl);
			for (pnf = NFSTART(nx); NOT_CURVED_EDGE(NBFACE(pnf)); pnf = NEXT(pnf));
			fy = NBFACE(pnf);
			d = (4*fy->newvertex->x[1] - 3*nx->newvertex->x[1] - nz->newvertex->x[1])/
				(4*fy->newvertex->x[0] - 3*nx->newvertex->x[0] - nz->newvertex->x[0]);
			a = atan(d);
			if (nx->newvertex->x[0] < nz->newvertex->x[0])
				a += QPI2;
			else
				a = QPI2 - a;
			a /= QPI180;
			/*printf("\nnz = [%g, %g]; fy = [%g, %g]\n", nz->myvertex->x[0], nz->myvertex->x[1],
													 fy->myvertex->x[0], fy->myvertex->x[1]);*/
			printf("angle in node [%g, %g] is: %g (DEG)\n", nx->newvertex->x[0],
														    nx->newvertex->x[1], a);
		}
	}
	printf("\n");
}

/*----------------------------------------------------------------*
void solve_lap_t(mg,tGrid,Z,g,to,tn)
MULTIGRID *mg;
GRID *tGrid;
INT Z;
FLOAT (*g)(), to, tn;
{ 
	sbound_val1_t(mg,Z,g,to,tn);
	sset_value(tGrid,0.,F,T_FOR_U);
	sset_mat_value(tGrid,A,0.);
	add_Laplace_matr(tGrid,1.,A,P1C,A_STRUCT,SCALAR,NO);
	subtract_Acontribution_from_bc(tGrid,A,Z,F,Q,D,R,T_FOR_U,Q_SN,A_STRUCT);
	sset_value(tGrid,0.,Z,T_FOR_U);
	sILU(tGrid,A,A1,T_FOR_U);
	PCG(tGrid,F,Z,D,Q,R,0,500,1.e50,EPS_PCGMOV,0,mult_A,A,
		D,T_FOR_U,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,A1,0,0,0,ILU_PR);
}

FLOAT hr_fce(x, time)
FLOAT x[2], time;
{
	if (fabs(x[0]) < EPS || fabs(x[0]-1.) < EPS || fabs(x[1]) < EPS)
		return 0.;
	else
		return ( 0.1*sin(2*M_PI*x[0])*sin(time) );
}

void move_new(mg, t_old, t_new)
MULTIGRID *mg;
FLOAT t_old, t_new;
{
	GRID *tGrid = TOP_GRID(mg);
	NODE *n;
	FACE *f;

//printf("\nmoving ...\n");
	solve_lap_t(mg,tGrid,UU,hr_fce,t_old,t_new);
	for (n = tGrid->firstN; n; n = n->succ)	{
		n->newvertex->x[1] += NDS(n,UU);
	}
//printf("nodes done\n");
	for (f = tGrid->firstF; f != NULL; f = f->succ) {
		if (IS_CURVED_EDGE(f))
			f->newvertex->x[1] += hr_fce(f->newvertex->x, t_new) - hr_fce(f->newvertex->x, t_old);
	}
//printf("faces done\n\n");
}

/*----------------------------------------------------------------*/
INT elem_ok(x0, x1, x2)
FLOAT *x0, *x1, *x2;
{
   FLOAT det;
  
   det = (x0[0]-x2[0])*(x1[1]-x2[1]) - (x1[0]-x2[0])*(x0[1]-x2[1]);
  
   if ( ( (det>0)?det:-det ) < 1.e-7 ){
      return (0);
   }
   return (1);
}

INT iso_elem_ok(x0, x1, x2, x12)
FLOAT *x0, *x1, *x2, *x12;
{
   FLOAT s0, s1, s2, a[2][2], alpha[2], jac[3];

   alpha[0] = 4.*(x12[0] - (x1[0] + x2[0])/2.);
   alpha[1] = 4.*(x12[1] - (x1[1] + x2[1])/2.);

   a[0][0] = x1[0] - x0[0];
   a[0][1] = x2[0] - x0[0];
   a[1][0] = x1[1] - x0[1];
   a[1][1] = x2[1] - x0[1];

   jac[0] = a[0][0]*alpha[1] - a[1][0]*alpha[0];
   jac[1] = a[1][1]*alpha[0] - a[0][1]*alpha[1];
   jac[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
  
   s0 = jac[2];
   s1 = jac[0]+jac[2];
   s2 = jac[1]+jac[2];
   if ( fabs(s0) < 1.e-7 || fabs(s1) < 1.e-7 || fabs(s2) < 1.e-7 ||
        s0*s1 < 1.e-14 || s0*s2 < 1.e-14 || s1*s2 < 1.e-14)    
   {
      return (0);
   }
   return (1);
}

INT elements_ok(tGrid)
GRID *tGrid;
{
	ELEMENT *pel;

	for (pel = FIRSTELEMENT(tGrid); pel; pel = SUCC(pel)) {
		if (CFACE(pel->f[0])) {
			if ( !iso_elem_ok(pel->n[0]->newvertex->x, pel->n[1]->newvertex->x,
							  pel->n[2]->newvertex->x, pel->f[0]->newvertex->x) )
				return (0);
		} else
		if (CFACE(pel->f[1])) {
			if ( !iso_elem_ok(pel->n[1]->newvertex->x, pel->n[0]->newvertex->x,
							  pel->n[2]->newvertex->x, pel->f[1]->newvertex->x) )
				return (0);
		} else
		if (CFACE(pel->f[2])) {
			if ( !iso_elem_ok(pel->n[2]->newvertex->x, pel->n[1]->newvertex->x,
							  pel->n[0]->newvertex->x, pel->f[2]->newvertex->x) )
				return (0);
		} else {
			if ( !elem_ok(pel->n[0]->newvertex->x, pel->n[1]->newvertex->x,
												   pel->n[2]->newvertex->x) )
				return (0);
		}
	}
	return (1);
}

void solve_lap_nd(mg,tGrid,Z,n,i,dt,t)
MULTIGRID *mg;
GRID *tGrid;
INT Z, n, i, t;
FLOAT dt;
{
	sbound_val1_nd(mg,Z,n,i,dt,t);
	sset_value(tGrid,0.,W,GO_THROUGH_ALL);
	sset_mat_value(tGrid,A,0.);
	add_Laplace_matr(tGrid,1.,A,P1C,A_STRUCT,SCALAR,NO);
	sset_value(tGrid,0.,Z,t);
	mso_smultA(tGrid,A,Z,W,t | T_FOR_BC);
	sinv(tGrid,W,W,t);
	sILU(tGrid,A,A1,t);
	PCG(tGrid,W,Z,D,Q,R,0,500,1.e50,EPS_PCGMOV,0,mult_A,A,
		D,t,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,A1,0,0,0,ILU_PR);
}

/*	moves triang. according to normal component
	of given velocity u_n at time t_n
	.. and returns time step
*/
FLOAT move_normdirect(mg, tau)	
MULTIGRID *mg;
FLOAT tau;	// = t_n+1 - t_n
{
	GRID *pg, *tGrid = TOP_GRID(mg);
	ELEMENT *pel;
	NODE *n, *n1, *n2;
	FACE *f, *f1, *f2;

BEGINNING:
	for (f = FIRSTF(tGrid); f != NULL; f = f->succ) 
		if (IS_CURVED_EDGE(f) && IS_TOP_FACE(f)) {
			f->newvertex->x[0] = f->myvertex->x[0] + tau*FDV(f,N,0);
			f->newvertex->x[1] = f->myvertex->x[1] + tau*FDV(f,N,1);
		}

	solve_lap_nd(mg,tGrid,X,N,0,tau,T_FOR_LAP);
	solve_lap_nd(mg,tGrid,Y,N,1,tau,T_FOR_LAP);
	for (n = FIRSTN(tGrid); n != NULL; n = n->succ)	{
		n->newvertex->x[0] = n->myvertex->x[0] + NDS(n,X);
		n->newvertex->x[1] = n->myvertex->x[1] + NDS(n,Y);
	}

#if STEP_CONTROL == YES
//	for (pg = tGrid; pg->coarser != NULL; pg = pg->coarser)
		if (!elements_ok(tGrid)) {
			if (tau > TAU_MIN) {
				tau /= 2.;
				printf("\nrefining time step.. tau = %e\n\n", tau);				
				goto BEGINNING;
			} else {
				eprintf("ERROR: Time step < TAU_MIN and still need refinement!\n");
			}
		}
#endif

#if CFACE_ON_COARSE == YES
	for (pg = tGrid->coarser; pg != NULL; pg = pg->coarser)
		for (f = FIRSTF(pg); f != NULL; f = f->succ)
			if (IS_CURVED_EDGE(f)) {
				f->newvertex->x[0] = f->myvertex->topnode->newvertex->x[0];
				f->newvertex->x[1] = f->myvertex->topnode->newvertex->x[1];
			}
#else
	for (pg = tGrid->coarser; pg != NULL; pg = pg->coarser)
		for (pel = FIRSTELEMENT(pg); pel; pel = pel->succ) {
			NODES_OF_ELEMENT(n,n1,n2,pel);
			FACES_OF_ELEMENT(f,f1,f2,pel);
			if (IS_CURVED_EDGE(f)) {
				f->newvertex->x[0] = (n1->newvertex->x[0] + n2->newvertex->x[0])/2.;
				f->newvertex->x[1] = (n1->newvertex->x[1] + n2->newvertex->x[1])/2.;
				f->myvertex->x[0] = (n1->myvertex->x[0] + n2->myvertex->x[0])/2.;
				f->myvertex->x[1] = (n1->myvertex->x[1] + n2->myvertex->x[1])/2.;
			} else
			if (IS_CURVED_EDGE(f1)) {
				f1->newvertex->x[0] = (n->newvertex->x[0] + n2->newvertex->x[0])/2.;
				f1->newvertex->x[1] = (n->newvertex->x[1] + n2->newvertex->x[1])/2.;
				f1->myvertex->x[0] = (n->myvertex->x[0] + n2->myvertex->x[0])/2.;
				f1->myvertex->x[1] = (n->myvertex->x[1] + n2->myvertex->x[1])/2.;
			} else
			if (IS_CURVED_EDGE(f2)) {
				f2->newvertex->x[0] = (n1->newvertex->x[0] + n->newvertex->x[0])/2.;
				f2->newvertex->x[1] = (n1->newvertex->x[1] + n->newvertex->x[1])/2.;
				f2->myvertex->x[0] = (n1->myvertex->x[0] + n->myvertex->x[0])/2.;
				f2->myvertex->x[1] = (n1->myvertex->x[1] + n->myvertex->x[1])/2.;
			}
		}
#endif

#if USE_WS_ANGLE == YES
	print_volumes(mg);
	print_angles(mg);
	surface_angle(mg, WS_ANGLE);
	print_volumes(mg);
	print_angles(mg);
#endif

	return ( tau );
}

/*	calculates normal vectors in boundary faces
	and bnd. nodes (as average value)
*/
void calc_normals(tGrid)
GRID *tGrid;
{
	ELEMENT *pel;
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	FLOAT vn, v0;
	NFLINK *nfl;
	int i;

	for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ) 
		if (IS_TOP_ELEMENT(pel)) {
			NODES_OF_ELEMENT(n0,n1,n2,pel);
			FACES_OF_ELEMENT(f0,f1,f2,pel);
			if (IS_CURVED_EDGE(f0) || IS_ZNF(f0)) {
				normal_vector(n1->myvertex->x, n2->myvertex->x, f0, FDVP(f0,N));
			} else
			if (IS_CURVED_EDGE(f1) || IS_ZNF(f1)) {
				normal_vector(n0->myvertex->x, n2->myvertex->x, f1, FDVP(f1,N));
			} else
			if (IS_CURVED_EDGE(f2) || IS_ZNF(f2)) {
				normal_vector(n0->myvertex->x, n1->myvertex->x, f2, FDVP(f2,N));
			}
		}

	for (n0 = FIRSTN(tGrid); n0 != NULL; n0 = n0->succ) {
		if ( IS_TOP_NODE(n0) && IS_FBN(n0) && NOT_ZNN(n0) ) {
			i = 0;
			ND(n0,N,0) = ND(n0,N,1) = 0.;
			for (nfl = TNFSTART(n0); nfl != NULL; nfl = nfl->next) {
				if (IS_CURVED_EDGE(nfl->nbface)) {
					i++;
					ND(n0,N,0) += FDV(nfl->nbface,N,0);
					ND(n0,N,1) += FDV(nfl->nbface,N,1);
				}
			}
			if (i > 1) {
				vn = sqrt(ND(n0,N,0)*ND(n0,N,0) + ND(n0,N,1)*ND(n0,N,1));
				if (fabs(vn) < 1e-15) eprintf("Zero normal in node !!\n");
				ND(n0,N,0) /= vn;
				ND(n0,N,1) /= vn;
			} else {
				ND(n0,N,0) = 0.;
				if (i == 1) ND(n0,N,1) = 1.;
				else ND(n0,N,1) = 0.;
			}
		}
		if ( IS_TOP_NODE(n0) && IS_ZNN(n0) ) {
			ND(n0,N,0) = ND(n0,N,1) = 0.;
			for (nfl = TNFSTART(n0); nfl != NULL; nfl = nfl->next)
				if (IS_ZNF(nfl->nbface)) {
					ND(n0,N,0) += FDV(nfl->nbface,N,0);
					ND(n0,N,1) += FDV(nfl->nbface,N,1);
				}
			/* normovani + otoceni doprava o Pi/2 na tecnu: */
			vn = sqrt(ND(n0,N,0)*ND(n0,N,0) + ND(n0,N,1)*ND(n0,N,1));
			if (fabs(vn) < 1e-15) eprintf("Zero normal in node !!\n");
			v0 = ND(n0,N,0) / vn;
			ND(n0,N,0) = ND(n0,N,1) / vn;
			ND(n0,N,1) = -v0;
		}
	}
}

/*	sets normal velocity to boundary nodes,
	according to starting velocity 'u0'
*/
void set_norm_u_t0(tGrid)
GRID *tGrid;
{
	NODE *n0;
	FACE *f0;
	FLOAT dp;

	for (n0 = FDBN(tGrid); n0 != NULL; n0 = n0->succ)
		if (IS_TOP_NODE(n0)) {
			if (IS_FBN(n0) || IS_ZNN(n0)) {
				dp = ND(n0,N,0)*u01(n0->myvertex->x) + ND(n0,N,1)*u02(n0->myvertex->x);
				ND(n0,N,0) *= dp;
				ND(n0,N,1) *= dp;
			}
			if (NOT_FN(n0)) {
				ND(n0,N,0) = ND(n0,N,1) = 0.;
			}
		}

	for (f0 = FDBF(tGrid); f0 != NULL; f0 = f0->succ)
		if (IS_TOP_FACE(f0) && IS_CURVED_EDGE(f0)) {
			dp = FDV(f0,N,0)*u01(f0->myvertex->x)+FDV(f0,N,1)*u02(f0->myvertex->x);
			FDV(f0,N,0) *= dp;
			FDV(f0,N,1) *= dp;
		}
	/* v FDV uz je samotna normalova rychlost, nikoli jeji koeficient */
}

/*	sets normal velocity to boundary nodes
*/
void set_norm_u(tGrid, Z)
GRID *tGrid;
int Z;
{
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	ELEMENT *pel;
	FLOAT dp;

	for (n0 = FDBN(tGrid); n0 != NULL; n0 = n0->succ)
		if (IS_TOP_NODE(n0)) {
			if (IS_FBN(n0) || IS_ZNN(n0)) {
				dp = ND(n0,N,0)*ND(n0,Z,0) + ND(n0,N,1)*ND(n0,Z,1);
				ND(n0,N,0) *= dp;
				ND(n0,N,1) *= dp;
			}
			if (NOT_FN(n0)) {
				ND(n0,N,0) = ND(n0,N,1) = 0.;
			}
		}

	for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
		if (IS_TOP_ELEMENT(pel)) {
			NODES_OF_ELEMENT(n0, n1, n2, pel);
			FACES_OF_ELEMENT(f0, f1, f2, pel);
			if (IS_CURVED_EDGE(f0)) {
				dp = (FDV(f0,N,0)*( FDV(f0,Z,0)/2. + ND(n1,Z,0) + ND(n2,Z,0) ) +
					  FDV(f0,N,1)*( FDV(f0,Z,1)/2. + ND(n1,Z,1) + ND(n2,Z,1) ))/2.;
				FDV(f0,N,0) *= dp;
				FDV(f0,N,1) *= dp;
			} else
			if (IS_CURVED_EDGE(f1)) {
				dp = (FDV(f1,N,0)*( FDV(f1,Z,0)/2. + ND(n0,Z,0) + ND(n2,Z,0) ) +
					  FDV(f1,N,1)*( FDV(f1,Z,1)/2. + ND(n0,Z,1) + ND(n2,Z,1) ))/2.;
				FDV(f1,N,0) *= dp;
				FDV(f1,N,1) *= dp;
			} else
			if (IS_CURVED_EDGE(f2)) {
				dp = (FDV(f2,N,0)*( FDV(f2,Z,0)/2. + ND(n1,Z,0) + ND(n0,Z,0) ) +
					  FDV(f2,N,1)*( FDV(f2,Z,1)/2. + ND(n1,Z,1) + ND(n0,Z,1) ))/2.;
				FDV(f2,N,0) *= dp;
				FDV(f2,N,1) *= dp;
			}
		}
	/* v FDV uz je samotna normalova rychlost, nikoli jeji koeficient */
}

/*	sets modified velocity to boundary nodes
*/
void set_modif_u(tGrid, Z)
GRID *tGrid;
int Z;
{
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	ELEMENT *pel;

	for (n0 = FDBN(tGrid); n0 != NULL; n0 = n0->succ)
		if (IS_TOP_NODE(n0)) {
			if (IS_FBN(n0) || IS_ZNN(n0)) {
				ND(n0,N,0) = ND(n0,Z,0);
				ND(n0,N,1) = ND(n0,Z,1);
			}
			if (NOT_FN(n0)) {
				ND(n0,N,0) = ND(n0,N,1) = 0.;
			}
		}

	for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ)
		if (IS_TOP_ELEMENT(pel)) {
			NODES_OF_ELEMENT(n0, n1, n2, pel);
			FACES_OF_ELEMENT(f0, f1, f2, pel);
			if (IS_CURVED_EDGE(f0)) {
				FDV(f0,N,0) = FDV(f0,Z,0)/4. + (ND(n1,Z,0) + ND(n2,Z,0))/2.;
				FDV(f0,N,1) = FDV(f0,Z,1)/4. + (ND(n1,Z,1) + ND(n2,Z,1))/2.;
			} else
			if (IS_CURVED_EDGE(f1)) {
				FDV(f1,N,0) = FDV(f1,Z,0)/4. + (ND(n0,Z,0) + ND(n2,Z,0))/2.;
				FDV(f1,N,1) = FDV(f1,Z,1)/4. + (ND(n0,Z,1) + ND(n2,Z,1))/2.;
			} else
			if (IS_CURVED_EDGE(f2)) {
				FDV(f2,N,0) = FDV(f2,Z,0)/4. + (ND(n0,Z,0) + ND(n1,Z,0))/2.;
				FDV(f2,N,1) = FDV(f2,Z,1)/4. + (ND(n0,Z,1) + ND(n1,Z,1))/2.;
			}
		}
	/* v FDV uz je samotna rychlost, nikoli jeji koeficient */
}

/*	sets initial velocity 'u0' to all nodes and faces
*/
void set_u0(tGrid, Z)
GRID *tGrid;
INT Z;
{
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;
	ELEMENT *pel;
	FLOAT xm[DIM];

	for (n0 = FIRSTN(tGrid); n0 != NULL; n0 = n0->succ) {
		ND(n0,Z,0) = u01(n0->myvertex->x);
		ND(n0,Z,1) = u02(n0->myvertex->x);
	}

	for (pel = FIRSTELEMENT(tGrid); pel != NULL; pel = pel->succ) {
		NODES_OF_ELEMENT(n0,n1,n2,pel);
		FACES_OF_ELEMENT(f0,f1,f2,pel);
		if (IS_CURVED_EDGE(f0)) {
			xm[0] = f0->myvertex->x[0];
			xm[1] = f0->myvertex->x[1];
		} else
			AVERAGE(n1->myvertex->x,n2->myvertex->x,xm);
		FDV(f0,Z,0) = 4.*u01(xm) - 2.*(ND(n1,Z,0)+ND(n2,Z,0));
		FDV(f0,Z,1) = 4.*u02(xm) - 2.*(ND(n1,Z,1)+ND(n2,Z,1));

		if (IS_CURVED_EDGE(f1)) {
			xm[0] = f1->myvertex->x[0];
			xm[1] = f1->myvertex->x[1];
		} else
			AVERAGE(n0->myvertex->x,n2->myvertex->x,xm);
		FDV(f1,Z,0) = 4.*u01(xm) - 2.*(ND(n0,Z,0)+ND(n2,Z,0));
		FDV(f1,Z,1) = 4.*u02(xm) - 2.*(ND(n0,Z,1)+ND(n2,Z,1));

		if (IS_CURVED_EDGE(f2)) {
			xm[0] = f2->myvertex->x[0];
			xm[1] = f2->myvertex->x[1];
		} else 
			AVERAGE(n0->myvertex->x,n1->myvertex->x,xm);
		FDV(f2,Z,0) = 4.*u01(xm) - 2.*(ND(n0,Z,0)+ND(n1,Z,0));
		FDV(f2,Z,1) = 4.*u02(xm) - 2.*(ND(n0,Z,1)+ND(n1,Z,1));
	}
}

/*
void solve_Lap_system(topGrid,tol,Z)
GRID *topGrid;
FLOAT tol;
INT Z;
{
//printf("Solving Laplace system...\n");
	sILU(topGrid,A,A1,T_FOR_U);

	PCG(topGrid,F,Z,D,Q,R,0,500,1.e50,tol,0,mult_A,A,
		D,T_FOR_U,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,A1,0,0,0,ILU_PR);
	//defect(topGrid,A,F,Z,D,Q,T_FOR_U,Q_SN,0);
//printf("defect: %e\n\n", sqrt(dot(topGrid,D,D,T_FOR_U,Q_SN)));
}

#define SOLVE_LAP(Z,g){	\
	sbound_val1(mg,Z,g);\
	sset_value(topGrid,0.,F,T_FOR_U);\
	sset_mat_value(topGrid,A,0.);\
	add_Laplace_matr(topGrid,1.,A,P1C,A_STRUCT,SCALAR,NO);\
	subtract_Acontribution_from_bc(topGrid,A,Z,F,Q,D,R,T_FOR_U,Q_SN,A_STRUCT);\
	sset_value(topGrid,0.,Z,T_FOR_U);\
	sILU(topGrid,A,A1,T_FOR_U);\
	PCG(topGrid,F,Z,D,Q,R,0,500,1.e50,EPS_PCGMOV,0,mult_A,A,	\
		D,T_FOR_U,Q_SN,0,0,1,2,3,4,5,6,7,8,9,0.,0.,0.,A1,0,0,0,ILU_PR);\
}

void move(mg)	// vicemene zkusebni //
MULTIGRID *mg;
{
	GRID *topGrid = TOP_GRID(mg);
	NODE *np;

printf("\nMoving triangulation...\n");

printf("first component change...\n");
	SOLVE_LAP(U,glap01);

	scopy(topGrid,U,UU,T_FOR_U);	// store first component change to UU //

printf("second component change...\n");
	SOLVE_LAP(U,glap02);

printf("changing positions of vertices...\n\n");
	for (np = topGrid->firstN; np != NULL; np = np->succ) {
		np->myvertex->x[0] += NDS(np,UU);
		np->myvertex->x[1] += NDS(np,U);
	}
}*/

/*--------------------------------------*/
/* toto neni primo o POHYBU triangulace */
/*--------------------------------------*/


void node_order(pgrid, mask)
GRID *pgrid;
INT mask;
{
   NODE *pn, *pnodeif, *pnodeil, *pnodebf, *pnodebl;
 
   pn = FIRSTN(pgrid);

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
   pnodeil->succ = NULL;
   pgrid->firstNode = pnodeif;
   pgrid->lastNode = pnodeil;
   pgrid->fdbn = pnodebf;
   pgrid->ldbn = pnodebl;
}

void check_links(mg)
MULTIGRID *mg;
{
	GRID *pg;
	NODE *pn;
	LINK *pli;

	printf("\nChecking node to node LINKS..\n");
	for (pg = TOP_GRID(mg); pg != NULL; pg = pg->coarser) {
		for (pn = FDBN(pg); pn != FIRSTNODE(pg); pn = pn->succ) {
			if (IS_FN(pn)) eprintf("\nWrong (dir)node order!\n");
		}
		for (pn = FIRSTNODE(pg); pn != NULL; pn = pn->succ) {
			if (NOT_FN(pn)) eprintf("\nWrong (free)node order!\n");
			for (pli = pn->tstart; pli != pn->start && pli != NULL; pli = pli->next)
				if (IS_FN(pli->nbnode)) {
					eprintf("\nError in LINK to dir. node !!\n");
					printf("node: %e %e\n",pli->nbnode->myvertex->x[0],pli->nbnode->myvertex->x[1]);
				}
			if (pli == NULL && pn->start != NULL) {
				eprintf("\nWrong START LINK in node: %e %e\n",
					pn->myvertex->x[0], pn->myvertex->x[1]);
			}
			for (pli = pn->start; pli != NULL; pli = pli->next)
				if (NOT_FN(pli->nbnode)) {
					eprintf("\nError in LINK to free node !!\n");
					printf("node: %e %e\n",pli->nbnode->myvertex->x[0],pli->nbnode->myvertex->x[1]);
				}
		}
	}
	printf(".. checking completed\n\n");
}

#if (DATA_S & N_LINK_TO_NODES) && (DATA_S & N_LINK_TO_FACES) && (DATA_S & F_LINK_TO_NODES)
/* clears NMASK bit from free boundary nodes 
*/
void unmark_fbnodes(mg)
MULTIGRID *mg;
{
	GRID *pg, *tGrid = FIRSTGRID(mg);
	NODE *n, *nsucc, *n_pg, *nbn;
	FACE *nbf;
	LINK *pli, *nbpli, *ppli;
	NFLINK *nfl;
	FNLINK *fnl, *pfnl;

	printf("\nUnmarking free boundary nodes..\n");
	check_links(mg);
	for (n = FDBN(tGrid); n != FIRSTNODE(tGrid); ) {
		if (IS_FBN(n)) {
			NTYPE(n) &= ~NMASK;
			n_pg = n;
			nsucc = n->succ;
			for (pg = tGrid; pg != NULL; pg = pg->finer, n_pg = n_pg->son) {
				for (pli = TSTART(n_pg); pli != NULL; pli = NEXT(pli))
					if (IS_FN(pli->nbnode)) {
						START(n_pg) = pli;
						break;
					}
				printf("Correction of START done.\n");
				for ( ; pli != NULL; pli = NEXT(pli)) {
					nbn = pli->nbnode;
					for (nbpli = TSTART(nbn); nbpli->nbnode != n_pg; nbpli = NEXT(nbpli));
					if (START(nbn) == NEXT(nbpli)) {
						START(nbn) = nbpli;
					} else {
						if (TSTART(nbn) == nbpli) {
							TSTART(nbn) = NEXT(nbpli);
						} else {
							for (ppli = TSTART(nbn); NEXT(ppli) != nbpli; ppli = NEXT(ppli));
							NEXT(ppli) = NEXT(nbpli);
						}
						NEXT(nbpli) = NEXT(START(nbn));
						NEXT(START(nbn)) = nbpli;
					}
				}
				printf("LINKS done.\n");
				for (nfl = TNFSTART(n_pg); nfl != NULL; nfl = NEXT(nfl)) {
					nbf = nfl->nbface;
					for (fnl = TFNSTART(nbf); fnl->nbnode != n_pg; fnl = NEXT(fnl));
					if (FNSTART(nbf) == NEXT(fnl)) {
						FNSTART(nbf) = fnl;
					} else {
						if (TFNSTART(nbf) == fnl) {
							TFNSTART(nbf) = NEXT(fnl);
						} else {
							for (pfnl = TFNSTART(nbf); NEXT(pfnl) != fnl; pfnl = NEXT(pfnl));
							NEXT(pfnl) = NEXT(fnl);
						}
						NEXT(fnl) = NEXT(FNSTART(nbf));
						NEXT(FNSTART(nbf)) = fnl;
					}
				}
				printf("FNLINKS done.\n");
			//zmena seznamu dirichletovskych uzlu:
				if (FIRSTNODE(pg) == n_pg->succ) {
					FIRSTNODE(pg) = n_pg;
				} else {
					if (FDBN(pg) == n_pg) {
						FDBN(pg) = n_pg->succ;
					} else {
						for (nbn = FDBN(pg); nbn->succ != n_pg; nbn = nbn->succ);
						nbn->succ = n_pg->succ;
					}
					n_pg->succ = FIRSTNODE(pg)->succ;
					FIRSTNODE(pg)->succ = n_pg;
				}
				printf("New order of nodes done.\n\n");
			}
			if (n != FIRSTNODE(tGrid)) n = nsucc;
		} else {
			n = n->succ;
		}
	}
	check_links(mg);
	printf(".. unmark done\n");
}
#else
void unmark_fbnodes(mg)
MULTIGRID *mg;
{
	eprintf("unmark_fbnodes not available!\n");
}
#endif

/**/
void set_newvertex(tGrid)
GRID *tGrid;
{
	NODE *n;
	FACE *f;
	GRID *pg;

	for (n = FIRSTN(tGrid); n != NULL; n = n->succ) {
		n->newvertex->x[0] = n->myvertex->x[0];
		n->newvertex->x[1] = n->myvertex->x[1];
	}
	for (pg = tGrid; pg != NULL; pg = pg->coarser)
		for (f = FIRSTF(pg); f != NULL; f = f->succ)
			if (IS_CURVED_EDGE(f)) {
				f->newvertex->x[0] = f->myvertex->x[0];
				f->newvertex->x[1] = f->myvertex->x[1];
			}
}

FLOAT c_volume_old(f0,n0,n1,n2)
FACE *f0;
NODE *n0, *n1, *n2;
{
	FLOAT a, b, c, d, e, g, x12[2];

	a = n1->myvertex->x[0] - n0->myvertex->x[0];
	e = n1->myvertex->x[1] - n0->myvertex->x[1];
	g = n2->myvertex->x[0] - n0->myvertex->x[0];
	c = n2->myvertex->x[1] - n0->myvertex->x[1];
	AVERAGE(n1->myvertex->x,n2->myvertex->x,x12);
	b = 4*(f0->myvertex->x[0] - x12[0])*(c-e);
	d = 4*(f0->myvertex->x[1] - x12[1])*(a-g);
	return ( 0.5*fabs(a*c - e*g + (b + d)/3.) );
}

FLOAT c_volume(f0,n0,n1,n2)
FACE *f0;
NODE *n0, *n1, *n2;
{
	FLOAT a, b, c, d, e, g;

	a = n1->myvertex->x[0] - n0->myvertex->x[0];
	e = n1->myvertex->x[1] - n0->myvertex->x[1];
	g = n2->myvertex->x[0] - n0->myvertex->x[0];
	c = n2->myvertex->x[1] - n0->myvertex->x[1];
	b = n2->myvertex->x[1]*(f0->myvertex->x[0] - n1->myvertex->x[0]) 
		- n1->myvertex->x[1]*(f0->myvertex->x[0] - n2->myvertex->x[0]);
	d = f0->myvertex->x[1]*(a-g);
	return ( 0.5*fabs(a*c - e*g + 4.*(b + d)/3.) );
}

FLOAT volume_of_omega(tGrid)
GRID *tGrid;
{
	FLOAT vol = 0.0;
	ELEMENT *pel;
	NODE *n0, *n1, *n2;
	FACE *f0, *f1, *f2;

	for (pel = FIRSTELEMENT(tGrid); pel; pel = pel->succ) {
		NODES_OF_ELEMENT(n0,n1,n2,pel);
		FACES_OF_ELEMENT(f0,f1,f2,pel);
		if (IS_CURVED_EDGE(f0))
			vol += c_volume(f0,n0,n1,n2);
		else
		if (IS_CURVED_EDGE(f1))
			vol += c_volume(f1,n1,n2,n0);
		else
		if (IS_CURVED_EDGE(f2))
			vol += c_volume(f2,n2,n0,n1);
		else
			vol += VOLUME(pel);
	}
	return ( vol );
}

#if DATA_STR & BD_DATA

typedef FLOAT DVEC[DIM];

void BDdata_alloc(mg)
MULTIGRID *mg;
{
	GRID *tGrid;
	NODE *n;
	FACE *f;
	LINK *pli;
	NFLINK *pnf;
	FNLINK *pfn;

	for (tGrid = FIRSTGRID(mg); tGrid != NULL; tGrid = tGrid->finer) {
		for (n = FIRSTN(tGrid); n != NULL; n = n->succ)
			if (IS_ZNN(n)) {
				n->sbddata = (FLOAT*)calloc(NS_VECT, sizeof(FLOAT));
				n->ddiag = (DVEC*)calloc(B_MATR, sizeof(DVEC));
				for (pli = TSTART(n); pli != NULL; pli = pli->next)
					if (IS_ZNN(pli->nbnode))
						pli->ddiag = (DVEC*)calloc(B_MATR, sizeof(DVEC));
					else
						pli->ddiag = NULL;
				for (pnf = TNFSTART(n); pnf != NULL; pnf = pnf->next)
					if (IS_ZNF(pnf->nbface))
						pnf->ddiag = (DVEC*)calloc(B_MATR, sizeof(DVEC));
					else
						pnf->ddiag = NULL;
			} else {
				n->sbddata = NULL;
				n->ddiag = NULL;
				for (pli = TSTART(n); pli != NULL; pli = pli->next)
					pli->ddiag = NULL;
				for (pnf = TNFSTART(n); pnf != NULL; pnf = pnf->next)
					pnf->ddiag = NULL;
			}
		for (f = FIRSTF(tGrid); f != NULL; f = f->succ)
			if (IS_ZNF(f)) {
				f->sbddata = (FLOAT*)calloc(FS_VECT, sizeof(FLOAT));
				f->ddiag = (DVEC*)calloc(B_MATR, sizeof(DVEC));
				for (pfn = TFNSTART(f); pfn != NULL; pfn = pfn->next)
					if (IS_ZNN(pfn->nbnode))
						pfn->ddiag = (DVEC*)calloc(B_MATR, sizeof(DVEC));
					else
						pfn->ddiag = NULL;
			} else {
				f->sbddata = NULL;
				f->ddiag = NULL;
				for (pfn = TFNSTART(f); pfn != NULL; pfn = pfn->next)
					pfn->ddiag = NULL;
			}
	}
}

void BDdata_free(mg)
MULTIGRID *mg;
{
	GRID *tGrid;
	NODE *n;
	FACE *f;
	LINK *pli;
	NFLINK *pnf;
	FNLINK *pfn;

	for (tGrid = FIRSTGRID(mg); tGrid != NULL; tGrid = tGrid->finer) {
		for (n = FIRSTNODE(tGrid); n != NULL; n = n->succ)
			if (IS_ZNN(n)) {
				free(n->sbddata);
				free(n->ddiag);
				for (pli = TSTART(n); pli != NULL; pli = pli->next)
					if (IS_ZNN(pli->nbnode))
						free(pli->ddiag);
				for (pnf = TNFSTART(n); pnf != NULL; pnf = pnf->next)
					if (IS_ZNF(pnf->nbface))
						free(pnf->ddiag);
			}
		for (f = FIRSTFACE(tGrid); f != NULL; f = f->succ)
			if (IS_ZNF(f)) {
				free(f->sbddata);
				free(f->ddiag);
				for (pfn = TFNSTART(f); pfn != NULL; pfn = pfn->next)
					if (IS_ZNN(pfn->nbnode))
						free(pfn->ddiag);
			}
	}
}
#else

void BDdata_alloc(mg)
MULTIGRID *mg;
{ eprintf("Error: BDdata_alloc not available.\n"); }

void BDdata_free(mg)
MULTIGRID *mg;
{ eprintf("Error: BDdata_free not available.\n"); }
#endif

#endif	//DIM == 2
