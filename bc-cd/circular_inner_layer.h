// ======================================================================
// circular inner layer
// John, Maubach, Tobiska, Num. Math. 1997
// ======================================================================
#define __TWO_INTERIOR_LAYERS__
#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: JohnMaubachTobiska.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{
 double d,r0,x0,y0,rxy,cxy,hxy,g1xy,g1xy_x,g1xy_xx,g1xy_y,g1xy_yy;
 double scale,g2xy,g2xy_x,g2xy_xx,g2xy_y,g2xy_yy,gxy,gxy_x,gxy_xx,gxy_y,gxy_yy;
 double gxy_l,v_x,v_y,c_m,para;

  para = TDatabase::ParamDB->RE_NR;
  d   =  2.0 * sqrt(para);
  r0  =  0.25;
  x0  =  0.5;
  y0  =  0.5;
  rxy =  (x-x0)*(x-x0) + (y-y0)*(y-y0);
  cxy =  d*(r0*r0 - rxy);
  hxy =  2.0*d / (1 + cxy*cxy);
  g1xy    =  atan(cxy) + Pi/2.0;
  g1xy_x  = -hxy*(x-x0);
  g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
  g1xy_y  = -hxy*(y-y0);
  g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
  g1xy    =  g1xy / Pi; 
  g1xy_x  =  g1xy_x / Pi;
  g1xy_xx =  g1xy_xx / Pi;
  g1xy_y  =  g1xy_y / Pi;
  g1xy_yy =  g1xy_yy / Pi;

  scale   =  16.0;
  g2xy    =  scale*x*(x - 1)*y*(y - 1);
  g2xy_x  =  scale*(2*x - 1)*y*(y - 1);
  g2xy_xx =  scale*2*y*(y - 1);
  g2xy_y  =  scale*(2*y - 1)*x*(x - 1);
  g2xy_yy =  scale*2*x*(x - 1);

  values[0] =  g1xy*g2xy;
  values[1] =  g1xy_x*g2xy + g1xy*g2xy_x;
  gxy_xx  =  g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx;
  values[2] =  g1xy_y*g2xy + g1xy*g2xy_y;
  gxy_yy  =  g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy;

  values[3]  =  gxy_xx + gxy_yy;
 }

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;
  double d,r0,x0,y0,rxy,cxy,hxy,g1xy,g1xy_x,g1xy_xx,g1xy_y,g1xy_yy;
  double scale,g2xy,g2xy_x,g2xy_xx,g2xy_y,g2xy_yy,gxy,gxy_x,gxy_xx,gxy_y,gxy_yy;
  double gxy_l,v_x,v_y,c_m;

  d   =  2.0 * sqrt(1.0/eps);
  r0  =  0.25;
  x0  =  0.5;
  y0  =  0.5;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];

    rxy =  (x-x0)*(x-x0) + (y-y0)*(y-y0);
    cxy =  d*(r0*r0 - rxy);
    hxy =  2.0*d / (1 + cxy*cxy);
    g1xy    =  atan(cxy) + Pi/2.0;
    g1xy_x  = -hxy*(x-x0);
    g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
    g1xy_y  = -hxy*(y-y0);
    g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
    g1xy    =  g1xy / Pi; 
    g1xy_x  =  g1xy_x / Pi;
    g1xy_xx =  g1xy_xx / Pi;
    g1xy_y  =  g1xy_y / Pi;
    g1xy_yy =  g1xy_yy / Pi;
    
    scale   =  16.0;
    g2xy    =  scale*x*(x - 1)*y*(y - 1);
    g2xy_x  =  scale*(2*x - 1)*y*(y - 1);
    g2xy_xx =  scale*2*y*(y - 1);
    g2xy_y  =  scale*(2*y - 1)*x*(x - 1);
    g2xy_yy =  scale*2*x*(x - 1);

    gxy     =  g1xy*g2xy;
    gxy_x   =  g1xy_x*g2xy + g1xy*g2xy_x;
    gxy_xx  =  g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx;
    gxy_y   =  g1xy_y*g2xy + g1xy*g2xy_y;
    gxy_yy  =  g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy;

    gxy_l  =  gxy_xx + gxy_yy;
  
    coeff[0] = eps;
    coeff[1] = 2;
    coeff[2] = 3;
    coeff[3] = 0;
    coeff[4] = -coeff[0]*gxy_l+coeff[1]*gxy_x+coeff[2]*gxy_y+coeff[3]*gxy; 
  }
}
void ComputeExtremalValues(int N, double *sol, double  *values)
{
   int i;
   double max, min;

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i] > max)
         max = sol[i];
      if(sol[i] < min)
         min = sol[i];
   }

   values[0] = min;
   values[1] = max;
}
   
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(int level, TFEFunction2D *ufct)
{
    /*const int max_bound_points = 10000;
  int i,j,k, N_Cells;
  double xi, eta;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  TFESpace2D *FESpace2D;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, *u1;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  
  int *Numbers, comp, found, N_Edges;
  double u, ux, uy, x, y, *Values;
  double val, x_min, val_min,  x_max, val_max, eps=1e-10;
  double y_coord[max_bound_points], uval[max_bound_points], min;
  int *GlobalNumbers, *BeginIndex, bound_points, index;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                   // find edges on boundary part 3  
      joint=cell->GetJoint(j);          // this is x=0   
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==3)
        {
          found = 1;
          break;
        }
      }
    }

    if (!found) continue;

    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    for (j=0;j<N_Edges;j++)
    {
      // check all vertices if they are on the boundary 
      x = cell->GetVertex(j)->GetX();
      // point not on boundary
      if (fabs(x)>eps)
         continue;    
      y = cell->GetVertex(j)->GetY();
      //check if this node is already in the array of boundpoints
      found = 0;
      for (k=bound_points-1;k>=0; k--)
      {
         if (fabs(y- y_coord[k]) < eps)
         {
            found = 1;
            break;
         }
      }
      if (found)
         continue;
      // new node
      y_coord[bound_points] = y;
      bound_points++;      
      if ( bound_points > max_bound_points)
      {
         OutPut("TwoInteriorLeyers.h: maximal number of boundary points reached !!!" << endl);
         exit(4711);
      }
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, (TGridCell *)cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);

      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      uval[bound_points-1] = u;
    } // endfor 
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
  // order the arrays

  for (i=0;i<bound_points; i++)
  {
     min = 1e6;
     for(j=i;j<bound_points; j++)
     {
        if (y_coord[j]< min)
        {
           min = y_coord[j];
           index = j;
        }
     }
     // change the entries
     y_coord[index] = y_coord[i];
     y_coord[i] = min;
     val = uval[i];
     uval[i] = uval[index];
     uval[index] = val;
  }
  for (i=0;i<bound_points; i++)
  {
     OutPut("outflow " << level << " " <<  y_coord[i] << 
            " " <<  uval[i] << endl);
  }
    */
    return;
  double h, x=1,values[3],y, max = 0;
  int i, bound_points = 200;
  h = 1.0/bound_points;
  for (i=0;i<=bound_points; i++)
  {
      y = i*h;
      ufct->FindGradient(x,y,values);
      if (values[0] > max)
	  max = values[0];
  }
      OutPut("max on x = 1 " << level << " " << 
            max << endl);
}
/** compute curve of the outflow boundary */
void ComputeDiagonal(int level, TFEFunction2D *ufct)
{
  const int max_bound_points = 10000;
  int i,j,k, N_Cells;
  double xi, eta;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  TFESpace2D *FESpace2D;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, *u1;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  return;
  int *Numbers, comp, found, N_Edges, N_add, diag =0;
  double u, ux, uy, x, y, *Values, x_add;
  double val, x_min, val_min,  x_max, val_max, eps=1e-10;
  double y_coord[max_bound_points], uval[max_bound_points], min;
  int *GlobalNumbers, *BeginIndex, bound_points, index;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;

    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    x_add = 0;
    N_add = 0;
    diag = 0;
    for (j=0;j<=N_Edges;j++)
    {
      // check all vertices if they are on the boundary 
      if (j <  N_Edges)
      {
         x = cell->GetVertex(j)->GetX();
         y = cell->GetVertex(j)->GetY();
         //check if this node is already in the array of boundpoints
         if (fabs(x-y)>eps)
            continue;
         x_add += x;
         N_add++;
      }
      else
      {
         // no mesh cell on the diagonal
         if (!diag)
            continue;
         else
         {
            x = x_add/N_add;
            y = x;
         }
      }
      diag++;
      found = 0;
      for (k=bound_points-1;k>=0; k--)
      {
         if (fabs(y- y_coord[k]) < eps)
         {
            found = 1;
            break;
         }
      }
      if (found)
         continue;
      // new node
      y_coord[bound_points] = y;
      bound_points++;      
      if ( bound_points > max_bound_points)
      {
         OutPut("TwoInteriorLeyers.h: maximal number of boundary points reached !!!" << endl);
         exit(4711);
      }
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, (TGridCell *)cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);

      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      uval[bound_points-1] = u;
    } // endfor 
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
  // order the arrays

  for (i=0;i<bound_points; i++)
  {
     min = 1e6;
     for(j=i;j<bound_points; j++)
     {
        if (y_coord[j]< min)
        {
           min = y_coord[j];
           index = j;
        }
     }
     // change the entries
     y_coord[index] = y_coord[i];
     y_coord[i] = min;
     val = uval[i];
     uval[i] = uval[index];
     uval[index] = val;
  }
  for (i=0;i<bound_points; i++)
  {
     OutPut("diagonal " << level << " " <<  y_coord[i]*sqrt(2.0) << 
            " " <<  uval[i] << endl);
  }
}
