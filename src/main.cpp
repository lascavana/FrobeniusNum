#include <vector>
#include <cassert>
#include <sstream>
#include <iostream>

// gurobi
#include "gurobi_c++.h"

// NTL library
#include <NTL/ZZ.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>

using namespace std;
using namespace NTL;

string itos(int i) {stringstream s; s << i; return s.str(); }

bool check_basis(
  const vec_ZZ &a,
  const mat_ZZ &Q
)
{
  int n = a.length();
  bool result = true;

  for (int i=0; i<n-1; i++)
  {
    ZZ prod = -a[0]*Q[0][i];
    for (int j=1; j<n; j++)
      prod += a[j]*Q[j][i];
    if (prod != to_ZZ(0)) result = false;
  }

  return result;
}

/* lattice reduction */
void reduce_ahl(
  const vec_ZZ  &a,
  mat_ZZ        &Q,
  vec_ZZ        &x0
)
{
  mat_ZZ L, U;
  ZZ N1=to_ZZ(10000000);
  ZZ N2=to_ZZ(1000000000);

  int n = a.length();

  /* create L matrix */
  L.SetDims(n+1,n+2);
  for (int i=0;i<n+1;i++)
  {
     for(int j=0;j<n+2;j++)
         L[i][j]=0;
     L[i][i]=to_ZZ(1);
  }
  L[n][n]= N1;
  for (int i=1;i<n;i++){
    L[i][n+1] = a[i]*N2;
  }
  L[0][n+1] = -a[0]*N2;
  L[n][n+1] = -N2;

  /* lattice reduction */
  ZZ determ;
  LLL(determ, L, U, 99, 100, 0);
  

  /* find N1 in L matrix */
    /* Note: it should appear in position [n-p][n], where
       p is the number of linearly independent rows (therefore
       smaller or equal to m). */
  if (L[n-1][n] != N1 && L[n-1][n] != -N1)
        throw std::runtime_error("ERROR: N1 not found. Instance infeasible");

  /* get kernel lattice basis */
  Q.SetDims(n,n-1);
  for (int i=0;i<n-1; i++){
    for (int j=0;j<n;j++)
      Q[j][i]=L[i][j];
   }

  /* check basis */
  assert(check_basis(a, Q));

  x0.SetLength(n);
  for (int i=0; i<n; i++) x0[i] = L[n-1][i];

}

vec_ZZ get_objfun(
  const mat_ZZ &Q,
  const vec_ZZ &xl,
  const vec_ZZ &a
)
{
  int n = a.length();

  vec_ZZ c;
  c.SetLength(n);

  for (int i=0; i<n-1; i++)
    c[i] = a[0]*Q[0][i];

  c[n-1] = to_ZZ(0);
  for (int j=1; j<n; j++) 
    c[n-1] += xl[j]*a[j];

  return c;
}


double Solve(
  const mat_ZZ &Q,
  const vec_ZZ &xl,
  const vec_ZZ &c
)
{
    int nvars = Q.NumCols();

    GRBEnv env = GRBEnv();
    vector<GRBVar> vars;

    /* disable console output */
    env.set(GRB_IntParam_OutputFlag, 0);

    /* create model */
    GRBModel* model = new GRBModel(env);

    /* create variables */
    for (int j=0; j<nvars; j++)
    {
      double obj = conv<double>(c[j]);
      vars.push_back( model->addVar(-1e20, 1e20, obj, 'I', "Lamb"+itos(j+1)) );
    }

    /* create constraints */
    for (int i=0; i<nvars+1; i++)
    {
        GRBLinExpr expr = 0;
        for (int j=0; j<nvars; j++)
        {
            double coeff = -conv<double>(Q[i][j]);
            if (coeff != 0) expr += coeff * vars[j];
        }
        double rhs = conv<double>(xl[i]);
        model->addConstr(expr, '<', rhs, "CX_"+itos(i+1));
    }
    
    /* solve */
    model->optimize();

    /* get status */
    int status = model->get(GRB_IntAttr_Status);
    assert(status == 2); // solution should be optimal

    /* get solution */
    double bestsol = model->get(GRB_DoubleAttr_ObjVal);
    bestsol += conv<double>(c[nvars]);

    return bestsol;

}


int main(
    int                        argc,          /**< number of arguments from the shell */
    char**                     argv           /**< array of shell arguments */
)
{
    /* print usage */
    std::cout << " usage: frobenius <filename> " << std::endl;

    /* check usage */
    if( argc < 2 )
    {
        std::cout << " No file provided. Please enter a filename. " << std::endl;
        return 1;
    }

    /* get row */
    vec_ZZ a;
    a.SetLength(5);
    a[0] = to_ZZ(12223); a[1] = to_ZZ(12224);
    a[2] = to_ZZ(36674); a[3] = to_ZZ(61119);
    a[4] = to_ZZ(85569);

    /* get reduced basis for the problem */
    mat_ZZ Q; vec_ZZ x0;
    reduce_ahl(a, Q, x0);
    
    ZZ z_star = to_ZZ(0);
    for(ZZ l = to_ZZ(1); l < a[0] ; l++)
    {
      vec_ZZ xl = x0*l;

      /* get new objective function */
      vec_ZZ cl = get_objfun(Q, xl, a);

      /* solve problem */
      double z = Solve(Q, xl, cl);
      if (to_ZZ(z) > z_star)
      {
        z_star = z;
      }

      cout << z << endl;

    }

    ZZ F = z_star - a[0];
    cout << "Frobenius number is " << F << endl;
   }