#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include <fstream>
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

/* read vector a from file */
vec_ZZ read_vec(
  const char *filename
)
{
  vec_ZZ a; a.SetLength(100000);

  /* open file for reading */
  ifstream input_file(filename);

  /* read data */
  if (input_file.is_open())
  {
    /* read a vector */
    ZZ ai; int idx = 0;
    while (input_file >> ai) {
      a[idx] = ai;
      idx++;

      if (idx == 100000)
      {
        cerr << "Vector a too long, must have less than 100000 components" << endl;
        break;
      }
    } 
    input_file.close();

    a.SetLength(idx);
  }
  else {cerr << "Unable to open file" << endl;}

  return a;
}


/* check that Q is basis for ker(a) */
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

/* get new objective function */
vec_ZZ get_objfun(
  const mat_ZZ &Q,
  const vec_ZZ &a
)
{
  int n = a.length();

  vec_ZZ c;
  c.SetLength(n-1);

  for (int i=0; i<n-1; i++)
    c[i] = a[0]*Q[0][i];

  return c;
}

/* lattice reduction */
void reduce_ahl(
  const vec_ZZ  &a,
  mat_ZZ        &Q,
  vec_ZZ        &x0,
  vec_ZZ        &c
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

  /* get shift vector */
  x0.SetLength(n);
  for (int i=0; i<n; i++) x0[i] = L[n-1][i];

  /* get new objective function */
  c = get_objfun(Q, a);

}

GRBModel* create_model(
  const mat_ZZ &Q,
  const vec_ZZ &c,
  vector<GRBConstr> &conss
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
      conss.push_back( model->addConstr(expr, '<', 0.0, "CX_"+itos(i+1)) );
  }

  return model;
}

double solve(
  GRBModel* model,
  vector<GRBConstr> &conss,
  const vec_ZZ &xl,
  const vec_ZZ &a
)
{
  /* change constraint rhs */
  for (int i=0; i<conss.size(); i++)
    conss[i].set(GRB_DoubleAttr_RHS, conv<double>(xl[i]));
    
  /* solve */
  model->optimize();

  /* get status */
  int status = model->get(GRB_IntAttr_Status);
  assert(status == 2); // solution should be optimal

  /* get solution */
  double bestsol = model->get(GRB_DoubleAttr_ObjVal);

  /* add objective offset */
  ZZ offset = to_ZZ(0);
  for (int j=1; j<a.length(); j++) 
    offset += xl[j]*a[j];
  bestsol += conv<double>(offset);

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

    vec_ZZ aa = read_vec(argv[1]);


    /* get row */
    vec_ZZ a = read_vec(argv[1]);

    /* get reduced basis for the problem */
    mat_ZZ Q; vec_ZZ x0; vec_ZZ c;
    reduce_ahl(a, Q, x0, c);

    /* create gurobi model */
    vector<GRBConstr> conss;
    GRBModel* model = create_model(Q, c, conss);
    
    ZZ z_star = to_ZZ(0);
    for(ZZ l = to_ZZ(1); l < a[0] ; l++)
    {
      vec_ZZ xl = x0*l;

      /* solve problem */
      double z = solve(model, conss, xl, a);

      /* save largest */
      if (to_ZZ(z) > z_star)
        z_star = z;

      /* print info */
      if (l % to_ZZ(100) == to_ZZ(0))
        cout << "Iteration " << l << ": best z "<< z_star - a[0] << endl;

    }

    ZZ F = z_star - a[0];
    cout << "Frobenius number is " << F << endl;
   }