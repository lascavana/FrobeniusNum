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

namespace path
{
  template<class T>
  T base_name(T const & path, T const & delims = "/")
  {
    return path.substr(path.find_last_of(delims) + 1);
  }
  template<class T>
  T remove_extension(T const & filename)
  {
    typename T::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
  }
} // namespace path


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
    if (prod != to_ZZ(0))
    {
      result = false; break;
    }
  }

  return result;
}

/* lattice reduction */
void reduce_ahl(
  ZZ            l,
  const vec_ZZ  &a,
  mat_ZZ        &Q,
  vec_ZZ        &x0
)
{
  mat_ZZ L;
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
  L[n][n+1] = -l*N2;

  /* lattice reduction */
  ZZ determ;
  LLL(determ, L, 99, 100, 0);
  
  /* find N1 in L matrix */
  if (L[n-1][n] != N1 && L[n-1][n] != -N1)
        throw std::runtime_error("ERROR: N1 not found. gcd(a) is not equal to 1!");

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

}

GRBModel* create_model(
  const mat_ZZ &Q,
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
    double obj = conv<double>(Q[0][j]);
    vars.push_back( model->addVar(-GRB_INFINITY, GRB_INFINITY, obj, GRB_INTEGER, "Lamb"+itos(j+1)) );
  }
  /* NOTE: the obj coefficient of variable j is a[0]*Q[0][j]. For numerical stability
  purposes, omit a[0] here and multiply it with the final result. **/

  /* create constraints */
  for (int i=0; i<nvars+1; i++)
  {
      GRBLinExpr expr = 0;
      for (int j=0; j<nvars; j++)
      {
          double coeff = -conv<double>(Q[i][j]);
          if (coeff != 0) expr += coeff * vars[j];
      }
      conss.push_back( model->addConstr(expr, GRB_LESS_EQUAL, 0.0, "CX_"+itos(i+1)) );
  }

  /* set tolerance */
  model->set(GRB_DoubleParam_OptimalityTol, 1e-9);

  return model;
}

ZZ solve(
  GRBModel*         model,
  vector<GRBConstr> &conss,
  const vec_ZZ      &xl,
  const vec_ZZ      &a
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

  /* objective shift */
  bestsol *= conv<double>(a[0]);

  /* add objective offset */
  ZZ offset = to_ZZ(0);
  for (int j=1; j<a.length(); j++) 
    offset += xl[j]*a[j];
  bestsol += conv<double>(offset);

  return to_ZZ(bestsol);

}

void write_prob(
  string        filename,
  const vec_ZZ  &a,
  ZZ            F
)
{
  int n = a.length();

  ofstream output_file(filename);
  output_file << "minimize +1 \n";
  output_file << "subject to" << "\n";
  output_file << "C1:";

  for (int j=0; j<n; j++)
    output_file << " +" << a[j] << " x" << j+1;

  output_file << " = " << F << endl;

  output_file << "Generals " << "\n";
  for (int j=0; j<n; j++)
    output_file << " x" << j+1;
  output_file << endl;

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

    ZZ z_star = to_ZZ(0);

    /* get row */
    vec_ZZ a = read_vec(argv[1]);

    /* get initial reduced basis */
    mat_ZZ Q; vec_ZZ x0;
    reduce_ahl(to_ZZ(1), a, Q, x0);

    /* create gurobi model */
    vector<GRBConstr> conss;
    GRBModel* model = create_model(Q, conss);
    
    /* calculate Frobenius number */
    for(ZZ l = to_ZZ(1); l < a[0] ; l++)
    {
      /* recalculate shift vector xl */
      vec_ZZ xl;
      reduce_ahl(l, a, Q, xl);

      /* solve problem */
      ZZ z = solve(model, conss, xl, a);

      /* save largest */
      if (z > z_star)
        z_star = z;

      /* print info */
      if (l % to_ZZ(100) == to_ZZ(0))
        cout << "Iteration " << l << ": best z "<< z_star - a[0] << endl;

    }
    ZZ F = z_star - a[0];
    cout << "Frobenius number is " << F << endl;

    /* write into an LP file */
    string filename {argv[1]};
    filename = path::remove_extension(filename) + ".lp";
    write_prob(filename, a, F);
   }