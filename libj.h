#define PLANCKBAR 1.054571818e-34
#define ECHARGE 1.602176634e-19

/* typedef struct _josephson_par {
  double C;
  double Ic;
  double R;
} josephson_par; */

typedef struct _j_ode_par {
  double B_c;
  double (* i_fun)(double, const void *);
  const void * i_fun_context;
} j_ode_par;

typedef struct _sim_par {
  double start;
  double end;
  int steps;
  double phi_0;
  double phi_p;
  double omega_0;
} sim_par;

typedef struct _j_ode_return {
  double * phis;
  double * omegas;
} j_ode_return;

j_ode_return run_josephson_ode(j_ode_par * jop, const sim_par *s) ;
/* 
double * run_josephson_phi_SM(double (* I_fun)(double, const void *), const void * context, const josephson_par * j, const sim_par * s) ;

double * run_josephson_phi_PF(double (*I_fun)(double, const void *), const void *context, const josephson_par *j, const sim_par *s) ;
 */