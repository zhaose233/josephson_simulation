#define PLANCKBAR 1.054571818e-34
#define ECHARGE 1.602176634e-19

typedef struct _josephson_par {
  double C;
  double Ic;
  double R;
} josephson_par;

typedef struct _sim_par {
  double t_end;
  int steps;
  double phi_0;
  double phi_p;
} sim_par;

double * cal_phi(double (* I_fun)(double, const void *), const void * context, const josephson_par * j, const sim_par * s) ;
