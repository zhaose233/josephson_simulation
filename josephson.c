#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libj.h"

#define PROJECT_NAME "josephson"

double t_end = 2e-9;
int steps = 20000;

double freq = 10e9;
double amp = 0.0005e-3;
double cap = 20e-12;
double Ic = 0.1e-3;
double R = 1;

//double Ib = 0.12e-3;

typedef struct _ac_par {
  double amp;
  double freq;
} ac_par;

double I_fun_ac_impl(double t, const void * context) {
  ac_par * a = (ac_par *)context;
  double amp = a->amp;
  double freq = a->freq;

  return amp * sin(freq * t * 2 * M_PI);
}

double I_fun_dc_impl(double t, const void * context) {
  double Ib = *(double *)context;

  return Ib;
}

int main(int argc, char **argv) {

  josephson_par jp = {cap, Ic, R};

  FILE * fp = fopen("out.csv", "w");
  if (fp == NULL) {
    fprintf(stderr,"ERROR: UNABLE TO OPEN OUTPUT FILE, TERMINATING.");
    return 1;
  }

  fprintf(fp,"V, Ib\n");

  sim_par dc_sim_par = {t_end, steps, 0, 0};

  for(int i = 0; i < 401; i ++) {
    double Ib;
    if (i < 201) 
      Ib = i * 1e-6;
    else 
      Ib = (- i + 400) * 1e-6;

    double * phis = cal_phi(I_fun_dc_impl, &Ib, &jp, &dc_sim_par);
    double tc = PLANCKBAR / (2 * ECHARGE * jp.Ic * jp.R);
    double tao_end = dc_sim_par.t_end / tc;
    double h = tao_end / dc_sim_par.steps;

    double s = 0;
    double n = 0;
    for(int i = steps/2; i < steps - 1; i += 3) {
      s = s + phis[i+1] - phis[i-1];
      n++;
    }
    double V = s / n / (2 * h * tc) * PLANCKBAR / 2 / ECHARGE;

    dc_sim_par.phi_0 = phis[dc_sim_par.steps - 1];
    dc_sim_par.phi_p = phis[dc_sim_par.steps - 2];

    fprintf(stdout, "Ib = %.32f\n", Ib);
    fprintf(fp, "%.32f, %.32f\n", V, Ib);
    free(phis);
  }

/*   fprintf(stdout, "CALCULATING BETA_C = %.2f", beta_c);
  double * phis = cal_phi_change_ib(ibf, steps, tao_end, beta_c);
  FILE * fp2 = fopen("out_driving.csv", "w");
  fprintf(fp2, "V, t\n0, 0\n");
  for(int i = 1; i < steps - 2; i++) {
    fprintf(fp2, "%.32f, %.32f\n", (phis[i+1] - phis[i-1]) / (2 * h * tc) * PLANCKBAR / 2 / ECHARGE , i * h * tc);
  }
  free(phis); */

  fclose(fp);

  return 0;
}
