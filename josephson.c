#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libj.h"

#include <gsl/gsl_rng.h>

#define PROJECT_NAME "josephson"

//double Ib = 0.12e-3;

typedef struct _ac_par {
  double amp;
  double freq;
} ac_par;

void run_noise(void);

void run_IV_curve(double B_c) ;

void run_chaos_sim(sim_par sp, j_ode_par jop, ac_par ap, char * csv_file) ;

double i_fun_ac_impl(double t, const void * context) {
  ac_par * a = (ac_par *)context;
  double amp = a->amp;
  double w = a->freq;

  return amp * sin(w * t);
}

int main(void) {
  //run_IV_curve();
  double tau_end = 800;
  double tau_start = 700;
  double steps = 10000;
  double h = (tau_end - tau_start) / steps;
  double B_c = 0.5;
  double w = 0.66;
  // double amp = 1.212;
  
  ac_par ap = {0, w};
  j_ode_par jop = {B_c,i_fun_ac_impl, &ap};
  sim_par sp = {tau_start, tau_end, steps, 0, 0, 0};

  double amps[] = {0.8, 1.0636, 1.0652, 1.08, 1.212, 1.32};
  for(int i = 0; i < 6; i++){
    ap.amp = amps[i];
    char * csv_file_suffix;
    asprintf(& csv_file_suffix, "%.4f", amps[i]);
    run_chaos_sim(sp, jop, ap, csv_file_suffix);
    free(csv_file_suffix);
  }

  run_IV_curve(0.5);
  run_IV_curve(1.5);

  return 0;
}

void run_chaos_sim(sim_par sp, j_ode_par jop, ac_par ap, char * csv_file_suffix) {

  double tau_end = sp.end;
  double tau_start = sp.start;
  double steps = sp.steps;
  double h = (tau_end - tau_start) / steps;

  j_ode_return r = run_josephson_ode(&jop, &sp);

  char * main_csv_file;
  asprintf(&main_csv_file, "chaos_%s.csv", csv_file_suffix);
  FILE * main_fp = fopen(main_csv_file, "w");
  free(main_csv_file);

  char * poincare_csv_file;
  asprintf(&poincare_csv_file, "chaos_poincare_%s.csv", csv_file_suffix);
  FILE * poincare_fp = fopen(poincare_csv_file, "w");
  free(poincare_csv_file);

  fprintf(main_fp, "phi, omega, tau\n");
  fprintf(poincare_fp, "phi, omega\n");

  for(int i = 0; i <= steps; i++) {
    double tau = tau_start + i * h;

    fprintf(main_fp, "%.17g, %.17g, %.17g\n", r.phis[i], r.omegas[i], tau);
  }

  int n_previous = (tau_start * ap.freq) / (2 * M_PI);
  int n_current;
  for(int i = 1; i <= steps; i++) {
    double tau = tau_start + i * h;
    n_current = (tau * ap.freq) / (2 * M_PI);

    if (n_current != n_previous) {
      fprintf(poincare_fp, "%.17g, %.17g\n", r.phis[i], r.omegas[i]);
    }

    n_previous = n_current;
  }

  free(r.omegas);
  free(r.phis);

  fclose(main_fp);
  fclose(poincare_fp);

}

double i_fun_dc_impl(double _, const void * context) {
  double ib = *(double *)context;
  return ib;
}

void run_IV_curve(double B_c) {

  double tau_end = 400;
  int steps = 1000;

  double ib;
  j_ode_par jop = {B_c, i_fun_dc_impl, &ib};

  char * csv_file;
  asprintf(&csv_file, "iv_%.3f.csv", B_c);
  FILE * fp = fopen(csv_file, "w");
  free(csv_file);

  if (fp == NULL) {
    fprintf(stderr,"ERROR: UNABLE TO OPEN OUTPUT FILE, TERMINATING.");
    return;
  }

  fprintf(fp,"omega, ib\n");

  sim_par dc_sim_par = {0,tau_end, steps, 0, 0, 0};

  const double ib_max = 2;
  const int NUM_POINTS= 200;
  const int TOTAL_POINTS = 4 * NUM_POINTS + 2;
  for (int i = 0; i < TOTAL_POINTS; i++) {
    if (i <= NUM_POINTS) {
      ib = i * ib_max / NUM_POINTS;
    } else if (i <= (NUM_POINTS * 2 + 1)) {
      ib = ib_max / NUM_POINTS * (NUM_POINTS * 2 + 1 - i);
    } else if (i <= (NUM_POINTS * 3 + 1)) {
      ib = - ib_max / NUM_POINTS * (i - NUM_POINTS * 2 - 1);
    } else {
      ib = - ib_max / NUM_POINTS * (NUM_POINTS * 4 + 2 - i);
    }

    fprintf(stdout, "ib = %.17e\n", ib);

    j_ode_return r = run_josephson_ode(&jop, &dc_sim_par);

    double s = 0;
    const int OMEGA_POINTS = steps / 2;
    for(int i = steps; i > steps - OMEGA_POINTS; i--) {
      s = s + r.omegas[i];
    }
    double V = s / OMEGA_POINTS;

    dc_sim_par.omega_0 = r.omegas[steps];
    dc_sim_par.phi_0 = r.phis[steps];

    fprintf(fp, "%.17e, %.17e\n", V, ib);
    free(r.phis);
    free(r.omegas);
  }

  fclose(fp);
}
