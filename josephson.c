#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PROJECT_NAME "josephson"

#define ENABLE_L

double PLANCKBAR = 1.054571818e-34;
double ECHARGE = 1.602176634e-19;

double t_end = 2e-9;
int steps = 20000;

double freq = 10e9;
double amp = 0.0005e-3;
double cap = 20e-12;
double Ic = 0.1e-3;
double R = 1;

double Ib = 0.12e-3;

double ibf(double tao) {
  double tc = PLANCKBAR / (2 * ECHARGE * Ic * R);
  if(tao * tc > t_end/2) return 0;
  return (amp / Ic * sin(freq * tao * tc * 2 * M_PI) + 0.9035) ;
}

double * cal_phi(double ib, int steps, double tao_end, double beta_c) {

  double h = tao_end / steps;
  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
  phi_current = 0;
  phi_previous = 0;
  phi_next = 0;

  double phi_next_coeff = beta_c / (h * h) + 1.0 / (2.0 * h);
  for (int i = 0; i < steps; i++) {

    double ib_modified = (1 - exp(- i * h / tao_end * 10)) * ib; 

    phi_next = (ib_modified - sin(phi_current) + (1/(2 * h) - beta_c / h / h) * phi_previous + 2 * beta_c / h / h * phi_current) / phi_next_coeff;
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;

}

double * cal_phi_change_ib(double (* ibf)(double), int steps, double tao_end, double beta_c) {

  double h = tao_end / steps;
  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
  phi_current = 0;
  phi_previous = 0;
  phi_next = 0;

  double phi_next_coeff = beta_c / (h * h) + 1.0 / (2.0 * h);
  for (int i = 0; i < steps; i++) {
    phi_next = (ibf(i * h) - sin(phi_current) + (1/(2 * h) - beta_c / h / h) * phi_previous + 2 * beta_c / h / h * phi_current) / phi_next_coeff;
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;

}

int main(int argc, char **argv) {

  double tc = PLANCKBAR / (2 * ECHARGE * Ic * R);
  double tao_end = t_end / tc;
  double h = tao_end / steps;
  double beta_c = 2 * ECHARGE / PLANCKBAR * Ic * R * R * cap;

  FILE * fp = fopen("out.csv", "w");
  if (fp == NULL) {
    fprintf(stderr,"ERROR: UNABLE TO OPEN OUTPUT FILE, TERMINATING.");
    return 1;
  }

  fprintf(fp,"V, Ib\n");

  for(double Ib = 0; Ib < 2e-4; Ib += 0.1e-6) {
    double * phis = cal_phi(Ib / Ic, steps, tao_end, beta_c);
    // double V = (phis[19999] - phis[19997]) / (2 * h * tc) * PLANCKBAR / 2 / ECHARGE;

    double s = 0;
    double n = 0;
    for(int i = steps/2; i < steps - 1; i += 3) {
      s = s + phis[i+1] - phis[i-1];
      n++;
    }
    double V = s / n / (2 * h * tc) * PLANCKBAR / 2 / ECHARGE;

    fprintf(stdout, "Ib = %.32f\n", Ib);
    fprintf(fp, "%.32f, %.32f\n", V, Ib);
    free(phis);
  }

  fprintf(stdout, "CALCULATING BETA_C = %.2f", beta_c);
  double * phis = cal_phi_change_ib(ibf, steps, tao_end, beta_c);
  FILE * fp2 = fopen("out_driving.csv", "w");
  fprintf(fp2, "V, t\n0, 0\n");
  for(int i = 1; i < steps - 2; i++) {
    fprintf(fp2, "%.32f, %.32f\n", (phis[i+1] - phis[i-1]) / (2 * h * tc) * PLANCKBAR / 2 / ECHARGE , i * h * tc);
  }
  free(phis);

  fclose(fp);

  return 0;
}
