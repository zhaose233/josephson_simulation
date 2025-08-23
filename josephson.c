#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PROJECT_NAME "josephson"

#define ENABLE_L

double PLANCKBAR = 1.054571818e-34;
double ECHARGE = 1.602176634e-19;

double dt = 1e-13;

double freq = 10e9;
double amp = 0.005e-3;
double cap = 20e-12;
double Ic = 0.1e-3;
double R = 1;

double Ib = 0.12e-3;

double driving_I(double t) {
  //return amp * sin(freq * t * 2 * M_PI);
    return 0.12e-3;
}

double * cal_phi(double Ib, int steps) {

  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
    phi_current = 0;
    phi_previous = 0;
    phi_next = 0;

  double a_C = cap / dt / dt;
  double a_R = 1.0 / 2.0 / dt / R;
  double a_h = 2.0 * ECHARGE / PLANCKBAR;

  for (int i = 0; i < steps; i++) {
    phi_next = (a_h * (Ib - Ic * sin(phi_current)) + (a_R - a_C) * phi_previous + a_C * 2 * phi_current) / (a_C + a_R);
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;

}

double * cal_phi_change_Ib(double (* Ibf)(double), int steps) {

  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
    phi_current = 0;
    phi_previous = 0;
    phi_next = 0;

  double a_C = cap / dt / dt;
  double a_R = 1.0 / 2.0 / dt / R;
  double a_h = 2.0 * ECHARGE / PLANCKBAR;

  for (int i = 0; i < steps; i++) {
    phi_next = (a_h * (Ibf(i * dt) - Ic * sin(phi_current)) + (a_R - a_C) * phi_previous + a_C * 2 * phi_current) / (a_C + a_R);
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;

}

double Ibf(double t) {
  return 0.1e-3 + amp * sin(freq * t * 2 * M_PI);
}

int main(int argc, char **argv) {
    FILE * fp = fopen("out.csv", "w");
    if (fp == NULL) {
      fprintf(stderr,"ERROR: UNABLE TO OPEN OUTPUT FILE, TERMINATING.");
      return 1;
    }

    fprintf(fp,"V, Ib\n");

    for(double Ib = 0; Ib < 2e-4; Ib += 0.1e-6) {
      double * phis = cal_phi(Ib, 20000);
      double V = (phis[19999] - phis[19997]) / 2 / dt * PLANCKBAR / 2 / ECHARGE;
      fprintf(stdout, "Ib = %.32f\n", Ib);
      fprintf(fp, "%.32f, %.32f\n", V, Ib);
      free(phis);
    }

    double * phis = cal_phi_change_Ib(Ibf, 20000);
    FILE * fp2 = fopen("out_driving.csv", "w");
    fprintf(fp2, "V, t\n0, 0\n");
    for(int i = 1; i < 19999; i++) {
      fprintf(fp2, "%.32f, %.32f\n", (phis[i] - phis[i-1]) * PLANCKBAR / 2 / ECHARGE / dt , i * dt);
    }
    
    
    

    fclose(fp);

    return 0;
}
