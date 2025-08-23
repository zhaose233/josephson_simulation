#include <stdio.h>
#include <math.h>

#define PROJECT_NAME "josephson"

//#define ENABLE_L

double PLANCKBAR = 1.054571818e-34;
double ECHARGE = 1.602176634e-19;

double dt = 1e-14;

double freq = 10e9;
double amp = 0.105e-3;
double cap = 20e-12;
double Rs = 1;
double Ic = 0.1e-3;
double ind = 8.6e-12;

double josephoson_r(double v) {
  return 1;
}

double driving_I(double t) {
  return amp * sin(freq * t * 2 * M_PI);
  // return 0.12e-3;
}

int main(int argc, char **argv) {
    FILE * fp = fopen("out.csv", "w");
    if (fp == NULL) {
      fprintf(stderr,"ERROR: UNABLE TO OPEN OUTPUT FILE, TERMINATING.");
      return 1;
    }

    fprintf(fp,"V, Is, phi, t, I\n");

    double V,Is,phi,t;
    V = 0;
    Is = 0;
    phi = 0;
    t = 0;

    for(int i = 0; i < 200000; i++) {
      fprintf(fp, "%.32f, %.32f, %.32f, %.32f, %.32f\n", V, Is, phi, t, driving_I(t));
      V += (driving_I(t) - Is - Ic * sin(phi) - V / josephoson_r(V)) / cap * dt;
      phi += V / PLANCKBAR * 2 * ECHARGE * dt;
      #ifdef ENABLE_L
      Is += (V - Is * Rs) / ind * dt;
      #endif
      t += dt;
    }

    fprintf(fp, "%.32f, %.32f, %.32f, %.32f, %.32f\n", V, Is, phi, t, driving_I(t));
    fclose(fp);

    return 0;
}
