#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libj.h"

double * cal_phi(double (* I_fun)(double, const void *), const void * context, const josephson_par * j, const sim_par * s) {
  double steps = s->steps;
  double t_end = s->t_end;

  double beta_c = 2 * ECHARGE / PLANCKBAR * j->Ic * j->R * j->R * j->C;
  double tc = PLANCKBAR / (2 * ECHARGE * j->Ic * j->R);
  double tao_end = t_end / tc;
  double h = tao_end / steps;
  double Ic = j->Ic;

  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
  phi_current = s->phi_0;
  phi_previous = s->phi_p;
  phi_next = 0;

  double phi_next_coeff = beta_c / (h * h) + 1.0 / (2.0 * h);
  for (int i = 0; i < steps; i++) {
    phi_next = (I_fun(i * h * tc, context) / Ic - sin(phi_current) + (1/(2 * h) - beta_c / h / h) * phi_previous + 2 * beta_c / h / h * phi_current) / phi_next_coeff;
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;

}