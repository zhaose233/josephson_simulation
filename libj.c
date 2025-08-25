#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "libj.h"

#define ODE_HSTART 1e-6
#define ODE_EPSABS 1e-6
#define ODE_ESPREL 1e-6

int j_ode_system(double tau, const double * y, double * f, void * par) {
  j_ode_par * p = (j_ode_par *)par;
  
  double phi = y[0];
  double omega = y[1];

  double current = p->i_fun(tau, p->i_fun_context);
  
  f[0] = omega;
  f[1] = current - sin(phi) - p->B_c * omega;

  return GSL_SUCCESS;
}

j_ode_return run_josephson_ode(j_ode_par * jop, const sim_par *s) {
  double steps = s->steps;
  double tau_start = s->start;
  double tau_end = s->end;
  double h = (tau_end - tau_start) / steps;
  
  gsl_odeiv2_system sys = {j_ode_system, NULL, 2, jop};

  gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
    ODE_HSTART, ODE_EPSABS, ODE_ESPREL);

  double * phis = malloc(sizeof(double) * (steps + 1));
  double * omegas = malloc(sizeof(double) * (steps + 1));

  double y[2] = {s->phi_0, s->omega_0};
  double tau = 0;
  j_ode_return r = {NULL, NULL};

  for(int i = 0; i <= steps; i++) {
    const double tau_i = i * h + tau_start;
    int status = gsl_odeiv2_driver_apply(driver, &tau, tau_i, y);
    if(status != GSL_SUCCESS) {
      fprintf(stderr, "SIMULATION ERROR\n");
      return r;
    }
    phis[i] = y[0];
    omegas[i] = y[1];
  }

  r.phis = phis;
  r.omegas = omegas;
  return r;

}

/* double * run_josephson_phi_SM(double (* I_fun)(double, const void *), const void * context,
                              const josephson_par * j, const sim_par * s) {
  double steps = s->steps;
  double t_end = s->end;

  double beta_c = 2 * ECHARGE / PLANCKBAR * j->Ic * j->R * j->R * j->C;
  double tc = PLANCKBAR / (2 * ECHARGE * j->Ic * j->R);
  double tau_end = t_end / tc;
  double h = tau_end / steps;
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

double * run_josephson_phi_PF(double (*I_fun)(double, const void *), const void *context, const josephson_par *j, const sim_par *s) {
  double steps = s->steps;
  double t_end = s->end;

  double inv_sqrt_beta_c =  1 / sqrt(2 * ECHARGE / PLANCKBAR * j->Ic * j->R * j->R * j->C);
  double Ic = j->Ic;
  double omega_p = sqrt(2 * ECHARGE * Ic / (PLANCKBAR * j->C));
  double tau_end = t_end * omega_p;
  double h = tau_end / steps;

  double * phis = malloc(sizeof(double) * steps);

  double phi_current, phi_next, phi_previous;
  phi_current = s->phi_0;
  phi_previous = s->phi_p;
  phi_next = 0;

  double phi_next_coeff = 1 / (h * h) + inv_sqrt_beta_c / (2.0 * h);
  for (int i = 0; i < steps; i++) {
    phi_next = (I_fun(i * h / omega_p, context) / Ic - sin(phi_current) + (inv_sqrt_beta_c /(2 * h) - 1 / h / h) * phi_previous + 2 / h / h * phi_current) / phi_next_coeff;
    phi_previous = phi_current;
    phi_current = phi_next;
    phis[i] = phi_current;
  }
  
  return phis;
} */
