#ifndef CONNECTOR_H
#define CONNECTOR_H
#include <stdio.h>
#include <stdlib.h>
#include "../_aux.h"

typedef struct Connector Connector;

struct Connector {
  double* pdf;
  double* cdf;
  double mu;
  double sigma;
  int max_len;
};

void parse_con(Connector *con, double* matrix, double mu, double sigma, int max_len);
void print_con(Connector* con);

void print_gap(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed, double num, double den);
double score_con(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed);

#endif