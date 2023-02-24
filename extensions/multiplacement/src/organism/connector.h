#include <stdio.h>
#include <stdlib.h>
#include "../_aux.h"
#ifndef CONNECTOR_H
#define CONNECTOR_H

typedef struct Connector Connector;

struct Connector {
  float* pdf;
  float* cdf;
  float mu;
  float sigma;
  int max_len;
};

void parse_con(Connector *con, float* matrix, float mu, float sigma, int max_len, bool precomputed){
  if (precomputed) {
    con->pdf = matrix;
    con->cdf = matrix + max_len;
    con->mu = 0.0;
    con->sigma = 0.0;
    con->max_len = max_len;
  } else {
    con->pdf = NULL;
    con->cdf = NULL;
    con->mu = mu;
    con->sigma = sigma;
    con->max_len = 0;
  }
}

float score_con(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed) {
  if (precomputed){
    float num = con->pdf[gap];
    float auc = con->cdf[s_len - 1] - con->cdf[0];
    if (auc > 1E-10) {
      num /= auc;
    } else {
      num /= 1E-5;
    }

    float den = (DEN_EXPANSION[(eff_len - (gap + 1)) - 1] -
                 DEN_EXPANSION[(n_recs - 1) - 1] -
                 DEN_EXPANSION[eff_len - (gap + 1) - (n_recs) - 1]) -
                (DEN_EXPANSION[eff_len - 1] -
                 DEN_EXPANSION[n_recs - 1] -
                 DEN_EXPANSION[eff_len - n_recs - 1]);
    //printf("num: %f, den: %f\n", num, den);
    if (num - den > -1E10){
      return num - den;
    }else{
      return -1E10;
    }
  }else{
    float num = get_numerator(s_len, gap, con->mu, con->sigma);
    float den = get_denominator(gap + 1, n_recs, eff_len);
    return log2f(num/den);
  }
}

#endif