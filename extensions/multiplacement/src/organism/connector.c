#include "connector.h"
void parse_con(Connector *con, double* matrix, double mu, double sigma, int max_len){
  con->pdf = matrix;
  con->cdf = matrix + max_len;
  con->mu = mu;
  con->sigma = sigma;
  con->max_len = max_len;
}

void print_con(Connector* con){
  printf("Connector:\n |mu: %f\n |sigma: %f\n |max_len: %i\n\n", con->mu, con->sigma, con->max_len);
  if (con->max_len > -1){
    printf(" |Precomputed pdfs:\n");
    print_matrixf(con->pdf, 1, con->max_len);
    printf("\n |Precomputed cdfs:\n");
    print_matrixf(con->cdf, 1, con->max_len);
  }
}

void print_gap(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed, double num, double den){
  printf("Connector score for gap of %i\n", gap);
  printf("gap length: %i\n", gap);
  printf("effective length: %i\n", eff_len);
  printf("sequence length: %i\n", s_len);
  printf("number of recognizers: %i\n", n_recs);
  printf("is precomputed: %s\n", (precomputed) ? "True" : "False");
  printf("numerator: %f\n", num);
  printf("denominator %f\n", den);
}

double score_con(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed) {
  double num;
  double auc;
  double den;
  double res;
  if (precomputed){
    num = con->pdf[gap];
    auc = con->cdf[s_len - 1] - con->cdf[0];
    if (num < SMALL_POSITIVE) {
      num = SMALL_POSITIVE;
    }
    if (auc > SMALL_POSITIVE) {
      num /= auc;
    } else {
      num /= SMALL_POSITIVE;
    }
    den = 0.0;
    if (gap > -1 && gap <= eff_len)
      den = (DEN_EXPANSION[(eff_len - (gap + 1)) - 1] -
             DEN_EXPANSION[(n_recs - 1) - 1] -
             DEN_EXPANSION[eff_len - (gap + 1) - (n_recs - 1) - 1]) -
            (DEN_EXPANSION[eff_len - 1] -
             DEN_EXPANSION[n_recs - 1] -
             DEN_EXPANSION[eff_len - n_recs - 1]);
    res = log2f(num) - den;
    if (res < BIG_NEGATIVE)
      return BIG_NEGATIVE;
    return res;
  }



  num = get_numerator(s_len, gap, con->mu, con->sigma);
  den = get_denominator(gap + 1, n_recs, eff_len);

  return log2f(num/den);
}

