#include "connector.h"
void parse_con(Connector *con, double* matrix, double mu, double sigma, int max_len){
  //if (precomputed) {
    con->pdf = matrix;
    con->cdf = matrix + max_len;
    con->mu = mu;
    con->sigma = sigma;
    con->max_len = max_len;
    printf("%i\n", max_len);
  /*
  } else {
    con->pdf = NULL;
    con->cdf = NULL;
    con->mu = mu;
    con->sigma = sigma;
    con->max_len = 0;
  }
  */
}

double score_con(Connector* con, int gap, int s_len, int eff_len, int n_recs, bool precomputed) {
  if (precomputed){
    double num = con->pdf[gap];
    double auc = con->cdf[s_len - 1] - con->cdf[0];
    if (auc > SMALL_POSITIVE) {
      num /= auc;
    } else {
      num /= SMALL_POSITIVE;
    }

    double den = 0.0;
    if (gap > -1 && gap < eff_len){
    den = (DEN_EXPANSION[(eff_len - (gap + 1)) - 1] -
                 DEN_EXPANSION[(n_recs - 1) - 1] -
                 DEN_EXPANSION[eff_len - (gap + 1) - (n_recs) - 1]) -
                (DEN_EXPANSION[eff_len - 1] -
                 DEN_EXPANSION[n_recs - 1] -
                 DEN_EXPANSION[eff_len - n_recs - 1]);
    } else {
      den = BIG_NEGATIVE;
    }

    if (num - den > BIG_POSITIVE){
      printf("%f, %f\n\n", num, den);
      printf("%i, %i, %i, %i\n\n", s_len, eff_len, gap, n_recs);
    }
    if (num - den > BIG_NEGATIVE){
      return num - den;
    }else{
      return BIG_NEGATIVE;
    }
  }

  double num = get_numerator(s_len, gap, con->mu, con->sigma);
  double den = get_denominator(gap + 1, n_recs, eff_len);
  return log2f(num/den);
}

