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

double score_con(Connector* con, int gap, int s_len, int eff_len, int n_recs, double auc, bool precomputed) {
  double num;
  double den;
  double res;
  if (precomputed){
    num = con->pdf[gap];

    //sometimes the length of the sequence is greater than the max precomputed length
    //if this is the case then we need to compute the cdf to used for the placement

    if (num < BIG_NEGATIVE) {
      num = BIG_NEGATIVE;
    }
    num -= auc;
    den = 0.0;

    //if (gap + n_recs + 1) > eff_len then we get negative indices for the 3rd term
    //model probablility for this case should be very unlikely so I return a large
    //negative value
    if (gap > -1 && (gap + n_recs + 1) <= eff_len)
      den = (DEN_EXPANSION[(eff_len - (gap + 1)) - 1] -
             DEN_EXPANSION[(n_recs - 1) - 1] -
             DEN_EXPANSION[eff_len - (gap + 1) - (n_recs - 1) - 1]) -
            (DEN_EXPANSION[eff_len - 1] -
             DEN_EXPANSION[n_recs - 1] -
             DEN_EXPANSION[eff_len - n_recs - 1]);
    else
      return BIG_NEGATIVE;
    res = num - den;
    if (res < BIG_NEGATIVE)
      return BIG_NEGATIVE;
    if (res > BIG_POSITIVE){
      return BIG_POSITIVE;
    }
    return res;
  }


  num = get_numerator(s_len, gap, con->mu, con->sigma);
  den = get_denominator(gap + 1, n_recs, eff_len);
  return log2f(num/den);
}

