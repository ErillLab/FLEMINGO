#include "_aux.h"

/* 
name:            print_matrixi
pre-conditions:  populated matrix of integers
post-conditions: .
parameters:      matrix: array of integers
                 row:    number of rows
                 col:    number columns
notes:           .
*/
void print_matrixi(int* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%3i", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

/* 
name:            print_matrixf
pre-conditions:  populated matrix of doubles
post-conditions: .
parameters:      matrix: array of doubles
                 row:    number of rows
                 col:    number of columns
notes:           .
*/
void print_matrixf(double* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%8.1e", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

/* 
name:            bin
pre-conditions:  n and k are positive integers
post-conditions: returns the binomial coefficient of n and k
parameters:      n: some positive integer
                 k: some positive integer
notes:           overflow can happen, but it would be pretty unlikely
*/
unsigned long long bin(unsigned long long n, unsigned long long k){
  if (n == k)
    return 1;
  unsigned long long c = 1;
  for (unsigned long long i = 1; i <= k; i++, n--) {

    if (c/i > ULONG_MAX/n) // return 0 on potential overflow
      return 0;

    c = c / i * n + c % i * n / i; // split c * n / i into (c / i * i + c % i) * n / i
  }

  return c; 
}

/* 
name:            norm_cdf
pre-conditions:  mu and sigma describe a normal distribution, x is a value
post-conditions: returns the cumulitive distribution around x in a 
                 normal distribution
parameters:      x:     a value in a normal distribution
                 mu:    the mean of the normal distribution
                 sigma: the standard deviation of the normal distribution
notes:           .
*/
double norm_cdf(double x, double mu, double sigma){
  double z = (x - mu) / fabs(sigma);
  return (1 + erff(z / sqrtf(2.00))) / 2.00;
}

/* 
name:            norm_pf
pre-conditions:  mu and sigma describe a normal distribution, x is a value
post-conditions: returns the probability function for a value x
                 returns 0.00 if sigma if 
parameters:      x:     a value in the normal distribution
                 mu:    the mean of the normal distribution
                 sigma: the standard deviation of the normal distribution
notes:           .
*/
double norm_pf(double x, double mu, double sigma){
  if (sigma != 0)
    return norm_cdf(x + 0.5, mu, sigma) - norm_cdf(x - 0.5, mu, sigma);
  if (x == mu)
    return 1.00;
  return 0.00;
}

/* 
name:            get_numerator
pre-conditions:  dna_length > 0, distance > -1, mu is real, sigma > 0.00
post-conditions: returns the numerator value for computing the probability
                 given the model
parameters:      dna_length: the length of the DNA
                 distance:   the length of the gap
                 mu:         the average of the normal distribution
                             of the connector
                 sigma:      the standard deviation of the nomral distribution
                             of the connector
notes:           .
*/
double get_numerator(int dna_length, int distance, double mu, double sigma){
  double numerator = norm_pf(distance, mu, sigma);
  if (sigma == 0.00)
    return numerator;

  double auc = norm_cdf(dna_length - 1, mu, sigma) - norm_cdf(0, mu, sigma);
  if (auc < SMALL_POSITIVE)
    auc = 1E-10;

  if (numerator < SMALL_POSITIVE)
    numerator = 1E-10;

  return (numerator / auc);

}

/* 
name:            get_denominator
pre-conditions:  attributes about the gap have been calculated
post-conditions: returns the denominator for calculating the
                 probability given d, N, and L
parameters:      d: the length of the gap plus 1
                 N: the number of recognizers
                 L: the length of the DNA strand minus the 
                    sum of the lengths of the recognizers plus 1
notes:           .
*/
double get_denominator(int d, int N, int L){
  if (1 <= d && d <= L - N + 1)
    return (double) bin(L - d, N - 1) / (double) bin(L, N);
  return SMALL_POSITIVE;
}

/* 
name:            max_index
pre-conditions:  array is populated with doubles, size is > 0
post-conditions: returns the index of the max element in arr
parameters:      arr:  array of doubles
                 size: the number of elements in arr
notes:           .
*/
int max_index(double *arr, int size) {
  int max_index = 0;
  for (int i = 0; i < size; i++) {

    if (arr[i] > arr[max_index]) {
      max_index = i;
    }
  }
  return max_index;
}



bool compare_doubles(double a, double b, double epsilon){
  if (fabs(a - b) < epsilon){
    return true;
  }

  return false;
}
/* 
name:            get_bin_frequency
pre-conditions:  score is calcualted from constant arrays of of shape scores
                 described in _constants.h
                 bin frequencies and edges are populated with appropriate 
                 information for null and alternative models
post-conditions: returns the index of the bin that the calculated score fits
                 in the range of 
parameters:      score:           the calculated score from sliding the pentamers
                                  over the sequence
                 bin_frequencies: the frequencies of the scores in each bin
                 bin_edges:       the range of each bin
                 num_bins:        the number of bins for the model
notes:           .
*/
double get_bin_frequency(double score, double bin_frequencies[], double bin_edges[], int num_bins){
  for (int i = 1; i < num_bins; i++){
    if (score < bin_edges[i] && fabs(score - bin_edges[i]) > 0.001){
      //printf("bin edge %i is %f, score is: %f, freq is: %f\n", i, bin_edges[i], score, bin_frequencies[i - 1]);
      return bin_frequencies[i - 1];
    }
  } 

  return bin_frequencies[num_bins - 1];
}

/* 
name:            shape_average
pre-conditions:  array of pentamer scores is populated
post-conditions: returns the average of all elements in the array
parameters:      pentamer_scores: the scores for each pentamer in the
                                  shape recognizer
                 num_elements:    the number of pentamers in the shape
                                  recognizer
notes:           the terminal elements are counted twice to weight
                 terminal scores for Roll and HelT the same as internal ones
*/
double shape_average(double pentamer_scores[], int num_elements){
  double sum = pentamer_scores[0] + pentamer_scores[num_elements - 1]; 
  for (int i = 0; i < num_elements; i++){
    sum += pentamer_scores[i];
  }
  return sum/((double)(num_elements + 2));
}

double shape_average_mgw_prot(double pentamer_scores[], int num_elements){
  double sum = 0.00;
  for (int i = 0; i < num_elements; i++){
    sum += pentamer_scores[i];
  }
  return sum/((double)(num_elements));  
  
}