#include "_aux.h"
void print_matrixi(int* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%3i", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

void print_matrixf(double* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%8.1e", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

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

double norm_cdf(double x, double mu, double sigma){
  double z = (x - mu) / fabs(sigma);
  return (1 + erff(z / sqrtf(2.0))) / 2.00;
}

double norm_pf(double x, double mu, double sigma){
  if (sigma != 0)
    return norm_cdf(x + 0.5, mu, sigma) - norm_cdf(x - 0.5, mu, sigma);
  if (x == mu)
    return 1.00;
  return 0.00;
}

double get_numerator(int dna_length, int distance, double mu, double sigma){
  double numerator = norm_pf(distance, mu, sigma);
  if (sigma == 0.00)
    return numerator;

  double auc = norm_cdf(dna_length - 1, mu, sigma) - norm_cdf(0, mu, sigma);
  if (auc < 0.000001)
    auc = 0.000001;

  if (numerator < 0.00001)
    numerator = 0.00001;

  return numerator /= auc;

}

double get_denominator(int d, int N, int L){
  if (1 <= d && d <= L - N + 1)
    return (double) bin(L - d, N - 1) / (double) bin(L, N);
  return 0.0001;
}

int max_index(double *arr, int size) {
  int max_index = 0;
  for (int i = 0; i < size; i++) {

    if (arr[i] > arr[max_index]) {
      max_index = i;
    }
  }
  return max_index;
}

double get_bin_frequency(double score, double bin_frequencies[], double bin_edges[], int num_bins){
  for (int i = 1; i < num_bins; i++){
    if (score < bin_edges[i]){
      if (bin_frequencies[num_bins - 1] == 0.00) {
        return bin_frequencies[i];
      }else{
        return bin_frequencies[i - 1];
      }
    }
  } 

  return bin_frequencies[num_bins - 1];
}

double shape_average(double array[], int num_elements){
  double sum = array[0] + array[num_elements - 1]; 
  for (int i = 0; i < num_elements; i++){
    sum += array[i];
  }
  return sum/(num_elements + 2);
}
