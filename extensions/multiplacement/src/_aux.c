#include "_aux.h"
void print_matrixi(int* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%3i", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

void print_matrixf(float* matrix, int row, int col){
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      printf("|%8.1e", matrix[i * col + j]);
    }
    printf("|\n");
  }
}

int int_arr_sum(int* int_arr, int stop){
  int sum = 0; 
  for (int i = 0; i < stop; i++){
    sum += int_arr[i];
  }
  return sum;
}

int min(int n, int k){
  if (n > k)
    return k;
  else
    return n;
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

float norm_cdf(float x, float mu, float sigma){
  float z = (x - mu) / fabs(sigma);
  return (1 + erff(z / sqrtf(2.0))) / 2.00;
}

float norm_pf(float x, float mu, float sigma){
  if (sigma != 0)
    return norm_cdf(x + 0.5, mu, sigma) - norm_cdf(x - 0.5, mu, sigma);
  if (x == mu)
    return 1.00;
  return 0.00;
}

float get_numerator(int dna_length, int distance, float mu, float sigma){
  float numerator = norm_pf(distance, mu, sigma);
  if (sigma == 0.00)
    return numerator;

  float auc = norm_cdf(dna_length - 1, mu, sigma) - norm_cdf(0, mu, sigma);
  if (auc < 0.000001)
    auc = 0.000001;

  if (numerator < 0.00001)
    numerator = 0.00001;

  return numerator /= auc;

}

float get_denominator(int d, int N, int L){
  if (1 <= d && d <= L - N + 1)
    return (float) bin(L - d, N - 1) / (float) bin(L, N);
  return 0.0001;
}

float get_score(float *arr, int dna_length, int effective_length, int num_rec, 
                int gap_size, int max_length, int curr_conn, bool is_precomputed){
  if (is_precomputed == false){
    return log2f(get_numerator(dna_length, gap_size, arr[curr_conn * 2], arr[curr_conn * 2 + 1]) / 
                 get_denominator(gap_size + 1, num_rec, effective_length));
  }
  /*
  return arr[curr_conn * max_length + gap_size] - 
          (
            (
              NUMERATORS[(effective_length - (gap_size + 1)) - 1] - 
              NUMERATORS[(num_rec - 1) - 1] - 
              NUMERATORS[effective_length - (gap_size + 1) - (num_rec - 1) - 1]
            ) -

            (
              NUMERATORS[effective_length - 1] - 
              NUMERATORS[num_rec - 1] - 
              NUMERATORS[effective_length - num_rec - 1]
            )
          );
  */
}

int get_forward_offset(int index, int cols[], int num_rec) {
  // finds the first possible possition for a pssm
  // based on number of columns of preceding pssms

  int offset = 0;
  for (int i = 0; i < index; i++) {
    offset += cols[i];
  }
  return offset;
}

int get_reverse_offset(int index, int cols[], int num_rec) {
  // finds the last possible possition for a pssm
  // based on number of columns of subsequenct pssms

  int offset = 0;
  for (int i = num_rec - 1; i >= index; i--) {
    offset += cols[i];
  }

  // this is important to get the right number of possible alignments
  return offset - 1;
}

int max_index(float *arr, int size) {
  int max_index = 0;
  for (int i = 0; i < size; i++) {

    if (arr[i] > arr[max_index]) {
      max_index = i;
    }
  }
  return max_index;
}

float get_bin_frequency(float score, float bin_frequencies[], float bin_edges[], int num_bins){
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

float shape_average(float array[], int num_elements){
  float sum = array[0] + array[num_elements - 1]; 
  for (int i = 0; i < num_elements; i++){
    sum += array[i];
  }
  return sum/(num_elements + 2);
}
