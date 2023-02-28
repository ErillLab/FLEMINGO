#ifndef _AUX_H
#define _AUX_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include "_constants.h"

unsigned long long bin(unsigned long long n, unsigned long long k);
double norm_cdf(double x, double mu, double sigma);
double norm_pf(double x, double mu, double sigma);
double get_numerator(int dna_length, int distance, double mu, double sigma);
double get_denominator(int d, int N, int L);
int max_index(double *arr, int size);
double get_bin_frequency(double score, double bin_frequencies[], double bin_edges[], int num_bins);
double shape_average(double array[], int num_elements);
#endif
