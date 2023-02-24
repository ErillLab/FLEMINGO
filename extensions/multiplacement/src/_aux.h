#ifndef _AUX_H
#define _AUX_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include "_constants.h"


int int_arr_sum(int*, int);

int min(int n, int k);

unsigned long long bin(unsigned long long n, unsigned long long k);

float norm_cdf(float x, float mu, float sigma);

float norm_pf(float x, float mu, float sigma);

float get_numerator(int dna_length, int distance, float mu, float sigma);

float get_denominator(int d, int N, int L);

float get_score(float *arr, int dna_length, int effective_length, int num_rec, int gap_size, int max_length, int curr_conn, bool is_precomputed);

int get_forward_offset(int index, int cols[], int num_rec);

int get_reverse_offset(int index, int cols[], int num_rec);

int max_index(float *arr, int size);

float get_bin_frequency(float score, float bin_frequencies[], float bin_edges[], int num_bins);

float shape_average(float array[], int num_elements);

#endif
