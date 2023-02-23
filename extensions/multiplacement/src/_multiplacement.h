#ifndef _MULTIPLACEMENT_H
#define _MULTIPLACEMENT_H
/*
#define PY_SSIZE_T_CLEAN
#include <python3.11/Python.h>
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "_aux.h"

void traceback(int num_rec, int len_seq, float* con_matrices, float* rec_score_matrix, int num_alignments, float* rec_alignments, int* con_alignments, float* rec_scores, float* con_scores, int* con_lengths, int max_seq_len, int effective_len, bool is_precomputed);

void fill_traceback_matrix(float *score_matrix, int num_alignments, float *gapMatrix, int *cols, int num_rec, int len_seq, float *con_scores, float *rec_scores, int *con_lengths, int max_length, bool is_precomputed);

void fill_matrix(const char seq[], int len_seq, float pssm[], int cols[], const char rec_types[], int num_rec, float score_matrix[], int num_alignments, float bin_freqs[], float bin_edges[], int num_bins[]);

void fill_row_pssm(const char* seq, int len_seq, int row, int curr_rec, float rec_matrices[], int pssm_score_offset, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[]);

void fill_row_mgw(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins);

void fill_row_prot(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins);

void fill_row_helt(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins);

void fill_row_roll(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins);

void place(char* seq, int s_len, float*** recs, int* r_lens, char* r_types, int n_rec, int* b_lens, float** cons, int m_len, float* r_scores, float* c_scores, int* c_lens);
#endif
