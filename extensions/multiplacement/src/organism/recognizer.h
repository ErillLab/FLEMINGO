#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../_constants.h"
#include "../_aux.h"
#ifndef RECOGNIZER_H
#define RECOGNIZER_H

typedef struct Recognizer Recognizer;
struct Recognizer {
  char feat;
  int bin_s;
  int len;
  float* matrix;
  float* alt_f;
  float* edges;
  float* null_f;
};

void parse_pssm(Recognizer *rec, float* matrix, int len){
  rec->matrix = matrix;
  rec->alt_f = NULL;
  rec->edges = NULL;
  rec->null_f = NULL;
  rec->bin_s = 0;
  rec->len = len;
  rec->feat = 'p';
}

void parse_shape(Recognizer* rec, float* null_f, float* alt_f, float* edges, int bin_s, int len, char feat) {
  rec->matrix = NULL;
  rec->null_f = null_f;
  rec->alt_f = alt_f;
  rec->edges = edges;
  rec->bin_s = bin_s;
  rec->len = len;
  rec->feat = feat;
}

void pssm_row( Recognizer* rec,  char* seq,  int len, float* row){
  float score = 0.0;
  float* matrix = rec->matrix;
  for (int i = 0; i < len; i++) {
    score = 0.0;
    for (int j = 0; j < rec->len; j++) {
      switch(seq[i + j]) {
        case 'a':
        case 'A':
          score += matrix[j * 4 + 0];
          break;
        
        case 'g':
        case 'G':
          score += matrix[j * 4 + 1];
          break;
        
        case 'c':
        case 'C':
          score += matrix[j * 4 + 2];
          break;
        
        case 't':
        case 'T':
          score += matrix[j * 4 + 3];
          break;
      }
    }
    row[i] = score;
  }

}

void mgw_row( Recognizer* rec,  char* seq,  int len, float* row){
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  float alt_f = 0.0;
  float null_f = 0.0;
  float* alt = rec->alt_f;
  float* null = rec->null_f;
  float* edges = rec->edges;
  float* pent_s = (float*)malloc(n_pent * sizeof(float));
  float score = 0.0;
  int idx = 0;
  for (int i = 0; i < len; i++){
    for (int j = 0; j < n_pent; j++){
      score = 0.0;
      idx = 0;
      for (int k = j; k < j + 5; k++){
        switch(seq[k]){
          case 'a':
          case 'A':
            idx += pow(4, 4 - (k - j)) * 0;
            break;

          case 'g':
          case 'G':
            idx += pow(4, 4 - (k - j)) * 1;
            break;

          case 'c':
          case 'C':
            idx += pow(4, 4 - (k - j)) * 2;
            break;

          case 't':
          case 'T':
            idx += pow(4, 4 - (k - j)) * 3;
            break;
        }
      }
      pent_s[j] = MGW_SCORES[idx];

    }
    score = shape_average(pent_s, n_pent);
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = log2f(alt_f / null_f);
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    if (score > BIG_POSITIVE)
      score = BIG_POSITIVE;

    row[i] = score;
  }
  free(pent_s);
}

void prot_row( Recognizer* rec,  char* seq,  int len, float* row){
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  float alt_f = 0.0;
  float null_f = 0.0;
  float* alt = rec->alt_f;
  float* null = rec->null_f;
  float* edges = rec->edges;
  float* pent_s = (float*)malloc(n_pent * sizeof(float));
  float score = 0.0;
  int idx = 0;

  for (int i = 0; i < len; i++){
    for (int j = 0; j < n_pent; j++){
      score = 0.0;
      idx = 0;
      for (int k = j; k < j + 5; k++){
        switch(seq[k]){
          case 'a':
          case 'A':
            idx += pow(4, 4 - (k - j)) * 0;
            break;

          case 'g':
          case 'G':
            idx += pow(4, 4 - (k - j)) * 1;
            break;

          case 'c':
          case 'C':
            idx += pow(4, 4 - (k - j)) * 2;
            break;

          case 't':
          case 'T':
            idx += pow(4, 4 - (k - j)) * 3;
            break;
        }
      }
      pent_s[j] = PROT_SCORES[idx];

    }
    score = shape_average(pent_s, n_pent);
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = log2f(alt_f / null_f);
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    if (score > BIG_POSITIVE)
      score = BIG_POSITIVE;

    row[i] = score;
  }
  free(pent_s);

}

void roll_row( Recognizer* rec,  char* seq,  int len, float* row){
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  float alt_f = 0.0;
  float null_f = 0.0;
  float* alt = rec->alt_f;
  float* null = rec->null_f;
  float* edges = rec->edges;
  float* pent_s = (float*)malloc((n_pent * 2) * sizeof(float));
  float score = 0.0;
  int idx = 0;

  for (int i = 0; i < len; i++){
    for (int j = 0; j < n_pent; j++){
      score = 0.0;
      idx = 0;
      for (int k = j; k < j + 5; k++){
        switch(seq[k]){
          case 'a':
          case 'A':
            idx += pow(4, 4 - (k - j)) * 0;
            break;

          case 'g':
          case 'G':
            idx += pow(4, 4 - (k - j)) * 1;
            break;

          case 'c':
          case 'C':
            idx += pow(4, 4 - (k - j)) * 2;
            break;

          case 't':
          case 'T':
            idx += pow(4, 4 - (k - j)) * 3;
            break;
        }
      }
      pent_s[j] = ROLL_SCORES[idx];
      pent_s[j + n_pent] = ROLL_SCORES[idx + 1024];

    }
    score = shape_average(pent_s, n_pent);
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = log2f(alt_f / null_f);
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    if (score > BIG_POSITIVE)
      score = BIG_POSITIVE;

    row[i] = score;
  }
  free(pent_s);

}

void helt_row( Recognizer* rec,  char* seq,  int len, float* row){
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  float alt_f = 0.0;
  float null_f = 0.0;
  float* alt = rec->alt_f;
  float* null = rec->null_f;
  float* edges = rec->edges;
  float* pent_s = (float*)malloc((n_pent * 2) * sizeof(float));
  float score = 0.0;
  int idx = 0;

  for (int i = 0; i < len; i++){
    for (int j = 0; j < n_pent; j++){
      score = 0.0;
      idx = 0;
      for (int k = j; k < j + 5; k++){
        switch(seq[k]){
          case 'a':
          case 'A':
            idx += pow(4, 4 - (k - j)) * 0;
            break;

          case 'g':
          case 'G':
            idx += pow(4, 4 - (k - j)) * 1;
            break;

          case 'c':
          case 'C':
            idx += pow(4, 4 - (k - j)) * 2;
            break;

          case 't':
          case 'T':
            idx += pow(4, 4 - (k - j)) * 3;
            break;
        }
      }
      pent_s[j] = HELT_SCORES[idx];
      pent_s[j + n_pent] = HELT_SCORES[idx + 1024];

    }
    score = shape_average(pent_s, n_pent);
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = log2f(alt_f / null_f);
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    if (score > BIG_POSITIVE)
      score = BIG_POSITIVE;

    row[i] = score;
  }
  free(pent_s);

}

void shape_row( Recognizer* rec,  char* seq,  int len, float* row,  char feat){
  switch(feat) {
    case 'm':
    case 'M':
      mgw_row(rec, seq, len, row);
      break;

    case 't':
    case 'T':
      prot_row(rec, seq, len, row);
      break;

    case 'r':
    case 'R':
      roll_row(rec, seq, len, row);
      break;

    case 'h':
    case 'H':
      helt_row(rec, seq, len, row);
      break;
  }

}

void score_row( Recognizer* rec,  char* seq,  int len, float* row){
  if (rec->feat == 'p'){
    pssm_row(rec, seq, len, row);
  }else{
    shape_row(rec, seq, len, row, rec->feat);
  }
}

#endif