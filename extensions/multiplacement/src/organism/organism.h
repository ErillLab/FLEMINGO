#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "recognizer.h"
#include "connector.h"
#ifndef ORGANISM_H
#define ORGANISM_H

typedef struct Organism Organism;

struct Organism {
  bool precomputed;
  int len;
  Recognizer* recs;
  Connector* cons;
};

void parse_org(Organism *org, float* matrix, int* rec_lengths, float* models, float* edges, int* model_lengths, char* rec_types, int num_recs, float* con_matrix, int max_len, bool precomputed){
  org->len = num_recs;
  org->precomputed = precomputed;
  int m_offset = 0;
  int a_offset = 0;
  int c_offset = 0;
  int pssm_i = 0;
  int shape_i = 0;

  for (int i = 0; i < num_recs; i++){
    if (rec_types[i] == 'p') {
      parse_pssm(&org->recs[i], matrix + m_offset, rec_lengths[i]);
      m_offset += rec_lengths[pssm_i] * 4;
      pssm_i += 1;
    } else { 
      parse_shape(&org->recs[i], models + a_offset, models + a_offset + model_lengths[shape_i], edges + a_offset + shape_i, model_lengths[shape_i], rec_lengths[i], rec_types[i]);
      printf("null_model_offset: %i, alt_model_offset: %i\n", a_offset, a_offset + model_lengths[shape_i]);
      a_offset += model_lengths[shape_i] * 2;
      shape_i += 1;
    }

    if (i != num_recs - 1){
      if (precomputed) {
        parse_con(&org->cons[i], con_matrix + c_offset, 0.0, 0.0, max_len, precomputed);      
        c_offset += 2 * max_len;
      }else {
        parse_con(&org->cons[i], NULL, con_matrix[c_offset], con_matrix[c_offset + 1], max_len, precomputed);      
        c_offset += 2; 
      }
    }

  }
}

void print_org(Organism *org) {
  for (int i = 0; i < org->len; i++) {
    if (org->recs[i].feat == 'p') {
      for (int j = 0; j < 4; j++){
        for (int k = 0; k < org->recs[i].len; k++){
          printf("%3.2f ", org->recs[i].matrix[k * 4 + j]);
        }
        printf("\n");
      }
      printf("\n");
    } else {
      Recognizer* c = &org->recs[i];
      printf("Null frequencies:\n");
      for (int j = 0; j < org->recs[i].bin_s - 1; j++){
        printf("[%5.2f, %5.2f]: %3.2f\n", c->edges[j], c->edges[j+1], c->null_f[j]);
      }


      printf("Alt frequencies:\n");
      for (int j = 0; j < org->recs[i].bin_s - 1; j++){
        printf("[%5.2f, %5.2f]: %3.2f\n", c->edges[j], c->edges[j+1], c->alt_f[j]);
      }
    }
  }
  
}


void place_org( Organism* org,  char* seq,  int s_len, float* rec_scores, float* con_scores, int* con_lengths) {
  int n_recs = org->len;
  int* r_lens = (int*)malloc(n_recs * sizeof(int));
  int m_len = 0;
  for (int i = 0; i < n_recs; i++){
    r_lens[i] = org->recs[i].len;
    m_len += r_lens[i];
  }

  int eff_len = s_len - m_len - n_recs;
  int n_align = s_len - m_len + 1;

  int f_offset = 0;
  //int r_offset = s_len - m_len;

  float* t_row = (float*)malloc(n_align * sizeof(float));
  for (int i = 0; i < n_align; i++){
    t_row[i] = -INFINITY;
  }
  //memset(t_row, -100000.0, n_align * sizeof(float));
  float* c_row = (float*)calloc(n_align, sizeof(float));

  float* rs_matrix = (float*)calloc(n_align * n_recs, sizeof(float));
  float* gs_matrix = (float*)calloc(n_align * (n_recs - 1), sizeof(float));
  int*   tr_matrix = (int*)calloc(n_align * (n_recs - 1), sizeof(int));
  int    gap       = 0;
  float  g_score   = 0.0;

  Recognizer* rec = NULL;
  Connector* con = NULL;
  printf("top\n");
  for (int i = 0; i < n_recs; i++) {
    rec = &org->recs[i];
    printf("rec: %i of type %c, of len: %i\n", i, rec->feat, rec->len);
    //r_offset += r_lens[i];
    score_row(rec, seq + f_offset, n_align, rs_matrix + (i * n_align));
    if (i > 0) {
      con = &org->cons[i - 1];

      for (int j = 0; j < n_align; j++){
        for (int k = 0; k < j + 1; k++){
          gap = j - k;
          g_score = score_con(con, gap, s_len, eff_len, n_recs);
          if (t_row[j] < c_row[k] + g_score + rs_matrix[i * n_align + j]){
            t_row[j] = c_row[k] + g_score + rs_matrix[i * n_align + j];
            tr_matrix[(i - 1) * n_align + j] = gap;
            gs_matrix[(i - 1) * n_align + j] = g_score; 
          }

        }

      }

      for (int j = 0; j < n_align; j++){
        c_row[j] = t_row[j];
        t_row[j] = -INFINITY;
      }
    } else {
      for (int j = 0; j < n_align; j++){
        c_row[j] = rs_matrix[(i * n_align) + j];
      }
    }

    f_offset += r_lens[i];
  }

  int m_idx = max_index(c_row, n_align);
  rec_scores[n_recs] = c_row[m_idx];
  for (int i = n_recs - 1; i >= 0; i--) {
    rec_scores[i] = rs_matrix[i * n_align + m_idx];
    if (i > 0){
      con_scores[i - 1] = gs_matrix[(i - 1) * n_align + m_idx];
      con_lengths[i - 1] = tr_matrix[(i - 1) * n_align + m_idx];
      m_idx -= tr_matrix[(i - 1) * n_align + m_idx];
    }
  }
  printf("end\n");
  printf("lens\n");
  free(r_lens);
  printf("t_row\n");
  free(t_row);
  printf("c_row\n");
  free(c_row);
  printf("rs_matrix\n");
  free(rs_matrix);
  printf("gs_matrix\n");
  free(gs_matrix);
  printf("tr_matrix\n");
  free(tr_matrix);
}
#endif
