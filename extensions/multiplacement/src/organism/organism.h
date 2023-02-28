#ifndef ORGANISM_H
#define ORGANISM_H
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "../_aux.h"
#include "recognizer.h"
#include "connector.h"

typedef struct Organism Organism;

struct Organism {
  bool precomputed;
  int len;
  Recognizer* recs;
  Connector* cons;
};
void print_scores(Organism* org, double* r_scores, double* c_scores, int* c_lens, int n_rec);
void print_placement(Organism* org, const char* seq, int s_len, int* pos);
void parse_org(Organism *org, double* matrix, int* rec_lengths, double* models, double* edges, int* model_lengths, const char* rec_types, int num_recs, double* con_matrix, int max_len, bool precomputed);
void print_org(Organism *org);
void place_org(Organism* org, const char* seq,  int s_len, double* r_scores, double* c_scores, int* c_lens, bool precomputed);
#endif
