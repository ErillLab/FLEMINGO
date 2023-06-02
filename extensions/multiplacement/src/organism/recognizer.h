#ifndef RECOGNIZER_H
#define RECOGNIZER_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../_constants.h"
#include "../_aux.h"

typedef struct Recognizer Recognizer;
struct Recognizer {
  char feat;
  int bin_s;
  int len;
  double* matrix;
  double* alt_f;
  double* edges;
  double* null_f;
};

void print_rec(Recognizer* rec);
void parse_pssm(Recognizer* rec, double* matrix, int len);
void print_pssm(Recognizer* rec);
void parse_shape(Recognizer* rec, double* null_f, double* alt_f, double* edges, int bin_s, int len, char feat) ;
void print_shape(Recognizer* rec);
void pssm_row( Recognizer* rec,  const char* seq,  int len, double* row);
void mgw_row( Recognizer* rec,  const char* seq,  int len, double* row);
void prot_row( Recognizer* rec,  const char* seq,  int len, double* row);
void roll_row( Recognizer* rec,  const char* seq,  int len, double* row);
void helt_row( Recognizer* rec,  const char* seq,  int len, double* row); 
void shape_row( Recognizer* rec,  const char* seq,  int len, double* row);
void score_row( Recognizer* rec,  const char* seq,  int len, double* row);

#endif