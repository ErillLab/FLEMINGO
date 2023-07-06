#include "recognizer.h"
/*
name:            parse_pssm
pre-conditions:  rec is not null, matrix is not null, len > 0
post-conditions: recognizer has pssm attributes set
parameters:      rec:    memory address of recognizer
                 matrix: pointer to the start of the recognizers pssm
                 len:    length of the pssm
notes:           sets the attributes of the pssm to whatever is passed
*/
void parse_pssm(Recognizer *rec, double* matrix, int len){
  rec->matrix = matrix;
  rec->alt_f = NULL;
  rec->edges = NULL;
  rec->null_f = NULL;
  rec->bin_s = 0;
  rec->len = len;
  rec->feat = 'p';
}

/*
name:            parse_shape
pre-conditions:  rec is not null, null_f is not null, alt_f is not null
                 edges is not null, bin_s is > 0, len > 4, feat is 
                 {m, t, r, or h}. number of elements in null_f and alt_f 
                 are both exactly 1 less than the number of elements in edges.
post-conditions: recognizer has shape attributes set
parameters:      rec:    memory address of recognizer struct
                 null_f: array of doubles containing log2 of frequencies in alt
                         model
                 alt_f:  array of doubles containing log2 of frequencies in null
                         model
                 edges:  array of doubles containing interval bounds for each
                         interval for the models
                 bin_s:  the number of bin edges for the models
                 len:    integer length of the shape recognizer
                 feat:   char representing dna shape feature
notes:           .
*/
void parse_shape(Recognizer* rec, double* null_f, double* alt_f, double* edges, int bin_s, int len, char feat) {
  rec->matrix = NULL;
  rec->null_f = null_f;
  rec->alt_f = alt_f;
  rec->edges = edges;
  rec->bin_s = bin_s;
  rec->len = len;
  rec->feat = feat;
}

/*
name:            print_rec
pre-conditions:  recognizer is populated with pssm or shape information
post-conditions: none
parameters:      rec: memory address of recognizer struct
notes:           calls the correct display function based on recognizer type
*/
void print_rec(Recognizer *rec){
  printf("Recognizer Feature: %c\n", rec->feat);
  printf("Recognizer Length %i\n", rec->len);
  if (rec->feat != 'p'){
    print_shape(rec);
  }else{
    print_pssm(rec);
  }
}

/*
name:            print_pssm
pre-conditions:  recognizer struct is populated with pssm information
post-conditions: none
parameters:      rec: memory address of the recognizer struct
notes:           prints all of the values in the pssm (rows are bases)
*/
void print_pssm(Recognizer *rec){
  for (int j = 0; j < 4; j++){
    if (j == 0)
      printf(" A");
    if (j == 1)
      printf(" G");
    if (j == 2)
      printf(" C");
    if (j == 3)
      printf(" T");
    for (int k = 0; k < rec->len; k++){
      printf("%5.2f ", rec->matrix[k * 4 + j]);
    }
    printf("\n");
  }
}

/*
name:            print_shape
pre-conditions:  recognizer struct is populated with shape information
post-conditions: none
parameters:      rec: memory address of the recognizer struct
notes:           prints all of the frequencies (log2 of frequenceies)
                 and the intervals
*/
void print_shape(Recognizer *rec){
  
  printf("Null frequencies:\n");
  for (int j = 0; j < rec->bin_s - 1; j++){
    printf(" [%5.2f, %5.2f]: %8.5f\n", rec->edges[j], rec->edges[j+1], rec->null_f[j]);
  }


  printf("Alt frequencies:\n");
  for (int j = 0; j < rec->bin_s - 1; j++){
    printf(" [%5.2f, %5.2f]: %8.5f\n", rec->edges[j], rec->edges[j+1], rec->alt_f[j]);
  }
  
}

/*
name:            pssm_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placement 
notes:           .
*/
void pssm_row( Recognizer* rec,  const char* seq,  int len, double* row){
  // score:  calculated score a single placement of the recognizer
  // matrix: pointer to the pssm of the recognizer
  double score = 0.0;
  double* matrix = rec->matrix;

  // for each possible starting position, i, we record the scores into row
  for (int i = 0; i < len; i++) {
    score = 0.0;

    // match the base for each column of the recognizer and add the
    // corresponding score to the current score
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

/*
name:            mgw_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placements

notes:           .
*/
void mgw_row( Recognizer* rec,  const char* seq,  int len, double* row){
  // n_pent: number of pentamers that the recognizer is placed on for a single 
  //         placement
  // bin_s:  the number of bin edges for the recogninzers models
  // alt_s:  frequency obtained based on which interval the placement score falls
  //         in for the alternative model
  // null_s: frequency obtained based on which interval the placement score falls
  //         in for the null model
  // edges:  boundaries of each interval for the models
  // pent_s: array to hold the placement scores for each pentamer that the
  //         recognizer is placed on for a single placement.
  // score:  score for the placement, obtained by averaging pentamer scores
  // idx:    index in the constant array of scores based on pentamer sequence
  int n_pent = rec->len - 4;
  double alt_s = 0.0;
  double null_s = 0.0;
  double* pent_s = (double*)malloc(n_pent * sizeof(double));
  double score = 0.0;
  int idx = 0;

  // for each possible starting position, i, we record the scores into row
  for (int i = 0; i < len; i++){

    // we have to obtain the score for each pentamer contained in the placement
    // of the recognizer, so we find the appropriate index based on the sequence
    // for each pentamer j
    for (int j = i; j < i + n_pent; j++){
      score = 0.0;
      idx = 0;

      // index the pentamer starting at position k, until position k + 5
      // indexing is done in base four, where A = 0, G = 1, C = 2, T = 3
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
      pent_s[j - i] = MGW_SCORES[idx];
    }
    score = shape_average_mgw_prot(pent_s, n_pent);

    // frequencies are obtained by observing which interval the score falls in
    // for the alternative and null models
    // NOTE, the logs of the frequencies are actually what is stored, so we just
    // subtract the obtained "frequencies" to obtain our log-liklihood ratio
    alt_s = get_bin_frequency(score, rec->alt_f, rec->edges, rec->bin_s);
    null_s = get_bin_frequency(score, rec->null_f, rec->edges, rec->bin_s);
    score = alt_s - null_s;
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    // this can only occur if a problem occurs when finding where the score
    // belongs in the null model, se the program exits.
    if (isnan(null_s)){
      printf("Something went wrong, impossible score selected from null model...\n");
      print_rec(rec);
      exit(1);
    }
    row[i] = score;
  }
  free(pent_s);
}

/*
name:            prot_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placement
notes:           .
*/
void prot_row( Recognizer* rec,  const char* seq,  int len, double* row){
  // n_pent: number of pentamers that the recognizer is placed on for a single 
  //         placement
  // n_bins: the number of bin edges for the recogninzers models
  // alt_f:  frequency obtained based on which interval the placement score falls
  //         in for the alternative model
  // null_f: frequency obtained based on which interval the placement score falls
  //         in for the null model
  // edges:  boundaries of each interval for the models
  // pent_s: array to hold the placement scores for each pentamer that the
  //         recognizer is placed on for a single placement.
  // score:  score for the placement, obtained by averaging pentamer scores
  // idx:    index in the constant array of scores based on pentamer sequence
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  double alt_f = 0.0;
  double null_f = 0.0;
  double* alt = rec->alt_f;
  double* null = rec->null_f;
  double* edges = rec->edges;
  double* pent_s = (double*)malloc(n_pent * sizeof(double));
  double score = 0.0;
  int idx = 0;

  // for each possible starting position, i, we record the scores into row
  for (int i = 0; i < len; i++){

    // we have to obtain the score for each pentamer contained in the placement
    // of the recognizer, so we find the appropriate index based on the sequence
    // for each pentamer j
    for (int j = i; j < i + n_pent; j++){
      score = 0.0;
      idx = 0;

      // index the pentamer starting at position k, until position k + 5
      // indexing is done in base four, where A = 0, G = 1, C = 2, T = 3
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
      pent_s[j - i] = PROT_SCORES[idx];

    }
    score = shape_average_mgw_prot(pent_s, n_pent);

    // frequencies are obtained by observing which interval the score falls in
    // for the alternative and null models
    // NOTE, the logs of the frequencies are actually what is stored, so we just
    // subtract the obtained "frequencies" to obtain our log-liklihood ratio
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = alt_f - null_f;
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    // this can only occur if a problem occurs when finding where the score
    // belongs in the null model, se the program exits.
    if (isnan(null_f)){
      printf("Something went wrong, impossible score selected from null model...\n");
      print_rec(rec);
      exit(1);
    }
    row[i] = score;
  }
  free(pent_s);

}


/*
name:            roll_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placements
notes:           .
*/
void roll_row( Recognizer* rec,  const char* seq,  int len, double* row){
  // n_pent: number of pentamers that the recognizer is placed on for a single 
  //         placement
  // n_bins: the number of bin edges for the recogninzers models
  // alt_f:  frequency obtained based on which interval the placement score falls
  //         in for the alternative model
  // null_f: frequency obtained based on which interval the placement score falls
  //         in for the null model
  // edges:  boundaries of each interval for the models
  // pent_s: array to hold the placement scores for each pentamer that the
  //         recognizer is placed on for a single placement. roll is a step
  //         based shape feature so it actually has 2 * n_pent indices
  // score:  score for the placement, obtained by averaging pentamer scores
  // idx:    index in the constant array of scores based on pentamer sequence
  int n_pent = rec->len - 4;
  int n_bins = rec->bin_s;
  double alt_f = 0.0;
  double null_f = 0.0;
  double* alt = rec->alt_f;
  double* null = rec->null_f;
  double* edges = rec->edges;
  double* pent_s = (double*)malloc((n_pent * 2) * sizeof(double));
  double score = 0.0;
  int idx = 0;

  // for each possible starting position, i, we record the scores into row
  for (int i = 0; i < len; i++){

    // we have to obtain the score for each pentamer contained in the placement
    // of the recognizer, so we find the appropriate index based on the sequence
    // for each pentamer j
    for (int j = i; j < i + n_pent; j++){
      score = 0.0;
      idx = 0;

      // index the pentamer starting at position k, until position k + 5
      // indexing is done in base four, where A = 0, G = 1, C = 2, T = 3
      for (int k = j; k < (j + 5); k++){
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
      pent_s[j - i] = ROLL_SCORES[idx];
      pent_s[j - i + n_pent] = ROLL_SCORES[idx + 1024];
    }
    score = shape_average_helt_roll(pent_s, n_pent * 2);

    // frequencies are obtained by observing which interval the score falls in
    // for the alternative and null models
    // NOTE, the logs of the frequencies are actually what is stored, so we just
    // subtract the obtained "frequencies" to obtain our log-liklihood ratio
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = alt_f - null_f;
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    // this can only occur if a problem occurs when finding where the score
    // belongs in the null model, se the program exits.
    if (isnan(null_f)){
      printf("Something went wrong, impossible score selected from null model...\n");
      print_rec(rec);
      exit(1);
    }
    row[i] = score;
  }
  free(pent_s);

}

/*
name:            helt_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placements

notes:           .
*/
void helt_row( Recognizer* rec,  const char* seq,  int len, double* row){

  // n_pent: number of pentamers that the recognizer is placed on for a single 
  //         placement
  // n_bins: the number of bin edges for the recogninzers models
  // alt_f:  frequency obtained based on which interval the placement score falls
  //         in for the alternative model
  // null_f: frequency obtained based on which interval the placement score falls
  //         in for the null model
  // edges:  boundaries of each interval for the models
  // pent_s: array to hold the placement scores for each pentamer that the
  //         recognizer is placed on for a single placement. helt is a step
  //         based shape feature so it actually has 2 * n_pent indices
  // score:  score for the placement, obtained by averaging pentamer scores
  // idx:    index in the constant array of scores based on pentamer sequence
  int n_pent = rec->len - 4; 
  int n_bins = rec->bin_s;     
  double alt_f = 0.0;
  double null_f = 0.0;
  double* alt = rec->alt_f;
  double* null = rec->null_f;
  double* edges = rec->edges;
  double* pent_s = (double*)malloc((n_pent * 2) * sizeof(double));
  double score = 0.0;
  int idx = 0;

  // for each possible starting position, i, we record the scores into row
  for (int i = 0; i < len; i++){

    // we have to obtain the score for each pentamer contained in the placement
    // of the recognizer, so we find the appropriate index based on the sequence
    // for each pentamer j
    for (int j = i; j < i + n_pent; j++){
      score = 0.0;
      idx = 0;

      // index the pentamer starting at position k, until position k + 5
      // indexing is done in base four, where A = 0, G = 1, C = 2, T = 3
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

      // helt features are step based and have 2 scores per pentamer, so both
      // are recorded
      pent_s[j - i] = HELT_SCORES[idx];
      pent_s[j - i + n_pent] = HELT_SCORES[idx + 1024];
    }

    // the scores for step based shape features are an average of all recorded
    // scores, plus the first and last score again, and divided by the number of
    // elements in the pent_s array + 2
    score = shape_average_helt_roll(pent_s, n_pent * 2);

    // frequencies are obtained by observing which interval the score falls in
    // for the alternative and null models
    // NOTE, the logs of the frequencies are actually what is stored, so we just
    // subtract the obtained "frequencies" to obtain our log-liklihood ratio
    alt_f = get_bin_frequency(score, alt, edges, n_bins);
    null_f = get_bin_frequency(score, null, edges, n_bins);
    score = alt_f - null_f;
    if (score < BIG_NEGATIVE)
      score = BIG_NEGATIVE;

    // this can only occur if a problem occurs when finding where the score
    // belongs in the null model, se the program exits.
    if (isnan(null_f)){
      printf("Something went wrong, impossible score selected from null model...\n");
      print_rec(rec);
      exit(1);
    }

    // score is recorded into row at index i
    row[i] = score;
  }
  free(pent_s);

}

/*
name:            shape_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placementsname:            

notes:           calls appropriate function for scoring a row based on shape
                 recognizer type
*/
void shape_row( Recognizer* rec,  const char* seq,  int len, double* row){
  switch(rec->feat) {
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

/*
name:            score_row

pre-conditions:  placement has been initiated, rec is populated, seq is composed
                 of only chars A, a, G, g, C, c, T or t, row is allocated space
                 to record placement scores of the recognizer

post-conditions: row is populated with the scores of all possible placements of
                 the recognizer along the sequence

parameters:      rec: pointer to recognizer struct in the organism
                 seq: DNA sequence, starting at first possible placement of rec
                 len: number of placements of rec along the dna sequence
                 row: allocated memory for scores of placements

notes:           calls the appropriate function depending on recognizer type
*/
void score_row( Recognizer* rec,  const char* seq,  int len, double* row){
  if (rec->feat == 'p' || rec->feat == 'P'){
    pssm_row(rec, seq, len, row);
  }else{
    shape_row(rec, seq, len, row);
  }
}