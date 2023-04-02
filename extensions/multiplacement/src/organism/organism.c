#include "organism.h"

/*
name:            print_scores
pre-conditions:  organism is populated organism information
                 r_scores is populated with computed placement scores
                 c_scores is populated with computed placement scores
                 c_lens is computed with placement lengths
post-conditions: .
parameters:      org:      organism with recognizers and connector arrays populated
                 r_scores: scores for each recognzier from place_org
                 c_scores: scores for each connector from place_org
                 c_lens:   length of each gap 
notes:           .
*/
void print_scores(Organism* org, double* r_scores, double* c_scores, int* c_lens){
  int n_rec = org->len;
  for (int i = 0; i < n_rec; i++){
    printf("|%c[%5.2f]|", org->recs[i].feat, r_scores[i]);
    if (i < n_rec - 1){
      printf("+|(%i)", c_lens[i + 1]);
      printf("c[%5.2f]|+", c_scores[i]);
    }
  }
  printf("=%5.2f\n", r_scores[n_rec]);

}

/*
name:            print_placement
pre-conditions:  organism is populated organism information
                 positions for placement have been computed in place_org
post-conditions: .
parameters:      org:   organism with recognziers and connector arrays popualted
                 seq:   are of characters representing DNA sequence
                 s_len: the length of seq
                 pos:   the positions of each recognizer placement
notes:           .
*/
void print_placement(Organism* org, const char* seq, int s_len, int* pos){
  int n_rec = org->len; 
  for (int i = 0; i < s_len; i++){
    printf("%c", seq[i]);
  }
  printf("\n");

  for (int i = 0; i < pos[0]; i++){
    printf(" ");
  }

  int s = 0;
  int e = 0;
  for (int i = 0; i < n_rec; i++){
    s = e + pos[i];
    e = s + org->recs[i].len;

    printf("%c", org->recs[i].feat);
    for (int j = s + 1; j < e - 1; j++)
      printf("#");

    if (org->recs[i].len > 1)
      printf("%c", org->recs[i].feat);

    if (i < n_rec - 1)
      for (int j = e; j < e + pos[i + 1]; j++)
        printf("-");
  }

  for (int i = pos[n_rec - 1] + org->recs[n_rec - 1].len; i < s_len; i++){
    printf(" ");
  }
  printf("\n\n");

}

/*
name:            parse_org
pre-conditions:  all information for recognizers and connectors has been loaded
                 number of recognziers is 1 greater than number of recognizers
post-conditions: organism has been loaded with parsed information
parameters:      org:           pointer to organism so it can be modified
                                in this function
                 matrix:        1D array concatenation of all PSSM scores in 
                                column first ordering
                 rec_lengths:   array of length of each recognizer
                 models:        1D array concatenation of null and alternative
                                model frequencies for all shape recognizers
                 edges:         1D array concatenation of the lower bound of 
                                each range for all bins of all shape 
                                recognizers
                 model_lengths: an array containing the lengths of all models
                                of all shape recognzers
                 rec_types:     an array of chars denoting the type of all
                                recognziers
                 num_recs:      the total number of recognizers
                 con_matrix:    1D array concatenation of mu, sigma,
                                pdfs, and cdfs, for all connectors
                 max_len:       the max length that has been precomputed up to
                                for connectors
notes:           .
*/
void parse_org(Organism *org, double* matrix, int* rec_lengths, double* models,
               double* edges, int* model_lengths, const char* rec_types, 
               int num_recs, double* con_matrix, int max_len){

  org->len  = num_recs;
  org->recs = (Recognizer*)malloc(num_recs * sizeof(Recognizer));
  org->cons = (Connector*)malloc((num_recs - 1) * sizeof(Connector));

  //m_offset: keeps track of the offset in the 1D concatenation of PSSM scores
  //a_offset: keeps track of the offset in the 1D concatenation of model
  //          frequencies
  //c_offset: keeps track of the offset in the 1D concatenation of mu's,
  //          sigma's, pdf's, and cdf's
  //e_offset: keeps track of the offset in the 1D concatenation of lower bounds
  //          of each bin for each model
  int m_offset = 0;
  int a_offset = 0;
  int c_offset = 0;
  int e_offset = 0;
  int shape_i  = 0;

  for (int i = 0; i < num_recs; i++){

    //parse the information for recognizers.
    //pssms are passed the matrix and the offset is update.
    //shape recognizers are passed null models, alternative models,
    //edges, length of their models, and the type of recognizer
    if (rec_types[i] == 'p') {
      parse_pssm(&org->recs[i], 
                 matrix + m_offset, 
                 rec_lengths[i]);
      m_offset += rec_lengths[i] * 4;
    } else { 
      parse_shape(&org->recs[i], 
                  models + a_offset, 
                  models + a_offset + model_lengths[shape_i], 
                  edges + e_offset, 
                  model_lengths[shape_i], 
                  rec_lengths[i], 
                  rec_types[i]);
      a_offset += model_lengths[shape_i] * 2;
      e_offset += model_lengths[shape_i];
      shape_i  += 1;
    }

    //connectors are passed pdfs and cdfs,
    //mu, sigma, and the length that is precomputed up until
    if (i < num_recs - 1){
        parse_con(&org->cons[i], 
                  con_matrix + 2 + c_offset, 
                  *(con_matrix + c_offset), 
                  *(con_matrix + c_offset + 1), 
                  max_len);      
        c_offset += 2 * (max_len + 1) + 2;
    }

  }

}

/*
name:            print_org
pre-conditions:  org has been populated
post-conditions: .
parameters:      org: organism populated with recognziers and connectors
notes:           .
*/
void print_org(Organism *org) {
  bool has_shape = false;
  for (int i = 0; i < org->len; i++){
    if (org->recs[i].feat != 'p')
      has_shape = true;
  }

  if (has_shape == false)
    return;
  printf("*************************************************\n");
  for (int i = 0; i < org->len; i++){
    print_rec(&org->recs[i]);
    printf("\n");
    if (i < org->len - 1){
      print_con(&org->cons[i]);
      printf("\n");
    }
  }
  printf("************************************************\n");
}

/*
name:            place_org
pre-conditions:  org has been populated
post-conditions: score arrays for recognizers and connectors
                 have been filled
parameters:      org:      a populated organism
                 seq:      an array of characters representing the dna sequence
                 s_len:    the length of the dna sequence
                 r_scores: array that holds the score for each recognizer
                           and the cumulitive score in the last positon
                 c_scores: array that holds the score for each connector
                 c_lens:   array that holds the start position of the
                           first recognzier and the length of each gap
notes:           .
*/
void place_org(Organism* org,  const char* seq,  int s_len, double* r_scores, 
               double* c_scores, int* c_lens) {

  //n_rec and r_lens are instantiated as variables to make code more readable
  //m_len is the sum of the lengths of all recognizers
  int n_rec = org->len;
  int* r_lens = (int*)malloc(n_rec * sizeof(int));
  int m_len = 0;

  //iteration of each recognizer to calculate the minimum length of a sequence
  //and to write the length of each recognizer to the r_lens array
  for (int i = 0; i < n_rec; i++){
    r_lens[i] = org->recs[i].len;
    m_len += r_lens[i];
  }

  //calcualtion of the effective length and number of alignments
  int eff_len = s_len - m_len + n_rec;
  int n_align = s_len - m_len + 1;

  //keeps track of first valid placement for the current recognizer
  int f_offset = 0;

  //t_row is an array for writing temporary scores to during calculation
  //of cumulative scores
  double* t_row = (double*)malloc(n_align * sizeof(double));
  for (int i = 0; i < n_align; i++){
    t_row[i] = -INFINITY;
  }

  //c_row keeps track of the best cumulative score for each position that
  //the current recognizer can be placed
  double* c_row = (double*)calloc(n_align, sizeof(double));

  //rs_matrix: is an n_align by n_rec array that holds scores for each placement
  //           placement of each recognizer. 
  //           row n corresponds to recognizer n
  //           column n corresponds to the nth valid placement along the sequence
  //gs_matrix: is an n_align by n_rec - 1 array that holds the scores for each
  //           gap length of each connector. 
  //           row n corresponds to connector n
  //           column n corresponds to the score for the gap
  //           that was taken to get to that position for that connector
  //tr_matrix: is a matrix that holds the gap length from the previous placement
  //           row n corresponds to 
  //           column n corresponds to the length of the gap
  //           that was taken to get to that position
  double* rs_matrix = (double*)calloc(n_align * n_rec, sizeof(double));
  double* gs_matrix = (double*)calloc(n_align * (n_rec - 1), sizeof(double));
  int*    tr_matrix = (int*)calloc(n_align * (n_rec - 1), sizeof(int));

  //gap:     keep tracks of the length of the current gap
  //max_len: keeps track of the max length that has been precomputed up until
  //         for a given connector
  //g_score: keeps track of the score for the current gap for the current
  //         connector
  int     gap;
  int     max_len;
  double  g_score;

  //rec: points to the current recognzier
  //con: points to the current connector
  Recognizer* rec = NULL;
  Connector* con = NULL;

  for (int i = 0; i < n_rec; i++) {
    rec = &org->recs[i];

    //populates the row for the current recognizer with the score for
    //each possible placement along the sequence
    score_row(rec, seq + f_offset, n_align, rs_matrix + (i * n_align));

    //on the first iteration we dont need to do any cumulative score
    //calculations because there are no connectors involved yet
    if (i > 0) {
      con = &org->cons[i - 1];
      max_len = con->max_len;

      //j current position in the temp score row
      //k is the current position in the cumulative score array
      //we compare the score currently stored in our array
      //of temporary scores at j to each possible placement of the previous
      //recognizer, which are held in our array for cumulative scores at k
      for (int j = 0; j < n_align; j++){
        for (int k = 0; k < j + 1; k++){
          gap = j - k;

          //this determines whether we can compute the score using passed
          //pdfs and cdfs, or if we need to compute it right now
          if (gap > max_len)
            g_score = score_con(con, gap, s_len, eff_len, n_rec, false);
          else 
            g_score = score_con(con, gap, s_len, eff_len, n_rec, true);

          //if the score at some position k in the cumulative scores +
          //the score from the gap of lenght j - k is greater than the 
          //current best score at j in our temp scores array, we overwrite it
          if (t_row[j] < c_row[k] + g_score + rs_matrix[i * n_align + j]){
            t_row[j] = c_row[k] + g_score + rs_matrix[i * n_align + j];
            tr_matrix[(i - 1) * n_align + j] = gap;
            gs_matrix[(i - 1) * n_align + j] = g_score; 
          }

        }

      }

      //each position in the cumulative array is overwritten with the 
      //corresponding position in the temp array since the best possible
      //score for each posiiton was obtained by the end of the current
      //iteration
      //the temp row must be overwritten with the first valid placement
      //for each j so it is all reset to negative infinity
      for (int j = 0; j < n_align; j++){
        c_row[j] = t_row[j];
        t_row[j] = -INFINITY;
      }

    } else {

      //the best cumulative scores for the first iteration is just each
      //possible placement of the first recognizer
      for (int j = 0; j < n_align; j++){
        c_row[j] = rs_matrix[(i * n_align) + j];
      }
    }

    f_offset += r_lens[i];
  }

  //m_idx:    index that keeps track of the current position in the traceback
  //r_scores: is the array that we write the cumulative score to in the last
  //          position and the score for the corresponding recognizer
  int m_idx = max_index(c_row, n_align);
  r_scores[n_rec] = c_row[m_idx];

  //starting from the last recognizer, iterate backwards and get the get the
  //score and gap length for each recognizer and connector
  for (int i = n_rec - 1; i >= 0; i--) {
    r_scores[i] = rs_matrix[i * n_align + m_idx];
    if (i > 0){
      c_scores[i - 1] = gs_matrix[(i - 1) * n_align + m_idx];
      c_lens[i] = tr_matrix[(i - 1) * n_align + m_idx];
      m_idx -= tr_matrix[(i - 1) * n_align + m_idx];
    }
    c_lens[0] = m_idx;
  }

  free(r_lens);
  free(t_row);
  free(c_row);
  free(rs_matrix);
  free(gs_matrix);
  free(tr_matrix);
}
