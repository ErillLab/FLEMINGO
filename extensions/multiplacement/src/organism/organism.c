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
  int n_rec = org->n_rec;
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
  int n_rec = org->n_rec; 
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

void debug_placement(Organism* org, const char* seq, int s_len, int* pos, double* r_scores, double* c_scores){
  print_scores(org, r_scores, c_scores, pos);
  print_placement(org, seq, s_len, pos);
}

/*
name:            parse_org

pre-conditions:  all information for recognizers and connectors has been loaded
                 number of recognziers is 1 greater than number of recognizers.
                 all of this information is passed from python, recognizers and
                 connectors store pointers to the 1D array passed from python.

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
                 model_lengths: an array containing the number of interval edges
                                for all shape recognzers
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

  org->n_rec  = num_recs;
  org->sum_rec_lens = 0;
  // allocates space for arrays of recognizers and connectors
  // recognizers and connectors just have pointers to information 
  // in arrays that were passed by from python
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
    org->sum_rec_lens += rec_lengths[i];
    // determine which parsing function to call based on recognizer type
    if (rec_types[i] == 'p') {
      parse_pssm(&org->recs[i],     // memory for the recognizer
                 matrix + m_offset, // matrix pointer shifted by the offset
                                    // so that the correct point in memory is
                                    // passed
                 rec_lengths[i]);   // length of the recognizer


      // each time a pssm is loaded, our offset in the 1D concatenated array
      // is increased by 4 * length of that recognizer
      // which is essentially (number of based * number of columns)
      m_offset += rec_lengths[i] * 4;

    } else { 
      parse_shape(
        &org->recs[i],          // memory for shape recognizer
        models + a_offset,      // null model offset by the correct amount
        models + (a_offset + model_lengths[shape_i]  - 1), // alt model offset
                                // by the correct amount
        edges + e_offset,       //intervals for recognizers models
        model_lengths[shape_i], // number of interval edges for the shape
        rec_lengths[i],         // length of the recognizer
        rec_types[i]            // type shape for the recognizer
      );

      // model_lengths[shape_i] is the number of interval edges the shapes models
      // have. null and alternative model are both (model_lengths[shape_i] - 1)
      // long. So to go to the start of the null model for the next shape, the
      // offset is increased by the total length occupied by the null and alt
      // models of the current shape.
      a_offset += (model_lengths[shape_i] - 1) * 2;

      // the array for intervals has exactly model_lengths[shape_i] edges
      // so to get to the intervals for the next shape, the offset for edges
      // is increased by that amount
      e_offset += model_lengths[shape_i];

      // information for shapes are all stored in 1D arrays,  we only want
      // to move to the next shape after encountering a shape object,
      // and not when we encounter a pssm object
      shape_i  += 1;
    }

    //connectors are passed pdfs and cdfs,
    //mu, sigma, and the length that is precomputed up until
    if (i < num_recs - 1){
        parse_con(
          &org->cons[i],                // memory for the connector
          con_matrix + 2 + c_offset,    // start of the pf array
          *(con_matrix + c_offset),     // mu
          *(con_matrix + c_offset + 1), // sigma
          max_len                       // max precomputed gap size
        ); 

        // each connector has mu, sigma, pfs, and aucs
        // mu and sigma both occupy one index each
        // pf and auc both occupy (max_len + 1) indices
        // since they are for gap sizes [0, max_len]
        // so getting from one connector to the next involves offsetting by
        // 1 for mu, 1 for sigma, (max_len + 1) for pfs, and (max_len + 1) for aucs
        // so  2 + (max_len + 1) * 2 in all
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
  for (int i = 0; i < org->n_rec; i++){
    if (org->recs[i].feat != 'p')
      has_shape = true;
  }

  if (has_shape == false)
    return;
  printf("*************************************************\n");
  for (int i = 0; i < org->n_rec; i++){
    print_rec(&org->recs[i]);
    printf("\n");
    if (i < org->n_rec - 1){
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
int place_org(Organism* org,  const char* seq,  int s_len, double* r_scores, 
               double* c_scores, int* c_lens) {
  
  double auc;

  // checks to see if the organism is too large to be placed on the sequence
  if (org->sum_rec_lens > s_len)
    return -1;

  // calcualtion of the effective length
  int eff_len = s_len - org->sum_rec_lens + org->n_rec;

  // n_rec_placements represents the number of possible placements of a single
  // recognizer along the sequence
  int n_rec_placements = s_len - org->sum_rec_lens + 1;

  // keeps track of first valid placement for the current recognizer
  int forward_offset = 0;

  // temp_curr_cumulative_scores is an array for writing temporary scores to during calculation
  // of cumulative scores
  double* temp_curr_cumulative_scores = (double*)malloc(n_rec_placements * sizeof(double));
  for (int i = 0; i < n_rec_placements; i++){

    // first evaluation of calculated score must always succeed so that a real gap
    // score is recorded
    temp_curr_cumulative_scores[i] = -INFINITY;
  }

  // prev_cumulative_scores keeps track of the best cumulative score for each position
  // this means that it holds the sum of scores for the best placements of all 
  // previous recognizers and connectors for each position
  double* prev_cumulative_scores = (double*)calloc(n_rec_placements, sizeof(double));

  // rec_scores_matrix: is an n_rec_placements by n_rec array that holds scores for each placement
  //            placement of each recognizer. 
  //            row n corresponds to recognizer n
  //            column n corresponds to the nth valid placement along the sequence
  // gap_scores_matrix: is an n_rec_placements by n_rec - 1 array that holds the scores for each
  //            gap length of each connector. 
  //            row n corresponds to connector n
  //            column n corresponds to the score for the gap
  //            that was taken to get to that position for that connector
  // traceback_matrix: is a matrix that holds the gap length from the previous placement
  //            row n corresponds to 
  //            column n corresponds to the length of the gap
  //            that was taken to get to that position
  double* rec_scores_matrix = (double*)calloc(n_rec_placements * org->n_rec, sizeof(double));
  double* gap_scores_matrix = (double*)calloc(n_rec_placements * (org->n_rec - 1), sizeof(double));
  int*    traceback_matrix = (int*)calloc(n_rec_placements * (org->n_rec - 1), sizeof(int));

  // gap:     keep tracks of the length of the current gap
  // max_len: keeps track of the max length that has been precomputed up until
  //          for a given connector
  // g_score: keeps track of the score for the current gap for the current
  //          connector
  int     gap;
  int     max_len;
  double  g_score;

  // rec: points to the current recognzier
  // con: points to the current connector
  Recognizer* rec = NULL;
  Connector* con = NULL;

  for (int i = 0; i < org->n_rec; i++) {
    rec = &org->recs[i];

    //populates the row for the current recognizer with the score for
    //each possible placement along the sequence
    score_row(rec, seq + forward_offset, n_rec_placements, rec_scores_matrix + (i * n_rec_placements));

    //on the first iteration we dont need to do any cumulative score
    //calculations because there are no connectors involved yet
    if (i > 0) {
      con = &org->cons[i - 1];
      max_len = con->max_len;

      //sometimes the length of the sequence is greater than the max precomputed length
      //if this is the case then we need to compute the cdf to used for the placement
      if (s_len - org->sum_rec_lens > max_len)
        auc = log2f(norm_cdf(eff_len - org->n_rec + 0.5, con->mu, con->sigma) 
                  - norm_cdf(-0.5, con->mu, con->sigma));
      else
        auc = con->cdf[s_len - org->sum_rec_lens];      

      //j current position in the temp score row
      //k is the current position in the cumulative score array
      //we compare the score currently stored in our array
      //of temporary scores at j to each possible placement of the previous
      //recognizer, which are held in our array for cumulative scores at k
      for (int j = 0; j < n_rec_placements; j++)
        for (int k = 0; k < j + 1; k++){
          gap = j - k;

          //this determines whether we can compute the score using passed
          //pdfs and cdfs, or if we need to compute it right now
          if (gap > max_len)
            g_score = score_con(con, gap, s_len, eff_len, org->n_rec, auc, false);
          else 
            g_score = score_con(con, gap, s_len, eff_len, org->n_rec, auc, true);

          //if the score at some position k in the cumulative scores +
          //the score from the gap of lenght j - k is greater than the 
          //current best score at j in our temp scores array, we overwrite it

          // evaluating the placement at position j and all possible gap lengths
          // to get from the previous placed cumulative scores to the placement
          // at position j 
          if (temp_curr_cumulative_scores[j] < prev_cumulative_scores[k] + g_score + rec_scores_matrix[i * n_rec_placements + j]){
            temp_curr_cumulative_scores[j] = prev_cumulative_scores[k] + g_score + rec_scores_matrix[i * n_rec_placements + j];
            traceback_matrix[(i - 1) * n_rec_placements + j] = gap;
            gap_scores_matrix[(i - 1) * n_rec_placements + j] = g_score; 
          }

        }
      //each position in the cumulative array is overwritten with the 
      //corresponding position in the temp array since the best possible
      //score for each posiiton was obtained by the end of the current
      //iteration
      //the temp row must be overwritten with the first valid placement
      //for each j so it is all reset to negative infinity
      for (int j = 0; j < n_rec_placements; j++){
        prev_cumulative_scores[j] = temp_curr_cumulative_scores[j];
        temp_curr_cumulative_scores[j] = -INFINITY;
      }

    } else {

      //the best cumulative scores for the first iteration is just each
      //possible placement of the first recognizer
      for (int j = 0; j < n_rec_placements; j++){
        prev_cumulative_scores[j] = rec_scores_matrix[(i * n_rec_placements) + j];
      }
    }

    forward_offset += rec->len;
  }

  //m_idx:    index that keeps track of the current position in the traceback
  //r_scores: is the array that we write the cumulative score to in the last
  //          position and the score for the corresponding recognizer
  int m_idx = max_index(prev_cumulative_scores, n_rec_placements);
  r_scores[org->n_rec] = prev_cumulative_scores[m_idx];

  //starting from the last recognizer, iterate backwards and get the get the
  //score and gap length for each recognizer and connector
  for (int i = org->n_rec - 1; i >= 0; i--) {

    // i * n_rec_placements gets to the row for the current recognizer
    // adding m_idx yields the appropriate column
    r_scores[i] = rec_scores_matrix[i * n_rec_placements + m_idx];

    // there are no connectors before the first recognizer so we stop
    // dealing with connectors when we hit the first recognizer
    if (i > 0){

      // connectors are offset 1 from the current recognizer
      c_scores[i - 1] = gap_scores_matrix[(i - 1) * n_rec_placements + m_idx];

      // c_lens is not offset one since it also holds the starting position
      // of the first recognizer
      c_lens[i] = traceback_matrix[(i - 1) * n_rec_placements + m_idx];

      // m_idx steps back by the gap length value at that postion in the
      // traceback matrix
      m_idx -= traceback_matrix[(i - 1) * n_rec_placements + m_idx];

    }

  }
  c_lens[0] = m_idx;

  free(temp_curr_cumulative_scores);
  free(prev_cumulative_scores);
  free(rec_scores_matrix);
  free(gap_scores_matrix);
  free(traceback_matrix);

  //debug_placement(org, seq, s_len, c_lens, r_scores, c_scores);
  return 0;
}
