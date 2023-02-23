#include "_multiplacement.h"
#include "_aux.h"

void fill_row_pssm(const char* seq, int len_seq, int row, int curr_rec, float rec_matrices[], int pssm_score_offset, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[]){
  float score = 0.0; 
  for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
    score = 0.0;
    for (int k = 0; k < rec_length; k++) {
      switch (seq[j + k]) {
        case 'A':
        case 'a':
          score += rec_matrices[(pssm_score_offset + k) * 4 + 0];
          break;
        case 'G':
        case 'g':
          score += rec_matrices[(pssm_score_offset + k) * 4 + 1];
          break;
        case 'C':
        case 'c':
          score += rec_matrices[(pssm_score_offset + k) * 4 + 2];
          break;
        case 'T':
        case 't':
          score += rec_matrices[(pssm_score_offset + k) * 4 + 3];
          break;
      }
    }
    //printf("score: %f\n", score);
    //printf("%f\n", score);
    score_matrix[(row * num_alignments) + j - forward_offset] = score;

  }
}

void fill_row_mgw(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins){
  float score = 0.0; 
  float null_freq = 0.0; 
  float alt_freq = 0.0; 
  float sum = 0.0;
  int index = 0;
  int num_pentamers = rec_length - 4;
  float *pentamer_scores = malloc(num_pentamers * sizeof(float));
  
  for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
    score = 0.0;
    sum = 0.0;
    //printf("j = %i\n", j);
    //printf("num_pent = %i\n", num_pentamers);
    for (int k = 0; k < num_pentamers; k++) {
      index = 0;
      for (int l = k; l < k + 5; l++){
        //printf("hi\n");
        //printf("%c", seq[j + l]);
        switch (seq[j + l]) {
          case 'A':
          case 'a':
            index += pow(4, 4 - (l - k)) * 0;
            break;
          case 'G':
          case 'g':
            index += pow(4, 4 - (l - k)) * 1;
            break;
          case 'C':
          case 'c':
            index += pow(4, 4 - (l - k)) * 2;
            break;
          case 'T':
          case 't':
            index += pow(4, 4 - (l - k)) * 3;
            break;
        }
      }
      //printf(": index: %i | MGW_SCORES[%i]: %f\n", index, index, MGW_SCORES[index]); 
      //printf("k = %i\n", k);    
      pentamer_scores[k] = MGW_SCORES[index];
    }

    //still need to get actual score from llr
    score = shape_average(pentamer_scores, num_pentamers);
    null_freq = get_bin_frequency(score, bin_frequencies, bin_edges, num_bins);
    //printf("\ngetting null freq from score: %f\n", score);
    alt_freq = get_bin_frequency(score, bin_frequencies + num_bins, bin_edges, num_bins);
    //printf("calculated score of %f\n from alt freq: %f\n and null freq: %f\n", log2(alt_freq/null_freq), alt_freq, null_freq);
    score_matrix[(row * num_alignments) + j - forward_offset] = log2(alt_freq/null_freq);
  }
  free(pentamer_scores);
}

void fill_row_prot(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins){

  float score = 0.0; 
  float null_freq = 0.0; 
  float alt_freq = 0.0; 
  float sum = 0.0;
  int index = 0;
  int num_pentamers = rec_length - 4;
  float *pentamer_scores = malloc(num_pentamers * sizeof(float));
  
  for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
    score = 0.0;
    sum = 0.0;
    //printf("j = %i\n", j);    
    //printf("num_pent = %i\n", num_pentamers);
    for (int k = 0; k < num_pentamers; k++) {
      index = 0;
      for (int l = k; l < k + 5; l++){
        //printf("hi\n");
        //printf("%c", seq[j + l]);
        switch (seq[j + l]) {
          case 'A':
          case 'a':
            index += pow(4, 4 - (l - k)) * 0;
            break;
          case 'G':
          case 'g':
            index += pow(4, 4 - (l - k)) * 1;
            break;
          case 'C':
          case 'c':
            index += pow(4, 4 - (l - k)) * 2;
            break;
          case 'T':
          case 't':
            index += pow(4, 4 - (l - k)) * 3;
            break;
        }
      }
      //printf("k = %i\n", k);    
      pentamer_scores[k] = PROT_SCORES[index];
    }

    //still need to get actual score from llr

    score = shape_average(pentamer_scores, num_pentamers);
    null_freq = get_bin_frequency(score, bin_frequencies, bin_edges, num_bins);
    //printf("\ngetting null freq from score: %f\n", score);
    alt_freq = get_bin_frequency(score, bin_frequencies + num_bins, bin_edges, num_bins);
    //printf("calculated score of %f\n from alt freq: %f\n and null freq: %f\n", log2(alt_freq/null_freq), alt_freq, null_freq);
    score_matrix[(row * num_alignments) + j - forward_offset] = log2(alt_freq/null_freq);
  }
  free(pentamer_scores);
}

void fill_row_helt(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins){
  float  score                 = 0.0; 
  float null_freq = 0.0; 
  float alt_freq = 0.0; 
  float  sum                   = 0.0;
  int    index                 = 0;
  int    num_pentamers         = rec_length - 4;
  float *pentamer_scores       = malloc(2 * num_pentamers  * sizeof(float));
  
  for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
    score = 0.0;
    sum   = 0.0;
    //printf("j = %i\n", j);    
    //printf("num_pent = %i\n", num_pentamers);
    
    for (int k = 0; k < num_pentamers; k++) {
      index = 0;

      for (int l = k; l < k + 5; l++){
        //printf("hi\n");
        //printf("%c", seq[j + l]);
        switch (seq[j + l]) {
          case 'A':
          case 'a':
            index += pow(4, 4 - (l - k)) * 0;
            break;
          case 'G':
          case 'g':
            index += pow(4, 4 - (l - k)) * 1;
            break;
          case 'C':
          case 'c':
            index += pow(4, 4 - (l - k)) * 2;
            break;
          case 'T':
          case 't':
            index += pow(4, 4 - (l - k)) * 3;
            break;
        }

      }
      //printf("k = %i\n", k);    
      pentamer_scores[k]     = HELT_SCORES[index];
      pentamer_scores[k + num_pentamers] = HELT_SCORES[1024 + index];
    }
    score = shape_average(pentamer_scores, num_pentamers * 2);
    //still need to get actual score from llr
    null_freq = get_bin_frequency(score, bin_frequencies, bin_edges, num_bins);
    //printf("getting null freq from score: %f\n", score);
    alt_freq = get_bin_frequency(score, bin_frequencies + num_bins, bin_edges, num_bins);
   // printf("calculated score of %f\n from alt freq: %f\n and null freq: %f\n", log2(alt_freq/null_freq), alt_freq, null_freq);
    score_matrix[(row * num_alignments) + j - forward_offset] = log2(alt_freq/null_freq);

  }
  free(pentamer_scores);
}

void fill_row_roll(const char* seq, int len_seq, int row, int curr_rec, int rec_length, int num_alignments, int forward_offset, int reverse_offset, float score_matrix[], float bin_frequencies[], float bin_edges[], int num_bins){
  float  score                 = 0.0; 
  float null_freq = 0.0; 
  float alt_freq = 0.0; 
  float  sum                   = 0.0;
  int    index                 = 0;
  int    num_pentamers         = rec_length - 4;
  float *pentamer_scores       = malloc(2 * num_pentamers  * sizeof(float));
  
  for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
    score = 0.0;
    sum   = 0.0;
    //printf("j = %i\n", j);    
    //printf("num_pent = %i\n", num_pentamers);
    
    for (int k = 0; k < num_pentamers; k++) {
      index = 0;

      for (int l = k; l < k + 5; l++){
        //printf("hi\n");
        //printf("%c", seq[j + l]);
        switch (seq[j + l]) {
          case 'A':
          case 'a':
            index += pow(4, 4 - (l - k)) * 0;
            break;
          case 'G':
          case 'g':
            index += pow(4, 4 - (l - k)) * 1;
            break;
          case 'C':
          case 'c':
            index += pow(4, 4 - (l - k)) * 2;
            break;
          case 'T':
          case 't':
            index += pow(4, 4 - (l - k)) * 3;
            break;
        }

      }
      //printf("k = %i\n", k);    
      pentamer_scores[k]     = ROLL_SCORES[index];
      pentamer_scores[k + num_pentamers] = ROLL_SCORES[1024 + index];
    }
    score = shape_average(pentamer_scores, num_pentamers * 2);
    //still need to get actual score from llr

    null_freq = get_bin_frequency(score, bin_frequencies, bin_edges, num_bins);
    //printf("getting null freq from score: %f\n", score);
    alt_freq = get_bin_frequency(score, bin_frequencies + num_bins, bin_edges, num_bins);
    //printf("calculated score of %f\n from alt freq: %f\n and null freq: %f\n", log2(alt_freq/null_freq), alt_freq, null_freq);
    score_matrix[(row * num_alignments) + j - forward_offset] = log2(alt_freq/null_freq);
  }
  free(pentamer_scores);
}

void traceback(int num_rec, int len_seq, float* con_matrices, float* rec_score_matrix, int num_alignments, float* rec_alignments, int* con_alignments, float* rec_scores, float* con_scores, int* con_lengths, int max_seq_len, int effective_len, bool is_precomputed){
  //printf("made it traceback\n");
  // finding gap lengths to ref orma alignments
  // starting from the index of the greatest score in alignments
  // trace back of gap alignments is conducted by subtracting the value
  // at the max score index from the index and going up a row
  // this is repeated until we know all the gap lengths in the alignment
  int index = max_index(rec_alignments, num_alignments);
  con_lengths[num_rec - 1] = con_alignments[(num_rec - 2) * num_alignments + index];

  // the cumulative best score is written into the last index of the scores
  // array
  rec_scores[num_rec] = rec_alignments[index];

  //the start and end of this for loop are weird because the first and last
  //indices are filled manually, so these bounds only encompass the middle
  //scores
  for (int i = num_rec - 3; i >= 0; i--) {
    con_lengths[i + 1] = con_alignments[i * num_alignments + index - con_lengths[i + 2]];
    index -= con_lengths[i + 2];
  }
  con_lengths[0] = index - con_lengths[1];

  int gapOffset = 0;
  for (int i = 0; i < num_rec - 1; i++) {
    //con_scores[i] = get_score(con_matrices[i * max_seq_len + con_lengths[i + 1]], effective_len, num_rec, con_lengths[i + 1]);
    con_scores[i] = get_score(con_matrices, len_seq, effective_len, num_rec, con_lengths[i + 1], max_seq_len, i, is_precomputed);
  }

  // scores for each PSSM are filled by iterating over the score matrix
  // and using the appropriate gap lengths as a cumulative offset
  for (int i = 0; i < num_rec; i++) {
    gapOffset += con_lengths[i];
    rec_scores[i] = rec_score_matrix[i * num_alignments + gapOffset];
  }

}

void fill_traceback_matrix(float *score_matrix, int num_alignments, float *gapMatrix, int *cols, int num_rec, int len_seq, float *con_scores, float *rec_scores, int *con_lengths, int max_length, bool is_precomputed) {
  //printf("Made it traceback matrix\n"); 
  int gap_length = 0;
  int effective_length = len_seq;
  int sum_of_lengths = 0;

  for (int i = 0; i < num_rec; i++){
    sum_of_lengths += cols[i];
  }

  effective_length -= sum_of_lengths;
  //  number of total alignments by number of pssms
  //  first index in each column holds current max score for that index
  //  all other indices hold gap lengths that got that alignment
  float *alignments = calloc(num_alignments, sizeof(*score_matrix));
  int *gap_alignments = calloc(num_alignments * (num_rec - 1), sizeof(*con_lengths));
  float *temp_max_scores = calloc(num_alignments, sizeof(*score_matrix));
  int *temp_gap_lengths = calloc(num_alignments, sizeof(*con_lengths));
  float temp_gap_score = 0.0;

  // start with first row as our current max
  for (int i = 0; i < num_alignments; i++) {
    alignments[i] = score_matrix[i];
    temp_max_scores[i] = score_matrix[i];
  }

  // for each connector (number of pssms - 1) populate alignments with
  // current maximum score for that index

  // overview: the algorithm iterates over the next row in the scoresMatrix
  // (first row is our starting cumulative score), for each index, we compare
  // the scores obtained by summing the PSSM score at that index with each
  // previous cumulative alignment score for k <= j, when k == j, the gap length
  // is 0
  for (int i = 1; i < num_rec; i++) {

    for (int j = 0; j < num_alignments; j++) {

      for (int k = 0; k <= j; k++) {

        // every column before or equal to the current column is a valid
        // alignment. Scores for each alignment
        // gap_length = difference between column normalized j and k
        gap_length = j - k;

        //compute on fly
        //temp_gap_score = get_score(gapMatrix[(i - 1) * max_length + gap_length], effective_length, num_rec, gap_length);
        temp_gap_score = get_score(gapMatrix, len_seq, effective_length, num_rec, gap_length, max_length, i - 1, is_precomputed);
        if (k == 0) {
          temp_max_scores[j] = alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j];
          temp_gap_lengths[j] = gap_length;
        }else{
          if (temp_max_scores[j] < alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j]) {
              temp_max_scores[j] = alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j];
              temp_gap_lengths[j] = gap_length;
          }
        }
      }
    }

    // we must reset our temp arrays so that they will be overwritten
    // when doing the comparison of scores
    for (int l = 0; l < num_alignments; l++) {
      alignments[l] = temp_max_scores[l];
      gap_alignments[(i - 1) * num_alignments + l] = temp_gap_lengths[l];
      temp_max_scores[l] = -INFINITY;
      temp_gap_lengths[l] = 0;
    }
  }
  free(temp_max_scores);
  free(temp_gap_lengths);

  traceback(num_rec, len_seq, gapMatrix, score_matrix, num_alignments, alignments, gap_alignments, rec_scores, con_scores, con_lengths, max_length, effective_length, is_precomputed);
  free(alignments);
  free(gap_alignments);  
}

void fill_matrix(const char seq[], int len_seq, float rec_matrices[], int rec_lengths[], const char rec_types[], int num_rec,float score_matrix[], int num_alignments, float bin_frequencies[], float bin_edges[], int num_bins[]) {
  //printf("made it matrix\n");
  int forward_offset = 0;
  int reverse_offset = 0;
  int pssm_score_offset = 0;
  int curr_count_pssm_rec = 0;
  int curr_count_shape_rec = 0;
  int* rec_offsets = (int*)calloc(num_rec, sizeof(int));
  int* pssm_offsets = (int*)calloc(num_rec, sizeof(int));
  int* bin_offsets = (int*)calloc(num_rec, sizeof(int));
  //this is for getting the correct offset for each recognizer within its own concatenated array
  for (int i = 0; i < num_rec; i++){
    //printf("%c\n", rec_types[i]);
    if (rec_types[i] == 'P' || rec_types[i] == 'p'){
     rec_offsets[i] = curr_count_pssm_rec;
     curr_count_pssm_rec++; 
     pssm_offsets[i] = rec_lengths[i];
    }else{
      rec_offsets[i] = curr_count_shape_rec;
      bin_offsets[i] = int_arr_sum(num_bins, curr_count_shape_rec);
      curr_count_shape_rec++;
    }     
  }

  if (curr_count_shape_rec > 1){
    for (int i = 0; i < num_rec; i++){
      //printf("bin_offset at %i is %i\n", i, bin_offsets[i]);
    }    

    for (int i = 0; i < curr_count_shape_rec; i++){
      //printf("bin size:\n", num_bins[i]);
    }
  }
  // pre computes alignments of each pssm at each possible position
  // i = current recognizer
  for (int i = 0; i < num_rec; i++) {
    forward_offset = get_forward_offset(i, rec_lengths, num_rec);
    reverse_offset = get_reverse_offset(i, rec_lengths, num_rec);
    //printf("doing rec: %i: of type %c\n", i, rec_types[i]);
    switch(rec_types[i]){
    case 'P':
    case 'p':
      //printf("fill row pssm: pssm offset: %i\n", rec_offsets[i]);
      pssm_score_offset = int_arr_sum(pssm_offsets, i);
      fill_row_pssm(seq, len_seq, i, rec_offsets[i], rec_matrices, pssm_score_offset, rec_lengths[i], num_alignments, forward_offset, reverse_offset, score_matrix);
      break;
    case 'M':
    case 'm':
      fill_row_mgw(seq, len_seq, i, rec_offsets[i], rec_lengths[i], num_alignments, forward_offset, reverse_offset, score_matrix, bin_frequencies, bin_edges, num_bins[0]);
      break;
    case 'T':
    case 't':
      fill_row_prot(seq, len_seq, i, rec_offsets[i], rec_lengths[i], num_alignments, forward_offset, reverse_offset, score_matrix, bin_frequencies, bin_edges, num_bins[0]);
      break;
    case 'H':
    case 'h':
      fill_row_helt(seq, len_seq, i, rec_offsets[i], rec_lengths[i], num_alignments, forward_offset, reverse_offset, score_matrix, bin_frequencies, bin_edges, num_bins[0]);
      break;
    case 'R':
    case 'r':
      fill_row_roll(seq, len_seq, i, rec_offsets[i], rec_lengths[i], num_alignments, forward_offset, reverse_offset, score_matrix, bin_frequencies, bin_edges, num_bins[0]);
      break;
      break;
    }
  }
}

float con(float* con, int s_len, int eff_len, int n_rec, int gap){
  return log2f(get_numerator(s_len, gap, con[0], con[1])/get_denominator(gap + 1, n_rec, eff_len));
}

void pssm(char* seq, int s_len, float* matrix, int len, float* row){
  float score = 0.0;
  for (int i = 0; i < s_len + 1; i++) {
    score = 0.0;
    for (int j = 0; j < len; j++) {
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

void mgw(char* seq, int s_len, float** rec, int len, int n_bin, float* row){
  int n_pent = len - 4;
  float alt_f = 0.0;
  float null_f = 0.0;
  float* alt = rec[0];
  float* null = rec[1];
  float* edges = rec[2];
  float* pent_s = (float*)malloc(n_pent * sizeof(float));
  float score = 0.0;
  int idx = 0;
  for (int i = 0; i < len + 1; i++){
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
    alt_f = get_bin_frequency(score, alt, edges, n_bin);
    null_f = get_bin_frequency(score, null, edges, n_bin);
    row[i] = log2f(alt_f / null_f);
    row[i] = 1.2;

  }
  free(pent_s);
}

void prot(char* seq, int s_len, float** rec, int len, int n_bin, float* row){
  
}

void roll(char* seq, int s_len, float** rec, int len, int n_bin, float* row){
  
}

void helt(char* seq, int s_len, float** rec, int len, int n_bin, float* row){
  
}

void shape(char* seq, int s_len, float** rec, int len, char feat, int n_bin, float* row){
  switch(feat){
    case 'm':
    case 'M':
      mgw(seq, s_len, rec, len, n_bin, row);
      break;

    case 't':
    case 'T':
      prot(seq, s_len, rec, len, n_bin, row);
      break;

    case 'r':
    case 'R':
      roll(seq, s_len, rec, len, n_bin, row);
      break;

    case 'h':
    case 'H':
      helt(seq, s_len, rec, len, n_bin, row);
      break;
    
  }
}

void row(char* seq, int s_len, float** rec, int len, char feat, int n_bin, float* row){
  if (feat == 'p'){
    pssm(seq, s_len, rec[0], len, row);
  }else{
    shape(seq, s_len, rec, len, feat, n_bin, row);
  }
}

void place(char* seq, int s_len, float*** recs, int* r_lens, char* r_types, int n_rec, int* b_lens, float** cons, int m_len, float* r_scores, float* c_scores, int* c_lens){
  int o_len = 0;
  int i_shape = 0;
  for (int i = 0; i < n_rec; i++){
    o_len += r_lens[i];
  }

  int eff_len = s_len - o_len - n_rec;
  int n_align = s_len - o_len + 1;

  int f_offset = 0;
  float* t_row = (float*)malloc(n_align * sizeof(float));
  for (int i = 0; i < n_align; i++){
    t_row[i] = -INFINITY;
  }
  float* c_row = (float*)calloc(n_align, sizeof(float));

  float* rs_matrix = (float*)calloc(n_align * n_rec, sizeof(float));
  float* gs_matrix = (float*)calloc(n_align * (n_rec - 1), sizeof(float));
  int*   tr_matrix = (int*)calloc(n_align * (n_rec - 1), sizeof(int));
  int    gap       = 0;
  float  g_score   = 0.0;
  
  for (int i = 0; i < n_rec; i++) {
    row(seq + f_offset, n_align, recs[i], r_lens[i], r_types[i], b_lens[i_shape], rs_matrix + (i * n_align));
    if (i > 0) {

      for (int j = 0; j < n_align; j++){
        for (int k = 0; k < j + 1; k++){
          
          gap = j - k;
          g_score = con(cons[i - 1], s_len, eff_len, n_rec, gap);
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
  r_scores[n_rec] = c_row[m_idx];
  for (int i = n_rec - 1; i >= 0; i--) {
    r_scores[i] = rs_matrix[i * n_align + m_idx];
    if (i > 0){
      c_scores[i - 1] = gs_matrix[(i - 1) * n_align + m_idx];
      c_lens[i] = tr_matrix[(i - 1) * n_align + m_idx];
      m_idx -= tr_matrix[(i - 1) * n_align + m_idx];
    }
  }

  if (n_rec > 1){
    c_lens[0] = m_idx - c_lens[1];
  } else {
    c_lens[0] = m_idx;
  }

  /*
  printf("%f: \n", r_scores[n_rec]);
  for (int i = 0; i < n_rec; i++){
    printf("rec %i: %f \n", i, r_scores[i]);
    if (i < n_rec - 1){
      printf("con: %i: %i, %f\n", i, c_lens[i], c_scores[i]);
    }
  }
  */

  printf("a\n");
  free(c_row);
  printf("b\n");
  free(t_row);
  printf("c\n");
  free(rs_matrix);
  printf("d\n");
  free(gs_matrix);
  printf("e\n");
  free(tr_matrix);
  printf("f\n");
  fflush(stdout);
  c_row = NULL;
  t_row = NULL;
  rs_matrix = NULL;
  gs_matrix = NULL;
  tr_matrix = NULL;
}