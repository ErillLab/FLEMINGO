#include "organism.h"
void print_scores(Organism* org, double* r_scores, double* c_scores, int* c_lens, int n_rec){
  for (int i = 0; i < n_rec; i++){
    printf("|%c[%5.2f]|", org->recs[i].feat, r_scores[i]);
    if (i < n_rec - 1){
      printf("+|(%i)", c_lens[i + 1]);
      printf("c[%5.2f]|+", c_scores[i]);
    }
  }
  printf("=%5.2f\n", r_scores[n_rec]);

}

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

void parse_org(Organism *org, double* matrix, int* rec_lengths, double* models, double* edges, int* model_lengths, const char* rec_types, int num_recs, double* con_matrix, int max_len){
  org->len = num_recs;
  org->recs = (Recognizer*)malloc(num_recs * sizeof(Recognizer));
  org->cons = (Connector*)malloc((num_recs - 1) * sizeof(Connector));
  int m_offset = 0;
  int a_offset = 0;
  int c_offset = 0;
  int e_offset = 0;
  int shape_i = 0;

  for (int i = 0; i < num_recs; i++){
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
      e_offset += model_lengths[shape_i] + 2;
      shape_i += 1;
    }

    if (i != num_recs - 1){
        parse_con(&org->cons[i], 
                  con_matrix + 2 + c_offset, 
                  *(con_matrix + c_offset), 
                  *(con_matrix + c_offset + 1), 
                  max_len);      
        c_offset += 2 * max_len + 2;
    }

  }
}

void print_org(Organism *org) {
  for (int i = 0; i < org->len; i++){
    print_rec(&org->recs[i]);
    printf("\n");
    if (i < org->len - 1){
      print_con(&org->cons[i]);
      printf("\n");
    }
  }
}

void place_org( Organism* org,  const char* seq,  int s_len, double* r_scores, double* c_scores, int* c_lens) {

  int n_rec = org->len;
  int* r_lens = (int*)malloc(n_rec * sizeof(int));
  int m_len = 0;
  for (int i = 0; i < n_rec; i++){
    r_lens[i] = org->recs[i].len;
    m_len += r_lens[i];
  }

  int eff_len = s_len - m_len + n_rec;
  int n_align = s_len - m_len + 1;

  int f_offset = 0;

  double* t_row = (double*)malloc(n_align * sizeof(double));
  for (int i = 0; i < n_align; i++){
    t_row[i] = -INFINITY;
  }
  double* c_row = (double*)calloc(n_align, sizeof(double));

  double* rs_matrix = (double*)calloc(n_align * n_rec, sizeof(double));
  double* gs_matrix = (double*)calloc(n_align * (n_rec - 1), sizeof(double));
  int*   tr_matrix = (int*)calloc(n_align * (n_rec - 1), sizeof(int));
  int    gap       = 0;
  int    max_len = 0;
  double  g_score   = 0.0;

  Recognizer* rec = NULL;
  Connector* con = NULL;
  for (int i = 0; i < n_rec; i++) {
    rec = &org->recs[i];
    score_row(rec, seq + f_offset, n_align, rs_matrix + (i * n_align));
    if (i > 0) {
      con = &org->cons[i - 1];
      max_len = con->max_len;
      for (int j = 0; j < n_align; j++){
        for (int k = 0; k < j + 1; k++){
          gap = j - k;
          if (gap > max_len - 1)
            g_score = score_con(con, gap, s_len, eff_len, n_rec, false);
          else 
            g_score = score_con(con, gap, s_len, eff_len, n_rec, true);
          
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
    c_lens[0] = m_idx;
  }

  free(r_lens);
  free(t_row);
  free(c_row);
  free(rs_matrix);
  free(gs_matrix);
  free(tr_matrix);

  //print_scores(org, r_scores, c_scores, c_lens, n_rec);
  //print_placement(org, seq, s_len, c_lens);
}
