#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
const int NUM_BASES = 4;

int get_forward_offset(int index, int cols[], int num_rec) {
  // finds the first possible possition for a pssm
  // based on number of columns of preceding pssms

  int offset = 0;
  for (int i = 0; i < index; i++) {
    offset += cols[i];
  }
  return offset;
}

int get_reverse_offset(int index, int cols[], int num_rec) {
  // finds the last possible possition for a pssm
  // based on number of columns of subsequenct pssms

  int offset = 0;
  for (int i = num_rec - 1; i >= index; i--) {
    offset += cols[i];
  }

  // this is important to get the right number of possible alignments
  return offset - 1;
}

int max_index(float *arr, int size) {
  int max_index = 0;
  for (int i = 0; i < size; i++) {

    if (arr[i] > arr[max_index]) {
      max_index = i;
    }
  }
  return max_index;
}

void traceback(float *score_matrix, int num_alignments, float *gapMatrix,
               int *cols, int num_rec, int len_seq, float *gapScores,
               float *scores, int *gaps) {
  int gap_length = 0;

  // printf("last in traceback\n");
  //  number of total alignments by number of pssms
  //  first index in each column holds current max score for that index
  //  all other indices hold gap lengths that got that alignment
  float *alignments = PyMem_Calloc(num_alignments, sizeof(*score_matrix));
  int *gap_alignments =
      PyMem_Calloc(num_alignments * (num_rec - 1), sizeof(*gaps));
  float *temp_max_scores = PyMem_Calloc(num_alignments, sizeof(*score_matrix));
  int *temp_gap_lengths = PyMem_Calloc(num_alignments, sizeof(*gaps));

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
    // printf("Calculating alignments for PSSM %i\n", i);

    for (int j = 0; j < num_alignments; j++) {

      for (int k = 0; k <= j; k++) {

        // every column before or equal to the current column is a valid
        // alignment. Scores for each alignment
        // gap_length = difference between column normalized j and k
        gap_length = j - k;

        if (k == 0) {
          temp_max_scores[j] = alignments[k] +
                               gapMatrix[(i - 1) * len_seq + gap_length] +
                               score_matrix[i * num_alignments + j];
          temp_gap_lengths[j] = gap_length;
        } else {
          if (temp_max_scores[j] <
              alignments[k] + gapMatrix[(i - 1) * len_seq + gap_length] +
                  score_matrix[i * num_alignments + j]) {
            temp_max_scores[j] = alignments[k] +
                                 gapMatrix[(i - 1) * len_seq + gap_length] +
                                 score_matrix[i * num_alignments + j];
            temp_gap_lengths[j] = gap_length;
          }
        }
      }
    }

    for (int l = 0; l < num_alignments; l++) {
      alignments[l] = temp_max_scores[l];
      gap_alignments[(i - 1) * num_alignments + l] = temp_gap_lengths[l];
      temp_max_scores[l] = -INFINITY;
      temp_gap_lengths[l] = 0;
    }
  }

  // finding gap lengths to ref orma alignments
  // starting from the index of the greatest score in alignments
  // trace back of gap alignments is conducted by subtracting the value
  // at the max score index from the index and going up a row
  // this is repeated until we know all the gap lengths in the alignment
  int index = max_index(alignments, num_alignments);
  gaps[num_rec - 1] = gap_alignments[(num_rec - 2) * num_alignments + index];

  // the cumulative best score is written into the last index of the scores
  // array
  scores[num_rec] = alignments[index];

  for (int i = num_rec - 3; i >= 0; i--) {
    gaps[i + 1] = gap_alignments[i * num_alignments + index - gaps[i + 2]];
    index -= gaps[i + 2];
  }
  gaps[0] = index - gaps[1];

  int gapOffset = 0;
  for (int i = 0; i < num_rec - 1; i++) {
    gapScores[i] = gapMatrix[i * len_seq + gaps[i + 1]];
  }

  // scores for each PSSM are filled by iterating over the score matrix
  // and using the appropriate gap lengths as a cumulative offset
  for (int i = 0; i < num_rec; i++) {
    gapOffset += gaps[i];
    scores[i] = score_matrix[i * num_alignments + gapOffset];
  }
  PyMem_Free(alignments);
  PyMem_Free(gap_alignments);
  PyMem_Free(temp_max_scores);
  PyMem_Free(temp_gap_lengths);
}

void fill_matrix(const char seq[], int len_seq, float pssm[], int cols[],
                 int num_rec, float score_matrix[], int num_alignments) {
  // length of the seq by number of pssms

  // printf("last in fill_matrix\n");
  float score = 0;
  int forward_offset = 0;
  int reverse_offset = 0;

  // pre computes alignments of each pssm at each possible position
  // i = current pssm
  // j = starting position on seq for computing score
  // k = current column in pssm
  // time complexity should be (num_recs * numAignments * avgLengthPSSMs)

  for (int i = 0; i < num_rec; i++) {
    // printf("Calculating scores for PSSM %i\n", i);
    // these functions are unnessesary, change to using += and -= cols[i]
    forward_offset = get_forward_offset(i, cols, num_rec);
    reverse_offset = get_reverse_offset(i, cols, num_rec);

    for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
      score = 0;
      // printf("\nstarting positon %i\n", j);
      for (int k = 0; k < cols[i]; k++) {
        switch (seq[j + k]) {
        case 'A':
        case 'a':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forward_offset + k) * 4 + 0]);
          score += pssm[(forward_offset + k) * 4 + 0];
          break;
        case 'G':
        case 'g':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forward_offset + k) * 4 + 1]);
          score += pssm[(forward_offset + k) * 4 + 1];
          break;
        case 'C':
        case 'c':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forward_offset + k) * 4 + 2]);
          score += pssm[(forward_offset + k) * 4 + 2];
          break;
        case 'T':
        case 't':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forward_offset + k) * 4 + 3]);
          score += pssm[(forward_offset + k) * 4 + 3];
          break;
        }
      }
      score_matrix[(i * num_alignments) + j - forward_offset] = score;
    }
  }
}

static int matrix_converter(PyObject *object, void *address) {
  const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
  char datatype;
  Py_buffer *view = address;

  if (object == NULL)
    goto exit;
  if (PyObject_GetBuffer(object, view, flags) == -1) {
    PyErr_SetString(PyExc_RuntimeError,
                    "position-weight matrix is not an array");
    return 0;
  }
  datatype = view->format[0];
  switch (datatype) {
  case '@':
  case '=':
  case '<':
  case '>':
  case '!':
    datatype = view->format[1];
    break;
  default:
    break;
  }
  return Py_CLEANUP_SUPPORTED;

exit:
  PyBuffer_Release(view);
  return 0;
}

static char calculate__doc__[] =
    "    calculate(sequence, recognizers, recognizer_lengths, connectors, "
    "recognizer_scores, connector_scores, connector_legnths)\n"
    "\n"
    "This function computes optimal placement for a \n"
    "transcription factor composed of PSSM recogniers\n"
    "and variable length connectors.\n";

static PyObject *py_calculate(PyObject *self, PyObject *args,
                              PyObject *keywords) {
  const char *seq;
  static char *kwlist[] = {
      "sequence",   "rec_matrices", "rec_lengths", "con_matrices",
      "rec_scores", "con_scores",   "con_lengths", NULL};
  // printf("back in c\n");
  Py_ssize_t len_seq;
  Py_ssize_t num_rec;
  PyObject *result = Py_None;
  Py_buffer rec_matrices;
  Py_buffer con_matrices;
  Py_buffer rec_lengths;
  Py_buffer rec_scores;
  Py_buffer con_scores;
  Py_buffer con_lengths;

  // double* test = malloc(100000000000000000000000 * sizeof(double));
  // free(test);

  rec_matrices.obj = NULL;
  con_matrices.obj = NULL;
  rec_lengths.obj = NULL;
  rec_scores.obj = NULL;
  con_scores.obj = NULL;
  con_lengths.obj = NULL;
  // printf("now parsing\n");
  if (!PyArg_ParseTupleAndKeywords(
          args, keywords, "y#O&O&O&O&O&O&", kwlist, &seq, &len_seq,
          matrix_converter, &rec_matrices, matrix_converter, &rec_lengths,
          matrix_converter, &con_matrices, matrix_converter, &rec_scores,
          matrix_converter, &con_scores, matrix_converter, &con_lengths))
    return NULL;
  // printf("Parsed successfully\n");
  num_rec = rec_lengths.shape[0];
  float *rec_matrices_ptr = rec_matrices.buf;
  float *con_matrices_ptr = con_matrices.buf;
  float *con_scores_ptr = con_scores.buf;
  float *rec_scores_ptr = rec_scores.buf;
  int *rec_lengths_ptr = rec_lengths.buf;
  int *con_lengths_ptr = con_lengths.buf;
  // if passes as a long (mpi on windows does it for some reason)
  // make an int array anf copy information in and typcast each entry
  // just do .buf in the function calls

  int forward_offset = get_forward_offset(0, rec_lengths_ptr, num_rec);
  int reverse_offset = get_reverse_offset(0, rec_lengths_ptr, num_rec);
  int num_alignments = len_seq - forward_offset - reverse_offset;
  float *score_matrix =
      PyMem_Calloc(num_alignments * num_rec, sizeof(*rec_matrices_ptr));
  fill_matrix(seq, len_seq, rec_matrices_ptr, rec_lengths_ptr, num_rec,
              score_matrix, num_alignments);

  if (num_rec == 1) {
    con_lengths_ptr[0] =
        max_index(score_matrix, len_seq - forward_offset - reverse_offset);
    rec_scores_ptr[0] = score_matrix[con_lengths_ptr[0]];
    con_scores_ptr[0] = 0.00;
  } else {
    traceback(score_matrix, num_alignments, con_matrices_ptr, rec_lengths_ptr,
              num_rec, len_seq, con_scores_ptr, rec_scores_ptr,
              con_lengths_ptr);
  }

  // printf("last in main\n");
  PyMem_Free(score_matrix);
  Py_INCREF(Py_None);
  result = Py_None;
  matrix_converter(NULL, &rec_matrices);
  matrix_converter(NULL, &rec_lengths);
  matrix_converter(NULL, &con_matrices);
  matrix_converter(NULL, &rec_scores);
  matrix_converter(NULL, &con_scores);
  matrix_converter(NULL, &con_lengths);
  // printf("Made it to returning org of size %i\n", num_rec);
  return result;
}

static struct PyMethodDef methods[] = {
    {
        "calculate",
        (PyCFunction)py_calculate,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR(calculate__doc__),
    },
    {NULL, NULL, 0, NULL} // sentinel
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_multiplacement",
    PyDoc_STR("Fast calculations involving multiple connected PSSMs"),
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL};

PyObject *PyInit__multiplacement(void) { return PyModule_Create(&moduledef); }
