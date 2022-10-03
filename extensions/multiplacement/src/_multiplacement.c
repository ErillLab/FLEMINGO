#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
const int NUM_BASES = 4;

int getForwardOffset(int index, int cols[], int numPSSM) {
  // finds the first possible possition for a pssm
  // based on number of columns of preceding pssms

  int offset = 0;
  for (int i = 0; i < index; i++) {
    offset += cols[i];
  }
  return offset;
}

int getReverseOffset(int index, int cols[], int numPSSM) {
  // finds the last possible possition for a pssm
  // based on number of columns of subsequenct pssms

  int offset = 0;
  for (int i = numPSSM - 1; i >= index; i--) {
    offset += cols[i];
  }

  // this is important to get the right number of possible alignments
  return offset - 1;
}

int maxIndex(float *arr, int size) {
  int maxIndex = 0;
  for (int i = 0; i < size; i++) {

    if (arr[i] > arr[maxIndex]) {
      maxIndex = i;
    }
  }
  return maxIndex;
}

void calculatePlacements(float *scoreMatrix, float *gapMatrix, int *cols,
                         int numPSSM, int lenSequence, float *gapScores,
                         float *scores, int *gaps) {
  int numAlignments = 0;
  int gapLength = 0;
  int forwardOffset = getForwardOffset(0, cols, numPSSM);
  int reverseOffset = getReverseOffset(0, cols, numPSSM);
  numAlignments = lenSequence - forwardOffset - reverseOffset;

  // number of total alignments by number of pssms
  // first index in each column holds current max score for that index
  // all other indices hold gap lengths that got that alignment
  float *alignments = calloc(numAlignments, sizeof(*scoreMatrix));
  int *gapAlignments = calloc(numAlignments * (numPSSM - 1), sizeof(*gaps));
  float *tempMax = calloc(numAlignments, sizeof(*scoreMatrix));
  int *tempGaps = calloc(numAlignments, sizeof(*gaps));

  // start with first row as our current max
  for (int i = 0; i < numAlignments; i++) {
    alignments[i] = scoreMatrix[i];
    tempMax[i] = scoreMatrix[i];
  }

// time complexity should be (numConnectors * numAlignments * sum(0 1 2 ...
// numAlignments))

// for each connector (number of pssms - 1) populate alignments with
// current maximum score for that index

// overview: the algorithm iterates over the next row in the scoresMatrix (first
// row is our starting cumulative score), for each index, we compare the scores
// obtained by summing the PSSM score at that index with each previous
// cumulative alignment score for k <= j, when k == j, the gap length is 0
  for (int i = 1; i < numPSSM; i++) {
    // printf("Calculating alignments for PSSM %i\n", i);

    for (int j = 0; j < numAlignments; j++) {

      for (int k = 0; k <= j; k++) {

        // every column before or equal to the current column is a valid
        // alignment. Scores for each alignment
        // gapLength = difference between column normalized j and k
        gapLength = j - k;

        if (k == 0) {
          tempMax[j] = alignments[k] +
                       gapMatrix[(i - 1) * lenSequence + gapLength] +
                       scoreMatrix[i * numAlignments + j];
          tempGaps[j] = gapLength;
        } else {
          if (tempMax[j] < alignments[k] +
                               gapMatrix[(i - 1) * lenSequence + gapLength] +
                               scoreMatrix[i * numAlignments + j]) {
            tempMax[j] = alignments[k] +
                         gapMatrix[(i - 1) * lenSequence + gapLength] +
                         scoreMatrix[i * numAlignments + j];
            tempGaps[j] = gapLength;
          }
        }
      }
    }

    for (int l = 0; l < numAlignments; l++) {
      alignments[l] = tempMax[l];
      gapAlignments[(i - 1) * numAlignments + l] = tempGaps[l];
      tempMax[l] = -INFINITY;
      tempGaps[l] = 0;
    }
  }

  // finding gap lengths to ref orma alignments
  // starting from the index of the greatest score in alignments
  // trace back of gap alignments is conducted by subtracting the value
  // at the max score index from the index and going up a row
  // this is repeated until we know all the gap lengths in the alignment
  int index = maxIndex(alignments, numAlignments);
  gaps[numPSSM - 1] = gapAlignments[(numPSSM - 2) * numAlignments + index];

  // the cumulative best score is written into the last index of the scores
  // array
  scores[numPSSM] = alignments[index];

  for (int i = numPSSM - 3; i >= 0; i--) {
    gaps[i + 1] = gapAlignments[i * numAlignments + index - gaps[i + 2]];
    index -= gaps[i + 2];
  }
  gaps[0] = index - gaps[1];

  int gapOffset = 0;
  for (int i = 0; i < numPSSM - 1; i++) {
    gapScores[i] = gapMatrix[i * lenSequence + gaps[i + 1]];
  }

  // scores for each PSSM are filled by iterating over the score matrix
  // and using the appropriate gap lengths as a cumulative offset
  for (int i = 0; i < numPSSM; i++) {
    gapOffset += gaps[i];
    scores[i] = scoreMatrix[i * numAlignments + gapOffset];
  }

  free(alignments);
  free(gapAlignments);
  free(tempMax);
  free(tempGaps);
  free(scoreMatrix);
}

float *calculateScores(const char seq[], int lenSequence, float pssm[],
                       int cols[], int numPSSM) {
  // length of the sequence by number of pssms

  float score = 0;
  int forwardOffset = 0;
  int reverseOffset = 0;
  forwardOffset = getForwardOffset(0, cols, numPSSM);
  reverseOffset = getReverseOffset(0, cols, numPSSM);
  int numAlignments = lenSequence - forwardOffset - reverseOffset;
  float *PSSMScores = calloc(numAlignments * numPSSM, sizeof(*pssm));

  // pre computes alignments of each pssm at each possible position
  // i = current pssm
  // j = starting position on sequence for computing score
  // k = current column in pssm
  // time complexity should be (numPSSMs * numAignments * avgLengthPSSMs)

  for (int i = 0; i < numPSSM; i++) {
    // printf("Calculating scores for PSSM %i\n", i);
    // these functions are unnessesary, change to using += and -= cols[i]
    forwardOffset = getForwardOffset(i, cols, numPSSM);
    reverseOffset = getReverseOffset(i, cols, numPSSM);

    for (int j = forwardOffset; j < lenSequence - reverseOffset; j++) {
      score = 0;
      // printf("\nstarting positon %i\n", j);
      for (int k = 0; k < cols[i]; k++) {
        switch (seq[j + k]) {
        case 'A':
        case 'a':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forwardOffset + k) * 4 + 0]);
          score += pssm[(forwardOffset + k) * 4 + 0];
          break;
        case 'G':
        case 'g':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forwardOffset + k) * 4 + 1]);
          score += pssm[(forwardOffset + k) * 4 + 1];
          break;
        case 'C':
        case 'c':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forwardOffset + k) * 4 + 2]);
          score += pssm[(forwardOffset + k) * 4 + 2];
          break;
        case 'T':
        case 't':
          // printf("   Score for pssm %i column %i with %c: %0.2f\n", i, k,
          // seq[j + k], pssm[(forwardOffset + k) * 4 + 3]);
          score += pssm[(forwardOffset + k) * 4 + 3];
          break;
        }
      }
      PSSMScores[(i * numAlignments) + j - forwardOffset] = score;
    }
  }

  return PSSMScores;
}

static int float_array_converter(PyObject *object, void *address) {
  const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
  char datatype;
  Py_buffer *view = address;

  if (object == NULL)
    goto exit;
  if (PyObject_GetBuffer(object, view, flags) == -1) {
    PyErr_SetString(PyExc_RuntimeError, "PSSMs is not an array");
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

  if (datatype != 'f') {
    PyErr_Format(
        PyExc_RuntimeError,
        "PSSMs matrix data format incorrect data format ('%c', expected 'f')",
        datatype);
    goto exit;
  }

  if (view->ndim != 1) {
    PyErr_Format(PyExc_RuntimeError,
                 "PSSMs matrix has incorrect rank (%d expected 1)", view->ndim);
    goto exit;
  }

  return Py_CLEANUP_SUPPORTED;

exit:
  PyBuffer_Release(view);
  return 0;
}

static int int_array_converter(PyObject *object, void *address) {
  const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
  char datatype;
  Py_buffer *view = address;

  if (object == NULL)
    goto exit;
  if (PyObject_GetBuffer(object, view, flags) == -1)
    return 0;
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

  if (datatype != 'i' && datatype != 'l') {
    PyErr_Format(PyExc_RuntimeError,
                 "columns array has incorrect data format ('%c', expected 'i' or 'l')",
                 datatype);
    goto exit;
  }

  if (view->ndim != 1) {
    PyErr_Format(PyExc_ValueError,
                 "columns array has incorrect rank (%d expected 1)",
                 view->ndim);
    goto exit;
  }
  return Py_CLEANUP_SUPPORTED;

exit:
  PyBuffer_Release(view);
  return 0;
}

static char calculate__doc__[] =
    "    _multiplacement(sequence, PSSMs[], columns[], gap_score_array, "
    "gap_scores, PSSM_scores, gap_lengths)\n"
    "\n"
    "This function calculates the optimal placement for a complex multi PSSM "
    "model\n"
    "along a given sequence and writes information to gap_scores, PSSM_scores "
    "and gap_lengths\n";

static PyObject *py_calculate(PyObject *self, PyObject *args,
                              PyObject *keywords) {
  const char *sequence;
  static char *kwlist[] = {"sequence",  "PSSMs",  "columns", "gapMatrix",
                           "gapScores", "scores", "gaps",    NULL};
  Py_ssize_t numPSSM;
  Py_ssize_t lenSequence;
  PyObject *result = Py_None;

  Py_buffer matrix;
  Py_buffer columns;
  Py_buffer gapMatrix;
  Py_buffer gapScores;
  Py_buffer scores;
  Py_buffer gaps;

  // I don't know why this is done, but it was in _pwm.c so I'll
  // stick to what they did
  matrix.obj = NULL;
  columns.obj = NULL;
  gapMatrix.obj = NULL;
  gapScores.obj = NULL;
  scores.obj = NULL;
  gaps.obj = NULL;

  if (!PyArg_ParseTupleAndKeywords(
          args, keywords, "s#O&O&O&O&O&O&", kwlist, &sequence, &lenSequence,
          float_array_converter, &matrix, int_array_converter, &columns,
          float_array_converter, &gapMatrix, float_array_converter, &gapScores,
          float_array_converter, &scores, int_array_converter, &gaps))
    return NULL;
  numPSSM = columns.shape[0];
  float *PSSMs_ptr = matrix.buf;
  float *gaps_matrix_ptr = gapMatrix.buf;
  float *gap_scores_ptr = gapScores.buf;
  float *scores_ptr = scores.buf;
  int *column_ptr = calloc(numPSSM, sizeof(int)); 
  int *gaps_ptr = calloc(numPSSM, sizeof(int));


  if (columns.itemsize == sizeof(long)){
    long *tempColumns = columns.buf;
    for (int i = 0; i < numPSSM; i++){
      column_ptr[i] = (int)tempColumns[i];
    }
  }else{
    column_ptr = columns.buf;

  }

  if(gaps.itemsize == sizeof(long)){
      long *tempGaps = gaps.buf;
      for (int i = 0; i < numPSSM - 1; i++){
        gaps_ptr[i] = (int)tempGaps[i];
      }
  }else{
    gaps_ptr = gaps.buf;
  }
 
  // just do .buf in the function calls
  float *scoreMatrix =
      calculateScores(sequence, lenSequence, PSSMs_ptr, column_ptr, numPSSM);

  if (numPSSM == 1) {

    int forwardOffset = getForwardOffset(0, column_ptr, numPSSM);
    int reverseOffset = getReverseOffset(0, column_ptr, numPSSM);
    gaps_ptr[0] =
        maxIndex(scoreMatrix, lenSequence - forwardOffset - reverseOffset);
    scores_ptr[0] = scoreMatrix[gaps_ptr[0]];
    gap_scores_ptr[0] = 0.00;
    free(scoreMatrix);

  } else {

    calculatePlacements(scoreMatrix, gaps_matrix_ptr, column_ptr, numPSSM,
                        lenSequence, gap_scores_ptr, scores_ptr, gaps_ptr);
  }

  float_array_converter(NULL, &matrix);
  int_array_converter(NULL, &columns);
  float_array_converter(NULL, &gapMatrix);
  float_array_converter(NULL, &gapScores);
  int_array_converter(NULL, &gaps);

  Py_INCREF(Py_None);
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
