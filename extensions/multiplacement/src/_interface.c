#include "_interface.h"
static int matrix_converter(PyObject *object, void *address) {
  const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
  char datatype;
  int dim;
  Py_buffer *view = address;
  if (object == NULL)
    goto exit;
  if (PyObject_GetBuffer(object, view, flags) == -1) {
    PyErr_SetString(PyExc_RuntimeError,
                    "position-weight matrix is not an array");
    return 0;
  }

  dim = view->shape[0]; 
  datatype = view->format[0];
  //printf("%c : %i \n", datatype, dim);
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
  const char *rec_types;
  static char *kwlist[] = {
      "sequence", "rec_types", "rec_matrices", "rec_lengths", "con_matrices",
      "rec_scores", "con_scores",   "con_lengths", "max_length", "bin_freqs", 
      "bin_edges", "num_bins", NULL};
  Py_ssize_t len_seq;
  Py_ssize_t num_rec;
  Py_ssize_t max_length;
  PyObject *result = Py_None;
  Py_buffer rec_matrices;
  Py_buffer con_matrices;
  Py_buffer rec_lengths;
  Py_buffer rec_scores;

  Py_buffer con_scores;
  Py_buffer con_lengths;

  Py_buffer bin_freqs;
  Py_buffer bin_edges;
  Py_buffer num_bins;

  rec_matrices.obj = NULL;
  con_matrices.obj = NULL;
  rec_lengths.obj = NULL;
  rec_scores.obj = NULL; 

  con_scores.obj = NULL;
  con_lengths.obj = NULL;

  bin_freqs.obj = NULL;
  bin_edges.obj = NULL;
  num_bins.obj = NULL;
  if (!PyArg_ParseTupleAndKeywords(
          args, keywords, "y#y#O&O&O&O&O&O&iO&O&O&", kwlist, 
          &seq, &len_seq,
          &rec_types, &num_rec,
          matrix_converter, &rec_matrices, 
          matrix_converter, &rec_lengths,
          matrix_converter, &con_matrices, 
          matrix_converter, &rec_scores,
          matrix_converter, &con_scores, 
          matrix_converter, &con_lengths, 
          &max_length,
          matrix_converter, &bin_freqs,
          matrix_converter, &bin_edges,
          matrix_converter, &num_bins))
    return NULL;

  // sequence:     DNA sequence used for placement
  // rec_matrices: one dimensional flattened representation of the
  //               scoring matrices for all recognizers
  // rec_lengths:  length of each recognizer (the number of columns)
  // con_matrices: one dimensional flattened representation of pre-computed
  //               scores for each connector for each gap length
  // rec_scores:   buffer used to store our calculated scores for recognizers
  // con_scores:   buffer used to score our calculated scores for connectors
  // con_lengths:  buffer used to store the length of each connector for our
  // placment
  double *rec_matrices_ptr = rec_matrices.buf;
  double *con_matrices_ptr = con_matrices.buf;

  double *con_scores_ptr = con_scores.buf;
  double *rec_scores_ptr = rec_scores.buf;
  int *rec_lengths_ptr = rec_lengths.buf;
  int *con_lengths_ptr = con_lengths.buf;

  double *bin_freqs_ptr = bin_freqs.buf;
  double *bin_edges_ptr = bin_edges.buf;
  int *num_bins_ptr = num_bins.buf;

  if (con_matrices.shape[0] == (num_rec - 1) * 2){
    max_length = 0;
  }

  Organism org;
  parse_org(&org, 
            rec_matrices_ptr, 
            rec_lengths_ptr, 
            bin_freqs_ptr, 
            bin_edges_ptr, 
            num_bins_ptr,
            rec_types, 
            num_rec, 
            con_matrices_ptr, 
            max_length);
  //print_org(&org);
  place_org(&org, 
             seq, len_seq, 
             rec_scores_ptr, 
             con_scores_ptr, 
             con_lengths_ptr);

  free(org.recs);
  free(org.cons);
  Py_INCREF(Py_None);
  result = Py_None;
  matrix_converter(NULL, &rec_matrices);
  matrix_converter(NULL, &rec_lengths);
  matrix_converter(NULL, &con_matrices);

  matrix_converter(NULL, &rec_scores);
  matrix_converter(NULL, &con_scores);
  matrix_converter(NULL, &con_lengths);

  matrix_converter(NULL, &bin_freqs);
  matrix_converter(NULL, &bin_edges);
  matrix_converter(NULL, &num_bins);
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
