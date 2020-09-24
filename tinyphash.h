#pragma once

// phash implementation re-written to pure C from
// https://github.com/aetilius/pHash

#include <stddef.h>
#include <stdint.h>

#define TINYPHASH_MATRIX_DIM 32
#define TINYPHASH_KERNEL_DIM 7
#define TINYPHASH_INPUT_DIM (TINYPHASH_MATRIX_DIM + TINYPHASH_KERNEL_DIM - 1)
#define TINYPHASH_SUBSEC_DIM 8

typedef struct {
  float *buf;
  size_t dim;
} tinyphash_smatrixf_t;

typedef struct {
  uint8_t *buf;
  size_t dim;
} tinyphash_smatrix_t;

// Compute the DCT matrix for a given dimensionality.
tinyphash_smatrixf_t tinyphash_dct_matrix(const size_t N);

// Transpose a matrix.
tinyphash_smatrixf_t tinyphash_transpose(tinyphash_smatrixf_t matrix);

// Compute the tinyphash fast
//
// @param data a single byte per pixel matrix in grayscale with
// TINYPHASH_INPUT_DIM*TINYPHASH_INPUT_DIM length.
// @param dct_matrix the result of `tinyphash_dct_matrix(TINYPHASH_MATRIX_DIM)`.
// @param dct_transpose the value of dct_matrix passed through
// `tinyphash_transpose`.
uint64_t tinyphash_dct_unchecked(uint8_t *data, tinyphash_smatrixf_t dct_matrix,
                                 tinyphash_smatrixf_t dct_transpose);

// Compute the tinyphash fast and simply.
//
// *Note* You require to pass a width and height
// equal to TINYPHASH_INPUT_DIM for the current implementation.
//
// For more effiency use `tinyphash_dct_unchecked` and re-use the matrix and
// transposed matrix.
//
// @param data a single byte per pixel matrix in grayscale.
// @param width the width of the matrix, must be TINYPHASH_INPUT_DIM.
// @param width the height of the matrix, must be TINYPHASH_INPUT_DIM.
uint64_t tinyphash_dct_easy(uint8_t *data, size_t width, size_t height);

// Compute the absolute hamming distance, with a value between 0 and 64.
uint8_t tinyphash_hamming_distance(uint64_t x, uint64_t y);

// Compute the relative hamming distance, with a value between 0.0 and 1.0.
float tinyphash_hamming_distancef(uint64_t x, uint64_t y);
