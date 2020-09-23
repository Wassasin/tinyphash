#pragma once

// phash implementation re-written to pure C from
// https://github.com/aetilius/pHash

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef TINYPHASH_DEBUG
#include "png.h"
#endif

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

#ifdef TINYPHASH_DEBUG
static tinyphash_smatrix_t tinyphash_to_discrete(tinyphash_smatrixf_t src) {
  tinyphash_smatrix_t dest = {
      .buf = calloc(src.dim * src.dim, sizeof(uint8_t)),
      .dim = src.dim,
  };

  float min = INFINITY, max = -INFINITY;
  for (size_t i = 0; i < src.dim * src.dim; ++i) {
    min = fminf(min, src.buf[i]);
    max = fmaxf(max, src.buf[i]);
  }

  float range = max - min;
  for (size_t i = 0; i < src.dim * src.dim; ++i) {
    dest.buf[i] = (src.buf[i] - min) / range * 255.;
  }

  return dest;
}

void tinyphash_write_debug(tinyphash_smatrix_t matrix, const char *path) {
  write_png(matrix.buf, matrix.dim, matrix.dim, path);
}

void tinyphash_write_debugf(tinyphash_smatrixf_t matrix, const char *path) {
  tinyphash_smatrix_t projected = tinyphash_to_discrete(matrix);
  tinyphash_write_debug(projected, path);
  free(projected.buf);
}
#else
#define UNUSED(x) (void)(x)
void tinyphash_write_debug(tinyphash_smatrix_t matrix, const char *path) {
  UNUSED(matrix);
  UNUSED(path);
}
void tinyphash_write_debugf(tinyphash_smatrixf_t matrix, const char *path) {
  UNUSED(matrix);
  UNUSED(path);
}
#endif

// Compute the DCT matrix for a given dimensionality.
static tinyphash_smatrixf_t tinyphash_dct_matrix(const size_t N) {
  float *buf = (float *)calloc(N * N, sizeof(float));
  float dv = 1 / sqrt((float)N);
  const float c1 = sqrt(2.0 / (float)N);

  // For y = 0, which is actually *really* significant. Especially for the
  // transpose later on.
  for (size_t x = 0; x < N; x++) {
    buf[x] = dv;
  }

  for (size_t y = 1; y < N; y++) {
    for (size_t x = 0; x < N; x++) {
      buf[y * N + x] =
          c1 * cosf((M_PI / 2. / (float)N) * (float)y * (2. * (float)x + 1.));
    }
  }

  tinyphash_smatrixf_t result = {
      .buf = buf,
      .dim = N,
  };

  return result;
}

// Transpose a matrix.
static tinyphash_smatrixf_t tinyphash_transpose(tinyphash_smatrixf_t matrix) {
  float *buf = (float *)calloc(matrix.dim * matrix.dim, sizeof(float));

  for (size_t y = 0; y < matrix.dim; ++y) {
    for (size_t x = 0; x < matrix.dim; ++x) {
      buf[x * matrix.dim + y] = matrix.buf[y * matrix.dim + x];
    }
  }

  tinyphash_smatrixf_t result = {
      .buf = buf,
      .dim = matrix.dim,
  };

  return result;
}

// Perform correlation filtering.
static tinyphash_smatrixf_t tinyphash_correlate_mean(tinyphash_smatrix_t image,
                                                     size_t mean_kernel_dim) {
  size_t dim = image.dim - mean_kernel_dim + 1;

  float *buf = (float *)calloc(dim * dim, sizeof(float));

  for (size_t y = 0; y < dim; ++y) {
    for (size_t x = 0; x < dim; ++x) {
      float v = 0;

      for (size_t ky = 0; ky < mean_kernel_dim; ++ky) {
        size_t ciy = (y + ky) * image.dim;
        size_t cix = ciy + x;

        for (size_t kx = 0; kx < mean_kernel_dim; ++kx) {
          v += image.buf[cix + kx];
        }
      }

      buf[y * dim + x] = v;
    }
  }

  tinyphash_smatrixf_t result = {
      .buf = buf,
      .dim = dim,
  };

  return result;
}

// Perform matrix multiplication, whilst simultaneously cropping the result to a
// pre-allocated dest. (resulting in less multiplications)
static void tinyphash_multiply_cropped(tinyphash_smatrixf_t a,
                                       tinyphash_smatrixf_t b,
                                       tinyphash_smatrixf_t dest) {
  assert(a.dim == b.dim);

  for (size_t i = 0; i < dest.dim; ++i) {
    for (size_t j = 0; j < dest.dim; ++j) {
      float v = 0.;
      for (size_t k = 0; k < a.dim; ++k) {
        v += a.buf[i * a.dim + k] * b.buf[k * a.dim + j];
      }
      dest.buf[i * dest.dim + j] = v;
    }
  }
}

static int tinyphash_float_cmp(const void *x, const void *y) {
  return *(float *)x < *(float *)y ? -1 : 1;
}

// Compute the median of a sequence
static float tinyphash_median(const float *src, size_t size) {
  float *buf = (float *)calloc(size, sizeof(float));
  memcpy(buf, src, size);
  qsort(buf, size, sizeof(float), tinyphash_float_cmp);
  float res = buf[size / 2 - 1];
  free(buf);
  return res;
}

uint64_t tinyphash_dct_unchecked(uint8_t *data, tinyphash_smatrixf_t dct_matrix,
                                 tinyphash_smatrixf_t dct_transpose) {
  tinyphash_smatrix_t image = {
      .buf = data,
      .dim = TINYPHASH_INPUT_DIM,
  };
  tinyphash_write_debug(image, "image_src.png");
  tinyphash_write_debugf(dct_matrix, "image_matrix.png");
  tinyphash_write_debugf(dct_transpose, "image_transpose.png");

  // Should do convolution, but our kernel is point-symmetric, thus this
  // is OK.
  tinyphash_smatrixf_t image_mean =
      tinyphash_correlate_mean(image, TINYPHASH_KERNEL_DIM);
  tinyphash_write_debugf(image_mean, "image_mean.png");

  tinyphash_smatrixf_t mult = {
      .buf = calloc(TINYPHASH_MATRIX_DIM * TINYPHASH_MATRIX_DIM, sizeof(float)),
      .dim = TINYPHASH_MATRIX_DIM,
  };
  tinyphash_smatrixf_t subsec = {
      .buf = calloc(TINYPHASH_SUBSEC_DIM * TINYPHASH_SUBSEC_DIM, sizeof(float)),
      .dim = TINYPHASH_SUBSEC_DIM,
  };
  tinyphash_multiply_cropped(dct_matrix, image_mean, mult);
  free(image_mean.buf);
  tinyphash_write_debugf(mult, "image_mult.png");
  tinyphash_multiply_cropped(mult, dct_transpose, subsec);
  free(mult.buf);
  tinyphash_write_debugf(subsec, "image_subsec.png");

  float median = tinyphash_median(subsec.buf, subsec.dim * subsec.dim);
  uint64_t result = 0;
  for (size_t i = 0; i < subsec.dim * subsec.dim; i++, result <<= 1) {
    float current = subsec.buf[i];
    if (current > median)
      result |= 0x01;
  }

  free(subsec.buf);

  return result;
}

uint64_t tinyphash_dct_easy(uint8_t *data, size_t width, size_t height) {
  assert(width == TINYPHASH_INPUT_DIM && height == TINYPHASH_INPUT_DIM);
  tinyphash_smatrixf_t dct_matrix = tinyphash_dct_matrix(TINYPHASH_MATRIX_DIM);
  tinyphash_smatrixf_t dct_transpose = tinyphash_transpose(dct_matrix);

  uint64_t result = tinyphash_dct_unchecked(data, dct_matrix, dct_transpose);
  free(dct_matrix.buf);
  free(dct_transpose.buf);

  return result;
}

uint8_t tinyphash_hamming_distance(uint64_t x, uint64_t y) {
  uint8_t hamming = 0;
  for (size_t i = 0; i < 64; ++i) {
    if (((x >> i) & 1) != ((y >> i) & 1)) {
      hamming += 1;
    }
  }
  return hamming;
}

float tinyphash_hamming_distancef(uint64_t x, uint64_t y) {
  return ((float)tinyphash_hamming_distance(x, y)) / 64.;
}
