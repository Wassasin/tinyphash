#pragma once

#include <pnglite.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint8_t *read_png(const char *path, size_t width, size_t height) {
  png_t png;
  bool valid = true;
  int err = png_open_file_read(&png, path);
  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to read file %s\n", path);
    return NULL;
  }

  png_print_info(&png);
  if (png.color_type != PNG_GREYSCALE) {
    fprintf(stderr, "File is not greyscale\n");
    valid = false;
  }

  if (png.width != width || png.height != height) {
    fprintf(stderr, "File is not %zux%zu\n", width, height);
    valid = false;
  }

  if (png.bpp != 1) {
    fprintf(stderr, "File is not grayscale with bytes per pixel = 1\n");
    valid = false;
  }

  uint8_t *buf = NULL;

  if (valid) {
    buf = (uint8_t *)calloc(png.width * png.height * png.bpp + png.height,
                            sizeof(uint8_t));
    err = png_get_data(&png, buf);
    if (err != PNG_NO_ERROR) {
      fprintf(stderr, "Failed to process file buffer\n");
      valid = false;
    }
  }

  err = png_close_file(&png);
  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to close file %s\n", path);
    valid = false;
    ;
  }

  if (valid) {
    return buf;
  } else {
    free(buf);
    return NULL;
  }
}

static bool write_png(const uint8_t *buf, size_t width, size_t height,
                      const char *path) __attribute__((unused));
static bool write_png(const uint8_t *buf, size_t width, size_t height,
                      const char *path) {
  png_t png;
  bool valid = true;
  int err = png_open_file_write(&png, path);
  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to open file for write %s\n", path);
    return false;
  }

  err =
      png_set_data(&png, width, height, 8, PNG_GREYSCALE, (unsigned char *)buf);

  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to write file %s\n", path);
    valid = false;
  }

  err = png_close_file(&png);
  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to close file %s\n", path);
    return false;
  }

  return valid;
}
