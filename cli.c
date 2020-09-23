#include "tinyphash.h"

#include "png.h"
#include "pnglite.h"
int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Please provide a (single) filename\n");
    return 1;
  }

  const char *path = argv[1];

  int err = png_init(0, 0);
  if (err != PNG_NO_ERROR) {
    fprintf(stderr, "Failed to init libpnglite\n");
    return 2;
  }

  uint8_t *buf = read_png(path, TINY_PHASH_BUF_DIM, TINY_PHASH_BUF_DIM);

  if (buf == NULL) {
    return 1;
  }

  uint64_t phash =
      tinyphash_dct_easy(buf, TINY_PHASH_BUF_DIM, TINY_PHASH_BUF_DIM);

  printf("phash: 0x%lx\n", phash);

  free(buf);
}
