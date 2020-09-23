#include "tinyphash.h"

#include "png.h"
#include "pnglite.h"
int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    fprintf(stderr, "Please provide a (single) filename, and optionally a "
                    "phash in hex (i.e. 0xa7a783838383b3b6)\n");
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

  if (argc == 3) {
    char *end;
    uint64_t original_phash = strtoul(argv[2], &end, 0);
    if (original_phash == 0) {
      fprintf(stderr, "Failed to parse comparison phash\n");
    } else {
      float distance = tinyphash_hamming_distancef(original_phash, phash);
      printf("distance: %.2f\n", distance);
    }
  }

  free(buf);
}
