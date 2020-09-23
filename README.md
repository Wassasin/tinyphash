# tiny-phash
A header-only pure-C implementation of (DCT) phash inspired by [aetilius/pHash](https://github.com/aetilius/pHash). Confirmed to be able to run on an ESP32 (with external RAM).

## Dependencies
This library has no dependencies. However, to run the CLI you require `libpnglite` and `zlib`.

For Ubuntu 18.04 this entails:
```bash
apt install libpnglite0 build-essential clang
```

## How to build
*Note:* with the libpnglite0 on Ubuntu 18.04 reading of files did not work. I rebuilt the library and linked statically, and hence it just worked. Generally you should:

```bash
clang ./cli.c -O2 -g -o tinyphash -lpnglite -lz -lm
```

But I do:
```bash
clang ./cli.c ./libpnglite.a -O2 -g -o tinyphash -lz -lm
```

## How to run
For now this library **requires** grayscale `uint8_t` images resampled to 38x38 pixels.
The commandline utility requires grayscale PNG images. You can quickly invoke these for a folder of JPEG files using Imagemagick:
```bash
convert -resize 38x38! -colorspace Gray '*.jpg' %d.png
```

And then run:

```bash
./tinyphash 0.png
```

## How to use the library
The library itself only uses libc as a dependency. You can invoke the easy API using: (note that any other size than 38x38 will not be accepted at this time)
```c
uint64_t phash = tinyphash_dct_easy(buf, 38, 38);
```
I recommend checking out the [imageresampler](https://github.com/rwohleb/imageresampler) library to sample your image down. Probably you would like to reuse the DCT transposed matrix between runs, by allocating it statically:

```c
static tinyphash_smatrixf_t *dct_matrix_transposed = NULL;
if (dct_matrix_transposed == NULL) {
    dct_matrix_transposed = malloc(sizeof(tinyphash_smatrixf_t));
    *dct_matrix_transposed = tinyphash_dct_precompute();
}

uint64_t result = tinyphash_dct_unchecked(data, *dct_matrix_transposed);
```
