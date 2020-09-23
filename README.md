# tiny-phash
A pure-C implementation of (DCT) phash inspired by [aetilius/pHash](https://github.com/aetilius/pHash). Confirmed to be able to run on an ESP32 (with external RAM).

## Dependencies
This library has no dependencies. However, to run the CLI you require `libpnglite` and `zlib`.

For Ubuntu 18.04 this entails:
```bash
apt install libpnglite0 build-essential clang
```

## How to build

```bash
clang ./cli.c -O2 -g -o tinyphash -lpnglite -lz -lm
```

## How to run

```bash
./tinyphash <filename>
```
