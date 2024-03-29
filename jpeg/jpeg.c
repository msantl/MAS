#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define BLOCK_SIZE 8
#define BUFF_LEN 256

const double MOJ_PI = 3.1415926535897932384626433832795028841971693993751058209;

typedef struct rgb_pixel_t { unsigned char rgb[3]; } rgb_pixel;
typedef struct rgb_image_t {
    rgb_pixel* colors;
    int w, h;
} image;

typedef struct ycbcr_pixel_t { double ycbcr[3]; } ycbcr_pixel;
typedef struct ycbcr_block_t {
    ycbcr_pixel colors[BLOCK_SIZE][BLOCK_SIZE];
} ycbcr_block;

typedef struct ycbcr_pixel_int_t { int ycbcr[3]; } ycbcr_pixel_int;
typedef struct ycbcr_block_int_t {
    ycbcr_pixel_int colors[BLOCK_SIZE][BLOCK_SIZE];
} ycbcr_block_int;

typedef struct ycbcr_array_t {
    ycbcr_pixel_int colors[BLOCK_SIZE * BLOCK_SIZE];
} ycbcr_array;

image* img_allocate(int w, int h) {
    image* img = (image *)malloc(sizeof(image));
    assert(img != NULL);
    img->colors = (rgb_pixel *)malloc(w * h * sizeof(rgb_pixel));
    assert(img->colors != NULL);
    img->w = w;
    img->h = h;
    return img;
}

void img_free(image* img) {
    free(img->colors);
    free(img);
    return;
}

/* usage */
const char usage[] = "Multimedijske arihtekture i sustavi\n"
                     "+---------------------------------+\n"
                     "|         1. domaca zadaca        |\n"
                     "+---------------------------------+\n"
                     "./jpeg <slika.ppm <oznaka_bloka>   \n";

/* luminanace quantization table K.1 */
int k1_table[BLOCK_SIZE][BLOCK_SIZE] = {{16, 11, 10, 16, 24, 40, 51, 61},
                                        {12, 12, 14, 19, 26, 58, 60, 55},
                                        {14, 13, 16, 24, 40, 57, 69, 56},
                                        {14, 17, 22, 29, 51, 87, 80, 62},
                                        {18, 22, 37, 56, 68, 109, 103, 77},
                                        {24, 35, 55, 64, 81, 104, 113, 92},
                                        {49, 64, 78, 87, 103, 121, 120, 101},
                                        {72, 92, 95, 98, 112, 100, 103, 99}};

/* chrominance quantization table K.2 */
int k2_table[BLOCK_SIZE][BLOCK_SIZE] = {{17, 18, 24, 47, 99, 99, 99, 99},
                                        {18, 21, 26, 66, 99, 99, 99, 99},
                                        {24, 26, 56, 99, 99, 99, 99, 99},
                                        {47, 66, 99, 99, 99, 99, 99, 99},
                                        {99, 99, 99, 99, 99, 99, 99, 99},
                                        {99, 99, 99, 99, 99, 99, 99, 99},
                                        {99, 99, 99, 99, 99, 99, 99, 99},
                                        {99, 99, 99, 99, 99, 99, 99, 99}};

/* helper wrapper functions */
inline int Round(double x) {
    return (x < 0) ? (int)(x - 0.5) : (int)(x + 0.5);
}

inline double Cos(int i, int u) {
    return cos(((2. * i +  1.) * u * MOJ_PI) / (2. * BLOCK_SIZE) );
}

inline void ConvertRGBtoYCbCr(int R, int G, int B,
                              double* Y, double *Cb, double *Cr) {

    *Y = 0.299 * R + 0.587 * G + 0.114 * B;
    *Cb = -0.1687 * R + -0.3313 * G + 0.5 * B + 128;
    *Cr = 0.5 * R + -0.4187 * G + -0.0813 * B + 128;
    return;
}

/* load ppm image from file */
image* load_ppm_image(FILE *f){
    char buf[BUFF_LEN], *t;
    unsigned int w, h, maxval;
    int r;
    image* img;

    assert(f != NULL);

    t = fgets(buf, BUFF_LEN, f);
    assert(t != NULL);
    assert(strncmp(buf, "P6\n", 3) == 0);

    do {
        t = fgets(buf, BUFF_LEN, f);
        assert(t != NULL);
    } while (strncmp(buf, "#", 1) == 0 );

    r = sscanf(buf, "%u %u", &w, &h);
    assert(r == 2);

    r = fscanf(f, "%u", &maxval);
    assert(r == 1);     assert(maxval < 256);

    fseek(f, 1, SEEK_CUR);

    img = img_allocate(w, h);

    r = fread(img->colors, sizeof(rgb_pixel) , w * h, f);
    assert(r == w * h);

    return img;
}

/* RGB to YCbCr conversion */
void RGBtoYCbCr(image* img, int blk_id, ycbcr_block* blk) {
    int i, j, k;
    double Y, Cb, Cr;

    int row = BLOCK_SIZE * (blk_id / (img->w / BLOCK_SIZE));
    int col = BLOCK_SIZE * (blk_id % (img->w / BLOCK_SIZE));

    for (i = 0; i < BLOCK_SIZE; ++i) {
        for (j = 0; j < BLOCK_SIZE; ++j) {
            k = (row + i) * img->w + (col + j);
            ConvertRGBtoYCbCr(img->colors[k].rgb[0],
                              img->colors[k].rgb[1],
                              img->colors[k].rgb[2],
                              &Y, &Cb, &Cr);

            blk->colors[i][j].ycbcr[0] = Y;
            blk->colors[i][j].ycbcr[1] = Cb;
            blk->colors[i][j].ycbcr[2] = Cr;
        }
    }
    return;
}

/* YCbCr color shift */
void YCbCr_color_shift(ycbcr_block* blk) {
    int i, j;

    for (i = 0; i < BLOCK_SIZE; ++i) {
        for (j = 0; j < BLOCK_SIZE; ++j) {
            blk->colors[i][j].ycbcr[0] -= 128.;
            blk->colors[i][j].ycbcr[1] -= 128.;
            blk->colors[i][j].ycbcr[2] -= 128.;
        }
    }
    return;
}

/* discrete cosine transformation */
void DCT(ycbcr_block* in_blk,
         ycbcr_block* out_blk) {
    int i, j, k, u, v;
    double sum, Cu, Cv;

    for (k = 0; k < 3; ++k) {
        for (u = 0; u < BLOCK_SIZE; ++u) {
            for (v = 0; v < BLOCK_SIZE; ++v) {
                sum = 0.0;

                if (u == 0) { Cu = sqrt(0.5); } else { Cu = 1.0; }
                if (v == 0) { Cv = sqrt(0.5); } else { Cv = 1.0; }

                for (i = 0; i < BLOCK_SIZE; ++i) {
                    for (j = 0; j < BLOCK_SIZE; ++j) {
                        sum += in_blk->colors[i][j].ycbcr[k] *
                               Cos(i, u) *
                               Cos(j, v);
                    }
                }

                out_blk->colors[u][v].ycbcr[k] = (2. / BLOCK_SIZE) *
                                                 Cu * Cv *
                                                 sum;
            }
        }
    }
    return;
}

/* quantize YCbCr block */
void quantize_YCbCr_block(ycbcr_block* blk, ycbcr_block_int *res) {
    int i, j, k;

    for (k = 0; k < 3; ++k) {
        for (i = 0; i < BLOCK_SIZE; ++i) {
            for (j = 0; j < BLOCK_SIZE; ++j) {
                int coef;
                /* use table K1 for Y */
                if (k == 0) { coef = k1_table[i][j]; }
                /* use table K2 for Cb and Cr */
                else { coef = k2_table[i][j]; }

                res->colors[i][j].ycbcr[k] =
                    Round(blk->colors[i][j].ycbcr[k] / coef);
            }
        }
    }

    return;
}

/* reorganize blocks into arrays (zig-zag) */
void zig_zag(ycbcr_block_int* blk, ycbcr_array* arr) {
    int i, j, k, direction, start, end;

    for (k = 0; k < 3; ++k) {
        i = 0, j = 0;
        direction = 1;
        start = 0, end = BLOCK_SIZE * BLOCK_SIZE;

        /* direction == 1 -> moving up-right */
        /* direction == -1 -> moving down-left */

        while (start < end) {
            arr->colors[start++].ycbcr[k] =
                blk->colors[i][j].ycbcr[k];
            arr->colors[--end].ycbcr[k] =
                blk->colors[BLOCK_SIZE - i - 1][BLOCK_SIZE - j - 1].ycbcr[k];

            i -= direction;
            j += direction;

            if (direction == 1) {
                if (i < 0) { i = 0; direction = -1; }
                if (j >= BLOCK_SIZE) { j = BLOCK_SIZE - 1; direction = -1; }

            } else if (direction == -1) {
                if (j < 0) { j = 0; direction = 1; }
                if (i >= BLOCK_SIZE) { i = BLOCK_SIZE - 1; direction = 1; }
            }
        }
    }

    return;
}

/* helper print block */
void __print_ycbcr_block(ycbcr_block_int* blk) {
    int i, j, k;
    for (k = 0; k < 3; ++k) {
        for (i = 0; i < BLOCK_SIZE; ++i) {
            for (j = 0; j < BLOCK_SIZE; ++j) {
                printf("%d ", blk->colors[i][j].ycbcr[k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return;
}

/* helper print array */
void __print_ycbcr_array(ycbcr_array* arr) {
    int i, k;
    for (k = 0; k < 3; ++k) {
        printf("Array %d: ", k);
        for (i = 0; i < BLOCK_SIZE * BLOCK_SIZE; ++i) {
            printf("%d ", arr->colors[i].ycbcr[k]);
        }
        printf("\n");
    }
    return;
}

int main(int argc, char** argv) {
    int blk_id;
    FILE* f;
    image* img;
    ycbcr_block blk, dct;
    ycbcr_block_int res;
    ycbcr_array arr;

    if (argc != 3 ) {
        printf("%s", usage);
        exit(1);
    } else {
        sscanf(argv[2], "%d", &blk_id);
    }

    /* load .ppm file */
    f = fopen(argv[1], "rb");
    img = load_ppm_image(f);
    fclose(f);

    /* convert from RGB to YCbCr */
    RGBtoYCbCr(img, blk_id, &blk);

    /* shift YCbCr values */
    YCbCr_color_shift(&blk);

    /* transform using 2D-DCT */
    DCT(&blk, &dct);

    /* quantization (K.1 for Y and K.2 for Cb and Cr) */
    quantize_YCbCr_block(&dct, &res);

    /* zig-zag reorganization */
    zig_zag(&res, &arr);

    __print_ycbcr_array(&arr);

    img_free(img);
    return 0;
}
