#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#define BUFF_LEN    256
#define BLOCK_SIZE  16

typedef struct bw_pixel_t { unsigned char bw; } bw_pixel;
typedef struct bw_image_t {
    bw_pixel* colors;
    int w, h;
} image;

typedef struct bw_block_t {
    bw_pixel colors[BLOCK_SIZE][BLOCK_SIZE];
} block;

unsigned int group_counter[16];

/* usage */
const char usage[] = "Multimedijske arihtekture i sustavi\n"
                     "+---------------------------------+\n"
                     "|         2. domaca zadaca        |\n"
                     "+---------------------------------+\n"
                     "./matija_santl_dz2.exe <id bloka>  \n";

image* img_allocate(int w, int h) {
    image* img = (image *)malloc(sizeof(image));
    assert(img != NULL);
    img->colors = (bw_pixel *)malloc(w * h * sizeof(bw_pixel));
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

/* load pgm image from file */
image* load_pgm_image(FILE *f){
    char buf[BUFF_LEN], *t;
    unsigned int w, h, maxval;
    int r;
    image* img;

    assert(f != NULL);

    t = fgets(buf, BUFF_LEN, f);
    assert(t != NULL);
    assert(strncmp(buf, "P5\n", 3) == 0);

    while (getc(f) == '#') {
        t = fgets(buf, BUFF_LEN, f);
        assert(t != NULL);
    }

    fseek(f, -1, SEEK_CUR);

    r = fscanf(f, "%u", &w);
    assert(r == 1);

    r = fscanf(f, "%u", &h);
    assert(r == 1);

    r = fscanf(f, "%u", &maxval);
    assert(r == 1);     assert(maxval < 256);

    fseek(f, 1, SEEK_CUR);

    img = img_allocate(w, h);

    r = fread(img->colors, sizeof(bw_pixel) , w * h, f);
    assert(r == w * h);

    return img;
}

void get_pgm_image(image** img, const char* filename) {
    FILE* f;

    f = fopen(filename, "rb");
    *img = load_pgm_image(f);
    fclose(f);

    return;
}

double get_mda(block* curr, block* prev) {
    int i, j;
    double ret = 0;

    for (i = 0; i < BLOCK_SIZE; ++i) {
        for (j = 0; j < BLOCK_SIZE; ++j) {
            ret += fabs(curr->colors[i][j].bw - prev->colors[i][j].bw);
        }
    }

    return ret / (BLOCK_SIZE * BLOCK_SIZE);
}

void get_block_by_x_y(image* img, block* blk, int x, int y) {
    int i, j, k;

    for (i = 0; i < BLOCK_SIZE; ++i) {
        for (j = 0; j < BLOCK_SIZE; ++j) {
            k = (x + i) * img->w + (y + j);

            blk->colors[i][j].bw = img->colors[k].bw;
        }
    }
    return;
}

void print_me_vector(image* curr, image* prev, int blk_id, int off) {
    block blk_curr, blk_prev;
    int i, j;
    int mx, my;
    double best = DBL_MAX, temp;

    int row = BLOCK_SIZE * (blk_id / (curr->w / BLOCK_SIZE));
    int col = BLOCK_SIZE * (blk_id % (curr->w / BLOCK_SIZE));

    get_block_by_x_y(curr, &blk_curr, row, col);

    for (i = -off; i <= off; ++i) {
        for (j = -off; j <= off; ++j) {
            if (row + i < 0 || row + i + BLOCK_SIZE > prev->h) continue;
            if (col + j < 0 || col + j + BLOCK_SIZE > prev->w) continue;

            get_block_by_x_y(prev, &blk_prev, row + i, col + j);
            temp = get_mda(&blk_curr, &blk_prev);
            if (temp < best) {
                best = temp;
                mx = j;
                my = i;
            }
        }
    }

    printf("(%d,%d)\n", mx, my);

    return;
}

int main(int argc, char** argv) {
    int blk_id;
    image *img_prev, *img_curr;

    if (argc != 2) {
        printf("%s", usage);
        exit(1);
    } else {
        sscanf(argv[1], "%d", &blk_id);
    }

    get_pgm_image(&img_prev, "lenna.pgm");
    get_pgm_image(&img_curr, "lenna1.pgm");

    print_me_vector(img_curr, img_prev, blk_id, 16);

    img_free(img_prev);
    img_free(img_curr);
    return 0;
}
