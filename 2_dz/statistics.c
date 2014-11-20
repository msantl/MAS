#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define BUFF_LEN    256

typedef struct bw_pixel_t { unsigned char bw; } bw_pixel;
typedef struct bw_image_t {
    bw_pixel* colors;
    int w, h;
} image;

unsigned int group_counter[16];

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

int main(int argc, char** argv) {
    FILE* f;
    image *img;
    int i, j, k;

    if (argc != 2) {
        printf("Please specify only the input file!");
        exit(1);
    }

    f = fopen(argv[1], "rb");
    img = load_pgm_image(f);
    fclose(f);

    memset(group_counter, 0, sizeof group_counter);

    for (i = 0; i < img->h; ++i) {
        for (j = 0; j < img->w; ++j) {
            k = img->colors[i * img->w + j].bw;

            group_counter[(k >> 4) & 0xf]++;
        }
    }

    for (i = 0; i < 16; ++i) {
        printf("[%2d] %5d = %.5lf\n",
                i,
                group_counter[i],
                group_counter[i] / (img->w * img->h * 1.));
    }

    img_free(img);

    return 0;
}
