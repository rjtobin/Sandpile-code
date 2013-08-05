#ifndef PNG_DRAW_H
#define PNG_DRAW_H

#include <stdlib.h>
#include <png.h>

#define PNG_SETJMP_NOT_SUPPORTED

struct 
{
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytepp row_pointers;
  FILE *f;
  int width;
  int height;
} typedef pngd_image;

/*
Initialises a pngd_image structure, and performs the necessary file opening
procedure.
 */
pngd_image* pngd_make_image(const char *fpath, int width, int height);

/*
Draws a pixel to an opened image structure.
 */
void pngd_draw_pixel(pngd_image *img, int x, int y, char r, char g, char b);

/*
Close the image and free the memory.
 */
void pngd_finalise(pngd_image *img);

#endif
