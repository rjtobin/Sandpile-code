#include "png_draw.h"

pngd_image* pngd_make_image(const char *fpath, int width, int height)
{
  pngd_image *img = (pngd_image*) malloc(sizeof(pngd_image)); 

  img->height = height;
  img->width = width;

  img->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0,0,0);
  if (!img->png_ptr)
    return 0;
  img->info_ptr = png_create_info_struct(img->png_ptr);
  if (!img->info_ptr)
  {
    png_destroy_write_struct(&(img->png_ptr), (png_infopp)NULL);
    return 0;
  }

  FILE *fp = fopen(fpath, "wb");
  png_init_io(img->png_ptr, fp);

  png_set_IHDR(img->png_ptr, img->info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  img->row_pointers = (png_bytepp) calloc(width, sizeof(png_bytep));
  int i,j;
  for(i=0; i<width; i++)
    img->row_pointers[i] = (png_bytep) malloc(img->info_ptr->rowbytes);

  for(i=0; i<width; i++)
    for(j=0; j < img->info_ptr->rowbytes; j++)
      img->row_pointers[i][j] = (png_byte) 0;

  return img;
}

void pngd_draw_pixel(pngd_image *img, int y, int x, char r, char g, char b)
{
  (img->row_pointers)[x][y*3] = r;
  (img->row_pointers)[x][y*3+1] = g;
  (img->row_pointers)[x][y*3+2] = b;
}

void pngd_finalise(pngd_image *img)
{
  png_write_info(img->png_ptr, img->info_ptr);
  png_write_image(img->png_ptr, img->row_pointers);
  png_write_end(img->png_ptr, img->info_ptr);

  int i,j;
  for(i=0; i<img->width; i++)
    free(img->row_pointers[i]);
  free(img->row_pointers);

  png_destroy_write_struct(&(img->png_ptr), &(img->info_ptr));
  free(img);
}
