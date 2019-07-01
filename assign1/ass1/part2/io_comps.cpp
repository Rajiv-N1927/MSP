#include "image_comps.h"
/* ========================================================================= */
/*                 Implementation of `filter_manager' functions              */
/* ========================================================================= */
filter_manager::filter_manager()
{ extent = 0;  taps = nullptr; }

filter_manager::~filter_manager()
{ if (taps != nullptr) delete[] taps; }

void filter_manager::init(float * filter, int extent) {
    int dim = 2*extent+1;
    int num_taps = dim*dim;
    taps = new float[num_taps];
    int r, c;
    for ( r = 0; r < dim; r++ ) {
        float *curpos = filter + dim*r;
        float *taps_pos = taps + dim*r;
        for ( c = 0; c < dim; c++ ) {
            taps_pos[c] = curpos[c];
        }
    }
    this->extent = extent;
}
//Returns a normalized filter
void filter_manager::normalize_filter() {
    //Sum up all the taps then divide each tap by the sum
    float sum = 0;
    int dim = 2*extent + 1;
    int num_taps = dim*dim;
    int r, c;
    for ( r = 0; r < dim; r++ ) {
        float *curpos = taps + dim*r;
        for ( c = 0; c < dim; c++ ) {
            sum += curpos[c];
        }
    }

    for ( r = 0; r < dim; r++ ) {
        float *curpos = taps + dim*r;
        for ( c = 0; c < dim; c++ ) {
            curpos[c] /= sum;
        }
    }

}

void filter_manager::mirror_filter() {
    //Go to the centre of the filter and mirror it
    float sum = 0;
    int dim = 2*extent + 1;
    int num_taps = dim*dim;
    int r, c;
    float * mirr_filt = new float[num_taps];
    int r1, r2, c1, c2;

    for ( r1 = 0, r2 = dim; r1 < dim && r2 >= 1; r1++, r2-- ) {
        float *filt_pos = taps + r1*dim;
        float *mirr_pos = mirr_filt + r2*dim;
        for ( c1 = 0, c2 = dim-1; c1 < dim && c2 > 0; c1++, c2-- ) {
            mirr_pos[c2] = filt_pos[c1];
        }
    }
    delete[] taps;
    taps = mirr_filt;
}
/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

my_image_comp::my_image_comp()
  { width = height = stride = border = 0;  handle = buf = nullptr; }
my_image_comp::~my_image_comp()
  { if (handle != nullptr) delete[] handle; }
void my_image_comp::init(int height, int width, int border)
  {
    this->width = width;  this->height = height;  this->border = border;
    stride = width + 2*border;
    if (handle != nullptr)
      delete[] handle; // Delete mem allocated by any previous `init' call
    handle = new float[stride*(height+2*border)];
    buf = handle + (border*stride) + border;
  }

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
}

void my_image_comp::perform_billinear_interpolation(my_image_comp *im, float sf)
{
    int row, col;
    for ( row = 0; row < height; row++ ) {
        float *out_pos = this->buf + row*this->stride;
        for ( col = 0; col < width; col++ ) {
            int posx = col/sf;
            int posy = row/sf;
            float delx = col/sf - posx;
            float dely = row/sf - posy;
            float *p1 = im->buf + posy*im->stride + posx;
            float *p2 = p1 + 1;
            float *p3 = p1 + im->stride;
            float *p4 = p3 + 1;
            //Do the 3 step interpolation
            float a = *p1 + delx*(*p2 - *p1);
            float b = *p3 + delx*(*p4 - *p3);
            float c = a + dely*(b-a);
            out_pos[col] = (int)(c+0.5);
        }
    }
}
