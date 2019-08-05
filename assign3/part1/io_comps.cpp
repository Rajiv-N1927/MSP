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
    handle = new int[stride*(height+2*border)];
    buf = handle + (border*stride) + border;
  }

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  int *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  int *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  int *left_edge = buf-border*stride;
  int *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
}


void my_image_comp::perform_symmetric_extension()
{
    /*
        LHS: x[-n1, n2] = x[n1, n2]
        RHS: x[N1-1+n1, n2] = x[N1 - n1, n2]
    */
    int *upper_bound = buf;
    int *lower_bound = buf + (height-1)*stride;
    int *left_edge = buf-border*stride;
    int *right_edge = left_edge + width - 1;
    int r, c;
    //Extending up
    for ( r = 1; r <= border; r++ )
        for ( c = 0; c < width; c++ )
            upper_bound[-r*stride+c] = buf[r*stride+c];
    //Extending down
    for ( r = 1; r <= border; r++ )
        for ( c = 0; c < width; c++ )
            lower_bound[r*stride+c] = lower_bound[-r*stride+c];

    for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
        for (c=1; c <= border; c++)
        {
          left_edge[-c] = left_edge[c];
          right_edge[c] = right_edge[-c];
        }
}

void my_image_comp::perform_point_symmetric_extension() {
  int *upper_bound = buf;
  int *lower_bound = buf + (height-1)*stride;
  int *left_edge = buf-border*stride;
  int *right_edge = left_edge + width - 1;
  int r, c;
  //Extending up
  for ( r = 1; r <= border; r++ )
      for ( c = 0; c < width; c++ )
          upper_bound[-r*stride+c] = buf[r*stride+c];
  //Extending down
  for ( r = 1; r <= border; r++ )
      for ( c = 0; c < width; c++ )
          lower_bound[r*stride+c] = lower_bound[-r*stride+c];

  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
      for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[c];
        right_edge[c] = right_edge[-c];
      }
}
/* Keypoint generator stuff */
kp_generator::kp_generator() {
  tgt = nullptr;
  kp_offsets = nullptr;
  local_vector = nullptr;
  b_extent = sa_extent = delta = 0;
  max_height = max_width = base_offset = num_points = 0;
}
//Image already gets destroyed later on
kp_generator::~kp_generator() {
  if (kp_offsets != nullptr) delete[] kp_offsets;
  if (local_vector != nullptr) delete[] local_vector;
}

void kp_generator::init(my_image_comp * tgt, int b_ext, int sa_ext, int del) {
  this->tgt = tgt;
  printf("Image props: w %d, h %d\n", tgt->width, tgt->height);
  this->b_extent = b_ext;
  this->sa_extent = sa_ext;
  this->delta = del;
  //The base offset is just the search area extent in the x and y direction
  //plus the block extent in the x and y direction
  base_offset = b_ext + sa_ext;
  max_width  = (int)((float)((tgt->width-2*(base_offset))/del));
  max_height = (int)((float)((tgt->height-2*(base_offset))/del));
  //Max number of keypoints
  float w_pts = max_width-base_offset;
  float h_pts = max_height-base_offset;
  float n_pts = max_width * max_height; //To ensure that all points fit in
  kp_offsets   = new mvector[(int)n_pts];
  local_vector = new mvector[(int)n_pts];
  printf("max_w: %d, max_h: %d, base_offset %d, w_pts: %0.2f, h_pts: %0.2f, n_pts %0.2f\n",
    max_width, max_height, base_offset, w_pts, h_pts, n_pts);
}
/*
  Formula: base_offset + delta*(n1, n2)
*/
void kp_generator::generate_kp(my_image_comp * ref) {
  int n1, n2, pos = 0;
  for ( n2 = 0; n2 < max_height; n2++ ) {
    for ( n1 = 0; n1 < max_width; n1++ ) {
      kp_offsets[pos].x = base_offset + n1*delta;
      kp_offsets[pos].y = base_offset + n2*delta;
      sad_vector(ref, kp_offsets[pos].y - b_extent,
        kp_offsets[pos].x - b_extent, pos);
      pos++;
    }
  }
  num_points = pos;
  produce_gv();
  printf("Global vector: %0.5f, %0.5f\n", gx, gy);
}

/*
  Technically a 'private' function. This function runs inside kp_generator
  to produce the vectors as the key points are supplied
*/
void kp_generator::sad_vector(my_image_comp *ref, int start_row, int start_col,
  int pos)
{
  mvector vec, best_vec;
  int block_width  = 2*b_extent + 1;
  int block_height = block_width;
  int sad, best_sad=256*block_width*block_height; //Max size
  //The size of the vec.y represents the area in which the search occurs in the y direction
  for (vec.y=-sa_extent; vec.y <= sa_extent; vec.y++)   {
  //The size of the vec.x represents the area in which the search occurs in the x direction
    for (vec.x=-sa_extent; vec.x <= sa_extent; vec.x++)
      {
        int ref_row = start_row-vec.y;
        int ref_col = start_col-vec.x;
        if ((ref_row < 0) || (ref_col < 0) ||
            ((ref_row+block_height) > ref->height) ||
            ((ref_col+block_width) > ref->width))
          continue; // Translated block not containe within reference frame
        int r, c;
        int *rp = ref->buf + ref_row*ref->stride + ref_col;
        int *tp = tgt->buf + start_row*tgt->stride + start_col;
        for (sad=0, r=block_height; r > 0; r--,
             rp+=ref->stride, tp+=tgt->stride) {
          for (c=0; c < block_width; c++)
            {
              int diff = tp[c] - rp[c];
              sad += (diff < 0)?(-diff):diff;
            }
        }
        if (sad < best_sad)
          {
            best_sad = sad;
            best_vec = vec;
          }
      }
  }

  local_vector[pos] = best_vec;
}

void kp_generator::produce_gv() {
  float vecx = 0;
  float vecy = 0;
  for ( int i = 0; i < num_points; i++ ) {
    vecx += local_vector[i].x;
    vecy += local_vector[i].y;
  }
  gx = vecx/num_points;
  gy = vecy/num_points;
}
