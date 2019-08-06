/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#define PI 3.141592F

/*****************************************************************************/
/* STATIC                         find_motion                                */
/*****************************************************************************/

static mvector
  find_motion(my_image_comp *ref, my_image_comp *tgt, int SA_extent,
              int start_row, int start_col, int block_width, int block_height)
  /* This function finds the motion vector which best describes the motion
     between the `ref' and `tgt' frames, over a specified block in the
     `tgt' frame.  Specifically, the block in the `tgt' frame commences
     at the coordinates given by `start_row' and `start_col' and extends
     over `block_width' columns and `block_height' rows.  The function finds
     the translational offset (the returned vector) which describes the
     best matching block of the same size in the `ref' frame, where
     the "best match" is interpreted as the one which minimizes the sum of
     absolute differences (SAD) metric. */
{
  mvector vec, best_vec;
  int sad, best_sad=256*block_width*block_height; //Max size
  //The size of the vec.y represents the area in which the search occurs in the y direction
  for (vec.y=-SA_extent; vec.y <= SA_extent; vec.y++)   {
  //The size of the vec.x represents the area in which the search occurs in the x direction
    for (vec.x=-SA_extent; vec.x <= SA_extent; vec.x++)
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

  return best_vec;
}

/*****************************************************************************/
/* STATIC                         motion_comp                                */
/*****************************************************************************/

static void
  motion_comp(my_image_comp *ref, my_image_comp *tgt, mvector vec,
              int start_row, int start_col, int block_width, int block_height)
  /* This function transfers data from the `ref' frame to a block within the
     `tgt' frame, thereby realizing motion compensation.  The motion in
     question has already been found by `find_motion' and is captured by
     the `vec' argument.  The block in the `tgt' frame commences
     at the coordinates given by `start_row' and `start_col' and extends
     for `block_width' columns and `block_height' rows. */
{
  int r, c;
  int ref_row = start_row - vec.y;
  int ref_col = start_col - vec.x;
  int *rp = ref->buf + ref_row*ref->stride + ref_col;
  int *tp = tgt->buf + start_row*tgt->stride + start_col;
  for (r=block_height; r > 0; r--,
       rp+=ref->stride, tp+=tgt->stride)
    for (c=0; c < block_width; c++)
      tp[c] = rp[c];
}

float * gaussian_filter(float var, int extent) {
  if ( var == 0.0F ) var = 1.0F;
  int dim = 2*extent + 1;
  float * gauss = new float[dim];
  float * cpos = gauss + extent;
  float var_sq = (float)(var*var);
  float cons = 1.0F/sqrt((2.0F*M_PI*var_sq));
  for ( int i = -extent; i <= extent; i++ ) {
    float num = (float)(i*i);
    cpos[i] = cons*((float)exp( -num/(2.0F*var_sq)) );
  }
  //Normalize
  float sum = 0.0F;
  for ( int i = -extent; i <= extent; i++ ) {
    sum += cpos[i];
  }
  for ( int i = -extent; i <= extent; i++ ) {
    cpos[i]/=sum;
  }
  return gauss;
}

float * filter_2d(float * gauss_filt, int extent) {
  int dim = 2*extent + 1;
  int num_taps = dim * dim;
  float * filt_2d = new float[num_taps];
  float * filt2dpos = filt_2d + extent + extent*dim;
  float * filtpos = gauss_filt + extent;
  for ( int i = -extent; i <= extent; i++ ) {
    float * c2dpos = filt2dpos + i*dim;
    float * cpos = filtpos + i;
    for ( int j = -extent; j <= extent; j++ ) {
      c2dpos[j] = (*cpos)*(filtpos[j]);
    }
  }
  return filt_2d;
}
//Sep convolve
void sep_convolve(my_image_comp * in, my_image_comp * out,
  float * filt, int extent) {
  my_image_comp med_comp;
  int r, c;
  float * filtpos = filt + extent;
  int * output, * input;
  med_comp.init( in->height, in->width, extent );
  for ( r = 0; r < med_comp.height; r++ ) {
      for ( c = 0; c < med_comp.width; c++ ) {
          output = med_comp.buf + r*med_comp.stride + c;
          input = in->buf + r*in->stride + c;
          float sum = 0.0F;
          for ( int i = -extent; i <= extent; i++ ) {
              sum += (float)(input[i])*filtpos[i];
          }
          *output = sum;
      }
  }
  //Perform boundary extension on intermediate filter
  med_comp.perform_symmetric_extension();
  //Perform convolution vertically
  for ( c = 0; c < out->width; c++ ) {
      for ( r = 0; r < out->height; r++ ) {
          output = out->buf + r*out->stride + c;
          input = med_comp.buf + med_comp.stride*r + c;
          float sum = 0.0F;
          for ( int i = -extent; i <= extent; i++ ) {
              sum += (float)(input[i*med_comp.stride])*filtpos[i];
          }
          *output = sum;
      }
  }

}

//Convolve
void convolve(my_image_comp * in, my_image_comp * out,
  float * filt, int extent) {
  int dim = 2*extent + 1;
  float * filtpos = filt + extent;
  //Filter horizontallly
  int r, c;
  for (r=0; r < out->height; r++)
    for (c=0; c < out->width; c++)
      {
        float *ip = (float *)(in->buf + r*in->stride + c);
        float *op = (float *)(out->buf + r*out->stride + c);
        float sum = 0.0F;
        for (int y=-extent; y <= extent; y++)
          for (int x=-extent; x <= extent; x++)
            sum += ip[y*in->stride+x] * filtpos[y*dim+x];
        *op = sum;
      }
}

//Get the gradient vector
mvector * DOG_vector(my_image_comp * tgt, int &length) {
  int npts = tgt->width * tgt->height;
  mvector * grad_vec = new mvector[npts];
  int row, col, pos = 0;
  for ( row = 0; row < tgt->height; row++ ) {
    int * tgt_pos = tgt->buf + row*tgt->stride;
    for ( col = 0; col < tgt->width; col++ ) {
      grad_vec[pos].x = (tgt_pos[col+1] - tgt_pos[col-1])/2;
      grad_vec[pos].y = (tgt_pos[col+tgt->stride] - tgt_pos[col-tgt->stride])/2;
      pos++;
    }
  }
  length = pos;
  return grad_vec;
}

/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 6)
    {
      fprintf(stderr,
              "Usage: %s <bmp img> <bmp out> <extent>"
              "<sigma> <numkp>\n",
              argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Inputs
      int extent = atoi(argv[3]);
      int sigma  = atoi(argv[4]);
      int noKp   = atoi(argv[5]);
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      my_image_comp mono;
      int height = in.rows;
      int width  = in.cols;
      mono.init(height,width,extent);
      int n, r, c;
      int num_comps = in.num_components;
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image
          // is stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&(in),line)) != 0)
            throw err_code;
          io_byte *src = line; // Points to first sample of component n
          int *dst = mono.buf + r * mono.stride;
          for (c=0; c < width; c++, src+=num_comps)
            dst[c] = *src;
        }
      bmp_in__close(&(in));

      // Allocate storage for the motion compensated output
      my_image_comp output;
      output.init(height,width,extent); // Don't need a border for output
      //Filter the input image here

      float * gauss_filt = gaussian_filter(sigma, extent);
      //Perform extension
      mono.perform_symmetric_extension();
      //Perform low pass filter with gaussian
      sep_convolve(&mono, &output, gauss_filt, extent);
      //Get the gradient vector on the target frame
      output.perform_symmetric_extension();
      int length_DOG = 0; //Should be img width * height
      mvector * dog_vecs = DOG_vector(&output, length_DOG);
      printf("%d\n", length_DOG);


      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,1)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          io_byte *dst = line; // Points to first sample of component n
          int *src = output.buf + r * output.stride;
          for (int c=0; c < width; c++, dst++)
            *dst = (io_byte) src[c];
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] gauss_filt;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  return 0;
}
