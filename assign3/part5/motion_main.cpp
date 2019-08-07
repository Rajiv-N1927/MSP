/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"

#define TRACE(a,b,c,d)     a + d
#define DET(a,b,c,d)       a*d-b*c
#define MAX_V              2<<16


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
/*
 * Generates the tensor struct
 * Each 'point' is comprised of 4 values
*/
int * tensor_struct(mvector * dogs, int num_dogs) {
  int * tensor = new int[4 * num_dogs];
  for ( int pos = 0, cpos = 0; pos < num_dogs; pos++, cpos += 4 ) {
    tensor[cpos]   = dogs[pos].x * dogs[pos].x;
    tensor[cpos+1] = dogs[pos].x * dogs[pos].y;
    tensor[cpos+2] = dogs[pos].x * dogs[pos].y;
    tensor[cpos+3] = dogs[pos].y * dogs[pos].y;
  }
  return tensor;
}

/*
 * At this point decide where to start the struct
 * Aggregation of tensor struct over a block
*/
int * int_tensor(int * tensor_struct, int height,
    int width, int extent, int &pts) {
  int dim = 2*extent+1;
  int num_pts = floor((height-(2*dim))*(width-(2*dim))*4);
  int * tensor_int = new int[num_pts];
  int row, col, block_row, block_col, pos = 0;
  for ( row = dim; row < height - dim; row++ ) {
    for ( col = dim; col < width - dim; col++, pos+=4 ) {
      int * ts_pos = tensor_struct + 4*(row*width + col);
      for ( block_row = 0; block_row < dim; block_row++ ) {
        int * ts_int_pos = ts_pos + 4*(block_row*width);
        for ( block_col = 0; block_col < dim; block_col++ ) {
          int * curpos = ts_int_pos + 4*block_col;
          tensor_int[pos] += curpos[0];
          tensor_int[pos+1] += curpos[1];
          tensor_int[pos+2] += curpos[2];
          tensor_int[pos+3] += curpos[3];
        }
      }
    }
  }
  pts = pos;
  return tensor_int;
}

/*
 * Generate the figure of merits being careful of division by 0
*/

float * FOM_calc(int * block_tensor, int block_pts) {
  int fom_pts = block_pts/4;
  float * fom = new float[fom_pts];
  int bt_pos = 0, fom_pos = 0;
  for ( int i = 0; i < fom_pts; i++, bt_pos+=4 ) {
    float a = (float)block_tensor[bt_pos];
    float b = (float)block_tensor[bt_pos+1];
    float c = (float)block_tensor[bt_pos+2];
    float d = (float)block_tensor[bt_pos+3];
    float num = DET(a, b, c, d);
    float den = TRACE(a, b, c, d);
    fom[i] = den != 0 ? num/den : MAX_V;
  }
  return fom;
}

/*
 * Generate the keypoints given a threshold tau, and number of keypoints
 * The vectors are all shifted by 2*extent + 1 in the x and y direction
 * This is because the FOM are within the bounds dim -> height - dim &
 * dim -> width - dim.
 * the low and high start off as -1
*/

mvector * keyPoints(float &thresh, float &thresh_low, float &thresh_high, int num_kp,
  float * fom, int sa_extent, int offset, int height, int width) {
  //Max number of keypoints possible
  mvector * kp = new mvector[(height-1)*(width-1)];
  int kp_pos = 0;
  int row, col;
  for ( row = sa_extent + 1; row < height-1-sa_extent; row++ ) {
    float * fom_pos = fom + row*width;
    for ( col = sa_extent + 1; col < width-1-sa_extent; col++ ) {
      float top = fom_pos[(row-1)*width];
      float left = fom_pos[col-1];
      float right = fom_pos[col+1];
      float bot = fom_pos[(row+1)*width];
      float cur = fom_pos[col];
      if ( cur >= top && cur >= left && cur >= right
        && cur >= bot && cur >= thresh) {
          kp[kp_pos].x = col;
          kp[kp_pos].y = row;
          kp_pos++;
      }
    }
  }
  if ( thresh_low == -1 || thresh_high == -1 ) {
    if ( kp_pos > num_kp ) {
      thresh_low = thresh;
      thresh = 2.0F * thresh;
    } else if ( kp_pos < num_kp ) {
      thresh_high = thresh;
      thresh = 0.5F * thresh;
    } else { //Success
      return kp;
    }
    //Otherwise stackoverflow!
    delete[] kp;
    return keyPoints(thresh, thresh_low,
      thresh_high, num_kp, fom, sa_extent, offset, height, width);
  } else { //Both thresholds have been set
    if ( kp_pos == num_kp ) { //Success
     printf("Found %d keypoints!\n", kp_pos);
     return kp;
    } else if ( kp_pos > num_kp ) { //Threshold Too low
      thresh_low = thresh;
    } else if ( kp_pos < num_kp ) { //Threshold Too high
      thresh_high = thresh;
    }
    thresh = 0.5F * (thresh_low + thresh_high);
    delete[] kp;
    return keyPoints(thresh, thresh_low,
      thresh_high, num_kp, fom, sa_extent, offset, height, width);
  }
  return kp;
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
  if (argc != 8)
    {
      fprintf(stderr,
              "Usage: %s <bmp ref> <bmp target> <bmp out> <block extent> "
              "<search area extent> <sigma> <number of keypoints>\n",
              argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      //Inputs
      int extent       = atoi(argv[4]);
      int dim          = 2*extent+1;
      int sa_extent    = atoi(argv[5]);
      int sigma        = atoi(argv[6]);
      int noKp         = atoi(argv[7]);

      // Read the input image
      bmp_in in[2];
      if ((err_code = bmp_in__open(&in[0],argv[1])) != 0)
        throw err_code;
      if ((err_code = bmp_in__open(&in[1],argv[2])) != 0)
        throw err_code;

      int width = in[0].cols, height = in[0].rows;
      if ((width != in[1].cols) || (height != in[1].rows))
        {
          fprintf(stderr,"The two input frames have different dimensions.\n");
          return -1;
        }

      my_image_comp mono[2];
      //Maximum extent is sa_extent, waste of memory but easiest way
      mono[0].init(height,width, sa_extent + extent);
      mono[1].init(height,width, sa_extent + extent);

      int n, r, c;
      int num_comps = in[0].num_components;
      io_byte *line = new io_byte[width*num_comps];
      for (n=0; n < 2; n++)
        {
          for (r=height-1; r >= 0; r--)
            { // "r" holds the true row index we are reading, since the image
              // is stored upside down in the BMP file.
              if ((err_code = bmp_in__get_line(&(in[n]),line)) != 0)
                throw err_code;
              io_byte *src = line; // Points to first sample of component n
              int *dst = mono[n].buf + r * mono[n].stride;
              for (c=0; c < width; c++, src+=num_comps)
                dst[c] = *src;
            }
          bmp_in__close(&(in[n]));
        }

      my_image_comp output;
      output.init(height, width, extent);
      //Filter the input image here
      float * gauss_filt = gaussian_filter(sigma, extent);
      //Perform extension
      mono[1].perform_symmetric_extension();
      //Perform low pass filter with gaussian
      sep_convolve(&(mono[1]), &output, gauss_filt, extent);
      //Get the gradient vector on the target frame
      output.perform_symmetric_extension();
      int length_DOG = 0, len_tensor = 0; //Should be img width * height
      mvector * dog_vecs = DOG_vector(&output, length_DOG);
      int * tensor   = tensor_struct(dog_vecs, length_DOG);
      int * block_tensor = int_tensor(tensor, output.height,
        output.width, extent, len_tensor);
      float * fom_pts = FOM_calc(block_tensor, len_tensor);
      //Find the keypoints
      float thresh = 10.0F;
      float thresh_low = -1.0F;
      float thresh_high = -1.0F;
      int n_height = output.height-(2*dim);
      int n_width  = output.width-(2*dim);
      mvector * kp = keyPoints(thresh, thresh_low, thresh_high, noKp,
        fom_pts, sa_extent, dim, n_height, n_width);

      //Setup the kp struct
      kp_generator kp_gen;
      kp_gen.tgt = &(mono[1]);
      kp_gen.ref = &(mono[0]);
      kp_gen.kp_offsets = kp;
      kp_gen.local_vector = new mvector[noKp];
      kp_gen.b_extent = extent;
      kp_gen.sa_extent = sa_extent;
      //Iterate over the keypoints


      // Write the image back out again
      width = output.width;
      height = output.height;
      io_byte *out_line = new io_byte[width*num_comps];
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[3],width,height,1)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
            io_byte *dst = out_line+n; // Points to first sample of component n
            int *src = output.buf + r * output.stride;
            for (int c=0; c < width; c++, dst+=num_comps) {
                *dst = (io_byte) (src[c]);
            }
          bmp_out__put_line(&out, out_line);
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
