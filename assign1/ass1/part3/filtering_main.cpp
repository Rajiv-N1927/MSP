/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include "filters.h"
#include <math.h>

#define PI 3.141592F
#define SCALE 2.5F
/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

//Generates a hanning windowed sinc
float *win_sinc(float shift, int extent) {
    int dim = 2*extent + 1;
    float *sinc = new float[dim]; // sin(Pi*x)/Pi*x
    float *hanning = new float[dim]; // 0.5(1+cos(Pi*x/t))
    float *sinc_filt = sinc + extent;
    for ( int i = -extent; i <= extent; i++ ) {
        sinc_filt[i] = sin(( /*(2.0F / 5.0F) */ PI * ( (float)i - (float)shift )) ) / (PI * (float)(i-shift));
        hanning[i+extent] = 0.5F - 0.5F*cos(( 2.0F * PI * ((float)(i + extent - shift)) / ((float)dim) ));
        sinc_filt[i] = sinc_filt[i]*hanning[i+extent];
    }
    //Normalization
    if ( (float)shift == 0.0F ) *sinc_filt = 1.0F;
    float sum = 0;
    for ( int i = -extent; i <= extent; i++ ) {
        sum+=sinc_filt[i];
    }
    for ( int i = -extent; i <= extent; i++ ) {
        sinc_filt[i]/=sum;
    }
    delete[] hanning;
    return sinc;
}
//The filter is a 1d filter
void sep_convolve(my_image_comp *in, my_image_comp *out, float ** filter,
    int extent, int no_filters) {
    //Only need the centre of the filter since its 1d;
    int r, c;
    float **filtpos = new float*[no_filters];
    for ( int i = 0; i < no_filters; i++ ) {
        filtpos[i] = filter[i] + extent;
    }
    //Need to follow the samples on the output
    /* First doing horizontal convolution from height - extent to
       height + extent. Take it from the input and use nearest sample
       Create an intermediate image
    */
    my_image_comp med_comp;
    float *input, *output, *filt;
    med_comp.init( in->height, ceil( (float)in->width * SCALE ), extent );
    for ( r = 0; r < med_comp.height; r++ ) {
        for ( c = 0; c < med_comp.width; c++ ) {
            output = med_comp.buf + r*med_comp.stride + c;
            int filt_num = c % no_filters;
            filt = filtpos[filt_num];
            if ( filt_num == 0 ) {
                input = in->buf + r*in->stride + 2*(c/5);
                *output = *input;
                continue;
            } else if ( filt_num == 1) {
                input = in->buf + r*in->stride + 2*((c-1)/5);
            } else if ( filt_num == 2) {
                input = in->buf + r*in->stride + 2*((c-2)/5)+1;
            } else if ( filt_num == 3) {
                input = in->buf + r*in->stride + 2*((c-3)/5)+1;
            } else if ( filt_num == 4) {
                input = in->buf + r*in->stride + 2*((c-4)/5)+2;
            }
            float sum = 0.0F;
            for ( int i = -extent; i <= extent; i++ ) {
                sum += input[i]*filt[i];
            }
            *output = sum;
        }
    }
    //Perform boundary extension on intermediate filter
    med_comp.perform_boundary_extension();
    //Perform convolution vertically
    for ( c = 0; c < out->width; c++ ) {
        for ( r = 0; r < out->height-1; r++ ) {
            output = out->buf + r*out->stride + c;
            int filt_num = r % no_filters;
            filt = filtpos[filt_num];
            if ( filt_num == 0 ) {
                input = med_comp.buf + med_comp.stride*(2*r/5) + c;
                *output = *input;
                continue;
            } else if ( filt_num == 1) {
                input = med_comp.buf + med_comp.stride*(2*(r-1)/5) + c;
            } else if ( filt_num == 2) {
                input = med_comp.buf + med_comp.stride*(2*(r-2)/5+1) + c;
            } else if ( filt_num == 3) {
                input = med_comp.buf + med_comp.stride*(2*(r-3)/5 + 1) + c;
            } else if ( filt_num == 4) {
                input = med_comp.buf + med_comp.stride*(2 * (r-4)/5 + 2) + c;
            }
            float sum = 0.0F;
            for ( int i = -extent; i <= extent; i++ ) {
                sum += input[i*med_comp.stride]*filt[i];
            }
            *output = sum;
        }
    }
    delete[] filtpos;
}



/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <extent>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];

      //Get the extent
      int in_extent = atoi(argv[3]);

      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width, in_extent);

      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init( ceil( (float)height * SCALE ),
            ceil( (float)width * SCALE ), 0 ); // Don't need a border for output

      // Process the image, all in floating point (easy)
      for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();

      //Get the sinc filters
      float *sinc_none = win_sinc(0.0F, in_extent);
      float *sinc_one = win_sinc(0.4F, in_extent);
      float *sinc_two = win_sinc(-0.2F, in_extent);
      float *sinc_three = win_sinc(0.2F, in_extent);
      float *sinc_four = win_sinc(-0.4F, in_extent);
      int num_filt = 5;
      float **filters = new float*[num_filt];
      filters[0] = sinc_none;
      filters[1] = sinc_one;
      filters[2] = sinc_two;
      filters[3] = sinc_three;
      filters[4] = sinc_four;

      //Time this
      for ( n=0; n < num_comps; n++ ) {
          sep_convolve(input_comps + n, output_comps + n, filters, in_extent, num_filt);
      }

      // Write the image back out again
      width = output_comps[0].width;
      height = output_comps[0].height;
      io_byte *out_line = new io_byte[width*num_comps];
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,num_comps)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = out_line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < width; c++, dst+=num_comps) {
                  *dst = (io_byte) (src[c]); // The cast to type "io_byte" is
                        // required here, since floats cannot generally be
                        // converted to bytes without loss of information.  The
                        // compiler will warn you of this if you remove the cast.
                        // There is in fact not the best way to do the
                        // conversion.  You should fix it up in the lab.
              }
            }
          bmp_out__put_line(&out,out_line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] out_line;
      delete[] input_comps;
      delete[] output_comps;
      delete[] sinc_none;
      delete[] sinc_one;
      delete[] sinc_two;
      delete[] sinc_three;
      delete[] sinc_four;
      delete[] filters;
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
