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


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

//Convolves based on the output perspective
void convolve(my_image_comp *in, my_image_comp *out, float * filter, int extent)
{
    int dim = 2*extent + 1;
    int taps = dim*dim;
    // Create the filter kernel as a local array on the stack, which can accept
    // row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
    float *mirror_psf = filter+(dim*extent)+extent;
          // `mirror_psf' points to the central tap in the filter

    // Check for consistent dimensions
    assert(in->border >= extent);
    assert((out->height <= in->height) && (out->width <= in->width));
    int r, c;
    // Perform the convolution
    for (r=0; r < out->height; r++)
    for (c=0; c < out->width; c++)
      {
        float *ip = in->buf + r*in->stride + c;
        float *op = out->buf + r*out->stride + c;
        float sum = 0.0F;
        for (int y=-extent; y <= extent; y++)
          for (int x=-extent; x <= extent; x++)
            sum += ip[y*in->stride+x] * mirror_psf[y*dim+x];
        *op = sum;
      }
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[])
{
  if (argc != 3)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
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
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width, 2); // No border required for bilin

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

      //Define the scaling factor
      float sf = 2.5;
      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++) {
          input_comps[n].perform_boundary_extension();
          output_comps[n].init( (int)((float)height*sf),
            (int)((float)width*sf), 0 ); // Don't need a border for output
        //Perform a billinear interpolation
          output_comps[n].perform_billinear_interpolation(input_comps + n, sf);
      }


      // Write the image back out again

      width = output_comps[0].width;
      height = output_comps[0].height;
      io_byte *out_line = new io_byte[width*num_comps];
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2], width, height,num_comps)) != 0) {
          throw err_code;
      }
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = out_line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < width; c++, dst+=num_comps)
                *dst = (io_byte) (src[c]+0.5); // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,out_line);
        }
      bmp_out__close(&out);
      delete[] out_line;
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
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
