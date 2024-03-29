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
#include <unistd.h>

#define ARGS_START 5


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

void get_intensity(my_image_comp *in, my_image_comp *out, int threshold) {
    int r, c;
    // Perform the convolution
    for (r=0; r < out->height; r++)
        for (c=0; c < out->width; c++)
          {
            float *ip = in->buf + r*in->stride + c;
            float *op = out->buf + r*out->stride + c;
            if ( *ip < threshold ) *op = 0;
            else *op = 0xFF;
          }
}

void erode_image(my_image_comp *in, my_image_comp *out,
    struct_set * in_set, int num_coord ) {

    int r, c;
    // Perform the convolution
    for (r=0; r < out->height; r++)
        for (c=0; c < out->width; c++) {
            float *ip = in->buf + r*in->stride + c;
            float *op = out->buf + r*out->stride + c;
            *op = 255;
            for ( int pos = 0; pos < num_coord; pos++ ) {
                float val = *( ip + in_set[pos].x + in->stride*in_set[pos].y );
                if ( val != 0xFF ) {
                    *op = 0;
                    break;
                }
            }
        }
}

void open_image(my_image_comp *in, my_image_comp *out,
    struct_set * in_set, int num_coord ) {
    int r, c;
    // Perform the convolution
    for (r=0; r < out->height; r++)
        for (c=0; c < out->width; c++) {
            float *ip = in->buf + r*in->stride + c;
            float *op = out->buf + r*out->stride + c;
            int set = 1;
            for ( int pos = 0; pos < num_coord; pos++ ) {
                int struct_pos = in_set[pos].x + in->stride*in_set[pos].y;
                float val = *( ip + struct_pos );
                if ( (int)val == 0x0 ) {
                    set = 0;
                    break;
                }
            }
            if ( set ) {
                for ( int pos = 0; pos < num_coord; pos++ ) {
                    int struct_pos = in_set[pos].x + in->stride*in_set[pos].y;
                    *(op + struct_pos) = 0xFF;
                }
            }
        }
}

void close_image(my_image_comp *in, my_image_comp *out,
    struct_set * in_set, int num_coord ) {
    int r, c;
    // Perform the convolution
    for (r=0; r < out->height; r++)
        for (c=0; c < out->width; c++) {
            float *ip = in->buf + r*in->stride + c;
            float *op = out->buf + r*out->stride + c;
            int set = 0;
            for ( int pos = 0; pos < num_coord; pos++ ) {
                int struct_pos = in_set[pos].x + in->stride*in_set[pos].y;
                float val = *( ip + struct_pos );
                if ( (int)val == 0xFF ) {
                    set = 1;
                    break;
                }
            }
            if ( set ) {
                for ( int pos = 0; pos < num_coord; pos++ ) {
                    int struct_pos = in_set[pos].x + in->stride*in_set[pos].y;
                    *(op + struct_pos) = 0xFF;
                }
            }
        }
}

int getExtent(struct_set * coords, int num_coords) {
    int maxVal = 0;
    for ( int i = 0; i < num_coords; i++ ) {
        if ( abs(coords[i].x) > maxVal ) maxVal = abs(coords[i].x);
        if ( abs(coords[i].y) > maxVal ) maxVal = abs(coords[i].y);
    }
    return maxVal;
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[])
{
  if (argc < ARGS_START)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <threshold>\n",argv[0]);
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

      //For morphology!!
      int thresh = atoi(argv[3]);
      if ( thresh < 0 || thresh > 255 ) {
          char * msg = new char[100];
          msg = "Threshold out of range please input between 0 and 255\n";
          write(2, msg, 100);
          exit(0);
      }
      //Setting up circular structure
      int radius = atoi(argv[4]);
      struct circle_set myset;
      myset.init(radius);

      struct_set * s_set = new struct_set[argc-ARGS_START];
      int num_coords = argc-ARGS_START;
      for ( int i = ARGS_START; i < argc; i++ ) {
          s_set[i-ARGS_START].setCoord(argv[i]);
          // printf("xval: %d, yval: %d\n",
          //   s_set[i-ARGS_START].x, s_set[i-ARGS_START].y);
      }

      //Get the extent
      int in_extent = getExtent(s_set, num_coords);
      if ( radius > in_extent ) in_extent = radius;

      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width, 0); // h1 extent is 2

       /* #################################################
        * --------------SETTING UP THE INPUTS--------------
        * #################################################
        */
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
      my_image_comp *int_output_comps = new my_image_comp[num_comps];
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++) {
        int_output_comps[n].init(height,width, in_extent);
        output_comps[n].init(height,width, in_extent);
    }

      // Process the image, all in floating point (easy)
      for (n=0; n < num_comps; n++) {
        get_intensity(input_comps+n,int_output_comps+n, thresh);
        int_output_comps[n].perform_symmetric_extension();
        open_image(int_output_comps+n,output_comps+n, myset.c_set, myset.no_components);
      }

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,num_comps)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < width; c++, dst+=num_comps)
                *dst = (io_byte) (src[c]+0.5); // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
      delete[] int_output_comps;
      delete[] s_set;
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
