/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"

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


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 7)
    {
      fprintf(stderr,
              "Usage: %s <bmp frame 1> <bmp frame 2> <bmp out> <Block extent>"
              " <Search area extent> <keypoint separation distance>\n",
              argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      //Inputs
      int block_extent = atoi(argv[4]);
      int sa_extent    = atoi(argv[5]);
      int kp_sep       = atoi(argv[6]);

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
      mono[0].init(height,width, sa_extent);
      mono[1].init(height,width, sa_extent);

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

      // Allocate storage for the motion compensated output
      my_image_comp output;
      output.init(height,width,0); // Don't need a border for output

      //Keypoint generator
      kp_generator kp;
      //Find the key points for the target image with specified inputs
      //First image is the target the next is the reference
      kp.init(&(mono[1]), &(mono[0]), block_extent, sa_extent, kp_sep);
      //Generate the vectors for image based on target image
      kp.generate_kp();
      kp.generate_shifted_img(&output);
      // Now perform simple motion estimation and compensation
      int nominal_block_width = 32;
      int nominal_block_height = 32;
      int block_width, block_height;
      /*
      * Modify to keypoints
      * K = H + S + delta*(N)
      */
      // for (r=0; r < height; r+=block_height)
      //   {
      //     block_height = nominal_block_height;
      //     if ((r+block_height) > height)
      //       block_height = height-r;
      //     for (c=0; c < width; c+=block_width)
      //       {
      //         block_width = nominal_block_width;
      //         if ((c+block_width) > width)
      //           block_width = width-c;
      //         mvector vec = find_motion(&(mono[0]),&(mono[1]), sa_extent,
      //                                   r,c,block_width,block_height);
      //         motion_comp(&(mono[0]),&output,vec,
      //                     r,c,block_width,block_height);
      //       }
      //   }

      //Write the motion compensated image out
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[3],width,height,1)) != 0)
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
