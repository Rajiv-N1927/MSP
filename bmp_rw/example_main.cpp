/*****************************************************************************/
// File: io_examples.cpp
// Author: David Taubman
// Last Revised: 19 August, 2002
/*****************************************************************************/
// Copyright 2001, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "io_custom.h"

/* ========================================================================= */
/*                             External Functions                            */
/* ========================================================================= */

/*****************************************************************************/
/* EXTERNAL                           main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 3)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      //Write the image back out again
      io_image my_img = io_image(argv[1], 0);

      //my_img.read(data);
      my_img.mod_comp(BLUE, 200);
      io_byte * new_data = my_img.write(), *dp;

      bmp_out out;
      int n;
      int height = my_img.comps[0].height;
      int width = my_img.comps[0].width;
      int planes = my_img.num_components;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,planes)) != 0)
        throw err_code;
      for (dp=new_data, n=height; n > 0; n--, dp+=width*planes)
        bmp_out__put_line(&out,dp);
      bmp_out__close(&out);

      delete[] new_data;
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
