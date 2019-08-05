#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*****************************************************************************/
// File: image_comps.h
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

/*****************************************************************************/
/* STRUCT                        my_image_comp                               */
/*****************************************************************************/
struct mvector {
    mvector() { x=y=0; }
    int x, y;
};

struct filter_manager {
    int extent;
    float *taps;
    filter_manager();
    ~filter_manager();
    void init(float *, int);
    //Returns a normalized filter
    void normalize_filter();
    void mirror_filter();
};

struct my_image_comp {
    // Data members: (these occupy space in the structure's block of memory)
    int width;
    int height;
    int stride;
    int border; // Extra rows/cols to leave around the boundary
    int *handle; // Points to start of allocated memory buffer
    int *buf; // Points to the first real image sample
    // Function members: (these do not occupy any space in memory)
    my_image_comp();
    ~my_image_comp();
    //height, width, border
    void init(int, int, int);
    void perform_boundary_extension();
       // This function is implemented in "io_comps.cpp".
    void perform_symmetric_extension();
    void perform_point_symmetric_extension();
  };

  /*
  * Keypoint generator generates the keypoints for the target image
  * This will be used in conjunction with an SAD vector generator with
  * the reference image to find the best match
  */

  struct kp_generator {
    my_image_comp * tgt;
    mvector * kp_offsets, * local_vector;
    float gx, gy;
    int b_extent, sa_extent, delta;
    //The bounds that can be travelled for the keypoints
    int max_height, max_width, base_offset, num_points;
    kp_generator();
    ~kp_generator();
    //img, block extent, search area extent, delta
    void init(my_image_comp *, int, int, int);
    void generate_kp(my_image_comp *);
    void sad_vector(my_image_comp *, int, int, int);
    void produce_gv();
  };
