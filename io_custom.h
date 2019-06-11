// Adding files from lab 1 Part 5
/*
    Creating a class for handling image components
    Each component corresponds to one of Red, Green or Blue (for 24 bit planes_bits)
    Otherwise it is just greyscale ( 8 bit )
*/
#include "io_bmp.h"

#define BLUE 0
#define GREEN 1
#define RED 2

struct io_comp {
    int width, height, border; // Border height and width are the same
    //Handlers will be integers when written to by typecasted to iobyte when written from
    int *buf, *handle;
    io_comp(int, int, int);
    ~io_comp();
    //Other functions to perform alterations on
    void modify(int); //This merely adds a number
};

struct io_image {
    int num_components;
    int image_size;
    io_comp *comps;
    io_image(int, int, int, int); //Num components, height, width, border
    ~io_image();
    void read(io_byte *); // Take in the image data to be broken into components
    void mod_comp(int, int);
    io_byte *write(void);
};

struct io_image_seq {
    io_image * images;
};
