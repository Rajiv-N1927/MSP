#include "io_custom.h"
#include <math.h>


//image composition functions
io_comp::io_comp() {}
io_comp::io_comp(int width, int height, int border) {
    this->width = width;
    this->height = height;
    this->border = border % 2 == 0 ? border : border + 1;
    this->handle = new int[(width+border)*(height+border)];
    this->buf = this->handle + (int)(border/2)*(width+border) + (int)border/2;
}

io_comp::~io_comp() {
    //delete [] this->handle;
}

void io_comp::modify(int toadd) {
    int *max = this->buf + (this->width+this->border)*(this->height-1);
    for (int *ptr = this->buf; ptr < max; ptr += this->width+this->border ) {
        for ( int *pos = ptr; pos < ptr + this->width; pos++ ) {
            *pos += toadd;
        }
    }
}

//image struct functions

io_image::io_image(int num_components, int width, int height, int border) {
    this->num_components = num_components;
    this->comps = new io_comp[num_components];
    for (int i = 0; i < num_components; i++) {
        this->comps[i] = io_comp(width, height, border);
    }
    this->image_size = this->comps[0].width*this->comps[0].height*this->num_components;
}

io_image::~io_image() {
    delete [] this->comps->handle;
    delete [] this->comps;
}

/*
    Traverse data and fill the components based on color
    Traverse through height, width then components
    Read image from the bottom left corner to the top left corner
*/

void io_image::read(io_byte * data) {
    int **cur_comp_pos = new int*[this->num_components]; // Tracks current component position
    for ( int i = 0; i < this->num_components; i++ ) {
        cur_comp_pos[i] = this->comps[i].buf;
    }
    io_byte *cur_data_pos, *row_pos; // Tracks current position in data
    int stride = this->comps[0].width*this->num_components;
    cur_data_pos = data+this->image_size - stride;
    for (; cur_data_pos != data; cur_data_pos -= stride) {
        for ( int i = 0; i < this->num_components; i++ ) {
            for ( row_pos=cur_data_pos; row_pos < cur_data_pos+stride;
                row_pos+=this->num_components) {
                *cur_comp_pos[i] = (int)(*(row_pos+i));
                cur_comp_pos[i]++;
            }
            //Add the border here and add w/e else you want to within the border
            cur_comp_pos[i] += this->comps[i].border;
        }
    }
    delete [] cur_comp_pos;
}

io_byte * io_image::write(void) {
    int **cur_comp_pos = new int*[this->num_components];
    io_byte *data = new io_byte[this->image_size];
    for ( int i = 0; i < this->num_components; i++ ) {
        cur_comp_pos[i] = this->comps[i].buf;
    }
    io_byte *cur_data_pos, *row_pos; // Tracks current position in data
    int stride_comp = (this->comps[0].width+this->comps[0].border)*this->num_components;
    int stride = this->comps[0].width*this->num_components;
    cur_data_pos = data+this->image_size-stride;
    for (; cur_data_pos != data; cur_data_pos -= stride) {
        for ( int i = 0; i < this->num_components; i++ ) {
            for ( row_pos=cur_data_pos; row_pos < cur_data_pos+stride;
                row_pos+=this->num_components) {
                *(row_pos+i) = (io_byte)(*cur_comp_pos[i]);
                cur_comp_pos[i]++;
            }
            //Add the border here and add w/e else you want to within the border
            cur_comp_pos[i] += this->comps[i].border;
        }
    }
    delete [] cur_comp_pos;
    return data;
}

void io_image::mod_comp(int color, int toadd) {
    if ( this->num_components != 1 ) {
        if ( color >= 0 && color <= 2 ) {
            this->comps[color].modify(toadd);
        }
    } else {
        printf("Greyscale image\n");
        this->comps[0].modify(toadd);
    }
}
