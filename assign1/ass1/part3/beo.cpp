/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include <math.h>

#define FILTER_EXTENT atoi(argv[3])
#define FILTER_LENGTH (2*FILTER_EXTENT+1)
/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
}


/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter_rows(my_image_comp *in, my_image_comp *out,
    float sinc_one[], float sinc_two[], float sinc_three[], float sinc_four[], float hanning[], int H)
{

	// Create the filter as a local array on the stack, which can accept
	// indices in the range -H to +H.
	float* filter_buf_one = new float[2 * H + 1];
	float* filter_buf_two = new float[2 * H + 1];
	float* filter_buf_three = new float[2 * H + 1];
	float* filter_buf_four = new float[2 * H + 1];

	float *mirror_psf_one = filter_buf_one + H;
	float *mirror_psf_two = filter_buf_two + H;
	float *mirror_psf_three = filter_buf_three + H;
	float *mirror_psf_four = filter_buf_four + H;
	// `mirror_psf' points to the central tap in the filter

	//Convolution performed on every row
	int r, c, i;
	for (r = 0; r < out->height; r++)
		for (c = 0; c < out->width; c++)
		{
			float *op = out->buf + r * out->stride + c;
			float sum = 0.0F;
			float normalization = 0.0F;
			if (c % 5 == 0)
			{
				float *ip = in->buf + r * in->stride + (2 * (c / 5));
				*op = *ip;
			}
			else if (c % 5 == 1)
			{
				float *ip = in->buf + r * in->stride + (2 * ((c - 1) / 5));
				for (i = -H; i <= H; i++)
				{
					mirror_psf_one[i] = hanning[i + H] * sinc_one[H + i];
					normalization = normalization + mirror_psf_one[i];
				}
				for (int x = -H; x <= H; x++)
					sum += ip[x] * mirror_psf_one[x];
				*op = sum / normalization;
			}
			else if (c % 5 == 2)
			{
				float *ip = in->buf + r * in->stride + (2 * ((c - 2) / 5) + 1);
				for (i = -H; i <= H; i++)
				{
					mirror_psf_two[i] = hanning[i + H] * sinc_two[H + i];
					normalization = normalization + mirror_psf_two[i];
				}
				for (int x = -H; x <= H; x++)
					sum += ip[x] * mirror_psf_two[x];
				*op = sum / normalization;
			}
			else if (c % 5 == 3)
			{
				float *ip = in->buf + r * in->stride + (2 * ((c - 3) / 5) + 1);
				for (i = -H; i <= H; i++)
				{
					mirror_psf_three[i] = hanning[i + H] * sinc_three[H + i];
					normalization = normalization + mirror_psf_three[i];
				}
				for (int x = -H; x <= H; x++)
					sum += ip[x] * mirror_psf_three[x];
				*op = sum / normalization;
			}
			else
			{
				float *ip = in->buf + r * in->stride + (2 * ((c - 4) / 5) + 2);
				for (i = -H; i <= H; i++)
				{
					mirror_psf_four[i] = hanning[i + H] * sinc_four[H + i];
					normalization = normalization + mirror_psf_four[i];
				}
				for (int x = -H; x <= H; x++)
					sum += ip[x] * mirror_psf_four[x];
				*op = sum / normalization;
			}
		}

	delete[] filter_buf_one;
	delete[] filter_buf_two;
	delete[] filter_buf_three;
	delete[] filter_buf_four;
}

void apply_filter_columns(my_image_comp *in, my_image_comp *out,
    float sinc_one[], float sinc_two[], float sinc_three[], float sinc_four[], float hanning[], int H)
{

	// Create the filter as a local array on the stack, which can accept
	// indices in the range -H to +H.
	float* filter_buf_one = new float[2 * H + 1];
	float* filter_buf_two = new float[2 * H + 1];
	float* filter_buf_three = new float[2 * H + 1];
	float* filter_buf_four = new float[2 * H + 1];

	float *mirror_psf_one = filter_buf_one + H;
	float *mirror_psf_two = filter_buf_two + H;
	float *mirror_psf_three = filter_buf_three + H;
	float *mirror_psf_four = filter_buf_four + H;
	// `mirror_psf' points to the central tap in the filter

	//Convolution performed on every column
	int r, c, i;
	for (c = 0; c < out->width; c++)
		for (r = 0; r < out->height; r++)
		{
			float *op = out->buf + r * out->stride + c;
			float sum = 0.0F;
			float normalization = 0.0F;
			if (r % 5 == 0)
			{
				float *ip = in->buf + (2 * (r / 5)) * in->stride + c;
				*op = *ip;
			}
			else if (r % 5 == 1)
			{
				float *ip = in->buf + (2 * ((r - 1) / 5)) * in->stride + c;
				for (i = -H; i <= H; i++)
				{
					mirror_psf_one[i] = hanning[i + H] * sinc_one[H + i];
					normalization = normalization + mirror_psf_one[i];
				}
				for (int y = -H; y <= H; y++)
					sum += ip[y * in->stride] * mirror_psf_one[y];
				*op = sum / normalization;
			}
			else if (r % 5 == 2)
			{
				float *ip = in->buf + (2 * ((r - 2) / 5) + 1) * in->stride + c;
				for (i = -H; i <= H; i++)
				{
					mirror_psf_two[i] = hanning[i + H] * sinc_two[H + i];
					normalization = normalization + mirror_psf_two[i];
				}
				for (int y = -H; y <= H; y++)
					sum += ip[y * in->stride] * mirror_psf_two[y];
				*op = sum / normalization;
			}
			else if (r % 5 == 3)
			{
				float *ip = in->buf + (2 * ((r - 3) / 5) + 1) * in->stride + c;
				for (i = -H; i <= H; i++)
				{
					mirror_psf_three[i] = hanning[i + H] * sinc_three[H + i];
					normalization = normalization + mirror_psf_three[i];
				}
				for (int y = -H; y <= H; y++)
					sum += ip[y * in->stride] * mirror_psf_three[y];
				*op = sum / normalization;
			}
			else
			{
				float *ip = in->buf + (2 * ((r - 4) / 5) + 2) * in->stride + c;
				for (i = -H; i <= H; i++)
				{
					mirror_psf_four[i] = hanning[i + H] * sinc_four[H + i];
					normalization = normalization + mirror_psf_four[i];
				}
				for (int y = -H; y <= H; y++)
					sum += ip[y * in->stride] * mirror_psf_four[y];
				*op = sum / normalization;
			}
		}

	delete[] filter_buf_one;
	delete[] filter_buf_two;
	delete[] filter_buf_three;
	delete[] filter_buf_four;
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{

	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> H\n", argv[0]);
		return -1;
	}

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;
		my_image_comp *input_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, FILTER_EXTENT); // Leave a border

		int r; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps)
					dst[c] = (float)*src; // The cast to type "float" is not
										  // strictly required here, since bytes can always be
										  // converted to floats without any loss of information.
			}
		}
		bmp_in__close(&in);

		// Create filters
#define PI 3.141592653589793F

		int H = FILTER_EXTENT;

		float *sinc_one = new float[2 * H + 1];
		float *sinc_two = new float[2 * H + 1];
		float *sinc_three = new float[2 * H + 1];
		float *sinc_four = new float[2 * H + 1];
		float i = 0.0F;

		while (i <= 2.0F * (float)H)
		{
			sinc_one[(int)i] = (sin(PI * (float)(i - (float)H - 0.4F))) / (PI * (float)(i - (float)H - 0.4F));
			sinc_two[(int)i] = (sin(PI * (float)(i - (float)H + 0.2F))) / (PI * (float)(i - (float)H + 0.2F));
			sinc_three[(int)i] = (sin(PI * (float)(i - (float)H - 0.2F))) / (PI * (float)(i - (float)H - 0.2F));
			sinc_four[(int)i] = (sin(PI * (float)(i - (float)H + 0.4F))) / (PI * (float)(i - (float)H + 0.4F));
			i = i + 1;
		}

		float m = 0.0F;
		float *hanning = new float[2 * H + 1];

		while (m < (float)(2 * H + 1))
		{
			hanning[(int)m] = 0.5F - 0.5F * cos((float)(2.0F * PI * m) / (float)(2 * H + 1));
			if (H == 0) {
				hanning[0] = 1.0F;
			}
			m = m + 1.0F;
		}

		// Allocate storage for the first filtered output
		my_image_comp *output_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			output_comps[n].init(height, floor((float)(width - 1) * 2.5F) + 1, FILTER_EXTENT);

		// Allocate storage for the second filtered output
		my_image_comp *output_comps_2 = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			output_comps_2[n].init(floor((float)(height - 1) * 2.5F) + 1, floor((float)(width - 1) * 2.5F) + 1, 0);

		// Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++)
			input_comps[n].perform_boundary_extension();
		for (n = 0; n < num_comps; n++)
			apply_filter_rows(input_comps + n, output_comps + n, sinc_one, sinc_two, sinc_three, sinc_four, hanning, H);
		for (n = 0; n < num_comps; n++)
			output_comps[n].perform_boundary_extension();
		for (n = 0; n < num_comps; n++)
			apply_filter_columns(output_comps + n, output_comps_2 + n, sinc_one, sinc_two, sinc_three, sinc_four, hanning, H);

		// Write the image back out again
		height = output_comps[0].height;
		width = output_comps[0].width;
        io_byte *line_new = new io_byte[width*num_comps];
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], width, height, num_comps)) != 0)
			throw err_code;
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = line_new + n; // Points to first sample of component n
				float *src = output_comps_2[n].buf + r * output_comps_2[n].stride;
				for (int c = 0; c < width; c++, dst += num_comps)
					*dst = (io_byte)src[c]; // The cast to type "io_byte" is
											// required here, since floats cannot generally be
											// converted to bytes without loss of information.  The
											// compiler will warn you of this if you remove the cast.
											// There is in fact not the best way to do the
											// conversion.  You should fix it up in the lab.
			}
			bmp_out__put_line(&out, line_new);
		}
		bmp_out__close(&out);
		delete[] line;
        delete[] line_new;
		delete[] input_comps;
		delete[] output_comps;
		delete[] output_comps_2;
		delete[] sinc_one;
		delete[] sinc_two;
		delete[] sinc_three;
		delete[] sinc_four;
		delete[] hanning;
	}
	catch (int exc) {
		if (exc == IO_ERR_NO_FILE)
			fprintf(stderr, "Cannot open supplied input or output file.\n");
		else if (exc == IO_ERR_FILE_HEADER)
			fprintf(stderr, "Error encountered while parsing BMP file header.\n");
		else if (exc == IO_ERR_UNSUPPORTED)
			fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
				"simple example supports only 8-bit and 24-bit data.\n");
		else if (exc == IO_ERR_FILE_TRUNC)
			fprintf(stderr, "Input or output file truncated unexpectedly.\n");
		else if (exc == IO_ERR_FILE_NOT_OPEN)
			fprintf(stderr, "Trying to access a file which is not open!(?)\n");
		return -1;
	}
	return 0;
}
