/* uvid_decompress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 

     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).

   B. Bird - 07/15/2020
*/

#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include "input_stream.hpp"
#include "uvid_common.hpp"
#include "yuv_stream.hpp"

std::vector<std::vector<double>> dctMatrix = create_2d_vector<double>(8, 8);
std::vector<std::vector<double>> dctTranspose = create_2d_vector<double>(8, 8);

/* Last frame's plane values. */
std::vector<std::vector<unsigned char>> lastYPlane;
std::vector<std::vector<unsigned char>> lastCbPlane;
std::vector<std::vector<unsigned char>> lastCrPlane;


/*
    delta_decode function - Decodes the mixed Improved Row-Delta compression (on values) 
                             and Variable-Length RLE (on repetition counts) into a padded
                             colour plane.
    @param stream - InputBitStream&
    @param height - int
    @param width - int
*/
void delta_decode (InputBitStream &stream, int height, int width, std::vector<std::vector<std::int16_t>> &data)
{
    std::int16_t last = 0;
    std::int16_t index = 0;

    /* Reconstruct the padded image with 8 by 8 matrices. */
    for (auto block_y = 0; block_y < height; block_y += 8)
    {
        for (auto block_x = 0; block_x < width; block_x += 8)
        {
            /* Follow the traversal order. */
            while (index < 64)
            {
                int i = travel_order.at(index).first + block_y;
                int j = travel_order.at(index).second + block_x;
        
                /* Encode first value of every block in 16 bits. */
                if (index < 1)
                {
                    data.at(i).at(j) = stream.read_u16();
                    last = data.at(i).at(j);
                    index++;
                    continue;
                }

                unsigned int bit = stream.read_bit();
                
                /* Case where delta is 0. */
                if (bit == 0)
                    data.at(i).at(j) = last;
                
                /* Decode the delta value.  */
                else
                {
                    bit = stream.read_bit();

                    int counter = 1;
                    while (stream.read_bit() != 0)
                        counter++;

                    if (bit == 0)
                        data.at(i).at(j) = last + counter;
                    else
                        data.at(i).at(j) = last - counter;

                    last = data.at(i).at(j);
                }

                /* RLE begins to decode the repetitions. */
                index++; 
                int bits = 0;
                int reps = 0;

                /* Count how many bits of information are used. */
                while (stream.read_bit())
                    bits++;
                
                /* If no bits used, no repetitions. */
                if (bits == 0)
                    continue;

                /* If one bit used, must be one or two repetitions. */
                else if (bits == 1) 
                {
                    if(!stream.read_bit())
                        reps = 1;
                    else
                        reps = 2;
                }   

                /* Else read in number of bits to find amount of repetitions. */
                else
                {
                    reps = stream.read_bits(bits - 1);
                    
                    if(reps == 0)
                        reps |= 1 << (bits - 1);
                }
                
                /* Add the last value a given amount of times in traversal order. */
                while(reps-- > 0)
                {
                    i = travel_order.at(index).first + block_y;
                    j = travel_order.at(index).second + block_x;
                    data.at(i).at(j) = last;
                    index++;
                }
            }
            index = 0;
        }
    }
}


/*
    recover_mB function - Takes a macroblock of 4 8 by 8 Y blocks and one each of Cb and Cr 8 by 8
                            blocks and recovers the original YCbCr values through reverse dct and
                            de-quantization. These are then transferred to their actual spots in the
                            colour plane for frame reconstruction.
    @param data - std::vector<std::vector<std::int16_t>>&
    @param recovered - std::vector<std::vector<unsigned char>>&
    @param iFrame - bool    // (If macroblock is an iFrame)
    @param h_size - int
    @param w_size - int
    @param quality - string
    @param plane - string
    @param lastPlane - std::vector<std::vector<unsigned char>>&
    @param mX - int         // (x position on the frame)
    @param mY - int         // (y position on the frame)
    @param mB - mBlock&
*/
void recover_mB (std::vector<std::vector<std::int16_t>> &data, std::vector<std::vector<unsigned char>> &recovered, \
            bool iFrame, int h_size, int w_size, std::string quality, std::string plane, \
            std::vector<std::vector<unsigned char>> &lastPlane, int mX, int mY, mBlock &mB)
{
    auto miniMatrix1 = create_2d_vector<std::int16_t>(8, 8);

    for (auto block_y = 0; block_y < h_size; block_y += 8)
    {
        for (auto block_x = 0; block_x < w_size; block_x += 8)
        {
            /* Fill the 8 by 8 matrix with Y values. */
            for (auto i = 0; i < 8; i++)
            {
                for (auto j = 0; j < 8; j++)
                {
                    miniMatrix1.at(i).at(j) = data.at(block_y + i).at(block_x + j);
                }
            }

            auto resultant = create_2d_vector<double>(8, 8);

            /* Undo Quantization. */
            dequantize(miniMatrix1, plane, quality, resultant);

            /* Undo DCT Matrix Multiplication. */
            resultant = inverseMultiply(dctTranspose, resultant, dctMatrix);
            int counter = 0;

            /* Reconstruct the padded image with 8 by 8 matrices. */
            for (auto i = 0; i < 8; i++)
            {
                for (auto j = 0; j < 8; j++)
                {
                    /* Intra-block clamping. */
                    if (iFrame)
                    {
                        if (resultant.at(i).at(j) + (128) > 255)
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)255;
                        else if (resultant.at(i).at(j) + (128) < 0)
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)0;
                        else
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = ((unsigned char)resultant.at(i).at(j) + (128));
                    }

                    /* Inter-block clamping. */
                    else
                    {
                        int vX = mX + mB.mVector.first;
                        int vY = mY + mB.mVector.second;
                        
                        /* If Cb/Cr plane, adjust to account for halved values. */
                        if(plane == "Cb" || plane == "Cr")
                            vX = mX + (mB.mVector.first / 2), vY = mY + (mB.mVector.second / 2);

                        if (resultant.at(i).at(j) + lastPlane.at(vY + block_y + i).at(vX + block_x + j) > 255)
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)255;
                        else if (resultant.at(i).at(j) + lastPlane.at(vY + block_y + i).at(vX + block_x + j) < 0)
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)0;
                        else
                            recovered.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)(resultant.at(i).at(j) + lastPlane.at(vY + block_y + i).at(vX + block_x + j));
                    }

                    lastPlane.at(mY + block_y + i).at(mX + block_x + j) = recovered.at(mY + block_y + i).at(mX + block_x + j);
                }
            }
        }
    }
}


int main(int argc, char** argv)
{    
    InputBitStream input_stream {std::cin};

    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};
    u8 qualint {input_stream.read_byte()};

    YUVStreamWriter writer {std::cout, width, height};

    /* Calculate DCT and DCT Transpose. */
    dct8Calc(dctMatrix);
    transpose(dctMatrix, dctTranspose);

    std::string quality = "";

    /* Interpret quality meta-data byte. */
    if (qualint == 1)
        quality = "low";
    else if (qualint == 2)
        quality = "medium";
    else
        quality = "high";

    /* Padded size for the Y planes. */
    unsigned int h_size = height;
    unsigned int w_size = width;
    unsigned int c_h_size = (height + 1) / 2;
    unsigned int c_w_size = (width + 1) / 2;

    /* Calculates the necessary padding height and width of the Y Plane. */
    if(height % 8 != 0)
        h_size += (8 - (height % 8));
    if(width % 8 != 0)
        w_size += (8 - (width % 8));

    /* Calculates the padded dimensions for the Cb and Cr planes after scaledown. */
    if(c_h_size % 8 != 0)
        c_h_size += (8 - (c_h_size % 8));
    if(c_w_size % 8 != 0)
        c_w_size += (8 - (c_w_size % 8));

    /* Create padded planes for each frame. */
    auto Y = create_2d_vector<unsigned char>(h_size, width);
    auto Cb = create_2d_vector<unsigned char>(c_h_size, c_w_size);
    auto Cr = create_2d_vector<unsigned char>(c_h_size, c_w_size);

    /* Create temporary macroblock data blocks. */
    auto Y_mB = create_2d_vector<std::int16_t>(16, 16);
    auto Cb_mB = create_2d_vector<std::int16_t>(8, 8);
    auto Cr_mB = create_2d_vector<std::int16_t>(8, 8);

    /* Initialize last frame's planes' values. */
    lastYPlane = create_2d_vector<unsigned char>(h_size, w_size);
    lastCbPlane = create_2d_vector<unsigned char>(c_h_size, c_w_size);
    lastCrPlane = create_2d_vector<unsigned char>(c_h_size, c_w_size);

    int frameType = 0;

    while (input_stream.read_byte())
    {
        YUVFrame420& frame = writer.frame();

        bool iFrame;
        
        /* Decodes and I-Frame ever 3 frames. */
        if(frameType % 3 == 0)
            iFrame = true;

        /* Decodes macroblocks and stores them in the planes. */
        /* Header Format: 1/0 I/P-Frame -- Motion Vector -- 4 Y blocks -- 1 Cb block -- 1 Cr block */
        for (auto mY = 0; mY < h_size; mY += 16)
        {
            for (auto mX = 0; mX < w_size; mX += 16)
            {
                mBlock mB;
                mB.mVector.first = 0;
                mB.mVector.second = 0;
                
                /* Check if I-Frame or P-Frame. */
                if(input_stream.read_byte() == 1)
                    iFrame = true;
                else
                {
                    /* Read in P-Frame meta-data. */
                    iFrame = false;
                    mB.mVector.first = input_stream.read_byte();
                    mB.mVector.second = input_stream.read_byte();
                }
                
                /* Decode the macroblock. */
                delta_decode(input_stream, 16, 16, Y_mB);
                delta_decode(input_stream, 8, 8, Cb_mB);
                delta_decode(input_stream, 8, 8, Cr_mB);
                
                int cp_mX = (mX + 1) / 2;
                int cp_mY = (mY + 1) / 2;

                /* Undo quantization, temporal compression, DCT Transformation, etc. */
                recover_mB(Y_mB, Y, iFrame, 16, 16, quality, "Y", lastYPlane, mX, mY, mB);
                recover_mB(Cb_mB, Cb, iFrame, 8, 8, quality, "Cb", lastCbPlane, cp_mX, cp_mY, mB);
                recover_mB(Cr_mB, Cr, iFrame, 8, 8, quality, "Cr", lastCrPlane, cp_mX, cp_mY, mB);
            }
        }

        /* Reassemble the frame for outputting.  */
        for (unsigned int y = 0; y < height; y++)
        {
            for (unsigned int x = 0; x < width; x++)
            {
                try
                {
                    frame.Y(x, y) = Y.at(y).at(x);

                    if(y < height / 2 && x < width / 2)
                    {
                        frame.Cb(x, y) = Cb.at(y).at(x);
                        frame.Cr(x, y) = Cr.at(y).at(x);
                    }
                }
                catch (const std::exception &e)
                {
                    std::cerr << "height : " << y << "and width: " << x << std::endl;
                }
            }
        }

        writer.write_frame();

        frameType++;
    }


    return 0;
}