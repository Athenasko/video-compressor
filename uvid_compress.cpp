/* uvid_compress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5

   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 

     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>

   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.

   B. Bird - 07/15/2020
*/

#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include "output_stream.hpp"
#include "uvid_common.hpp"
#include "yuv_stream.hpp"


std::vector<std::vector<double>> dctMatrix = create_2d_vector<double>(8, 8);
std::vector<std::vector<double>> dctTranspose = create_2d_vector<double>(8, 8);

std::vector<std::vector<unsigned char>> lastYPlane;
std::vector<std::vector<unsigned char>> lastCbPlane;
std::vector<std::vector<unsigned char>> lastCrPlane;


/*
    delta_encode function - Encodes an 8 by 8 block using Improved Row-Delta compression
                             to compress the values, and Variable-Length RLE to compress
                             any repetitions.
    @param stream - OutputBitStream&
    @param data - std::vector<std::vector<double>>&
*/
void delta_encode (OutputBitStream &stream, std::vector<std::vector<std::int16_t>> &data)
{
    std::int16_t last = 0;
    std::int16_t delta = 0;
    int index = 0;

    while(index < 64)
    {
        int reps = 0;
        std::int16_t current = data.at(travel_order.at(index).first).at(travel_order.at(index).second);
        
        /* Pushes first two values. */
        if (index < 1)
        {
            index++;
            stream.push_u16(current);
            last = current;
            continue;
        }

        while(index < 63 && data.at(travel_order.at(index).first).at(travel_order.at(index).second) == \
                data.at(travel_order.at(index + 1).first).at(travel_order.at(index + 1).second))
        {
            index++;
            reps++;
        }

        delta = current - last;
        //std::cerr << "index: " << index << "  delta: " << delta << std::endl;
        if (delta == 0)
            stream.push_bit(0);

        else
        {
            stream.push_bit(1);

            if (delta > 0)
            {
                stream.push_bit(0);
                
                for (int i = delta; i > 1; i--)
                    stream.push_bit(1);
            }
            else
            {
                stream.push_bit(1);

                if (delta != -1)
                    for (int i = delta; i < -1; i++)
                        stream.push_bit(1);
            }
            stream.push_bit(0);
        }

        last = current;

        /* Variable Length RLE for Repetitions. */

        /* If not repetitions push 0. */
        if(reps == 0)
        {
            stream.push_bit(0);
            index++;
            continue;
        }

        /* Edge cases for RLE repetitions. */
        else if (reps == 1 || reps == 2)
        {
            stream.push_bit(1);

            if(reps == 1)
            {
                stream.push_bit(0);
                stream.push_bit(0);
            }

            else
            {
                stream.push_bit(0);
                stream.push_bit(1);
            }
        }

        /* Push the number of repetitions using Variable-Length RLE. */
        else
        {
            int bits = ceil(log2(reps));

            if (reps && (reps - 1))
                bits++;

            for (int i = 0; i < bits; i++)
                stream.push_bit(1);

            stream.push_bit(0);

            stream.push_bits(reps, bits - 1);
        }
        index++;
    }
}


/*
    stream_out_mb function - Takes a macroblock and puts each 8 by 8 block through
                             dct, quantization, and finally delta compresses it. Each block
                             is also decompressed and saved in a past frame.
    @param output_stream - OutputBitStream&
    @param mB - mBlock&
    @param iFrame - bool
    @param quality - string
    @param mX - int
    @param mY - int
*/
void stream_out_mb(OutputBitStream &output_stream, mBlock &mB, bool iFrame, std::string quality, int mX, int mY)
{
    auto miniMatrix1 = create_2d_vector<std::int16_t>(8, 8);
    auto miniMatrix2 = create_2d_vector<std::int16_t>(8, 8);
    auto resultant = create_2d_vector<double>(8, 8);
    auto resultant2 = create_2d_vector<double>(8, 8);
    auto rounded_res = create_2d_vector<std::int16_t>(8, 8);

    /* Intra-block Encoding flag. */
    if (iFrame)
        output_stream.push_byte(1);
    else
    {
        /* Inter-block Encoding flag. */
        output_stream.push_byte(0);

        /* Macroblock vector deltas. */
        output_stream.push_byte(mB.mVector.first);
        output_stream.push_byte(mB.mVector.second);
    }

    /* Creating a Y Macroblock chunk. */
    for (auto block_y = 0; block_y < 16; block_y += 8)
    {
        for (auto block_x = 0; block_x < 16; block_x += 8)
        {
            for (auto i = 0; i < 8; i++)
            {
                for (auto j = 0; j < 8; j++)
                {
                    /* Intra-block Y encoding. */
                    if (iFrame)
                        miniMatrix1.at(i).at(j) = (std::int16_t)mB.Y.at(i + block_y).at(j + block_x) - 128;
                    /* Inter-block Y encoding. */
                    else
                        miniMatrix1.at(i).at(j) = (std::int16_t)mB.Y.at(i + block_y).at(j + block_x) - (std::int16_t)lastYPlane.at(mY + mB.mVector.second + block_y + i).at(mX + mB.mVector.first + block_x + j);
                }
            }

            /* Transform each 8 by 8 Y matrix with DCT. */
            resultant = multiply(dctMatrix, miniMatrix1, dctTranspose);

            /* Quantize each 8 by 8 Y Matrix. */
            quantize(resultant, "Y", quality, rounded_res);

            /* Stream out each Y block in delta-encoded traversal order. */
            delta_encode(output_stream, rounded_res);

            /* Undo Quantization. */
            dequantize(rounded_res, "Y", quality, resultant);

            /* Undo DCT Matrix Multiplication. */
            resultant = inverseMultiply(dctTranspose, resultant, dctMatrix);

            /* Reconstruct the Y padded image with 8 by 8 matrices. */
            for (auto i = 0; i < 8; i++)
            {
                for (auto j = 0; j < 8; j++)
                {
                    /* Intra-block Y clamping. */
                    if (iFrame)
                    {
                        if (resultant.at(i).at(j) + (128) > 255)
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)255;
                        else if (resultant.at(i).at(j) + (128) < 0)
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)0;
                        else
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)(resultant.at(i).at(j) + (128));
                    }
                    /* Inter-block Y clamping. */
                    else
                    {
                        if (resultant.at(i).at(j) + lastYPlane.at(mY + mB.mVector.second + block_y + i).at(mX + mB.mVector.first + block_x + j) > 255)
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)255;
                        else if (resultant.at(i).at(j) + lastYPlane.at(mY + mB.mVector.second + block_y + i).at(mX + mB.mVector.first + block_x + j) < 0)
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)0;
                        else
                            lastYPlane.at(mY + block_y + i).at(mX + block_x + j) = (unsigned char)(resultant.at(i).at(j) + lastYPlane.at(mY + mB.mVector.second + block_y + i).at(mX + mB.mVector.first + block_x + j));
                    }
                }
            }
        }
    }


    /* Halve the values for use on Cb and Cr planes. */
    int cpX = (mX + 1) / 2;
    int cpY = (mY + 1) / 2;    
    int mvX, mvY;

    if(iFrame)
        mvX = cpX, mvY = cpY;
    else
    {
        mvX = (mY + 1 + mB.mVector.second) / 2;
        mvY = (mY + 1 + mB.mVector.second) / 2;
    }
    
    /* Creating the Cb and Cr macroblock chunks. */
    for (auto i = 0; i < 8; i++)
    {
        for (auto j = 0; j < 8; j++)
        {
            /* Intra-block Cb/Cr encoding. */
            if (iFrame)
            {
                miniMatrix1.at(i).at(j) = (std::int16_t)mB.Cb.at(i).at(j) - 128;
                miniMatrix2.at(i).at(j) = (std::int16_t)mB.Cr.at(i).at(j) - 128;
            }

            /* Inter-block Cb/Cr encoding. */
            else
            {
                miniMatrix1.at(i).at(j) = (std::int16_t)mB.Cb.at(i).at(j) - (std::int16_t)lastCbPlane.at(mvY + i).at(mvX + j);
                miniMatrix2.at(i).at(j) = (std::int16_t)mB.Cr.at(i).at(j) - (std::int16_t)lastCrPlane.at(mvY + i).at(mvX + j);
            }
        }
    }
    /* Transform the 8 by 8 Cb and Cr matrices with DCT. */
    resultant = multiply(dctMatrix, miniMatrix1, dctTranspose);
    resultant2 = multiply(dctMatrix, miniMatrix2, dctTranspose);

    /* Quantize the 8 by 8 Cb Matrix. */
    quantize(resultant, "Cb", quality, miniMatrix1);
    quantize(resultant2, "Cr", quality, miniMatrix2);

    /* Stream out the Cb/Cr blocks in delta-encoded traversal order. */
    delta_encode(output_stream, miniMatrix1);
    delta_encode(output_stream, miniMatrix2);

    /* Undo Quantization. */
    dequantize(miniMatrix1, "Cb", quality, resultant);
    dequantize(miniMatrix2, "Cr", quality, resultant2);


    /* Undo DCT Matrix Multiplication. */
    resultant = inverseMultiply(dctTranspose, resultant, dctMatrix);
    resultant2 = inverseMultiply(dctTranspose, resultant2, dctMatrix);

    /* Reconstruct the Cb padded image with 8 by 8 matrices. */
    for (auto i = 0; i < 8; i++)
    {
        for (auto j = 0; j < 8; j++)
        {
            /* Intra-block encoding (I-Frame). */
            if (iFrame)
            {
                /* Cb clamping. */
                if (resultant.at(i).at(j) + (128) > 255)
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)255;
                else if (resultant.at(i).at(j) + (128) < 0)
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)0;
                else
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)(resultant.at(i).at(j) + (128));
                /* Cr clamping. */
                if (resultant2.at(i).at(j) + (128) > 255)
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)255;
                else if (resultant2.at(i).at(j) + (128) < 0)
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)0;
                else
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)(resultant2.at(i).at(j) + (128));
            }

            /* Inter-block encoding (P-Frame). */
            else
            {
                /* Cb clamping. */
                if (resultant.at(i).at(j) + lastCbPlane.at(mvY + i).at(mvX + j) > 255)
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)255;
                else if (resultant.at(i).at(j) + lastCbPlane.at(mvY + i).at(mvX + j) < 0)
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)0;
                else
                    lastCbPlane.at(cpY + i).at(cpX + j) = (unsigned char)(resultant.at(i).at(j) + lastCbPlane.at(mvY + i).at(mvX + j));

                /* Cr clamping. */
                if (resultant2.at(i).at(j) + lastCrPlane.at(mvY + i).at(mvX + j) > 255)
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)255;
                else if (resultant2.at(i).at(j) + lastCrPlane.at(mvY + i).at(mvX + j) < 0)
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)0;
                else
                    lastCrPlane.at(cpY + i).at(cpX + j) = (unsigned char)(resultant2.at(i).at(j) + lastCrPlane.at(mvY + i).at(mvX + j));
            }
        }
    }
}


/*
    vectorFinder function - Checks a macroblock's Y plane data if it would benefit from P-Frames.
                             Uses AAD and PSNR to check.
    @param mB - mBlock&
    @param mX - int
    @param mY - int
    @param height - int
    @param width - int
    @param quality - std::string
*/
bool vectorFinder (mBlock &mB, int mX, int mY, int height, int width, std::string quality)
{
    double best_psnr = 30;

    int startingX, startingY, endingX, endingY;
    int max_i;

    startingX = -8, startingY = -8, endingX = 8, endingY = 8, max_i = 64;

    if(mX + startingX < 0)
        startingX = 0;
    if(mY + startingY < 0)
        startingY = 0;
    if(mX + endingX > width - 16)
        endingX = width - mX - 16;
    if(mY + endingY > height - 16)
        endingY = height - mY - 16;   

    /* Initialized search values. */
    double aad = 500, best_aad = 500;
    double mse = 0;

    /* Local search with a distance of 8 pixels. */
    for(auto search_Y = startingY; search_Y < endingY; search_Y++)
    {
        for(auto search_X = startingX; search_X < endingX; search_X++)
        {
            /* Ensure the motion vector is actually searched for and isn't overrided. */
            if(search_Y == 0 && search_X == 0)
                continue;

            /* Calculates the AAD value of a given 16 by 16 area. */
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 16; j++)
                {
                    double p0 = mB.Y.at(i).at(j);
                    double p1 = lastYPlane.at(mY + i).at(mX + j);

                    double diff = abs(p0 - p1);
                    aad += diff;
                    mse += diff * diff;
                }
            }
            mse /= 256;

            /* If the AAD is the new lowest calculate the PSNR. */
            if (aad < best_aad)
            {
                /* Peak-Signal-to-Noise-Ratio */
                double psnr = (20 * log10((max_i) - 10 * log10(mse)));

                /* If the PSNR is the new best, save the motion vector. */
                if (psnr > best_psnr)
                {
                    mB.mVector.first = search_X, mB.mVector.second = search_Y;
                    best_psnr = psnr;
                }
            }
        }
    }

    /* If a motion vector is found, save it. */
    if(best_psnr > 30)
    {
        //std::cerr << "PFrame" << std::endl;
        return false;
    }

    /* If no motion vector found, encode using intra-block encoding (I-Frame). */
    return true;
}


/*
    main function - Takes in a stream of video data and compresses it.
    @param argc - int
    @param argv - char**
*/
int main(int argc, char** argv)
{
    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};

    /* Calculate DCT and DCT Transpose. */
    dct8Calc(dctMatrix);
    transpose(dctMatrix, dctTranspose);

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

    /* Convert quality into a one byte integer to pass along. */
    std::uint8_t qualint;
    if (quality == "low")
        qualint = 1;
    else if (quality == "medium")
        qualint = 2;
    else
        qualint = 3;

    /* Push meta-data. */
    output_stream.push_u32(height);
    output_stream.push_u32(width);
    output_stream.push_byte(qualint);

    /* Create planes to store last frame's plane values. */
    lastYPlane = create_2d_vector<unsigned char>(h_size, w_size);
    lastCbPlane = create_2d_vector<unsigned char>(c_h_size, c_w_size);
    lastCrPlane = create_2d_vector<unsigned char>(c_h_size, c_w_size);

    int frameType = 0;

    while (reader.read_next_frame())
    {
        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        YUVFrame420& frame = reader.frame();

        /* 
            Creates Macroblock chunks.
            Header Format: 1/0 I/P-Frame -- Motion Vector -- 4 Y blocks -- 1 Cb block -- 1 Cr block 
        */
        for (auto mY = 0; mY < h_size; mY += 16)
        {
            for (auto mX = 0; mX < w_size; mX += 16)
            {
                mBlock mB;
                /* Initialize the motion vector to point to the start of the macroblock. */
                mB.mVector.first = 0;
                mB.mVector.second = 0;

                for (auto y = 0; y < 16; y++)
                {
                    for (auto x = 0; x < 16; x++)
                    {
                        mB.Y.at(y).at(x) = frame.Y(x + mX, y + mY);

                        /* Fill the Cb and Cr blocks. */
                        if (y < 8 && x < 8)
                        {
                            mB.Cb.at(y).at(x) = frame.Cb(((mX + 1) / 2) + x, ((mY + 1) / 2) + y);
                            mB.Cr.at(y).at(x) = frame.Cr(((mX + 1) / 2) + x, ((mY + 1) / 2) + y);
                        }
                    }
                }

                /* Dedicated I-Frame. (Every 3). */
                if (frameType % 3 == 0)
                {
                    /* Send Macroblock to Stream_out. */
                    stream_out_mb(output_stream, mB, true, quality, mX, mY);
                }

                /* Decide if I-frame or P-frame works better for the Macroblock. */
                else
                {
                    bool iBlock = false;

                    /* Uses local search to find good motion vector match. */
                    iBlock = vectorFinder(mB, mX, mY, h_size, w_size, quality);

                    stream_out_mb(output_stream, mB, iBlock, quality, mX, mY);
                }
            }
        }

        frameType++;
    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();

    return 0;
}