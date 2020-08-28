/* uvg_common.hpp

   B. Bird - 07/02/2020
*/

#ifndef UVG_COMMON_HPP
#define UVG_COMMON_HPP

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

double luminance_matrix [8][8] = {{16, 11, 10, 16, 24, 40, 51, 61},\
                               {12, 12, 14, 19, 26, 58, 60, 55},\
                               {14, 13, 16, 24, 40, 57, 69, 56},\
                               {14, 17, 22, 29, 51, 87, 80, 62},\
                               {18, 22, 37, 56, 68, 109, 103, 77},\
                               {24, 35, 55, 64, 81, 104, 113, 92},\
                               {49, 64, 78, 87, 103, 121, 120, 101},\
                               {72, 92, 95, 98, 112, 100, 103, 99}};

double low_luminance_matrix [8][8] = {{25, 20, 17, 25, 37, 50, 65, 75},\
                               {23, 23, 25, 30, 37, 71, 75, 80},\
                               {25, 35, 25, 37, 50, 80, 85, 75},\
                               {25, 24, 28, 38, 80, 92, 90, 75},\
                               {24, 28, 45, 80, 85, 114, 108, 90},\
                               {32, 42, 80, 85, 81, 104, 130, 100},\
                               {49, 90, 95, 100, 108, 130, 125, 110},\
                               {90, 110, 120, 106, 130, 110, 115, 105}};


double chrominance_matrix [8][8] = {{17, 18, 24, 47, 99, 99, 99, 99},\
                                 {18, 21, 26, 66, 99, 99, 99, 99},\
                                 {24, 26, 56, 99, 99, 99, 99, 99},\
                                 {47, 66, 99, 99, 99, 99, 99, 99},\
                                 {99, 99, 99, 99, 99, 99, 99, 99},\
                                 {99, 99, 99, 99, 99, 99, 99, 99},\
                                 {99, 99, 99, 99, 99, 99, 99, 99},\
                                 {99, 99, 99, 99, 99, 99, 99, 99}};
                                 
double low_chrominance_matrix [8][8] = {{25, 26, 32, 53, 105, 120, 120, 120},\
                                 {26, 28, 35, 72, 105, 120, 120, 120},\
                                 {31, 35, 62, 105, 120, 120, 120, 120},\
                                 {55, 73, 105, 120, 120, 120, 120, 120},\
                                 {105, 120, 120, 120, 120, 120, 120, 120},\
                                 {120, 120, 120, 120, 120, 120, 120, 120},\
                                 {120, 120, 120, 120, 120, 120, 120, 120},\
                                 {120, 120, 120, 120, 120, 120, 120, 120}};

//Convenience function to wrap around the nasty notation for 2d vectors
template<typename T>
std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner)
{
    std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
    return V;
}

std::vector<std::pair<int, int>> travel_order = {
    {0,0}, {0,1}, {1,0}, {2,0}, {1,1}, {0,2}, {0,3}, {1,2}, 
    {2,1}, {3,0}, {4,0}, {3,1}, {2,2}, {1,3}, {0,4}, {0,5}, 
    {1,4}, {2,3}, {3,2}, {4,1}, {5,0}, {6,0}, {5,1}, {4,2}, 
    {3,3}, {2,4}, {1,5}, {0,6}, {0,7}, {1,6}, {2,5}, {3,4}, 
    {4,3}, {5,2}, {6,1}, {7,0}, {7,1}, {6,2}, {5,3}, {4,4}, 
    {3,5}, {2,6}, {1,7}, {2,7}, {3,6}, {4,5}, {5,4}, {6,3}, 
    {7,2}, {7,3}, {6,4}, {5,5}, {4,6}, {3,7}, {4,7}, {5,6}, 
    {6,5}, {7,4}, {7,5}, {6,6}, {5,7}, {6,7}, {7,6}, {7,7}
};

struct mBlock {
    std::pair <std::int8_t, std::int8_t> mVector;
    std::vector<std::vector<unsigned char>> Y = create_2d_vector<unsigned char>(16, 16);
    std::vector<std::vector<unsigned char>> Cb = create_2d_vector<unsigned char>(8, 8);
    std::vector<std::vector<unsigned char>> Cr = create_2d_vector<unsigned char>(8, 8);
};


/* 
    dct8Calc function - Computes and saves an 8 by 8 DCT matrix.
    param @dctMatrix - double &
*/
void dct8Calc (std::vector<std::vector<double>>(&dctMatrix))
{
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            if (i == 0)
                dctMatrix.at(i).at(j) = sqrt((double)1/(double)8);
            else
                dctMatrix.at(i).at(j) = (sqrt((double)1/(double)4)) * cos(((double)(2*j + 1) * (double)i * M_PI)/(16));
        }
    }
}

/* 
    transpose function - Computes and saves the transpose of a matrix.
    param @dctMatrix - double &
*/
void transpose (std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> (&resultant))
{
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            resultant.at(i).at(j) = matrix.at(j).at(i);
        }
    }
}


/* 
    multiply function - Computes matrix multiplication for two equal sized 8 by 8 matrices. Must be
                        non-DCT matrix as left and DCT matrix as right due to typing.
    @param leftMat - std::vector<std::vector<unsigned char>>
    @param rightMat - std::vector<std::vector<double>> 
*/
std::vector<std::vector<double>> multiply (std::vector<std::vector<double>> &dct, std::vector<std::vector<std::int16_t>> &dataMat,\
                                            std::vector<std::vector<double>> &dctTp)
{
    /* Left Multiply the DCT Matrix. */
    auto temp = create_2d_vector<double>(8, 8);
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            temp.at(i).at(j) = 0;
            for(auto k = 0; k < 8; k++)
                temp.at(i).at(j) += dct.at(i).at(k) * (double)dataMat.at(k).at(j);
        }
    }
    
    /* Right Multiply the DCT Transpose. */
    auto resultant = create_2d_vector<double>(8, 8);
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            resultant.at(i).at(j) = 0;
            for(auto k = 0; k < 8; k++)
                resultant.at(i).at(j) += temp.at(i).at(k) * dctTp.at(k).at(j);
        }
    }

    return resultant;
}

/* 
    inverseMultiply function - Computes matrix multiplication for two equal sized 8 by 8 matrices. Must
                               be DCT matrix as left and non-DCT matrix as right due to typing.
    @param leftMat - std::vector<std::vector<double>>
    @param rightMat - std::vector<std::vector<unsigned char>> 
*/
std::vector<std::vector<double>> inverseMultiply (std::vector<std::vector<double>> &dctTp, std::vector<std::vector<double>> &dataMat,\
                                                    std::vector<std::vector<double>> &dct)
{
    /* Left Multiply the DCT Transpose. */
    auto temp = create_2d_vector<double>(8, 8);
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            temp.at(i).at(j) = 0;
            for(auto k = 0; k < 8; k++)
                temp.at(i).at(j) += dctTp.at(i).at(k) * dataMat.at(k).at(j);
        }
    }

    /* Right Multiply the DCT Matrix. */
    auto resultant = create_2d_vector<double>(8, 8);
    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            resultant.at(i).at(j) = 0;
            for(auto k = 0; k < 8; k++)
                resultant.at(i).at(j) += temp.at(i).at(k) * dct.at(k).at(j);
        }
    }

    return resultant;
}


/* 
    quantize function - Computes matrix quantization for color planes. Must pass flag to denote
                        which plane. 0 for Y, 1 for Cb/Cr.
    @param dataMat - std::vector<std::vector<double>>&
    @param plane - std::string
    @param quality - std::string
    @param - resultant - std::vector<std::vector<std::int16_t>>&
*/
void quantize (std::vector<std::vector<double>> &dataMat, std::string plane, \
                std::string quality, std::vector<std::vector<std::int16_t>> &resultant)
{
    // double modifier = 0.5;
    // if(quality == "high")
    //     modifier = 2;
    // else if (quality == "medium")
    //     modifier = 1;

    double modifier = 0.3;
    if(quality == "high")
        modifier = 1.5;
    else if (quality == "medium")
        modifier = 0.7;

    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            /* Case where the YPlane is being quantized. */
            if (plane == "Y")
            {
                if (quality == "low")
                    resultant.at(i).at(j) = round(dataMat.at(i).at(j) / (low_luminance_matrix[i][j]/ modifier));
                else
                    resultant.at(i).at(j) = round(dataMat.at(i).at(j) / (luminance_matrix[i][j]/ modifier));

            }

            /* Case where Cb/Cr Planes are being quantized. */
            else if (plane == "Cb" || plane == "Cr")
            {
                if (quality == "low")
                    resultant.at(i).at(j) = round(dataMat.at(i).at(j) / (low_chrominance_matrix[i][j]/ modifier));
                else
                    resultant.at(i).at(j) = round(dataMat.at(i).at(j) / (chrominance_matrix[i][j]/ modifier));
            }

        }
    }
}


/* 
    dequantize function - Computes matrix dequantization for color planes. Must pass flag to denote
                        which type of plane. 0 for Y, 1 for Cb/Cr.
    @param dataMat - std::vector<std::vector<std::int16_t>>&
    @param plane - std::string
    @param quality - std::string
    @param resultant - std::vector<std::vector<double>>&
*/
void dequantize (std::vector<std::vector<std::int16_t>> &dataMat, std::string plane, \
                    std::string quality, std::vector<std::vector<double>> &resultant)
{
    // double modifier = 0.5;
    // if(quality == "high")
    //     modifier = 2;
    // else if (quality == "medium")
    //     modifier = 1;

    double modifier = 0.3;
    if(quality == "high")
        modifier = 1.5;
    else if (quality == "medium")
        modifier = 0.7;

    for(auto i = 0; i < 8; i++)
    {
        for(auto j = 0; j < 8; j++)
        {
            /* Case where YPlane is being de-quantized. */
            if (plane == "Y")
            {
                if (quality == "low")
                    resultant.at(i).at(j) = dataMat.at(i).at(j) * (low_luminance_matrix[i][j]/ modifier);
                else
                    resultant.at(i).at(j) = dataMat.at(i).at(j) * (luminance_matrix[i][j]/ modifier);
            }            
            /* Case where Cb/Cr Planes are being de-quantized. */
            else if (plane == "Cb" || plane == "Cr")
            {                
                if (quality == "low")
                    resultant.at(i).at(j) = dataMat.at(i).at(j) * (low_chrominance_matrix[i][j]/ modifier);   
                else
                    resultant.at(i).at(j) = dataMat.at(i).at(j) * (chrominance_matrix[i][j]/ modifier);   
            }
        }
    }
}



/*
    print function - Prints a 2D vector of type double.
    @param matrix - std::vector<std::vector<double>>
*/
void print(std::vector<std::vector<double>> matrix)
{
    /* Print 8 by 8 block */
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.at(i).size(); j++)
        {
            //std::cerr << (unsigned int) miniMatrix1.at(i).at(j) << ' ';
            std::cerr << std::setprecision(3) << matrix.at(i).at(j) << ' ';
        }
        std::cerr << std::endl;
    }
    std::cerr << "-------------------------------------------" << std::endl;
}


//The floating point calculations we use while converting between 
//RGB and YCbCr can occasionally yield values slightly out of range
//for an unsigned char (e.g. -1 or 255.9).
//Furthermore, we want to ensure that any conversion uses rounding
//and not truncation (to improve accuracy).
inline unsigned char round_and_clamp_to_char(double v)
{
    //Round to int 
    int i = (int)(v+0.5);
    //Clamp to the range [0,255]
    if (i < 0)
        return 0;
    else if (i > 255)
        return 255;
    return i;
}

/* The exact RGB <-> YCbCr conversion formula used here is the "JPEG style"
   conversion (there is some debate over the best conversion formula) */
struct PixelYCbCr;
struct PixelRGB
{
    unsigned char r, g, b;
    PixelYCbCr to_ycbcr(); //Implementation is below (since the PixelYCbCr type has to exist before we can fully define this function)
};

struct PixelYCbCr
{
    unsigned char Y, Cb, Cr;
    inline PixelRGB to_rgb()
    {
        return {
            round_and_clamp_to_char(Y + 1.402*(Cr-128.0)),
            round_and_clamp_to_char(Y-0.344136*(Cb-128.0)-0.714136*(Cr-128.0)),
            round_and_clamp_to_char(Y+1.772*(Cb-128.0))
        };
    }
};


inline PixelYCbCr PixelRGB::to_ycbcr()
{
    return {
        round_and_clamp_to_char(0.299*r + 0.587*g + 0.114*b),
        round_and_clamp_to_char(128 + -0.168736*r + -0.331264*g + 0.5*b),
        round_and_clamp_to_char(128 + 0.5*r + -0.418688*g + -0.081312*b)
    };
}


#endif