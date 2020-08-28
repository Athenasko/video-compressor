/* yuv_stream.hpp

   Reader/writer for raw frames of YUV data with 4:2:0 subsampling
   Streams of this type can be played or encoded/decoded with ffplay/ffmpeg
   with parameters like

   ./some_program_generating_this-data | ffplay -f rawvideo -framerate 24 -pixel_format yuv420p -video_size 854x480

   Note that the stream does not have the resolution encoded into it in any way (so that must be
   provided externally somehow).

   Recall that YUV is a misnomer (it's actually YCbCr since YUV is an analog encoding)

   B. Bird - 04/07/2020
*/

#ifndef YUV_STREAM_HPP
#define YUV_STREAM_HPP

#include <iostream>
#include <vector>
#include <array>

class YUVStreamReader;
class YUVStreamWriter;

template<int chroma_ratio_x, int chroma_ratio_y> //chroma_ratio is the number of Y pixels per chroma pixel in each direction (for 4:2:0 or 4:1:1 these numbers are both 2)
class YUVFrame{
public:
    YUVFrame(unsigned int width, unsigned int height): width{width}, height{height} {
        assert(width%chroma_ratio_x == 0);
        assert(height%chroma_ratio_y == 0);
        Y_data.resize(width*height);
        Cb_data.resize(width*height/(chroma_ratio_x*chroma_ratio_y));
        Cr_data.resize(width*height/(chroma_ratio_x*chroma_ratio_y));
    }
    unsigned char& Y(unsigned int x, unsigned int y){
        if (y >= height)
            y = height-1;
        if (x >= width)
            x = width-1;
        return Y_data.at(y*width+x);
    }
    //Note that the coordinate systems for Cb and Cr are distinct from Y (e.g. in 4:2:0, the chroma values for Y pixel (10,10) are at Cb/Cr coordinates (5,5))
    unsigned char& Cb(unsigned int x, unsigned int y){
        if (y >= height/chroma_ratio_y)
            y = height/chroma_ratio_y-1;
        if (x >= width/chroma_ratio_x)
            x = width/chroma_ratio_x-1;
        return Cb_data.at(y*width/chroma_ratio_x + x);
    }
    unsigned char& Cr(unsigned int x, unsigned int y){
        if (y >= height/chroma_ratio_y)
            y = height/chroma_ratio_y-1;
        if (x >= width/chroma_ratio_x)
            x = width/chroma_ratio_x-1;
        return Cr_data.at(y*width/chroma_ratio_x + x);
    }
private:
    unsigned int width, height;
    std::vector<unsigned char> Y_data, Cb_data, Cr_data;
    friend class YUVStreamReader;
    friend class YUVStreamWriter;
};

using YUVFrame420 = YUVFrame<2,2>;

class YUVStreamReader{
public:
    YUVStreamReader(std::istream& stream, unsigned int width, unsigned int height): input_stream{stream}, active_frame{width,height}, frame_counter{0U} {
        
    }

    YUVFrame420& frame(){
        return active_frame;
    }

    bool read_next_frame(){
        frame_counter++;
        return read_into_array(active_frame.Y_data) && 
               read_into_array(active_frame.Cb_data) && 
               read_into_array(active_frame.Cr_data);
    }

private:
    bool read_into_array(std::vector<unsigned char>& A){
        //input_stream.read(&A.at(0),A.size());
        //return (bool)input_stream;
        char c;
        unsigned int i = 0;
        while(i < A.size() && input_stream.get(c))
            A.at(i++) = (unsigned char)c;
        return i == A.size();
    }
    std::istream& input_stream;
    YUVFrame420 active_frame;
    unsigned int frame_counter;
};

class YUVStreamWriter{
public:
    YUVStreamWriter(std::ostream& stream, unsigned int width, unsigned int height): output_stream{stream}, active_frame{width,height}, frame_counter{0U} {
        
    }

    YUVFrame420& frame(){
        return active_frame;
    }

    bool write_frame(){
        return write_from_array(active_frame.Y_data) && 
               write_from_array(active_frame.Cb_data) && 
               write_from_array(active_frame.Cr_data);
    }

private:
    bool write_from_array(std::vector<unsigned char> const& A){
        //output_stream.write(&A.at(0), A.size());
        //return output_stream
        for(unsigned char c: A)
            output_stream.put(c);
        return true;
    }
    std::ostream& output_stream;
    YUVFrame420 active_frame;
    unsigned int frame_counter;
};

#endif