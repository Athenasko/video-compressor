----------------------------------
README SENG480B Assignment 5

Coded by: Austin Lee - V00878834
Starter code provided by: B. Bird

Date: 08/10/2020
----------------------------------

-----------------------------------------------------------------------------------
-                                     NOTES                                       -
-----------------------------------------------------------------------------------
This was an interesting assignment to take on. Building up to it with Assignment 4
was definitely a 100% great idea because although it takes a little to work out
how to do it, being able to canibalize 80% of the previous assignment's code is 
something that many classes don't do but should since it's good to have a 
functioning starting point and in reality we reuse code a lot. 
Since I was in all project courses, being able to reuse it was also a godsend.

I thought about trying to go for B-Frames but that mind-fuckery (pardon my french)
was too much for my time. On the plus side, Temporal and Motion Compensation make
sense and overall, implementing them (while tedious sometimes to find the bugs), 
wasn't too bad overall. Definitely (although very time consuming) enjoyed this
assignment way more than A2 or A3 and although I didn't quite manage to get real-
time streaming on medium or high like I wanted, I feel that I did well (at least
personally) on this assignment.
-----------------------------------------------------------------------------------
-                                 IMPLEMENTATION                                  -
-----------------------------------------------------------------------------------
    This implementation of a video compressor contains the following pipeline:

    The starter code, provided by B. Bird, takes raw video in YUV420p format and 
    splits each fram into the YUVFrame format.
    -> The program send metadata of 3 pieces of information. The original height, 
        original width, and the quantization quality.
    -> For each frame:
        --> A single byte is sent to denote another frame exists.
        --> The frame is split into macroblocks of 4 Y 8 by 8 blocks and one 8 by 8
             block each of Cb and Cr.
        --> Each macroblock is encoded with a header containing a one byte flag
             denoting it as an intra-encoded macroblock.
        --> The macroblock is then sent to the to be encoded in 8 by 8 chunks where
             as the chunks are split, the values are also subtracted 128, similarly
             done by the JPEG to shift the average to 0.
        --> Each 8 by 8 block for the plane is transformed by the DCT, of which
             the DCT and DCT Transpose matrices are generated once at the start.
        --> Each 8 by 8 block also quantized by the JPEG schematic derived 
             luminance and chrominance matrices. Depending on the quality, these
             matrices are more than doubled (low), slightly increased (medium) or
             slightly decreased (high). These quantized blocks are also rounded.
        --> The 8 by 8 block is then delta compressed using Improved Row-Delta
             compression based on slides provided by B. Bird, and encoded in
             traversal order. After a value is output in delta compression, if 
             there are repetitions, the repetitions are encoded using Variable-
             Length RLE.
        --> Each 8 by 8 is also decoded and saved into its location in plane to
             act as a "past frame" reference for P-Frames.
        --> This repeats for each other plane in the frame.
    -> The compressor uses P-Frames for up to the next two frames.
        --> For each macroblock, the compressor checks if it would benefit 
             from temporal compression using local search of 16 pixels for high 
             quality or 8 for medium or low, in either direction checking for a 
             Peak-Signal-to-Noise-Ratio of at least 30.
             The algorithm checks also Average Absolute Difference initially,
             and only the top left 8 by 8 block of each macroblock to improve
             runtime. 
                ---> If temporal compression would help, it encodes the 1 byte flag
                      as a 0.
                ---> It then calculates a motion vector based on a local search
                      algorithm and sends two bytes containing the delta values of
                      the x and y for the calculated motion vector.
        --> The macroblock is then encoded similarly as to the steps in intra-block
             encoding. If inter-block encoding, instead of subtracting 128 like JPEG,
             motion compression subtracts the best last frame instead.
    
    These steps in reverse are then followed to reconstruct and decompress the 
    original image.
-----------------------------------------------------------------------------------
-                              WHAT WAS ACCOMPLISHED                              -
-----------------------------------------------------------------------------------
    The following has been accomplished:

    1. Frames are read in and split into their respective YCbCr planes.
    
    2. Each frame is decomposed into macroblocks of size 16 x 16 for Y and 8 x 8
        for Cb and Cr.
        
    3. Each macroblock is transformed using the DCT and quantized before being
        encoded using a mix of Improved-Row Delta and Variable-Length Run Length.

    4. For P-Frames, each macroblock is checked using Average Absolute Difference 
        and confirmed with Peak-Signal-to-Noise-Ratio. If determined to be good for
        inter-block encoding, the best motion vector in an 8 pixel square is added.
    
    5. All together, a compression ratio of roughly 6 - 12 on low is achieved
        depending on the video. By manipulating the chrominance and luminance
        matrices and quality settings especially for low, the low quality is quite
        a bit worse than medium or high due to the matrices being more harsh.
        As such, low quality is not great, but it is watchable, just with lots of
        artifacting.

    6. The stream-ability is not quite real-time. For 352 x 288 resolutions, it
        can go around 20+ frames a second on high, and kinda real-time on low, 
        though this depends on the video.
-----------------------------------------------------------------------------------
-                                      BUGS                                       -
-----------------------------------------------------------------------------------
    This code has no known bugs.
-----------------------------------------------------------------------------------
-                                   REFERENCES                                    -
-----------------------------------------------------------------------------------
    1. Basic Pseudocode provided by B. Bird.
    2. Cooper Mountford - Student in the class. Discorded a lot (verbally only).
    3. https://media.xiph.org/video/derf/ - Videos for testing.
    4. https://github.com/ruofeidu/ImageQualityCompare/blob/master/Compare - Used 
            to understand and guide AAD and MSE calculations.
            